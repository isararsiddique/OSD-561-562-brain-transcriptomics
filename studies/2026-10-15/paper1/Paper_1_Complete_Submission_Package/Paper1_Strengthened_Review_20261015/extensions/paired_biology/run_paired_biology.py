#!/usr/bin/env python3
"""Matched-animal regional interactions; every tested panel is exported.

Gene expression is log2(median-ratio-normalized counts + 1), normalized within
each original deposit. Regional subtraction is within the same Source Name.
Permutation inference requires animal exchangeability within age (and cohort
for the pooled flight contrast). This is not causal identification of flight.
"""
import argparse
import itertools
import json
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / "scripts"))
from run_revision import load_dataset, indices, PANELS, check_input_integrity
from statistics_core import Contrast, animal_permutations, bh, permutation_p

SEED = 56156215


def paired_data(cb_number, hp_number):
    cb, hp = load_dataset(cb_number), load_dataset(hp_number)
    cmeta, hmeta = cb["meta"].copy(), hp["meta"].copy()
    cmeta.index = cmeta.animal.to_numpy()
    hmeta.index = hmeta.animal.to_numpy()
    common = sorted(set(cmeta.index) & set(hmeta.index))
    if not common:
        raise ValueError("No matched source animals")
    for column in ["flight", "cohort", "age"]:
        if not (cmeta.loc[common, column].to_numpy() == hmeta.loc[common, column].to_numpy()).all():
            raise ValueError(f"Cross-region metadata conflict for {column}")
    genes = cb["y"].index.intersection(hp["y"].index)
    c = cb["y"].loc[genes, [f"OSD-{cb_number}:{a}" for a in common]].to_numpy()
    h = hp["y"].loc[genes, [f"OSD-{hp_number}:{a}" for a in common]].to_numpy()
    meta = cmeta.loc[common].copy().reset_index(drop=True)
    meta["cb_accession"], meta["hip_accession"] = f"OSD-{cb_number}", f"OSD-{hp_number}"
    meta["cerebellum_sample"] = cmeta.loc[common, "sample"].to_numpy()
    meta["hippocampus_sample"] = hmeta.loc[common, "sample"].to_numpy()
    meta["location_conflict"] = (cmeta.loc[common, "location_conflict"].to_numpy()
                                  | hmeta.loc[common, "location_conflict"].to_numpy())
    accounting = {"cb_accession": f"OSD-{cb_number}", "hip_accession": f"OSD-{hp_number}",
                  "n_cb_animals": len(cmeta), "n_hip_animals": len(hmeta), "n_pairs": len(common),
                  "unmatched_cb": sorted(set(cmeta.index) - set(hmeta.index)),
                  "unmatched_hip": sorted(set(hmeta.index) - set(cmeta.index)),
                  "n_common_filtered_genes": len(genes)}
    return {"genes": genes, "delta": h - c, "cb": c, "hip": h, "meta": meta, "accounting": accounting}


def design(meta, field, positive, pooled=False):
    labels = (meta[field] == positive).to_numpy(float)
    zs = []
    if meta.age.nunique() > 1:
        zs.append((meta.age == "29").to_numpy(float))
    if pooled:
        zs.append((meta.cohort == "ISS_T").to_numpy(float))
    return labels, (np.column_stack(zs) if zs else None)


def ols_hc3(y, x, coefficient, allow_undefined=False):
    """Matrix OLS and independently implemented HC3 marginal standard errors."""
    n, k = x.shape
    if np.linalg.matrix_rank(x) != k:
        raise ValueError("Non-identifiable design")
    inv = np.linalg.inv(x.T @ x)
    coef = (inv @ x.T @ y.T).T
    residual = y - coef @ x.T
    leverage = np.einsum("ij,ji->i", x, inv @ x.T)
    if (leverage > 1 - 1e-9).any():
        if not allow_undefined:
            raise ValueError("HC3 undefined for unit-leverage animal")
        # A singleton nuisance stratum still permits the contrast coefficient,
        # but its 0/0 HC3 residual must never be silently treated as zero variance.
        missing = np.full(y.shape[0], np.nan)
        return coef[:, coefficient], missing, missing, missing, coef, residual, leverage
    weights = (x @ inv[:, coefficient]) / (1 - leverage)
    se = np.sqrt((residual**2) @ (weights**2))
    effect = coef[:, coefficient]
    t = np.divide(effect, se, out=np.zeros_like(effect), where=se > 0)
    p = 2 * stats.t.sf(abs(t), n - k)
    return effect, se, t, p, coef, residual, leverage


def full_design(meta):
    f = (meta.flight == "FLT").to_numpy(float)
    c = (meta.cohort == "ISS_T").to_numpy(float)
    xs = [np.ones(len(meta))]
    if meta.age.nunique() > 1:
        xs.append((meta.age == "29").to_numpy(float))
    xs.extend([f, c, f*c])
    return np.column_stack(xs)


def bootstrap_panel(y, meta, group, z, count, rng):
    frame = meta.copy()
    frame["label"] = group
    # Resample full animals within actual exposure × collection × age strata.
    strata = [f.index.to_numpy() for _, f in frame.groupby(["label", "cohort", "age"], sort=True)]
    values = []
    for _ in range(count):
        chosen = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
        try:
            values.append(Contrast(y[:, chosen], None if z is None else z[chosen]).fit(group[chosen])[0])
        except ValueError:
            pass
    if len(values) < count * .95:
        raise ValueError("Too many non-identifiable bootstrap designs")
    return np.asarray(values)


def triple_interaction(y, meta, bootstraps, wild_draws, rng):
    """Region × flight × cohort in already paired Hip-minus-Cb expression.

    HC3 studentized restricted-model wild bootstrap, Rademacher animal weights.
    Approximate fixed-design inference; no claim of exact randomization testing.
    """
    x = full_design(meta)
    effect, se, t, p, _, _, _ = ols_hc3(y, x, x.shape[1]-1)
    x0 = x[:, :-1]
    coef0 = np.linalg.lstsq(x0, y.T, rcond=None)[0].T
    fitted0 = coef0 @ x0.T
    r0 = y - fitted0
    h0 = np.einsum("ij,ji->i", x0, np.linalg.pinv(x0))
    adjusted_residual = r0 / np.sqrt(1-h0)
    # Batch draws without constructing gene × animal × all-draw arrays.
    exceeded = np.zeros(y.shape[0], dtype=int)
    for _ in range(wild_draws):
        simulated = fitted0 + adjusted_residual * rng.choice([-1., 1.], len(meta))
        _, _, tstar, _, _, _, _ = ols_hc3(simulated, x, x.shape[1]-1)
        exceeded += abs(tstar) >= abs(t)-1e-12
    wildp = (exceeded+1)/(wild_draws+1)
    groups = [f.index.to_numpy() for _, f in meta.groupby(["flight", "cohort", "age"], sort=True)]
    boots = []
    for _ in range(bootstraps):
        chosen = np.concatenate([rng.choice(s, len(s), replace=True) for s in groups])
        xb = x[chosen]
        if np.linalg.matrix_rank(xb) == x.shape[1]:
            boots.append(np.linalg.lstsq(xb, y[:, chosen].T, rcond=None)[0][-1])
    if len(boots) < bootstraps * .95:
        raise ValueError("Too many non-identifiable triple-interaction bootstrap designs")
    return effect, se, p, wildp, np.asarray(boots), x


def label_tests(name, data, mask, field, positive, annotation, sensitivity, rows, genes, influence, gene_influence,
                bootstraps=2000):
    meta = data["meta"].loc[mask].reset_index(drop=True)
    y = data["delta"][:, mask]
    pooled = name == "RRRM2_region_by_flight_cohort_adjusted"
    g, z = design(meta, field, positive, pooled)
    blocks = tuple(v for v in ["cohort" if pooled else None, "age" if meta.age.nunique()>1 else None] if v)
    # The largest paired label family contains only 105,000 assignments;
    # enumerate it fully rather than sample almost as many Monte Carlo draws.
    perms, method, total = animal_permutations(meta, g, blocks, max_exact=200000, seed=SEED)
    panel_index = indices(data["genes"], annotation)
    names = list(panel_index)
    yp = np.vstack([y[panel_index[p]].mean(axis=0) for p in names])
    observed = Contrast(yp, z).fit(g)[0]
    null = Contrast(yp, z).permutation_betas(perms)
    p = permutation_p(observed[:, None], null, method)
    boots = bootstrap_panel(yp, meta, g, z, bootstraps, np.random.default_rng(SEED))
    # Gene-level t-tests are exploratory; their BH family is distinct from panels.
    x = np.column_stack([np.ones(len(meta)), g]) if z is None else np.column_stack([np.ones(len(meta)), z, g])
    gb, gs, gt, gp, _, _, leverage = ols_hc3(y, x, x.shape[1]-1, allow_undefined=True)
    hc3_status = ("undefined: singleton nuisance stratum / unit leverage"
                  if (leverage > 1 - 1e-9).any() else "estimable; approximate inference")
    genes.append(pd.DataFrame({"contrast": name, "sensitivity": sensitivity,
                               "ensembl_id": data["genes"], "region_interaction_log2": gb,
                               "HC3_se": gs, "p_HC3_approximate": gp, "HC3_status": hc3_status}))
    for i, panel in enumerate(names):
        rows.append({"contrast": name, "sensitivity": sensitivity, "panel": panel, "n_genes": len(panel_index[panel]),
                     "n_paired_animals": len(meta), "n_positive": int(g.sum()), "n_negative": int(len(g)-g.sum()),
                     "n_age_12": int((meta.age=="12").sum()), "n_age_29": int((meta.age=="29").sum()),
                     "n_positive_age_12": int(((g==1)&(meta.age=="12")).sum()),
                     "n_positive_age_29": int(((g==1)&(meta.age=="29")).sum()),
                     "n_negative_age_12": int(((g==0)&(meta.age=="12")).sum()),
                     "n_negative_age_29": int(((g==0)&(meta.age=="29")).sum()),
                     "paired_interaction_log2": observed[i], "bootstrap_lower95": np.quantile(boots[:, i], .025),
                     "bootstrap_upper95": np.quantile(boots[:, i], .975), "p_inference": p[i],
                     "inference_method": method, "n_assignments": perms.shape[1], "total_assignments": total,
                     "permutation_strata": "|".join(blocks), "n_bootstraps": len(boots),
                     "interpretation_status": "exploratory; conditional label exchangeability",
                     "approximate_gene_HC3_status": hc3_status,
                     "definition": "(Hip-Cb) positive group minus (Hip-Cb) negative group"})
        for j, gene in enumerate(panel_index[panel]):
            left = np.delete(panel_index[panel], j)
            effect = Contrast(y[left].mean(axis=0)[None, :], z).fit(g)[0][0]
            gene_influence.append({"contrast":name, "sensitivity":sensitivity, "panel":panel,
                                   "omitted_ensembl_id": data["genes"][gene], "effect_without_gene":effect,
                                   "full_effect":observed[i]})
    for animal in meta.animal:
        keep = (meta.animal != animal).to_numpy()
        effect = Contrast(yp[:, keep], None if z is None else z[keep]).fit(g[keep])[0]
        for i, panel in enumerate(names):
            influence.append({"contrast": name, "sensitivity": sensitivity, "panel": panel,
                              "omitted_animal": animal, "effect_without_animal": effect[i], "full_effect": observed[i]})


def make_figure(results, score_frame, accounting):
    plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "axes.titlesize": 13,
                         "xtick.labelsize": 11, "ytick.labelsize": 11, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(17, 12), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 1.25])
    a, b, c, d = [fig.add_subplot(grid[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
    groups = ["GC LAR", "FLT LAR", "GC ISS-T", "FLT ISS-T"]
    colors = ["#4c78a8", "#f28e2b", "#4c78a8", "#f28e2b"]
    sizes = [len(score_frame[(score_frame.mission=="RRRM-2") & (score_frame.panel=="SRG") & (score_frame.group==g)]) for g in groups]
    a.bar(np.arange(4), sizes, color=colors, width=.6)
    for i, n in enumerate(sizes):
        a.text(i, n+.1, str(n), ha="center", fontsize=13)
    a.set_xticks(range(4), groups, rotation=15)
    a.set_ylabel("Matched cerebellum–hippocampus animals")
    a.set_ylim(0, max(sizes)+2)
    a.set_title("A  Pairing preserves the animal as the unit", loc="left", fontweight="bold")
    rr, rr10 = accounting
    union = rr["n_pairs"] + len(rr["unmatched_cb"]) + len(rr["unmatched_hip"])
    a.text(.02,.93,f"RRRM-2: {rr['n_pairs']} of {union} distinct source animals\nRR-10: {rr10['n_pairs']} matched animals", transform=a.transAxes,
           va="top", fontsize=11, bbox={"facecolor":"white", "edgecolor":"none", "alpha":.85})
    labels = {"RRRM2_region_by_flight_ISST":"Flight in ISS-T", "RRRM2_region_by_flight_LAR":"Flight in LAR",
              "RRRM2_region_by_cohort_GC":"Collection cohort in GC", "RRRM2_region_by_cohort_FLT":"Collection cohort in FLT",
              "RRRM2_region_by_flight_cohort_adjusted":"Flight, cohort-adjusted", "RR10_region_by_flight":"Flight in RR-10",
              "RRRM2_region_flight_cohort_triple":"Flight × collection cohort (approx.)"}
    main = results[results.sensitivity=="all_matched_animals"]
    matrix = main.pivot(index="contrast", columns="panel", values="paired_interaction_log2").loc[list(labels), list(PANELS)]
    q = main.pivot(index="contrast", columns="panel", values="q_BH_all_panel_tests").loc[matrix.index,matrix.columns]
    v = max(.1, abs(matrix.to_numpy()).max())
    im = b.imshow(matrix.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-v,vcenter=0,vmax=v), aspect="auto")
    b.set_xticks(range(len(matrix.columns)), [p.replace("Glia_Complement","Glia/complement").replace("BBB_Endothelium","BBB/endothelium") for p in matrix.columns], rotation=40, ha="right")
    b.set_yticks(range(len(matrix)), [labels[x] for x in matrix.index])
    for i in range(len(matrix)):
        for j in range(len(matrix.columns)):
            value=matrix.iloc[i,j]
            star="*" if q.iloc[i,j] <= .05 else ""
            b.text(j,i,f"{value:.2f}{star}",ha="center",va="center",color="white" if abs(value)>v*.55 else "black",fontsize=11)
    b.set_title(f"B  Hip-minus-Cb interactions; * BH q ≤ 0.05\nComplete family: {len(results)} panel tests",loc="left",fontweight="bold")
    fig.colorbar(im,ax=b,shrink=.7,label="Interaction effect on normalized log2 scale")
    rng=np.random.default_rng(SEED)
    sr=score_frame[(score_frame.mission=="RRRM-2") & (score_frame.panel=="SRG")]
    for i,g in enumerate(groups):
        vals=sr[sr.group==g]
        c.scatter(i+rng.uniform(-.13,.13,len(vals)), vals.paired_score_hip_minus_cb,
                  color=colors[i],s=65,edgecolors=np.where(vals.location_conflict,"black",colors[i]),linewidths=1.4,zorder=3)
        c.plot([i-.23,i+.23],[vals.paired_score_hip_minus_cb.mean()]*2,color="black",lw=2.2)
    c.set_xticks(range(4),groups,rotation=15)
    c.set_ylabel("Animal SRG score: hippocampus minus cerebellum")
    c.set_title("C  SRG variation across paired animals",loc="left",fontweight="bold")
    c.set_ylim(top=sr.paired_score_hip_minus_cb.max()+.65)
    c.text(.02,.98,"Black outlines: location metadata disagreement\nHorizontal lines: descriptive group means",transform=c.transAxes,va="top",fontsize=11)
    contrasts=["RRRM2_region_by_flight_ISST","RRRM2_region_by_flight_LAR","RR10_region_by_flight"]
    j=0
    for name in contrasts:
        for sensitivity, offset, color in [("all_matched_animals",0,"#4c78a8"),("exclude_metadata_conflicts",.23,"#b279a2")]:
            row=results[(results.contrast==name)&(results.panel=="SRG")&(results.sensitivity==sensitivity)]
            if row.empty: continue
            row=row.iloc[0]
            x=row.paired_interaction_log2
            d.hlines(j+offset,row.bootstrap_lower95,row.bootstrap_upper95,color=color,lw=1.8)
            d.scatter(x,j+offset,color=color,s=49,
                      label="All matched" if j==0 and offset==0 else ("Metadata sensitivity" if j==0 and offset else None))
            d.text(max(.65,results[results.panel=="SRG"].bootstrap_upper95.max()+.1), j+offset, f"q={row.q_BH_all_panel_tests:.3f}",fontsize=11,va="center")
        j+=1
    d.axvline(0,color="#777777",ls="--",lw=1)
    d.set_yticks(range(len(contrasts)),["RRRM-2 ISS-T","RRRM-2 LAR","RR-10"])
    d.invert_yaxis()
    d.set_xlabel("SRG regional flight interaction (Hip minus Cb), 95% bootstrap CI")
    d.set_title("D  SRG flight contrasts and metadata sensitivity",loc="left",fontweight="bold")
    d.legend(frameon=False,loc="upper left",bbox_to_anchor=(0,-.19),ncol=2)
    # Reserve room for q values without clipping labels.
    limits=d.get_xlim();d.set_xlim(limits[0],max(limits[1],results[results.panel=="SRG"].bootstrap_upper95.max()+.6))
    for ax in [a,c,d]:
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y" if ax!=d else "x",alpha=.15)
    fig.savefig(OUT/"figure_paired_regional_interactions.png",dpi=600,bbox_inches="tight")
    fig.savefig(OUT/"figure_paired_regional_interactions.svg",bbox_inches="tight")
    plt.close(fig)


def write_manuscript_summary(result, summary):
    main = result[result.sensitivity == "all_matched_animals"]
    def srg(name):
        return main[(main.contrast == name) & (main.panel == "SRG")].iloc[0]
    iss = srg("RRRM2_region_by_flight_ISST")
    lar = srg("RRRM2_region_by_flight_LAR")
    rr10 = srg("RR10_region_by_flight")
    cf = srg("RRRM2_region_by_cohort_FLT")
    cg = srg("RRRM2_region_by_cohort_GC")
    triple = srg("RRRM2_region_flight_cohort_triple")
    text = f"""# Paired regional biological follow-up

Source Name matching identified 24 cerebellum–hippocampus pairs in RRRM-2 (OSD-561/562; 24 of 29 distinct animals across the two deposits) and 10 pairs in RR-10 (OSD-563/564). For each animal, the mean panel expression in cerebellum was subtracted from that in hippocampus after per-deposit median-ratio normalization and log2 transformation. Eight inherited exploratory panels were tested for regional flight and ISS-T/LAR collection-cohort contrasts using animal-level label permutations restricted by age; the pooled RRRM-2 flight test additionally restricted exchange within collection cohort. The region × flight × collection-cohort coefficient used an additive-age model and approximate HC3-studentized restricted-model wild bootstrap. All 104 panel tests, including metadata-exclusion sensitivities and the approximate interactions, entered one BH family.

No panel test survived BH correction (minimum q = {result.q_BH_all_panel_tests.min():.3f}). SRG regional flight effects were {iss.paired_interaction_log2:+.3f} in ISS-T, {lar.paired_interaction_log2:+.3f} in LAR, and {rr10.paired_interaction_log2:+.3f} in RR-10 on the normalized log2 scale (q = {iss.q_BH_all_panel_tests:.3f}, {lar.q_BH_all_panel_tests:.3f}, and {rr10.q_BH_all_panel_tests:.3f}, respectively). Descriptive SRG ISS-T-minus-LAR regional cohort effects were larger in both FLT ({cf.paired_interaction_log2:+.3f}; permutation p = {cf.p_inference:.4f}, q = {cf.q_BH_all_panel_tests:.3f}) and GC ({cg.paired_interaction_log2:+.3f}; p = {cg.p_inference:.3f}, q = {cg.q_BH_all_panel_tests:.3f}). The approximate SRG three-way coefficient was {triple.paired_interaction_log2:+.3f} (wild-bootstrap p = {triple.p_inference:.3f}, q = {triple.q_BH_all_panel_tests:.3f}). These observations motivate further study of regional collection-cohort differences but do not support a multiplicity-adjusted biological discovery or establish an absence of flight effects.

Excluding six animals with location metadata disagreement in either region did not produce a significant panel result. This exclusion removed all four older FLT ISS-T animals; consequently, the ISS-T flight sensitivity is supported only by younger animals and has just 20 age-restricted assignments. Its approximate gene-level HC3 uncertainties are undefined because the remaining older GC animal has unit leverage, and those uncertainties are explicitly exported as unavailable. Percentile intervals use 2,000 animal bootstraps within exposure × cohort × age strata; the three-way test uses 9,999 wild draws. Effects describe relative bulk expression on separately normalized deposit scales and cannot identify causal flight responses, absolute tissue abundance, or cell-type-specific regulation.

Figure legend: A, matched RRRM-2 group sizes and pair accounting. B, effects for all eight panels in the primary contrasts; asterisks denote q ≤ 0.05 across the complete 104-test family (none observed). C, individual animal SRG hippocampus-minus-cerebellum scores; horizontal lines show descriptive group means and black outlines identify metadata disagreement. D, SRG regional flight effects with percentile 95% bootstrap intervals, including the metadata-exclusion sensitivity where applicable; annotations report the complete-family BH q values. The three-way heatmap row is exploratory approximate inference.
"""
    (OUT / "paired_manuscript_summary.md").write_text(text)


def run(args):
    if args.bootstraps < 1 or args.wild_draws < 1:
        raise ValueError("Bootstrap and wild-bootstrap draw counts must be positive")
    check_input_integrity()
    annotation=pd.read_csv(BASE/"inputs/gene_annotation.csv.gz")
    rr=paired_data("561","562")
    r10=paired_data("563","564")
    pd.concat([rr["meta"],r10["meta"]],ignore_index=True).to_csv(OUT/"paired_animal_manifest.csv",index=False)
    rows=[];genes=[];influence=[];gene_influence=[];scores=[]
    for dataset in [rr,r10]:
        pidx=indices(dataset["genes"],annotation)
        if any(len(idx) < 2 for idx in pidx.values()):
            raise ValueError("Insufficient retained genes for a panel and its leave-one-gene-out check")
        for panel,idx in pidx.items():
            for j,meta in dataset["meta"].iterrows():
                scores.append({"animal":meta.animal,"mission":meta.mission,"flight":meta.flight,"cohort":meta.cohort,
                               "group":f"{meta.flight} {meta.cohort.replace('ISS_T','ISS-T')}","age":meta.age,
                               "panel":panel,"n_genes":len(idx),"cb_score":dataset["cb"][idx,j].mean(),
                               "hip_score":dataset["hip"][idx,j].mean(),"paired_score_hip_minus_cb":dataset["delta"][idx,j].mean(),
                               "location_conflict":bool(meta.location_conflict)})
    m=rr["meta"]
    specs=[("RRRM2_region_by_flight_ISST",(m.cohort=="ISS_T").to_numpy(),"flight","FLT"),
           ("RRRM2_region_by_flight_LAR",(m.cohort=="LAR").to_numpy(),"flight","FLT"),
           ("RRRM2_region_by_cohort_GC",(m.flight=="GC").to_numpy(),"cohort","ISS_T"),
           ("RRRM2_region_by_cohort_FLT",(m.flight=="FLT").to_numpy(),"cohort","ISS_T"),
           ("RRRM2_region_by_flight_cohort_adjusted",np.ones(len(m),bool),"flight","FLT")]
    for sensitivity in ["all_matched_animals","exclude_metadata_conflicts"]:
        keep=np.ones(len(m),bool) if sensitivity=="all_matched_animals" else ~m.location_conflict.to_numpy()
        for name,mask,field,positive in specs:
            label_tests(name,rr,mask&keep,field,positive,annotation,sensitivity,rows,genes,influence,gene_influence,args.bootstraps)
        if sensitivity=="all_matched_animals":
            label_tests("RR10_region_by_flight",r10,np.ones(len(r10["meta"]),bool),"flight","FLT",annotation,sensitivity,rows,genes,influence,gene_influence,args.bootstraps)
        # Triple interaction uses separate null model but enters the same panel family.
        meta=m.loc[keep].reset_index(drop=True)
        idx=indices(rr["genes"],annotation)
        y=rr["delta"][:,keep]
        yp=np.vstack([y[idx[p]].mean(axis=0) for p in idx])
        effect,se,hc3p,p,boots,x=triple_interaction(yp,meta,args.bootstraps,args.wild_draws,np.random.default_rng(SEED))
        for i,panel in enumerate(idx):
            rows.append({"contrast":"RRRM2_region_flight_cohort_triple","sensitivity":sensitivity,"panel":panel,
                         "n_genes":len(idx[panel]),"n_paired_animals":len(meta),"n_positive":int((meta.flight=="FLT").sum()),
                         "n_negative":int((meta.flight=="GC").sum()),"paired_interaction_log2":effect[i],
                         "n_age_12":int((meta.age=="12").sum()),"n_age_29":int((meta.age=="29").sum()),
                         "bootstrap_lower95":np.quantile(boots[:,i],.025),"bootstrap_upper95":np.quantile(boots[:,i],.975),
                         "p_inference":p[i],"p_HC3_approximate":hc3p[i],"HC3_se":se[i],
                         "inference_method":"Approximate HC3-studentized restricted-model wild bootstrap",
                         "n_assignments":args.wild_draws,"total_assignments":"not a randomization enumeration",
                         "permutation_strata":"fixed design / animal wild weights","n_bootstraps":len(boots),
                         "interpretation_status":"exploratory; approximate additive-age interaction model",
                         "definition":"Regional flight interaction in ISS-T minus regional flight interaction in LAR"})
        gb,gs,_,gp,_,_,_=ols_hc3(y,x,x.shape[1]-1)
        genes.append(pd.DataFrame({"contrast":"RRRM2_region_flight_cohort_triple","sensitivity":sensitivity,
                                   "ensembl_id":rr["genes"],"region_interaction_log2":gb,"HC3_se":gs,"p_HC3_approximate":gp,
                                   "HC3_status":"estimable; approximate inference"}))
    result=pd.DataFrame(rows)
    result["q_BH_all_panel_tests"]=bh(result.p_inference)
    result["BH_family_size"]=len(result)
    assert result.groupby(["contrast", "sensitivity"]).panel.nunique().eq(len(PANELS)).all()
    assert np.isfinite(result.p_inference).all()
    result.to_csv(OUT/"paired_panel_interactions.csv",index=False)
    gframe=pd.concat(genes,ignore_index=True)
    gframe=gframe.merge(annotation[["ensembl_id","gene_symbol"]],on="ensembl_id",how="left",validate="many_to_one")
    gframe["q_BH_all_gene_interaction_tests"]=bh(gframe.p_HC3_approximate)
    gframe["gene_BH_family_size"]=int(gframe.p_HC3_approximate.notna().sum())
    gframe["n_planned_gene_effects"]=len(gframe)
    gframe.to_csv(OUT/"paired_gene_interactions_exploratory.csv.gz",index=False)
    pd.DataFrame(influence).to_csv(OUT/"leave_one_animal_out_effects.csv",index=False)
    pd.DataFrame(gene_influence).to_csv(OUT/"leave_one_gene_out_panel_effects.csv",index=False)
    score_frame=pd.DataFrame(scores)
    score_frame.to_csv(OUT/"paired_animal_panel_scores.csv",index=False)
    significant=result[result.q_BH_all_panel_tests<=.05]
    summary={"seed":SEED,"n_bootstraps_requested":args.bootstraps,"n_wild_bootstrap_draws":args.wild_draws,
             "dataset_pairing":[rr["accounting"],r10["accounting"]],"n_panel_tests":len(result),
             "n_significant_panel_tests_BH05":len(significant),"n_gene_interaction_tests":len(gframe),
             "n_gene_tests_with_undefined_HC3":int(gframe.p_HC3_approximate.isna().sum()),
             "n_evaluable_gene_HC3_tests":int(gframe.p_HC3_approximate.notna().sum()),
             "n_exploratory_gene_tests_BH05":int((gframe.q_BH_all_gene_interaction_tests<=.05).sum()),
             "normalization":"Per-deposit positive-gene median-ratio; log2(normalized counts + 1); common filtered gene universe per pair",
             "panel_status":"Inherited exploratory gene panels, not externally validated ground truth",
             "inference":"Animal-restricted label exchangeability tests; triple interaction uses approximate fixed-design wild bootstrap",
             "metadata_sensitivity":"Exclude an animal if either region disagrees with Source Name-derived location; 6 RRRM-2 animals excluded, including all 4 older FLT ISS-T animals. The sensitivity ISS-T flight coefficient is supported only by the younger age stratum; its gene-level HC3 SEs are undefined due to the single older GC animal.",
             "scale_limitation":"Effects compare relative bulk expression on separately normalized deposit scales; they do not estimate absolute tissue abundance, cell-type-specific regulation, or causal flight effects.",
             "significant_panels":significant.to_dict(orient="records"),"python":platform.python_version()}
    (OUT/"paired_analysis_summary.json").write_text(json.dumps(summary,indent=2,default=str)+"\n")
    write_manuscript_summary(result, summary)
    make_figure(result,score_frame,summary["dataset_pairing"])
    print(result[result.panel.isin(["SRG","Circadian","Glia_Complement"])][["contrast","sensitivity","panel","paired_interaction_log2","p_inference","q_BH_all_panel_tests"]].to_string(index=False))
    print(f"Wrote {len(result)} panel tests and {len(gframe)} exploratory gene tests to {OUT}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--bootstraps",type=int,default=2000)
    parser.add_argument("--wild-draws",type=int,default=9999)
    run(parser.parse_args())
