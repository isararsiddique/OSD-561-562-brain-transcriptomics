"""Paper 1 / Worker D (run inline by the parent): collection protocol versus spaceflight.

Question: in ISS rodent brain bulk RNA-seq, is the immediate-early-gene (IEG) panel signal
explained by the collection protocol rather than by spaceflight?

Design facts (deposited metadata, RRRM-2 OSD-561 cerebellum / OSD-562 hippocampus):
  FLT ISS-T: ketamine/xylazine overdose, cryochiller      GC ISS-T: ketamine/xylazine, liquid N2
  FLT LAR  : gas anaesthesia, LN2 (bulk)                  GC LAR  : gas anaesthesia, LN2 (bulk)
So: protocol effect in ground   P_G = GC ISS-T - GC LAR  (no flight; euthanasia + schedule differ)
    matched flight effect       F_L = FLT LAR  - GC LAR  (collection fully matched)
    flight, freezing differs    F_I = FLT ISS-T - GC ISS-T (euthanasia matched)
    protocol effect in flight   P_F = FLT ISS-T - FLT LAR
Dominance test (new): GC ISS-T versus FLT LAR. These two groups differ in BOTH factors in
opposite directions relative to the shared GC LAR reference, so a positive contrast means the
protocol effect exceeds the matched flight effect. Exact whole-animal permutation within age.

The package's own validated loaders/estimators are reused (read-only import) so that every
reproduced number is identical to results/panel_results.csv. Nothing in the package is written.
"""
import itertools
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import pandas as pd

PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(PKG / "scripts"))
import run_revision as rr  # noqa: E402  (module-level code only defines constants/functions)
from statistics_core import Contrast, animal_permutations, permutation_p, bh  # noqa: E402

SEED = 15102026
N_BOOT = 2000
PANELS = list(rr.PANELS)
(OUT / "D_tables").mkdir(exist_ok=True)


def panel_index(genes, annotation):
    return rr.indices(genes, annotation)


def subset_effect(data, mask, field, positive, idx):
    """Panel-mean OLS effect (age/tissue-adjusted) exactly as rr.analyse_contrast computes it."""
    meta = data["meta"].loc[mask]
    y = data["y"].loc[:, meta.index].to_numpy()
    z = rr.nuisance(meta)
    g = (meta[field] == positive).to_numpy(int)
    beta, _ = Contrast(y[idx], z).fit(g)
    return float(beta.mean())


def exact_test(data, mask, field, positive, idx_by_panel):
    meta = data["meta"].loc[mask].copy()
    y = data["y"].loc[:, meta.index].to_numpy()
    z = rr.nuisance(meta)
    g = (meta[field] == positive).to_numpy(int)
    blocks = ("age",) if meta.age.nunique() > 1 else ()
    labels, method, total = animal_permutations(meta, g, blocks)
    out = {}
    for p, idx in idx_by_panel.items():
        c = Contrast(y[idx], z)
        obs = float(c.fit(g)[0].mean())
        null = c.permutation_betas(labels).mean(axis=0)
        out[p] = {"effect": obs, "p_exact": float(permutation_p(obs, null, method)),
                  "n_assignments": int(labels.shape[1]), "method": method,
                  "n_animals": int(meta.animal.nunique()),
                  "n_pos_animals": int(meta.loc[g == 1].animal.nunique())}
    return out


def main():
    annotation = pd.read_csv(PKG / "inputs/gene_annotation.csv.gz")
    data = {n: rr.load_dataset(n) for n in rr.DATASETS}
    pool = rr.pooled(data["561"], data["562"])
    regions = {"cerebellum": data["561"], "hippocampus": data["562"], "pooled": pool}
    ref = pd.read_csv(PKG / "results/panel_results.csv")
    ref = ref[ref.sensitivity == "design_cohort"].set_index(["contrast", "panel"])

    # ------------------------------------------------------------------ 1. reproduction
    repro, rows = [], []
    contrast_specs = {
        "P_G": ("GC", "cohort", "ISS_T", lambda m: m.flight == "GC", "{lab}_GC_ISST_vs_LAR"),
        "P_F": ("FLT", "cohort", "ISS_T", lambda m: m.flight == "FLT", "{lab}_FLT_ISST_vs_LAR"),
        "F_I": (None, "flight", "FLT", lambda m: m.cohort == "ISS_T", "{lab}_FLT_vs_GC_within_ISS_T"),
        "F_L": (None, "flight", "FLT", lambda m: m.cohort == "LAR", "{lab}_FLT_vs_GC_within_LAR"),
    }
    labels = {"cerebellum": "OSD561_cerebellum", "hippocampus": "OSD562_hippocampus", "pooled": "RRRM2_pooled"}
    exact = {}
    for region, d in regions.items():
        idx = panel_index(d["y"].index, annotation)
        for key, (_, field, pos, sel, name) in contrast_specs.items():
            if region == "pooled" and key in ("F_I", "F_L"):
                continue  # the package fits within-cohort flight effects per region only
            res = exact_test(d, sel(d["meta"]), field, pos, idx)
            exact[(region, key)] = res
            cname = name.format(lab=labels[region])
            for p in PANELS:
                pkg = ref.loc[(cname, p)]
                diff = res[p]["effect"] - float(pkg.mean_log2_effect)
                repro.append({"region": region, "contrast": key, "package_contrast": cname, "panel": p,
                              "recomputed": res[p]["effect"], "package": float(pkg.mean_log2_effect),
                              "abs_diff": abs(diff), "p_recomputed": res[p]["p_exact"],
                              "p_package": float(pkg.p_animal_permutation)})
    repro = pd.DataFrame(repro)
    repro.to_csv(OUT / "D_tables/D1_reproduction_vs_package.csv", index=False)
    max_diff = float(repro.abs_diff.max())
    max_pdiff = float((repro.p_recomputed - repro.p_package).abs().max())
    print(f"reproduction: max |effect diff| = {max_diff:.2e}; max |p diff| = {max_pdiff:.2e}", flush=True)
    if max_diff > 1e-9 or max_pdiff > 1e-12:
        raise SystemExit("Reproduction failed; stop before extending.")

    # ------------------------------------------------------------------ 2. dominance test
    dom_rows = []
    for region, d in regions.items():
        m = d["meta"]
        mask = ((m.flight == "GC") & (m.cohort == "ISS_T")) | ((m.flight == "FLT") & (m.cohort == "LAR"))
        dd = {"y": d["y"], "meta": m.assign(dom=np.where(m.flight == "GC", "GC_ISST", "FLT_LAR"))}
        idx = panel_index(d["y"].index, annotation)
        res = exact_test(dd, mask, "dom", "GC_ISST", idx)
        for p in PANELS:
            r = res[p]
            dom_rows.append({"region": region, "panel": p, "GC_ISST_minus_FLT_LAR": r["effect"],
                             "p_exact": r["p_exact"], "n_assignments": r["n_assignments"],
                             "n_animals": r["n_animals"], "method": r["method"]})
    dom = pd.DataFrame(dom_rows)
    dom["q_BH_family24"] = bh(dom.p_exact)  # complete family: 8 panels x 3 analyses (not preregistered)
    dom.to_csv(OUT / "D_tables/D2_dominance_GC_ISST_vs_FLT_LAR.csv", index=False)

    # ------------------------------------------------------------------ 3. bootstrap of the 2x2
    rng = np.random.default_rng(SEED)
    boot_rows, summ_rows = [], []
    for region in ("cerebellum", "hippocampus", "pooled"):
        d = regions[region]
        m = d["meta"]
        idx = panel_index(d["y"].index, annotation)
        avg = pd.DataFrame({p: d["y"].to_numpy()[ix].mean(axis=0) for p, ix in idx.items()}, index=m.index)
        per = m.drop_duplicates("animal")
        cells = [g.animal.tolist() for _, g in per.groupby(["flight", "cohort", "age"], sort=True)]

        def effects(meta):
            out = {}
            for key, (_, field, pos, sel, _) in contrast_specs.items():
                sub = meta[sel(meta)]
                g = (sub[field] == pos).to_numpy(int)
                out[key] = Contrast(avg.loc[sub.index].T.to_numpy(), rr.nuisance(sub)).fit(g)[0]
            return out

        obs = effects(m)
        reps = {k: [] for k in obs}
        for _ in range(N_BOOT):
            chosen = list(itertools.chain.from_iterable(rng.choice(c, len(c), replace=True) for c in cells))
            parts = []
            for j, a in enumerate(chosen):
                rows_a = m[m.animal == a].copy()
                rows_a["animal"] = f"{a}#b{j}"
                rows_a.index = [f"{i}#b{j}" for i in rows_a.index]
                parts.append(rows_a)
            bm = pd.concat(parts)
            avg_b = avg.loc[[i.split("#b")[0] for i in bm.index]].copy()
            avg_b.index = bm.index
            try:
                e = {}
                for key, (_, field, pos, sel, _) in contrast_specs.items():
                    sub = bm[sel(bm)]
                    z = rr.nuisance(sub)
                    g = (sub[field] == pos).to_numpy(int)
                    e[key] = Contrast(avg_b.loc[sub.index].T.to_numpy(), z).fit(g)[0]
            except ValueError:
                continue
            for k in reps:
                reps[k].append(e[k])
        reps = {k: np.array(v) for k, v in reps.items()}
        n_ok = len(reps["P_G"])
        for j, p in enumerate(avg.columns):
            pg, fl, fi, pf = (reps[k][:, j] for k in ("P_G", "F_L", "F_I", "P_F"))
            delta = pg - fl
            row = {"region": region, "panel": p, "n_boot_ok": n_ok,
                   **{f"{k}": float(obs[k][j]) for k in obs},
                   **{f"{k}_lo95": float(np.quantile(reps[k][:, j], .025)) for k in obs},
                   **{f"{k}_hi95": float(np.quantile(reps[k][:, j], .975)) for k in obs},
                   "delta_PG_minus_FL": float(obs["P_G"][j] - obs["F_L"][j]),
                   "delta_lo95": float(np.quantile(delta, .025)),
                   "delta_hi95": float(np.quantile(delta, .975)),
                   "boot_frac_PG_gt_FL": float(np.mean(pg > fl)),
                   "boot_frac_absPG_gt_absFL": float(np.mean(np.abs(pg) > np.abs(fl))),
                   "ratio_PG_over_FL": float(obs["P_G"][j] / obs["F_L"][j]) if obs["F_L"][j] != 0 else np.nan}
            summ_rows.append(row)
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(OUT / "D_tables/D3_two_by_two_bootstrap.csv", index=False)

    # ------------------------------------------------------------------ 4. gene-level view
    gene_rows = []
    for region in ("cerebellum", "hippocampus"):
        d = regions[region]
        m = d["meta"]
        idx = panel_index(d["y"].index, annotation)
        union = np.unique(np.concatenate([idx[p] for p in ("rPRG", "dPRG", "HSP", "SRG")]))
        sym = annotation.set_index("ensembl_id").gene_symbol.reindex(d["y"].index).to_numpy()
        fits = {}
        for key in ("P_G", "F_L", "F_I", "P_F"):
            _, field, pos, sel, _ = contrast_specs[key]
            sub = m[sel(m)]
            c = Contrast(d["y"].loc[:, sub.index].to_numpy()[union], rr.nuisance(sub))
            b, se, t, pv, df = c.cluster_fit((sub[field] == pos).to_numpy(int), sub.animal.to_numpy())
            fits[key] = (b, se)
        for j, gi in enumerate(union):
            panels = [p for p in ("rPRG", "dPRG", "HSP", "SRG") if gi in set(idx[p])]
            gene_rows.append({"region": region, "ensembl_id": d["y"].index[gi], "symbol": sym[gi],
                              "panels": "|".join(panels),
                              **{f"{k}": float(fits[k][0][j]) for k in fits},
                              **{f"{k}_se": float(fits[k][1][j]) for k in fits}})
    genes = pd.DataFrame(gene_rows)
    genes.to_csv(OUT / "D_tables/D4_gene_level_effects.csv", index=False)

    # ------------------------------------------------------------------ 5. cross-mission classes
    allref = pd.read_csv(PKG / "results/panel_results.csv")
    allref = allref[allref.sensitivity == "design_cohort"]
    classes = [
        ("RRRM2_pooled_GC_ISST_vs_LAR", "RRRM-2 pooled", "Protocol differs, no flight"),
        ("OSD561_cerebellum_GC_ISST_vs_LAR", "RRRM-2 cerebellum", "Protocol differs, no flight"),
        ("OSD562_hippocampus_GC_ISST_vs_LAR", "RRRM-2 hippocampus", "Protocol differs, no flight"),
        ("RRRM2_pooled_FLT_ISST_vs_LAR", "RRRM-2 pooled", "Protocol differs, all flown"),
        ("OSD561_cerebellum_FLT_ISST_vs_LAR", "RRRM-2 cerebellum", "Protocol differs, all flown"),
        ("OSD562_hippocampus_FLT_ISST_vs_LAR", "RRRM-2 hippocampus", "Protocol differs, all flown"),
        ("OSD561_cerebellum_FLT_vs_GC_within_ISS_T", "RRRM-2 cerebellum", "Flight, euthanasia matched, freezing differs"),
        ("OSD562_hippocampus_FLT_vs_GC_within_ISS_T", "RRRM-2 hippocampus", "Flight, euthanasia matched, freezing differs"),
        ("OSD613_FLT_vs_GC", "RRRM-2 left hemisphere", "Flight, euthanasia matched, freezing differs"),
        ("OSD561_cerebellum_FLT_vs_GC_within_LAR", "RRRM-2 cerebellum", "Flight, protocol matched"),
        ("OSD562_hippocampus_FLT_vs_GC_within_LAR", "RRRM-2 hippocampus", "Flight, protocol matched"),
        ("OSD563_FLT_vs_GC", "RR-10 cerebellum", "Flight, protocol matched"),
        ("OSD564_FLT_vs_GC", "RR-10 hippocampus", "Flight, protocol matched"),
        ("OSD612_FLT_vs_GC", "RR-10 left hemisphere", "Flight, protocol matched"),
        ("OSD525_FLT_vs_HGC", "RR-23 cerebellum (vs HGC)", "Flight, protocol matched"),
        ("OSD525_FLT_vs_VGC", "RR-23 cerebellum (vs VGC)", "Flight, protocol matched"),
        ("OSD525_HGC_vs_VGC", "RR-23 cerebellum", "Housing only, protocol matched"),
    ]
    cm = []
    for cname, label, cls in classes:
        for p in PANELS:
            r = allref[(allref.contrast == cname) & (allref.panel == p)]
            if len(r) != 1:
                raise SystemExit(f"missing {cname} {p}")
            r = r.iloc[0]
            cm.append({"contrast": cname, "label": label, "class": cls, "panel": p,
                       "effect": r.mean_log2_effect, "lo95": r.bootstrap_lower95, "hi95": r.bootstrap_upper95,
                       "p_exact": r.p_animal_permutation, "q_BH_232": r.q_BH_all_tests_in_family,
                       "n_animals": int(r.n_animals), "n_pos": int(r.n_positive_animals)})
    cm = pd.DataFrame(cm)
    cm.to_csv(OUT / "D_tables/D5_cross_mission_by_protocol_class.csv", index=False)
    cls_summary = (cm.assign(abs_effect=cm.effect.abs(), sig=cm.q_BH_232 <= 0.05)
                   .groupby(["panel", "class"]).agg(n=("effect", "size"), mean_effect=("effect", "mean"),
                                                    mean_abs=("abs_effect", "mean"),
                                                    max_abs=("abs_effect", "max"), n_q05=("sig", "sum"))
                   .reset_index())
    cls_summary.to_csv(OUT / "D_tables/D5b_class_summary.csv", index=False)

    out = {
        "seed": SEED, "n_boot": N_BOOT,
        "reproduction": {"max_abs_effect_diff": max_diff, "max_abs_p_diff": max_pdiff,
                         "n_rows": int(len(repro))},
        "dominance": dom.to_dict(orient="records"),
        "two_by_two": summ.to_dict(orient="records"),
        "cross_mission_class_summary": cls_summary.to_dict(orient="records"),
        "design_note": ("RR-10 OSD-563/564/612: ketamine/xylazine and cryochiller carcass freezing for every "
                        "animal; RR-23 OSD-525: inhalant sedation, cardiac puncture, RNAlater for every animal; "
                        "OSD-613: cryochiller (FLT) vs liquid nitrogen (GC)."),
    }
    (OUT / "D_results.json").write_text(json.dumps(out, indent=1, default=float))
    print("wrote D_results.json and D_tables/", flush=True)


if __name__ == "__main__":
    main()
