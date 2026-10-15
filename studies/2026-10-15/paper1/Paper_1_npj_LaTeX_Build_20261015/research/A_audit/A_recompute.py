"""Worker A, step 2: independent recomputation of Paper 1 headline numbers from raw inputs.

Own code throughout. The package scripts were READ to learn definitions only; nothing is
imported from them. Inputs are the bundled NASA count matrices and deposited metadata in
the read-only source package. Each section appends rows of
  {id, claim, expected, recomputed, status, note}
to A_recompute_results.json (merged across runs, so sections can be run separately):

  python A_recompute.py a b c d e f g h i      # or any subset

Sections
  a  accounting: 129 libraries, 126 animal-region columns, 83 animals, 24 Cb+Hip pairs, 6 conflict animals
  b  RRRM-2 panel effects (rPRG/HSP/SRG + all 8 panels) with exact whole-animal permutation p in age blocks
  c  score-operator spectral gains ||gamma|| for all 15 anchor-target geometries
  d  ground holdout (81 splits each, OSD-561 and OSD-562): sensor AUC, contrast cosine, score MSE ratio
  e  cross-cohort intersection-union screen (54,989 tests) from component Wald p; plus an
     estimator-independent OLS re-screen for one region
  f  cerebellar Hbb-bt / Hbb-bs flight x collection interaction (own NB GLM + own OLS)
  g  hippocampal Circadian interaction: OLS, stratified bootstrap, HC3 restricted wild bootstrap
  h  external units 64 = 46 individual + 18 pooled (from count columns + NASA SampleTables)
  i  gain-constrained mean MSE-ratio improvements from per-holdout rows
  j  NEW own nested gain-2 constrained regression for the 4 external datasets from raw STAR counts
  k  NEW random-sensor null for the score-operator gain
  l  NEW margin of the constrained estimator over matched clipped ridge, per holdout
  m  NEW own nested gain-2 constrained regression for the 81+81 original ground holdouts

Status: MATCH (equal at the written precision), ROUNDING (differs only in the last written
digit), DIRECTION (same sign / qualitative conclusion, magnitude differs because the
resampling seed or a nuisance estimator differs), MISMATCH, NOT_CHECKED.
"""
import itertools
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
OUT = Path(__file__).resolve().parent
RESULTS = OUT / "A_recompute_results.json"
CORE = {"525": ("RR-23", "cerebellum"), "561": ("RRRM-2", "cerebellum"), "562": ("RRRM-2", "hippocampus"),
        "613": ("RRRM-2", "left hemisphere"), "563": ("RR-10", "cerebellum"), "564": ("RR-10", "hippocampus"),
        "612": ("RR-10", "left hemisphere")}
# Panel gene lists as printed in results/panel_definitions.csv (read at run time, not copied).
PANEL_DEF = pd.read_csv(PKG / "results/panel_definitions.csv")
PANELS = {p: g.symbol.tolist() for p, g in PANEL_DEF.groupby("panel", sort=False)}
ANNOT = pd.read_csv(PKG / "inputs/gene_annotation.csv.gz").set_index("ensembl_id").gene_symbol.to_dict()
ROWS = []


def record(rid, claim, expected, recomputed, status, note=""):
    def clean(x):
        if isinstance(x, (np.floating, np.integer)):
            return x.item()
        return x
    ROWS.append(dict(id=rid, claim=claim, expected=clean(expected), recomputed=clean(recomputed),
                     status=status, note=note))
    print(f"[{status:9s}] {rid}: expected {expected} | recomputed {recomputed}  {note}", flush=True)


def at_precision(written, value):
    """Return MATCH if value rounds to the written string, ROUNDING if off by one last digit."""
    s = str(written).replace("−", "-").replace(",", "")
    dec = len(s.split(".")[1]) if "." in s else 0
    w = float(s)
    r = round(float(value), dec)
    if abs(r - w) < 10 ** (-dec) * 0.01:
        return "MATCH"
    if abs(r - w) <= 10 ** (-dec) * 1.01:
        return "ROUNDING"
    return "MISMATCH"


# ----------------------------------------------------------------------------- data loading

def source_groups(n, sn):
    mission = CORE[n][0]
    if mission == "RRRM-2":
        return ("FLT" if "_FL_" in sn else "GC"), ("ISS_T" if "_ISS_" in sn else "LAR")
    if mission == "RR-23":
        return {"F": "FLT", "G": "HGC", "V": "VGC"}[sn.split("_")[-1][0]], "RR23"
    return ("FLT" if "_FL-" in sn else "GC"), "matched"


def load(n):
    counts = pd.read_csv(PKG / f"inputs/OSD-{n}_counts.csv.gz", index_col=0)
    md = pd.read_csv(PKG / f"inputs/OSD-{n}_metadata.csv.gz", dtype=str).set_index("sample").loc[counts.columns]
    animal = md["Source Name"].str.strip()
    collapsed = counts.T.groupby(animal.values).sum().T          # sum technical libraries per animal
    meta = pd.DataFrame(index=collapsed.columns)
    meta["animal"] = meta.index
    fl = [source_groups(n, a) for a in meta.index]
    meta["flight"] = [f for f, _ in fl]
    meta["cohort"] = [c for _, c in fl]
    agecol = "Factor Value[Age]" if "Factor Value[Age]" in md.columns else None
    first = md.assign(animal=animal.values).drop_duplicates("animal").set_index("animal")
    meta["age"] = first.loc[meta.index, agecol].values if agecol else "single"
    loc = "Factor Value[Euthanasia Location]"
    if CORE[n][0] == "RRRM-2" and loc in md.columns:
        exp = np.where((meta.flight == "FLT") & (meta.cohort == "ISS_T"), "On ISS", "On Earth")
        dep = first.loc[meta.index, loc].fillna("").values
        meta["deposited_location"] = dep
        meta["conflict"] = (dep != "") & (dep != exp)
    else:
        meta["deposited_location"] = ""
        meta["conflict"] = False
    meta["n_lib"] = animal.value_counts().reindex(meta.index).values
    return dict(raw=counts, counts=collapsed, meta=meta, md=md)


def normalise(counts, min_mean=10):
    canon = counts.index.astype(str).str.fullmatch(r"ENSMUSG\d+")
    c = counts.loc[canon]
    c = c.loc[c.mean(axis=1) >= min_mean]
    a = c.to_numpy(float)
    pos = (a > 0).all(axis=1)
    la = np.log(a[pos])
    sf = np.exp(np.median(la - la.mean(axis=1, keepdims=True), axis=0))
    return pd.DataFrame(np.log2(a / sf + 1), index=c.index, columns=c.columns), sf


def panel_index(genes):
    sym = np.array([ANNOT.get(g, "") for g in genes])
    out = {p: np.flatnonzero(np.isin(sym, s)) for p, s in PANELS.items()}
    return out


def resid_proj(z, n):
    x = np.ones((n, 1)) if z is None else np.column_stack([np.ones(n), z])
    q, r = np.linalg.qr(x)
    keep = np.abs(np.diag(r)) > 1e-10
    return q[:, keep]


def nuisance_cols(meta):
    cols = []
    if meta.age.nunique() > 1:
        cols.append((meta.age == "29").to_numpy(float))
    if "tissue" in meta and meta.tissue.nunique() > 1:
        cols.append((meta.tissue == meta.tissue.unique()[1]).to_numpy(float))
    return np.column_stack(cols) if cols else None


def fwl_beta(y, g, q):
    """Group coefficient(s) after projecting out nuisance basis q. y genes x n, g n or n x B."""
    y0 = y - (y @ q) @ q.T
    g = np.asarray(g, float)
    gt = g - q @ (q.T @ g)
    return (y0 @ gt) / np.sum(gt * gt, axis=0)


def exact_assignments(meta, label, blocks):
    per = meta.assign(label=label).drop_duplicates("animal").set_index("animal")
    strata = [(sorted(b.index), int(b.label.sum())) for _, b in
              (per.groupby(list(blocks), sort=True) if blocks else [(None, per)])]
    choices = [list(itertools.combinations(a, k)) for a, k in strata]
    sets = [set(itertools.chain.from_iterable(c)) for c in itertools.product(*choices)]
    animals = meta.animal.to_numpy()
    return np.column_stack([np.isin(animals, list(s)).astype(float) for s in sets])


# ----------------------------------------------------------------------------- section a

def section_a():
    total_lib, total_cols, animals_by_mission = 0, 0, {}
    per = {}
    for n in CORE:
        d = load(n)
        y, _ = normalise(d["counts"])
        per[n] = d
        total_lib += d["raw"].shape[1]
        total_cols += d["counts"].shape[1]
        animals_by_mission.setdefault(CORE[n][0], set()).update(d["meta"].animal)
        acc = pd.read_csv(PKG / "results/dataset_accounting.csv").set_index("accession").loc[f"OSD-{n}"]
        record(f"a.libs.OSD-{n}", f"OSD-{n} libraries", int(acc.n_deposited_libraries_in_count_matrix), d["raw"].shape[1],
               "MATCH" if acc.n_deposited_libraries_in_count_matrix == d["raw"].shape[1] else "MISMATCH")
        record(f"a.cols.OSD-{n}", f"OSD-{n} animal-region columns", int(acc.n_animals_in_region), d["counts"].shape[1],
               "MATCH" if acc.n_animals_in_region == d["counts"].shape[1] else "MISMATCH")
        record(f"a.genes.OSD-{n}", f"OSD-{n} genes retained (canonical ENSMUSG, mean>=10)", int(acc.n_filtered_genes), len(y),
               "MATCH" if acc.n_filtered_genes == len(y) else "MISMATCH")
    n_animals = sum(len(s) for s in animals_by_mission.values())
    record("a.libraries", "Seven count matrices contain 129 libraries", 129, total_lib, "MATCH" if total_lib == 129 else "MISMATCH")
    record("a.columns", "126 animal-region columns", 126, total_cols, "MATCH" if total_cols == 126 else "MISMATCH")
    record("a.animals", "83 distinct animals across RRRM-2, RR-10 and RR-23", 83, n_animals,
           "MATCH" if n_animals == 83 else "MISMATCH",
           "distinct Source Name per mission: " + ", ".join(f"{k} {len(v)}" for k, v in animals_by_mission.items()))
    cb, hip = set(per["561"]["meta"].animal), set(per["562"]["meta"].animal)
    record("a.pairs_rrrm2", "24 animals contributed both cerebellum and hippocampus", 24, len(cb & hip),
           "MATCH" if len(cb & hip) == 24 else "MISMATCH")
    record("a.rrrm2_union", "OSD-561/562 represent 29 animals", 29, len(cb | hip), "MATCH" if len(cb | hip) == 29 else "MISMATCH")
    rr10 = set(per["563"]["meta"].animal) & set(per["564"]["meta"].animal)
    record("a.pairs_rr10", "OSD-563/564 provide ten matched pairs", 10, len(rr10), "MATCH" if len(rr10) == 10 else "MISMATCH")
    in612 = set(per["612"]["meta"].animal) <= (set(per["563"]["meta"].animal) | set(per["564"]["meta"].animal))
    record("a.612_subset", "OSD-612 contains nine of that RR-10 cohort's animals", "9, subset", f"{len(per['612']['meta'])}, subset={in612}",
           "MATCH" if in612 and len(per["612"]["meta"]) == 9 else "MISMATCH")
    other = set(per["613"]["meta"].animal) - (cb | hip)
    record("a.613_other", "OSD-613 contained 18 other RRRM-2 animals", 18, len(other), "MATCH" if len(other) == 18 else "MISMATCH")
    m525 = per["525"]["meta"]
    vc = m525.flight.value_counts().to_dict()
    record("a.rr23", "RR-23: 29 libraries, 26 animals: 9 FLT, 9 HGC, 8 VGC; 3 technical replicates summed",
           "29/26/9/9/8/3", f"{per['525']['raw'].shape[1]}/{len(m525)}/{vc.get('FLT')}/{vc.get('HGC')}/{vc.get('VGC')}/{int((m525.n_lib > 1).sum())}",
           "MATCH" if (per['525']['raw'].shape[1], len(m525), vc.get('FLT'), vc.get('HGC'), vc.get('VGC'), int((m525.n_lib > 1).sum())) == (29, 26, 9, 9, 8, 3) else "MISMATCH")
    # conflicts
    entries = sum(int(per[n]["meta"].conflict.sum()) for n in ("561", "562"))
    animals = set(per["561"]["meta"].animal[per["561"]["meta"].conflict]) | set(per["562"]["meta"].animal[per["562"]["meta"].conflict])
    record("a.conflicts", "Seven location entries from six core animals disagreed with Source Name cohort", "7 / 6",
           f"{entries} / {len(animals)}", "MATCH" if (entries, len(animals)) == (7, 6) else "MISMATCH",
           "conflict animals: " + ", ".join(sorted(animals)))
    # Older FLT ISS-T animals among conflicts
    for n in ("561", "562"):
        m = per[n]["meta"]
        old = m[(m.flight == "FLT") & (m.cohort == "ISS_T") & (m.age == "29")]
        kept = old[~old.animal.isin(animals)]
        record(f"a.old_flt_isst_after_exclusion.OSD-{n}", f"OSD-{n} older FLT ISS-T animals left after excluding conflict animals",
               {"561": 0, "562": 1}[n], len(kept), "MATCH" if len(kept) == {"561": 0, "562": 1}[n] else "MISMATCH",
               f"{len(old)} older FLT ISS-T before exclusion")
    pairs_old_flt = [a for a in (cb & hip) if a in set(per["561"]["meta"].query("flight=='FLT' and cohort=='ISS_T' and age=='29'").animal)]
    record("a.paired_old_flt_isst", "Metadata exclusion removed all four older FLT ISS-T paired animals", 4,
           f"{len(pairs_old_flt)} paired, {sum(a in animals for a in pairs_old_flt)} of them conflict animals",
           "MATCH" if len(pairs_old_flt) == 4 and all(a in animals for a in pairs_old_flt) else "MISMATCH")
    # ground animals used by the 362 holdouts: GC in 561/562 plus HGC/VGC in 525
    g = set(per["561"]["meta"].query("flight=='GC'").animal) | set(per["562"]["meta"].query("flight=='GC'").animal)
    g525 = set(m525.query("flight in ['HGC','VGC']").animal)
    record("a.ground_unique", "362 splits use 30 unique ground animals", 30, len(g) + len(g525),
           "MATCH" if len(g) + len(g525) == 30 else "MISMATCH", f"RRRM-2 GC {len(g)} + RR-23 HGC/VGC {len(g525)}")


# ----------------------------------------------------------------------------- section b

def rrrm2(which):
    d = load(which)
    y, _ = normalise(d["counts"])
    m = d["meta"].copy()
    m["tissue"] = CORE[which][1]
    m.index = [f"{which}:{a}" for a in m.animal]
    y.columns = m.index
    return y, m


def panel_contrast(y, m, mask, field, positive, panels, nboot=0, seed=0):
    mm = m.loc[mask]
    yy = y.loc[:, mm.index].to_numpy()
    lab = (mm[field] == positive).to_numpy(float)
    blocks = ("age",) if mm.age.nunique() > 1 else ()
    q = resid_proj(nuisance_cols(mm), len(mm))
    idx = panel_index(y.index)
    pm = np.vstack([yy[idx[p]].mean(axis=0) for p in panels])        # panel means (linear => mean of gene betas)
    obs = fwl_beta(pm, lab, q)
    L = exact_assignments(mm, lab, blocks)
    null = fwl_beta(pm, L, q)
    p = (np.abs(null) >= np.abs(obs)[:, None] - 1e-12).mean(axis=1)
    out = {pn: dict(effect=obs[i], p=p[i], n_genes=len(idx[pn]), n_assign=L.shape[1], n_cols=len(mm),
                    n_animals=mm.animal.nunique()) for i, pn in enumerate(panels)}
    if nboot:
        rng = np.random.default_rng(seed)
        per = mm.assign(label=lab).drop_duplicates("animal")
        strata = [f.animal.tolist() for _, f in per.groupby(["label", "age"], sort=True)]
        an = mm.animal.to_numpy()
        boots = []
        for _ in range(nboot):
            ch = list(itertools.chain.from_iterable(rng.choice(s, len(s), replace=True) for s in strata))
            ix = np.concatenate([np.flatnonzero(an == a) for a in ch])
            z = nuisance_cols(mm.iloc[ix])
            boots.append(fwl_beta(pm[:, ix], lab[ix], resid_proj(z, len(ix))))
        boots = np.array(boots)
        for i, pn in enumerate(panels):
            out[pn]["boot_lo"], out[pn]["boot_hi"] = np.quantile(boots[:, i], [.025, .975])
    return out


def section_b():
    y1, m1 = rrrm2("561")
    y2, m2 = rrrm2("562")
    common = y1.index.intersection(y2.index)
    yp = pd.concat([y1.loc[common], y2.loc[common]], axis=1)
    mp = pd.concat([m1, m2])
    pkg = pd.read_csv(PKG / "results/panel_results.csv")
    pkg = pkg[pkg.sensitivity == "design_cohort"].set_index(["contrast", "panel"])
    panels = list(PANELS)
    sets = {"RRRM2_pooled": (yp, mp), "OSD561_cerebellum": (y1, m1), "OSD562_hippocampus": (y2, m2)}
    table = []
    for label, (y, m) in sets.items():
        specs = [(f"{label}_GC_ISST_vs_LAR", m.flight == "GC", "cohort", "ISS_T"),
                 (f"{label}_FLT_ISST_vs_LAR", m.flight == "FLT", "cohort", "ISS_T"),
                 (f"{label}_FLT_vs_GC_within_ISS_T", m.cohort == "ISS_T", "flight", "FLT"),
                 (f"{label}_FLT_vs_GC_within_LAR", m.cohort == "LAR", "flight", "FLT")]
        for name, mask, field, pos in specs:
            boot = 1000 if name == "RRRM2_pooled_GC_ISST_vs_LAR" else 0
            res = panel_contrast(y, m, mask.to_numpy(), field, pos, panels, nboot=boot, seed=20261015)
            for pn, r in res.items():
                key = (name, pn)
                row = dict(contrast=name, panel=pn, **r)
                if key in pkg.index:
                    pr = pkg.loc[key]
                    row.update(pkg_effect=pr.mean_log2_effect, pkg_p=pr.p_animal_permutation, pkg_n_assign=pr.n_assignments_evaluated,
                               pkg_n_genes=pr.n_genes, pkg_boot_lo=pr.bootstrap_lower95, pkg_boot_hi=pr.bootstrap_upper95)
                table.append(row)
    # remaining design contrasts: OSD-613, RR-10 deposits, RR-23 pairs
    for n in ("563", "564", "612", "613"):
        y, m = rrrm2(n)
        name = f"OSD{n}_FLT_vs_GC"
        res = panel_contrast(y, m, np.ones(len(m), bool), "flight", "FLT", panels)
        for pn, r in res.items():
            row = dict(contrast=name, panel=pn, **r)
            if (name, pn) in pkg.index:
                pr = pkg.loc[(name, pn)]
                row.update(pkg_effect=pr.mean_log2_effect, pkg_p=pr.p_animal_permutation, pkg_n_assign=pr.n_assignments_evaluated)
            table.append(row)
    y, m = rrrm2("525")
    for a_, b_ in [("HGC", "VGC"), ("FLT", "HGC"), ("FLT", "VGC")]:
        name = f"OSD525_{a_}_vs_{b_}"
        res = panel_contrast(y, m, m.flight.isin([a_, b_]).to_numpy(), "flight", a_, panels)
        for pn, r in res.items():
            row = dict(contrast=name, panel=pn, **r)
            if (name, pn) in pkg.index:
                pr = pkg.loc[(name, pn)]
                row.update(pkg_effect=pr.mean_log2_effect, pkg_p=pr.p_animal_permutation, pkg_n_assign=pr.n_assignments_evaluated)
            table.append(row)
    t = pd.DataFrame(table)
    t.to_csv(OUT / "A_recompute_panel_effects.csv", index=False)
    have17 = t.dropna(subset=["pkg_effect"])
    record("b.all_17_design_contrasts", f"All {have17.contrast.nunique()} primary design contrasts x 8 panels ({len(have17)} rows) in panel_results.csv",
           "package panel_results.csv (design_cohort)", f"max |d effect|={(have17.effect - have17.pkg_effect).abs().max():.2e}; max |d p|={(have17.p - have17.pkg_p).abs().max():.2e}",
           "MATCH" if (have17.effect - have17.pkg_effect).abs().max() < 1e-9 and (have17.p - have17.pkg_p).abs().max() < 1e-12 else "MISMATCH")
    tt0 = t.set_index(["contrast", "panel"])
    ladder = []
    for c in ["RRRM2_pooled_GC_ISST_vs_LAR", "RRRM2_pooled_FLT_vs_GC_within_ISS_T", "RRRM2_pooled_FLT_vs_GC_within_LAR",
              "OSD613_FLT_vs_GC", "OSD563_FLT_vs_GC", "OSD564_FLT_vs_GC", "OSD612_FLT_vs_GC", "OSD525_FLT_vs_HGC", "OSD525_HGC_vs_VGC"]:
        r = tt0.loc[(c, "rPRG")]
        ladder.append(f"{c} {r.effect:+.3f} (p={r.p:.3g})")
    record("b.new.rprg_design_ladder", "NEW: rPRG effect across contrasts ordered by which protocol factors differ", "n/a", "; ".join(ladder), "NEW")
    have = t.dropna(subset=["pkg_effect"])
    have = have[have.contrast.str.startswith(("RRRM2", "OSD561", "OSD562"))]
    de = (have.effect - have.pkg_effect).abs().max()
    dp = (have.p - have.pkg_p).abs().max()
    record("b.all_panel_rows", f"All {len(have)} package design-cohort panel rows for these contrasts (effect, exact p)",
           "package panel_results.csv", f"max |d effect|={de:.2e}, max |d p|={dp:.2e}, assignments equal={bool((have.n_assign == have.pkg_n_assign).all())}",
           "MATCH" if de < 1e-9 and dp < 1e-12 else "MISMATCH")
    brief = [("RRRM2_pooled_GC_ISST_vs_LAR", "rPRG", "0.889", "0.001429"),
             ("OSD561_cerebellum_GC_ISST_vs_LAR", "rPRG", "0.704", "0.005"),
             ("OSD562_hippocampus_GC_ISST_vs_LAR", "rPRG", "1.118", "0.010"),
             ("RRRM2_pooled_FLT_ISST_vs_LAR", "rPRG", "0.656", None),
             ("RRRM2_pooled_FLT_ISST_vs_LAR", "HSP", "0.400", None),
             ("OSD561_cerebellum_FLT_ISST_vs_LAR", "rPRG", "0.332", "0.21"),
             ("OSD562_hippocampus_FLT_ISST_vs_LAR", "rPRG", "1.004", "0.0005"),
             ("OSD562_hippocampus_FLT_ISST_vs_LAR", "SRG", "0.629", None),
             ("OSD561_cerebellum_FLT_vs_GC_within_ISS_T", "rPRG", "-0.206", "0.36"),
             ("OSD562_hippocampus_FLT_vs_GC_within_ISS_T", "rPRG", "0.164", "0.49"),
             ("OSD561_cerebellum_FLT_vs_GC_within_LAR", "rPRG", "0.188", "0.38"),
             ("OSD562_hippocampus_FLT_vs_GC_within_LAR", "rPRG", "0.281", "0.06"),
             ("RRRM2_pooled_GC_ISST_vs_LAR", "SRG", "0.032", None)]
    tt = t.set_index(["contrast", "panel"])
    for c, p, e, pv in brief:
        r = tt.loc[(c, p)]
        st = at_precision(e, r.effect)
        if pv is not None:
            st2 = at_precision(pv, r.p)
            st = st if st2 == "MATCH" else (st2 if st == "MATCH" else st)
        record(f"b.{c}.{p}", f"{c} {p} effect (p)", f"{e} (p={pv})", f"{r.effect:.4f} (p={r.p:.6f}, {int(r.n_assign)} assignments, {int(r.n_genes)} genes)", st)
    r = tt.loc[("RRRM2_pooled_GC_ISST_vs_LAR", "rPRG")]
    record("b.pooled_gc_rprg_boot", "pooled GC ISS-T-LAR rPRG 95% animal-bootstrap interval 0.595-1.183 (1,000 draws)",
           "0.595-1.183", f"{r.boot_lo:.3f}-{r.boot_hi:.3f}", "DIRECTION" if abs(r.boot_lo - .595) < .06 and abs(r.boot_hi - 1.183) < .06 else "MISMATCH",
           "own seed 20261015; Monte Carlo interval, exact equality not expected")
    # new pooled flight-effect contrasts (not in package)
    for c in ("RRRM2_pooled_FLT_vs_GC_within_ISS_T", "RRRM2_pooled_FLT_vs_GC_within_LAR"):
        for p in ("rPRG", "HSP", "SRG"):
            r = tt.loc[(c, p)]
            record(f"b.new.{c}.{p}", f"NEW (no package row): {c} {p}", "n/a", f"{r.effect:.4f} (exact p={r.p:.4f}, {int(r.n_assign)} assignments)", "NEW")
    # NEW: bootstrap intervals for the pooled matched-protocol flight contrasts (2,000 whole-animal draws)
    for name, mask, field, pos in [("RRRM2_pooled_FLT_vs_GC_within_LAR", (mp.cohort == "LAR").to_numpy(), "flight", "FLT"),
                                   ("RRRM2_pooled_FLT_vs_GC_within_ISS_T", (mp.cohort == "ISS_T").to_numpy(), "flight", "FLT"),
                                   ("RRRM2_pooled_GC_ISST_vs_LAR", (mp.flight == "GC").to_numpy(), "cohort", "ISS_T")]:
        r = panel_contrast(yp, mp, mask, field, pos, ["rPRG", "HSP", "SRG"], nboot=2000, seed=15102026)
        record(f"b.new.boot.{name}", f"NEW: {name} 95% whole-animal bootstrap (2,000 draws within label x age)", "n/a",
               "; ".join(f"{k} {v['effect']:+.3f} [{v['boot_lo']:+.3f}, {v['boot_hi']:+.3f}]" for k, v in r.items()), "NEW")
    # ratio of protocol effect to flight effect, pooled rPRG
    g = tt.loc[("RRRM2_pooled_GC_ISST_vs_LAR", "rPRG")].effect
    f_lar = tt.loc[("RRRM2_pooled_FLT_vs_GC_within_LAR", "rPRG")].effect
    f_iss = tt.loc[("RRRM2_pooled_FLT_vs_GC_within_ISS_T", "rPRG")].effect
    record("b.new.protocol_vs_flight_rPRG", "NEW: pooled rPRG ground protocol effect vs matched flight effects",
           "n/a", f"GC ISS-T-LAR {g:+.3f}; FLT-GC within ISS-T {f_iss:+.3f}; within LAR {f_lar:+.3f}", "NEW")


# ----------------------------------------------------------------------------- section c

def anchor(yt, gt, z, sensor):
    q = resid_proj(z, yt.shape[1])
    y0 = yt - (yt @ q) @ q.T
    gtil = gt - q @ (q.T @ gt)
    beta = y0 @ gtil / (gtil @ gtil)
    v = beta / np.linalg.norm(beta)
    w = np.zeros_like(beta)
    w[sensor] = beta[sensor]
    w /= np.linalg.norm(w)
    s = w @ y0
    gamma = y0 @ s / (s @ s)
    return v, w, gamma


def section_c():
    data = {}
    for n in CORE:
        d = load(n)
        y, _ = normalise(d["counts"])
        data[n] = (y, d["meta"])
    geo = pd.read_csv(PKG / "results/operator_geometry.csv")
    gains = []
    for _, r in geo.iterrows():
        a, t = r.anchor_accession.split("-")[1], r.target_accession.split("-")[1]
        ya, ma = data[a]
        common = ya.index.intersection(data[t][0].index)
        if a == "525":
            mask = ma.flight.isin(["HGC", "VGC"]).to_numpy()
            lab = (ma.flight[mask] == "HGC").to_numpy(float)
            z = None
        else:
            mask = (ma.flight == "GC").to_numpy()
            lab = (ma.cohort[mask] == "ISS_T").to_numpy(float)
            z = (ma.age[mask] == "29").to_numpy(float)
        yt = ya.loc[common, ma.index[mask]].to_numpy()
        sensor = panel_index(common)["Sensor31"]
        v, w, gamma = anchor(yt, lab, z, sensor)
        # spectral norm of I - gamma w^T computed directly (2x2 block), not via the identity
        h = gamma - (w @ gamma) * w
        B = np.array([[1 - w @ gamma, 0.0], [-np.linalg.norm(h), 1.0]])
        snorm = np.linalg.svd(B, compute_uv=False)[0]
        gains.append(dict(anchor=r.anchor_accession, target=r.target_accession, n_genes=len(common), n_sensor=len(sensor),
                          w_dot_gamma=w @ gamma, gamma_norm=np.linalg.norm(gamma), spectral_norm=snorm, pkg=r.gamma_l2_norm))
        record(f"c.gain.{r.anchor_accession}->{r.target_accession}", "score-operator spectral gain ||gamma||",
               f"{r.gamma_l2_norm:.3f}", f"{np.linalg.norm(gamma):.3f} (||A||2={snorm:.3f}, w.gamma={w @ gamma:.6f}, genes {len(common)}/{r.n_genes}, sensors {len(sensor)}/{r.n_sensor})",
               at_precision(f"{r.gamma_l2_norm:.3f}", np.linalg.norm(gamma)))
    g = pd.DataFrame(gains)
    g.to_csv(OUT / "A_recompute_operator_gains.csv", index=False)
    record("c.range", "learned spectral gains ranged from 8.129 to 15.309 across all 15 geometries", "8.129-15.309",
           f"{g.gamma_norm.min():.3f}-{g.gamma_norm.max():.3f} (n={len(g)})",
           "MATCH" if at_precision("8.129", g.gamma_norm.min()) == "MATCH" and at_precision("15.309", g.gamma_norm.max()) == "MATCH" else "MISMATCH")
    # also bounded loading at L=2
    ya, ma = data["561"]
    mask = (ma.flight == "GC").to_numpy()
    sensor = panel_index(ya.index)["Sensor31"]
    v, w, gamma = anchor(ya.loc[:, ma.index[mask]].to_numpy(), (ma.cohort[mask] == "ISS_T").to_numpy(float),
                         (ma.age[mask] == "29").to_numpy(float), sensor)
    h = gamma - w
    c = min(1, math.sqrt(3) / np.linalg.norm(h))
    gl = w + c * h
    A = np.eye(3)
    record("c.bounded_561", "bounded loading at L=2 has norm 2 and keeps w.gamma=1 (OSD-561 anchor)", "2, 1",
           f"{np.linalg.norm(gl):.12f}, {w @ gl:.12f}", "MATCH" if abs(np.linalg.norm(gl) - 2) < 1e-9 and abs(w @ gl - 1) < 1e-9 else "MISMATCH")


# ----------------------------------------------------------------------------- section d

def train_ref_expression(c, train, test):
    canon = c.index.astype(str).str.fullmatch(r"ENSMUSG\d+")
    keep = canon & (c.iloc[:, train].mean(axis=1).to_numpy() >= 10)
    a = c.loc[keep].to_numpy(float)
    pos = (a[:, train] > 0).all(axis=1)
    logref = np.log(a[pos][:, train]).mean(axis=1)

    def expr(cols):
        f = []
        for j in cols:
            x = a[pos, j]
            nz = x > 0
            f.append(np.exp(np.median(np.log(x[nz]) - logref[nz])))
        f = np.array(f)
        return np.log2(a[:, cols] / f + 1)
    return c.index[keep], expr(train), expr(test)


def strat_auc(score, age, lab):
    wins = []
    for a in np.unique(age):
        s = age == a
        P, N = score[s & (lab == 1)], score[s & (lab == 0)]
        wins += [(p > q) + .5 * (p == q) for p in P for q in N]
    return float(np.mean(wins))


def section_d():
    summ = pd.read_csv(PKG / "extensions/ground_holdout/holdout_summary.csv").set_index(["accession", "method"])
    per_split_pkg = pd.read_csv(PKG / "extensions/ground_holdout/holdout_split_results.csv")
    for n in ("561", "562"):
        d = load(n)
        m = d["meta"][d["meta"].flight == "GC"].copy()
        c = d["counts"].loc[:, m.index]
        m["label"] = (m.cohort == "ISS_T").astype(int)
        mr = m.reset_index(drop=True)
        cells = [np.asarray(v) for _, v in mr.groupby(["label", "age"], sort=True).indices.items()]
        rows = []
        for sid, held in enumerate(itertools.product(*cells)):
            test = np.array(held)
            train = np.array([i for i in range(len(mr)) if i not in held])
            genes, yt, ye = train_ref_expression(c, train, test)
            mt, me = mr.iloc[train], mr.iloc[test]
            zt = (mt.age == "29").to_numpy(float)
            ze = (me.age == "29").to_numpy(float)
            sensor = panel_index(genes)["Sensor31"]
            outside = np.setdiff1d(np.arange(len(genes)), sensor)
            gt, ge = mt.label.to_numpy(float), me.label.to_numpy(float)
            v, w, gamma = anchor(yt, gt, zt, sensor)
            X = np.column_stack([np.ones(len(train)), zt])
            B = np.linalg.lstsq(X, yt.T, rcond=None)[0]
            res = ye - (np.column_stack([np.ones(len(test)), ze]) @ B).T
            score = w @ res
            wscore = v @ res
            qt, qe = resid_proj(zt, len(train)), resid_proj(ze, len(test))
            bt = fwl_beta(yt, gt, qt)
            be = fwl_beta(ye, ge, qe)
            cos = bt @ be / np.linalg.norm(bt) / np.linalg.norm(be)
            cos_s = bt[sensor] @ be[sensor] / np.linalg.norm(bt[sensor]) / np.linalg.norm(be[sensor])
            r = res[outside]
            mse0 = np.mean(r ** 2)
            mse_score = np.mean((r - np.outer(gamma[outside], score)) ** 2) / mse0
            corrected = ye - np.outer(gamma, w @ ye)
            bc = fwl_beta(corrected, ge, qe)
            energy = 100 * (1 - np.sum(bc[outside] ** 2) / np.sum(be[outside] ** 2))
            rows.append(dict(split=sid, auc=strat_auc(score, me.age.to_numpy(), me.label.to_numpy()),
                             auc_whole=strat_auc(wscore, me.age.to_numpy(), me.label.to_numpy()),
                             cosine=cos, cosine_sensor=cos_s, gain=np.linalg.norm(gamma), mse_ratio_score=mse_score,
                             outside_energy_removed=energy, n_genes=len(genes), n_sensor=len(sensor)))
        df = pd.DataFrame(rows)
        df.to_csv(OUT / f"A_recompute_holdout_OSD-{n}.csv", index=False)
        pk = per_split_pkg[(per_split_pkg.accession == f"OSD-{n}") & (per_split_pkg.method == "score")].sort_values("split_id")
        dcos = np.abs(df.cosine.to_numpy() - pk.train_test_contrast_cosine.to_numpy()).max()
        dmse = np.abs(df.mse_ratio_score.to_numpy() - pk.outside_prediction_MSE_ratio_to_nuisance_only.to_numpy()).max()
        s = summ.loc[(f"OSD-{n}", "score")]
        exp_cos = {"561": "0.228", "562": "0.521"}[n]
        exp_mse = {"561": "1.091", "562": "0.962"}[n]
        exp_en = {"561": "-23.26", "562": "16.85"}[n]
        record(f"d.auc.OSD-{n}", f"OSD-{n} mean age-stratified sensor AUC over 81 splits", "1.000",
               f"{df.auc.mean():.3f} (min {df.auc.min():.3f}; {len(df)} splits; {df.n_sensor.iloc[0]} sensors)",
               at_precision("1.000", df.auc.mean()))
        record(f"d.cos.OSD-{n}", f"OSD-{n} mean training-test whole-genome contrast cosine", exp_cos,
               f"{df.cosine.mean():.4f} (max per-split |diff| vs package {dcos:.1e})", at_precision(exp_cos, df.cosine.mean()))
        record(f"d.mse.OSD-{n}", f"OSD-{n} original-score outside-sensor MSE ratio", exp_mse,
               f"{df.mse_ratio_score.mean():.4f} (max per-split |diff| {dmse:.1e})", at_precision(exp_mse, df.mse_ratio_score.mean()))
        record(f"d.energy.OSD-{n}", f"OSD-{n} mean held-out outside-sensor contrast energy removal (%)", exp_en,
               f"{df.outside_energy_removed.mean():.2f}", at_precision(exp_en, df.outside_energy_removed.mean()))
        exp_wauc = {"561": "0.796", "562": "1.000"}[n]
        record(f"d.wauc.OSD-{n}", f"OSD-{n} mean whole-genome AUC (Supplement)", exp_wauc, f"{df.auc_whole.mean():.3f}",
               at_precision(exp_wauc, df.auc_whole.mean()))
        exp_sc = {"561": "0.683", "562": "0.843"}[n]
        record(f"d.scos.OSD-{n}", f"OSD-{n} mean sensor-vector cosine (Supplement)", exp_sc, f"{df.cosine_sensor.mean():.3f}",
               at_precision(exp_sc, df.cosine_sensor.mean()))


# ----------------------------------------------------------------------------- section e

def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def section_e():
    design = json.loads((PKG / "results/regional_gene_replication_design.json").read_text())["contrasts"]
    cm = {}
    for f in sorted((PKG / "results/count_models").glob("OSD-*_gene_results.csv.gz")):
        g = pd.read_csv(f, usecols=["contrast", "ensembl_id", "log2FoldChange", "pvalue"])
        for c, s in g.groupby("contrast"):
            cm[c] = s.set_index("ensembl_id")
    allp, per_region = [], {}
    for region, comps in design.items():
        genes = None
        for c in comps:
            genes = cm[c].index if genes is None else genes.intersection(cm[c].index)
        L = np.column_stack([cm[c].loc[genes, "log2FoldChange"].to_numpy() for c in comps])
        P = np.column_stack([cm[c].loc[genes, "pvalue"].to_numpy() for c in comps])
        signs = np.sign(L)
        same_all = (np.abs(signs.sum(axis=1)) == len(comps)) & np.isfinite(L).all(axis=1)
        valid = same_all & np.isfinite(P).all(axis=1)
        maxp = np.where(valid, np.nanmax(np.where(np.isfinite(P), P, np.nan), axis=1), 1.0)
        per_region[region] = dict(n=len(genes), concordant=int(same_all.sum()), nominal=int((maxp < .05).sum()), minp=float(maxp.min()))
        allp.append(maxp)
    allp = np.concatenate(allp)
    q = bh(allp)
    hsum = np.sum(1 / np.arange(1, len(allp) + 1))
    qby = np.minimum(q * hsum, 1)
    record("e.n_tests", "54,989 region-gene combinations", 54989, len(allp), "MATCH" if len(allp) == 54989 else "MISMATCH")
    record("e.survivors", "no gene survived global BH or BY", "0 / 0", f"{int((q <= .05).sum())} / {int((qby <= .05).sum())} (min BH q={q.min():.3f})",
           "MATCH" if (q <= .05).sum() == 0 and (qby <= .05).sum() == 0 else "MISMATCH")
    exp = {"cerebellum": (5047, 3), "hippocampus": (5001, 0), "left_cerebral_hemisphere": (9884, 554)}
    for r, v in per_region.items():
        e = exp[r]
        record(f"e.{r}", f"{r}: direction-concordant genes / nominal max-p<0.05", f"{e[0]} / {e[1]}",
               f"{v['concordant']} / {v['nominal']} (n={v['n']}, min max-p={v['minp']:.2e})",
               "MATCH" if (v["concordant"], v["nominal"]) == e else "MISMATCH")
    # Estimator-independent sensitivity: same IU screen on own log-expression OLS for the left hemisphere
    comps = {"OSD613": "613", "OSD612": "612"}
    effs, ps = [], []
    genes = None
    fits = {}
    for c, n in comps.items():
        d = load(n)
        y, _ = normalise(d["counts"])
        m = d["meta"]
        z = (m.age == "29").to_numpy(float) if m.age.nunique() > 1 else None
        q_ = resid_proj(z, len(m))
        g = (m.flight == "FLT").to_numpy(float)
        Y = y.to_numpy()
        y0 = Y - (Y @ q_) @ q_.T
        gt = g - q_ @ (q_.T @ g)
        beta = y0 @ gt / (gt @ gt)
        dfree = len(m) - q_.shape[1] - 1
        rss = np.sum(y0 ** 2, axis=1) - beta ** 2 * (gt @ gt)
        se = np.sqrt(np.maximum(rss, 0) / dfree / (gt @ gt))
        p = 2 * stats.t.sf(np.abs(beta / se), dfree)
        fits[c] = pd.DataFrame(dict(b=beta, p=p), index=y.index)
        genes = y.index if genes is None else genes.intersection(y.index)
    B = np.column_stack([fits[c].loc[genes, "b"] for c in comps])
    P = np.column_stack([fits[c].loc[genes, "p"] for c in comps])
    same = np.abs(np.sign(B).sum(axis=1)) == 2
    mp = np.where(same, P.max(axis=1), 1.0)
    qq = bh(mp)
    record("e.ols_left_hemisphere", "Estimator sensitivity (NEW): IU screen with own log-expression OLS t tests, left hemisphere (OSD-613 + OSD-612)",
           "package NB screen: 9,884 concordant, 554 nominal, 0 BH", f"{int(same.sum())} concordant, {int((mp < .05).sum())} nominal, {int((qq <= .05).sum())} BH within region (n={len(genes)})", "NEW")


# ----------------------------------------------------------------------------- section f

def deseq_size_factors(counts):
    a = counts.to_numpy(float)
    pos = (a > 0).all(axis=1)
    la = np.log(a[pos])
    return np.exp(np.median(la - la.mean(axis=1, keepdims=True), axis=0))


def nb_fit(yc, X, sf):
    """NB2 GLM MLE (log link, offset log sf) with Cox-Reid adjusted profile dispersion. Own code."""
    off = np.log(sf)

    def irls(alpha, b0=None):
        b = np.linalg.lstsq(X, np.log(yc / sf + .5), rcond=None)[0] if b0 is None else b0.copy()
        for _ in range(200):
            mu = np.exp(X @ b + off)
            wts = mu / (1 + alpha * mu)
            zz = X @ b + (yc - mu) / mu
            bn = np.linalg.solve(X.T @ (wts[:, None] * X), X.T @ (wts * zz))
            if np.max(np.abs(bn - b)) < 1e-10:
                b = bn
                break
            b = bn
        return b

    def negll(la):
        alpha = math.exp(la)
        b = irls(alpha)
        mu = np.exp(X @ b + off)
        r = 1 / alpha
        ll = np.sum(special.gammaln(yc + r) - special.gammaln(r) - special.gammaln(yc + 1)
                    + r * np.log(r / (r + mu)) + yc * np.log(mu / (r + mu)))
        wts = mu / (1 + alpha * mu)
        cr = 0.5 * np.linalg.slogdet(X.T @ (wts[:, None] * X))[1]
        return -(ll - cr)

    opt = optimize.minimize_scalar(negll, bounds=(-12, 3), method="bounded")
    alpha = math.exp(opt.x)
    b = irls(alpha)
    mu = np.exp(X @ b + off)
    wts = mu / (1 + alpha * mu)
    cov = np.linalg.inv(X.T @ (wts[:, None] * X))
    return b, np.sqrt(np.diag(cov)), alpha


def interaction_design(m):
    f = (m.flight == "FLT").to_numpy(float)
    i = (m.cohort == "ISS_T").to_numpy(float)
    return np.column_stack([np.ones(len(m)), (m.age == "29").to_numpy(float), f, i, f * i])


def section_f():
    d = load("561")
    m = d["meta"]
    c = d["counts"]
    canon = c.index.astype(str).str.fullmatch(r"ENSMUSG\d+")
    cf = c.loc[canon & (c.mean(axis=1) >= 10)]
    sf = deseq_size_factors(cf)
    fitjs = json.loads((PKG / "extensions/flight_collection_interaction/OSD-561_all_animals_NB_fit.json").read_text())
    record("f.size_factors", "PyDESeq2 size-factor range (OSD-561 fit json)", f"{fitjs['size_factor_range'][0]:.4f}-{fitjs['size_factor_range'][1]:.4f}",
           f"{sf.min():.4f}-{sf.max():.4f}", at_precision(f"{fitjs['size_factor_range'][1]:.4f}", sf.max()))
    X = interaction_design(m)
    y, _ = normalise(c)
    for sym, ens, exp, expq in [("Hbb-bt", "ENSMUSG00000073940", "-3.495", "0.000197"), ("Hbb-bs", "ENSMUSG00000052305", "-3.009", "0.001406")]:
        yc = cf.loc[ens].to_numpy(float)
        b, se, alpha = nb_fit(yc, X, sf)
        lfc, lse = b[4] / math.log(2), se[4] / math.log(2)
        p = 2 * stats.norm.sf(abs(b[4] / se[4]))
        ols = np.linalg.lstsq(X, y.loc[ens].to_numpy(), rcond=None)[0][4]
        record(f"f.{sym}", f"cerebellar {sym} flight x collection interaction (NB, log2)", f"{exp} (global q={expq})",
               f"NB {lfc:.3f} (SE {lse:.3f}, Wald p={p:.2e}, own CR dispersion {alpha:.3f}); own OLS log2 {ols:.3f}",
               at_precision(exp, lfc) if at_precision(exp, lfc) != "MISMATCH" else ("DIRECTION" if np.sign(lfc) == np.sign(float(exp)) else "MISMATCH"),
               "independent NB MLE with gene-wise Cox-Reid dispersion (no empirical-Bayes shrinkage); q not recomputed (needs all 72,023 fits)")


# ----------------------------------------------------------------------------- section g

def hc3_t(X, Y, k):
    """HC3 t for coefficient k; Y can be n x B. Returns beta_k, se_k, t_k."""
    XtXi = np.linalg.inv(X.T @ X)
    H = X @ XtXi
    h = np.sum(H * X, axis=1)
    B = XtXi @ X.T @ Y
    E = Y - X @ B
    a = XtXi[k] @ X.T                      # row giving beta_k = a @ y
    var = np.sum((a[:, None] ** 2) * (E ** 2) / (1 - h[:, None]) ** 2, axis=0)
    return B[k], np.sqrt(var), B[k] / np.sqrt(var)


def section_g():
    d = load("562")
    y, _ = normalise(d["counts"])
    m = d["meta"]
    idx = panel_index(y.index)["Circadian"]
    pm = y.to_numpy()[idx].mean(axis=0)
    X = interaction_design(m)
    beta, se, t = hc3_t(X, pm[:, None], 4)
    record("g.estimate", "hippocampal Circadian panel interaction (log2)", "-0.405", f"{beta[0]:.4f} ({len(idx)} genes, {len(m)} animals)",
           at_precision("-0.405", beta[0]))
    # stratified whole-animal bootstrap, 4000 draws, flight x cohort x age strata
    rng = np.random.default_rng(15102026)
    strata = [np.flatnonzero(((m.flight == f) & (m.cohort == c) & (m.age == a)).to_numpy())
              for f in ("FLT", "GC") for c in ("ISS_T", "LAR") for a in ("12", "29")]
    strata = [s for s in strata if len(s)]
    bs = []
    for _ in range(4000):
        ix = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
        Xb = X[ix]
        if np.linalg.matrix_rank(Xb) < 5:
            continue
        bs.append(np.linalg.lstsq(Xb, pm[ix], rcond=None)[0][4])
    lo, hi = np.quantile(bs, [.025, .975])
    record("g.interval", "95% stratified animal-bootstrap interval -0.591 to -0.215 (4,000 draws)", "-0.591 to -0.215",
           f"{lo:.3f} to {hi:.3f} ({len(bs)} usable draws, own seed)",
           "DIRECTION" if abs(lo + .591) < .04 and abs(hi + .215) < .04 else "MISMATCH", "Monte Carlo interval; exact equality not expected")
    # HC3-studentized restricted-model Rademacher wild bootstrap, 499,999 draws
    Xr = X[:, :4]
    br = np.linalg.lstsq(Xr, pm, rcond=None)[0]
    fit_r = Xr @ br
    er = pm - fit_r
    rng = np.random.default_rng(20261015)
    B = 499_999
    exceed = 0
    done = 0
    while done < B:
        k = min(50_000, B - done)
        s = rng.choice([-1.0, 1.0], size=(len(pm), k))
        Ys = fit_r[:, None] + er[:, None] * s
        _, _, ts = hc3_t(X, Ys, 4)
        exceed += int(np.sum(np.abs(ts) >= abs(t[0]) - 1e-12))
        done += k
    p = (exceed + 1) / (B + 1)
    record("g.wild_p", "restricted-model HC3 wild bootstrap p (499,999 draws)", "0.000838", f"{p:.6f} (HC3 t={t[0]:.3f})",
           "DIRECTION" if 0.0006 < p < 0.0011 else "MISMATCH", "Monte Carlo p; package MC 95% interval 0.000758-0.000920")
    # metadata exclusion
    conf = set()
    for n in ("561", "562"):
        mm = load(n)["meta"]
        conf |= set(mm.animal[mm.conflict])
    keep = ~m.animal.isin(conf).to_numpy()
    be = np.linalg.lstsq(X[keep], pm[keep], rcond=None)[0][4]
    record("g.exclusion", "metadata-exclusion hippocampal Circadian interaction", "-0.460", f"{be:.4f} ({keep.sum()} animals)", at_precision("-0.460", be))
    # leave-one-animal-out range
    loo = [np.linalg.lstsq(np.delete(X, i, 0), np.delete(pm, i), rcond=None)[0][4] for i in range(len(pm))
           if np.linalg.matrix_rank(np.delete(X, i, 0)) == 5]
    record("g.loo_animal", "one-animal omission range -0.459 to -0.365", "-0.459 to -0.365", f"{min(loo):.3f} to {max(loo):.3f}",
           "MATCH" if at_precision("-0.459", min(loo)) == "MATCH" and at_precision("-0.365", max(loo)) == "MATCH" else "MISMATCH")
    Yg = y.to_numpy()[idx]
    log_ = [np.linalg.lstsq(X, np.delete(Yg, j, 0).mean(axis=0), rcond=None)[0][4] for j in range(len(idx))]
    record("g.loo_gene", "one-gene omission range -0.452 to -0.360", "-0.452 to -0.360", f"{min(log_):.3f} to {max(log_):.3f}",
           "MATCH" if at_precision("-0.452", min(log_)) == "MATCH" and at_precision("-0.360", max(log_)) == "MATCH" else "MISMATCH")
    # cells
    cell = m.groupby(["flight", "cohort", "age"]).size().to_dict()
    record("g.cells", "hippocampal design-cell sizes (flight, cohort, age)", "n/a", json.dumps({"|".join(k): int(v) for k, v in cell.items()}), "NEW")


# ----------------------------------------------------------------------------- section h

def section_h():
    ext = PKG / "extensions/strengthening_20261015/external_data"
    tot_cols, indiv, pooled = 0, 0, 0
    for n, sample_table in [("202", "GLDS-202"), ("352", "GLDS-352"), ("374", "GLDS-374"), ("611", "GLDS-592")]:
        raw = pd.read_csv(ext / f"OSD-{n}/raw_sources/{sample_table}_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv", index_col=0, nrows=2)
        st = pd.read_csv(ext / f"OSD-{n}/raw_sources/{sample_table}_rna_seq_SampleTable_GLbulkRNAseq.csv")
        stid = set(st.iloc[:, 0].astype(str)) | set(st.iloc[:, 1].astype(str)) if st.shape[1] > 1 else set(st.iloc[:, 0].astype(str))
        cols = list(raw.columns)
        in_st = [c for c in cols if c in stid]
        md = pd.read_csv(ext / f"OSD-{n}/metadata.csv", dtype=str)
        pool = md.pool_size_protocol.astype(int).max()
        tot_cols += len(cols)
        units = len(in_st)
        if pool == 1:
            indiv += units
        else:
            pooled += units
        record(f"h.OSD-{n}", f"OSD-{n} count columns / columns in NASA SampleTable / pool size", "see Table 2",
               f"{len(cols)} / {units} / {pool}", "INFO", f"SampleTable rows {len(st)}")
    record("h.total_columns", "66 acquired count columns", 66, tot_cols, "MATCH" if tot_cols == 66 else "MISMATCH")
    record("h.units", "64 primary units = 46 individual + 18 pooled", "64 = 46 + 18", f"{indiv + pooled} = {indiv} + {pooled}",
           "MATCH" if (indiv, pooled) == (46, 18) else "MISMATCH")


# ----------------------------------------------------------------------------- section i

def section_i():
    mb = PKG / "extensions/strengthening_20261015/method_benchmark"
    o = pd.read_csv(mb / "strengthened_outer_results.csv")
    o = o[o.method == "spectral_constrained_ridge_CV"]
    exp = {"OSD-561": ("0.9878", "0.9884", "78/81"), "OSD-562": ("0.9484", "0.9492", "81/81"), "OSD-525": ("0.9628", "0.9628", "200/200"),
           "OSD-202": ("0.9284", "0.9080", "37/40"), "OSD-352": ("0.9857", "0.9843", "5/6"), "OSD-374": ("0.9501", "0.9744", "8/8"),
           "OSD-611": ("0.9719", "0.9836", "9/10")}
    red = {}
    for acc, (em, ep, ei) in exp.items():
        s = o[o.accession == acc]
        ratio = s.prediction_SSE / s.baseline_SSE
        consistent = np.allclose(ratio, s.MSE_ratio, rtol=1e-9)
        mean = ratio.mean()
        pooled = s.prediction_SSE.sum() / s.baseline_SSE.sum()
        imp = f"{int((ratio < 1).sum())}/{len(s)}"
        red[acc] = 100 * (1 - mean)
        st = "MATCH" if (at_precision(em, mean) == "MATCH" and at_precision(ep, pooled) == "MATCH" and imp == ei) else "MISMATCH"
        record(f"i.{acc}", f"{acc} constrained mean MSE ratio; pooled SSE ratio; improved holdouts", f"{em}; {ep}; {ei}",
               f"{mean:.4f}; {pooled:.4f}; {imp}; max gain {s.operator_gain.max():.6f}", st,
               f"per-row MSE_ratio==SSE ratio: {consistent}; n rows {len(s)}")
    orig = [red[a] for a in ("OSD-561", "OSD-562", "OSD-525")]
    ext = [red[a] for a in ("OSD-202", "OSD-352", "OSD-374", "OSD-611")]
    record("i.range_original", "mean held-out error reduction 1.2-5.2% (original settings)", "1.2-5.2",
           f"{min(orig):.2f}-{max(orig):.2f}", "MATCH" if (f"{min(orig):.1f}", f"{max(orig):.1f}") == ("1.2", "5.2") else "MISMATCH")
    record("i.range_external", "mean held-out error reduction 1.4-7.2% (external)", "1.4-7.2",
           f"{min(ext):.2f}-{max(ext):.2f}", "MATCH" if (f"{min(ext):.1f}", f"{max(ext):.1f}") == ("1.4", "7.2") else "MISMATCH")
    # equal-weight vs pooled for external: the paper reports equal-unit means
    nunits = o[o.accession.isin(["OSD-202", "OSD-352", "OSD-374", "OSD-611"])].shape[0]
    record("i.external_units", "external holdouts = 64 animal/pool units", 64, nunits, "MATCH" if nunits == 64 else "MISMATCH")
    # comparison with unconstrained ridge from the same file
    o2 = pd.read_csv(mb / "strengthened_outer_results.csv")
    for acc in ("OSD-202", "OSD-374"):
        for meth in ("ridge_CV", "original_score"):
            s = o2[(o2.accession == acc) & (o2.method == meth)]
            if len(s):
                record(f"i.{acc}.{meth}", f"{acc} {meth} mean MSE ratio (paper: ridge OSD-202 0.8028; original score OSD-374 0.7344)", "see text",
                       f"{(s.prediction_SSE / s.baseline_SSE).mean():.4f}; mean gain {s.operator_gain.mean():.2f}", "INFO")


# ----------------------------------------------------------------------------- section k (new)

def section_k():
    """NEW: is the 8-15 spectral gain specific to the IEG/HSP sensor?  Random 30-gene sensors.

    For the same anchor (ground ISS-T vs LAR, age-adjusted, n=12 animals) draw 500 random
    gene sets of the same size as Sensor31 from the retained genes, form w exactly as for
    Sensor31 (restricted ground contrast, unit length) and record ||gamma||.  Also report the
    fixed-sensor isotropic-noise expectation sqrt(1+(G-1)/(r-2)) from the paper's Eq. 11.
    """
    rng = np.random.default_rng(561562015)
    for n in ("561", "562", "525"):
        d = load(n)
        y, _ = normalise(d["counts"])
        m = d["meta"]
        if n == "525":
            mask = m.flight.isin(["HGC", "VGC"]).to_numpy()
            lab = (m.flight[mask] == "HGC").to_numpy(float)
            z = None
        else:
            mask = (m.flight == "GC").to_numpy()
            lab = (m.cohort[mask] == "ISS_T").to_numpy(float)
            z = (m.age[mask] == "29").to_numpy(float)
        Y = y.loc[:, m.index[mask]].to_numpy()
        sensor = panel_index(y.index)["Sensor31"]
        g_obs = np.linalg.norm(anchor(Y, lab, z, sensor)[2])
        q = resid_proj(z, Y.shape[1])
        y0 = Y - (Y @ q) @ q.T
        gtil = lab - q @ (q.T @ lab)
        beta = y0 @ gtil / (gtil @ gtil)
        null = []
        for _ in range(500):
            idx = rng.choice(len(beta), len(sensor), replace=False)
            w = np.zeros_like(beta)
            w[idx] = beta[idx]
            w /= np.linalg.norm(w)
            s = w @ y0
            null.append(np.linalg.norm(y0 @ s / (s @ s)))
        null = np.array(null)
        r = Y.shape[1] - q.shape[1]
        G = len(beta)
        iso = math.sqrt(1 + (G - 1) / (r - 2)) if r > 2 else float("nan")
        record(f"k.gain_null.OSD-{n}", f"NEW: Sensor31 gain vs 500 random {len(sensor)}-gene sensors (OSD-{n} anchor, own universe)",
               "n/a", f"Sensor31 {g_obs:.2f}; random median {np.median(null):.2f} (2.5-97.5%: {np.quantile(null, .025):.2f}-{np.quantile(null, .975):.2f}); "
               f"fraction random >= Sensor31: {np.mean(null >= g_obs):.3f}; isotropic fixed-w expectation (G={G}, r={r}) {iso:.1f}", "NEW")


# ----------------------------------------------------------------------------- section j (new)

EXT = PKG / "extensions/strengthening_20261015/external_data"
SENSOR31 = set(PANELS["Sensor31"])
ALPHAS = [None, 10., 1., .1, .01]


def ext_load(n):
    st_name = {"202": "GLDS-202", "352": "GLDS-352", "374": "GLDS-374", "611": "GLDS-592"}[n]
    raw = pd.read_csv(EXT / f"OSD-{n}/raw_sources/{st_name}_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv", index_col=0)
    st = pd.read_csv(EXT / f"OSD-{n}/raw_sources/{st_name}_rna_seq_SampleTable_GLbulkRNAseq.csv", index_col=0)
    raw = raw.loc[:, [c for c in raw.columns if c in st.index]]
    md = pd.read_csv(EXT / f"OSD-{n}/metadata.csv", dtype=str).set_index("sample_id").loc[raw.columns]
    X = [np.ones(raw.shape[1])]
    if n == "202":
        X += [(md.collection_after_treatment.astype(float) == 4).to_numpy(float), (md.hindlimb_unloading == "Hindlimb Unloaded").to_numpy(float)]
    return raw, np.column_stack(X)


def ext_build(c, genes, X, train):
    canon = genes.str.fullmatch(r"ENSMUSG\d+")
    keep = np.flatnonzero(canon & (c[:, train].mean(axis=1) >= 10))
    a = c[keep][:, train]
    pos = keep[(a > 0).all(axis=1)]
    logref = np.log(c[pos][:, train]).mean(axis=1)
    sf = np.exp(np.median(np.log(c[pos][:, train]) - logref[:, None], axis=0))
    y = np.log2(c[keep][:, train] / sf + 1).T
    B = np.linalg.lstsq(X[train], y, rcond=None)[0]
    res = y - X[train] @ B
    sym = np.array([ANNOT.get(g, "") for g in genes[keep]])
    sens = np.flatnonzero(np.isin(sym, list(SENSOR31)))
    outs = np.setdiff1d(np.arange(len(keep)), sens)
    S, Yo = res[:, sens], res[:, outs]
    scale = np.sqrt(np.mean(S * S, axis=0))
    scale = np.where(scale < 1e-10, 1.0, scale)
    return dict(keep=keep, pos=pos, logref=logref, B=B, sens=sens, outs=outs, S=S, Y=Yo, scale=scale)


def ext_target(c, X, mdl, test):
    out = []
    for j in test:
        x = c[mdl["pos"], j]
        nz = x > 0
        f = np.exp(np.median(np.log(x[nz]) - mdl["logref"][nz]))
        yv = np.log2(c[mdl["keep"], j] / f + 1) - X[j] @ mdl["B"]
        out.append(yv)
    r = np.vstack(out)
    return r[:, mdl["sens"]], r[:, mdl["outs"]]


def spectral_ridge(S, Y, scale, alpha, cap=math.sqrt(3), iters=200000, rel_gap=1e-10):
    """Own solver for min ||Y-SC||_F^2/(2n) + alpha||DC||_F^2/2 s.t. ||C||_2 <= cap.

    Written through the thin SVD of Y (exact: C = T V^T, T = p x r).  Unconstrained ridge
    solution is returned when feasible; otherwise accelerated projected gradient with
    gradient restart, stopped on a Frank-Wolfe duality gap (own implementation)."""
    n = S.shape[0]
    U, s, Vt = np.linalg.svd(Y, full_matrices=False)
    k = s > s.max() * 1e-10
    U, s, Vt = U[:, k], s[k], Vt[k]
    H = S.T @ S / n + alpha * np.diag(scale ** 2)
    F = S.T @ (U * s) / n
    T0 = np.linalg.solve(H, F)
    if np.linalg.svd(T0, compute_uv=False)[0] <= cap:
        return T0, Vt
    const = np.sum(s ** 2) / (2 * n)
    L = np.linalg.eigvalsh(H)[-1]

    def proj(A):
        u, sv, vt = np.linalg.svd(A, full_matrices=False)
        return (u * np.minimum(sv, cap)) @ vt
    T = proj(T0)
    Z, tk = T.copy(), 1.0
    for it in range(iters):
        Tn = proj(Z - (H @ Z - F) / L)
        G = H @ Tn - F
        gap = np.sum(G * Tn) + cap * np.linalg.svd(G, compute_uv=False).sum()
        if gap <= rel_gap * max(1.0, const):
            T = Tn
            break
        if np.sum((Z - Tn) * (Tn - T)) > 0:      # gradient restart
            tk = 1.0
        tn = (1 + math.sqrt(1 + 4 * tk * tk)) / 2
        Z = Tn + (tk - 1) / tn * (Tn - T)
        T, tk = Tn, tn
    return T, Vt


def section_j():
    pkg = pd.read_csv(PKG / "extensions/strengthening_20261015/method_benchmark/strengthened_outer_results.csv")
    pkg = pkg[pkg.method == "spectral_constrained_ridge_CV"]
    for n in ("352", "374", "611", "202"):
        raw, X = ext_load(n)
        c = raw.to_numpy(float)
        genes = raw.index.astype(str)
        N = c.shape[1]
        ratios, chosen, gains = [], [], []
        for o in range(N):
            train = np.array([i for i in range(N) if i != o])
            err = np.zeros(len(ALPHAS))
            for left in train:
                inner = np.array([i for i in train if i != left])
                mdl = ext_build(c, genes, X, inner)
                St, Rt = ext_target(c, X, mdl, [left])
                for a_i, a in enumerate(ALPHAS):
                    if a is None:
                        err[a_i] += np.sum(Rt ** 2)
                    else:
                        T, Vt = spectral_ridge(mdl["S"], mdl["Y"], mdl["scale"], a)
                        err[a_i] += np.sum((Rt - (St @ T) @ Vt) ** 2)
            a = ALPHAS[int(np.argmin(err))]
            mdl = ext_build(c, genes, X, train)
            St, Rt = ext_target(c, X, mdl, [o])
            if a is None:
                sse, g = np.sum(Rt ** 2), 1.0
            else:
                T, Vt = spectral_ridge(mdl["S"], mdl["Y"], mdl["scale"], a)
                sse = np.sum((Rt - (St @ T) @ Vt) ** 2)
                g = math.sqrt(1 + np.linalg.svd(T, compute_uv=False)[0] ** 2)
            ratios.append(sse / np.sum(Rt ** 2))
            chosen.append(a)
            gains.append(g)
        ratios = np.array(ratios)
        pk = pkg[pkg.accession == f"OSD-{n}"].sort_values("split_id")
        diff = np.abs(np.sort(ratios) - np.sort(pk.MSE_ratio.to_numpy())).max()
        exp = {"352": "0.9857", "374": "0.9501", "611": "0.9719", "202": "0.9284"}[n]
        record(f"j.OSD-{n}", f"OSD-{n} nested gain-2 constrained regression, own implementation from raw STAR counts + SampleTable",
               f"{exp} (package {len(pk)} rows; alphas {pk.parameter.tolist()})",
               f"mean MSE ratio {ratios.mean():.4f}; improved {int((ratios < 1).sum())}/{N}; max gain {max(gains):.6f}; alphas {chosen}; max |sorted diff| {diff:.1e}",
               at_precision(exp, ratios.mean()))


# ----------------------------------------------------------------------------- section m (new)

def section_m():
    """Own nested gain-2 constrained regression on the ORIGINAL ground holdouts (OSD-561, OSD-562; 81 splits
    each, enumerated independently), from collapsed raw counts, training-only preprocessing."""
    pkg = pd.read_csv(PKG / "extensions/strengthening_20261015/method_benchmark/strengthened_outer_results.csv")
    pkg = pkg[pkg.method == "spectral_constrained_ridge_CV"]
    for n, exp in (("561", "0.9878"), ("562", "0.9484")):
        d = load(n)
        m = d["meta"][d["meta"].flight == "GC"].copy()
        m["label"] = (m.cohort == "ISS_T").astype(int)
        mr = m.reset_index(drop=True)
        c = d["counts"].loc[:, m.index].to_numpy(float)
        genes = d["counts"].index.astype(str)
        X = np.column_stack([np.ones(len(mr)), (mr.age == "29").to_numpy(float)])
        cells = [np.asarray(v) for _, v in mr.groupby(["label", "age"], sort=True).indices.items()]
        ratios = []
        for held in itertools.product(*cells):
            test = np.array(held)
            train = np.array([i for i in range(len(mr)) if i not in held])
            err = np.zeros(len(ALPHAS))
            for left in train:
                inner = np.array([i for i in train if i != left])
                mdl = ext_build(c, genes, X, inner)
                St, Rt = ext_target(c, X, mdl, [left])
                for a_i, a in enumerate(ALPHAS):
                    if a is None:
                        err[a_i] += np.sum(Rt ** 2)
                    else:
                        T, Vt = spectral_ridge(mdl["S"], mdl["Y"], mdl["scale"], a)
                        err[a_i] += np.sum((Rt - (St @ T) @ Vt) ** 2)
            a = ALPHAS[int(np.argmin(err))]
            mdl = ext_build(c, genes, X, train)
            St, Rt = ext_target(c, X, mdl, list(test))
            if a is None:
                sse = np.sum(Rt ** 2)
            else:
                T, Vt = spectral_ridge(mdl["S"], mdl["Y"], mdl["scale"], a)
                sse = np.sum((Rt - (St @ T) @ Vt) ** 2)
            ratios.append(sse / np.sum(Rt ** 2))
        ratios = np.array(ratios)
        pk = pkg[pkg.accession == f"OSD-{n}"]
        diff = np.abs(np.sort(ratios) - np.sort(pk.MSE_ratio.to_numpy())).max()
        record(f"m.OSD-{n}", f"OSD-{n} nested gain-2 constrained regression over 81 ground holdouts, own implementation from raw counts",
               exp, f"mean MSE ratio {ratios.mean():.4f}; improved {int((ratios < 1).sum())}/{len(ratios)}; max |sorted diff| vs package {diff:.1e}",
               at_precision(exp, ratios.mean()))


# ----------------------------------------------------------------------------- section l (new)

def section_l():
    """NEW: margin of the convex gain-2 estimator over the simplest comparator (ridge then spectral
    clipping to gain 2, same 5-candidate grid), recomputed per holdout from the package rows."""
    mb = PKG / "extensions/strengthening_20261015/method_benchmark"
    mg = pd.read_csv(mb / "matched_grid_outer_results.csv")
    so = pd.read_csv(mb / "strengthened_outer_results.csv")
    con = so[so.method == "spectral_constrained_ridge_CV"].set_index(["accession", "split_id"])
    for acc in ["OSD-561", "OSD-562", "OSD-525", "OSD-202", "OSD-352", "OSD-374", "OSD-611"]:
        clip = mg[(mg.accession == acc) & (mg.method == "matched_grid_gain2_clipped_ridge_CV")].set_index(["accession", "split_id"])
        ridge = mg[(mg.accession == acc) & (mg.method == "matched_grid_ridge_CV")].set_index(["accession", "split_id"])
        c = con.loc[clip.index]
        rc_ = c.prediction_SSE / c.baseline_SSE
        rk = clip.prediction_SSE / clip.baseline_SSE
        rr = ridge.prediction_SSE / ridge.baseline_SSE
        d = rc_ - rk
        record(f"l.{acc}", f"NEW: {acc} constrained vs matched clipped ridge vs unbounded ridge (mean MSE ratio)",
               "n/a", f"constrained {rc_.mean():.4f}; clipped {rk.mean():.4f}; unbounded ridge {rr.mean():.4f} (mean gain {ridge.operator_gain.mean():.1f}); "
               f"constrained-minus-clipped mean {d.mean():+.4f}, constrained better in {int((d < 0).sum())}/{len(d)} holdouts", "NEW")


# ----------------------------------------------------------------------------- main

SECTIONS = dict(a=section_a, b=section_b, c=section_c, d=section_d, e=section_e, f=section_f, g=section_g, h=section_h, i=section_i,
                j=section_j, k=section_k, l=section_l, m=section_m)

if __name__ == "__main__":
    wanted = sys.argv[1:] or list(SECTIONS)
    old = json.loads(RESULTS.read_text()) if RESULTS.exists() else []
    for s in wanted:
        print(f"===== section {s}", flush=True)
        SECTIONS[s]()
    keep = [r for r in old if r["id"].split(".")[0] not in wanted]
    RESULTS.write_text(json.dumps(keep + ROWS, indent=1, default=str))
    print(f"wrote {len(keep) + len(ROWS)} rows to {RESULTS}")
