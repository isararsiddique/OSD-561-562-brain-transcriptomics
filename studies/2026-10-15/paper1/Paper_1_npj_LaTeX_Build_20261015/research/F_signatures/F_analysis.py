"""Paper 1 / F: published-signature test of the protocol-versus-flight finding.

Implements F_PLAN.md exactly. Refuses to run unless F_PLAN.md and signatures.json match the hashes recorded
in F_PLAN.sha256, and refuses to compute any signature number unless the rPRG reproduction gate passes.
Reuses the package's loaders and statistics (read-only import) and the D_protocol estimator.
"""
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
D_TABLES = HERE.parent / "D_protocol" / "D_tables"
sys.path.insert(0, str(PKG / "scripts"))
import run_revision as rr  # noqa: E402
from statistics_core import Contrast, animal_permutations, permutation_p, bh  # noqa: E402

CONFIRMATORY = ["S1_Hrvatin2018_ERG_up", "S3_Marsh2022_exAM_allCNS", "S4_Marsh2022_exAM_myeloid",
                "S5_vandenBrink2017_dissociation", "S6_Wu2017_IEG"]
CONTROL = "S2_Hrvatin2018_LRG_only"
MIN_GENES = 5
BASE_PANELS = ["rPRG", "dPRG", "SRG", "HSP", "Glia_Complement", "BBB_Endothelium", "Circadian"]


def check_hashes():
    recorded = {}
    for line in (HERE / "F_PLAN.sha256").read_text().splitlines():
        if line and not line.startswith("#"):
            digest, name = line.split(maxsplit=1)
            recorded[name.strip()] = digest
    for name in ("F_PLAN.md", "signatures.json"):
        now = hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        if recorded.get(name) != now:
            raise SystemExit(f"{name} changed after the plan was hashed; refusing to run.")
    return recorded


def exact(data, mask, field, positive, idx_sets):
    """Exact whole-animal permutation test for several gene-index sets on one contrast."""
    meta = data["meta"].loc[mask].copy()
    y = data["y"].loc[:, meta.index].to_numpy()
    z = rr.nuisance(meta)
    g = (meta[field] == positive).to_numpy(int)
    blocks = ("age",) if meta.age.nunique() > 1 else ()
    labels, method, total = animal_permutations(meta, g, blocks)
    out = {}
    for key, idx in idx_sets.items():
        if len(idx) < MIN_GENES:
            out[key] = None
            continue
        c = Contrast(y[idx], z)
        obs = float(c.fit(g)[0].mean())
        null = c.permutation_betas(labels).mean(axis=0)
        out[key] = {"effect": obs, "p_exact": float(permutation_p(obs, null, method)),
                    "n_assignments": int(labels.shape[1]), "method": method,
                    "n_animals": int(meta.animal.nunique()),
                    "n_pos_animals": int(meta.loc[g == 1].animal.nunique())}
    return out


def main():
    hashes = check_hashes()
    sigs = json.loads((HERE / "signatures.json").read_text())["signatures"]
    annotation = pd.read_csv(PKG / "inputs/gene_annotation.csv.gz")
    e2s = annotation.set_index("ensembl_id").gene_symbol.to_dict()
    panel_syms = set(pd.read_csv(PKG / "results/panel_definitions.csv")
                     .query("panel in @BASE_PANELS").symbol)

    data = {n: rr.load_dataset(n) for n in ("561", "562", "563", "564", "612", "525")}
    pool = rr.pooled(data["561"], data["562"])

    # RRRM-2 analyses and contrasts (D_protocol definitions)
    def dom_view(d):
        m = d["meta"]
        return {"y": d["y"], "meta": m.assign(dom=np.where(m.flight == "GC", "GC_ISST", "FLT_LAR"))}

    def is_dom(m):
        return ((m.flight == "GC") & (m.cohort == "ISS_T")) | ((m.flight == "FLT") & (m.cohort == "LAR"))

    rrrm2 = {"cerebellum": data["561"], "hippocampus": data["562"], "pooled": pool}
    specs = {
        "P_G": (lambda m: m.flight == "GC", "cohort", "ISS_T", False),
        "P_F": (lambda m: m.flight == "FLT", "cohort", "ISS_T", False),
        "F_L": (lambda m: m.cohort == "LAR", "flight", "FLT", False),
        "F_I": (lambda m: m.cohort == "ISS_T", "flight", "FLT", False),
        "D": (is_dom, "dom", "GC_ISST", True),
    }
    plan = []
    for region, d in rrrm2.items():
        for key, (sel, field, pos, needs_dom) in specs.items():
            if region == "pooled" and key in ("F_L", "F_I"):
                continue
            plan.append((f"RRRM-2 {region}", region, key, d, sel, field, pos, needs_dom))
    cross = [("RR-10 cerebellum", "563", "FLT_vs_GC", None, "FLT"), ("RR-10 hippocampus", "564", "FLT_vs_GC", None, "FLT"),
             ("RR-10 left hemisphere", "612", "FLT_vs_GC", None, "FLT"),
             ("RR-23 cerebellum", "525", "FLT_vs_HGC", ("FLT", "HGC"), "FLT"),
             ("RR-23 cerebellum", "525", "FLT_vs_VGC", ("FLT", "VGC"), "FLT"),
             ("RR-23 cerebellum", "525", "HGC_vs_VGC", ("HGC", "VGC"), "HGC")]
    for label, n, key, pair, pos in cross:
        d = data[n]
        sel = (lambda m: np.ones(len(m), bool)) if pair is None else (lambda m, pr=pair: m.flight.isin(pr))
        plan.append((label, f"OSD-{n}", key, d, sel, "flight", pos, False))

    # ---------------------------------------------------------------- reproduction gate
    gate_rows = []
    d1 = pd.read_csv(D_TABLES / "D1_reproduction_vs_package.csv").query("panel == 'rPRG'")
    d2 = pd.read_csv(D_TABLES / "D2_dominance_GC_ISST_vs_FLT_LAR.csv").query("panel == 'rPRG'")
    pkg = pd.read_csv(PKG / "results/panel_results.csv").query("sensitivity == 'design_cohort' and panel == 'rPRG'")
    pkg_names = {"563": "OSD563_FLT_vs_GC", "564": "OSD564_FLT_vs_GC", "612": "OSD612_FLT_vs_GC"}
    for label, analysis, key, d, sel, field, pos, needs_dom in plan:
        dd = dom_view(d) if needs_dom else d
        idx = rr.indices(d["y"].index, annotation)["rPRG"]
        r = exact(dd, sel(dd["meta"]), field, pos, {"rPRG": idx})["rPRG"]
        if label.startswith("RRRM-2"):
            if key == "D":
                ref = d2[d2.region == analysis].iloc[0]
                exp_e, exp_p = ref.GC_ISST_minus_FLT_LAR, ref.p_exact
            else:
                ref = d1[(d1.region == analysis) & (d1.contrast == key)].iloc[0]
                exp_e, exp_p = ref.recomputed, ref.p_recomputed
        else:
            name = pkg_names.get(analysis[-3:], f"OSD525_{key}")
            ref = pkg[pkg.contrast == name].iloc[0]
            exp_e, exp_p = ref.mean_log2_effect, ref.p_animal_permutation
        gate_rows.append({"analysis": label, "contrast": key, "effect": r["effect"], "expected": float(exp_e),
                          "p": r["p_exact"], "p_expected": float(exp_p),
                          "abs_diff": abs(r["effect"] - exp_e), "abs_pdiff": abs(r["p_exact"] - exp_p)})
    gate = pd.DataFrame(gate_rows)
    (HERE / "F_tables").mkdir(exist_ok=True)
    gate.to_csv(HERE / "F_tables/F0_reproduction_gate.csv", index=False)
    print(f"gate: {len(gate)} rPRG rows; max |effect diff| {gate.abs_diff.max():.1e}; max |p diff| {gate.abs_pdiff.max():.1e}")
    if gate.abs_diff.max() > 1e-12 or gate.abs_pdiff.max() > 1e-12:
        raise SystemExit("Reproduction gate failed; no signature number computed.")

    # ---------------------------------------------------------------- signature index sets
    def index_sets(d):
        genes = d["y"].index
        syms = np.array([e2s.get(g, "") for g in genes])
        out, mapping = {}, {}
        for sid, s in sigs.items():
            if s["map_by"] == "ensembl":
                want = {g["ensembl"] for g in s["genes"]}
                hit = np.isin(genes, list(want))
                src_sym = {g["ensembl"]: g["symbol"] for g in s["genes"]}
                in_panel = np.array([(sym in panel_syms) or (src_sym.get(e, "") in panel_syms)
                                     for e, sym in zip(genes, syms)])
                missing = sorted(src_sym[e] for e in want - set(genes[hit]))
            else:
                want = {g["symbol"] for g in s["genes"]}
                hit = np.isin(syms, list(want))
                in_panel = np.isin(syms, list(panel_syms))
                missing = sorted(want - set(syms[hit]))
            out[(sid, "full")] = np.flatnonzero(hit)
            out[(sid, "novel")] = np.flatnonzero(hit & ~in_panel)
            mapping[sid] = {"n_source": len(s["genes"]), "n_mapped": int(hit.sum()),
                            "n_mapped_novel": int((hit & ~in_panel).sum()), "unmapped": missing}
        return out, mapping

    rows, mapping_report = [], {}
    for label, analysis, key, d, sel, field, pos, needs_dom in plan:
        dd = dom_view(d) if needs_dom else d
        sets, mapping = index_sets(d)
        mapping_report[label] = mapping
        res = exact(dd, sel(dd["meta"]), field, pos, sets)
        for (sid, ver), r in res.items():
            rows.append({"signature": sid, "version": ver, "analysis": label, "contrast": key,
                         "n_genes": int(len(sets[(sid, ver)])), "testable": r is not None,
                         **(r or {})})
        print(f"done {label} {key}", flush=True)
    res = pd.DataFrame(rows)

    # ---------------------------------------------------------------- families and decision rules
    def family(contrast):
        f = res[(res.analysis == "RRRM-2 pooled") & (res.contrast == contrast) &
                res.signature.isin(CONFIRMATORY) & res.testable].copy()
        f["q_BH"] = bh(f.p_exact.to_numpy())
        return f

    prim, sec = family("D"), family("P_G")
    res = res.merge(prim[["signature", "version", "q_BH"]].rename(columns={"q_BH": "q_primary_D"}),
                    on=["signature", "version"], how="left")
    res.loc[~((res.analysis == "RRRM-2 pooled") & (res.contrast == "D")), "q_primary_D"] = np.nan
    res = res.merge(sec[["signature", "version", "q_BH"]].rename(columns={"q_BH": "q_secondary_PG"}),
                    on=["signature", "version"], how="left")
    res.loc[~((res.analysis == "RRRM-2 pooled") & (res.contrast == "P_G")), "q_secondary_PG"] = np.nan
    res.to_csv(HERE / "F_tables/F1_signature_results.csv", index=False)

    ok = lambda f, ver: f[(f.version == ver) & (f.effect > 0) & (f.q_BH <= 0.05)]
    r1 = len(ok(prim, "full"))
    novel_testable = int((prim.version == "novel").sum())
    r2 = len(ok(prim, "novel"))
    pg = res[(res.analysis == "RRRM-2 pooled") & (res.contrast == "P_G") & (res.version == "full")].set_index("signature")
    r3 = bool(pg.loc[CONTROL, "p_exact"] > 0.05 and pg.loc[CONTROL, "effect"] < pg.loc["S1_Hrvatin2018_ERG_up", "effect"])
    decisions = {
        "R1_published_signatures_reproduce_dominance": {"n_pass": r1, "of": 5, "needed": 3, "verdict": r1 >= 3},
        "R2_beyond_package_panel_genes": {"n_pass": r2, "of_testable": novel_testable, "needed": 3, "verdict": r2 >= 3},
        "R3_time_scale_specificity": {"S2_PG": float(pg.loc[CONTROL, "effect"]), "S2_p": float(pg.loc[CONTROL, "p_exact"]),
                                      "S1_PG": float(pg.loc["S1_Hrvatin2018_ERG_up", "effect"]), "verdict": r3},
    }
    out = {"plan_hashes": hashes, "gate": {"rows": len(gate), "max_abs_effect_diff": float(gate.abs_diff.max()),
                                            "max_abs_p_diff": float(gate.abs_pdiff.max())},
           "mapping": mapping_report, "primary_family": prim.to_dict(orient="records"),
           "secondary_family": sec.to_dict(orient="records"), "decisions": decisions}
    (HERE / "F_results.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(decisions, indent=1, default=float))


if __name__ == "__main__":
    main()
