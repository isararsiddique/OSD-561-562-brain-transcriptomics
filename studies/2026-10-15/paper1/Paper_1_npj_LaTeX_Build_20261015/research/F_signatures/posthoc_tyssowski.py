"""POST HOC (not in F_PLAN.md): Tyssowski et al. 2018's own published gene classes.

Run after the registered test failed R1. Uses the identical estimator and exact tests. Nothing here can change
the registered verdict in F_results.json; results go to posthoc_results.json and F_tables/PH1_*.csv only.

Question: the package rPRG panel is a 20-gene selection (15 shared with Tyssowski's rapid class). Do the
originating paper's complete classes (19 rapid PRGs, 116 delayed PRGs, 38 SRGs) show the protocol pattern? Under
Tyssowski's model, brief activity induces only the rapid wave, so a brief perimortem stimulus should move rapid
PRGs, much less delayed PRGs, and not SRGs.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import F_analysis as fa  # noqa: E402  (reuses rr, exact(), PKG; main() is not called)
from statistics_core import bh  # noqa: E402

rr = fa.rr


def main():
    cls = json.loads((HERE / "posthoc_tyssowski_classes.json").read_text())["classes"]
    annotation = pd.read_csv(fa.PKG / "inputs/gene_annotation.csv.gz")
    e2s = annotation.set_index("ensembl_id").gene_symbol.to_dict()
    data = {n: rr.load_dataset(n) for n in ("561", "562")}
    regions = {"cerebellum": data["561"], "hippocampus": data["562"], "pooled": rr.pooled(data["561"], data["562"])}

    def dom_view(d):
        m = d["meta"]
        return {"y": d["y"], "meta": m.assign(dom=np.where(m.flight == "GC", "GC_ISST", "FLT_LAR"))}

    specs = {
        "P_G": (lambda m: m.flight == "GC", "cohort", "ISS_T", False),
        "F_L": (lambda m: m.cohort == "LAR", "flight", "FLT", False),
        "D": (lambda m: ((m.flight == "GC") & (m.cohort == "ISS_T")) | ((m.flight == "FLT") & (m.cohort == "LAR")),
              "dom", "GC_ISST", True),
    }
    rows, mapped = [], {}
    for region, d in regions.items():
        syms = np.array([e2s.get(g, "") for g in d["y"].index])
        sets = {k: np.flatnonzero(np.isin(syms, v)) for k, v in cls.items()}
        sets["package_rPRG"] = rr.indices(d["y"].index, annotation)["rPRG"]
        mapped[region] = {k: int(len(v)) for k, v in sets.items()}
        for key, (sel, field, pos, needs_dom) in specs.items():
            if region == "pooled" and key == "F_L":
                continue
            dd = dom_view(d) if needs_dom else d
            for name, r in fa.exact(dd, sel(dd["meta"]), field, pos, sets).items():
                rows.append({"analysis": region, "contrast": key, "gene_class": name,
                             "n_genes": mapped[region][name], **(r or {})})
    res = pd.DataFrame(rows)
    fam = res[(res.analysis == "pooled") & (res.contrast == "D") & (res.gene_class != "package_rPRG")].copy()
    fam["q_BH_posthoc3"] = bh(fam.p_exact.to_numpy())
    res = res.merge(fam[["analysis", "contrast", "gene_class", "q_BH_posthoc3"]],
                    on=["analysis", "contrast", "gene_class"], how="left")
    res.to_csv(HERE / "F_tables/PH1_tyssowski_classes.csv", index=False)
    (HERE / "posthoc_results.json").write_text(json.dumps(
        {"label": "POST HOC - cannot change the registered verdict in F_results.json",
         "mapped_genes": mapped, "rows": res.to_dict(orient="records")}, indent=1, default=float))
    show = res.pivot_table(index=["gene_class", "analysis"], columns="contrast", values=["effect", "p_exact"])
    print(show.round(4).to_string())
    print(fam[["gene_class", "effect", "p_exact", "q_BH_posthoc3"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
