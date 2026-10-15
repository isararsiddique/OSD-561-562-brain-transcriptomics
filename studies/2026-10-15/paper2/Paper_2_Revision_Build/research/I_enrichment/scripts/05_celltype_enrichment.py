#!/usr/bin/env python3
"""Cerebellar cell-type over-representation against the MEASURED gene universe.

Marker source (downloaded at source, not recalled):
  Kozareva V., Martin C., Osorno T., Rudolph S., Pan C., Regehr W.G.,
  Macosko E.Z. et al. "A transcriptomic atlas of mouse cerebellar cortex
  comprehensively defines cell types." Nature 598, 214-219 (2021).
  DOI 10.1038/s41586-021-03220-z
  Supplementary Table 2 ("Highly differentially expressed genes across all
  clusters"), sheet Supplementary_Table2 of Supplementary Tables 1-4
  (41586_2021_3220_MOESM4_ESM.xlsx).

Marker definition: within each of the atlas's 46 annotated subtypes, genes with
padj <= 0.05 and logFC > 0, ranked by logFC, top TOP_PER_SUBTYPE retained; the
subtypes are then collapsed to major cerebellar classes. A fixed per-subtype cap
is used because the published table is itself truncated at 500 rows for many
subtypes, so an unbounded "all significant markers" rule would hand bigger
marker sets to the subtypes that happen not to be truncated.

Test: one-sided Fisher exact (over-representation) of each class's markers in
each query gene set, with BOTH marker set and query set intersected with the
measured background first, so the universe is the genes this junction assay
could see -- never the genome. BH correction across classes within each gene set.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parents[1]
TAB = HERE / "tables"
MRK = HERE / "markers"

TOP_PER_SUBTYPE = 200

SOURCE = {
    "citation": (
        "Kozareva V. et al. A transcriptomic atlas of mouse cerebellar cortex "
        "comprehensively defines cell types. Nature 598, 214-219 (2021)."
    ),
    "doi": "10.1038/s41586-021-03220-z",
    "table": "Supplementary Table 2 (sheet Supplementary_Table2), "
             "Highly differentially expressed genes across all clusters",
    "file": "41586_2021_3220_MOESM4_ESM.xlsx",
    "url": "https://media.springernature.com/original/springer-static/esm/"
           "art%3A10.1038%2Fs41586-021-03220-z/MediaObjects/"
           "41586_2021_3220_MOESM4_ESM.xlsx",
}

# Subtype -> major class. Names are the atlas's own cluster labels.
def to_class(subtype: str) -> str:
    s = subtype
    if s.startswith("Purkinje"):
        return "Purkinje"
    if s.startswith("Granule"):
        return "Granule"
    if s.startswith("MLI1"):
        return "MLI1"
    if s.startswith("MLI2"):
        return "MLI2"
    if s.startswith("Golgi"):
        return "Golgi"
    if s.startswith("PLI"):
        return "PLI"
    if s.startswith("UBC"):
        return "Unipolar brush"
    if s.startswith("Bergmann"):
        return "Bergmann glia"
    if s.startswith("Astrocyte"):
        return "Astrocyte"
    if s.startswith("ODC"):
        return "Oligodendrocyte"
    if s.startswith("OPC"):
        return "OPC"
    if s == "Microglia":
        return "Microglia"
    if s == "Macrophage":
        return "Macrophage"
    if s.startswith("Endothelial"):
        return "Endothelial/mural"
    if s == "Fibroblast":
        return "Fibroblast"
    if s == "Choroid":
        return "Choroid plexus"
    if s == "Ependymal":
        return "Ependymal"
    raise ValueError(f"unmapped subtype {subtype!r}")


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m, dtype=float)
    q[order] = np.minimum(ranked, 1.0)
    return q


def read_list(name: str) -> list[str]:
    return [g for g in (TAB / name).read_text().split("\n") if g.strip()]


def main() -> None:
    t2 = pd.read_csv(MRK / "Kozareva2021_Supplementary_Table2.csv")
    t2 = t2.dropna(subset=["gene", "subtype"])
    sig = t2[(t2.padj <= 0.05) & (t2.logFC > 0)].copy()
    sig["cell_class"] = sig.subtype.map(to_class)
    sig = sig.sort_values(["subtype", "logFC"], ascending=[True, False])
    top = sig.groupby("subtype", group_keys=False).head(TOP_PER_SUBTYPE)

    background = set(read_list("background_measured_genes.txt"))

    # Marker sets per class, before and after restriction to the measured universe.
    markers_all = {c: set(d.gene) for c, d in top.groupby("cell_class")}
    markers = {c: sorted(g & background) for c, g in markers_all.items()}

    # Class-unique markers (sensitivity): drop genes claimed by >1 class.
    counts: dict[str, int] = {}
    for g in markers_all.values():
        for x in g:
            counts[x] = counts.get(x, 0) + 1
    markers_unique = {
        c: sorted({x for x in g if counts[x] == 1} & background)
        for c, g in markers_all.items()
    }

    # Export the marker definitions actually used.
    rows = []
    for c in sorted(markers_all):
        for g in sorted(markers_all[c]):
            rows.append({
                "cell_class": c,
                "gene": g,
                "in_measured_background": g in background,
                "class_unique": counts[g] == 1,
            })
    pd.DataFrame(rows).to_csv(TAB / "cerebellar_marker_sets_used.csv", index=False)

    queries = {
        "FLT_GC_BY": read_list("set_FLT_GC_BY.txt"),
        "FLT_GC_top100": read_list("set_FLT_GC_top100.txt"),
        "GC_VIV_top100": read_list("set_GC_VIV_top100.txt"),
    }

    out = []
    for defn, mset in (("all_markers", markers), ("class_unique_markers", markers_unique)):
        for qlabel, qgenes in queries.items():
            q = set(qgenes) & background
            recs = []
            for cls in sorted(mset):
                M = set(mset[cls])
                if not M:
                    continue
                a = len(q & M)                       # query & marker
                b = len(q) - a                       # query, not marker
                c_ = len(M) - a                      # marker, not query
                d = len(background) - len(q) - c_    # neither
                orr, p = fisher_exact([[a, b], [c_, d]], alternative="greater")
                # Haldane-Anscombe OR so a zero cell still gives a finite estimate.
                orr_ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c_ + 0.5))
                recs.append({
                    "marker_definition": defn,
                    "gene_set": qlabel,
                    "cell_class": cls,
                    "n_background": len(background),
                    "n_query_in_background": len(q),
                    "n_markers_in_background": len(M),
                    "overlap": a,
                    "expected_overlap": len(q) * len(M) / len(background),
                    "odds_ratio_fisher": orr,
                    "odds_ratio_haldane": orr_ha,
                    "p_fisher_one_sided": p,
                    "overlapping_genes": ";".join(sorted(q & M)),
                })
            df = pd.DataFrame(recs)
            df["q_BH_within_set"] = bh(df.p_fisher_one_sided.to_numpy())
            df["significant_BH_05"] = df.q_BH_within_set <= 0.05
            out.append(df)

    res = pd.concat(out, ignore_index=True).sort_values(
        ["marker_definition", "gene_set", "p_fisher_one_sided"]
    )
    res.to_csv(TAB / "celltype_enrichment_results.csv", index=False)

    meta = {
        "marker_source": SOURCE,
        "marker_rule": {
            "filter": "padj <= 0.05 and logFC > 0 in the published table",
            "top_per_subtype_by_logFC": TOP_PER_SUBTYPE,
            "n_subtypes": int(top.subtype.nunique()),
            "n_classes": int(top.cell_class.nunique()),
        },
        "background_n": len(background),
        "marker_sizes_in_background": {c: len(g) for c, g in sorted(markers.items())},
        "class_unique_marker_sizes_in_background": {
            c: len(g) for c, g in sorted(markers_unique.items())
        },
        "test": "one-sided Fisher exact (over-representation), BH across classes "
                "within each gene set; universe = measured genes",
        "n_significant_BH_05": int(res.significant_BH_05.sum()),
    }
    (TAB / "celltype_enrichment_meta.json").write_text(json.dumps(meta, indent=2))

    print(json.dumps(meta["marker_sizes_in_background"], indent=1))
    for (defn, gs), sub in res.groupby(["marker_definition", "gene_set"]):
        print(f"\n== {defn} / {gs}  (min q={sub.q_BH_within_set.min():.3g}, "
              f"n sig={int(sub.significant_BH_05.sum())})")
        print(sub.head(5)[[
            "cell_class", "overlap", "n_markers_in_background", "expected_overlap",
            "odds_ratio_haldane", "p_fisher_one_sided", "q_BH_within_set"
        ]].to_string(index=False))


if __name__ == "__main__":
    main()
