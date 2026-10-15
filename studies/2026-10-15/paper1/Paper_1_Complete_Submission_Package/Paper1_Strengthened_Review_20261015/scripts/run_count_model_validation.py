#!/usr/bin/env python3
"""Regional PyDESeq2 robustness checks on actual integer counts."""
import argparse
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from run_revision import DATASETS, check_input_integrity, indices, load_dataset

HERE = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    check_input_integrity()
    output = HERE / "results/count_models"
    output.mkdir(parents=True, exist_ok=True)
    annotation = pd.read_csv(HERE / "inputs/gene_annotation.csv.gz")
    all_results, panels, diagnostics = [], [], []
    for n in DATASETS:
        d = load_dataset(n)
        counts = d["counts"].loc[d["y"].index].T.astype(np.int64)
        meta = d["meta"].copy()
        if n in ("561", "562"):
            meta["condition"] = meta.flight + "_" + meta.cohort
            label = "OSD561_cerebellum" if n == "561" else "OSD562_hippocampus"
            contrasts = [(f"{label}_{f}_ISST_vs_LAR", f"{f}_ISS_T", f"{f}_LAR") for f in ("GC", "FLT")]
            contrasts += [(f"{label}_FLT_vs_GC_within_{c}", f"FLT_{c}", f"GC_{c}") for c in ("ISS_T", "LAR")]
        else:
            meta["condition"] = meta.flight
            contrasts = [(f"OSD{n}_FLT_vs_GC", "FLT", "GC")]
            if n == "525":
                contrasts = [(f"OSD525_{a}_vs_{b}", a, b) for a, b in [("HGC", "VGC"), ("FLT", "HGC"), ("FLT", "VGC")]]
        meta["condition"] = pd.Categorical(meta.condition)
        meta["age"] = pd.Categorical(meta.age)
        design = "~ age + condition" if meta.age.nunique() > 1 else "~ condition"
        destination = output / f"OSD-{n}_gene_results.csv.gz"
        diag_path = output / f"OSD-{n}_fit.json"
        if destination.exists() and diag_path.exists() and not args.force:
            frame = pd.read_csv(destination)
            diag = json.loads(diag_path.read_text())
            print(f"OSD-{n}: cached fit", flush=True)
        else:
            print(f"OSD-{n}: fitting {counts.shape[0]} animals and {counts.shape[1]} genes; {design}", flush=True)
            dds = DeseqDataSet(counts=counts, metadata=meta[["age", "condition"]], design=design,
                              refit_cooks=True, n_cpus=1, quiet=True, low_memory=True)
            dds.deseq2()
            frames = []
            for name, positive, reference in contrasts:
                stat = DeseqStats(dds, contrast=["condition", positive, reference], alpha=.05,
                                  independent_filter=False, cooks_filter=True, quiet=True, n_cpus=1)
                stat.summary()
                result = stat.results_df.rename_axis("ensembl_id").reset_index()
                result.insert(0, "contrast", name)
                result.insert(1, "accession", f"OSD-{n}")
                frames.append(result)
            frame = pd.concat(frames, ignore_index=True)
            frame.to_csv(destination, index=False, compression={"method": "gzip", "mtime": 0})
            diag = dict(accession=f"OSD-{n}", design=design, n_animals=len(counts), n_genes=counts.shape[1],
                        condition_counts={str(k): int(v) for k, v in meta.condition.value_counts().items()},
                        n_contrasts=len(contrasts), pydeseq2_version=pydeseq2.__version__,
                        size_factor_min=float(np.min(dds.obs['size_factors'])),
                        size_factor_max=float(np.max(dds.obs['size_factors'])),
                        n_finite_dispersions=int(np.isfinite(dds.var['dispersions']).sum()),
                        fit_type=str(getattr(dds, 'fit_type', 'parametric')),
                        independent_filter=False, cooks_filter=True,
                        coefficient_shrinkage=False, python=platform.python_version())
            diag_path.write_text(json.dumps(diag, indent=2) + "\n")
            print(f"OSD-{n}: complete", flush=True)
        diagnostics.append(diag)
        all_results.append(frame)
        sets = indices(d["y"].index, annotation)
        for name, _, _ in contrasts:
            sub = frame[frame.contrast == name].set_index("ensembl_id").reindex(d["y"].index)
            for panel, idx in sets.items():
                part = sub.iloc[idx]
                panels.append(dict(contrast=name, accession=f"OSD-{n}", panel=panel, n_genes=len(idx),
                                   mean_NB_log2_effect=float(part.log2FoldChange.mean()),
                                   median_NB_log2_effect=float(part.log2FoldChange.median()),
                                   n_NB_FDR05=int((part.padj < .05).sum()),
                                   n_finite_NB_p=int(part.pvalue.notna().sum())))
    pd.concat(all_results, ignore_index=True).merge(annotation, on="ensembl_id", how="left").to_csv(
        HERE / "results/count_model_gene_results.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    panel_frame = pd.DataFrame(panels)
    linear = pd.read_csv(HERE / "results/panel_results.csv")
    linear = linear[linear.sensitivity == "design_cohort"]
    panel_frame = panel_frame.merge(linear[["contrast", "panel", "mean_log2_effect", "q_BH_all_tests_in_family"]], on=["contrast", "panel"], validate="one_to_one")
    panel_frame["same_effect_direction"] = np.sign(panel_frame.mean_NB_log2_effect) == np.sign(panel_frame.mean_log2_effect)
    panel_frame.to_csv(HERE / "results/count_model_panel_comparison.csv", index=False)
    (HERE / "results/count_model_fit_diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    print(f"Completed {len(diagnostics)} regional fits and {panel_frame.contrast.nunique()} contrasts.", flush=True)


if __name__ == "__main__":
    main()
