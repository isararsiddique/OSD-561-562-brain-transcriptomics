"""Descriptive, outcome-unselected ordination and single-animal influence."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]


def main():
    results = BASE / "results"
    meta = pd.read_csv(BASE / "inputs/animal_metadata.csv").set_index("sample")
    gates = pd.read_csv(results / "measurement_gates.csv.gz").set_index("event_id")
    fractions = pd.read_csv(results / "all26_animal_fractions.csv.gz", index_col=0)
    x = fractions.loc[gates.measured_primary_all26].T
    centered = x.to_numpy() - x.to_numpy().mean(axis=0)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    explained = singular**2 / np.sum(singular**2)
    scores = pd.DataFrame(u[:, :2] * singular[:2], index=x.index, columns=["PC1", "PC2"])
    scores["group"] = meta.loc[x.index, "group"]
    scores.to_csv(results / "unsupervised_animal_pca.csv")
    rows = []
    for ga, gb in [("FLT", "GC"), ("GC", "VIV"), ("FLT", "VIV")]:
        a, b = x.loc[meta.group == ga].to_numpy(), x.loc[meta.group == gb].to_numpy()
        full = 100 * np.sqrt(np.mean((a.mean(axis=0)-b.mean(axis=0))**2))
        for group, samples in [(ga, x.index[meta.loc[x.index, "group"] == ga]), (gb, x.index[meta.loc[x.index, "group"] == gb])]:
            for omitted in samples:
                reduced = x.drop(index=omitted)
                diff = reduced.loc[meta.loc[reduced.index, "group"] == ga].mean().to_numpy()-reduced.loc[meta.loc[reduced.index, "group"] == gb].mean().to_numpy()
                rms = 100 * np.sqrt(np.mean(diff**2))
                rows.append(dict(contrast=ga+"_minus_"+gb, omitted_animal=omitted, omitted_group=group, RMS_contrast_pp=rms, full_RMS_pp=full))
    loo = pd.DataFrame(rows)
    assert len(loo) == 52 and np.isfinite(loo.RMS_contrast_pp).all()
    loo.to_csv(results / "leave_one_animal_out_descriptive.csv", index=False)
    (results / "descriptive_summary.json").write_text(json.dumps(dict(PCA_PC1_variance_fraction=float(explained[0]), PCA_PC2_variance_fraction=float(explained[1]), PCA_feature_count=x.shape[1], PCA_animal_count=x.shape[0], leave_one_out_rows=len(loo), interpretation="PCA uses the complete primary measured panel without event significance selection. Leave-one-out ranges are descriptive and are not confidence intervals or new hypothesis tests."), indent=2))
    print(loo.groupby("contrast").RMS_contrast_pp.agg(["min", "max"]).to_string())


if __name__ == "__main__":
    main()
