#!/usr/bin/env python3
"""Conservative intersection-union screen of regional flight-associated genes.

Rejecting a replication null requires every listed contrast to be supported.
Max-p tests require no independence assumption between contrast p values.
Sign discordance and missing Cook-filtered probabilities are conservative gates.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from statistics_core import bh

HERE = Path(__file__).resolve().parents[1]
CONTRASTS = {
    "hippocampus": ["OSD562_hippocampus_FLT_vs_GC_within_ISS_T", "OSD562_hippocampus_FLT_vs_GC_within_LAR", "OSD564_FLT_vs_GC"],
    "cerebellum": ["OSD561_cerebellum_FLT_vs_GC_within_ISS_T", "OSD561_cerebellum_FLT_vs_GC_within_LAR", "OSD563_FLT_vs_GC"],
    "left_cerebral_hemisphere": ["OSD613_FLT_vs_GC", "OSD612_FLT_vs_GC"],
}


def main():
    d = pd.read_csv(HERE / "results/count_model_gene_results.csv.gz")
    rows = []
    for tissue, contrasts in CONTRASTS.items():
        parts = [d[d.contrast == c].set_index("ensembl_id") for c in contrasts]
        common = parts[0].index
        for part in parts[1:]:
            common = common.intersection(part.index)
        betas = np.vstack([part.reindex(common).log2FoldChange.to_numpy() for part in parts])
        probabilities = np.vstack([part.reindex(common).pvalue.fillna(1).to_numpy() for part in parts])
        sign = np.sign(betas)
        concordant = np.isfinite(betas).all(axis=0) & ((sign == 1).all(axis=0) | (sign == -1).all(axis=0))
        p_rep = probabilities.max(axis=0)
        p_rep[~concordant] = 1
        frame = pd.DataFrame(dict(tissue=tissue, ensembl_id=common,
                                  gene_symbol=parts[0].reindex(common).gene_symbol,
                                  n_contrasts=len(parts), same_direction=concordant,
                                  replication_max_p=p_rep,
                                  minimum_absolute_log2_effect=np.min(abs(betas), axis=0),
                                  signed_minimum_log2_effect=np.sign(betas.mean(axis=0)) * np.min(abs(betas), axis=0)))
        for i, contrast in enumerate(contrasts, 1):
            frame[f"contrast_{i}"] = contrast
            frame[f"log2_effect_{i}"] = betas[i - 1]
            frame[f"pvalue_{i}"] = probabilities[i - 1]
        rows.append(frame)
    combined = pd.concat(rows, ignore_index=True)
    combined["q_BH_global_replication_family"] = bh(combined.replication_max_p)
    harmonic = np.sum(1 / np.arange(1, len(combined) + 1))
    combined["q_BY_global_replication_family"] = np.minimum(1, combined.q_BH_global_replication_family * harmonic)
    combined.to_csv(HERE / "results/regional_gene_replication.csv.gz", index=False,
                    compression={"method": "gzip", "mtime": 0})
    summary = []
    for tissue, part in combined.groupby("tissue"):
        summary.append(dict(tissue=tissue, n_genes=len(part), n_contrasts=int(part.n_contrasts.iloc[0]),
                            direction_concordant=int(part.same_direction.sum()),
                            nominal_max_p05=int((part.replication_max_p < .05).sum()),
                            BH_FDR05=int((part.q_BH_global_replication_family < .05).sum()),
                            BY_FDR05=int((part.q_BY_global_replication_family < .05).sum()),
                            minimum_replication_p=float(part.replication_max_p.min())))
    pd.DataFrame(summary).to_csv(HERE / "results/regional_gene_replication_summary.csv", index=False)
    combined.sort_values("replication_max_p").groupby("tissue").head(20).to_csv(HERE / "results/regional_gene_replication_ranked_candidates.csv", index=False)
    (HERE / "results/regional_gene_replication_design.json").write_text(json.dumps(dict(
        contrasts=CONTRASTS, number_of_tests=len(combined), family="All region-gene replication tests together",
        inference="max of valid two-sided component Wald p values; discordant signs or missing values are conservative gates",
        dependence="max-p validity does not require independence of component contrasts; BH standard dependence assumptions, BY arbitrary dependence",
        limitations="Exploratory post-audit screen; shared estimated dispersion/age parameters are not independent biological replication. Source cohorts and collection protocols differ. No undiscovered validated gene is claimed from ranking alone."), indent=2) + "\n")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
