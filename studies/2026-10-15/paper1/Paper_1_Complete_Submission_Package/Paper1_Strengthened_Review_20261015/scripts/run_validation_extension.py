#!/usr/bin/env python3
"""Known-truth simulations, operator geometry and joint animal bootstraps."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from run_revision import DATASETS, check_input_integrity, indices, load_dataset, nuisance
from statistics_core import Contrast, bh, bounded_gamma, learn_anchor, transform

HERE = Path(__file__).resolve().parents[1]


def spectral_gain(w, gamma, strength=1.0):
    """Exact largest singular value of I - strength*gamma*w.T, with unit w."""
    a = float(w @ gamma)
    h = float(np.linalg.norm(gamma - a * w))
    # Identity on the remaining G-2 dimensions, this matrix on span(w, gamma).
    block = np.array([[1 - strength * a, 0.0], [-strength * h, 1.0]])
    return max(1.0, float(np.linalg.svd(block, compute_uv=False)[0]))


def benchmark(repeats):
    seed, genes, n_train, n_test = 20261014, 1200, 24, 24
    setup = np.random.default_rng(seed)
    sensor = np.arange(31)
    outside = np.arange(31, genes)
    u = np.zeros(genes)
    u[sensor] = 1.0
    u[31:151] = setup.choice([-1., 1.], 120) * .7
    beta_disjoint = np.zeros(genes)
    beta_disjoint[200:260] = setup.choice([-1., 1.], 60) * .8
    noise_sd = setup.uniform(.65, 1.15, genes)
    loadings = setup.normal(0, .12, (genes, 5))
    group_train = np.r_[np.zeros(n_train // 2), np.ones(n_train // 2)]
    group_test = np.r_[np.zeros(n_test // 2), np.ones(n_test // 2)]
    oracle_w = np.zeros(genes)
    oracle_w[sensor] = 1 / np.sqrt(len(sensor))
    oracle = (u / np.linalg.norm(u), oracle_w, u / (oracle_w @ u))
    scenarios = ["pure_null", "confounded_no_biology", "orthogonal_biology",
                 "overlapping_biology", "changed_coupling"]
    methods = ["raw", "projection", "score", "bounded_score", "half_score", "oracle_projection", "oracle_score"]
    result = []
    for sindex, scenario in enumerate(scenarios):
        biology = np.zeros(genes) if scenario in ("pure_null", "confounded_no_biology") else beta_disjoint.copy()
        if scenario == "overlapping_biology":
            biology += .5 * u
        target_u = u.copy()
        if scenario == "changed_coupling":
            target_u[31:151] *= -1
        artifact_strength = 0.0 if scenario == "pure_null" else 1.0
        true_positive = np.abs(biology) > 1e-12
        for repeat in range(repeats):
            rng = np.random.default_rng(np.random.SeedSequence([seed, sindex, repeat]))
            train_score = 1.5 * group_train + rng.normal(0, .4, n_train)
            test_score = artifact_strength * group_test + rng.normal(0, .4, n_test)
            y_train = 10 + np.outer(u, train_score) + noise_sd[:, None] * rng.normal(size=(genes, n_train))
            y_train += loadings @ rng.normal(size=(5, n_train))
            y_test = 10 + np.outer(target_u, test_score) + np.outer(biology, group_test)
            y_test += noise_sd[:, None] * rng.normal(size=(genes, n_test))
            y_test += loadings @ rng.normal(size=(5, n_test))
            learned = learn_anchor(y_train, group_train, None, sensor)
            raw_beta, _ = Contrast(y_test).fit(group_test)
            for method in methods:
                anchor = oracle if method.startswith("oracle") else learned
                operator = method.replace("oracle_", "")
                if operator == "half_score":
                    corrected = y_test - .5 * np.outer(anchor[2], anchor[1] @ y_test)
                    bio_corrected = biology - .5 * anchor[2] * (anchor[1] @ biology)
                else:
                    corrected = transform(y_test, anchor, operator)
                    bio_corrected = transform(biology[:, None], anchor, operator)[:, 0]
                fit = Contrast(corrected)
                estimate, tvalue = fit.fit(group_test)
                q = bh(2 * stats.t.sf(abs(tvalue), fit.df))
                called = q < .05
                discoveries = int(called.sum())
                false = int((called & ~true_positive).sum())
                power = float((called & true_positive).sum() / true_positive.sum()) if true_positive.any() else np.nan
                retention = float(np.sum(bio_corrected**2) / np.sum(biology**2)) if true_positive.any() else np.nan
                gamma_gain = bounded_gamma(anchor[1], anchor[2]) if operator == "bounded_score" else anchor[2]
                gain = 1. if operator in ("raw", "projection") else spectral_gain(anchor[1], gamma_gain, .5 if operator == "half_score" else 1.)
                result.append(dict(scenario=scenario, repeat=repeat, method=method, n_train=n_train, n_test=n_test,
                                   n_shared_samples=0, n_genes=genes, n_true_positive=int(true_positive.sum()),
                                   biological_effect_RMSE=float(np.sqrt(np.mean((estimate - biology)**2))),
                                   outside_sensor_effect_RMSE=float(np.sqrt(np.mean((estimate[outside] - biology[outside])**2))),
                                   biological_false_discovery_proportion=false / max(discoveries, 1),
                                   power=power, discoveries=discoveries,
                                   true_biology_energy_retained_fraction=retention,
                                   sensor_effect_energy_removed_pct=float(100 * (1 - np.sum(estimate[sensor]**2) / np.sum(raw_beta[sensor]**2))),
                                   spectral_gain=gain))
        print(f"Simulation: {scenario} / {repeats} independent repetitions complete", flush=True)
    frame = pd.DataFrame(result)
    frame.to_csv(HERE / "results/known_truth_benchmark_repetitions.csv", index=False)
    metrics = ["biological_effect_RMSE", "outside_sensor_effect_RMSE", "biological_false_discovery_proportion",
               "power", "discoveries", "true_biology_energy_retained_fraction", "sensor_effect_energy_removed_pct", "spectral_gain"]
    summary = []
    for (scenario, method), sub in frame.groupby(["scenario", "method"], sort=False):
        row = dict(scenario=scenario, method=method, repetitions=len(sub))
        for metric in metrics:
            values = sub[metric].dropna().to_numpy()
            row[metric + "_mean"] = float(values.mean()) if len(values) else np.nan
            row[metric + "_mc_se"] = float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else np.nan
            row[metric + "_lower95"] = float(np.quantile(values, .025)) if len(values) else np.nan
            row[metric + "_upper95"] = float(np.quantile(values, .975)) if len(values) else np.nan
        summary.append(row)
    pd.DataFrame(summary).to_csv(HERE / "results/known_truth_benchmark_summary.csv", index=False)
    config = dict(seed=seed, scenarios=scenarios, methods=methods, repetitions=repeats, genes=genes,
                  ground_samples=n_train, target_samples=n_test, sensor_genes=31,
                  nuisance_support_genes=151, disjoint_biology_genes=60, biological_effect_disjoint=.8,
                  overlapping_biology_nuisance_multiplier=.5, ground_cohort_nuisance_difference=1.5,
                  target_cohort_nuisance_difference=1.0, latent_nuisance_sd=.4,
                  noise_sd_range=[.65, 1.15], additional_latent_noise_factors=5,
                  gene_loading_sd=.12, oracle_description="Actual simulation nuisance; explanatory upper-information reference",
                  probability_interpretation="Biological false calls compare conditional OLS/BH calls with simulated biological effects, not with the full confounded group mean.",
                  limitation="Gaussian mechanism stress tests are not external experimental validation or a realistic complete RNA-seq simulator.")
    (HERE / "results/known_truth_benchmark_config.json").write_text(json.dumps(config, indent=2) + "\n")


def normalized_fixed_universe(counts):
    positive = (counts > 0).all(axis=1)
    log = np.log(counts[positive])
    factors = np.exp(np.median(log - log.mean(axis=1, keepdims=True), axis=0))
    return np.log2(counts / factors + 1)


def joint_bootstrap(datasets, annotation, repetitions):
    rows, summaries = [], []
    for n in ("561", "562", "525"):
        d = datasets[n]
        counts = d["counts"].loc[d["y"].index].to_numpy(float)
        meta = d["meta"].copy()
        sets = indices(d["y"].index, annotation)
        sensor = sets["Sensor31"]
        outside = np.setdiff1d(np.arange(len(counts)), sensor)
        strata = [np.array(pos) for pos in meta.reset_index(drop=True).groupby(["flight", "cohort", "age"], sort=True).indices.values()]
        rng = np.random.default_rng(20261014 + int(n))
        observed = {}
        for repetition in range(-1, repetitions):
            draw = np.arange(len(meta)) if repetition == -1 else np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
            m = meta.iloc[draw].reset_index(drop=True)
            y = normalized_fixed_universe(counts[:, draw])
            if n == "525":
                train = m.flight.isin(["HGC", "VGC"]).to_numpy()
                test = m.flight.isin(["FLT", "HGC"]).to_numpy()
                gt = (m.flight[train] == "HGC").to_numpy(int)
                ge = (m.flight[test] == "FLT").to_numpy(int)
            else:
                train = (m.flight == "GC").to_numpy()
                test = (m.flight == "FLT").to_numpy()
                gt = (m.cohort[train] == "ISS_T").to_numpy(int)
                ge = (m.cohort[test] == "ISS_T").to_numpy(int)
            anchor = learn_anchor(y[:, train], gt, nuisance(m.loc[train]), sensor)
            raw_beta, _ = Contrast(y[:, test], nuisance(m.loc[test])).fit(ge)
            for mode in ("raw", "projection", "score", "bounded_score"):
                effect, _ = Contrast(transform(y[:, test], anchor, mode), nuisance(m.loc[test])).fit(ge)
                common = dict(accession=f"OSD-{n}", repetition=repetition, mode=mode)
                values = {panel: float(effect[idx].mean()) for panel, idx in sets.items()}
                values["sensor_energy_removed_pct"] = float(100 * (1 - np.sum(effect[sensor]**2) / np.sum(raw_beta[sensor]**2)))
                values["outside_sensor_energy_removed_pct"] = float(100 * (1 - np.sum(effect[outside]**2) / np.sum(raw_beta[outside]**2)))
                for metric, value in values.items():
                    row = {**common, "metric": metric, "value": value}
                    if repetition == -1:
                        observed[(mode, metric)] = value
                    else:
                        rows.append(row)
        part = pd.DataFrame([row for row in rows if row["accession"] == f"OSD-{n}"])
        for (mode, metric), sub in part.groupby(["mode", "metric"]):
            values = sub.value.to_numpy()
            summaries.append(dict(accession=f"OSD-{n}", mode=mode, metric=metric,
                                  point_estimate=observed[(mode, metric)], bootstrap_lower95=float(np.quantile(values, .025)),
                                  bootstrap_upper95=float(np.quantile(values, .975)), repetitions=len(values),
                                  shared_HGC_resampled_once=n == "525", fixed_initial_gene_universe=True,
                                  normalization_refitted=True, ground_anchor_refitted=True, target_refitted=True))
        print(f"Joint normalization/anchor/target bootstrap: OSD-{n} / {repetitions} complete", flush=True)
    pd.DataFrame(rows).to_csv(HERE / "results/joint_bootstrap_repetitions.csv", index=False)
    pd.DataFrame(summaries).to_csv(HERE / "results/joint_bootstrap_summary.csv", index=False)


def geometry(datasets, annotation):
    records = []
    combinations = [(train, target) for train in ("561", "562") for target in ("561", "562", "563", "564", "612", "613", "525")]
    combinations.append(("525", "525"))
    for train_n, target_n in combinations:
        train, target = datasets[train_n], datasets[target_n]
        common = train["y"].index.intersection(target["y"].index)
        meta = train["meta"]
        mask = meta.flight.isin(["HGC", "VGC"]) if train_n == "525" else meta.flight == "GC"
        group = (meta.loc[mask, "flight"] == "HGC").to_numpy(int) if train_n == "525" else (meta.loc[mask, "cohort"] == "ISS_T").to_numpy(int)
        sensor = indices(common, annotation)["Sensor31"]
        v, w, gamma = learn_anchor(train["y"].loc[common, meta.index[mask]].to_numpy(), group, nuisance(meta.loc[mask]), sensor)
        # ||A-A.T||_F follows from the rank-two skew update, without allocating G x G.
        symmetry_departure = np.sqrt(max(0., 2 * (np.dot(gamma, gamma) - np.dot(w, gamma)**2)))
        records.append(dict(anchor_accession=f"OSD-{train_n}", target_accession=f"OSD-{target_n}", n_genes=len(common),
                            n_sensor=len(sensor), w_dot_gamma=float(w @ gamma), gamma_l2_norm=float(np.linalg.norm(gamma)),
                            score_spectral_gain=spectral_gain(w, gamma), half_score_spectral_gain=spectral_gain(w, gamma, .5),
                            bounded_score_spectral_gain=spectral_gain(w, bounded_gamma(w, gamma)),
                            projection_spectral_gain=1., score_symmetry_departure_frobenius=float(symmetry_departure),
                            sensor_identity_error=abs(float(w @ gamma) - 1)))
    pd.DataFrame(records).to_csv(HERE / "results/operator_geometry.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=200)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--skip-simulations", action="store_true")
    args = parser.parse_args()
    check_input_integrity()
    if not args.skip_simulations:
        benchmark(args.repeats)
    annotation = pd.read_csv(HERE / "inputs/gene_annotation.csv.gz")
    datasets = {n: load_dataset(n) for n in DATASETS}
    geometry(datasets, annotation)
    joint_bootstrap(datasets, annotation, args.bootstrap)


if __name__ == "__main__":
    main()
