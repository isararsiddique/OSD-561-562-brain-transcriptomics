"""Frozen native-composition directional transport; see retrospective plan.

This script reads only its extension-local inputs and writes script-relative outputs.
All allocation and bootstrap inference uses whole biological animals.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
from pathlib import Path
from collections import defaultdict
import csv
import gzip
import hashlib
import itertools
import json
import math
import platform
import numpy as np
import pandas as pd
import scipy

B = Path(__file__).resolve().parents[1]
I, T = B / "inputs", B / "tables"
TOL = 1e-12
SEED, NBOOT = 2026101503, 2000


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reconstruct(path, definitions, membership, expected_rows):
    """Stream the full actual count matrix, retaining all alternatives at needed boundaries."""
    keys = {bound: set(definitions.loc[definitions.boundary_type.eq(bound), "boundary_key"]) for bound in ["start", "end"]}
    junctions = set(membership.junction_id)
    kept, totals = {}, {"start": {}, "end": {}}
    rows, boundary_rows = 0, 0
    with gzip.open(path, "rt") as stream:
        reader = csv.reader(stream)
        samples = next(reader)[1:]
        n = len(samples)
        for row in reader:
            rows += 1
            chrom, strand, start, end = row[0].split(":")
            k = {"start": chrom + ":" + strand + ":" + start,
                 "end": chrom + ":" + strand + ":" + end}
            relevant = [bound for bound in ["start", "end"] if k[bound] in keys[bound]]
            if not relevant and row[0] not in junctions:
                continue
            vector = np.asarray(row[1:], dtype=np.int64)
            assert len(vector) == n and np.all(vector >= 0)
            boundary_rows += 1
            for bound in relevant:
                if k[bound] not in totals[bound]:
                    totals[bound][k[bound]] = np.zeros(n, dtype=np.int64)
                totals[bound][k[bound]] += vector
            if row[0] in junctions:
                assert row[0] not in kept
                kept[row[0]] = vector
    assert rows == expected_rows
    alternatives = np.stack([kept.get(j, np.zeros(n, dtype=np.int64)) for j in membership.junction_id])
    offsets = np.r_[0, definitions.n_alternatives.to_numpy().cumsum()]
    cluster_counts = np.add.reduceat(alternatives, offsets[:-1], axis=0)
    all_counts = np.stack([totals[row.boundary_type].get(row.boundary_key, np.zeros(n, dtype=np.int64)) for row in definitions.itertuples()])
    assert np.all(cluster_counts <= all_counts)
    dominance = np.divide(cluster_counts, all_counts, out=np.zeros_like(cluster_counts, dtype=float), where=all_counts > 0)
    fractions = np.divide(alternatives, np.repeat(cluster_counts, definitions.n_alternatives.to_numpy(), axis=0),
                          out=np.full_like(alternatives, np.nan, dtype=float),
                          where=np.repeat(cluster_counts, definitions.n_alternatives.to_numpy(), axis=0) > 0)
    return dict(samples=samples, counts=alternatives, selected_counts=cluster_counts, all_boundary_counts=all_counts,
                dominance=dominance, fractions=fractions, source_rows=rows, relevant_boundary_rows=boundary_rows)


def adjust(p, method):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind="stable")
    ranked, m = p[order], len(p)
    if method == "holm":
        q = np.maximum.accumulate(ranked * (m - np.arange(m)))
    else:
        factor = sum(1 / i for i in range(1, m + 1)) if method == "by" else 1
        q = np.minimum.accumulate((ranked * m * factor / (np.arange(m) + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1)
    return out


def stratified_weights(meta):
    weights, observed, design = [], np.zeros(len(meta)), []
    for st in sorted(meta.stratum.unique()):
        ids = np.flatnonzero(meta.stratum.eq(st))
        lab = meta.iloc[ids].exposure.eq("FLT").to_numpy()
        nf, ng = int(lab.sum()), int((~lab).sum())
        assert nf >= 2 and ng >= 2
        w = len(ids) / len(meta)
        combinations = list(itertools.combinations(ids, nf))
        d = np.zeros((len(combinations), len(meta)))
        d[:, ids] = -w / ng
        for j, selected in enumerate(combinations):
            d[j, list(selected)] = w / nf
        weights.append(d)
        observed[ids] = np.where(lab, w / nf, -w / ng)
        design.append(dict(stratum=st, n_animals=len(ids), n_FLT=nf, n_GC=ng, weight=w, assignments=len(d)))
    allocations = np.zeros((1, len(meta)))
    for d in weights:
        allocations = (allocations[:, None, :] + d[None, :, :]).reshape(-1, len(meta))
    assert len(allocations) == 11200 == math.prod(r["assignments"] for r in design)
    assert np.any(np.max(np.abs(allocations - observed), axis=1) <= TOL)
    return observed, allocations, pd.DataFrame(design)


def bootstrap_weights(meta, rng):
    result = np.zeros((NBOOT, len(meta)))
    for st in sorted(meta.stratum.unique()):
        ids = np.flatnonzero(meta.stratum.eq(st))
        w = len(ids) / len(meta)
        for exposure, sign in [("FLT", 1), ("GC", -1)]:
            ix = ids[meta.iloc[ids].exposure.eq(exposure).to_numpy()]
            result[:, ix] = sign * w * rng.multinomial(len(ix), np.full(len(ix), 1 / len(ix)), NBOOT) / len(ix)
    return result


def main():
    T.mkdir(exist_ok=True)
    frozen = pd.read_csv(I / "frozen_input_manifest.csv")
    for row in frozen.itertuples():
        assert digest(B / row.file) == row.sha256
    plan_hash = digest(B / "RETROSPECTIVE_PLAN_2026-10-15.md")
    definitions = pd.read_csv(I / "RR23_frozen_cluster_definitions.csv", dtype={"chrom": str})
    assert len(definitions) == 923 and definitions.cluster_id.is_unique and definitions.selected.all()
    definitions["boundary_key"] = definitions.chrom + ":" + definitions.strand.astype(str) + ":" + definitions.boundary.astype(str)
    long_rr = pd.read_csv(I / "RR23_frozen_alternative_compositions.csv.gz")
    membership = long_rr[["cluster_id", "alternative_index", "junction_id", "matching_gene_ids"]].drop_duplicates().sort_values(["cluster_id", "alternative_index"]).reset_index(drop=True)
    assert len(membership) == 3376 and not membership.duplicated(["cluster_id", "alternative_index"]).any()
    assert list(membership.cluster_id.drop_duplicates()) == list(definitions.cluster_id)
    assert np.array_equal(membership.groupby("cluster_id", sort=False).size().to_numpy(), definitions.n_alternatives)
    offsets = np.r_[0, definitions.n_alternatives.to_numpy().cumsum()]
    frozen_features = np.load(I / "RR23_frozen_features.npz")
    assert np.array_equal(frozen_features["offsets"], offsets)
    assert np.array_equal(frozen_features["cluster_ids"], definitions.cluster_id)
    rr = reconstruct(I / "RR23_animal_junction_counts.csv.gz", definitions, membership, 805562)
    val = reconstruct(I / "RRRM2_animal_junction_counts.csv.gz", definitions, membership, 668237)
    assert len(rr["samples"]) == 26 and len(val["samples"]) == 19
    assert np.array_equal(frozen_features["samples"], rr["samples"])
    np.testing.assert_allclose(rr["fractions"].T, frozen_features["H"] ** 2, atol=1e-14, rtol=0)
    np.testing.assert_array_equal(rr["selected_counts"].min(1), definitions.min_selected_reads)
    np.testing.assert_allclose(rr["dominance"].min(1), definitions.min_dominance, atol=1e-14, rtol=0)
    rr_meta = pd.read_csv(I / "RR23_animal_metadata.csv").set_index("sample", drop=False).loc[rr["samples"]]
    val_all = pd.read_csv(I / "RRRM2_animal_metadata.csv").set_index("sample", drop=False)
    expected_include = ~(val_all.cohort.eq("ISS_T") & val_all.age.eq(29))
    assert np.array_equal(expected_include, val_all.included)
    assert set(val["samples"]) == set(val_all.loc[expected_include].index)
    val_meta = val_all.loc[val["samples"]]
    assert rr_meta["Source Name"].is_unique and val_meta["Source Name"].is_unique
    assert rr_meta.group.value_counts().to_dict() == {"FLT": 9, "GC": 9, "VIV": 8}
    assert val_meta.exposure.value_counts().to_dict() == {"FLT": 10, "GC": 9}
    val_all.to_csv(T / "RRRM2_complete_animal_decisions.csv", index=False)
    rr_meta[["sample", "Source Name", "group", "Comment[Euthanasia Date]"]].to_csv(T / "RR23_animals.csv", index=False)
    membership.to_csv(T / "frozen_alternative_membership.csv", index=False)
    rr_measured = (rr["selected_counts"].min(1) >= 20) & (rr["dominance"].min(1) >= .9 - TOL)
    val_measured = (val["selected_counts"].min(1) >= 20) & (val["dominance"].min(1) >= .9 - TOL)
    assert rr_measured.all()
    common = rr_measured & val_measured
    assert common.any(), "No common measurable cluster: report measurement failure without fitting a score."
    gates = definitions.copy()
    gates["RR23_minimum_selected_reads"] = rr["selected_counts"].min(1)
    gates["RR23_minimum_dominance"] = rr["dominance"].min(1)
    gates["RRRM2_minimum_selected_reads"] = val["selected_counts"].min(1)
    gates["RRRM2_minimum_dominance"] = val["dominance"].min(1)
    gates["RRRM2_coverage_pass"] = val["selected_counts"].min(1) >= 20
    gates["RRRM2_dominance_pass"] = val["dominance"].min(1) >= .9 - TOL
    gates["common_measurable"] = common
    gates.to_csv(T / "all_923_cluster_measurement_gates.csv", index=False)
    for tag, cohort in [("RR23", rr), ("RRRM2", val)]:
        for name in ["counts", "fractions"]:
            frame = membership[["cluster_id", "alternative_index", "junction_id"]].copy()
            frame[cohort["samples"]] = cohort[name]
            frame.to_csv(T / f"{tag}_all_frozen_alternative_{name}.csv.gz", index=False)
        for name in ["selected_counts", "all_boundary_counts", "dominance"]:
            pd.DataFrame(cohort[name], index=definitions.cluster_id, columns=cohort["samples"]).to_csv(T / f"{tag}_all_cluster_{name}.csv.gz")
    rr_flt, rr_gc = rr_meta.group.eq("FLT").to_numpy(), rr_meta.group.eq("GC").to_numpy()
    rr_direction_raw = rr["fractions"][:, rr_flt].mean(1) - rr["fractions"][:, rr_gc].mean(1)
    training_norms = np.sqrt(np.add.reduceat(rr_direction_raw ** 2, offsets[:-1]))
    unit = np.divide(rr_direction_raw, np.repeat(training_norms, definitions.n_alternatives),
                     out=np.zeros_like(rr_direction_raw), where=np.repeat(training_norms, definitions.n_alternatives) > TOL)
    unit_sums = np.add.reduceat(unit, offsets[:-1])
    assert np.max(np.abs(unit_sums)) < 1e-10
    assert np.allclose(np.add.reduceat(unit ** 2, offsets[:-1]), np.where(training_norms > TOL, 1, 0), atol=1e-12)
    rr_gc_mean = rr["fractions"][:, rr_gc].mean(1)
    direction_table = membership.copy()
    direction_table["RR23_FLT_minus_GC_composition_difference"] = rr_direction_raw
    direction_table["RR23_unit_direction"] = unit
    direction_table["RR23_GC_mean_fraction"] = rr_gc_mean
    direction_table.to_csv(T / "RR23_frozen_training_vectors.csv.gz", index=False)
    cluster_scores_all = np.add.reduceat(unit[:, None] * (val["fractions"] - rr_gc_mean[:, None]), offsets[:-1], axis=0)
    cluster_scores_all[~common] = np.nan
    pd.DataFrame(cluster_scores_all, index=definitions.cluster_id, columns=val["samples"]).to_csv(T / "RRRM2_cluster_directional_scores.csv.gz")
    weights, allocations, design = stratified_weights(val_meta)
    design.to_csv(T / "RRRM2_stratified_design.csv", index=False)
    alternative_summary = direction_table.copy()
    for exposure in ["FLT", "GC", "VIV"]:
        alternative_summary["RR23_mean_" + exposure] = rr["fractions"][:, rr_meta.group.eq(exposure).to_numpy()].mean(1)
    feature_measurable = np.repeat(common, definitions.n_alternatives)
    for exposure, group_weights in [("FLT", np.maximum(weights, 0)), ("GC", -np.minimum(weights, 0))]:
        means = val["fractions"] @ group_weights
        alternative_summary["RRRM2_route_age_weighted_mean_" + exposure] = np.where(feature_measurable, means, np.nan)
    alternative_summary["RRRM2_weighted_native_FLT_minus_GC"] = alternative_summary.RRRM2_route_age_weighted_mean_FLT - alternative_summary.RRRM2_route_age_weighted_mean_GC
    alternative_summary["common_cluster_measurable"] = feature_measurable
    alternative_summary.to_csv(T / "all_3376_alternative_group_composition_summaries.csv.gz", index=False)
    cluster_scores = cluster_scores_all[common]
    local_effect = cluster_scores @ weights
    global_scores = cluster_scores.mean(0)
    observed = float(global_scores @ weights)
    null = allocations @ global_scores
    global_p = float(np.mean(null >= observed - TOL))
    global_p_two = float(np.mean(np.abs(null) >= abs(observed) - TOL))
    p_common = np.zeros(common.sum())
    local_two = np.zeros(common.sum())
    for start in range(0, common.sum(), 32):
        local_null = cluster_scores[start:start + 32] @ allocations.T
        delta = local_effect[start:start + 32, None]
        p_common[start:start + 32] = (local_null >= delta - TOL).mean(1)
        local_two[start:start + 32] = (np.abs(local_null) >= np.abs(delta) - TOL).mean(1)
    p_all, two_all = np.ones(923), np.ones(923)
    p_all[common], two_all[common] = p_common, local_two
    rng = np.random.default_rng(SEED)
    val_boot_weights = bootstrap_weights(val_meta, rng)
    boot_local = val_boot_weights @ cluster_scores.T
    boot_global = boot_local.mean(1)
    boot_positive_fraction = (boot_local > 0).mean(1)
    # Secondary training uncertainty: jointly resample independent cohort animals.
    rr_train = rr_meta.loc[rr_meta.group.isin(["FLT", "GC"])].copy()
    rr_train["exposure"], rr_train["stratum"] = rr_train.group, "all"
    rr_train_ix = rr_meta.index.get_indexer(rr_train.index)
    rr_boot_weights = bootstrap_weights(rr_train, rng)
    feature_common = np.repeat(common, definitions.n_alternatives)
    rr_boot_delta = rr_boot_weights @ rr["fractions"][feature_common][:, rr_train_ix].T
    val_boot_delta = val_boot_weights @ val["fractions"][feature_common].T
    common_sizes = definitions.n_alternatives.to_numpy()[common]
    common_offsets = np.r_[0, common_sizes.cumsum()]
    rr_boot_norm = np.sqrt(np.add.reduceat(rr_boot_delta ** 2, common_offsets[:-1], axis=1))
    numerator = np.add.reduceat(rr_boot_delta * val_boot_delta, common_offsets[:-1], axis=1)
    joint_local = np.divide(numerator, rr_boot_norm, out=np.zeros_like(numerator), where=rr_boot_norm > TOL)
    joint_global = joint_local.mean(1)
    val_native_effect = val["fractions"][feature_common] @ weights
    val_native_norm = np.sqrt(np.add.reduceat(val_native_effect ** 2, common_offsets[:-1]))
    alignment = np.divide(local_effect, val_native_norm, out=np.full_like(local_effect, np.nan), where=val_native_norm > TOL)
    out = gates.copy()
    out["RR23_training_native_L2_norm"] = training_norms
    out["RR23_zero_direction"] = training_norms <= TOL
    for name in ["RRRM2_directional_effect", "RRRM2_native_effect_L2_norm", "RRRM2_direction_cosine", "conditional_bootstrap_CI95_low", "conditional_bootstrap_CI95_high"]:
        out[name] = np.nan
    out.loc[common, "RRRM2_directional_effect"] = local_effect
    out.loc[common, "RRRM2_native_effect_L2_norm"] = val_native_norm
    out.loc[common, "RRRM2_direction_cosine"] = alignment
    out.loc[common, "conditional_bootstrap_CI95_low"] = np.quantile(boot_local, .025, axis=0)
    out.loc[common, "conditional_bootstrap_CI95_high"] = np.quantile(boot_local, .975, axis=0)
    out["exact_one_sided_directional_p"] = p_all
    out["exact_absolute_two_sided_p"] = two_all
    for method in ["bh", "by", "holm"]:
        out[method.upper() + "_full_923"] = adjust(p_all, method)
    out["test_family_includes_cluster"] = True
    out.to_csv(T / "all_923_cluster_directional_results.csv", index=False)
    all_score_rows = val_meta[["sample", "Source Name", "exposure", "stratum", "cohort", "age"]].copy()
    all_score_rows["equal_cluster_directional_score"] = global_scores
    all_score_rows["common_cluster_count"] = int(common.sum())
    all_score_rows.to_csv(T / "RRRM2_equal_cluster_animal_scores.csv", index=False)
    pd.DataFrame(dict(assignment=np.arange(len(null)), equal_cluster_directional_effect=null,
                      as_or_more_positive=null >= observed - TOL,
                      as_or_more_absolute=np.abs(null) >= abs(observed) - TOL)).to_csv(T / "exact_global_null_all_11200.csv.gz", index=False)
    pd.DataFrame(dict(draw=np.arange(NBOOT), conditional_fixed_direction_global_effect=boot_global,
                      conditional_positive_cluster_fraction=boot_positive_fraction,
                      joint_retrained_direction_global_effect=joint_global,
                      joint_retrained_positive_cluster_fraction=(joint_local > 0).mean(1))).to_csv(T / "animal_bootstrap_global_draws.csv.gz", index=False)
    np.savez_compressed(T / "animal_bootstrap_cluster_effects.npz", conditional_fixed_direction=boot_local,
                        joint_retrained_direction=joint_local, common_cluster_ids=definitions.cluster_id.to_numpy()[common],
                        RRRM2_bootstrap_contrast_weights=val_boot_weights, RR23_bootstrap_contrast_weights=rr_boot_weights)
    stratum_rows = []
    for row in design.itertuples():
        ix = np.flatnonzero(val_meta.stratum.eq(row.stratum))
        lab = val_meta.iloc[ix].exposure.eq("FLT").to_numpy()
        d = np.where(lab, 1 / lab.sum(), -1 / (~lab).sum())
        bs = val_boot_weights[:, ix] @ global_scores[ix] / row.weight
        stratum_rows.append(dict(stratum=row.stratum, n_FLT=row.n_FLT, n_GC=row.n_GC, weight=row.weight,
                                 effect=float(global_scores[ix] @ d), CI95_low=float(np.quantile(bs, .025)), CI95_high=float(np.quantile(bs, .975))))
    pd.DataFrame(stratum_rows).to_csv(T / "descriptive_stratum_effects.csv", index=False)
    sensitivity = []
    for name, depth, dom in [("primary", 20, .90), ("dominance95", 20, .95), ("depth50", 50, .90)]:
        keep = (rr["selected_counts"].min(1) >= depth) & (rr["dominance"].min(1) >= dom - TOL) & (val["selected_counts"].min(1) >= depth) & (val["dominance"].min(1) >= dom - TOL)
        if keep.any():
            score = cluster_scores_all[keep].mean(0)
            delta, perm = float(score @ weights), allocations @ score
            sensitivity.append(dict(panel=name, read_threshold=depth, dominance_threshold=dom,
                                    common_clusters=int(keep.sum()), global_effect=delta,
                                    exact_one_sided_p=float((perm >= delta - TOL).mean()),
                                    exact_absolute_two_sided_p=float((np.abs(perm) >= abs(delta) - TOL).mean())))
        else:
            sensitivity.append(dict(panel=name, read_threshold=depth, dominance_threshold=dom, common_clusters=0,
                                    global_effect=np.nan, exact_one_sided_p=1.0, exact_absolute_two_sided_p=1.0))
    sensitivities = pd.DataFrame(sensitivity)
    sensitivities["Holm_two_secondary_sensitivities"] = np.nan
    sensitivities.loc[1:, "Holm_two_secondary_sensitivities"] = adjust(sensitivities.loc[1:, "exact_one_sided_p"], "holm")
    sensitivities.to_csv(T / "measurement_sensitivities.csv", index=False)
    # Historically selected Tcaf1 cluster: measurement governs whether it is described.
    tcaf = definitions.gene_symbol.eq("Tcaf1").to_numpy()
    tcaf_results = out.loc[tcaf].copy()
    tcaf_results.to_csv(T / "Tcaf1_measurement_and_complete_family_result.csv", index=False)
    alternative_summary.loc[alternative_summary.cluster_id.isin(definitions.loc[tcaf & common, "cluster_id"])].to_csv(T / "Tcaf1_all_frozen_alternative_mean_contrasts.csv", index=False)
    case_rows = []
    if (tcaf & common).any():
        for c in np.flatnonzero(tcaf & common):
            for tag, cohort, meta in [("RR23", rr, rr_meta.assign(exposure=rr_meta.group)), ("RRRM2", val, val_meta)]:
                for j in range(offsets[c], offsets[c + 1]):
                    for animal, sample in enumerate(cohort["samples"]):
                        case_rows.append(dict(cohort=tag, cluster_id=definitions.iloc[c].cluster_id,
                                              alternative_index=int(membership.iloc[j].alternative_index),
                                              junction_id=membership.iloc[j].junction_id, sample=sample,
                                              exposure=meta.loc[sample, "exposure"],
                                              stratum=meta.loc[sample, "stratum"] if "stratum" in meta else "all",
                                              count=int(cohort["counts"][j, animal]), fraction=float(cohort["fractions"][j, animal]),
                                              selected_total=int(cohort["selected_counts"][c, animal]), dominance=float(cohort["dominance"][c, animal])))
        pd.DataFrame(case_rows).to_csv(T / "Tcaf1_all_measurable_alternatives_and_animals.csv", index=False)
    binary = json.loads((I / "prior_binary_transport_summary.json").read_text())
    transport = pd.DataFrame([
        dict(observable="binary local junction pairs", frozen_family=binary["frozen_definitions"], RR23_measured=binary["RR23_measured_all26"], common_measured=binary["common_measured_events"], fraction_of_RR23_measured=binary["common_transport_fraction"]),
        dict(observable="multiway fixed-alternative compositions", frozen_family=923, RR23_measured=923, common_measured=int(common.sum()), fraction_of_RR23_measured=float(common.mean()))])
    transport.to_csv(T / "binary_multiway_measurement_transport_comparison.csv", index=False)
    summary = dict(status="retrospective full-family multiway independent-launch directional transport", author="Isarar Siddique", affiliation="", email="",
                   analysis_date="2026-10-15", plan_sha256=plan_hash, frozen_clusters=923, fixed_alternative_occurrences=3376,
                   RR23_animals=26, RR23_training_FLT=9, RR23_training_GC=9, RRRM2_animals=19, RRRM2_FLT=10, RRRM2_GC=9,
                   exact_stratified_assignments=11200, common_measured_clusters=int(common.sum()), common_measurement_fraction=float(common.mean()),
                   RRRM2_coverage_fail=int((~gates.RRRM2_coverage_pass).sum()), RRRM2_dominance_fail=int((~gates.RRRM2_dominance_pass).sum()),
                   RRRM2_both_measurement_fail=int((~gates.RRRM2_coverage_pass & ~gates.RRRM2_dominance_pass).sum()),
                   zero_RR23_training_directions=int((training_norms <= TOL).sum()),
                   normalization="native RR23 FLT−GC composition contrast divided by Euclidean norm; zero norm≤1e−12 gives zero vector",
                   global_equal_cluster_effect=observed, global_equal_cluster_effect_pp=100 * observed,
                   global_exact_one_sided_p=global_p, global_exact_absolute_two_sided_p=global_p_two,
                   global_one_sided_extreme_allocations=int((null >= observed - TOL).sum()),
                   global_two_sided_extreme_allocations=int((np.abs(null) >= abs(observed) - TOL).sum()),
                   conditional_fixed_direction_bootstrap_CI95=np.quantile(boot_global, [.025, .975]).tolist(),
                   joint_retrained_direction_bootstrap_CI95=np.quantile(joint_global, [.025, .975]).tolist(),
                   positive_directional_cluster_fraction=float((local_effect > 0).mean()),
                   positive_directional_cluster_fraction_conditional_bootstrap_CI95=np.quantile(boot_positive_fraction, [.025, .975]).tolist(),
                   mean_native_cosine_alignment=float(np.nanmean(alignment)),
                   local_BH_full923_05=int((out.BH_full_923 < .05).sum()),
                   local_BY_full923_05=int((out.BY_full_923 < .05).sum()),
                   local_Holm_full923_05=int((out.HOLM_full_923 < .05).sum()),
                   smallest_local_one_sided_p=float(p_common.min()),
                   minimum_attainable_raw_one_sided_p=1 / 11200,
                   Holm_first_step_alpha05_threshold=.05 / 923,
                   minimum_attainable_first_Holm_adjusted_p=923 / 11200,
                   local_Holm05_attainable_with_this_design=False,
                   Tcaf1_frozen_clusters=int(tcaf.sum()), Tcaf1_common_measurable_clusters=int((tcaf & common).sum()),
                   bootstrap_draws=NBOOT, bootstrap_seed=SEED, bootstrap_unit="whole animal within exposure and route-age cells",
                   primary_global_positive_transfer_supported=global_p < .05,
                   no_RR23_significance_filter=True, no_gene_or_cluster_independence_p_values=True,
                   RR23_flight_and_collection_day_aliased=True, causal_flight_identification=False,
                   separate_retrospective_global_and_local_families=True)
    (T / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (T / "software_versions.json").write_text(json.dumps(dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__), indent=2) + "\n")
    after = []
    for row in frozen.itertuples():
        assert digest(B / row.file) == row.sha256
        after.append(dict(file=row.file, sha256_before=row.sha256, sha256_after=digest(B / row.file), source_unchanged=True))
    assert digest(B / "RETROSPECTIVE_PLAN_2026-10-15.md") == plan_hash
    pd.DataFrame(after).to_csv(T / "source_hashes_before_after.csv", index=False)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
