"""Independent chunk/groupby reconstruction and direct allocation checks.

Does not import the analysis implementation. Reads only extension-local inputs.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
from pathlib import Path
from collections import defaultdict
import hashlib
import itertools
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

B = Path(__file__).resolve().parents[1]
I, T = B / "inputs", B / "tables"


def raw_reconstruct(path, defs, membership, expected_rows):
    key_sets = {bound: set(defs.loc[defs.boundary_type.eq(bound), "boundary_key"]) for bound in ["start", "end"]}
    wanted = set(membership.junction_id)
    retained, totals, all_ids = {}, {"start": {}, "end": {}}, set()
    rows = 0
    for chunk in pd.read_csv(path, index_col=0, chunksize=40000):
        columns = list(chunk.columns)
        assert np.all(chunk.to_numpy() >= 0) and not chunk.index.duplicated().any()
        assert not all_ids.intersection(chunk.index)
        all_ids.update(chunk.index)
        rows += len(chunk)
        parts = chunk.index.to_series().str.split(":", expand=True)
        for bound, col in [("start", 2), ("end", 3)]:
            keys = parts[0] + ":" + parts[1] + ":" + parts[col]
            ok = keys.isin(key_sets[bound])
            summed = chunk.loc[ok].groupby(keys.loc[ok]).sum()
            for key, row in summed.iterrows():
                if key not in totals[bound]:
                    totals[bound][key] = np.zeros(len(columns), dtype=np.int64)
                totals[bound][key] += row.to_numpy(dtype=np.int64)
        for jid, row in chunk.loc[chunk.index.isin(wanted)].iterrows():
            assert jid not in retained
            retained[jid] = row.to_numpy(dtype=np.int64)
    assert rows == expected_rows
    arr = np.stack([retained.get(jid, np.zeros(len(columns), dtype=np.int64)) for jid in membership.junction_id])
    offsets = np.r_[0, defs.n_alternatives.to_numpy().cumsum()]
    counts, all_counts, proportions = [], [], np.full_like(arr, np.nan, dtype=float)
    for c, row in enumerate(defs.itertuples()):
        sl = slice(offsets[c], offsets[c + 1])
        total = arr[sl].sum(axis=0)
        full = totals[row.boundary_type].get(row.boundary_key, np.zeros(len(columns), dtype=int))
        assert np.all(total <= full)
        proportions[sl] = np.divide(arr[sl], total, out=np.full_like(arr[sl], np.nan, dtype=float), where=total > 0)
        counts.append(total)
        all_counts.append(full)
    counts, all_counts = np.array(counts), np.array(all_counts)
    dominance = np.divide(counts, all_counts, out=np.zeros_like(counts, dtype=float), where=all_counts > 0)
    return dict(samples=columns, alternative_counts=arr, counts=counts, all_counts=all_counts, fractions=proportions, dominance=dominance)


def main():
    checks = []
    frozen = pd.read_csv(I / "frozen_input_manifest.csv")
    for row in frozen.itertuples():
        assert hashlib.sha256((B / row.file).read_bytes()).hexdigest() == row.sha256
    summary = json.loads((T / "summary.json").read_text())
    assert hashlib.sha256((B / "RETROSPECTIVE_PLAN_2026-10-15.md").read_bytes()).hexdigest() == summary["plan_sha256"]
    checks.append(dict(check="Frozen input and retrospective-plan SHA256", passed=True, detail="All 15 original input snapshots and the dated plan unchanged"))
    definitions = pd.read_csv(I / "RR23_frozen_cluster_definitions.csv", dtype={"chrom": str})
    definitions["boundary_key"] = definitions.chrom + ":" + definitions.strand.astype(str) + ":" + definitions.boundary.astype(str)
    membership = pd.read_csv(T / "frozen_alternative_membership.csv")
    offsets = np.r_[0, definitions.n_alternatives.to_numpy().cumsum()]
    raw = {}
    for tag, rows in [("RR23", 805562), ("RRRM2", 668237)]:
        raw[tag] = raw_reconstruct(I / f"{tag}_animal_junction_counts.csv.gz", definitions, membership, rows)
        x = raw[tag]
        for file, array in [("all_frozen_alternative_counts", x["alternative_counts"]), ("all_frozen_alternative_fractions", x["fractions"])]:
            saved = pd.read_csv(T / f"{tag}_{file}.csv.gz")
            assert list(saved.cluster_id) == list(membership.cluster_id) and list(saved.junction_id) == list(membership.junction_id)
            np.testing.assert_allclose(saved[x["samples"]].to_numpy(), array, atol=1e-14, rtol=0, equal_nan=True)
        for file, array in [("selected_counts", x["counts"]), ("all_boundary_counts", x["all_counts"]), ("dominance", x["dominance"])]:
            saved = pd.read_csv(T / f"{tag}_all_cluster_{file}.csv.gz", index_col=0).loc[definitions.cluster_id, x["samples"]]
            np.testing.assert_allclose(saved.to_numpy(), array, atol=1e-14, rtol=0)
        checks.append(dict(check=f"{tag} complete counts/compositions/boundary dominance", passed=True,
                           detail=f"Independently read {rows:,} nonnegative unique raw junction rows in chunks; checked all 923 clusters and 3,376 alternative occurrences"))
    rr, val = raw["RR23"], raw["RRRM2"]
    rr_meta = pd.read_csv(I / "RR23_animal_metadata.csv").set_index("sample").loc[rr["samples"]]
    val_all = pd.read_csv(I / "RRRM2_animal_metadata.csv").set_index("sample")
    val_meta = val_all.loc[val["samples"]]
    assert rr_meta["Source Name"].is_unique and val_meta["Source Name"].is_unique
    assert np.array_equal(val_all.included, ~(val_all.cohort.eq("ISS_T") & val_all.age.eq(29)))
    assert set(val_meta.index) == set(val_all.loc[val_all.included].index)
    assert len(val_meta) == 19 and val_meta.exposure.value_counts().to_dict() == {"FLT": 10, "GC": 9}
    rr_pass = (rr["counts"].min(1) >= 20) & (rr["dominance"].min(1) >= .9 - 1e-12)
    val_pass = (val["counts"].min(1) >= 20) & (val["dominance"].min(1) >= .9 - 1e-12)
    common = rr_pass & val_pass
    result = pd.read_csv(T / "all_923_cluster_directional_results.csv")
    np.testing.assert_array_equal(common, result.common_measurable)
    assert len(result) == 923 and result.test_family_includes_cluster.all()
    assert np.all(result.loc[~common, "exact_one_sided_directional_p"] == 1)
    assert result.loc[~common, "RRRM2_directional_effect"].isna().all()
    checks.append(dict(check="Distinct animals, inclusion and all measurement decisions", passed=True, detail="19 held-out animals; entire ISS_T_29 exclusion; all 923 family members retained with failed-gate p=1"))
    flt, gc = rr_meta.group.eq("FLT").to_numpy(), rr_meta.group.eq("GC").to_numpy()
    direction, centering, norms, projections = [], [], [], []
    for c in range(923):
        sl = slice(offsets[c], offsets[c + 1])
        delta = np.array([float(row[flt].mean() - row[gc].mean()) for row in rr["fractions"][sl]])
        norm = math.sqrt(sum(float(d) ** 2 for d in delta))
        u = delta / norm if norm > 1e-12 else np.zeros_like(delta)
        center = np.array([float(row[gc].mean()) for row in rr["fractions"][sl]])
        direction.extend(u)
        centering.extend(center)
        norms.append(norm)
        score = np.array([sum(float(a) * float(p - z) for a, p, z in zip(u, val["fractions"][sl, i], center)) for i in range(19)])
        if not common[c]:
            score[:] = np.nan
        projections.append(score)
    projections = np.array(projections)
    saved_direction = pd.read_csv(T / "RR23_frozen_training_vectors.csv.gz")
    # Different mean-reduction orders can be amplified by unit normalization of
    # tiny training effects. Use the declared 1e-12 inferential arithmetic tolerance.
    direction_max_error = float(np.max(np.abs(np.asarray(direction) - saved_direction.RR23_unit_direction.to_numpy())))
    np.testing.assert_allclose(direction, saved_direction.RR23_unit_direction, atol=1e-12, rtol=0)
    np.testing.assert_allclose(centering, saved_direction.RR23_GC_mean_fraction, atol=1e-14, rtol=0)
    np.testing.assert_allclose(norms, result.RR23_training_native_L2_norm, atol=1e-14, rtol=0)
    saved_scores = pd.read_csv(T / "RRRM2_cluster_directional_scores.csv.gz", index_col=0).loc[definitions.cluster_id, val["samples"]].to_numpy()
    np.testing.assert_allclose(projections, saved_scores, atol=1e-12, rtol=0, equal_nan=True)
    all_scores = projections[common].mean(0)
    score_table = pd.read_csv(T / "RRRM2_equal_cluster_animal_scores.csv").set_index("sample").loc[val["samples"]]
    np.testing.assert_allclose(all_scores, score_table.equal_cluster_directional_score, atol=1e-12, rtol=0)
    checks.append(dict(check="RR23-only unit directions and every held-out projection", passed=True, detail="All 923 unit vectors, every 19-animal measurable projection and complete equal-cluster score independently reproduced"))
    null_parts, obs, local_obs, designs = [], 0.0, np.zeros(common.sum()), []
    selected_rows = np.unique(np.linspace(0, common.sum() - 1, 17, dtype=int))
    local_parts = []
    for st in sorted(val_meta.stratum.unique()):
        ids = np.flatnonzero(val_meta.stratum.eq(st))
        lab = val_meta.iloc[ids].exposure.eq("FLT").to_numpy()
        nf, ng = int(lab.sum()), int((~lab).sum())
        w = len(ids) / len(val_meta)
        obs += w * (all_scores[ids[lab]].mean() - all_scores[ids[~lab]].mean())
        local_obs += w * (projections[common][:, ids[lab]].mean(1) - projections[common][:, ids[~lab]].mean(1))
        combinations = list(itertools.combinations(ids, nf))
        nv, lv = [], []
        for chosen in combinations:
            others = sorted(set(ids) - set(chosen))
            nv.append(w * (all_scores[list(chosen)].mean() - all_scores[others].mean()))
            lv.append(w * (projections[common][selected_rows][:, list(chosen)].mean(1) - projections[common][selected_rows][:, others].mean(1)))
        null_parts.append(nv)
        local_parts.append(lv)
        designs.append((ids, lab, w))
    null = np.array([sum(values) for values in itertools.product(*null_parts)])
    assert len(null) == 11200 == math.prod(len(part) for part in null_parts)
    np.testing.assert_allclose(np.sort(null), np.sort(pd.read_csv(T / "exact_global_null_all_11200.csv.gz").equal_cluster_directional_effect), atol=1e-14, rtol=0)
    pu = float((null >= obs - 1e-12).mean())
    pt = float((np.abs(null) >= abs(obs) - 1e-12).mean())
    assert np.isclose(pu, summary["global_exact_one_sided_p"]) and np.isclose(pt, summary["global_exact_absolute_two_sided_p"])
    assert np.isclose(obs, summary["global_equal_cluster_effect"])
    np.testing.assert_allclose(local_obs, result.loc[common, "RRRM2_directional_effect"], atol=1e-14, rtol=0)
    checks.append(dict(check="Independent complete global exact allocations", passed=True,
                       detail=f"All 11,200 direct stratified combinations; one-sided p={pu}, absolute two-sided p={pt}; every local effect checked"))
    local_null = np.array([sum(values) for values in itertools.product(*local_parts)])
    lp = (local_null >= local_obs[selected_rows] - 1e-12).mean(0)
    np.testing.assert_allclose(lp, result.loc[common, "exact_one_sided_directional_p"].iloc[selected_rows], atol=1e-14, rtol=0)
    p_all = result.exact_one_sided_directional_p.to_numpy()
    np.testing.assert_allclose(false_discovery_control(p_all, method="bh"), result.BH_full_923, atol=1e-14, rtol=0)
    np.testing.assert_allclose(false_discovery_control(p_all, method="by"), result.BY_full_923, atol=1e-14, rtol=0)
    ranked, last, h = np.argsort(p_all), 0.0, np.zeros(923)
    for k, j in enumerate(ranked):
        last = max(last, p_all[j] * (923 - k))
        h[j] = min(1.0, last)
    np.testing.assert_allclose(h, result.HOLM_full_923, atol=1e-14, rtol=0)
    assert 1 / 11200 > .05 / 923 and 923 / 11200 > .05
    checks.append(dict(check="Local exact tests, all three corrections and resolution bound", passed=True, detail="17 row-order-selected local exact p values; all 923 BH/BY/Holm values; Holm05 unattainable under 11,200 assignments"))
    saved_boot = np.load(T / "animal_bootstrap_cluster_effects.npz")
    vw, rw = saved_boot["RRRM2_bootstrap_contrast_weights"], saved_boot["RR23_bootstrap_contrast_weights"]
    assert vw.shape == (2000, 19) and rw.shape == (2000, 18)
    for ids, lab, w in designs:
        for group_ids, sign in [(ids[lab], 1), (ids[~lab], -1)]:
            multiplicities = vw[:, group_ids] * sign * len(group_ids) / w
            np.testing.assert_allclose(multiplicities, np.rint(multiplicities), atol=1e-12, rtol=0)
            np.testing.assert_allclose(multiplicities.sum(1), len(group_ids), atol=1e-12, rtol=0)
    local_boot = vw @ projections[common].T
    np.testing.assert_allclose(local_boot, saved_boot["conditional_fixed_direction"], atol=1e-14, rtol=0)
    boot_table = pd.read_csv(T / "animal_bootstrap_global_draws.csv.gz")
    np.testing.assert_allclose(local_boot.mean(1), boot_table.conditional_fixed_direction_global_effect, atol=1e-14, rtol=0)
    np.testing.assert_allclose(np.quantile(local_boot.mean(1), [.025, .975]), summary["conditional_fixed_direction_bootstrap_CI95"], atol=1e-14, rtol=0)
    rr_train_ix = np.flatnonzero(rr_meta.group.isin(["FLT", "GC"]))
    joint, common_ids = [], np.flatnonzero(common)
    for c in common_ids:
        sl = slice(offsets[c], offsets[c + 1])
        a = rw @ rr["fractions"][sl][:, rr_train_ix].T
        v = vw @ val["fractions"][sl].T
        norm = np.linalg.norm(a, axis=1)
        joint.append(np.divide((a * v).sum(1), norm, out=np.zeros(2000), where=norm > 1e-12))
    joint = np.array(joint).T
    np.testing.assert_allclose(joint, saved_boot["joint_retrained_direction"], atol=1e-12, rtol=0)
    np.testing.assert_allclose(joint.mean(1), boot_table.joint_retrained_direction_global_effect, atol=1e-14, rtol=0)
    np.testing.assert_allclose(np.quantile(joint.mean(1), [.025, .975]), summary["joint_retrained_direction_bootstrap_CI95"], atol=1e-14, rtol=0)
    checks.append(dict(check="Whole-animal conditional and two-cohort bootstrap", passed=True, detail="Verified cell multiplicities, every common-cluster effect in 2,000 draws, both global draw series and percentile intervals"))
    if (T / "all_3376_alternative_group_composition_summaries.csv.gz").exists():
        component = pd.read_csv(T / "all_3376_alternative_group_composition_summaries.csv.gz")
        for exposure in ["FLT", "GC", "VIV"]:
            ix = np.flatnonzero(rr_meta.group.eq(exposure))
            expected = np.array([sum(row[ix]) / len(ix) for row in rr["fractions"]])
            np.testing.assert_allclose(expected, component["RR23_mean_" + exposure], atol=1e-14, rtol=0)
        feature_common = np.repeat(common, definitions.n_alternatives)
        means = {}
        for exposure in ["FLT", "GC"]:
            expected = np.zeros(len(membership))
            for ids, lab, weight in designs:
                ix = ids[lab] if exposure == "FLT" else ids[~lab]
                expected += weight * np.array([sum(row[ix]) / len(ix) for row in val["fractions"]])
            np.testing.assert_allclose(expected[feature_common], component.loc[feature_common, "RRRM2_route_age_weighted_mean_" + exposure], atol=1e-14, rtol=0)
            means[exposure] = expected
        np.testing.assert_allclose((means["FLT"] - means["GC"])[feature_common], component.loc[feature_common, "RRRM2_weighted_native_FLT_minus_GC"], atol=1e-14, rtol=0)
        checks.append(dict(check="Complete descriptive alternative-level mean contrasts", passed=True, detail="Every frozen RR23 group mean and common-panel RRRM2 route-age-weighted native component mean/difference independently reproduced"))
    checks.append(dict(check="All full-family conclusions retained", passed=True, detail=f"{common.sum()} common clusters; no significance subset; all family corrections and failed gates retained"))
    pd.DataFrame(checks).to_csv(T / "independent_numerical_checks.csv", index=False)
    out = dict(status="PASS", all_checks_passed=True, analysis_code_imported=False, independent_check_count=len(checks), checks=checks,
               common_clusters=int(common.sum()), exact_assignments=11200, source_and_plan_hashes_unchanged=True,
               maximum_independent_unit_direction_absolute_difference=direction_max_error,
               declared_inferential_arithmetic_tolerance=1e-12)
    (T / "independent_numerical_checks.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
