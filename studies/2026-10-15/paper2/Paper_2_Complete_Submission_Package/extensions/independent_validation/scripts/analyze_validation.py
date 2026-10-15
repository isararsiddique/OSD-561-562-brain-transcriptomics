"""Retrospective full-panel independent-launch validation; see frozen plan.

Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1.
All inferential resampling is at the biological animal level.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import rankdata

B = Path(__file__).resolve().parents[1]
T = B / "tables"
SEED, NBOOT = 20261015, 2000


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reconstruct(counts, definitions):
    a = counts.reindex(definitions.junction_a).fillna(0).to_numpy()
    b = counts.reindex(definitions.junction_b).fillna(0).to_numpy()
    pair = a + b
    frac = np.divide(a, pair, out=np.full_like(a, np.nan, dtype=float), where=pair > 0)
    parts = counts.index.str.split(":", expand=True)
    info = pd.DataFrame({"chr": parts.get_level_values(0), "strand": parts.get_level_values(1),
                         "start": parts.get_level_values(2), "end": parts.get_level_values(3)}, index=counts.index)
    all_pair = np.zeros_like(pair)
    for boundary in ["start", "end"]:
        sums = counts.groupby(info.chr + ":" + info.strand + ":" + info[boundary]).sum()
        ix = np.flatnonzero(definitions.boundary_type.eq(boundary))
        rows = definitions.iloc[ix]
        keys = rows.chr.astype(str) + ":" + rows.strand.astype(str) + ":" + rows.boundary.astype(str)
        all_pair[ix] = sums.reindex(keys).fillna(0).to_numpy()
    dominance = np.divide(pair, all_pair, out=np.zeros_like(pair, dtype=float), where=all_pair > 0)
    assert np.all(pair <= all_pair)
    return frac, pair, dominance


def stratified_weights(meta):
    groups, observed = [], np.zeros(len(meta))
    strata = []
    for st in sorted(meta.stratum.unique()):
        ix = np.flatnonzero(meta.stratum.eq(st))
        flt = meta.iloc[ix].exposure.eq("FLT").to_numpy()
        nf, ng = int(flt.sum()), int((~flt).sum())
        assert nf > 0 and ng > 0
        w = len(ix) / len(meta)
        combos = list(itertools.combinations(ix, nf))
        mat = np.full((len(combos), len(meta)), 0.0)
        mat[:, ix] = -w / ng
        for j, selected in enumerate(combos):
            mat[j, list(selected)] = w / nf
        groups.append(mat)
        observed[ix] = np.where(flt, w / nf, -w / ng)
        strata.append(dict(stratum=st, n_animals=len(ix), n_FLT=nf, n_GC=ng,
                           weight=w, label_assignments=len(combos)))
    perms = np.zeros((1, len(meta)))
    for mat in groups:
        perms = (perms[:, None, :] + mat[None, :, :]).reshape(-1, len(meta))
    assert len(perms) == math.prod(r["label_assignments"] for r in strata)
    assert np.allclose(perms.sum(1), 0)
    assert np.any(np.max(np.abs(perms - observed), axis=1) < 1e-12)
    return observed, perms, pd.DataFrame(strata)


def bootstrap_weights(meta, rng):
    mat = np.zeros((NBOOT, len(meta)))
    for st in sorted(meta.stratum.unique()):
        ix = np.flatnonzero(meta.stratum.eq(st))
        w = len(ix) / len(meta)
        for exposure, sign in [("FLT", 1), ("GC", -1)]:
            ids = ix[meta.iloc[ix].exposure.eq(exposure).to_numpy()]
            assert len(ids) >= 2
            mat[:, ids] = sign * w * rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), NBOOT) / len(ids)
    return mat


def rho(a, b):
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def holm(p):
    p = np.asarray(p)
    order = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate(p[order] * (len(p) - np.arange(len(p)))))
    out = np.empty_like(adj)
    out[order] = adj
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True, help="Beta analysis/package root containing inputs/ and results/")
    args = parser.parse_args()
    RR = args.base.resolve()
    VAL = B / "inputs"
    T.mkdir(exist_ok=True)
    paths = [RR / "inputs/primary_all26_event_definitions.csv.gz",
             RR / "results/all26_animal_fractions.csv.gz", RR / "results/measurement_gates.csv.gz",
             RR / "inputs/animal_metadata.csv", RR / "inputs/all_animal_junction_counts.csv.gz",
             VAL / "OSD-561_animal_junction_counts.csv.gz", VAL / "OSD-561_sample_decisions.csv",
             B / "RETROSPECTIVE_ANALYSIS_PLAN.md", VAL / "RRRM2_publisher_provenance.csv",
             VAL / "frozen_input_provenance.json", VAL / "OSD-561_publisher_inventory.json",
             VAL / "GLDS-556_rna_seq_processed_md5sum_GLbulkRNAseq.tsv"]
    provenance = json.loads((VAL / "frozen_input_provenance.json").read_text())
    assert provenance["collapsed_matrix_independently_matches_all_retained_STAR_files"]
    assert provenance["retained_junction_files"] == 19
    for record in provenance["frozen_inputs"]:
        assert sha(VAL / record["file"]) == record["sha256"]
    manifest = pd.DataFrame([dict(scope="extension" if p.is_relative_to(B) else "beta",
                                 path=str(p.relative_to(B if p.is_relative_to(B) else RR)), bytes=p.stat().st_size,
                                 sha256_before=sha(p)) for p in paths])
    definitions = pd.read_csv(paths[0], dtype={"chr": str}).set_index("event_id", drop=False)
    assert len(definitions) == 5724 and definitions.index.is_unique
    saved_frac = pd.read_csv(paths[1], index_col=0).loc[definitions.index]
    rr_gates = pd.read_csv(paths[2], index_col=0).loc[definitions.index]
    rr_counts = pd.read_csv(paths[4], index_col=0)
    rr_meta = pd.read_csv(paths[3]).set_index("sample", drop=False).loc[rr_counts.columns]
    val_counts = pd.read_csv(paths[5], index_col=0)
    val_all = pd.read_csv(paths[6]).set_index("sample", drop=False)
    exclusion = ~(val_all.cohort.eq("ISS_T") & val_all.age.eq(29))
    assert (exclusion == val_all.included).all()
    assert val_all.loc[~exclusion & val_all.exposure.eq("FLT"), "location_conflict"].all()
    assert set(val_counts.columns) == set(val_all.loc[exclusion].index)
    val_meta = val_all.loc[val_counts.columns].copy()
    assert len(val_meta) == 19 and val_meta["Source Name"].is_unique
    assert val_meta.exposure.value_counts().to_dict() == {"FLT": 10, "GC": 9}
    assert rr_meta.group.value_counts().to_dict() == {"FLT": 9, "GC": 9, "VIV": 8}
    assert rr_meta["Source Name"].is_unique
    val_all["analysis_decision"] = np.where(exclusion, "include", "exclude entire ISS_T age29 stratum: conflicting FLT location")
    val_all.to_csv(T / "RRRM2_animal_inclusion_decisions.csv", index=False)
    rr_meta[["sample", "Source Name", "group", "Comment[Euthanasia Date]",
             "Characteristics[Sex]", "Characteristics[Strain]"]].to_csv(T / "RR23_animal_metadata.csv", index=False)
    rr_frac, rr_pair, rr_dom = reconstruct(rr_counts, definitions)
    val_frac, val_pair, val_dom = reconstruct(val_counts, definitions)
    np.testing.assert_allclose(rr_frac, saved_frac[rr_counts.columns].to_numpy(), rtol=0, atol=1e-14)
    np.testing.assert_allclose(rr_pair.min(1), rr_gates.minimum_pair_reads_all26)
    np.testing.assert_allclose(rr_dom.min(1), rr_gates.minimum_dominance_all26, atol=1e-14)
    rr_measured = (rr_pair.min(1) >= 20) & (rr_dom.min(1) >= .9)
    assert np.array_equal(rr_measured, rr_gates.measured_primary_all26.to_numpy())
    assert rr_measured.sum() == 5249
    val_measured = (val_pair.min(1) >= 20) & (val_dom.min(1) >= .9)
    common = rr_measured & val_measured
    mapping = definitions.copy()
    mapping["RR23_minimum_pair_reads"] = rr_pair.min(1)
    mapping["RR23_minimum_dominance"] = rr_dom.min(1)
    mapping["RR23_measured_all26"] = rr_measured
    mapping["RRRM2_minimum_pair_reads"] = val_pair.min(1)
    mapping["RRRM2_minimum_dominance"] = val_dom.min(1)
    mapping["RRRM2_coverage20_pass"] = val_pair.min(1) >= 20
    mapping["RRRM2_dominance90_pass"] = val_dom.min(1) >= .9
    mapping["RRRM2_measured"] = val_measured
    mapping["common_primary_panel"] = common
    mapping.to_csv(T / "complete_feature_mapping_and_gates.csv.gz", index=False)
    for name, arr in [("fractions", val_frac), ("pair_counts", val_pair), ("pair_dominance", val_dom)]:
        pd.DataFrame(arr, index=definitions.index, columns=val_meta.index).to_csv(T / f"RRRM2_all26_{name}.csv.gz")
    rr_train_meta = rr_meta.loc[rr_meta.group.isin(["FLT", "GC"])].copy()
    rr_train_meta["exposure"] = rr_train_meta.group
    rr_train_meta["stratum"] = "all"
    rr_ix = rr_meta.index.get_indexer(rr_train_meta.index)
    rr_x = rr_frac[:, rr_ix]
    # No RR23 permutations needed: learned directions are trained once and fixed.
    rr_lab = rr_train_meta.exposure.eq("FLT").to_numpy()
    rr_w = np.where(rr_lab, 1 / rr_lab.sum(), -1 / (~rr_lab).sum())
    val_w, val_perm, strata = stratified_weights(val_meta)
    strata.to_csv(T / "RRRM2_stratum_design.csv", index=False)
    d_rr, d_val = rr_x @ rr_w, val_frac[common] @ val_w
    all_dval = np.full(len(definitions), np.nan)
    all_dval[common] = d_val
    rng = np.random.default_rng(SEED)
    rr_bw = bootstrap_weights(rr_train_meta, rng)
    val_bw = bootstrap_weights(val_meta, rng)
    # Reuse identical animal bootstrap draws across measurement sensitivities.
    rr_boot = rr_bw @ rr_x[common].T
    val_boot = val_bw @ val_frac[common].T
    rho_draws = np.array([rho(a, b) for a, b in zip(rr_boot, val_boot)])
    concordance_draws = (rr_boot * val_boot > 0).mean(1)
    directions = np.sign(d_rr[common])
    rr_gc_center = rr_x[common][:, ~rr_lab].mean(1)
    scores = directions @ (val_frac[common] - rr_gc_center[:, None]) / common.sum()
    score_delta = float(scores @ val_w)
    score_null = val_perm @ scores
    p_up = float(np.mean(score_null >= score_delta - 1e-12))
    p_two = float(np.mean(np.abs(score_null) >= abs(score_delta) - 1e-12))
    score_boot = val_bw @ scores
    assert np.isclose(score_delta, np.mean(directions * d_val), atol=1e-14)
    score_table = val_meta[["sample", "Source Name", "exposure", "stratum", "cohort", "age"]].copy()
    score_table["learned_direction_score"] = scores
    score_table["fixed_event_count"] = int(common.sum())
    score_table.to_csv(T / "RRRM2_heldout_animal_scores.csv", index=False)
    pd.DataFrame(dict(assignment=np.arange(len(score_null)), score_delta=score_null,
                      as_or_more_positive=score_null >= score_delta - 1e-12,
                      as_or_more_absolute=np.abs(score_null) >= abs(score_delta) - 1e-12)).to_csv(T / "RRRM2_exact_score_null.csv.gz", index=False)
    pd.DataFrame(dict(draw=np.arange(NBOOT), spearman_rho=rho_draws,
                      same_direction_fraction=concordance_draws, fixed_score_delta=score_boot)).to_csv(T / "animal_bootstrap_draws.csv.gz", index=False)
    # Secondary full common-panel event map; no RR23 p-value or effect-size filtering.
    event_p = np.zeros(common.sum())
    event_two = np.zeros(common.sum())
    common_x = val_frac[common]
    for start in range(0, len(common_x), 64):
        null = common_x[start:start + 64] @ val_perm.T
        effect = d_val[start:start + 64, None]
        sign = directions[start:start + 64, None]
        event_p[start:start + 64] = (sign * null >= sign * effect - 1e-12).mean(1)
        event_two[start:start + 64] = (np.abs(null) >= np.abs(effect) - 1e-12).mean(1)
    effects = definitions.loc[common, ["event_id", "gene_id", "gene_symbol", "boundary_type"]].copy()
    effects["RR23_delta_FLT_minus_GC"] = d_rr[common]
    effects["RRRM2_delta_stratified_FLT_minus_GC"] = d_val
    effects["RR23_learned_direction"] = directions
    effects["same_direction"] = d_rr[common] * d_val > 0
    effects["RRRM2_exact_directional_p"] = event_p
    effects["RRRM2_directional_Holm_full_common_panel"] = holm(event_p)
    effects["RRRM2_exact_two_sided_p"] = event_two
    effects["RR23_GC_center"] = rr_gc_center
    effects.to_csv(T / "complete_common_panel_effects.csv.gz", index=False)
    rows = []
    for st in strata.stratum:
        ix = np.flatnonzero(val_meta.stratum.eq(st))
        lab = val_meta.iloc[ix].exposure.eq("FLT").to_numpy()
        sw = np.where(lab, 1 / lab.sum(), -1 / (~lab).sum())
        # Bootstrap coefficients already carry global stratum weight.
        w = len(ix) / len(val_meta)
        bs = val_bw[:, ix] @ scores[ix] / w
        rows.append(dict(stratum=st, score_delta=float(scores[ix] @ sw),
                         CI95_low=float(np.quantile(bs, .025)), CI95_high=float(np.quantile(bs, .975)),
                         n_FLT=int(lab.sum()), n_GC=int((~lab).sum()), weight=w))
    pd.DataFrame(rows).to_csv(T / "descriptive_stratum_score_effects.csv", index=False)
    sensitivities = []
    for name, depth, dom in [("primary20_dominance90", 20, .90),
                             ("dominance95", 20, .95), ("depth50", 50, .90)]:
        keep = (rr_pair.min(1) >= depth) & (rr_dom.min(1) >= dom) & (val_pair.min(1) >= depth) & (val_dom.min(1) >= dom)
        a, z = d_rr[keep], val_frac[keep] @ val_w
        sign = np.sign(a)
        center = rr_x[keep][:, ~rr_lab].mean(1)
        s = sign @ (val_frac[keep] - center[:, None]) / keep.sum()
        obs, nul = float(s @ val_w), val_perm @ s
        sensitivities.append(dict(panel=name, pair_read_threshold=depth, dominance_threshold=dom,
                                  n_events=int(keep.sum()), spearman_rho=rho(a, z),
                                  same_direction_fraction=float((a * z > 0).mean()), score_delta=obs,
                                  exact_one_sided_p=float(np.mean(nul >= obs - 1e-12)),
                                  exact_two_sided_p=float(np.mean(np.abs(nul) >= abs(obs) - 1e-12))))
    sens = pd.DataFrame(sensitivities)
    sens["Holm_two_secondary_sensitivity_tests"] = np.nan
    sens.loc[1:, "Holm_two_secondary_sensitivity_tests"] = holm(sens.loc[1:, "exact_one_sided_p"])
    sens.to_csv(T / "measurement_sensitivities.csv", index=False)
    flow = [dict(step="Frozen RR23 all26 definitions", events=len(definitions)),
            dict(step="RR23 measured all26", events=int(rr_measured.sum())),
            dict(step="RRRM2 measured from definitions", events=int(val_measured.sum())),
            dict(step="Common primary measured panel", events=int(common.sum()))]
    pd.DataFrame(flow).to_csv(T / "measurement_transport_counts.csv", index=False)
    summary = dict(analysis_date="2026-10-15", status="retrospective full-panel independent-launch validation",
                   author="Isarar Siddique", affiliation="", email="", plan_sha256=sha(B / "RETROSPECTIVE_ANALYSIS_PLAN.md"),
                   RR23_animals_all26=len(rr_meta), RR23_training_FLT=9, RR23_training_GC=9,
                   RRRM2_animals=19, RRRM2_FLT=10, RRRM2_GC=9, RRRM2_strata=3,
                   RRRM2_exact_label_assignments=len(val_perm), frozen_definitions=len(definitions),
                   RR23_measured_all26=int(rr_measured.sum()), RRRM2_measured_definitions=int(val_measured.sum()),
                   common_measured_events=int(common.sum()), common_transport_fraction=float(common.sum() / rr_measured.sum()),
                   RR23_measured_RRRM2_coverage_fail=int((rr_measured & (val_pair.min(1) < 20)).sum()),
                   RR23_measured_RRRM2_dominance_fail=int((rr_measured & (val_dom.min(1) < .9)).sum()),
                   RR23_measured_RRRM2_both_fail=int((rr_measured & (val_pair.min(1) < 20) & (val_dom.min(1) < .9)).sum()),
                   spearman_rho=rho(d_rr[common], d_val),
                   spearman_animal_bootstrap_CI95=np.quantile(rho_draws, [.025, .975]).tolist(),
                   pearson_r=float(np.corrcoef(d_rr[common], d_val)[0, 1]),
                   same_direction_fraction=float((d_rr[common] * d_val > 0).mean()),
                   same_direction_animal_bootstrap_CI95=np.quantile(concordance_draws, [.025, .975]).tolist(),
                   learned_direction_score_delta=score_delta,
                   learned_direction_score_delta_percentage_points=100 * score_delta,
                   learned_direction_score_animal_bootstrap_CI95=np.quantile(score_boot, [.025, .975]).tolist(),
                   learned_direction_score_exact_one_sided_p=p_up,
                   learned_direction_score_exact_two_sided_p=p_two,
                   primary_positive_transfer_supported=p_up < .05,
                   full_common_panel_Holm_directional05=int((effects.RRRM2_directional_Holm_full_common_panel < .05).sum()),
                   minimum_event_directional_p=float(event_p.min()),
                   bootstrap_draws=NBOOT, bootstrap_seed=SEED,
                   bootstrap_unit="biological animal within exposure and route-age cells",
                   no_RR23_significance_filter=True, no_gene_independence_inference=True,
                   RR23_flight_collection_day_aliased=True, causal_flight_effect_identified=False)
    manifest["sha256_after"] = [sha(p) for p in paths]
    manifest["source_unchanged"] = manifest.sha256_before.eq(manifest.sha256_after)
    assert manifest.source_unchanged.all()
    manifest.to_csv(T / "input_manifest_sha256.csv", index=False)
    (T / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
