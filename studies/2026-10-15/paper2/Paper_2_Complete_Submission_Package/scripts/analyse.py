"""Exact animal-level three-arm characterization of frozen RR23 junctions."""
from pathlib import Path
from collections import Counter
import hashlib
import itertools
import json
import math

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t

BASE = Path(__file__).resolve().parents[1]
INPUTS, RESULTS = BASE / "inputs", BASE / "results"
CONTRASTS = [("FLT", "GC"), ("GC", "VIV"), ("FLT", "VIV")]
SEED = 20261014


def bh(p, by=False):
    p = np.asarray(p, float)
    order = np.argsort(p)
    factor = np.sum(1 / np.arange(1, len(p)+1)) if by else 1
    sorted_q = np.minimum.accumulate((p[order] * len(p) * factor / np.arange(1, len(p)+1))[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.minimum(sorted_q, 1)
    return q


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.empty_like(p)
    q[order] = np.minimum(np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)), 1)
    return q


def verify_inputs():
    for row in pd.read_csv(INPUTS / "checksums.csv").to_dict("records"):
        path = INPUTS / row["file"]
        assert path.stat().st_size == row["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], path


def binary_weights(n, na):
    weights = np.full((math.comb(n, na), n), -1/(n-na), dtype=float)
    for k, chosen in enumerate(itertools.combinations(range(n), na)):
        weights[k, list(chosen)] = 1/na
    return weights


def fractions(counts, definitions):
    ca = counts.reindex(definitions.junction_a).fillna(0).to_numpy(dtype=float)
    cb = counts.reindex(definitions.junction_b).fillna(0).to_numpy(dtype=float)
    total = ca + cb
    frac = np.divide(ca, total, out=np.full_like(total, np.nan), where=total > 0)
    return frac, total


def dominance(counts, definitions, total):
    ids = counts.index.str.split(":", expand=True)
    info = pd.DataFrame(dict(chr=ids.get_level_values(0), strand=ids.get_level_values(1), start=ids.get_level_values(2), end=ids.get_level_values(3)), index=counts.index)
    denominator = np.zeros_like(total)
    for boundary in ["start", "end"]:
        summed = counts.groupby(info.chr + ":" + info.strand + ":" + info[boundary]).sum()
        ix = np.flatnonzero(definitions.boundary_type.to_numpy() == boundary)
        sub = definitions.iloc[ix]
        keys = sub.chr.astype(str) + ":" + sub.strand.astype(str) + ":" + sub.boundary.astype(str)
        denominator[ix] = summed.reindex(keys).fillna(0).to_numpy()
    return np.divide(total, denominator, out=np.zeros_like(total), where=denominator > 0)


def global_test(x, weights, na, gene_weights=None):
    """Exact label test of mean squared difference; animals are the units."""
    assert np.isfinite(x).all() and x.shape[1] == weights.shape[1]
    feature_weights = np.full(len(x), 1/len(x)) if gene_weights is None else gene_weights / np.sum(gene_weights)
    x = x - x.mean(axis=1, keepdims=True)
    kernel = (x * feature_weights[:, None]).T @ x
    observed_weights = np.r_[np.full(na, 1/na), np.full(x.shape[1]-na, -1/(x.shape[1]-na))]
    observed = float(observed_weights @ kernel @ observed_weights)
    null = np.einsum("ij,jk,ik->i", weights, kernel, weights, optimize=True)
    p = float(np.mean(null >= observed - 1e-14))
    ka, kb, kab = kernel[:na, :na], kernel[na:, na:], kernel[:na, na:]
    nb = len(kb)
    unbiased = float((ka.sum()-np.trace(ka))/(na*(na-1)) + (kb.sum()-np.trace(kb))/(nb*(nb-1)) - 2*kab.mean())
    return dict(features=len(x), mean_squared_contrast=observed, RMS_contrast_pp=100*np.sqrt(max(observed, 0)), unbiased_squared_distance=unbiased, exact_p=p, allocations=len(weights))


def event_test(x, weights, na):
    d = x[:, :na].mean(axis=1) - x[:, na:].mean(axis=1)
    p = np.ones(len(x))
    for start in range(0, len(x), 96):
        null = x[start:start+96] @ weights.T
        p[start:start+96] = np.mean(np.abs(null) >= np.abs(d[start:start+96, None]) - 1e-12, axis=1)
    va = x[:, :na].var(axis=1, ddof=1)/na
    vb = x[:, na:].var(axis=1, ddof=1)/(x.shape[1]-na)
    se = np.sqrt(va+vb)
    df = np.divide((va+vb)**2, va**2/(na-1) + vb**2/(x.shape[1]-na-1), out=np.full(len(x), np.inf), where=(va+vb) > 0)
    return d, p, se, df


def rho(x, y):
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def main():
    verify_inputs()
    definitions = pd.read_csv(INPUTS / "primary_all26_event_definitions.csv.gz", dtype={"chr": str}).set_index("event_id", drop=False)
    counts = pd.read_csv(INPUTS / "all_animal_junction_counts.csv.gz", index_col=0)
    meta = pd.read_csv(INPUTS / "animal_metadata.csv").set_index("sample", drop=False).loc[counts.columns]
    assert len(meta) == 26 and meta["Source Name"].is_unique and len(definitions) > 1000
    assert meta.group.value_counts().to_dict() == {"FLT": 9, "GC": 9, "VIV": 8}
    frac, total = fractions(counts, definitions)
    dom = dominance(counts, definitions, total)
    measured = (total.min(axis=1) >= 20) & (dom.min(axis=1) >= .9)
    measured95 = (total.min(axis=1) >= 20) & (dom.min(axis=1) >= .95)
    assert measured.sum() > 1000
    gates = definitions[["event_id", "gene_id", "gene_symbol"]].copy()
    gates["minimum_pair_reads_all26"] = total.min(axis=1)
    gates["minimum_dominance_all26"] = dom.min(axis=1)
    gates["measured_primary_all26"] = measured
    gates["measured_dominance95"] = measured95
    gates.to_csv(RESULTS / "measurement_gates.csv.gz", index=False)
    pd.DataFrame(frac, index=definitions.index, columns=counts.columns).to_csv(RESULTS / "all26_animal_fractions.csv.gz")
    indices = {g: np.flatnonzero(meta.group.to_numpy() == g) for g in ["FLT", "GC", "VIV"]}
    global_rows, event_rows, cache, sensitivity_rows = [], [], {}, []
    for group_a, group_b in CONTRASTS:
        contrast = group_a + "_minus_" + group_b
        ix = np.r_[indices[group_a], indices[group_b]]
        na = len(indices[group_a])
        weights = binary_weights(len(ix), na)
        cache[contrast] = (ix, na, weights)
        data = frac[measured][:, ix]
        glob = global_test(data, weights, na)
        glob.update(contrast=contrast, n_group_a=na, n_group_b=len(ix)-na)
        global_rows.append(glob)
        delta, p, se, df = event_test(data, weights, na)
        table = definitions[["event_id", "gene_id", "gene_symbol"]].copy()
        table["contrast"] = contrast
        table["measured"] = measured
        table["delta_fraction"] = np.nan
        table["exact_p"] = 1.
        table["SE"] = np.nan
        table["Welch_df"] = np.nan
        table.loc[measured, ["delta_fraction", "exact_p", "SE", "Welch_df"]] = np.column_stack([delta, p, se, df])
        table["q_within_contrast_M"] = bh(table.exact_p)
        table["q_BY_within_contrast_M"] = bh(table.exact_p, by=True)
        event_rows.append(table)
        if measured95.sum():
            s = global_test(frac[measured95][:, ix], weights, na)
            s.update(contrast=contrast, sensitivity="95% dominance")
            sensitivity_rows.append(s)
        mapped = measured & definitions.gene_id.notna().to_numpy()
        genes = definitions.loc[mapped, "gene_id"]
        frequency = Counter(genes)
        gene_weights = np.array([1/frequency[g] for g in genes])
        s = global_test(frac[mapped][:, ix], weights, na, gene_weights=gene_weights)
        s.update(contrast=contrast, sensitivity="equal weight per uniquely mapped gene", mapped_genes=len(frequency))
        sensitivity_rows.append(s)
        print(f"{contrast}: {len(ix)} animals, {len(data)} common events, exact global p {glob['exact_p']:.6g}; event tests complete", flush=True)
    global_results = pd.DataFrame(global_rows)
    global_results["Holm_three_global_tests"] = holm(global_results.exact_p)
    global_results.to_csv(RESULTS / "global_contrasts.csv", index=False)
    events = pd.concat(event_rows, ignore_index=True)
    events["q_joint_3M_tests"] = bh(events.exact_p)
    events["q_BY_joint_3M_tests"] = bh(events.exact_p, by=True)
    events.to_csv(RESULTS / "three_arm_event_results.csv.gz", index=False)
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity["Holm_within_sensitivity_three_tests"] = sensitivity.groupby("sensitivity")["exact_p"].transform(lambda p: holm(p))
    sensitivity.to_csv(RESULTS / "global_measurement_sensitivities.csv", index=False)

    # Composite criterion: difference from BOTH controls + housing equivalence.
    tables = {g: events[events.contrast == g].set_index("event_id").loc[definitions.index] for g in cache}
    fg, gv, fv = [tables[g] for g in ["FLT_minus_GC", "GC_minus_VIV", "FLT_minus_VIV"]]
    eq = np.ones(len(definitions))
    se = gv.SE.to_numpy()
    delta = gv.delta_fraction.to_numpy()
    df = gv.Welch_df.to_numpy()
    usable = measured & (se > 0)
    eq[usable] = np.maximum(t.sf((delta[usable]+.05)/se[usable], df[usable]), t.cdf((delta[usable]-.05)/se[usable], df[usable]))
    zero_se = measured & (se == 0)
    eq[zero_se] = np.where(np.abs(delta[zero_se]) < .05, 0., 1.)
    conjunction = np.maximum.reduce([fg.exact_p.to_numpy(), fv.exact_p.to_numpy(), eq])
    same_direction = fg.delta_fraction.to_numpy() * fv.delta_fraction.to_numpy() > 0
    conjunction[~same_direction] = 1
    classification = definitions[["event_id", "gene_id", "gene_symbol"]].copy()
    classification["measured"] = measured
    classification["FLT_GC_delta"] = fg.delta_fraction.to_numpy()
    classification["GC_VIV_delta"] = gv.delta_fraction.to_numpy()
    classification["FLT_VIV_delta"] = fv.delta_fraction.to_numpy()
    classification["GC_VIV_equivalence_p_5pp"] = eq
    classification["conjunction_p"] = conjunction
    classification["q_conjunction_M"] = bh(conjunction)
    classification["same_direction_against_both_controls"] = same_direction
    classification["classification"] = np.where(classification.q_conjunction_M < .05, "flight-group-distinct within operational housing-equivalence criterion", "criterion not met")
    classification.to_csv(RESULTS / "exploratory_control_classification.csv.gz", index=False)

    # A GC term appears with opposite signs in the two effect profiles. Show
    # how removing shared animals changes the apparent profile relationship.
    use = frac[measured]
    f, g, v = [use[:, indices[k]] for k in ["FLT", "GC", "VIV"]]
    shared_rho = rho(f.mean(axis=1)-g.mean(axis=1), g.mean(axis=1)-v.mean(axis=1))
    split_rows = []
    for first in itertools.combinations(range(9), 4):
        second = sorted(set(range(9))-set(first))
        split_rows.append(dict(control_arm_A=";".join(meta.index[indices['GC'][list(first)]]), control_arm_B=";".join(meta.index[indices['GC'][second]]), rho=rho(f.mean(axis=1)-g[:, list(first)].mean(axis=1), g[:, second].mean(axis=1)-v.mean(axis=1))))
    splits = pd.DataFrame(split_rows)
    splits.to_csv(RESULTS / "disjoint_control_profile_correlations.csv", index=False)
    rng = np.random.default_rng(SEED)
    boot = []
    for _ in range(2000):
        bf = f[:, rng.choice(9, 9, replace=True)].mean(axis=1)
        bg = g[:, rng.choice(9, 9, replace=True)].mean(axis=1)
        bv = v[:, rng.choice(8, 8, replace=True)].mean(axis=1)
        boot.append(rho(bf-bg, bg-bv))
    profile = dict(common_events=int(measured.sum()), shared_control_rho=shared_rho, shared_control_animal_bootstrap_low=float(np.quantile(boot,.025)), shared_control_animal_bootstrap_high=float(np.quantile(boot,.975)), disjoint_control_partitions=len(splits), disjoint_rho_median=float(splits.rho.median()), disjoint_rho_min=float(splits.rho.min()), disjoint_rho_max=float(splits.rho.max()), interpretation="Descriptive profile relationship; shared control introduces negative covariance. No gene-wise binomial or independence p assigned.")
    (RESULTS / "profile_dependence.json").write_text(json.dumps(profile, indent=2))

    # Thin the FULL junction matrix once per replicate so repeated uses of a
    # junction preserve their count dependence. Undefined pairs are reported.
    raw = counts.to_numpy(dtype=np.int64)
    depths = raw.sum(axis=0)
    minimum = int(depths.min())
    probability = minimum / depths
    thinning_rows = []
    for replicate in range(10):
        rng = np.random.default_rng(SEED + 100 + replicate)
        thinned = pd.DataFrame(rng.binomial(raw, probability[None, :]), index=counts.index, columns=counts.columns)
        thin_fraction, thin_total = fractions(thinned, definitions)
        valid = measured & np.isfinite(thin_fraction).all(axis=1)
        for contrast, (ix, na, weights) in cache.items():
            s = global_test(thin_fraction[valid][:, ix], weights, na)
            s.update(replicate=replicate+1, contrast=contrast, target_expected_junction_depth=minimum, undefined_pairs_excluded=int((measured & ~valid).sum()))
            thinning_rows.append(s)
        print(f"Depth sensitivity {replicate+1}/10 complete", flush=True)
    thinning = pd.DataFrame(thinning_rows)
    thinning["Holm_three_global_tests"] = thinning.groupby("replicate")["exact_p"].transform(lambda p: holm(p))
    thinning.to_csv(RESULTS / "depth_thinning_global_sensitivity.csv", index=False)

    # Check preserved primary estimates without conflating historical and new n.
    old = pd.read_csv(INPUTS / "prior_OSD-525_frozen_event_results.csv.gz", index_col=0)
    overlap = definitions.index[measured].intersection(old.index)
    old_delta = pd.Series(f.mean(axis=1)-g.mean(axis=1), index=definitions.index[measured])
    np.testing.assert_allclose(old.loc[overlap, "delta"], old_delta.loc[overlap], atol=1e-12, rtol=0)
    primary_counts = events.groupby("contrast").apply(lambda z: pd.Series(dict(measured=int(z.measured.sum()), q05_joint3M=int((z.q_joint_3M_tests<.05).sum()), q05_joint_and5pp=int(((z.q_joint_3M_tests<.05)&z.delta_fraction.abs().ge(.05)).sum()), q05_BY_joint=int((z.q_BY_joint_3M_tests<.05).sum()), q05_within=int((z.q_within_contrast_M<.05).sum()))), include_groups=False)
    primary_counts.to_csv(RESULTS / "event_discovery_summary.csv")
    fields = ["Source Name", "sample", "group", "Comment[Euthanasia Date]", "Comment[BSP Dissection Date]", "Characteristics[Strain]", "Characteristics[Sex]", "Parameter Value[Sample Preservation Method]"]
    meta[fields].to_csv(RESULTS / "animal_design.csv", index=False)
    qc = pd.read_csv(INPUTS / "library_qc.csv")
    qc["sample"] = qc["sample"].str.replace(r"_techrep[0-9]+$", "_techrep1", regex=True)
    qca = qc.groupby("sample").agg(read_depth=("read_depth", "sum"), rin=("rin", "mean"), rrna_contamination=("rrna_contamination", "mean"), uniquely_mapped_percent=("uniquely_mapped_percent", "mean"))
    qca["group"] = meta.group
    qca["source_animal"] = meta["Source Name"]
    qca["retained_junction_depth"] = pd.Series(depths, index=meta.index)
    qca.to_csv(RESULTS / "animal_qc.csv")
    summary = dict(status="retrospective three-arm extension executed with all26 label-invariant features", animals=26, source_libraries=29, feature_definitions=len(definitions), common_measured_events=int(measured.sum()), dominance95_events=int(measured95.sum()), primary_global_significant_Holm05=int((global_results.Holm_three_global_tests<.05).sum()), event_family=len(events), event_discoveries=primary_counts.to_dict('index'), exploratory_control_criterion_q05=int((classification.q_conjunction_M<.05).sum()), collection_date_confounded_with_flight=True, no_causal_flight_identification=True, prior_FLT_GC_deltas_independently_reproduced=True, prior_definition_overlap_checked=len(overlap), seed=SEED)
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
