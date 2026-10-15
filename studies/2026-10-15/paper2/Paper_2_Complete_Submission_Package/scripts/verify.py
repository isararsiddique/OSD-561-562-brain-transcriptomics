"""Independent arithmetic and inferential-unit checks for the beta extension."""
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

import analyse as a


def main():
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(dict(check=name, status="PASS"))
    a.verify_inputs()
    check("all frozen input sizes and SHA256 hashes match", True)
    meta = pd.read_csv(a.INPUTS / "animal_metadata.csv").set_index("sample", drop=False)
    count = pd.read_csv(a.INPUTS / "all_animal_junction_counts.csv.gz", index_col=0)
    manifest = pd.read_csv(a.INPUTS / "junction_library_manifest.csv")
    defs = pd.read_csv(a.INPUTS / "primary_all26_event_definitions.csv.gz", dtype={"chr": str}).set_index("event_id", drop=False)
    gates = pd.read_csv(a.RESULTS / "measurement_gates.csv.gz").set_index("event_id").loc[defs.index]
    fraction = pd.read_csv(a.RESULTS / "all26_animal_fractions.csv.gz", index_col=0).loc[defs.index]
    events = pd.read_csv(a.RESULTS / "three_arm_event_results.csv.gz")
    glob = pd.read_csv(a.RESULTS / "global_contrasts.csv").set_index("contrast")
    summary = json.loads((a.RESULTS / "summary.json").read_text())
    check("26 distinct animal sources; 9 flight, 9 habitat, 8 vivarium", len(meta) == 26 and meta["Source Name"].is_unique and meta.group.value_counts().to_dict() == {"FLT": 9, "GC": 9, "VIV": 8})
    check("29 official-MD5-matching libraries collapse to 26 samples", len(manifest) == 29 and manifest["sample"].nunique() == 26 and manifest.official_md5_match.all())
    declared = manifest.groupby("sample").retained_junction_reads.sum().loc[count.columns]
    check("every collapsed animal count total matches its source-library sum", np.array_equal(count.sum(axis=0), declared))
    ca = count.reindex(defs.junction_a).to_numpy()
    cb = count.reindex(defs.junction_b).to_numpy()
    check("all saved fractions equal independent raw-count ratios", np.allclose(fraction.to_numpy(), ca/(ca+cb), atol=1e-14))
    check("pooled measurement gates require all 26 animals", np.all((ca+cb).min(axis=1) >= 20) and np.allclose(defs.pooled_mean_fraction, fraction.mean(axis=1)))
    check("all first junctions sort before their paired alternatives", all(first < second for first, second in zip(defs.junction_a, defs.junction_b)))
    check("both donor and acceptor boundaries follow saved coordinates", all((r.junction_a.split(':')[2 if r.boundary_type == 'start' else 3] == str(r.boundary)) and (r.junction_b.split(':')[2 if r.boundary_type == 'start' else 3] == str(r.boundary)) for r in defs.itertuples()))
    check("event multiplicity retains all definitions in all three contrasts", len(events) == 3*len(defs) and events.groupby("contrast").size().eq(len(defs)).all())
    check("joint BH matches an independent SciPy implementation", np.allclose(events.q_joint_3M_tests, false_discovery_control(events.exact_p.to_numpy())))
    check("joint BY matches an independent SciPy implementation", np.allclose(events.q_BY_joint_3M_tests, false_discovery_control(events.exact_p.to_numpy(), method="by")))
    check("unmeasured events retain p=1 in correction family", events.loc[~events.measured, "exact_p"].eq(1).all())
    measured = gates.measured_primary_all26.to_numpy()
    x = fraction.loc[measured]
    rng = np.random.default_rng(19)
    examples = rng.choice(len(x), 4, replace=False)
    for ga, gb in a.CONTRASTS:
        tag = ga + "_minus_" + gb
        samples = list(meta.index[meta.group == ga]) + list(meta.index[meta.group == gb])
        values = x[samples].to_numpy()
        na, nb = int((meta.group == ga).sum()), int((meta.group == gb).sum())
        choices = np.array(list(itertools.combinations(range(na+nb), na)))
        data = values[examples]
        sums_a = data[:, choices].sum(axis=2)
        sums_all = data.sum(axis=1, keepdims=True)
        null = sums_a/na - (sums_all-sums_a)/nb
        observed = data[:, :na].mean(axis=1)-data[:, na:].mean(axis=1)
        exact_p = np.mean(np.abs(null) >= np.abs(observed[:, None])-1e-12, axis=1)
        saved = events[events.contrast == tag].set_index("event_id").loc[x.index[examples]]
        check(tag + ": four exact event p values match separate subset enumeration", np.allclose(exact_p, saved.exact_p, atol=1e-14))
        # Independently reconstruct aggregate distances from kernel block sums,
        # rather than using the production weight-matrix statistic.
        kernel = values.T @ values / len(values)
        scalar_statistics = []
        for first in choices:
            second = np.setdiff1d(np.arange(na+nb), first)
            scalar_statistics.append(kernel[np.ix_(first, first)].sum()/na**2 + kernel[np.ix_(second, second)].sum()/nb**2 - 2*kernel[np.ix_(first, second)].sum()/(na*nb))
        scalar_statistics = np.array(scalar_statistics)
        obs = float(np.mean((values[:, :na].mean(axis=1)-values[:, na:].mean(axis=1))**2))
        p = float(np.mean(scalar_statistics >= obs-1e-13))
        check(tag + ": exact global p and squared distance match independent scalar calculation", abs(p-glob.loc[tag, "exact_p"]) < 1e-14 and abs(obs-glob.loc[tag, "mean_squared_contrast"]) < 1e-14)
        check(tag + ": exact label space agrees with group sizes", int(glob.loc[tag, "allocations"]) == math.comb(na+nb, na))
    check("global Holm values agree with a separately sorted calculation", np.allclose(glob.Holm_three_global_tests, a.holm(glob.exact_p)))
    dates = meta.groupby("group")["Comment[Euthanasia Date]"].agg(set)
    check("flight date is aliased while both ground groups share a date", dates["GC"] == dates["VIV"] and dates["FLT"].isdisjoint(dates["GC"]))
    split = pd.read_csv(a.RESULTS / "disjoint_control_profile_correlations.csv")
    check("126 control partitions contain disjoint 4/5 biological samples", len(split) == math.comb(9,4) and all(len(set(r.control_arm_A.split(';'))) == 4 and len(set(r.control_arm_B.split(';'))) == 5 and set(r.control_arm_A.split(';')).isdisjoint(r.control_arm_B.split(';')) for r in split.itertuples()))
    thinning = pd.read_csv(a.RESULTS / "depth_thinning_global_sensitivity.csv")
    check("all 10 thinning replicates report every global comparison", len(thinning) == 30 and thinning.groupby("replicate").size().eq(3).all())
    qc = pd.read_csv(a.RESULTS / "animal_qc.csv", index_col=0)
    check("animal QC depth joins are aligned by sample ID", qc.retained_junction_depth.equals(count.sum(axis=0).loc[qc.index]))
    check("summary matches observed family and global discovery counts", summary["feature_definitions"] == len(defs) and summary["event_family"] == len(events) and summary["primary_global_significant_Holm05"] == int((glob.Holm_three_global_tests < .05).sum()))
    # Exact finite null: labels are exchangeable, including features correlated
    # across the same six animals; the test's rejection probability is bounded.
    toy = np.array([[0,1,2,1,3,4],[1,2,3,0,1,2],[2,0,1,4,2,3]],float)
    weights = a.binary_weights(6,3)
    kernel = toy.T @ toy / len(toy)
    stats = np.einsum('ij,jk,ik->i',weights,kernel,weights)
    null_p = np.mean(stats[None,:] >= stats[:,None]-1e-12,axis=1)
    check("correlated-feature exact null controls nominal rejection frequency", np.mean(null_p <= .05) <= .05)
    result = dict(status="PASS", checks_passed=len(checks), checks=checks, limitation="Arithmetic, source provenance and randomization checks do not establish causality, functional injury, full isoform identity or cross-launch replication.")
    (a.RESULTS / "verification.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({k:v for k,v in result.items() if k != "checks"}, indent=2))


if __name__ == "__main__":
    main()
