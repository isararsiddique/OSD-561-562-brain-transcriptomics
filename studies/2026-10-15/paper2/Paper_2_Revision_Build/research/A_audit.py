#!/usr/bin/env python3
"""A_audit.py -- evidence audit of Paper 2 (cerebellar RNA junction profiles).

Stdlib only (no pandas/numpy needed; host memory is tight and pandas is absent).
Reads EXISTING result files only -- never reruns the pipeline, never writes
inside the submission package.

Two jobs:
  1. RECOMPUTE a set of headline numbers from the raw/complete tables.
  2. CHECK every numeric claim in the manuscript Abstract/Results/Discussion
     against a located source value, writing A_claims_check.csv.

Usage:  OPENBLAS_NUM_THREADS=1 python3 A_audit.py [--section recompute|claims|all]
"""
import csv
import gzip
import io
import json
import os
import sys
from collections import Counter, defaultdict

P = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "Paper_2_Complete_Submission_Package")
OUT = os.path.dirname(os.path.abspath(__file__))

csv.field_size_limit(1 << 24)


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------
def rows(rel):
    """Yield dict rows from a csv or csv.gz under the package root."""
    path = os.path.join(P, rel)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        for r in csv.DictReader(fh):
            yield r


def jload(rel):
    with open(os.path.join(P, rel)) as fh:
        return json.load(fh)


def fnum(x):
    """Float from a result-file cell."""
    if x is None or x == "":
        return None
    return float(x)


def quantile(sorted_vals, q):
    """Linear-interpolation quantile (numpy 'linear' convention)."""
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def median(vals):
    s = sorted(vals)
    return quantile(s, 0.5)


# --------------------------------------------------------------------------
# multiplicity corrections, implemented from the definitions
# --------------------------------------------------------------------------
def holm(pvals):
    """Holm step-down adjusted p values, monotone-enforced."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    out = [0.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        adj = (m - rank) * pvals[i]
        running = max(running, adj)
        out[i] = min(1.0, running)
    return out


def bh(pvals):
    """Benjamini-Hochberg adjusted q values (step-up, monotone-enforced)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    out = [0.0] * m
    running = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        running = min(running, m * pvals[i] / (rank + 1))
        out[i] = min(1.0, running)
    return out


def harmonic(m):
    return sum(1.0 / k for k in range(1, m + 1))


def by(pvals):
    """Benjamini-Yekutieli adjusted q values = BH scaled by H_m."""
    m = len(pvals)
    hm = harmonic(m)
    order = sorted(range(m), key=lambda i: pvals[i])
    out = [0.0] * m
    running = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        running = min(running, m * hm * pvals[i] / (rank + 1))
        out[i] = min(1.0, running)
    return out


# --------------------------------------------------------------------------
# RECOMPUTATIONS
# --------------------------------------------------------------------------
RECOMP = {}


def recompute():
    R = RECOMP
    line = lambda s: print(s)

    print("=" * 78)
    print("PART 2 -- INDEPENDENT RECOMPUTATION FROM RAW TABLES")
    print("=" * 78)

    # ---- (1) BY threshold arithmetic for the 2,769-test multiway family ----
    m = 2769
    hm = harmonic(m)
    p_floor = 2.0 / 48620.0
    R["H_2769"] = hm
    R["p_floor"] = p_floor
    R["inv_24310"] = 1.0 / 24310.0
    R["q_block20"] = m * hm * p_floor / 20
    R["q_block19"] = m * hm * p_floor / 19
    R["min_ties_for_BY05"] = next(k for k in range(1, m + 1)
                                  if m * hm * p_floor / k < 0.05)
    print("\n(1) BY threshold arithmetic, m = 2769 tests")
    print("    harmonic factor H_m            = %.8f" % hm)
    print("    2/48620                        = %.12g" % p_floor)
    print("    1/24310                        = %.12g" % R["inv_24310"])
    print("    equal?                           %s"
          % (abs(p_floor - R["inv_24310"]) < 1e-18))
    print("    q for a 20-test floor block    = %.8f" % R["q_block20"])
    print("    q for a 19-test floor block    = %.8f" % R["q_block19"])
    print("    smallest floor block crossing 0.05 = %d"
          % R["min_ties_for_BY05"])

    # ---- (2) Holm values in Table 1, from the exact p values ----
    print("\n(2) Table 1 Holm correction, recomputed from exact p")
    binary = [(r["contrast"], fnum(r["exact_p"]), fnum(r["RMS_contrast_pp"]),
               fnum(r["Holm_three_global_tests"]))
              for r in rows("results/global_contrasts.csv")]
    bh_adj = holm([b[1] for b in binary])
    R["holm_binary"] = {}
    for (name, p, mag, stored), mine in zip(binary, bh_adj):
        R["holm_binary"][name] = mine
        print("    binary %-14s p=%.8f  Holm mine=%.8f  stored=%.8f  %s"
              % (name, p, mine, stored,
                 "MATCH" if abs(mine - stored) < 1e-12 else "MISMATCH"))
    multi = [(r["contrast"], fnum(r["exact_p"]),
              fnum(r["hellinger_centroid_distance"]), fnum(r["holm_p_3"]))
             for r in rows("extensions/multiway_splicing/tables/"
                           "global_exact_results.csv")]
    mh_adj = holm([b[1] for b in multi])
    R["holm_multi"] = {}
    for (name, p, mag, stored), mine in zip(multi, mh_adj):
        R["holm_multi"][name] = mine
        print("    multi  %-14s p=%.8f  Holm mine=%.8f  stored=%.8f  %s"
              % (name, p, mine, stored,
                 "MATCH" if abs(mine - stored) < 1e-12 else "MISMATCH"))
    # magnitude identity: multialternative magnitude = sqrt(T/2)
    chk = []
    for r in rows("extensions/multiway_splicing/tables/global_exact_results.csv"):
        t = fnum(r["statistic"])
        chk.append((r["contrast"], (t / 2.0) ** 0.5,
                    fnum(r["hellinger_centroid_distance"])))
    R["sqrtT2"] = chk
    print("    sqrt(T/2) identity: " + ", ".join(
        "%s %.6f vs %.6f" % (c, a, b) for c, a, b in chk))

    # ---- (3) BY-positive counts, recomputed over the frozen families ----
    print("\n(3) BY-positive counts recomputed from complete families")
    # gene omnibus: 5,949 joint gene-contrast tests
    g = list(rows("extensions/gene_omnibus/results/gene_omnibus_complete.csv.gz"))
    gp = [fnum(r["p_exact"]) for r in g]
    gq = by(gp)
    gene_pos = Counter()
    rms_pos = []
    for r, q in zip(g, gq):
        if q < 0.05:
            gene_pos[r["contrast"]] += 1
            rms_pos.append(fnum(r["rms_delta_percentage_points"]))
    R["gene_family_size"] = len(g)
    R["gene_BY_by_contrast"] = dict(gene_pos)
    R["gene_BY_total"] = sum(gene_pos.values())
    R["gene_BY_minq"] = min(q for q in gq if q < 0.05)
    R["gene_BY_rms_min"] = min(rms_pos)
    R["gene_BY_rms_max"] = max(rms_pos)
    R["gene_BY_all_at_floor"] = all(
        abs(fnum(r["p_exact"]) - 2.0 / 48620.0) < 1e-15
        for r, q in zip(g, gq) if q < 0.05)
    print("    gene family size                = %d" % len(g))
    print("    BY q<0.05 by contrast           = %s" % dict(gene_pos))
    print("    minimum BY q among positives    = %.9f" % R["gene_BY_minq"])
    print("    all positives at p=2/48620?       %s" % R["gene_BY_all_at_floor"])
    print("    RMS pp range of positives       = %.2f to %.2f"
          % (R["gene_BY_rms_min"], R["gene_BY_rms_max"]))

    # multiway local: 2,769 joint cluster-contrast tests
    c = list(rows("extensions/multiway_splicing/tables/"
                  "local_exact_results_complete.csv.gz"))
    cp = [fnum(r["exact_p"]) for r in c]
    cq = by(cp)
    cl_pos = Counter()
    uniq = set()
    for r, q in zip(c, cq):
        if q < 0.05:
            cl_pos[r["contrast"]] += 1
            uniq.add(r["cluster_id"])
    R["cluster_family_size"] = len(c)
    R["cluster_BY_by_contrast"] = dict(cl_pos)
    R["cluster_BY_total"] = sum(cl_pos.values())
    R["cluster_BY_unique"] = len(uniq)
    R["cluster_n_at_floor"] = sum(
        1 for p in cp if abs(p - 2.0 / 48620.0) < 1e-15)
    R["cluster_BY_minq"] = min(cq)
    print("    cluster family size             = %d" % len(c))
    print("    BY q<0.05 by contrast           = %s" % dict(cl_pos))
    print("    total BY-positive tests         = %d" % R["cluster_BY_total"])
    print("    unique BY clusters              = %d" % R["cluster_BY_unique"])
    print("    tests at the exact p floor      = %d" % R["cluster_n_at_floor"])
    print("    minimum BY q in family          = %.9f" % R["cluster_BY_minq"])

    # ---- (4) 354 split-pair positive alignments ----
    print("\n(4) Complementary half-pair alignment on the original 43 tests")
    pos = Counter()
    meas = Counter()
    total = Counter()
    cos = defaultdict(list)
    proj = defaultdict(list)
    for r in rows("extensions/effect_stability/results/"
                  "all_complementary_half_concordance.csv.gz"):
        if r["primary_BY_test"] != "True":
            continue
        k = r["contrast"]
        total[k] += 1
        if r["both_halves_measured"] == "True":
            meas[k] += 1
            cv = fnum(r["cosine"])
            cos[k].append(cv)
            proj[k].append(fnum(r["relative_projection"]))
            if cv is not None and cv > 0:
                pos[k] += 1
    R["split_total"] = sum(total.values())
    R["split_measured"] = sum(meas.values())
    R["split_measured_by"] = dict(meas)
    R["split_positive"] = sum(pos.values())
    R["split_median_cos"] = {k: median(v) for k, v in cos.items()}
    R["split_median_proj"] = {k: median(v) for k, v in proj.items()}
    print("    possible instances              = %d" % R["split_total"])
    print("    jointly measured instances      = %d  %s"
          % (R["split_measured"], dict(meas)))
    print("    with positive dot product       = %d" % R["split_positive"])
    for k in sorted(cos):
        print("    %-12s median cosine %.6f   median rel. projection %.6f"
              % (k, R["split_median_cos"][k], R["split_median_proj"][k]))

    # ---- (5) RR10 / RRRM2 measurable clusters ----
    print("\n(5) External measurable clusters out of the frozen 923")
    for coh in ("RR10", "RRRM2"):
        f = ("extensions/comparator_transport/results/%s_all923_"
             "measurement_gates.csv" % coh)
        rr = list(rows(f))
        n_meas = sum(1 for r in rr if r["common_measurable"] == "True")
        R["%s_total" % coh] = len(rr)
        R["%s_measurable" % coh] = n_meas
        R["%s_fraction" % coh] = n_meas / len(rr)
        print("    %-6s %d of %d measurable (%.4f%%)"
              % (coh, n_meas, len(rr), 100 * n_meas / len(rr)))

    # ---- (6) B-half BY calls ----
    print("\n(6) Complementary B-half local BY calls")
    per_rep = Counter()
    for r in rows("extensions/effect_stability/results/"
                  "complement_B_local_exact_full_families.csv.gz"):
        if r["BY_positive"] == "True":
            per_rep[int(r["replicate"])] += 1
    counts = [per_rep.get(i, 0) for i in range(1, 11)]
    R["B_half_counts"] = counts
    R["B_half_total"] = sum(counts)
    R["B_half_with_calls"] = sum(1 for x in counts if x > 0)
    print("    per-replicate B-half calls      = %s" % counts)
    print("    total call instances            = %d" % R["B_half_total"])
    print("    B halves with >=1 call          = %d of 10"
          % R["B_half_with_calls"])
    envel = list(rows("extensions/effect_stability/results/"
                      "BY_resolution_envelopes.csv"))
    a_floor = [int(r["n_at_minimum_p"]) for r in envel
               if r["family"].endswith("_A")]
    b_floor = [int(r["n_at_minimum_p"]) for r in envel
               if r["family"].endswith("_B")]
    R["A_floor_counts"] = a_floor
    R["B_floor_counts"] = b_floor
    R["A_max_floor"] = max(a_floor)
    print("    A-half floor-level test counts  = %s (max %d)"
          % (a_floor, max(a_floor)))
    print("    B-half floor-level test counts  = %s" % b_floor)

    # ---- (7) Tcaf1 native fractions ----
    print("\n(7) Tcaf1 longest-intron fractions, recomputed per animal")
    byg = defaultdict(list)
    for r in rows("extensions/transcript_context/tables/"
                  "Tcaf1_all26_alternative_read_support.csv"):
        if int(r["intron_length_bp"]) == 23000:
            byg[r["group"]].append(fnum(r["proportion"]))
    means = {g: sum(v) / len(v) for g, v in byg.items()}
    R["Tcaf1_means"] = means
    R["Tcaf1_FLT_GC_pp"] = 100 * (means["FLT"] - means["GC"])
    R["Tcaf1_FLT_VIV_pp"] = 100 * (means["FLT"] - means["VIV"])
    for g in ("FLT", "GC", "VIV"):
        print("    mean %s long-intron fraction   = %.4f%% (n=%d)"
              % (g, 100 * means[g], len(byg[g])))
    print("    FLT-GC  difference              = %.4f pp" % R["Tcaf1_FLT_GC_pp"])
    print("    FLT-VIV difference              = %.4f pp" % R["Tcaf1_FLT_VIV_pp"])

    # ---- (8) DM overdispersion and bound hits ----
    print("\n(8) Dirichlet-multinomial fit diagnostics")
    rho_full = []
    bound_full = 0
    bound_half = 0
    gap_half = 0
    nfit = Counter()
    conv = Counter()
    for r in rows("extensions/count_model_benchmark/results/"
                  "all_model_fit_diagnostics.csv.gz"):
        rep = r["replicate"]
        nfit[rep] += 1
        if r["fit_valid"] == "True":
            conv[rep] += 1
        hit = r["alpha_bound_hit"] == "True"
        if rep == "full":
            rho_full.append(fnum(r["rho"]))
            bound_full += hit
        else:
            bound_half += hit
            if fnum(r["best_solution_loglik_gap"]) > 1e-3:
                gap_half += 1
    s = sorted(rho_full)
    R["rho_median"] = quantile(s, 0.5)
    R["rho_q25"] = quantile(s, 0.25)
    R["rho_q75"] = quantile(s, 0.75)
    R["full_fits"] = nfit["full"]
    R["full_converged"] = conv["full"]
    R["bound_full"] = bound_full
    R["bound_half"] = bound_half
    R["gap_half"] = gap_half
    print("    full-depth fits / converged     = %d / %d"
          % (nfit["full"], conv["full"]))
    print("    median rho                      = %.8f" % R["rho_median"])
    print("    rho IQR                         = %.8f to %.8f"
          % (R["rho_q25"], R["rho_q75"]))
    print("    full-depth alpha-bound hits     = %d" % bound_full)
    print("    A-half alpha-bound hits         = %d" % bound_half)
    print("    A-half fits w/ start-loglik gap = %d" % gap_half)

    # ---- (9) Holm attainability floors in the external cohorts ----
    print("\n(9) External Holm resolution floors")
    R["holm_floor_923"] = 923 / 11200
    R["holm_floor_1512"] = 1512 / 11200
    print("    923 / 11200                     = %.8f" % R["holm_floor_923"])
    print("    1512 / 11200                    = %.8f" % R["holm_floor_1512"])

    # ---- (10) design / QC descriptives ----
    print("\n(10) Design and QC descriptives")
    qc = list(rows("results/animal_qc.csv"))
    depth = [fnum(r["retained_junction_depth"]) for r in qc]
    rin = [fnum(r["rin"]) for r in qc]
    grp = Counter(r["group"] for r in qc)
    R["n_animals"] = len(qc)
    R["depth_min_M"] = min(depth) / 1e6
    R["depth_max_M"] = max(depth) / 1e6
    R["rin_min"] = min(rin)
    R["rin_max"] = max(rin)
    R["group_counts"] = dict(grp)
    man = list(rows("inputs/junction_library_manifest.csv"))
    R["n_libraries"] = len(man)
    R["n_md5_match"] = sum(1 for r in man if r["official_md5_match"] == "True")
    print("    animals / libraries             = %d / %d"
          % (len(qc), len(man)))
    print("    official MD5 matches            = %d" % R["n_md5_match"])
    print("    group counts                    = %s" % dict(grp))
    print("    retained depth range (M reads)  = %.2f to %.2f"
          % (R["depth_min_M"], R["depth_max_M"]))
    print("    RIN range                       = %.1f to %.1f"
          % (R["rin_min"], R["rin_max"]))
    R["total_animals_3_cohorts"] = (
        len(qc)
        + jload("extensions/comparator_transport/results/summary.json")
        ["RRRM2_animals"]
        + jload("extensions/comparator_transport/results/summary.json")
        ["RR10_animals"])
    print("    RR23 + RRRM2 + RR10 animals     = %d"
          % R["total_animals_3_cohorts"])
    return R


# --------------------------------------------------------------------------
# claim resolvers
# --------------------------------------------------------------------------
_CACHE = {}


def _tab(rel):
    if rel not in _CACHE:
        _CACHE[rel] = list(rows(rel))
    return _CACHE[rel]


def j(rel, *keys):
    """Nested JSON value."""
    v = jload(rel)
    for k in keys:
        v = v[k]
    return v


def jl(rel, listkey, key, val, field):
    """Field of the first dict in a JSON list whose key == val."""
    for d in jload(rel)[listkey]:
        if d[key] == val:
            return d[field]
    return None


def cv(rel, field, **filters):
    """First matching cell of a csv/csv.gz."""
    for r in _tab(rel):
        if all(r.get(k) == v for k, v in filters.items()):
            return fnum(r[field]) if _looks_num(r[field]) else r[field]
    return None


def cs(rel, field, **filters):
    """First matching cell, as a raw string."""
    for r in _tab(rel):
        if all(r.get(k) == v for k, v in filters.items()):
            return r[field]
    return None


def ca(rel, field, how, **filters):
    """Aggregate over matching rows: min/max/sum/count/nuniq."""
    vals = []
    for r in _tab(rel):
        if all(r.get(k) == v for k, v in filters.items()):
            vals.append(r[field] if how == "nuniq" else fnum(r[field]))
    if how == "count":
        return len(vals)
    if how == "nuniq":
        return len(set(vals))
    if not vals:
        return None
    return {"min": min, "max": max, "sum": sum}[how](vals)


def cn(rel, **filters):
    """Count matching rows."""
    return sum(1 for r in _tab(rel)
               if all(r.get(k) == v for k, v in filters.items()))


def _looks_num(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def R(key):
    """A value produced by the independent recomputation."""
    return RECOMP[key]


# short aliases for the heavily used files
F_GLOBAL = "results/global_contrasts.csv"
F_MULTI_G = "extensions/multiway_splicing/tables/global_exact_results.csv"
F_SEL = "extensions/multiway_splicing/tables/selection_summary.json"
F_INF = "extensions/multiway_splicing/tables/inference_summary.json"
F_SUM = "results/summary.json"
F_GENE = "extensions/gene_omnibus/results/discovery_summary.csv"
F_GENEM = "extensions/gene_omnibus/results/manifest.json"
F_GLOO = "extensions/gene_omnibus/results/descriptive_leave_one_out_summary.json"
F_DOM95 = "extensions/gene_omnibus/results/dominance95_sensitivity_summary.csv"
F_ENV = "extensions/effect_stability/results/BY_resolution_envelopes.csv"
F_ES = "extensions/effect_stability/results/effect_stability_summary.csv"
F_ESS = "extensions/effect_stability/results/summary.json"
F_ROB = "extensions/multiway_robustness/results/summary.json"
F_ROBR = "extensions/multiway_robustness/results/original_hit_retention_by_contrast.csv"
F_ROBS = "extensions/multiway_robustness/results/original_hit_stability_summary.csv"
F_BOOT = "extensions/multiway_robustness/results/global_animal_bootstrap_intervals.csv"
F_DISP = "extensions/multiway_robustness/results/dispersion_and_studentization_diagnostics.csv"
F_CM = "extensions/count_model_benchmark/results/all_global_results.csv"
F_CMC = "extensions/count_model_benchmark/results/model_concordance_summary.csv"
F_CMA = "extensions/count_model_benchmark/results/exact_resolution_and_BY_threshold_audit.csv"
F_CMT = "extensions/count_model_benchmark/results/Tcaf1_full_model_comparison.csv"
F_CMF = "extensions/count_model_benchmark/results/fit_diagnostic_summary.csv"
F_CMS = "extensions/count_model_benchmark/results/summary.json"
F_QC = "extensions/confounding_sensitivity/results/QC_adjustment_sensitivity.csv"
F_TIP = "extensions/confounding_sensitivity/results/bias_tipping_points.csv"
F_GEO = "extensions/confounding_sensitivity/results/comparator_geometry.csv"
F_NUL = "extensions/confounding_sensitivity/results/identification_nullspace.csv"
F_CFS = "extensions/confounding_sensitivity/results/summary.json"
F_CT = "extensions/comparator_transport/results/all_comparator_global_transport.csv"
F_CTS = "extensions/comparator_transport/results/summary.json"
F_CTT = "extensions/comparator_transport/results/Tcaf1_conditional_intervals.csv"
F_IV = "extensions/independent_validation/tables/summary.json"
F_MT = "extensions/multiway_transport/tables/summary.json"
F_AB = "extensions/abundance_integration/results/summary.json"
F_TC = "extensions/transcript_context/tables/summary.json"
F_TCS = "extensions/transcript_context/tables/Tcaf1_transcript_CDS_UTR_summary.csv"
F_TCJ = "extensions/transcript_context/tables/Tcaf1_measured_junction_transcripts.csv"
F_FV = "extensions/multiway_splicing/tables/three_flight_vivarium_BY_clusters.csv"
F_FUN = "extensions/multiway_splicing/tables/selection_funnel.csv"
F_PD = "results/profile_dependence.json"
F_AQ = "results/animal_qc.csv"
F_AD = "results/animal_design.csv"
F_LM = "inputs/junction_library_manifest.csv"
F_EDS = "results/event_discovery_summary.csv"


if __name__ == "__main__":
    sec = "all"
    if "--section" in sys.argv:
        sec = sys.argv[sys.argv.index("--section") + 1]
    if sec in ("recompute", "all"):
        recompute()
    if sec in ("claims", "all"):
        if not RECOMP:
            recompute()
        from A_claims import run_claims
        run_claims(globals())
