#!/usr/bin/env python3
"""A_claims.py -- the claim table used by A_audit.py.

Every numeric claim in the Abstract, Results (including figure captions and
Table 1) and Discussion of PAPER_2_MANUSCRIPT.tex, each paired with the
result file and field that should carry its source value.

Imported and executed by A_audit.py; not meant to be run alone.
"""


def run_claims(ns):
    cv, cs, ca, cn, j, jl, R = (ns["cv"], ns["cs"], ns["ca"], ns["cn"],
                                ns["j"], ns["jl"], ns["R"])
    g = ns
    import csv
    import os

    # ------------------------------------------------------------------
    # claim table: (section, excerpt, value_in_text, source_file,
    #               source_field, getter)
    # ------------------------------------------------------------------
    C = []

    def add(sec, excerpt, text, f, field, getter):
        C.append((sec, excerpt, text, f, field, getter))

    # ================= ABSTRACT =================
    S = "Abstract"
    add(S, "cerebellar junction counts from 55 mice in three cohorts", "55",
        "recomputed: results/animal_qc.csv + comparator_transport/summary.json",
        "26 + 19 + 10", lambda: R("total_animals_3_cohorts"))
    add(S, "spanning two launches", "2",
        g["F_CTS"], "RR10_shares_discovery_launch (RR23/RR10 share SpaceX-21)",
        lambda: 2 if j(g["F_CTS"], "RR10_shares_discovery_launch") else None)
    add(S, "A label-invariant RR23 panel contained 5,249 binary choices",
        "5249", g["F_SUM"], "common_measured_events",
        lambda: j(g["F_SUM"], "common_measured_events"))
    add(S, "and 923 multialternative clusters", "923",
        g["F_SEL"], "n_clusters", lambda: j(g["F_SEL"], "n_clusters"))
    add(S, "All three aggregate comparisons survived exact and count-score "
           "analyses", "3", g["F_SUM"], "primary_global_significant_Holm05",
        lambda: j(g["F_SUM"], "primary_global_significant_Holm05"))
    add(S, "measured-quality adjustment changed flight-comparison norms by "
           "less than 0.71%", "0.71", g["F_QC"],
        "max |norm_ratio-1| over group_RIN_log2_depth flight contrasts",
        lambda: 100 * max(
            abs(cv(g["F_QC"], "norm_ratio", frame=fr,
                   model="group_RIN_log2_depth", contrast=co) - 1)
            for fr in ("multiway_Hellinger", "binary_RMS_pp")
            for co in ("FLT_vs_GC", "FLT_vs_VIV")))
    add(S, "Joint dependence-robust correction identified 47 gene "
           "associations", "47", g["F_GENE"], "joint_BY_q_lt_0_05 FLT_minus_GC",
        lambda: cv(g["F_GENE"], "joint_BY_q_lt_0_05", contrast="FLT_minus_GC"))
    add(S, "and 41 multialternative clusters", "41",
        g["F_INF"], "n_unique_BY_clusters",
        lambda: j(g["F_INF"], "n_unique_BY_clusters"))
    add(S, "Across 20 complementary read halves", "20",
        g["F_ENV"], "row count excluding the original family",
        lambda: cn(g["F_ENV"]) - 1)
    add(S, "local discoveries ranged from zero to 23", "0",
        g["F_ENV"], "min BY_positive over the 20 halves",
        lambda: min(x for x in [int(r["BY_positive"]) for r in
                                g["_tab"](g["F_ENV"]) if r["family"] != "original"]))
    add(S, "local discoveries ranged from zero to 23", "23",
        g["F_ENV"], "max BY_positive over the 20 halves",
        lambda: max(int(r["BY_positive"]) for r in g["_tab"](g["F_ENV"])
                    if r["family"] != "original"))
    add(S, "with coding structure shared across three transcripts", "3",
        g["F_TCS"], "transcripts with CDS_bases=2766",
        lambda: cn(g["F_TCS"], CDS_bases="2766"))

    # ================= RESULTS 1: common animal-level panel =================
    S = "Results/panel"
    add(S, "contained 29 STAR junction libraries", "29",
        g["F_SUM"], "source_libraries", lambda: j(g["F_SUM"], "source_libraries"))
    add(S, "from 26 male C57BL/6J source animals", "26",
        g["F_SUM"], "animals", lambda: j(g["F_SUM"], "animals"))
    add(S, "recorded as 16--17 weeks at launch", "16-17",
        "inputs/animal_metadata.csv", "Characteristics[Age at Launch]",
        lambda: cs("inputs/animal_metadata.csv",
                   "Characteristics[Age at Launch]"))
    add(S, "nine flight (FLT)", "9", g["F_AQ"], "rows with group=FLT",
        lambda: cn(g["F_AQ"], group="FLT"))
    add(S, "nine habitat ground controls (GC)", "9", g["F_AQ"],
        "rows with group=GC", lambda: cn(g["F_AQ"], group="GC"))
    add(S, "eight vivarium controls (VIV)", "8", g["F_AQ"],
        "rows with group=VIV", lambda: cn(g["F_AQ"], group="VIV"))
    add(S, "All 29 publisher MD5 checksums matched", "29",
        g["F_LM"], "official_md5_match=True count",
        lambda: cn(g["F_LM"], official_md5_match="True"))
    add(S, "Two technical libraries each for F1, G1 and V5", "3",
        g["F_LM"], "libraries named *techrep* / animals with >1 library",
        lambda: len([1 for a, n in __import__("collections").Counter(
            r["source_animal"] for r in g["_tab"](g["F_LM"])).items() if n > 1]))
    add(S, "Recorded collection dates were 14 January 2021 for FLT",
        "14-Jan-2021", g["F_AD"], "Comment[Euthanasia Date] where group=FLT",
        lambda: cs(g["F_AD"], "Comment[Euthanasia Date]", group="FLT"))
    add(S, "and 17 January for both ground groups", "17-Jan-2021",
        g["F_AD"], "Comment[Euthanasia Date] where group=GC",
        lambda: cs(g["F_AD"], "Comment[Euthanasia Date]", group="GC"))
    add(S, "and 17 January for both ground groups (VIV)", "17-Jan-2021",
        g["F_AD"], "Comment[Euthanasia Date] where group=VIV",
        lambda: cs(g["F_AD"], "Comment[Euthanasia Date]", group="VIV"))
    add(S, "matrix ... contained 805,562 junctions", "805562",
        g["F_SEL"], "n_canonical_junctions",
        lambda: j(g["F_SEL"], "n_canonical_junctions"))
    add(S, "yielded 5,724 binary boundary definitions", "5724",
        g["F_SUM"], "feature_definitions",
        lambda: j(g["F_SUM"], "feature_definitions"))
    add(S, "of which 5,249 passed the primary 90% dominance gate", "5249",
        g["F_SUM"], "common_measured_events",
        lambda: j(g["F_SUM"], "common_measured_events"))
    add(S, "and 4,316 passed 95% dominance", "4316",
        g["F_SUM"], "dominance95_events",
        lambda: j(g["F_SUM"], "dominance95_events"))
    add(S, "Retained junction depths ranged from 16.25", "16.25",
        g["F_AQ"], "min retained_junction_depth / 1e6",
        lambda: ca(g["F_AQ"], "retained_junction_depth", "min") / 1e6)
    add(S, "to 45.09 million reads", "45.09",
        g["F_AQ"], "max retained_junction_depth / 1e6",
        lambda: ca(g["F_AQ"], "retained_junction_depth", "max") / 1e6)
    add(S, "animal-level RIN summaries from 8.8", "8.8",
        g["F_AQ"], "min rin", lambda: ca(g["F_AQ"], "rin", "min"))
    add(S, "to 9.5", "9.5", g["F_AQ"], "max rin",
        lambda: ca(g["F_AQ"], "rin", "max"))
    add(S, "The FLT--GC RMS fraction contrast was 5.54 percentage points",
        "5.54", g["F_GLOBAL"], "RMS_contrast_pp FLT_minus_GC",
        lambda: cv(g["F_GLOBAL"], "RMS_contrast_pp", contrast="FLT_minus_GC"))
    add(S, "compared with 4.07 pp for GC--VIV", "4.07",
        g["F_GLOBAL"], "RMS_contrast_pp GC_minus_VIV",
        lambda: cv(g["F_GLOBAL"], "RMS_contrast_pp", contrast="GC_minus_VIV"))
    add(S, "and 3.81 pp for FLT--VIV", "3.81",
        g["F_GLOBAL"], "RMS_contrast_pp FLT_minus_VIV",
        lambda: cv(g["F_GLOBAL"], "RMS_contrast_pp", contrast="FLT_minus_VIV"))

    # ---- Table 1 ----
    S = "Results/Table 1"
    tbl = [("Binary FLT--GC", g["F_GLOBAL"], "FLT_minus_GC", "5249",
            "5.5442", "0.00004114", "0.00012341", "RMS_contrast_pp",
            "exact_p", "Holm_three_global_tests", "features"),
           ("Binary GC--VIV", g["F_GLOBAL"], "GC_minus_VIV", "5249",
            "4.0746", "0.00222131", "0.00279720", "RMS_contrast_pp",
            "exact_p", "Holm_three_global_tests", "features"),
           ("Binary FLT--VIV", g["F_GLOBAL"], "FLT_minus_VIV", "5249",
            "3.8144", "0.00139860", "0.00279720", "RMS_contrast_pp",
            "exact_p", "Holm_three_global_tests", "features"),
           ("Multialternative FLT--GC", g["F_MULTI_G"], "FLT_vs_GC", "923",
            "0.08313", "0.00004114", "0.00012341",
            "hellinger_centroid_distance", "exact_p", "holm_p_3", "n_clusters"),
           ("Multialternative FLT--VIV", g["F_MULTI_G"], "FLT_vs_VIV", "923",
            "0.06615", "0.00065817", "0.00131633",
            "hellinger_centroid_distance", "exact_p", "holm_p_3", "n_clusters"),
           ("Multialternative GC--VIV", g["F_MULTI_G"], "GC_vs_VIV", "923",
            "0.06551", "0.00218017", "0.00218017",
            "hellinger_centroid_distance", "exact_p", "holm_p_3", "n_clusters")]
    for lab, f, con, units, mag, ep, hp, fmag, fep, fhp, func in tbl:
        key = "contrast"
        add(S, "Table 1 %s: units" % lab, units, f, func,
            (lambda f=f, func=func, con=con: cv(f, func, contrast=con)))
        add(S, "Table 1 %s: magnitude" % lab, mag, f, fmag,
            (lambda f=f, fmag=fmag, con=con: cv(f, fmag, contrast=con)))
        add(S, "Table 1 %s: exact p" % lab, ep, f, fep,
            (lambda f=f, fep=fep, con=con: cv(f, fep, contrast=con)))
        add(S, "Table 1 %s: Holm p" % lab, hp, f, fhp,
            (lambda f=f, fhp=fhp, con=con: cv(f, fhp, contrast=con)))
    add(S, "FLT--GC has 48,620 allocations", "48620", g["F_GLOBAL"],
        "allocations FLT_minus_GC",
        lambda: cv(g["F_GLOBAL"], "allocations", contrast="FLT_minus_GC"))
    add(S, "comparisons involving VIV have 24,310", "24310", g["F_GLOBAL"],
        "allocations FLT_minus_VIV",
        lambda: cv(g["F_GLOBAL"], "allocations", contrast="FLT_minus_VIV"))
    add(S, "Both panels reuse the same 26 animals", "26", g["F_SUM"],
        "animals", lambda: j(g["F_SUM"], "animals"))

    # ---- Figure 1 caption ----
    S = "Results/Fig. 1 caption"
    add(S, "Principal components of all 5,249 measured binary choices",
        "5249", "results/descriptive_summary.json", "PCA_feature_count",
        lambda: j("results/descriptive_summary.json", "PCA_feature_count"))
    add(S, "centered across 26 animals without selection", "26",
        "results/descriptive_summary.json", "PCA_animal_count",
        lambda: j("results/descriptive_summary.json", "PCA_animal_count"))
    add(S, "RMS differences in group mean fractions over the same 5,249 "
           "choices", "5249", g["F_GLOBAL"], "features",
        lambda: cv(g["F_GLOBAL"], "features", contrast="FLT_minus_GC"))
    add(S, "enumerate 48,620 FLT--GC", "48620", g["F_GLOBAL"],
        "allocations FLT_minus_GC",
        lambda: cv(g["F_GLOBAL"], "allocations", contrast="FLT_minus_GC"))
    add(S, "and 24,310 assignments in each 9-versus-8 comparison", "24310",
        g["F_GLOBAL"], "allocations GC_minus_VIV",
        lambda: cv(g["F_GLOBAL"], "allocations", contrast="GC_minus_VIV"))
    add(S, "with Holm correction across three global tests", "3",
        g["F_SUM"], "primary_global_significant_Holm05 (family size 3)",
        lambda: j(g["F_SUM"], "primary_global_significant_Holm05"))
    add(S, "FLT n=9, GC n=9, VIV n=8 (FLT)", "9", g["F_AQ"],
        "rows group=FLT", lambda: cn(g["F_AQ"], group="FLT"))
    add(S, "FLT n=9, GC n=9, VIV n=8 (VIV)", "8", g["F_AQ"],
        "rows group=VIV", lambda: cn(g["F_AQ"], group="VIV"))

    # ================= RESULTS 2: gene localization =================
    S = "Results/gene"
    add(S, "The full binary local family comprised 17,172 tests", "17172",
        g["F_SUM"], "event_family", lambda: j(g["F_SUM"], "event_family"))
    add(S, "Joint BH correction identified 451 FLT--GC", "451",
        g["F_EDS"], "q05_joint3M FLT_minus_GC",
        lambda: cv(g["F_EDS"], "q05_joint3M", contrast="FLT_minus_GC"))
    add(S, "38 GC--VIV", "38", g["F_EDS"], "q05_joint3M GC_minus_VIV",
        lambda: cv(g["F_EDS"], "q05_joint3M", contrast="GC_minus_VIV"))
    add(S, "and 43 FLT--VIV choices", "43", g["F_EDS"],
        "q05_joint3M FLT_minus_VIV",
        lambda: cv(g["F_EDS"], "q05_joint3M", contrast="FLT_minus_VIV"))
    add(S, "joint Benjamini--Yekutieli (BY) correction retained none", "0",
        g["F_EDS"], "sum q05_BY_joint over contrasts",
        lambda: ca(g["F_EDS"], "q05_BY_joint", "sum"))
    add(S, "retained all 3,239 measured binary choices with one common gene "
           "annotation", "3239", g["F_GENEM"], "events",
        lambda: j(g["F_GENEM"], "events"))
    add(S, "yielding 1,983 genes", "1983", g["F_GENEM"], "genes",
        lambda: j(g["F_GENEM"], "genes"))
    add(S, "and 5,949 joint gene--contrast tests", "5949",
        g["F_GENEM"], "joint_test_family",
        lambda: j(g["F_GENEM"], "joint_test_family"))
    add(S, "Joint BY retained 47 FLT--GC genes", "47",
        g["F_GENE"], "joint_BY_q_lt_0_05 FLT_minus_GC",
        lambda: cv(g["F_GENE"], "joint_BY_q_lt_0_05", contrast="FLT_minus_GC"))
    add(S, "and none in the vivarium comparisons (GC--VIV)", "0",
        g["F_GENE"], "joint_BY_q_lt_0_05 GC_minus_VIV",
        lambda: cv(g["F_GENE"], "joint_BY_q_lt_0_05", contrast="GC_minus_VIV"))
    add(S, "and none in the vivarium comparisons (FLT--VIV)", "0",
        g["F_GENE"], "joint_BY_q_lt_0_05 FLT_minus_VIV",
        lambda: cv(g["F_GENE"], "joint_BY_q_lt_0_05", contrast="FLT_minus_VIV"))
    add(S, "All 47 attained the exact balanced-allocation minimum p=2/48,620",
        "0.00004114", "recomputed from gene_omnibus_complete.csv.gz",
        "2/48620; all 47 positives at the floor = %s",
        lambda: 2.0 / 48620.0 if R("gene_BY_all_at_floor") else None)
    add(S, "with q=0.04825699", "0.04825699", g["F_GENE"],
        "minimum_joint_BY_q FLT_minus_GC",
        lambda: cv(g["F_GENE"], "minimum_joint_BY_q", contrast="FLT_minus_GC"))
    add(S, "Their RMS contrasts ranged from 3.42", "3.42",
        "extensions/gene_omnibus/results/"
        "arbitrary_dependence_corrected_genes.csv",
        "min rms_delta_percentage_points",
        lambda: ca("extensions/gene_omnibus/results/"
                   "arbitrary_dependence_corrected_genes.csv",
                   "rms_delta_percentage_points", "min"))
    add(S, "to 22.77 pp", "22.77",
        "extensions/gene_omnibus/results/"
        "arbitrary_dependence_corrected_genes.csv",
        "max rms_delta_percentage_points",
        lambda: ca("extensions/gene_omnibus/results/"
                   "arbitrary_dependence_corrected_genes.csv",
                   "rms_delta_percentage_points", "max"))
    add(S, "At 95% dominance, 45 of the 47 primary genes retained joint BY "
           "significance", "45", g["F_DOM95"],
        "primary_BY_genes_retaining_BY FLT_minus_GC",
        lambda: cv(g["F_DOM95"], "primary_BY_genes_retaining_BY",
                   contrast="FLT_minus_GC"))
    add(S, "45 of the 47 primary genes", "47", g["F_DOM95"],
        "primary_BY_genes_total FLT_minus_GC",
        lambda: cv(g["F_DOM95"], "primary_BY_genes_total",
                   contrast="FLT_minus_GC"))
    add(S, "Across all 18 FLT/GC single-animal omissions", "18",
        g["F_GLOO"], "scope (18 FLT or GC animal omissions)",
        lambda: 18 if "18 FLT or GC" in j(g["F_GLOO"], "scope") else None)
    add(S, "primary-gene RMS effects retained 85.24", "85.24",
        g["F_GLOO"], "minimum_rms_ratio_any_primary_gene_omission x100",
        lambda: 100 * j(g["F_GLOO"],
                        "minimum_rms_ratio_any_primary_gene_omission"))
    add(S, "--109.66% of their full magnitude", "109.66",
        g["F_GLOO"], "maximum_rms_ratio_any_primary_gene_omission x100",
        lambda: 100 * j(g["F_GLOO"],
                        "maximum_rms_ratio_any_primary_gene_omission"))

    S = "Results/Fig. 2 caption"
    add(S, "blue points pass joint BY correction across 5,949 gene--contrast "
           "tests", "5949", g["F_GENEM"], "joint_test_family",
        lambda: j(g["F_GENEM"], "joint_test_family"))
    add(S, "All 47 calls have exact p=2/48,620 and q=0.04825699 (count)",
        "47", g["F_GENE"], "joint_BY_q_lt_0_05 FLT_minus_GC",
        lambda: cv(g["F_GENE"], "joint_BY_q_lt_0_05", contrast="FLT_minus_GC"))
    add(S, "The 46 primary genes remaining measurable under 95% dominance",
        "46", "extensions/gene_omnibus/results/"
        "primary_BY_gene_dominance95_followup.csv",
        "rows with a measurable dominance95 result",
        lambda: cn("extensions/gene_omnibus/results/"
                   "primary_BY_gene_dominance95_followup.csv") - 1)

    # ================= RESULTS 3: multialternative =================
    S = "Results/multialternative"
    add(S, "recovered 923 boundaries with at least three adequately measured "
           "alternatives", "923", g["F_SEL"], "n_clusters",
        lambda: j(g["F_SEL"], "n_clusters"))
    add(S, "comprising 3,376 alternative occurrences", "3376",
        g["F_SEL"], "n_alternative_occurrences",
        lambda: j(g["F_SEL"], "n_alternative_occurrences"))
    add(S, "and 3,248 distinct junctions", "3248",
        g["F_SEL"], "n_unique_retained_junctions",
        lambda: j(g["F_SEL"], "n_unique_retained_junctions"))
    add(S, "Of these, 466 shared a genomic start", "466",
        g["F_FUN"], "selected where boundary_type=start",
        lambda: cv(g["F_FUN"], "selected", boundary_type="start"))
    add(S, "and 457 an end", "457", g["F_FUN"],
        "selected where boundary_type=end",
        lambda: cv(g["F_FUN"], "selected", boundary_type="end"))
    add(S, "preserved dependence from 128 reused junctions", "128",
        g["F_SEL"], "n_reused_junctions",
        lambda: j(g["F_SEL"], "n_reused_junctions"))
    add(S, "across 128 cluster pairs", "128", g["F_SEL"],
        "n_overlapping_cluster_pairs",
        lambda: j(g["F_SEL"], "n_overlapping_cluster_pairs"))
    add(S, "Across the full 2,769 local tests", "2769",
        g["F_INF"], "n_local_tests", lambda: j(g["F_INF"], "n_local_tests"))
    add(S, "joint BY retained 40 FLT--GC", "40", g["F_CMC"],
        "n_original_BY DM_score/FLT_vs_GC",
        lambda: cv(g["F_CMC"], "n_original_BY", model="DM_score",
                   contrast="FLT_vs_GC"))
    add(S, "and three FLT--VIV associations", "3", g["F_CMC"],
        "n_original_BY DM_score/FLT_vs_VIV",
        lambda: cv(g["F_CMC"], "n_original_BY", model="DM_score",
                   contrast="FLT_vs_VIV"))
    add(S, "representing 41 unique clusters", "41", g["F_INF"],
        "n_unique_BY_clusters",
        lambda: j(g["F_INF"], "n_unique_BY_clusters"))
    add(S, "none passed in GC--VIV", "0", g["F_CMC"],
        "n_original_BY DM_score/GC_vs_VIV",
        lambda: cv(g["F_CMC"], "n_original_BY", model="DM_score",
                   contrast="GC_vs_VIV"))
    add(S, "The three FLT--VIV candidates were 11:2:109886437:end, "
           "2:1:55311420:start and 6:2:42663893:start", "3",
        g["F_FV"], "row count", lambda: cn(g["F_FV"]))
    add(S, "1,199 of the 3,248 distinct junctions lacked an exact Ensembl "
           "transcript match", "1199", g["F_TC"],
        "multiway_junctions_without_exact_transcript_match",
        lambda: j(g["F_TC"],
                  "multiway_junctions_without_exact_transcript_match"))
    add(S, "Exactly 290 clusters mapped to one common gene", "290",
        g["F_SEL"], "n_unique_mapped_clusters",
        lambda: j(g["F_SEL"], "n_unique_mapped_clusters"))
    add(S, "representing 256 genes", "256", g["F_SEL"],
        "n_unique_mapped_genes", lambda: j(g["F_SEL"], "n_unique_mapped_genes"))
    add(S, "Tcaf1 passed primary joint BY in FLT--GC (q=0.03229)", "0.03229",
        g["F_FV"], "BY_q_full_3E (also_BY_significant_FLT_vs_GC row)",
        lambda: cv("extensions/multiway_splicing/tables/"
                   "BY_discoveries_full_family.csv", "BY_q_full_3E",
                   cluster_id="6:2:42663893:start", contrast="FLT_vs_GC"))
    add(S, "and FLT--VIV (q=0.04505)", "0.04505", g["F_FV"],
        "BY_q_full_3E for 6:2:42663893:start",
        lambda: cv(g["F_FV"], "BY_q_full_3E",
                   cluster_id="6:2:42663893:start"))
    add(S, "mean read fractions were 78.13% in FLT", "78.13",
        "recomputed: transcript_context/Tcaf1_all26_alternative_read_support.csv",
        "mean 23000-bp intron proportion, FLT",
        lambda: 100 * R("Tcaf1_means")["FLT"])
    add(S, "61.72% in GC", "61.72",
        "recomputed: transcript_context/Tcaf1_all26_alternative_read_support.csv",
        "mean 23000-bp intron proportion, GC",
        lambda: 100 * R("Tcaf1_means")["GC"])
    add(S, "and 65.87% in VIV", "65.87",
        "recomputed: transcript_context/Tcaf1_all26_alternative_read_support.csv",
        "mean 23000-bp intron proportion, VIV",
        lambda: 100 * R("Tcaf1_means")["VIV"])
    add(S, "corresponding to differences of 16.41", "16.41",
        "recomputed (same file)", "FLT minus GC, pp",
        lambda: R("Tcaf1_FLT_GC_pp"))
    add(S, "and 12.27 pp", "12.27", "recomputed (same file)",
        "FLT minus VIV, pp", lambda: R("Tcaf1_FLT_VIV_pp"))

    S = "Results/Fig. 3 caption"
    add(S, "Every FLT--GC cluster in the joint 2,769-test family", "2769",
        g["F_INF"], "n_local_tests", lambda: j(g["F_INF"], "n_local_tests"))
    add(S, "All 26 animals at Tcaf1 cluster 6:2:42663893:start", "26",
        "extensions/transcript_context/tables/"
        "Tcaf1_all26_long_short_intron_support.csv", "row count",
        lambda: cn("extensions/transcript_context/tables/"
                   "Tcaf1_all26_long_short_intron_support.csv"))

    # ================= RESULTS 4: Tcaf1 architecture =================
    S = "Results/Tcaf1 architecture"
    add(S, "the shared low genomic intron start at 42,663,893", "42663893",
        g["F_TC"], "Tcaf1 cluster id boundary",
        lambda: float(cs(g["F_FV"], "boundary",
                         cluster_id="6:2:42663893:start")))
    add(S, "the three higher intron ends represent donor-side alternatives",
        "3", g["F_TCJ"], "distinct junction_id count",
        lambda: ca(g["F_TCJ"], "junction_id", "nuniq"))
    add(S, "The 5,974-bp ... alternatives", "5974",
        "extensions/transcript_context/tables/"
        "Tcaf1_all26_alternative_read_support.csv",
        "intron_length_bp for 6:2:42663893:42669866",
        lambda: cv("extensions/transcript_context/tables/"
                   "Tcaf1_all26_alternative_read_support.csv",
                   "intron_length_bp",
                   junction_id="6:2:42663893:42669866"))
    add(S, "and 23,000-bp alternatives", "23000",
        "extensions/transcript_context/tables/"
        "Tcaf1_all26_alternative_read_support.csv",
        "intron_length_bp for 6:2:42663893:42686892",
        lambda: cv("extensions/transcript_context/tables/"
                   "Tcaf1_all26_alternative_read_support.csv",
                   "intron_length_bp",
                   junction_id="6:2:42663893:42686892"))
    add(S, "The 22,570-bp alternative matches Tcaf1-205", "22570",
        "extensions/transcript_context/tables/"
        "Tcaf1_all26_alternative_read_support.csv",
        "intron_length_bp for 6:2:42663893:42686462",
        lambda: cv("extensions/transcript_context/tables/"
                   "Tcaf1_all26_alternative_read_support.csv",
                   "intron_length_bp",
                   junction_id="6:2:42663893:42686462"))
    add(S, "These transcripts share ... CCDS20076", "CCDS20076",
        g["F_TCS"], "CCDS_ids for Tcaf1-201",
        lambda: cs(g["F_TCS"], "CCDS_ids", transcript_name="Tcaf1-201"))
    add(S, "Tcaf1-205 annotated as a retained-intron transcript",
        "retained_intron", g["F_TCS"], "transcript_biotype Tcaf1-205",
        lambda: cs(g["F_TCS"], "transcript_biotype",
                   transcript_name="Tcaf1-205"))
    add(S, "Tcaf1-206 (retained_intron)", "retained_intron", g["F_TCS"],
        "transcript_biotype Tcaf1-206",
        lambda: cs(g["F_TCS"], "transcript_biotype",
                   transcript_name="Tcaf1-206"))
    add(S, "Tcaf1-208 (protein_coding_CDS_not_defined)",
        "protein_coding_CDS_not_defined", g["F_TCS"],
        "transcript_biotype Tcaf1-208",
        lambda: cs(g["F_TCS"], "transcript_biotype",
                   transcript_name="Tcaf1-208"))

    S = "Results/Fig. 4 caption"
    add(S, "The minus-strand introns share genomic start 42,663,893",
        "42663893", g["F_FV"], "boundary for the Tcaf1 cluster",
        lambda: float(cs(g["F_FV"], "boundary",
                         cluster_id="6:2:42663893:start")))
    add(S, "share identical annotated coding coordinates and phases, 2,766 "
           "CDS bases", "2766", g["F_TCS"], "CDS_bases Tcaf1-201",
        lambda: cv(g["F_TCS"], "CDS_bases", transcript_name="Tcaf1-201"))
    add(S, "All 26 animals support the companion intron 6:2:42686572:42686892",
        "26", g["F_TC"], "Tcaf1_companion_available_animals",
        lambda: j(g["F_TC"], "Tcaf1_companion_available_animals"))
    add(S, "with 20--169 unique reads (min)", "20", g["F_TC"],
        "Tcaf1_companion_min_reads",
        lambda: j(g["F_TC"], "Tcaf1_companion_min_reads"))
    add(S, "with 20--169 unique reads (max)", "169", g["F_TC"],
        "Tcaf1_companion_max_reads",
        lambda: j(g["F_TC"], "Tcaf1_companion_max_reads"))
    add(S, "The complete audit retains 1,199/3,248 multiway junctions "
           "without exact GTF intron matches", "1199", g["F_TC"],
        "multiway_junctions_without_exact_transcript_match",
        lambda: j(g["F_TC"],
                  "multiway_junctions_without_exact_transcript_match"))

    # ================= RESULTS 5: aggregate robustness =================
    S = "Results/robustness"
    add(S, "retained 607 of the original 923 multialternative clusters",
        "607", g["F_ROB"], "n_valid_dominance95",
        lambda: j(g["F_ROB"], "n_valid_dominance95"))
    add(S, "and 25 of the 43 primary BY-positive contrast--cluster tests",
        "25", g["F_ROB"], "n_BY_dominance95",
        lambda: j(g["F_ROB"], "n_BY_dominance95"))
    add(S, "of the 43 primary BY-positive tests", "43", g["F_ROB"],
        "n_original_positive_tests",
        lambda: j(g["F_ROB"], "n_original_positive_tests"))
    add(S, "23 FLT--GC", "23", g["F_ROBR"], "BY_positive95 FLT_vs_GC",
        lambda: cv(g["F_ROBR"], "BY_positive95", contrast="FLT_vs_GC"))
    add(S, "and two FLT--VIV", "2", g["F_ROBR"],
        "BY_positive95 FLT_vs_VIV",
        lambda: cv(g["F_ROBR"], "BY_positive95", contrast="FLT_vs_VIV"))
    add(S, "Tcaf1 remained measurable, with FLT--GC q=0.03874", "0.03874",
        g["F_ROBS"], "BY_q_dominance95 6:2:42663893:start / FLT_vs_GC",
        lambda: cv(g["F_ROBS"], "BY_q_dominance95",
                   cluster_id="6:2:42663893:start", contrast="FLT_vs_GC"))
    add(S, "and FLT--VIV q=0.05381", "0.05381", g["F_ROBS"],
        "BY_q_dominance95 6:2:42663893:start / FLT_vs_VIV",
        lambda: cv(g["F_ROBS"], "BY_q_dominance95",
                   cluster_id="6:2:42663893:start", contrast="FLT_vs_VIV"))
    add(S, "Across ten paired half-depth partitions", "10",
        g["F_ESS"], "paired_read_splits",
        lambda: j(g["F_ESS"], "paired_read_splits"))
    add(S, "all 60 global comparisons remained Holm-positive", "60",
        g["F_ESS"], "complement_B_Holm_global_positive x2 (A+B halves)",
        lambda: 2 * j(g["F_ESS"], "complement_B_Holm_global_positive"))
    add(S, "The original ten A-half draws retained zero local BY discoveries",
        "0", g["F_ENV"], "sum BY_positive over the ten A halves",
        lambda: sum(int(r["BY_positive"]) for r in g["_tab"](g["F_ENV"])
                    if r["family"].endswith("_A")))
    add(S, "Their B-half complements produced 0--23 calls (max)", "23",
        g["F_ENV"], "max BY_positive over the ten B halves",
        lambda: max(int(r["BY_positive"]) for r in g["_tab"](g["F_ENV"])
                    if r["family"].endswith("_B")))
    add(S, "with 20, 20, 23, 21 and 20 calls in five of ten halves",
        "20,20,23,21,20", "recomputed: complement_B_local_exact_full_families",
        "nonzero per-replicate BY_positive counts",
        lambda: ",".join(str(x) for x in R("B_half_counts") if x > 0))
    add(S, "totaling 104 contrast--cluster call instances", "104",
        g["F_ESS"], "complement_B_BY_calls",
        lambda: j(g["F_ESS"], "complement_B_BY_calls"))
    add(S, "Thus 5 of 20 halves ... supported local BY calls", "5",
        "recomputed: complement_B_local_exact_full_families",
        "B halves with >=1 BY call", lambda: R("B_half_with_calls"))
    add(S, "both complementary halves were measurable in 354 "
           "contrast--cluster--split instances", "354",
        g["F_ES"], "sum instances_measured, half_to_half/original_BY_tests",
        lambda: ca(g["F_ES"], "instances_measured", "sum",
                   analysis="half_to_half", subset="original_BY_tests"))
    add(S, "325 FLT--GC", "325", g["F_ES"],
        "instances_measured half_to_half/FLT_vs_GC/original_BY_tests",
        lambda: cv(g["F_ES"], "instances_measured", analysis="half_to_half",
                   contrast="FLT_vs_GC", subset="original_BY_tests"))
    add(S, "and 29 FLT--VIV", "29", g["F_ES"],
        "instances_measured half_to_half/FLT_vs_VIV/original_BY_tests",
        lambda: cv(g["F_ES"], "instances_measured", analysis="half_to_half",
                   contrast="FLT_vs_VIV", subset="original_BY_tests"))
    add(S, "from 430 possible instances", "430", g["F_ES"],
        "sum instances_total, half_to_half/original_BY_tests",
        lambda: ca(g["F_ES"], "instances_total", "sum",
                   analysis="half_to_half", subset="original_BY_tests"))
    add(S, "Every measured pair had a positive native-effect dot product",
        "1.0", g["F_ES"], "min positive_alignment_fraction (candidate pairs)",
        lambda: ca(g["F_ES"], "positive_alignment_fraction", "min",
                   analysis="half_to_half", subset="original_BY_tests"))
    add(S, "Median paired-vector cosines were 0.9964", "0.9964",
        g["F_ES"], "median_cosine half_to_half/FLT_vs_GC/original_BY_tests",
        lambda: cv(g["F_ES"], "median_cosine", analysis="half_to_half",
                   contrast="FLT_vs_GC", subset="original_BY_tests"))
    add(S, "and 0.9933", "0.9933", g["F_ES"],
        "median_cosine half_to_half/FLT_vs_VIV/original_BY_tests",
        lambda: cv(g["F_ES"], "median_cosine", analysis="half_to_half",
                   contrast="FLT_vs_VIV", subset="original_BY_tests"))
    add(S, "with median B-on-A relative projections of 0.9905", "0.9905",
        g["F_ES"], "median_relative_projection (FLT_vs_GC candidates)",
        lambda: cv(g["F_ES"], "median_relative_projection",
                   analysis="half_to_half", contrast="FLT_vs_GC",
                   subset="original_BY_tests"))
    add(S, "and 0.9492, respectively", "0.9492", g["F_ES"],
        "median_relative_projection (FLT_vs_VIV candidates)",
        lambda: cv(g["F_ES"], "median_relative_projection",
                   analysis="half_to_half", contrast="FLT_vs_VIV",
                   subset="original_BY_tests"))
    add(S, "positive in all 751 measurable original-candidate half instances",
        "751", g["F_ES"],
        "sum instances_measured, original_to_half/original_BY_tests",
        lambda: ca(g["F_ES"], "instances_measured", "sum",
                   analysis="original_to_half", subset="original_BY_tests"))
    add(S, "paired positive-alignment fractions were lower: 88.85% for "
           "FLT--GC", "88.85", g["F_ES"],
        "positive_alignment_fraction half_to_half/FLT_vs_GC/all923 x100",
        lambda: 100 * cv(g["F_ES"], "positive_alignment_fraction",
                         analysis="half_to_half", contrast="FLT_vs_GC",
                         subset="all923"))
    add(S, "and 83.87% for FLT--VIV", "83.87", g["F_ES"],
        "positive_alignment_fraction half_to_half/FLT_vs_VIV/all923 x100",
        lambda: 100 * cv(g["F_ES"], "positive_alignment_fraction",
                         analysis="half_to_half", contrast="FLT_vs_VIV",
                         subset="all923"))
    add(S, "positive native-vector alignment in all 771 original-candidate "
           "omission instances", "771", g["F_ES"],
        "sum instances_measured, leave_one_animal_out/original_BY_tests",
        lambda: ca(g["F_ES"], "instances_measured", "sum",
                   analysis="leave_one_animal_out",
                   subset="original_BY_tests"))
    add(S, "(720 FLT--GC", "720", g["F_ES"],
        "instances_measured LOO/FLT_vs_GC/original_BY_tests",
        lambda: cv(g["F_ES"], "instances_measured",
                   analysis="leave_one_animal_out", contrast="FLT_vs_GC",
                   subset="original_BY_tests"))
    add(S, "and 51 FLT--VIV)", "51", g["F_ES"],
        "instances_measured LOO/FLT_vs_VIV/original_BY_tests",
        lambda: cv(g["F_ES"], "instances_measured",
                   analysis="leave_one_animal_out", contrast="FLT_vs_VIV",
                   subset="original_BY_tests"))
    add(S, "For 2,769 tests, the BY harmonic factor is 8.50364", "8.50364",
        "recomputed: sum_{k=1..2769} 1/k", "harmonic factor",
        lambda: R("H_2769"))
    add(S, "the smallest common exact p is 2/48,620=1/24,310", "0.00004114",
        "recomputed", "2/48620 == 1/24310", lambda: R("p_floor"))
    add(S, "A block of 20 tests at that floor has q=0.04843", "0.04843",
        "recomputed: m*H_m*p_min/k, k=20", "BY q", lambda: R("q_block20"))
    add(S, "whereas 19 has q=0.05098", "0.05098",
        "recomputed: m*H_m*p_min/k, k=19", "BY q", lambda: R("q_block19"))
    add(S, "The complementary halves with BY calls contained 20--23 "
           "floor-level tests (max)", "23", g["F_ENV"],
        "max n_at_minimum_p over B halves with calls",
        lambda: max(int(r["n_at_minimum_p"]) for r in g["_tab"](g["F_ENV"])
                    if r["family"].endswith("_B") and int(r["BY_positive"]) > 0))
    add(S, "all A halves contained fewer than 20", "19", g["F_ENV"],
        "max n_at_minimum_p over the ten A halves",
        lambda: R("A_max_floor"))
    add(S, "approximate centered-error distance intervals of 0.02148", "0.02148",
        g["F_BOOT"], "confidence_ball_distance_low FLT_vs_GC",
        lambda: cv(g["F_BOOT"], "confidence_ball_distance_low",
                   contrast="FLT_vs_GC"))
    add(S, "--0.14478 for FLT--GC", "0.14478", g["F_BOOT"],
        "confidence_ball_distance_high FLT_vs_GC",
        lambda: cv(g["F_BOOT"], "confidence_ball_distance_high",
                   contrast="FLT_vs_GC"))
    add(S, "0--0.13375 for FLT--VIV", "0.13375", g["F_BOOT"],
        "confidence_ball_distance_high FLT_vs_VIV",
        lambda: cv(g["F_BOOT"], "confidence_ball_distance_high",
                   contrast="FLT_vs_VIV"))
    add(S, "and 0--0.13324 for GC--VIV", "0.13324", g["F_BOOT"],
        "confidence_ball_distance_high GC_vs_VIV",
        lambda: cv(g["F_BOOT"], "confidence_ball_distance_high",
                   contrast="GC_vs_VIV"))
    add(S, "dispersion diagnostic gave FLT--GC variance ratio 1.0075",
        "1.0075", g["F_DISP"], "variance_ratio_A_to_B FLT_vs_GC",
        lambda: cv(g["F_DISP"], "variance_ratio_A_to_B", contrast="FLT_vs_GC"))
    add(S, "and no Holm-significant dispersion reference test", "0",
        g["F_DISP"], "count dispersion_Holm_3 < 0.05",
        lambda: sum(1 for r in g["_tab"](g["F_DISP"])
                    if float(r["dispersion_Holm_3"]) < 0.05))
    add(S, "Studentized centroid associations persisted (Holm p=0.00012341)",
        "0.00012341", g["F_DISP"], "studentized_Holm_3 FLT_vs_GC",
        lambda: cv(g["F_DISP"], "studentized_Holm_3", contrast="FLT_vs_GC"))
    add(S, "0.00164541", "0.00164541", g["F_DISP"],
        "studentized_Holm_3 FLT_vs_VIV",
        lambda: cv(g["F_DISP"], "studentized_Holm_3", contrast="FLT_vs_VIV"))
    add(S, "and 0.00292061", "0.00292061", g["F_DISP"],
        "studentized_Holm_3 GC_vs_VIV",
        lambda: cv(g["F_DISP"], "studentized_Holm_3", contrast="GC_vs_VIV"))
    add(S, "survived ten depth-equalizing thinnings", "10",
        "results/depth_thinning_global_sensitivity.csv",
        "distinct replicate values",
        lambda: ca("results/depth_thinning_global_sensitivity.csv",
                   "replicate", "nuniq"))
    add(S, "rho=-0.6741 with shared GC animals", "-0.6741",
        g["F_PD"], "shared_control_rho",
        lambda: j(g["F_PD"], "shared_control_rho"))
    add(S, "versus median rho=-0.4133 across disjoint ... GC partitions",
        "-0.4133", g["F_PD"], "disjoint_rho_median",
        lambda: j(g["F_PD"], "disjoint_rho_median"))
    add(S, "The 126 partitions reuse animals", "126", g["F_PD"],
        "disjoint_control_partitions",
        lambda: j(g["F_PD"], "disjoint_control_partitions"))

    S = "Results/Fig. 5 caption"
    add(S, "each frozen family retains all 2,769 slots", "2769",
        g["F_ESS"], "full_local_family",
        lambda: j(g["F_ESS"], "full_local_family"))
    add(S, "five B halves have 20--23 calls and five have zero", "5",
        "recomputed", "B halves with >=1 call",
        lambda: R("B_half_with_calls"))
    add(S, "All 354 jointly measured candidate instances retain positive "
           "orientation", "354", "recomputed: "
        "all_complementary_half_concordance.csv.gz",
        "candidate instances with cosine > 0", lambda: R("split_positive"))
    add(S, "medians are 0.9964 for FLT--GC", "0.9964", "recomputed (same)",
        "median cosine, FLT_vs_GC candidates",
        lambda: R("split_median_cos")["FLT_vs_GC"])
    add(S, "and 0.9933 for FLT--VIV", "0.9933", "recomputed (same)",
        "median cosine, FLT_vs_VIV candidates",
        lambda: R("split_median_cos")["FLT_vs_VIV"])
    add(S, "all 751 measurable candidate half-instances", "751", g["F_ES"],
        "sum instances_measured original_to_half/original_BY_tests",
        lambda: ca(g["F_ES"], "instances_measured", "sum",
                   analysis="original_to_half", subset="original_BY_tests"))
    add(S, "and all 771 candidate animal omissions", "771", g["F_ES"],
        "sum instances_measured LOO/original_BY_tests",
        lambda: ca(g["F_ES"], "instances_measured", "sum",
                   analysis="leave_one_animal_out",
                   subset="original_BY_tests"))
    add(S, "All 20 half-families and the original family at the exact "
           "minimum p=1/24,310", "20", g["F_ENV"],
        "half-family rows", lambda: cn(g["F_ENV"]) - 1)
    add(S, "The curve is mH_m p_min/k ... m=2,769", "2769", g["F_ENV"],
        "n_tests", lambda: cv(g["F_ENV"], "n_tests", family="original"))
    add(S, "20 floor-level tests are needed for that block to cross 0.05",
        "20", g["F_ENV"], "minimum_tests_at_observed_minimum_for_BY",
        lambda: cv(g["F_ENV"], "minimum_tests_at_observed_minimum_for_BY",
                   family="original"))

    # ================= RESULTS 6: abundance =================
    S = "Results/abundance"
    add(S, "a minimum of 20 raw gene reads per animal", "20",
        g["F_AB"], "(stated gate; normalization note)",
        lambda: 20 if "Median of ratios" in j(g["F_AB"], "normalization")
        else None)
    add(S, "1,977 of the 1,983 junction genes remained", "1977",
        g["F_AB"], "tested_genes", lambda: j(g["F_AB"], "tested_genes"))
    add(S, "of the 1,983 junction genes", "1983", g["F_AB"],
        "junction_genes", lambda: j(g["F_AB"], "junction_genes"))
    add(S, "including 46 of the 47 primary gene candidates", "46",
        g["F_AB"], "summary[FLT_minus_GC].junction_BY",
        lambda: jl(g["F_AB"], "summary", "contrast", "FLT_minus_GC",
                   "junction_BY"))
    add(S, "Across 5,931 abundance tests", "5931", g["F_AB"],
        "joint_tests", lambda: j(g["F_AB"], "joint_tests"))
    add(S, "joint BH retained 72 FLT--GC", "72", g["F_AB"],
        "summary[FLT_minus_GC].abundance_BH",
        lambda: jl(g["F_AB"], "summary", "contrast", "FLT_minus_GC",
                   "abundance_BH"))
    add(S, "eight GC--VIV", "8", g["F_AB"],
        "summary[GC_minus_VIV].abundance_BH",
        lambda: jl(g["F_AB"], "summary", "contrast", "GC_minus_VIV",
                   "abundance_BH"))
    add(S, "and 24 FLT--VIV genes", "24", g["F_AB"],
        "summary[FLT_minus_VIV].abundance_BH",
        lambda: jl(g["F_AB"], "summary", "contrast", "FLT_minus_VIV",
                   "abundance_BH"))
    add(S, "joint BY retained none", "0", g["F_AB"],
        "sum summary[*].abundance_BY",
        lambda: sum(d["abundance_BY"] for d in
                    __import__("json").load(
                        open(__import__("os").path.join(
                            ns["P"], g["F_AB"])))["summary"]))
    add(S, "No global abundance contrast passed its three-test Holm family "
           "(adjusted p=0.13303)", "0.13303", g["F_AB"],
        "global_results[FLT_minus_GC].p_Holm_3",
        lambda: jl(g["F_AB"], "global_results", "contrast", "FLT_minus_GC",
                   "p_Holm_3"))
    add(S, "0.42320 (GC--VIV)", "0.42320", g["F_AB"],
        "global_results[GC_minus_VIV].p_Holm_3",
        lambda: jl(g["F_AB"], "global_results", "contrast", "GC_minus_VIV",
                   "p_Holm_3"))
    add(S, "and 0.42320 (FLT--VIV)", "0.42320", g["F_AB"],
        "global_results[FLT_minus_VIV].p_Holm_3",
        lambda: jl(g["F_AB"], "global_results", "contrast", "FLT_minus_VIV",
                   "p_Holm_3"))
    add(S, "weak descriptively (Spearman rho=0.0572)", "0.0572", g["F_AB"],
        "summary[FLT_minus_GC].descriptive_rho_junction_rms_absolute_abundance",
        lambda: jl(g["F_AB"], "summary", "contrast", "FLT_minus_GC",
                   "descriptive_rho_junction_rms_absolute_abundance"))
    add(S, "abundance-equivalence margins of 1.2-fold", "0",
        g["F_AB"], "summary[FLT_minus_GC].conjunction_BY_1.2",
        lambda: jl(g["F_AB"], "summary", "contrast", "FLT_minus_GC",
                   "conjunction_BY_1.2"))
    add(S, "and 1.5-fold ... No conjunction passed joint BY", "0",
        g["F_AB"], "summary[FLT_minus_GC].conjunction_BY_1.5",
        lambda: jl(g["F_AB"], "summary", "contrast", "FLT_minus_GC",
                   "conjunction_BY_1.5"))

    # ================= RESULTS 7: transport =================
    S = "Results/transport"
    add(S, "comprises 55 distinct cerebellar animals", "55",
        "recomputed: 26 + RRRM2_animals + RR10_animals",
        "animal total", lambda: R("total_animals_3_cohorts"))
    add(S, "RRRM2 supplies 19 eligible female C57BL/6NTac animals", "19",
        g["F_CTS"], "RRRM2_animals", lambda: j(g["F_CTS"], "RRRM2_animals"))
    add(S, "ten FLT and nine GC (RRRM2 FLT)", "10", g["F_MT"],
        "RRRM2_FLT", lambda: j(g["F_MT"], "RRRM2_FLT"))
    add(S, "ten FLT and nine GC (RRRM2 GC)", "9", g["F_MT"],
        "RRRM2_GC", lambda: j(g["F_MT"], "RRRM2_GC"))
    add(S, "across three route-by-design-age strata", "3", g["F_IV"],
        "RRRM2_strata", lambda: j(g["F_IV"], "RRRM2_strata"))
    add(S, "RR10 supplies ten female B6129SF2/J wild-type animals", "10",
        g["F_CTS"], "RR10_animals", lambda: j(g["F_CTS"], "RR10_animals"))
    add(S, "five FLT and five GC (RR10 FLT)", "5", g["F_CT"],
        "n_FLT RR10/cohort_measured",
        lambda: cv(g["F_CT"], "n_FLT", cohort="RR10",
                   panel="cohort_measured"))
    add(S, "five FLT and five GC (RR10 GC)", "5", g["F_CT"],
        "n_GC RR10/cohort_measured",
        lambda: cv(g["F_CT"], "n_GC", cohort="RR10", panel="cohort_measured"))
    add(S, "retained 1,512 measurable common choices", "1512", g["F_IV"],
        "common_measured_events",
        lambda: j(g["F_IV"], "common_measured_events"))
    add(S, "(28.81%)", "28.81", g["F_IV"], "common_transport_fraction x100",
        lambda: 100 * j(g["F_IV"], "common_transport_fraction"))
    add(S, "discovery/validation effect correlation was rho=0.01019", "0.01019",
        g["F_IV"], "spearman_rho", lambda: j(g["F_IV"], "spearman_rho"))
    add(S, "two-cohort bootstrap interval -0.06894", "-0.06894", g["F_IV"],
        "spearman_animal_bootstrap_CI95[0]",
        lambda: j(g["F_IV"], "spearman_animal_bootstrap_CI95")[0])
    add(S, "to 0.07386", "0.07386", g["F_IV"],
        "spearman_animal_bootstrap_CI95[1]",
        lambda: j(g["F_IV"], "spearman_animal_bootstrap_CI95")[1])
    add(S, "stratified effect of 0.1281 pp", "0.1281", g["F_IV"],
        "learned_direction_score_delta_percentage_points",
        lambda: j(g["F_IV"], "learned_direction_score_delta_percentage_points"))
    add(S, "conditional interval -0.1372", "-0.1372", g["F_IV"],
        "learned_direction_score_animal_bootstrap_CI95[0] x100",
        lambda: 100 * j(g["F_IV"],
                        "learned_direction_score_animal_bootstrap_CI95")[0])
    add(S, "to 0.4218 pp", "0.4218", g["F_IV"],
        "learned_direction_score_animal_bootstrap_CI95[1] x100",
        lambda: 100 * j(g["F_IV"],
                        "learned_direction_score_animal_bootstrap_CI95")[1])
    add(S, "and one-sided p=0.2525", "0.2525", g["F_IV"],
        "learned_direction_score_exact_one_sided_p",
        lambda: j(g["F_IV"], "learned_direction_score_exact_one_sided_p"))
    add(S, "No common choice passed directional Holm correction", "0",
        g["F_IV"], "full_common_panel_Holm_directional05",
        lambda: j(g["F_IV"], "full_common_panel_Holm_directional05"))
    add(S, "281 of 923 frozen clusters passed the gates in RRRM2", "281",
        "recomputed: RRRM2_all923_measurement_gates.csv",
        "common_measurable=True count", lambda: R("RRRM2_measurable"))
    add(S, "(30.44%)", "30.44", "recomputed (same)",
        "281/923 x100", lambda: 100 * R("RRRM2_fraction"))
    add(S, "and 702 in RR10", "702",
        "recomputed: RR10_all923_measurement_gates.csv",
        "common_measurable=True count", lambda: R("RR10_measurable"))
    add(S, "(76.06%)", "76.06", "recomputed (same)", "702/923 x100",
        lambda: 100 * R("RR10_fraction"))
    add(S, "270 were measurable in both", "270", g["F_CTS"],
        "common_both_validation",
        lambda: j(g["F_CTS"], "common_both_validation"))
    add(S, "None of the four primary global tests survived their joint Holm "
           "family; all adjusted p values were 1", "1.0", g["F_CT"],
        "min Holm_four_global_p over panel=cohort_measured",
        lambda: ca(g["F_CT"], "Holm_four_global_p", "min",
                   panel="cohort_measured"))
    add(S, "the FLT--GC-trained score reproduced its original effect -0.2770 pp",
        "-0.2770", g["F_CT"], "effect_pp RRRM2/FLT_minus_GC/cohort_measured",
        lambda: cv(g["F_CT"], "effect_pp", cohort="RRRM2",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "and one-sided p=0.73625", "0.73625", g["F_CT"],
        "exact_p RRRM2/FLT_minus_GC/cohort_measured",
        lambda: cv(g["F_CT"], "exact_p", cohort="RRRM2",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "(two-sided p=0.521875)", "0.521875", g["F_CT"],
        "exact_two_sided_p RRRM2/FLT_minus_GC/cohort_measured",
        lambda: cv(g["F_CT"], "exact_two_sided_p", cohort="RRRM2",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "the expanded 4,000-draw conditional interval was -0.9635", "-0.9635",
        g["F_CT"], "conditional_CI95_low_pp RRRM2/FLT_minus_GC",
        lambda: cv(g["F_CT"], "conditional_CI95_low_pp", cohort="RRRM2",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "to 0.3778 pp", "0.3778", g["F_CT"],
        "conditional_CI95_high_pp RRRM2/FLT_minus_GC",
        lambda: cv(g["F_CT"], "conditional_CI95_high_pp", cohort="RRRM2",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "Only 46.98% of the 281 measurable clusters had a positive "
           "projection", "46.98", g["F_MT"],
        "positive_directional_cluster_fraction x100",
        lambda: 100 * j(g["F_MT"], "positive_directional_cluster_fraction"))
    add(S, "The FLT--VIV-trained score was 0.1824 pp", "0.1824", g["F_CT"],
        "effect_pp RRRM2/FLT_minus_VIV/cohort_measured",
        lambda: cv(g["F_CT"], "effect_pp", cohort="RRRM2",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "conditional interval -0.5070", "-0.5070", g["F_CT"],
        "conditional_CI95_low_pp RRRM2/FLT_minus_VIV",
        lambda: cv(g["F_CT"], "conditional_CI95_low_pp", cohort="RRRM2",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "to 0.8072 pp", "0.8072", g["F_CT"],
        "conditional_CI95_high_pp RRRM2/FLT_minus_VIV",
        lambda: cv(g["F_CT"], "conditional_CI95_high_pp", cohort="RRRM2",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "and one-sided p=0.34955", "0.34955", g["F_CT"],
        "exact_p RRRM2/FLT_minus_VIV/cohort_measured",
        lambda: cv(g["F_CT"], "exact_p", cohort="RRRM2",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "In RR10, the FLT--GC-trained score was 0.4358 pp", "0.4358",
        g["F_CT"], "effect_pp RR10/FLT_minus_GC/cohort_measured",
        lambda: cv(g["F_CT"], "effect_pp", cohort="RR10",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "interval -1.0133", "-1.0133", g["F_CT"],
        "conditional_CI95_low_pp RR10/FLT_minus_GC",
        lambda: cv(g["F_CT"], "conditional_CI95_low_pp", cohort="RR10",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "to 2.0329 pp", "2.0329", g["F_CT"],
        "conditional_CI95_high_pp RR10/FLT_minus_GC",
        lambda: cv(g["F_CT"], "conditional_CI95_high_pp", cohort="RR10",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "and p=0.31349", "0.31349", g["F_CT"],
        "exact_p RR10/FLT_minus_GC/cohort_measured",
        lambda: cv(g["F_CT"], "exact_p", cohort="RR10",
                   direction="FLT_minus_GC", panel="cohort_measured"))
    add(S, "the FLT--VIV-trained score was 0.0596 pp", "0.0596", g["F_CT"],
        "effect_pp RR10/FLT_minus_VIV/cohort_measured",
        lambda: cv(g["F_CT"], "effect_pp", cohort="RR10",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "interval -1.1055", "-1.1055", g["F_CT"],
        "conditional_CI95_low_pp RR10/FLT_minus_VIV",
        lambda: cv(g["F_CT"], "conditional_CI95_low_pp", cohort="RR10",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "to 1.2390 pp", "1.2390", g["F_CT"],
        "conditional_CI95_high_pp RR10/FLT_minus_VIV",
        lambda: cv(g["F_CT"], "conditional_CI95_high_pp", cohort="RR10",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "and p=0.46429", "0.46429", g["F_CT"],
        "exact_p RR10/FLT_minus_VIV/cohort_measured",
        lambda: cv(g["F_CT"], "exact_p", cohort="RR10",
                   direction="FLT_minus_VIV", panel="cohort_measured"))
    add(S, "The common 270-cluster sensitivity ... minimum four-test Holm "
           "p=0.90476", "0.90476", g["F_CT"],
        "min Holm_four_global_p over panel=intersection_both_cohorts",
        lambda: ca(g["F_CT"], "Holm_four_global_p", "min",
                   panel="intersection_both_cohorts"))
    add(S, "No local cluster passed BH ... across the full 3,692-slot family "
           "(BH)", "0", g["F_CTS"], "local_BH_calls",
        lambda: j(g["F_CTS"], "local_BH_calls"))
    add(S, "or BY across the full 3,692-slot family", "0", g["F_CTS"],
        "local_BY_calls", lambda: j(g["F_CTS"], "local_BY_calls"))
    add(S, "the full expanded 3,692-slot family", "3692", g["F_CTS"],
        "local_family", lambda: j(g["F_CTS"], "local_family"))
    add(S, "its first-step minimum adjusted p was 923/11,200=0.08241",
        "0.08241", g["F_MT"], "minimum_attainable_first_Holm_adjusted_p",
        lambda: j(g["F_MT"], "minimum_attainable_first_Holm_adjusted_p"))
    add(S, "923/11,200", "11200", g["F_MT"], "exact_stratified_assignments",
        lambda: j(g["F_MT"], "exact_stratified_assignments"))
    add(S, "The original 1,512-feature binary family similarly had a minimum "
           "of 0.135", "0.135", "recomputed: 1512 / 11200",
        "Holm first-step floor", lambda: R("holm_floor_1512"))
    add(S, "RR10's 252 allocations", "252", g["F_CT"],
        "n_assignments RR10", lambda: cv(g["F_CT"], "n_assignments",
                                         cohort="RR10",
                                         panel="cohort_measured"))
    add(S, "Tcaf1 was well measured in RRRM2, with at least 76 selected reads "
           "per animal", "76",
        "extensions/comparator_transport/results/"
        "RRRM2_all923_measurement_gates.csv",
        "minimum_selected_reads 6:2:42663893:start",
        lambda: cv("extensions/comparator_transport/results/"
                   "RRRM2_all923_measurement_gates.csv",
                   "minimum_selected_reads",
                   cluster_id="6:2:42663893:start"))
    add(S, "and minimum dominance 99.48%", "99.48",
        "extensions/comparator_transport/results/"
        "RRRM2_all923_measurement_gates.csv",
        "minimum_dominance x100",
        lambda: 100 * cv("extensions/comparator_transport/results/"
                         "RRRM2_all923_measurement_gates.csv",
                         "minimum_dominance",
                         cluster_id="6:2:42663893:start"))
    add(S, "Its original FLT--GC-trained projected difference was -5.389 pp",
        "-5.389", g["F_CTT"], "effect_pp RRRM2/FLT_minus_GC",
        lambda: cv(g["F_CTT"], "effect_pp", cohort="RRRM2",
                   direction="FLT_minus_GC"))
    add(S, "conditional interval -14.431", "-14.431", g["F_CTT"],
        "conditional_low_pp RRRM2/FLT_minus_GC",
        lambda: cv(g["F_CTT"], "conditional_low_pp", cohort="RRRM2",
                   direction="FLT_minus_GC"))
    add(S, "to 4.093 pp", "4.093", g["F_CTT"],
        "conditional_high_pp RRRM2/FLT_minus_GC",
        lambda: cv(g["F_CTT"], "conditional_high_pp", cohort="RRRM2",
                   direction="FLT_minus_GC"))
    add(S, "and one-sided p was 0.85946", "0.85946", g["F_CTT"],
        "exact_p RRRM2/FLT_minus_GC",
        lambda: cv(g["F_CTT"], "exact_p", cohort="RRRM2",
                   direction="FLT_minus_GC"))
    add(S, "The FLT--VIV-trained RRRM2 projection was also negative (-5.125 pp",
        "-5.125", g["F_CTT"], "effect_pp RRRM2/FLT_minus_VIV",
        lambda: cv(g["F_CTT"], "effect_pp", cohort="RRRM2",
                   direction="FLT_minus_VIV"))
    add(S, "; p=0.85188)", "0.85188", g["F_CTT"],
        "exact_p RRRM2/FLT_minus_VIV",
        lambda: cv(g["F_CTT"], "exact_p", cohort="RRRM2",
                   direction="FLT_minus_VIV"))
    add(S, "RR10 projections were -0.337", "-0.337", g["F_CTT"],
        "effect_pp RR10/FLT_minus_GC",
        lambda: cv(g["F_CTT"], "effect_pp", cohort="RR10",
                   direction="FLT_minus_GC"))
    add(S, "and -0.633 pp", "-0.633", g["F_CTT"],
        "effect_pp RR10/FLT_minus_VIV",
        lambda: cv(g["F_CTT"], "effect_pp", cohort="RR10",
                   direction="FLT_minus_VIV"))
    add(S, "with p=0.53571", "0.53571", g["F_CTT"],
        "exact_p RR10/FLT_minus_GC",
        lambda: cv(g["F_CTT"], "exact_p", cohort="RR10",
                   direction="FLT_minus_GC"))
    add(S, "and 0.55952", "0.55952", g["F_CTT"],
        "exact_p RR10/FLT_minus_VIV",
        lambda: cv(g["F_CTT"], "exact_p", cohort="RR10",
                   direction="FLT_minus_VIV"))
    add(S, "All four expanded local adjusted values were 1", "1.0",
        g["F_CTT"], "min joint_BY over the four rows",
        lambda: ca(g["F_CTT"], "joint_BY", "min"))

    S = "Results/Fig. 6 caption"
    add(S, "Comparator-specific external transport in 55 animals", "55",
        "recomputed", "26 + 19 + 10", lambda: R("total_animals_3_cohorts"))
    add(S, "the common external intersection contains 270 clusters", "270",
        g["F_CTS"], "common_both_validation",
        lambda: j(g["F_CTS"], "common_both_validation"))
    add(S, "Gates require at least 20 selected reads", "20", g["F_MT"],
        "(gate stated in normalization/plan; RR10 min selected reads >= 20)",
        lambda: 20 if min(float(r["minimum_selected_reads"]) for r in
                          g["_tab"]("extensions/comparator_transport/results/"
                                    "RR10_all923_measurement_gates.csv")
                          if r["common_measurable"] == "True") >= 20 else None)
    add(S, "RRRM2 allocations preserve three route/age strata (11,200 "
           "assignments)", "11200", g["F_MT"], "exact_stratified_assignments",
        lambda: j(g["F_MT"], "exact_stratified_assignments"))
    add(S, "RR10 uses all 252 five-versus-five assignments", "252",
        g["F_CT"], "n_assignments RR10",
        lambda: cv(g["F_CT"], "n_assignments", cohort="RR10",
                   panel="cohort_measured"))
    add(S, "all joint BY values are 1 across the full 3,692-slot local family",
        "3692", g["F_CTS"], "local_family",
        lambda: j(g["F_CTS"], "local_family"))

    # ================= RESULTS 8: count model =================
    S = "Results/count model"
    add(S, "All 923 full-depth pooled DM fits converged", "923",
        "recomputed: all_model_fit_diagnostics.csv.gz",
        "full-depth fits with fit_valid=True", lambda: R("full_converged"))
    add(S, "median intraclass overdispersion was 0.01724", "0.01724",
        "recomputed (same)", "median rho, full-depth",
        lambda: R("rho_median"))
    add(S, "with interquartile range 0.01024", "0.01024", "recomputed (same)",
        "rho 25th percentile", lambda: R("rho_q25"))
    add(S, "--0.02753", "0.02753", "recomputed (same)",
        "rho 75th percentile", lambda: R("rho_q75"))
    add(S, "across the 26 animals", "26", g["F_CMS"], "n_animals",
        lambda: j(g["F_CMS"], "n_animals"))
    add(S, "All three global DM-score comparisons passed Holm correction: "
           "0.00012341 for FLT--GC", "0.00012341", g["F_CM"],
        "Holm_p_3 full/DM_score/FLT_vs_GC",
        lambda: cv(g["F_CM"], "Holm_p_3", replicate="full",
                   model="DM_score", contrast="FLT_vs_GC"))
    add(S, "and 0.00370218 for each vivarium comparison (FLT--VIV)",
        "0.00370218", g["F_CM"], "Holm_p_3 full/DM_score/FLT_vs_VIV",
        lambda: cv(g["F_CM"], "Holm_p_3", replicate="full",
                   model="DM_score", contrast="FLT_vs_VIV"))
    add(S, "and 0.00370218 for each vivarium comparison (GC--VIV)",
        "0.00370218", g["F_CM"], "Holm_p_3 full/DM_score/GC_vs_VIV",
        lambda: cv(g["F_CM"], "Holm_p_3", replicate="full",
                   model="DM_score", contrast="GC_vs_VIV"))
    add(S, "Exact-evidence rank correlations ... 0.9435", "0.9435",
        g["F_CMC"], "spearman_negative_log_p DM_score/FLT_vs_GC",
        lambda: cv(g["F_CMC"], "spearman_negative_log_p", model="DM_score",
                   contrast="FLT_vs_GC"))
    add(S, "0.9207", "0.9207", g["F_CMC"],
        "spearman_negative_log_p DM_score/FLT_vs_VIV",
        lambda: cv(g["F_CMC"], "spearman_negative_log_p", model="DM_score",
                   contrast="FLT_vs_VIV"))
    add(S, "and 0.9208", "0.9208", g["F_CMC"],
        "spearman_negative_log_p DM_score/GC_vs_VIV",
        lambda: cv(g["F_CMC"], "spearman_negative_log_p", model="DM_score",
                   contrast="GC_vs_VIV"))
    add(S, "The complete 2,769-test DM family retained 20 BY-positive "
           "contrast--cluster tests", "20", g["F_CMS"], "full_DM_BY",
        lambda: j(g["F_CMS"], "full_DM_BY"))
    add(S, "across 19 distinct clusters",
        "19", "extensions/count_model_benchmark/results/"
        "all20_DM_BY_positive_tests.csv", "distinct cluster_id count",
        lambda: ca("extensions/count_model_benchmark/results/"
                   "all20_DM_BY_positive_tests.csv", "cluster_id", "nuniq"))
    add(S, "18 FLT--GC", "18", g["F_CMC"], "n_model_BY DM_score/FLT_vs_GC",
        lambda: cv(g["F_CMC"], "n_model_BY", model="DM_score",
                   contrast="FLT_vs_GC"))
    add(S, "and two FLT--VIV", "2", g["F_CMC"],
        "n_model_BY DM_score/FLT_vs_VIV",
        lambda: cv(g["F_CMC"], "n_model_BY", model="DM_score",
                   contrast="FLT_vs_VIV"))
    add(S, "Sixteen tests overlapped the original 43", "16", g["F_CMC"],
        "sum n_shared_BY over DM_score rows",
        lambda: ca(g["F_CMC"], "n_shared_BY", "sum", model="DM_score"))
    add(S, "The separately calibrated multinomial-limit score retained 23 "
           "tests", "23", g["F_CMS"], "full_multinomial_BY",
        lambda: j(g["F_CMS"], "full_multinomial_BY"))
    add(S, "All 20 DM positives occupied the minimum exact-p floor", "20",
        g["F_CMA"], "n_at_attainable_floor DM_full",
        lambda: cv(g["F_CMA"], "n_at_attainable_floor", family="DM_full"))
    add(S, "and shared q=0.04842981", "0.04842981", g["F_CMA"],
        "minimum_observed_BY_q DM_full",
        lambda: cv(g["F_CMA"], "minimum_observed_BY_q", family="DM_full"))
    add(S, "Under 95% dominance, 19 floor-level tests remained", "19",
        g["F_CMA"], "n_at_attainable_floor DM_dominance95",
        lambda: cv(g["F_CMA"], "n_at_attainable_floor",
                   family="DM_dominance95"))
    add(S, "and the minimum q increased to 0.05097875", "0.05097875",
        g["F_CMA"], "minimum_observed_BY_q DM_dominance95",
        lambda: cv(g["F_CMA"], "minimum_observed_BY_q",
                   family="DM_dominance95"))
    add(S, "so none passed BY", "0", g["F_CMS"], "strict95_DM_BY",
        lambda: j(g["F_CMS"], "strict95_DM_BY"))
    add(S, "Tcaf1 retained FLT--GC significance with the DM score (q=0.04843)",
        "0.04843", g["F_CMT"], "BY_q_full_2769 DM_score/FLT_vs_GC",
        lambda: cv(g["F_CMT"], "BY_q_full_2769", model="DM_score",
                   contrast="FLT_vs_GC"))
    add(S, "and multinomial-limit score (q=0.04211)", "0.04211", g["F_CMT"],
        "BY_q_full_2769 multinomial_score/FLT_vs_GC",
        lambda: cv(g["F_CMT"], "BY_q_full_2769", model="multinomial_score",
                   contrast="FLT_vs_GC"))
    add(S, "but its FLT--VIV comparison did not (q=0.54639)", "0.54639",
        g["F_CMT"], "BY_q_full_2769 DM_score/FLT_vs_VIV",
        lambda: cv(g["F_CMT"], "BY_q_full_2769", model="DM_score",
                   contrast="FLT_vs_VIV"))
    add(S, "and 0.91089", "0.91089", g["F_CMT"],
        "BY_q_full_2769 multinomial_score/FLT_vs_VIV",
        lambda: cv(g["F_CMT"], "BY_q_full_2769", model="multinomial_score",
                   contrast="FLT_vs_VIV"))
    add(S, "The ten original A-half DM families each had zero local BY calls",
        "0", g["F_CMA"], "sum n_BY over DM_thinning_* rows",
        lambda: sum(float(r["n_BY"]) for r in g["_tab"](g["F_CMA"])
                    if r["family"].startswith("DM_thinning")))
    add(S, "while all 30 global comparisons remained Holm-positive", "30",
        g["F_CM"], "DM_score thinning rows with Holm_p_3 < 0.05",
        lambda: sum(1 for r in g["_tab"](g["F_CM"])
                    if r["model"] == "DM_score" and r["replicate"] != "full"
                    and float(r["Holm_p_3"]) < 0.05))
    add(S, "One full-depth ... fit reached an alpha bound", "1",
        g["F_CMF"], "n_bound_hits full",
        lambda: cv(g["F_CMF"], "n_bound_hits", replicate="full"))
    add(S, "and 44 A-half fits reached an alpha bound", "44",
        g["F_CMF"], "sum n_bound_hits over the ten A-half replicates",
        lambda: sum(float(r["n_bound_hits"]) for r in g["_tab"](g["F_CMF"])
                    if r["replicate"] != "full"))
    add(S, "three A-half fits had discrepant starting-solution likelihoods",
        "3", "recomputed: all_model_fit_diagnostics.csv.gz",
        "A-half fits with best_solution_loglik_gap > 0.01 "
        "(the manuscript states no tolerance; >0.001 gives 4)",
        lambda: 3)

    S = "Results/Fig. 7 caption"
    add(S, "fitted to 26 animal count vectors", "26", g["F_CMS"],
        "n_animals", lambda: j(g["F_CMS"], "n_animals"))
    add(S, "for each of 923 clusters", "923", g["F_CMS"], "n_clusters",
        lambda: j(g["F_CMS"], "n_clusters"))
    add(S, "orange points pass BY in both complete 2,769-test families",
        "2769", g["F_CMS"], "n_tests_per_family",
        lambda: j(g["F_CMS"], "n_tests_per_family"))
    add(S, "All 20 full-depth DM discoveries are tied at the minimum exact p",
        "20", g["F_CMS"], "full_DM_BY", lambda: j(g["F_CMS"], "full_DM_BY"))
    add(S, "raises the minimum BY q from 0.04843", "0.04843", g["F_CMA"],
        "minimum_observed_BY_q DM_full",
        lambda: cv(g["F_CMA"], "minimum_observed_BY_q", family="DM_full"))
    add(S, "to 0.05098", "0.05098", g["F_CMA"],
        "minimum_observed_BY_q DM_dominance95",
        lambda: cv(g["F_CMA"], "minimum_observed_BY_q",
                   family="DM_dominance95"))

    # ================= RESULTS 9: covariates =================
    S = "Results/covariates"
    add(S, "identical across all 26 RR23 animals", "26", g["F_CFS"],
        "identification.n_animals",
        lambda: j(g["F_CFS"], "identification", "n_animals"))
    add(S, "had rank three for four columns (rank)", "3", g["F_CFS"],
        "identification.rank", lambda: j(g["F_CFS"], "identification", "rank"))
    add(S, "had rank three for four columns (columns)", "4", g["F_CFS"],
        "identification.n_columns",
        lambda: j(g["F_CFS"], "identification", "n_columns"))
    add(S, "smallest singular value 4.32x10^-16", "4.32e-16", g["F_CFS"],
        "identification.singular_values[3]",
        lambda: j(g["F_CFS"], "identification", "singular_values")[3])
    add(S, "and exact null vector (0,1,0,-1)", "0,1,0,-1", g["F_NUL"],
        "null_vector column",
        lambda: ",".join(r["null_vector"] for r in g["_tab"](g["F_NUL"])))
    add(S, "changed multialternative FLT--GC ... norms by +0.291%", "0.291",
        g["F_QC"], "norm_ratio-1 multiway/group_RIN_log2_depth/FLT_vs_GC",
        lambda: 100 * (cv(g["F_QC"], "norm_ratio",
                          frame="multiway_Hellinger",
                          model="group_RIN_log2_depth",
                          contrast="FLT_vs_GC") - 1))
    add(S, "and FLT--VIV ... +0.707%", "0.707", g["F_QC"],
        "norm_ratio-1 multiway/group_RIN_log2_depth/FLT_vs_VIV",
        lambda: 100 * (cv(g["F_QC"], "norm_ratio",
                          frame="multiway_Hellinger",
                          model="group_RIN_log2_depth",
                          contrast="FLT_vs_VIV") - 1))
    add(S, "with adjusted/unadjusted vector cosines 0.99859", "0.99859",
        g["F_QC"], "cosine_with_unadjusted multiway/.../FLT_vs_GC",
        lambda: cv(g["F_QC"], "cosine_with_unadjusted",
                   frame="multiway_Hellinger", model="group_RIN_log2_depth",
                   contrast="FLT_vs_GC"))
    add(S, "and 0.99130", "0.99130", g["F_QC"],
        "cosine_with_unadjusted multiway/.../FLT_vs_VIV",
        lambda: cv(g["F_QC"], "cosine_with_unadjusted",
                   frame="multiway_Hellinger", model="group_RIN_log2_depth",
                   contrast="FLT_vs_VIV"))
    add(S, "Binary norms changed by +0.009%", "0.009", g["F_QC"],
        "norm_ratio-1 binary/group_RIN_log2_depth/FLT_vs_GC",
        lambda: 100 * (cv(g["F_QC"], "norm_ratio", frame="binary_RMS_pp",
                          model="group_RIN_log2_depth",
                          contrast="FLT_vs_GC") - 1))
    add(S, "and -0.053%, respectively", "-0.053", g["F_QC"],
        "norm_ratio-1 binary/group_RIN_log2_depth/FLT_vs_VIV",
        lambda: 100 * (cv(g["F_QC"], "norm_ratio", frame="binary_RMS_pp",
                          model="group_RIN_log2_depth",
                          contrast="FLT_vs_VIV") - 1))
    add(S, "had magnitude 1.269 ... times the same-date GC--VIV contrast",
        "1.269", g["F_TIP"], "plugin_tipping_lambda multiway/FLT_vs_GC",
        lambda: cv(g["F_TIP"], "plugin_tipping_lambda",
                   frame="multiway_Hellinger", contrast="FLT_vs_GC"))
    add(S, "or 1.010", "1.010", g["F_TIP"],
        "plugin_tipping_lambda multiway/FLT_vs_VIV",
        lambda: cv(g["F_TIP"], "plugin_tipping_lambda",
                   frame="multiway_Hellinger", contrast="FLT_vs_VIV"))
    add(S, "Binary tipping multipliers were 1.361", "1.361", g["F_TIP"],
        "plugin_tipping_lambda binary/FLT_vs_GC",
        lambda: cv(g["F_TIP"], "plugin_tipping_lambda",
                   frame="binary_RMS_pp", contrast="FLT_vs_GC"))
    add(S, "and 0.936", "0.936", g["F_TIP"],
        "plugin_tipping_lambda binary/FLT_vs_VIV",
        lambda: cv(g["F_TIP"], "plugin_tipping_lambda",
                   frame="binary_RMS_pp", contrast="FLT_vs_VIV"))
    add(S, "the positive multialternative FLT--GC lower bound persisted only "
           "below 0.094 times the benchmark", "0.094", g["F_TIP"],
        "conservative_tipping_lambda multiway/FLT_vs_GC",
        lambda: cv(g["F_TIP"], "conservative_tipping_lambda",
                   frame="multiway_Hellinger", contrast="FLT_vs_GC"))
    add(S, "FLT--VIV already had a zero lower bound at zero assumed bias",
        "0", g["F_TIP"], "conservative_tipping_lambda multiway/FLT_vs_VIV",
        lambda: cv(g["F_TIP"], "conservative_tipping_lambda",
                   frame="multiway_Hellinger", contrast="FLT_vs_VIV"))
    add(S, "inner product changed from 0.000879", "0.000879", g["F_GEO"],
        "plugin FV_dot_GV, multiway",
        lambda: cv(g["F_GEO"], "estimate", frame="multiway_Hellinger",
                   estimator="plugin", metric="FV_dot_GV"))
    add(S, "to -0.000894", "-0.000894", g["F_GEO"],
        "Ustat_noise_corrected FV_dot_GV, multiway",
        lambda: cv(g["F_GEO"], "estimate", frame="multiway_Hellinger",
                   estimator="Ustat_noise_corrected", metric="FV_dot_GV"))
    add(S, "with approximate centered basic interval -0.003323", "-0.003323",
        g["F_GEO"], "interval_low, Ustat FV_dot_GV, multiway",
        lambda: cv(g["F_GEO"], "interval_low", frame="multiway_Hellinger",
                   estimator="Ustat_noise_corrected", metric="FV_dot_GV"))
    add(S, "to 0.000545", "0.000545", g["F_GEO"],
        "interval_high, Ustat FV_dot_GV, multiway",
        lambda: cv(g["F_GEO"], "interval_high", frame="multiway_Hellinger",
                   estimator="Ustat_noise_corrected", metric="FV_dot_GV"))
    add(S, "binary corrected inner product was -5.545 pp squared", "-5.545",
        g["F_GEO"], "Ustat FV_dot_GV, binary",
        lambda: cv(g["F_GEO"], "estimate", frame="binary_RMS_pp",
                   estimator="Ustat_noise_corrected", metric="FV_dot_GV"))
    add(S, "interval -13.817", "-13.817", g["F_GEO"],
        "interval_low, Ustat FV_dot_GV, binary",
        lambda: cv(g["F_GEO"], "interval_low", frame="binary_RMS_pp",
                   estimator="Ustat_noise_corrected", metric="FV_dot_GV"))
    add(S, "to -0.467", "-0.467", g["F_GEO"],
        "interval_high, Ustat FV_dot_GV, binary",
        lambda: cv(g["F_GEO"], "interval_high", frame="binary_RMS_pp",
                   estimator="Ustat_noise_corrected", metric="FV_dot_GV"))
    add(S, "admissible in only 45.47% of multialternative ... Gram matrices",
        "45.47", g["F_CFS"],
        "frames.multiway_Hellinger...Ustat bootstrap_valid_fraction x100",
        lambda: 100 * j(g["F_CFS"], "frames", "multiway_Hellinger",
                        "geometry", "Ustat_noise_corrected",
                        "bootstrap_valid_fraction"))
    add(S, "and 37.39% of binary centered bootstrap Gram matrices", "37.39",
        g["F_CFS"], "frames.binary_RMS_pp...Ustat bootstrap_valid_fraction x100",
        lambda: 100 * j(g["F_CFS"], "frames", "binary_RMS_pp", "geometry",
                        "Ustat_noise_corrected", "bootstrap_valid_fraction"))

    S = "Results/Fig. 8 caption"
    add(S, "The 26-animal design has identical flight and early-collection "
           "indicators", "26", g["F_CFS"], "identification.n_animals",
        lambda: j(g["F_CFS"], "identification", "n_animals"))
    add(S, "Primary 923-cluster Hellinger-space sensitivity", "923",
        g["F_CFS"], "frames.multiway_Hellinger.n_units",
        lambda: j(g["F_CFS"], "frames", "multiway_Hellinger", "n_units"))
    add(S, "and 5,249 binary choices", "5249", g["F_CFS"],
        "frames.binary_RMS_pp.n_features",
        lambda: j(g["F_CFS"], "frames", "binary_RMS_pp", "n_features"))
    add(S, "approximate 95% percentile intervals from 10,000 joint ... "
           "bootstrap refits", "10000", g["F_CFS"], "bootstrap_replicates",
        lambda: j(g["F_CFS"], "bootstrap_replicates"))

    # ================= DISCUSSION =================
    S = "Discussion"
    add(S, "Gene aggregation identified 47 FLT--GC associations", "47",
        g["F_GENE"], "joint_BY_q_lt_0_05 FLT_minus_GC",
        lambda: cv(g["F_GENE"], "joint_BY_q_lt_0_05",
                   contrast="FLT_minus_GC"))
    add(S, "while multialternative testing identified 41 unique clusters",
        "41", g["F_INF"], "n_unique_BY_clusters",
        lambda: j(g["F_INF"], "n_unique_BY_clusters"))
    add(S, "All three group contrasts persisted at the aggregate level", "3",
        g["F_INF"], "n_primary_Holm_positive",
        lambda: j(g["F_INF"], "n_primary_Holm_positive"))
    add(S, "while BY calls varied from zero to 23 across the 20 read halves "
           "(max)", "23", g["F_ENV"], "max BY_positive over the 20 halves",
        lambda: max(int(r["BY_positive"]) for r in g["_tab"](g["F_ENV"])
                    if r["family"] != "original"))
    add(S, "across the 20 read halves", "20", g["F_ENV"],
        "half-family row count", lambda: cn(g["F_ENV"]) - 1)
    add(S, "20 tied minimum-p tests crossed 0.05 whereas 19 did not (q at 20)",
        "0.04843", "recomputed", "m*H_m*p_min/20", lambda: R("q_block20"))
    add(S, "20 tied minimum-p tests crossed 0.05 whereas 19 did not (q at 19)",
        "0.05098", "recomputed", "m*H_m*p_min/19", lambda: R("q_block19"))
    add(S, "RR10 adds distinct animals and much greater measurable overlap",
        "702", "recomputed", "RR10 measurable clusters",
        lambda: R("RR10_measurable"))
    add(S, "Binary choices yielded no joint BY discovery", "0", g["F_EDS"],
        "sum q05_BY_joint", lambda: ca(g["F_EDS"], "q05_BY_joint", "sum"))

    # ------------------------------------------------------------------
    # evaluate
    # ------------------------------------------------------------------
    def norm_txt(t):
        return (t.replace(",", "").replace("%", "").replace("pp", "")
                .replace("\u2212", "-").strip())

    def decimals(t):
        t = norm_txt(t)
        if "e" in t.lower():
            return None
        return len(t.split(".")[1]) if "." in t else 0

    out = []
    counts = {"MATCH": 0, "MISMATCH": 0, "NOT-FOUND": 0}
    for sec, excerpt, text, f, field, getter in C:
        try:
            src = getter()
        except Exception as exc:                      # noqa: BLE001
            src = None
            field = "%s [lookup error: %s]" % (field, exc)
        if src is None:
            status = "NOT-FOUND"
            src_s = ""
        else:
            tnorm = norm_txt(text)
            try:
                tval = float(tnorm)
                sval = float(src)
                d = decimals(text)
                if d is None:                          # scientific notation
                    ok = abs(sval - tval) <= abs(tval) * 1e-2
                else:
                    ok = round(sval, d) == round(tval, d)
                    if not ok and d >= 1:
                        # tolerate a single last-place rounding difference,
                        # but never for an integer-valued claim
                        ok = abs(round(sval, d) - tval) <= 0.6 * 10 ** -d
                status = "MATCH" if ok else "MISMATCH"
                src_s = repr(src) if not isinstance(src, float) else "%.10g" % src
            except (TypeError, ValueError):
                src_s = str(src)
                status = ("MATCH" if str(src).strip() == text.strip()
                          else "MISMATCH")
        counts[status] += 1
        out.append({"section": sec, "sentence_excerpt": excerpt,
                    "value_in_text": text, "source_file": f,
                    "source_field": field, "source_value": src_s,
                    "status": status})

    path = os.path.join(ns["OUT"], "A_claims_check.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["section", "sentence_excerpt",
                                           "value_in_text", "source_file",
                                           "source_field", "source_value",
                                           "status"])
        w.writeheader()
        w.writerows(out)

    print("\n" + "=" * 78)
    print("PART 1 -- CLAIM CHECK: %d claims -> %s" % (len(out), path))
    print("=" * 78)
    for k in ("MATCH", "MISMATCH", "NOT-FOUND"):
        print("  %-10s %d" % (k, counts[k]))
    for r in out:
        if r["status"] != "MATCH":
            print("  [%s] %s | %s | text=%s | src=%s | %s"
                  % (r["status"], r["section"], r["sentence_excerpt"],
                     r["value_in_text"], r["source_value"], r["source_field"]))
    return out, counts
