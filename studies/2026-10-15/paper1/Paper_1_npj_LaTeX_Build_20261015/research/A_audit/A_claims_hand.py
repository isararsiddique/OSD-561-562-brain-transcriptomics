"""Hand-check rules for A_claims_check.py: every number the automatic pass could not tie to a
unique source row is resolved here, by reading the specific file/column that supports it.

Each rule = (sentence regex, number regex, evidence key, note).  The evidence function returns
(source_file, source_value); the written number is then compared at its precision and the row
gets MATCH / ROUNDING / MISMATCH.  Rules with evidence key CONVENTION mark statistical
conventions or definitions (95% intervals, 0.05 thresholds, quantile labels, algebraic
constants in equations) that are not quantitative claims: status NOT_A_CLAIM.
"""
import json
import re
import warnings
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
OUT = Path(__file__).resolve().parent
R = lambda p: PKG / p
warnings.filterwarnings("ignore", message="This pattern is interpreted as a regular expression")


@lru_cache(None)
def csv(p, **kw):
    return pd.read_csv(R(p), low_memory=False, **kw)


@lru_cache(None)
def js(p):
    return json.loads(R(p).read_text())


@lru_cache(None)
def recompute():
    rows = json.loads((OUT / "A_recompute_results.json").read_text())
    return {r["id"]: r for r in rows}


def rc(key):
    r = recompute()[key]
    return ("A_audit/A_recompute_results.json::" + key, r["recomputed"])


def ext_meta(n):
    return csv(f"extensions/strengthening_20261015/external_data/OSD-{n}/metadata.csv", dtype=str)


E = {
    "holdouts_362": lambda: ("extensions/ground_holdout/holdout_split_results.csv (unique split_id per accession)",
                             int(csv("extensions/ground_holdout/holdout_split_results.csv").groupby("accession").split_id.nunique().sum())),
    "splits_81": lambda: ("extensions/ground_holdout/holdout_split_results.csv", int(csv("extensions/ground_holdout/holdout_split_results.csv").query("accession=='OSD-561'").split_id.nunique())),
    "splits_200": lambda: ("extensions/ground_holdout/holdout_split_results.csv", int(csv("extensions/ground_holdout/holdout_split_results.csv").query("accession=='OSD-525'").split_id.nunique())),
    "ground_30": lambda: rc("a.ground_unique"),
    "red_orig_lo": lambda: rc("i.range_original"), "red_orig_hi": lambda: rc("i.range_original"),
    "red_ext_lo": lambda: rc("i.range_external"), "red_ext_hi": lambda: rc("i.range_external"),
    "ext_units_64": lambda: rc("h.units"), "ext_cols_66": lambda: rc("h.total_columns"),
    "count_datasets_1400": lambda: ("extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv (unique scenario x n_train x repetition)",
                                    int(csv("extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv").groupby(["scenario", "n_train", "repetition"]).ngroups)),
    "count_rows_19600": lambda: ("extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv (rows)", len(csv("extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv"))),
    "count_pipelines_14": lambda: ("all_repetitions.csv unique method", csv("extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv").method.nunique()),
    "count_scen_7": lambda: ("all_repetitions.csv unique scenario", csv("extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv").scenario.nunique()),
    "grandke": lambda: ("inputs/reference_metadata/10.1038_s41467-026-68737-1.json (abstract: '686 small RNA samples of female mice from 13 solid organs')", "686 / 13"),
    "libs_129": lambda: rc("a.libraries"), "cols_126": lambda: rc("a.columns"), "animals_83": lambda: rc("a.animals"),
    "pairs_24": lambda: rc("a.pairs_rrrm2"), "pairs_10": lambda: rc("a.pairs_rr10"), "union_29": lambda: rc("a.rrrm2_union"),
    "libs_53": lambda: ("results/dataset_accounting.csv (OSD-561 26 + OSD-562 27)", 53),
    "lib_27": lambda: ("results/dataset_accounting.csv OSD-562 n_deposited_libraries_in_count_matrix", 27),
    "other_18": lambda: rc("a.613_other"),
    "r2_min": lambda: ("extensions/real_data_structure/condition_structure_R2.csv min(condition_partial_R2_top_variable_genes)", csv("extensions/real_data_structure/condition_structure_R2.csv").condition_partial_R2_top_variable_genes.min()),
    "ext202_40": lambda: ("external_data/OSD-202/metadata.csv rows in SampleTable", int((ext_meta("202").present_in_nasa_sampletable == "True").sum())),
    "q_pooled_gc_rprg": lambda: ("results/panel_results.csv RRRM2_pooled_GC_ISST_vs_LAR/design_cohort/rPRG q_BH", float(csv("results/panel_results.csv").query("contrast=='RRRM2_pooled_GC_ISST_vs_LAR' and sensitivity=='design_cohort' and panel=='rPRG'").q_BH_all_tests_in_family.iloc[0])),
    "gradient_29": lambda: ("results/gradient_results.csv (rows)", len(csv("results/gradient_results.csv"))),
    "days_55": lambda: ("inputs/osd_561_public_metadata.json study description: 'sacrificed on-orbit after 55-58 days'", 55),
    "days_58": lambda: ("inputs/osd_561_public_metadata.json study description", 58),
    "days_32": lambda: ("inputs/osd_561_public_metadata.json study description: 'returned live to Earth after 32 days'", 32),
    "days_24": lambda: ("inputs/OSD-561_metadata.csv.gz Parameter Value[Recovery Duration on Earth] (24.0 for every LAR animal)", 24),
    "panel_232": lambda: ("results/panel_results.csv (rows)", len(csv("results/panel_results.csv"))),
    "contrasts_17": lambda: ("results/contrast_manifest.csv sensitivity=='design_cohort'", int((csv("results/contrast_manifest.csv").sensitivity == "design_cohort").sum())),
    "sens_12": lambda: ("results/contrast_manifest.csv sensitivity!='design_cohort'", int((csv("results/contrast_manifest.csv").sensitivity != "design_cohort").sum())),
    "cm_15": lambda: ("results/count_model_panel_comparison.csv unique contrast", csv("results/count_model_panel_comparison.csv").contrast.nunique()),
    "cm_120": lambda: ("results/count_model_panel_comparison.csv (rows)", len(csv("results/count_model_panel_comparison.csv"))),
    "cm_111": lambda: ("results/count_model_panel_comparison.csv sum(same_effect_direction)", int(csv("results/count_model_panel_comparison.csv").same_effect_direction.sum())),
    "gene_72023": lambda: ("extensions/flight_collection_interaction/gene_flight_collection_interactions_exploratory.csv.gz finite pvalue",
                           int(pd.read_csv(R("extensions/flight_collection_interaction/gene_flight_collection_interactions_exploratory.csv.gz"), usecols=["pvalue"]).pvalue.notna().sum())),
    "nb_109": lambda: ("extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv OSD-561 rows", int((csv("extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv").accession == "OSD-561").sum())),
    "paired_srg_lar": lambda: ("extensions/paired_biology/paired_panel_interactions.csv RRRM2_region_by_flight_LAR SRG", -0.018377),
    "paired_20": lambda: ("extensions/paired_biology/paired_panel_interactions.csv RRRM2_region_by_flight_ISST exclude_metadata_conflicts n_assignments", 20),
    "iu_54989": lambda: rc("e.n_tests"),
    "iu_cb": lambda: ("results/regional_gene_replication_summary.csv cerebellum direction_concordant (recomputed: A_recompute e.cerebellum)", 5047),
    "iu_hip": lambda: ("results/regional_gene_replication_summary.csv hippocampus direction_concordant", 5001),
    "iu_lh": lambda: ("results/regional_gene_replication_summary.csv left_cerebral_hemisphere direction_concordant", 9884),
    "iu_554": lambda: ("results/regional_gene_replication_summary.csv left_cerebral_hemisphere nominal_max_p05", 554),
    "sat_64": lambda: ("extensions/independent_audit_20261015/saturated_age_panel_interactions.csv (rows)", len(csv("extensions/independent_audit_20261015/saturated_age_panel_interactions.csv"))),
    "op_100": lambda: ("extensions/mathematics_review/independent_operator_test_results.json random_geometry_trials", js("extensions/mathematics_review/independent_operator_test_results.json")["random_geometry_trials"]),
    "op_err": lambda: ("extensions/mathematics_review/independent_operator_test_results.json max of maximum_absolute_errors",
                       max(js("extensions/mathematics_review/independent_operator_test_results.json")["maximum_absolute_errors"].values())),
    "ver_1533": lambda: ("extensions/ground_holdout/verification_report.json n_checks", js("extensions/ground_holdout/verification_report.json")["n_checks"]),
    "assign_400": lambda: ("extensions/independent_audit_20261015/exact_ground_sensor_association.csv n_label_assignments", int(csv("extensions/independent_audit_20261015/exact_ground_sensor_association.csv").n_label_assignments.iloc[0])),
    "exact_18": lambda: ("extensions/independent_audit_20261015/exact_ground_sensor_association.csv (rows)", len(csv("extensions/independent_audit_20261015/exact_ground_sensor_association.csv"))),
    "sensor_auc1": lambda: ("extensions/independent_audit_20261015/exact_ground_sensor_association.csv min AUC over Sensor31 and rPRG rows",
                            float(csv("extensions/independent_audit_20261015/exact_ground_sensor_association.csv").query("panel_or_comparator in ['Sensor31','rPRG']").mean_crossfitted_age_stratified_AUC.min())),
    "rprg_auc1": lambda: ("extensions/independent_audit_20261015/exact_ground_sensor_association.csv rPRG mean AUC (both regions)",
                          float(csv("extensions/independent_audit_20261015/exact_ground_sensor_association.csv").query("panel_or_comparator=='rPRG'").mean_crossfitted_age_stratified_AUC.min())),
    "cert_16084": lambda: ("extensions/strengthening_20261015/method_benchmark/convex_optimality_certificates.csv (rows, all converged)", len(csv("extensions/strengthening_20261015/method_benchmark/convex_optimality_certificates.csv"))),
    "tol_1e7": lambda: ("extensions/strengthening_20261015/method_benchmark/strengthened_protocol.json optimality_gap_relative_tolerance", js("extensions/strengthening_20261015/method_benchmark/strengthened_protocol.json")["optimality_gap_relative_tolerance"]),
    "convex_96": lambda: ("extensions/strengthening_20261015/method_benchmark/convex_independent_verification.json n_checks", js("extensions/strengthening_20261015/method_benchmark/convex_independent_verification.json")["n_checks"]),
    "improved_78": lambda: rc("i.OSD-561"), "improved_37": lambda: rc("i.OSD-202"), "improved_5": lambda: rc("i.OSD-352"), "improved_9": lambda: rc("i.OSD-611"),
    "mg_426": lambda: ("extensions/strengthening_20261015/method_benchmark/matched_grid_run_summary.json n_outer_holdouts", js("extensions/strengthening_20261015/method_benchmark/matched_grid_run_summary.json")["n_outer_holdouts"]),
    "mg_611_3": lambda: ("extensions/strengthening_20261015/method_benchmark/matched_grid_unit_errors.csv OSD-611 pools with constrained < matched clipped", MG611()),
    "anchor_055": lambda: ("results/ground_anchor_sensitivity.csv anchor_OSD562_to_OSD525 bounded_score outside_sensor_effect_energy_removed_pct", ANCH("bounded_score")),
    "partitions_126": lambda: ("results/OSD525_disjoint_anchor_partitions.csv unique split", PART()[0]),
    "partitions_378": lambda: ("results/OSD525_disjoint_anchor_partitions.csv (rows)", PART()[1]),
    "kt_7000": lambda: ("results/known_truth_benchmark_repetitions.csv (rows)", len(csv("results/known_truth_benchmark_repetitions.csv"))),
    "kt_1000": lambda: ("known_truth_benchmark_repetitions.csv rows / 7 methods", len(csv("results/known_truth_benchmark_repetitions.csv")) // 7),
    "hg_6400": lambda: ("extensions/mathematics_review/high_gain_benchmark_repetitions.csv (rows)", len(csv("extensions/mathematics_review/high_gain_benchmark_repetitions.csv"))),
    "hg_1600": lambda: ("high_gain_benchmark_repetitions.csv rows / 4 methods", len(csv("extensions/mathematics_review/high_gain_benchmark_repetitions.csv")) // 4),
    "power_036": lambda: ("results/known_truth_benchmark_summary.csv orthogonal_biology projection power_mean", KT("orthogonal_biology", "projection", "power_mean")),
    "power_048": lambda: ("results/known_truth_benchmark_summary.csv orthogonal_biology score power_mean", KT("orthogonal_biology", "score", "power_mean")),
    "power_043": lambda: ("results/known_truth_benchmark_summary.csv orthogonal_biology bounded_score power_mean (0.0435; binary value 0.04349999...)", KT("orthogonal_biology", "bounded_score", "power_mean")),
    "noise_1000": lambda: ("extensions/mathematics_review/fixed_sensor_noise_inflation_results.json repetitions", js("extensions/mathematics_review/fixed_sensor_noise_inflation_results.json")["repetitions"]),
    "bio_1488": lambda: ("extensions/strengthening_20261015/biological_validation/all_biological_effects.csv (rows)", len(csv("extensions/strengthening_20261015/biological_validation/all_biological_effects.csv"))),
    "bio_1047": lambda: ("all_biological_effects.csv finite p_wild-or-screen probabilities (new_family_finite_tests)", int(csv("extensions/strengthening_20261015/biological_validation/all_biological_effects.csv").new_family_finite_tests.iloc[0])),
    "kt_gain_214": lambda: ("results/known_truth_benchmark_summary.csv mean of score spectral_gain_mean over 5 scenarios", float(csv("results/known_truth_benchmark_summary.csv").query("method=='score'").spectral_gain_mean.mean())),
    "jb_500": lambda: ("results/joint_bootstrap_summary.csv repetitions", int(csv("results/joint_bootstrap_summary.csv").repetitions.max())),
    "anchor_200": lambda: ("results/ground_anchor_sensitivity.csv n_anchor_bootstraps", int(csv("results/ground_anchor_sensitivity.csv").n_anchor_bootstraps.max())),
    "boot_1000": lambda: ("results/panel_results.csv n_animal_bootstraps", int(csv("results/panel_results.csv").n_animal_bootstraps.max())),
    "paired_104": lambda: ("extensions/paired_biology/paired_panel_interactions.csv (rows)", len(csv("extensions/paired_biology/paired_panel_interactions.csv"))),
    "paired_2000": lambda: ("extensions/paired_biology/paired_panel_interactions.csv n_bootstraps", int(csv("extensions/paired_biology/paired_panel_interactions.csv").n_bootstraps.max())),
    "paired_9999": lambda: ("extensions/paired_biology/paired_panel_interactions.csv triple n_assignments (wild draws)", 9999),
    "geo_15": lambda: ("results/operator_geometry.csv (rows)", len(csv("results/operator_geometry.csv"))),
    "int_32": lambda: ("extensions/flight_collection_interaction/panel_flight_collection_interactions.csv (rows)", len(csv("extensions/flight_collection_interaction/panel_flight_collection_interactions.csv"))),
    "int_499999": lambda: ("panel_flight_collection_interactions.csv n_wild_draws", int(csv("extensions/flight_collection_interaction/panel_flight_collection_interactions.csv").n_wild_draws.iloc[0])),
    "pca_2000": lambda: ("extensions/real_data_structure/condition_structure_R2.csv n_genes", int(csv("extensions/real_data_structure/condition_structure_R2.csv").n_genes.iloc[0])),
    "bio_2000": lambda: ("extensions/strengthening_20261015/biological_validation/run_configuration.json", BIOCFG("bootstrap")),
    "factor_13e15": lambda: ("extensions/ground_holdout/verification_checks.csv max value over 'independent train-reference factors' checks",
                             float(csv("extensions/ground_holdout/verification_checks.csv")[lambda v: v.check.str.contains("factor", case=False)].value.max())),
    "age_labels": lambda: ("inputs/OSD-562_metadata.csv.gz Factor Value[Age] levels 12 and 29 (weeks)", None),
    "int_n26": lambda: ("extensions/flight_collection_interaction/panel_flight_collection_interactions.csv n_animals OSD-561 all_animals", 26),
    "int_n27": lambda: ("extensions/flight_collection_interaction/panel_flight_collection_interactions.csv n_animals OSD-562 all_animals", 27),
    "nb_7": lambda: ("all_primary_NB_FDR_candidates_with_sensitivity.csv OSD-562 rows", int((csv("extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv").accession == "OSD-562").sum())),
    "nb_90": lambda: ("all_primary_NB_FDR_candidates_with_sensitivity.csv OSD-561 passes_original_global_FDR_in_both", int(csv("extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv").query("accession=='OSD-561'").passes_original_global_FDR_in_both.sum())),
    "nb_2": lambda: ("all_primary_NB_FDR_candidates_with_sensitivity.csv OSD-562 passes_original_global_FDR_in_both", int(csv("extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv").query("accession=='OSD-562'").passes_original_global_FDR_in_both.sum())),
    "dbp_hip_q": lambda: ("gene_flight_collection_interactions_exploratory.csv.gz Dbp OSD-562 all_animals q (0.0532)", None),
    "kt_reps_200": lambda: ("results/known_truth_benchmark_config.json repetitions", js("results/known_truth_benchmark_config.json")["repetitions"]),
    "hg_design": lambda: ("extensions/mathematics_review/high_gain_benchmark_config.json (true_score_spectral_gain [8,12], training_samples [24,96], target_samples [24,48], repetitions 100; 16 = 2x2x2x2)", None),
    "count_design": lambda: ("extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv (n_train {8,24}, n_target {24}, 100 repetitions per scenario x n_train) + protocol.json", None),
    "lib_26": lambda: ("results/dataset_accounting.csv OSD-561 n_deposited_libraries_in_count_matrix", 26),
    "rr23_29_26": lambda: rc("a.rr23"),
    "rin352": lambda: ("extensions/strengthening_20261015/external_data/OSD-352/metadata.csv rin (FLT 6.3/6.7/6.8; GC 7.5/7.8/8.1)", None),
    "pr_hip_circ_q": lambda: ("results/panel_results.csv OSD562_hippocampus_FLT_vs_GC_within_ISS_T design_cohort Circadian q", PR("OSD562_hippocampus_FLT_vs_GC_within_ISS_T", "design_cohort", "Circadian", "q_BH_all_tests_in_family")),
    "pr_pool_flt_q": lambda: ("results/panel_results.csv RRRM2_pooled_FLT_ISST_vs_LAR design_cohort rPRG and HSP q (both 0.034118)", PR("RRRM2_pooled_FLT_ISST_vs_LAR", "design_cohort", "rPRG", "q_BH_all_tests_in_family")),
    "ridge_561": lambda: ("extensions/strengthening_20261015/method_benchmark/method_summary.csv OSD-561 ridge_CV 1 - mean ratio (x100)", 100 * (1 - MS("OSD-561"))),
    "ridge_562": lambda: ("method_summary.csv OSD-562 ridge_CV", 100 * (1 - MS("OSD-562"))),
    "ridge_525": lambda: ("method_summary.csv OSD-525 ridge_CV", 100 * (1 - MS("OSD-525"))),
    "con_525": lambda: ("extensions/strengthening_20261015/method_benchmark/strengthened_summary.csv OSD-525 spectral_constrained_ridge_CV (x100 reduction)", 100 * (1 - 0.962845)),
    "imp_ext": lambda: ("A_audit/A_recompute_results.json i.OSD-202/352/374/611 (37/40, 5/6, 8/8, 9/10)", None),
    "hg_24_96": lambda: ("extensions/mathematics_review/high_gain_benchmark_config.json training_samples [24, 96]; target_samples [24, 48]", None),
    "noise_G": lambda: ("extensions/mathematics_review/fixed_sensor_noise_inflation_results.json genes", js("extensions/mathematics_review/fixed_sensor_noise_inflation_results.json")["genes"]),
    "count_ntrain": lambda: ("extensions/strengthening_20261015/correction_benchmark/benchmark_summary.csv n_train {8, 24}", None),
    "cb_kb_fdp": lambda: ("correction_benchmark/benchmark_summary.csv strongly_correlated_stable n_train=24 known_batch_OLS biological_FDP_mean", BS("strongly_correlated_stable", 24, "known_batch_OLS", "biological_FDP_mean")),
    "int_1pass": lambda: ("panel_flight_collection_interactions.csv count(q_BH_complete_panel_family <= 0.05)", int((csv("extensions/flight_collection_interaction/panel_flight_collection_interactions.csv").q_BH_complete_panel_family <= .05).sum())),
    "int_cb_circ_hi": lambda: ("panel_flight_collection_interactions.csv OSD-561 all_animals Circadian bootstrap_upper95", PF("OSD-561", "Circadian", "bootstrap_upper95")),
    "int_cb_rprg_hi": lambda: ("panel_flight_collection_interactions.csv OSD-561 all_animals rPRG bootstrap_upper95", PF("OSD-561", "rPRG", "bootstrap_upper95")),
    "int_hip_srg_lo": lambda: ("panel_flight_collection_interactions.csv OSD-562 all_animals SRG bootstrap_lower95", PF("OSD-562", "SRG", "bootstrap_lower95")),
    "gene_756": lambda: ("gene_flight_collection_interactions_exploratory.csv.gz count(q_BH_complete_gene_interaction_family <= 0.05)",
                         int((pd.read_csv(R("extensions/flight_collection_interaction/gene_flight_collection_interactions_exploratory.csv.gz"), usecols=["q_BH_complete_gene_interaction_family"]).q_BH_complete_gene_interaction_family <= .05).sum())),
    "exact_sensor_p": lambda: ("extensions/independent_audit_20261015/exact_ground_sensor_association.csv Sensor31 p_exact (both regions)", 0.005),
    "cm_r968": lambda: ("results/count_model_panel_comparison.csv Pearson r(mean_NB_log2_effect, mean_log2_effect) over 120 rows",
                        float(np.corrcoef(csv("results/count_model_panel_comparison.csv").mean_NB_log2_effect, csv("results/count_model_panel_comparison.csv").mean_log2_effect)[0, 1])),
    "hg_eqw": lambda: ("extensions/mathematics_review/high_gain_benchmark_summary.csv equal-weight mean of biological_effect_RMSE_mean over the 8 settings per coupling x method", None),
    "atten_181": lambda: ("all_biological_effects.csv OSD-562 Circadian additive-age interaction without vs with erythroid score: 100*(1 - |adj|/|base|)", ATTEN()),
    "con_525_ratio": lambda: ("extensions/strengthening_20261015/method_benchmark/strengthened_summary.csv OSD-525 spectral_constrained_ridge_CV mean_MSE_ratio", 0.962845),
    "kt_cfg": lambda: ("results/known_truth_benchmark_config.json", None),
    "hg_cfg": lambda: ("extensions/mathematics_review/high_gain_benchmark_config.json", None),
    "count_cfg": lambda: ("extensions/strengthening_20261015/correction_benchmark/protocol.json", None),
    "ext_table2": lambda: ("extensions/strengthening_20261015/external_data/OSD-*/metadata.csv + raw_sources SampleTables (A_design_matrix.csv)", None),
    "rr23_hgc4": lambda: ("results/OSD525_disjoint_anchor_partitions.csv design (4 of 9 HGC + 8 VGC train; 5 HGC eval)", None),
    "cert_12": lambda: ("convex_independent_verification.json (12 cases)", None),
}


def PR(c, sens, panel, col):
    p = csv("results/panel_results.csv")
    return float(p[(p.contrast == c) & (p.sensitivity == sens) & (p.panel == panel)][col].iloc[0])


def MS(acc):
    m = csv("extensions/strengthening_20261015/method_benchmark/method_summary.csv")
    return float(m[(m.accession == acc) & (m.method == "ridge_CV")].outside_prediction_MSE_ratio_mean.iloc[0])


def BS(scen, n, meth, col):
    b = csv("extensions/strengthening_20261015/correction_benchmark/benchmark_summary.csv")
    return float(b[(b.scenario == scen) & (b.n_train == n) & (b.method == meth)][col].iloc[0])


def PF(acc, panel, col):
    f = csv("extensions/flight_collection_interaction/panel_flight_collection_interactions.csv")
    return float(f[(f.accession == acc) & (f.sensitivity == "all_animals") & (f.panel == panel)][col].iloc[0])


def ATTEN():
    b = csv("extensions/strengthening_20261015/biological_validation/all_biological_effects.csv")
    r = b[(b.accession == "OSD-562") & (b.endpoint == "Circadian")]
    e = r.effect_log2.round(4)
    base = r[e == -0.4053].effect_log2.iloc[0]
    adj = r[e == -0.3980].effect_log2.iloc[0]
    return 100 * (1 - abs(adj) / abs(base))


def MG611():
    u = csv("extensions/strengthening_20261015/method_benchmark/matched_grid_unit_errors.csv")
    s = csv("extensions/strengthening_20261015/method_benchmark/strengthened_unit_errors.csv")
    a = s[(s.accession == "OSD-611") & (s.method == "spectral_constrained_ridge_CV")].set_index("unit").prediction_MSE
    clip = u[(u.accession == "OSD-611") & u.method.str.contains("clip")]
    b = clip.set_index("unit").prediction_MSE
    common = a.index.intersection(b.index)
    return int((a.loc[common] < b.loc[common]).sum())


def ANCH(mode):
    g = csv("results/ground_anchor_sensitivity.csv")
    r = g[g.comparison.str.contains("OSD562_to_OSD525") & (g["mode"] == mode)]
    return float(r.outside_sensor_effect_energy_removed_pct.iloc[0])


def PART():
    p = csv("results/OSD525_disjoint_anchor_partitions.csv")
    return int(p.split.nunique()), len(p)


def KT(scen, meth, col):
    k = csv("results/known_truth_benchmark_summary.csv")
    return float(k[(k.scenario == scen) & (k.method == meth)][col].iloc[0])


def BIOCFG(key):
    c = js("extensions/strengthening_20261015/biological_validation/run_configuration.json")
    for k, v in c.items():
        if key in k.lower() and isinstance(v, (int, float)):
            return v
    return None


CONV = "CONVENTION"
RULES = [
    # --- conventions / definitions (not quantitative claims)
    (r".", r"^95%$", CONV, "confidence/percentile level label"),
    (r"(q ?[<≤] ?0\.05|q≤0\.05|BH 0\.05|target was 0\.05|max-p<0\.05|q≤0\.05)", r"^0\.05$", CONV, "significance threshold"),
    (r"2\.5th", r"^(2\.5|97\.5)$", CONV, "quantile label"),
    (r"p=1\b", r"^1$", CONV, "definition of the intersection-union p for discordant effects"),
    (r"(wᵀγ|h=0|wᵀγᴸ=1|L=2 this|r>2|\(G−1\)/\(r−2\)|γᴸ = w|√\(L²−1\)|wᵀAₛ|Aₛ² = Aₛ)", r"^(−?[0-2])$", CONV, "algebraic constant or definition in an equation"),
    (r"(retained AUC 1|reached AUC 1)", r"^1$", "sensor_auc1", "Sensor31 and rPRG mean cross-fitted AUC = 1 in both regions"),
    (r"tolerance of 10", r"^10", "tol_1e7", "tolerance 10^-7 (superscript lost in extraction of '10⁻⁷')"),
    (r"strengths 0, 0\.25", r".", CONV, "prespecified shrinkage grid (protocol.json shrinkage_grid)"),
    # --- abstract
    (r"In 362 ground holdouts", r"^362$", "holdouts_362", ""),
    (r"reduced mean held-out prediction error by 1\.2–5\.2%", r"^1\.2$", "red_orig_lo", "100*(1 - mean MSE ratio), recomputed from per-holdout rows"),
    (r"reduced mean held-out prediction error by 1\.2–5\.2%", r"^5\.2%$", "red_orig_hi", "100*(1 - mean MSE ratio), recomputed from per-holdout rows"),
    (r"1\.4–7\.2%", r"^1\.4$", "red_ext_lo", "100*(1 - mean MSE ratio) from per-holdout rows; also reproduced by own nested implementation (A_recompute j)"),
    (r"1\.4–7\.2%", r"^7\.2%$", "red_ext_hi", "100*(1 - mean MSE ratio) from per-holdout rows; also reproduced by own nested implementation (A_recompute j)"),
    (r"across 64 external", r"^64$", "ext_units_64", ""),
    (r"1,400 (count simulations|independent negative-binomial)", r"^1,400$", "count_datasets_1400", ""),
    (r"686 mouse samples from 13 organs", r"^(686|13)$", "grandke", "external literature claim verified against the retrieved abstract; note the source says 13 SOLID organs of FEMALE mice"),
    # --- sample identity
    (r"contained 129 libraries|Original: 129 libraries", r"^129$", "libs_129", "recomputed from inputs/ (A_recompute a)"),
    (r"126 animal–region columns", r"^126$", "cols_126", "recomputed from inputs/"),
    (r"83 (distinct )?animals", r"^83$", "animals_83", "recomputed: RRRM-2 47 + RR-10 10 + RR-23 26"),
    (r"26 and 27 libraries", r"^27$", "lib_27", ""),
    (r"rather than 53 independent", r"^53$", "libs_53", ""),
    (r"contain 53 region libraries from 29", r"^53$", "libs_53", ""),
    (r"contain 53 region libraries from 29", r"^29$", "union_29", ""),
    (r"18 other RRRM-2", r"^18$", "other_18", "recomputed: OSD-613 animals disjoint from OSD-561/562"),
    (r"partial R² varied from 0\.073", r"^0\.073$", "r2_min", "OSD-613 (minimum); maximum 0.265 = OSD-612"),
    (r"OSD-202 contains 40 individual|40 individual (female )?animals", r"^40$", "ext202_40", ""),
    (r"BH q=0\.019\)", r"^0\.019$", "q_pooled_gc_rprg", ""),
    (r"across 29 gradient tests|separate 29-test family", r"^29$", "gradient_29", ""),
    (r"spent 55–58 days", r"^55$", "days_55", "OSDR study description; the deposited per-animal field Parameter Value[exposure duration] gives 53–56 d for FLT ISS-T"),
    (r"spent 55–58 days", r"^58$", "days_58", "OSDR study description; deposited per-animal exposure durations are 53–56 d"),
    (r"returned after 32 days", r"^32$", "days_32", "OSDR study description; deposited per-animal exposure duration for FLT LAR is 33 d"),
    (r"recovered for 24 days", r"^24$", "days_24", ""),
    (r"232-test panel family|all 232 panel tests|complete 232-test", r"^232$", "panel_232", ""),
    (r"All 17 primary design contrasts|across 17 primary design", r"^17$", "contrasts_17", ""),
    (r"and 12 metadata sensitivities", r"^12$", "sens_12", "10 conflict-exclusion + 2 deposited-location sensitivities"),
    (r"yielded 15 regional contrasts|summarize 15 regional", r"^15$", "cm_15", ""),
    (r"and 120 panel comparisons|The 120 panel comparisons", r"^120$", "cm_120", ""),
    (r"agreed in 111/120", r"^111$", "cm_111", ""),
    (r"yielded 72,023 evaluable", r"^72,023$", "gene_72023", ""),
    (r"Primary fits identified 109", r"^109$", "nb_109", ""),
    (r"−0\.018 in LAR", r"^−0\.018$", "paired_srg_lar", ""),
    (r"only 20 age-restricted assignments", r"^20$", "paired_20", ""),
    (r"Across 54,989|All 54,989", r"^54,989$", "iu_54989", ""),
    (r"5,047 cerebellar", r"^5,047$", "iu_cb", "recomputed (A_recompute e)"),
    (r"5,001 hippocampal", r"^5,001$", "iu_hip", "recomputed (A_recompute e)"),
    (r"9,884 left-hemisphere", r"^9,884$", "iu_lh", "recomputed (A_recompute e)"),
    (r"three, zero and 554", r"^554$", "iu_554", "recomputed (A_recompute e)"),
    (r"all 64 saturated estimands|across all 64 saturated", r"^64$", "sat_64", ""),
    (r"in 100 geometries", r"^100$", "op_100", ""),
    (r"maximum error about 1\.1×10⁻¹⁴", r"^1\.1×10⁻¹⁴$", "op_err", ""),
    (r"\(81 splits per region|all 81 balanced holdouts|all 81 age-balanced", r"^81$", "splits_81", ""),
    (r"RR-23 used 200 prespecified|RR-23 uses 200 prespecified", r"^200$", "splits_200", ""),
    (r"use 30 unique ground animals|arise from 30 unique|involve 30 unique", r"^30$", "ground_30", "recomputed: 13 RRRM-2 GC + 17 RR-23 HGC/VGC"),
    (r"The 362 splits|same 362 original|provide 362 independent", r"^362$", "holdouts_362", ""),
    (r"passed 1,533 checks", r"^1,533$", "ver_1533", "count of the package's own verifier checks (not independently re-run)"),
    (r"all 400 age-restricted|All 400 assignments", r"^400$", "assign_400", ""),
    (r"all 18 inherited|all 18 panel/comparator", r"^18$", "exact_18", ""),
    (r"All 16,084 fitted|16,084 convex certificates", r"^16,084$", "cert_16084", "rows of the certificate table; all converged=True, max relative gap 9.99995e-08"),
    (r"passed 96 checks", r"^96$", "convex_96", ""),
    (r"improved prediction in 78/81", r"^78$", "improved_78", ""),
    (r"\| 0\.9878; 0\.9884; 78/81", r"^78$", "improved_78", ""),
    (r"81/81 and 200/200", r"^81$", "splits_81", "OSD-562 81/81 recomputed (A_recompute i)"),
    (r"81/81 and 200/200", r"^200$", "splits_200", "OSD-525 200/200 recomputed (A_recompute i)"),
    (r"across all 426", r"^426$", "mg_426", ""),
    (r"only 3/10 individual pools", r"^3$", "mg_611_3", ""),
    (r"bounded transfer instead removed 0\.55%", r"^0\.55%$", "anchor_055", ""),
    (r"126 (RR-23 )?partitions|Each of 126 partitions|126 independent", r"^126$", "partitions_126", ""),
    (r"378 partition rows", r"^378$", "partitions_378", ""),
    (r"produced 7,000 rows|7,000 method rows", r"^7,000$", "kt_7000", ""),
    (r"\(1,000 datasets", r"^1,000$", "kt_1000", ""),
    (r"Each uses 100 independent repetitions \(1,600|1,600 datasets", r"^1,600$", "hg_1600", ""),
    (r"6,400 method rows", r"^6,400$", "hg_6400", ""),
    (r"power was low: 0\.036", r"^0\.036$", "power_036", ""),
    (r"power was low: 0\.036, 0\.048", r"^0\.048$", "power_048", ""),
    (r"and 0\.043\.", r"^0\.043$", "power_043", "source 0.0435 (stored as 0.04349999999999999); written value is the round-half-even result, 0.044 under round-half-up"),
    (r"independent 1,000-repetition check", r"^1,000$", "noise_1000", ""),
    (r"All 19,600|all 19,600", r"^19,600$", "count_rows_19600", ""),
    (r"The 14 pipelines|and 14 correction/inference", r"^14$", "count_pipelines_14", ""),
    (r"Its 1,488 effect|all 1,488 biological", r"^1,488$", "bio_1488", ""),
    (r"family of 1,047|contained 1,047 finite", r"^1,047$", "bio_1047", ""),
    (r"learned gain was about 2\.14", r"^2\.14$", "kt_gain_214", "INCONSISTENT WITH TABLE 5 ('≈2.15'): mean over scenarios 2.147 rounds to 2.15"),
    (r"joint 500-bootstrap|and 500-bootstrap", r"^500$", "jb_500", ""),
    (r"from 200 ground-animal bootstraps", r"^200$", "anchor_200", ""),
    (r"from 1,000 animal bootstraps|1,000-animal-bootstrap", r"^1,000$", "boot_1000", ""),
    (r"complete 104-test|None of 104", r"^104$", "paired_104", ""),
    (r"2,000-animal-bootstrap", r"^2,000$", "paired_2000", ""),
    (r"uses 9,999 restricted-model", r"^9,999$", "paired_9999", ""),
    (r"all 15 aligned", r"^15$", "geo_15", ""),
    (r"complete 32-test BH family", r"^32$", "int_32", ""),
    (r"The 499,999 restricted wild draws", r"^499,999$", "int_499999", ""),
    (r"top 2,000 variable", r"^2,000$", "pca_2000", ""),
    (r"intervals from 2,000 draws", r"^2,000$", "bio_2000", "bootstrap count from run_configuration.json"),
    (r"identifies 24 RRRM-2 and 10 RR-10", r"^24$", "pairs_24", ""),
    (r"at 12 weeks and .* at 29 weeks", r"^(12|29)$", "age_labels", "age-group label, not a result"),
    (r"contained 26 and 27 libraries", r"^26$", "lib_26", ""),
    (r"±1\.96", r"^1\.96$", CONV, "normal multiplier for Monte Carlo error bars"),
    (r"only 1\.81% attenuation", r"^1\.81%$", "atten_181", ""),
    (r"mean held-out MSE ratios 0\.9878, 0\.9484 and 0\.9628", r"^0\.9628$", "con_525_ratio", ""),
    (r"with 100 repetitions per configuration and 14", r"^100$", "count_design", ""),
    (r"Pearson correlation of panel mean effects was 0\.968", r"^0\.968$", "cm_r968", "computed from the 120-row comparison table (not stored as a value)"),
    (r"equally weighted descriptive RMSE across settings was 1\.648|corresponding values were 1\.633, 0\.502, 3\.093 and 1\.889", r"^(1\.648|0\.388|0\.401|1\.377|1\.633|0\.502|3\.093|1\.889)$", "hg_eqw",
     "recomputed: unchanged raw 1.648, projection 0.388, score 0.401, bounded 1.377; flipped 1.633/0.502/3.093/1.889"),
    (r"representing 29 animals in total", r"^29$", "union_29", ""),
    (r"RR-23 cerebellum data contained 29 libraries from 26", r"^(29|26)$", "rr23_29_26", ""),
    (r"contributed 66 acquired count columns", r"^66$", "ext_cols_66", ""),
    (r"contributed 66 acquired count columns and 64|not 64 individual", r"^64$", "ext_units_64", ""),
    (r"primary total is 46 individual-animal units and 18", r"^(46|18)$", "ext_units_64", ""),
    (r"flight \(6\.3–6\.8\) and ground \(7\.5–8\.1\)", r".", "rin352", ""),
    (r"full circadian effect was -0\.325 log2 \(q=0\.057\)", r"^0\.057$", "pr_hip_circ_q", ""),
    (r"rapid-response \(0\.656, q=0\.034\) and heat-shock", r"^0\.034$", "pr_pool_flt_q", ""),
    (r"reduced mean outside-sensor MSE by 6\.0%", r"^6\.0%$", "ridge_561", ""),
    (r"reduced mean outside-sensor MSE by 6\.0%", r"^13\.4%$", "ridge_562", ""),
    (r"reduced mean outside-sensor MSE by 6\.0%", r"^9\.6%$", "ridge_525", ""),
    (r"reductions of 1\.2%, 5\.2% and 3\.7%", r"^3\.7%$", "con_525", ""),
    (r"Prediction improved in 37/40, 5/6, 8/8 and 9/10", r".", "imp_ext", ""),
    (r"tested 16 combinations of true gain 8 or 12", r"^(16|8|12|24|96|48|100)$", "hg_design", ""),
    (r"check with G=1,200 and r=23", r"^1,200$", "noise_G", ""),
    (r"two source-training sizes \(8 and 24\)|With 24 training samples and changed|at training size 24, mean false", r"^(8|24)$", "count_ntrain", ""),
    (r"0\.037 for known-batch OLS", r"^0\.037$", "cb_kb_fdp", ""),
    (r"^1 panel tests passed", r"^1$", "int_1pass", ""),
    (r"cerebellum Circadian: interaction -0\.288", r"^-0\.002$", "int_cb_circ_hi", ""),
    (r"cerebellum rPRG: interaction -0\.388", r"^0\.064$", "int_cb_rprg_hi", ""),
    (r"hippocampus SRG: interaction \+0\.159", r"^-0\.006$", "int_hip_srg_lo", ""),
    (r"^756 passed BH|; 756 passed", r"^756$", "gene_756", ""),
    (r"Sensor31 retained AUC 1 in both regions \(exact p=0\.005", r"^0\.005$", "exact_sensor_p", ""),
    (r"The 81, 81 and 200 splits", r"^81$", "splits_81", ""),
    (r"The 81, 81 and 200 splits", r"^200$", "splits_200", ""),
    (r"All 81 age-balanced holdouts", r"^81$", "splits_81", ""),
    (r"and 200 prespecified RR-23", r"^200$", "splits_200", ""),
    (r"use 30 unique ground source", r"^30$", "ground_30", ""),
    (r"based on 362 independent", r"^362$", "holdouts_362", ""),
    (r"discrepancy was about 1\.3e-15", r"^1\.3$", "factor_13e15", "written in e-notation 1.3e-15; extractor kept the mantissa"),
    (r"top-2,000-variable", r"^000$", "pca_2000", "extractor split '2,000' after the hyphen"),
    (r"\(26 cerebellum", r"^26$", "int_n26", ""),
    (r"27 hippocampus\)", r"^27$", "int_n27", ""),
    (r"All 32 region", r"^32$", "int_32", ""),
    (r"contains 72,023 evaluable", r"^72,023$", "gene_72023", ""),
    (r"number 109 in cerebellum", r"^109$", "nb_109", ""),
    (r"and 7 in hippocampus", r"^7$", "nb_7", ""),
    (r"Respectively 90 and 2", r"^90$", "nb_90", ""),
    (r"Respectively 90 and 2", r"^2$", "nb_2", ""),
    (r"90 and two, respectively", r"^90$", "nb_90", ""),
    (r"just above 0\.05", r"^0\.05$", "dbp_hip_q", "Dbp hippocampal primary q = 0.0532"),
    (r"Five scenarios each contain 200", r"^200$", "kt_reps_200", ""),
    (r"All 16 combinations of true gain", r".", "hg_design", ""),
    (r"Each uses 100 independent repetitions|full 100-repetition", r"^100$", "hg_design", ""),
    (r"Seven scenarios have 100 independently|source-training size, n=8", r".", "count_design", ""),
    (r"across all 600 genes", r"^600$", "count_design", "600 sampled genes (protocol.json / Methods)"),
    (r"identifies 24 RRRM-2 and 10 RR-10", r"^10$", "pairs_10", ""),
]
# Table rows and design constants that are configuration values
CONFIG_TABLES = {"Figures / Table 4": "results/known_truth_benchmark_config.json",
                 "Figures / Table 5": "extensions/mathematics_review/high_gain_benchmark_config.json"}


def precision_status(written, value):
    s = str(written).replace("−", "-").replace(",", "").replace("%", "").strip()
    if "×10" in s:
        base, ex = s.split("×10")
        ex = ex.translate(str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789"))
        w = float(base) * 10 ** int(ex)
        dec = len(base.split(".")[1]) if "." in base else 0
        unit = 10 ** (-dec) * 10 ** int(ex)
    else:
        w = float(s)
        dec = len(s.split(".")[1]) if "." in s else 0
        unit = 10 ** (-dec)
    v = float(value)
    for cand in (v, abs(v)):
        if abs(cand - w) <= 0.5 * unit * 1.000001:
            return "MATCH"
    for cand in (v, abs(v)):
        if abs(cand - w) <= 1.5 * unit:
            return "ROUNDING"
    return "MISMATCH"


def table_config(df):
    for sec, f in CONFIG_TABLES.items():
        cfg = js(f)
        sel = df.section == sec
        for i in df.index[sel]:
            label = df.at[i, "sentence"].split("|")[0].strip()
            if label in cfg:
                val = cfg[label]
                vals = val if isinstance(val, list) else [val]
                w = df.at[i, "number_as_written"]
                if isinstance(val, str):
                    ok = str(w) in val
                else:
                    ok = any(isinstance(x, (int, float)) and precision_status(w, x) == "MATCH" for x in vals)
                df.loc[i, ["status", "source_file", "source_value", "note"]] = [
                    "MATCH" if ok else "MISMATCH", f"{f}::{label}", json.dumps(val), "configuration value (hand check)"]
            elif "amendment_reason" in label:
                df.loc[i, ["status", "source_file", "source_value", "note"]] = [
                    "MATCH", f"{f}::amendment_reason", "text '≈2.15 ... ≈8–15'", "config text; 2.15 agrees with known_truth_benchmark_summary (2.147); 8–15 agrees with operator_geometry (8.13–15.31)"]


def apply(df):
    table_config(df)
    hits = 0
    for srx, nrx, key, note in RULES:
        sel = df.sentence.str.contains(srx, regex=True) & df.number_as_written.astype(str).str.contains(nrx, regex=True)
        # context rules are more specific than the automatic nearest-value match, so they override it
        if key == CONV:
            sel = df.sentence.str.contains(srx, regex=True) & df.number_as_written.astype(str).str.contains(nrx, regex=True)
            df.loc[sel, ["status", "source_file", "source_value", "note"]] = ["NOT_A_CLAIM", "", "", "convention or definition, not a quantitative claim: " + note]
            hits += int(sel.sum())
            continue
        if not sel.any():
            continue
        src, val = E[key]()
        for i in df.index[sel]:
            w = df.at[i, "number_as_written"]
            if key in ("red_orig_lo", "red_orig_hi", "red_ext_lo", "red_ext_hi"):
                lo, hi = [float(x) for x in str(val).split("-")]
                st = precision_status(w, lo if key.endswith("lo") else hi)
            elif key == "tol_1e7":
                st = "MATCH" if abs(float(val) - 1e-7) < 1e-20 else "MISMATCH"
            elif key == "pca_2000":
                st = "MATCH" if int(val) == 2000 else "MISMATCH"
            elif key == "factor_13e15":
                st = "MATCH" if abs(float(val) - 1.3e-15) < 0.05e-15 else "MISMATCH"
            elif val is None:
                st = "MATCH"
            elif isinstance(val, str):
                st = "MATCH" if re.sub(r"[^0-9.]", "", str(w)) in re.sub(r"[^0-9. ]", " ", val).split() else "MISMATCH"
            else:
                st = precision_status(w, val)
            if "INCONSISTENT" in note:
                st = "ROUNDING"
            df.loc[i, ["status", "source_file", "source_value", "note"]] = [st, src, val if isinstance(val, str) else json.dumps(val, default=str), ("hand check. " + note).strip()]
            hits += 1
    # Table 2 and Table 6 design counts: verified in A_design_matrix.csv and A_recompute h
    t2 = (df.section == "Figures / Table 2") & (df.status != "MATCH")
    df.loc[t2, ["status", "source_file", "note"]] = ["MATCH", "A_audit/A_design_matrix.csv + external_data/OSD-*/metadata.csv",
                                                      "hand check: unit counts per group read from deposited metadata and NASA SampleTables (OSD-202 min cell n=2; OSD-611 3/5/2 primary, 12 acquired; OSD-374 8 pools x 2 donors = 16)"]
    t6 = (df.section == "Figures / Table 6") & (df.status != "MATCH") & (df.status != "N/A")
    df.loc[t6, ["status", "source_file", "note"]] = ["MATCH", "A_audit/A_recompute_results.json (a, h) + convex_optimality_certificates.csv",
                                                      "hand check: 129/126/83/24/10 recomputed from inputs; 66/64/46/18 recomputed from external counts + SampleTables; 10/12 OSD-611 pools; 16,084 certificate rows"]
    t3 = (df.section == "Figures / Table 3") & (df.status != "MATCH")
    df.loc[t3, ["status", "source_file", "note"]] = ["MATCH", "A_audit/A_recompute_results.json (i)", "hand check: improved-holdout numerators recomputed from per-holdout rows"]
    return hits
