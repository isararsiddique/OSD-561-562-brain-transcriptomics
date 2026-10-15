#!/usr/bin/env python3
"""Export whole electronic result tables under numbered Supplementary Data names.

Copies preserve every source byte, including missing values and failed tests.
The full original result trees accompany these convenient journal-upload names.
"""
import hashlib
import json
import shutil
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'supplementary_data'
TABLES = [
 ('results/sample_manifest.csv', 'sample_provenance', '129 deposited libraries; animal and design labels'),
 ('results/animal_region_manifest.csv', 'animal_region_provenance', '126 collapsed animal-region observations'),
 ('results/metadata_conflicts.csv', 'metadata_conflicts', 'All recorded deposited location conflicts'),
 ('results/panel_results.csv', 'primary_panel_tests', 'Complete 232-test primary and sensitivity panel family'),
 ('results/gradient_results.csv', 'gradient_tests', 'Complete 29-test gradient family'),
 ('results/count_model_gene_results.csv.gz', 'regional_count_gene_tests', 'All regional negative-binomial gene results'),
 ('results/count_model_panel_comparison.csv', 'count_panel_agreement', 'All 120 panel estimator comparisons'),
 ('extensions/paired_biology/paired_panel_interactions.csv', 'paired_panel_tests', 'Complete 104-test paired panel family'),
 ('extensions/paired_biology/paired_gene_interactions_exploratory.csv.gz', 'paired_gene_tests', 'Complete exported exploratory paired gene tests'),
 ('results/regional_gene_replication.csv.gz', 'cross_cohort_gene_screen', 'Complete 54989-entry intersection-union family'),
 ('extensions/ground_holdout/holdout_split_results.csv', 'ground_holdout_results', 'All 362 overlapping splits and five methods'),
 ('results/ground_anchor_sensitivity.csv', 'ground_anchor_transfer', 'All 15 anchor applications and four methods'),
 ('results/operator_geometry.csv', 'operator_geometry', 'All 15 aligned real-data correction geometries'),
 ('results/joint_bootstrap_summary.csv', 'joint_refit_uncertainty', 'Complete joint refit summary'),
 ('results/joint_bootstrap_repetitions.csv', 'joint_refit_repetitions', 'All retained joint bootstrap repetitions'),
 ('results/known_truth_benchmark_repetitions.csv', 'primary_simulation_repetitions', 'All 1000 independent primary datasets and seven methods'),
 ('extensions/mathematics_review/high_gain_benchmark_repetitions.csv', 'high_gain_simulation_repetitions', 'All 1600 independent stress datasets and four methods'),
 ('extensions/flight_collection_interaction/panel_flight_collection_interactions.csv', 'flight_collection_panel_tests', 'Complete 32-test exploratory interaction family'),
 ('extensions/flight_collection_interaction/gene_flight_collection_interactions_exploratory.csv.gz', 'flight_collection_gene_tests', 'All gene interactions including unavailable p-values; 72023 evaluable tests'),
 ('extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv', 'primary_gene_candidate_sensitivity', 'Every primary FDR candidate including sensitivity failures'),
 ('extensions/independent_audit_20261015/exact_ground_sensor_association.csv', 'exact_ground_complete_family', 'All 18 inherited panel and comparator by region tests'),
 ('extensions/independent_audit_20261015/OSD-561_exact_label_null.csv', 'cerebellum_exact_label_null', 'All 400 assignments and all nine comparators'),
 ('extensions/independent_audit_20261015/OSD-562_exact_label_null.csv', 'hippocampus_exact_label_null', 'All 400 assignments and all nine comparators'),
 ('extensions/independent_audit_20261015/saturated_age_panel_interactions.csv', 'saturated_age_complete_family', 'All 64 saturated-age panel estimands'),
 ('extensions/independent_audit_20261015/heldout_animal_prediction_errors.csv', 'heldout_animal_prediction_errors', 'Every retained animal and overlapping split prediction error'),
 ('extensions/independent_audit_20261015/training_expression_tertile_prediction_errors.csv', 'expression_tertile_prediction_errors', 'Every training-expression tertile prediction diagnostic'),
 ('extensions/real_data_structure/animal_count_qc.csv', 'all_deposit_count_diagnostics', 'All 126 animal-region observations in seven deposits'),
 ('extensions/real_data_structure/animal_pca_scores.csv', 'all_deposit_pca_scores', 'Every animal-region PCA score'),
 ('results/contrast_manifest.csv', 'complete_contrast_manifest', 'Full primary and sensitivity contrast definitions'),
('extensions/strengthening_20261015/external_data/cohort_inventory.csv', 'external_inventory', 'All acquired count matrices and provenance'),
('extensions/strengthening_20261015/external_data/harmonized_metadata.csv', 'external_sample_metadata', 'All 66 acquired units with primary mask'),
('extensions/strengthening_20261015/external_data/sample_to_source_membership.csv', 'external_pool_membership', 'All recorded donor and library memberships'),
('extensions/strengthening_20261015/method_benchmark/outer_split_results.csv', 'expanded_original_methods', 'All 3258 initial method by split rows'),
('extensions/strengthening_20261015/method_benchmark/outer_animal_prediction_errors.csv', 'expanded_original_unit_errors', 'All individual animal errors across overlapping splits'),
('extensions/strengthening_20261015/method_benchmark/strengthened_outer_results.csv', 'constrained_external_predictions', 'All 1002 new original and external method rows'),
('extensions/strengthening_20261015/method_benchmark/strengthened_unit_errors.csv', 'constrained_external_unit_errors', 'All 2488 individual retained errors'),
('extensions/strengthening_20261015/method_benchmark/convex_optimality_certificates.csv', 'convex_certificates', 'All 16084 fitted convex optimization certificates'),
('extensions/strengthening_20261015/method_benchmark/strengthened_inner_CV.csv', 'constrained_nested_selection', 'All nested tuning candidates'),
('extensions/strengthening_20261015/method_benchmark/strengthened_split_roles.csv', 'new_split_membership', 'All original and external training/test roles'),
('extensions/strengthening_20261015/method_benchmark/matched_grid_outer_results.csv', 'matched_grid_predictions', 'All 852 matchedgrid comparator rows'),
('extensions/strengthening_20261015/method_benchmark/matched_grid_fixed_alpha_objectives.csv', 'fixed_penalty_objective_checks', 'All 28 representative full-space objective comparisons'),
('extensions/strengthening_20261015/correction_benchmark/all_repetitions.csv', 'count_benchmark_repetitions', 'All 19600 method runs over 1400 independent datasets'),
('extensions/strengthening_20261015/correction_benchmark/real_count_calibration.csv', 'count_benchmark_calibration', 'All 600 empirical count calibration genes'),
('extensions/strengthening_20261015/biological_validation/all_biological_effects.csv', 'expanded_biological_family', 'All 1488 effect and screen rows;1047 finite probabilities'),
('extensions/strengthening_20261015/biological_validation/all_animal_omissions.csv', 'expanded_animal_omissions', 'All animal omission effects including failures'),
('extensions/strengthening_20261015/biological_validation/all_gene_omissions.csv', 'expanded_gene_omissions', 'All gene omission effects including failures'),
('extensions/strengthening_20261015/biological_validation/design_cell_support.csv', 'expanded_design_support', 'All recorded design cells'),
('extensions/strengthening_20261015/biological_validation/endpoint_availability.csv', 'expanded_endpoint_availability', 'All panel and marker gene availability'),
('extensions/strengthening_20261015/biological_validation/harmonized_panel_direction_screen.csv', 'expanded_direction_screen', 'All common-gene direction screens'),
('extensions/strengthening_20261015/biological_validation/all_cross_context_effect_differences.csv', 'expanded_context_differences', 'All cross-context differences'),

]

def main():
    OUT.mkdir(exist_ok=True)
    records = []
    for number, (source, label, scope) in enumerate(TABLES, 1):
        path = ROOT / source
        suffix = '.csv.gz' if source.endswith('.csv.gz') else '.csv'
        name = f'Supplementary_Data_{number:02d}_{label}{suffix}'
        target = OUT / name
        shutil.copy2(path, target)
        original = hashlib.sha256(path.read_bytes()).hexdigest()
        assert hashlib.sha256(target.read_bytes()).hexdigest() == original
        data = pd.read_csv(target)
        records.append({'number':number, 'file':name, 'original_analysis_path':source,
                        'scope':scope, 'data_rows':len(data), 'columns':len(data.columns),
                        'sha256':original})
    pd.DataFrame(records).to_csv(OUT/'INDEX.csv', index=False)
    (OUT/'README.md').write_text(
        '# Electronic Supplementary Data\n\n'
        f'INDEX.csv maps {len(records)} numbered complete-table exports to the original analysis paths, '
        'row/column dimensions and SHA-256 hashes. Compressed CSV files contain large gene tables; '
        'decompress them for a journal portal that requires plain CSV. These are byte-identical copies, '
        'with no discovery-only filtering or revised multiplicity families. The full results/ and '
        'extensions/ trees also retain constituent genes, normalization, design support, omitted-case '
        'diagnostics, software protocols and additional intermediate tables.\n\n'
        'Overlapping holdouts are not independent animal samples. Simulated repetitions are generated '
        'datasets, not real animals. Exact-label families, primary panel families, paired families and '
        'exploratory interaction families are separately declared in Methods. Missing probabilities '
        'remain unavailable. This folder is an upload convenience for the electronic data referenced '
        'in the Article and Supplementary Table S1; reproduction uses the original paths.\n')
    (ROOT/'results/supplementary_data_validation.json').write_text(json.dumps(
        {'complete_table_exports':len(records), 'source_copy_hashes_match':True,
         'all_source_rows_retained':True, 'records':records}, indent=2)+'\n')
    print(f'Exported {len(records)} complete numbered Supplementary Data tables.')

if __name__ == '__main__':
    main()
