"""The complete study in dependency order. Each stage is a separate process."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Step:
    name: str
    script: str
    args: tuple = ()


def step(name, group, script, *args):
    return Step(name, f"stages/{group}/{script}", args)


STEPS = [
    step("core_inputs", "core", "prepare_inputs.py"),
    step("orthology", "core", "refine_mapping.py"),
    Step("check_panel", "check_panel.py"),
    step("primary_source_identity", "core", "verify_deposited_inputs.py"),
    step("bulk_deseq2", "core", "fit_bulk.py"),
    step("bulk_rsem_sensitivity", "core", "fit_bulk.py", "--mode", "rounded_RSEM"),
    step("geomx", "core", "fit_spatial.py"),
    step("core_sensitivities", "core", "secondary_checks.py"),
    step("category_exact", "core", "category_tests.py"),
    step("gene_exact", "core", "gene_exact_sensitivity.py"),
    step("nucleus_context", "core", "single_nucleus_context.py"),
    step("independent_inputs", "extended", "prepare_validation.py"),
    step("multiome", "extended", "aggregate_multiome.py"),
    step("bulk_edger_voom", "extended", "compare_expression_methods.R"),
    step("reference_mixtures", "extended", "composition_models.py"),
    step("composition_sensitivity", "extended", "adjustment_sensitivity.py"),
    step("promoter", "extended", "promoter_context.py"),
    step("pathway_camera", "extended", "pathway_camera.R"),
    step("pathway_exact", "extended", "pathway_animal_tests.py"),
    step("rr9_ncounter", "extended", "analyze_RR9_panel.py"),
    step("protein_source_audit", "extended", "audit_protein_sources.py"),
    step("visium_inputs", "extended", "extract_visium_inputs.py"),
    step("visium", "extended", "analyze_visium.py"),
    step("sorted_identity", "extended", "audit_additional_cell_labels.py"),
    step("sorted_qc", "extended", "prepare_sorted_cells.py"),
    step("sorted_accessions", "extended", "audit_RR8_NASA_accessions.py"),
    step("sorted_marker_qc", "extended", "sorted_marker_sensitivity.py"),
    step("sorted_models", "extended", "analyze_sorted_cells.R"),
    step("sorted_exact_pairing", "extended", "sorted_exact_and_pairing.py"),
    step("sorted_camera", "extended", "sorted_pathway_camera.R"),
    step("sorted_pathway_exact", "extended", "sorted_pathway_exact.py"),
    step("sorted_mouse_deletions", "extended", "sorted_leave_one_mouse.R"),
    step("sorted_depth", "extended", "sorted_cell_depth_endpoints.py"),
    step("sorted_relative", "extended", "infer_additional_validation.py"),
    step("camera_assumptions", "extended", "audit_camera_assumptions.R"),
    step("continuous_evidence", "extended", "analyze_evidence_distributions.py"),
    step("additional_inputs", "deep", "prepare_additional_cohorts.py"),
    step("additional_models", "deep", "analyze_additional_cohorts.R"),
    step("additional_mouse_deletions", "deep", "validate_additional_mouse_deletions.R"),
    step("regional_inputs", "deep", "prepare_terminal_regions.py"),
    step("regional_models", "deep", "analyze_terminal_regions.R"),
    step("expanded_families", "deep", "expand_validation_family.py"),
    step("expression_states", "deep", "prepare_expression_states.py"),
    step("state_models", "deep", "analyze_expression_states.R"),
    step("state_composition", "deep", "analyze_mouse_state_composition.py"),
    step("technical_states", "deep", "prepare_expression_states.py", "--depth-adjusted"),
    step("technical_state_models", "deep", "analyze_expression_states.R", "{deep}/depth_sensitivity"),
    step("technical_composition", "deep", "analyze_mouse_state_composition.py", "--depth-adjusted"),
    step("article_qc", "deep", "prepare_article_QC.py"),
    step("article_models", "deep", "analyze_article_QC.R"),
    Step("result_tables", "build_result_tables.py"),
    step("core_figures", "core", "make_figures.py"),
    step("integrated_figures", "extended", "make_figures.py"),
    step("immune_figure", "extended", "make_mechanistic_figure.py"),
    step("mouse_validation_figure", "extended", "make_validation_figure.py"),
    step("deep_figures", "deep", "make_deep_figures.py"),
    step("verify_core", "core", "verify_results.py"),
    step("verify_extension", "extended", "verify_extension.py"),
    step("verify_figure_values", "extended", "verify_figure_values.py"),
    step("verify_deep", "deep", "verify_deep_analysis.py"),
    Step("report", "build_summary.py"),
    Step("figure_files", "check_figures.py"),
]
FIGURES = [s for s in STEPS if s.name in {"result_tables", "core_figures", "integrated_figures",
           "immune_figure", "mouse_validation_figure", "deep_figures", "verify_figure_values", "figure_files"}]
VERIFICATION = [s for s in STEPS if s.name.startswith("verify_") or s.name == "figure_files"]
