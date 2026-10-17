# Figure captions and source tables

## fig01_design_and_QC

Counts are sequenced animals within region. PCA uses centered log2 normalized endogenous expression of the 2,000 highest-variance genes; age 29 uses larger markers. Terminal flight preservation is confounded with exposure.

Sources: bulk_sample_registry.csv;verified_panel_mapping.csv;bulk_panel_complete.csv;spatial_panel_complete.csv

## fig02_complete_bulk_panel

All 151 corrected source analytes are retained. Rows use current human symbols; original labels and orthology ambiguity remain in the tables. Color saturation is display clipping, not changed estimates.

Sources: verified_panel_mapping.csv;bulk_panel_complete.csv

## fig03_AIF1_PYCARD

Individual animal points and full-model Wald intervals. These intervals are nominal, not simultaneous. Results do not establish age or region specificity. Rounded RSEM and leave-one-out are shared-sample sensitivities; exact-gene sensitivity has no full-family q<0.05.

Sources: bulk_panel_complete.csv;bulk_count_source_sensitivity.csv;bulk_flight_leave_one_out.csv;bulk_gene_exact_sensitivity.csv

## fig04_collection_context_and_pairing

Fos and Junb are measured in the same hippocampal animals. Paired analysis is within RRRM-2, not independent replication; per-gene paired q-values and the 127-gene omnibus are nonsignificant.

Sources: bulk_panel_complete.csv;paired_region_flight_associations.csv;paired_region_omnibus.json;paired_region_omnibus_null.csv

## fig05_spatial_factorial_uncertainty

The full interaction is (FLT_BuOE−GC_BuOE)−(FLT_saline−GC_saline). Five fixed-df HC3 calls do not survive the Welch family correction. No model separates slide from group.

Sources: spatial_panel_complete.csv;spatial_sample_registry.csv;spatial_welch_sensitivity.csv

## figS01_complete_spatial_panel

All 151 corrected source analytes are retained. Rows use current human symbols; original labels and orthology ambiguity remain in the tables. Color saturation is display clipping, not changed estimates.

Sources: verified_panel_mapping.csv;spatial_panel_complete.csv

## fig06_RR3_celltype_context

21,178 author-annotated nuclei aggregate to five animals. Neurons, CNS glia and immune groups pass the fixed cell-count gate; vascular and neural crest-like groups do not. AIF1/PYCARD effects do not replicate the negative cerebellum direction.

Sources: OSD-352_panel_celltype_context.csv;OSD-352_pseudobulk_registry.csv;OSD-352_panel_celltype_exact_tests.csv

## figS02_five_views

Scatter, hexbin, smoothed density contours, marginal distributions and 2D binned density use exactly the same gene estimates. Density is descriptive; shared-sample agreement is not independent validation.

Sources: five_view_source_data.csv;bulk_count_source_sensitivity.csv
