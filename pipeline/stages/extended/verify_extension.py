from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Independent numerical and identity checks; no implementation-mirroring unit tests."""
from pathlib import Path
import hashlib,json,re
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
from statsmodels.stats.multitest import multipletests
B=stage_root('extended');V1=stage_root('core');checks=[]
def check(name,condition):
    checks.append({'check':name,'passes':bool(condition)})
    assert condition,name
for name,expected in [('bulk_method_comparison_panel',4228),('independent_bulk_panel',453),('celltype_pseudobulk_panel',755)]:
    a=pd.read_csv(B/f'tables/{name}.csv')
    check(name+'_family_size',a.groupby('method').size().eq(expected).all())
    for method,s in a.groupby('method'):
        q=multipletests(s.pvalue.fillna(1).to_numpy(),method='fdr_bh')[1];check(name+'_'+method+'_independent_BH',np.allclose(q,s.q_full_family,rtol=1e-10,atol=1e-12))
        check(name+'_'+method+'_unique_keys',not s[['study','contrast','analyte']].duplicated().any())
map=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('')
check('151_corrected_targets_and_143_primary_symbols',len(map)==151 and len(set(map.mouse_symbol_primary)-{''})==143)
ann=pd.read_csv(B/'prepared/OSD352_RNA_annotation.csv').set_index('ENSEMBL')
rna=pd.read_csv(B/'prepared/OSD352_RNA_pseudobulk.csv.gz',index_col=0);rna.index=ann.loc[rna.index,'SYMBOL'];rna=rna.groupby(level=0).sum()
old=pd.read_csv(V1/'prepared/OSD-352_panel_pseudobulk_counts.csv',index_col=0)
check('full_RNA_aggregation_matches_frozen_143_gene_25_donor_type_counts',np.array_equal(old.to_numpy(),rna.loc[old.index,old.columns].to_numpy()))
reg=pd.read_csv(B/'tables/OSD352_l2_registry.csv');prev=pd.read_csv(V1/'tables/OSD-352_pseudobulk_registry.csv').set_index('column')
check('full_RNA_library_totals_match_frozen_sparse_aggregation',np.array_equal(reg.set_index('column').RNA_total_UMIs,prev.loc[reg.column].total_RNA_UMIs))
check('21178_nuclei_five_donors_three_flight_two_ground',reg.n_nuclei.sum()==21178 and reg.bio_origin.nunique()==5 and reg.drop_duplicates('bio_origin').flight_status.value_counts().to_dict()=={'flight':3,'ground':2})
check('all_RNA_feature_IDs_unique',ann.index.is_unique)
cv=pd.read_csv(B/'tables/deconvolution_held_donor_validation.csv')
check('five_held_donors_and_500_synthetic_mixtures_each_method',cv.groupby('method').apply(lambda x:len(x[['held_out_donor','simulation']].drop_duplicates()),include_groups=False).eq(500).all())
check('no_held_donor_in_training_reference',all(r.held_out_donor not in r.training_ground_donors.split(';') for r in cv.itertuples()))
signature=pd.read_csv(B/'tables/deconvolution_signature_genes.csv').signature_gene
check('400_signature_genes_exclude_entire_CNS_panel',len(signature)==400 and not set(signature)&set(map.mouse_symbol_primary))
fractions=pd.read_csv(B/'tables/bulk_RNA_fraction_estimates.csv')
check('RNA_fraction_simplexes',np.allclose(fractions.groupby(['study','sample','method']).estimated_RNA_fraction.sum(),1))
check('same_reference_all_five_bulk_studies',fractions.study.nunique()==5)
pro=pd.read_csv(B/'tables/RNA_ATAC_promoter_context.csv')
check('453_promoter_tests_full_family',len(pro)==151*3)
check('promoter_exact_p_grid',np.allclose(pro.ATAC_exact_p*10,np.round(pro.ATAC_exact_p*10)))
check('promoter_BH_independent',np.allclose(pro.q_panel_three_celltypes,multipletests(pro.ATAC_exact_p,method='fdr_bh')[1]))
sec=pd.read_csv(B/'tables/Visium_section_registry.csv');sp=pd.read_csv(B/'prepared/Visium_lead_gene_spots.csv.gz')
check('Visium_six_animals_twelve_technical_sections',len(sec)==12 and sec.donor.nunique()==6 and sec.groupby('donor').size().eq(2).all())
check('Visium_29839_spots_section_counts',len(sp)==29839 and np.array_equal(sp.groupby('section').size(),sec.set_index('section').n_spots.sort_index()))
check('Visium_unique_section_barcodes_and_count_nonnegative',not sp[['section','barcode']].duplicated().any() and sp.filter(regex='_raw$').ge(0).all().all())
vs=pd.read_csv(B/'tables/Visium_spatial_structure_tests.csv');vf=pd.read_csv(B/'tables/Visium_donor_flight_tests.csv')
check('Visium_full151gene12section_and50set12section_families',vs.groupby('kind').size().to_dict()=={'Hallmark_transcript_sum':600,'panel_gene':1812})
check('Visium_full151gene50set_donor_flight_families',vf.groupby('kind').size().to_dict()=={'Hallmark_transcript_sum':50,'panel_gene':151})
check('Visium_no_spot_pseudoreplication_exposure',vf.n_exact_labels.eq(20).all() and np.allclose(vf.exact_animal_p*20,np.round(vf.exact_animal_p*20)))
dv=pd.read_csv(B/'tables/Visium_donor_mean_scores.csv')
for gene in ['AIF1','PYCARD']:
    s=dv[dv.analyte.eq(gene)];v=s.groupby('flight').mean_spatial_score.mean();effect=v['Space Flight']-v['Ground Control'];r=vf[vf.analyte.eq(gene)].iloc[0]
    check('Visium_'+gene+'_recomputed_donor_effect',np.isclose(effect,r.mean_score_flight_minus_ground))
ca=pd.read_csv(B/'tables/pathway_CAMERA_complete.csv')
check('all50sets_allcontrast_correlation_families',ca.groupby(['family','correlation_model']).size().isin([250,1400,150]).all() and len(ca)==5400)
for (fam,cor),s in ca.groupby(['family','correlation_model']):check('CAMERA_'+fam+'_'+cor+'_BH',np.allclose(s.q_full_family,multipletests(s.PValue,method='fdr_bh')[1]))
check('CAMERA_inadequate_coverage_p1',ca.loc[ca.NGenes.lt(10)].PValue.eq(1).all())
ex=pd.read_csv(B/'tables/pathway_animal_exact_tests.csv');check('exact_singscore_fixed_families',ex.groupby('family').size().to_dict()=={'cells_250':250,'core_flight_400':400,'independent_150':150})
for fam,s in ex.groupby('family'):check('singscore_'+fam+'_BH',np.allclose(s.q_full_family,multipletests(s.exact_p,method='fdr_bh')[1]))
rr=pd.read_csv(B/'tables/OSD536_panel_results.csv');check('RR9_full151panel_family',len(rr)==151);check('RR9_AIF1_absent_PYCARD_below_gate',rr.set_index('analyte').loc['AIF1','status']=='not_in_assay_or_unmapped' and rr.set_index('analyte').loc['PYCARD','status']=='below_fixed_detection_gate')
q=pd.read_csv(B/'tables/OSD536_nCounter_QC.csv');check('RR9_six_flight_six_ground_unique_sources',len(q)==12 and q.animal.is_unique and q.flight.value_counts().to_dict()=={'Space Flight':6,'Ground Control':6})
checksums=pd.read_csv(B/'tables/validation_checksum_verification.csv');check('three_fresh_STAR_deposited_MD5_checks',len(checksums)==3 and checksums.matches.all())
for sid,glds in [(352,352),(563,568),(564,569)]:
    side=B/f'sources/OSD-{sid}/GLDS-{glds}_rna_seq_processed_md5sum_GLbulkRNAseq.tsv';d=dict(line.split('\t')[:2] for line in side.read_text().splitlines() if '\t' in line)
    p=inputs_root()/f'research/brain_integrative_2026/data/OSD-{sid}/GLDS-{glds}_rna_seq_differential_expression_GLbulkRNAseq.csv'
    check('OSD-'+str(sid)+'_annotation_source_full_table_deposited_MD5',hashlib.md5(p.read_bytes()).hexdigest()==d['differential_expression_GLbulkRNAseq.csv'])
pdfs=list((B/'figures').glob('*.pdf'));check('eight_figures_each_PDF_PNG_SVG',len(pdfs)==8 and all(p.with_suffix('.png').exists() and p.with_suffix('.svg').exists() for p in pdfs))
sorted_registry=pd.read_csv(B/'tables/GSE295428_current_sample_registry.csv')
sorted_qc=pd.read_csv(B/'tables/GSE295428_sorted_cell_QC.csv')
labels=pd.read_csv(B/'tables/GSE295428_sample_vs_countfile_labels.csv')
chars=pd.read_csv(B/'tables/GSE295428_current_characteristics_audit.csv')
check('RR8_47_current_GEO_brain_library_labels_match_files',len(labels)==47 and labels.all_labels_match.all() and set(labels.accession)==set(sorted_registry['sample']))
check('RR8_47_known_current_characteristics_after_series_correction',len(chars)==47 and chars.all_characteristics_known_and_match_analysis.all() and not chars.filter(regex='characteristic$').astype(str).eq('unknown').any().any())
cr=chars.set_index('sample').loc[sorted_registry['sample']]
check('RR8_current_tissue_age_treatment_match_fit_registry',np.array_equal(cr.tissue_characteristic,sorted_registry.tissue) and np.array_equal(cr.age_characteristic.astype(int),sorted_registry.age) and np.array_equal(cr.treatment_characteristic.eq('spaceflight'),sorted_registry.flight.eq('FLT')))
check('RR8_eight_shared_group_qualified_exposure_donors',sorted_registry.donor.nunique()==8 and not sorted_registry[['tissue','donor']].duplicated().any() and sorted_registry.groupby('donor').flight.nunique().eq(1).all() and sorted_registry.groupby('donor').age.nunique().eq(1).all())
check('RR8_six_populations_47_libraries',sorted_registry.groupby('tissue').size().sort_values().tolist()==[7,8,8,8,8,8])
check('RR8_156746_barcodes_122084_primary_QC_cells',sorted_qc.n_input_barcodes.sum()==156746 and sorted_qc.n_retained_cells.sum()==122084 and sorted_qc.n_retained_cells.ge(50).all())
check('RR8_marker_cell_gate_nested_primary_QC',sorted_qc.n_marker_positive_cells.le(sorted_qc.n_retained_cells).all() and sorted_qc.loc[sorted_qc.tissue.str.startswith('Myeloid')].n_marker_positive_cells.ge(50).all())
check('RR8_marker_gate_excludes_entire_CNS_target_panel',not set(json.loads((B/'tables/sorted_marker_gate_definition.json').read_text())['marker_genes'])&set(map.mouse_symbol_primary))
provenance=json.loads((B/'sources/sorted_cell_download_provenance.json').read_text())
check('RR8_141_primary_count_barcode_feature_files_CRC_and_SHA',len(provenance)==141 and all(r['gzip_CRC_checked'] and resolve_record(r['file']).stat().st_size==r['bytes'] and hashlib.sha256(resolve_record(r['file']).read_bytes()).hexdigest()==r['SHA256'] for r in provenance))
sa=pd.read_csv(B/'prepared/GSE295428_sorted_annotation.csv')
for name,cell_column in [('donor','n_retained_cells'),('sensitivity','n_sensitivity_cells'),('marker_positive','n_marker_positive_cells')]:
    namefile='GSE295428_sorted_'+('donor' if name=='donor' else name)+'_pseudobulk.csv.gz'
    z=pd.read_csv(B/'prepared'/namefile,index_col=0)
    check('RR8_'+name+'_counts_unique_nonnegative_integer_features',z.index.is_unique and z.columns.is_unique and z.ge(0).all().all() and np.array_equal(z.to_numpy(),np.floor(z.to_numpy())) and np.array_equal(z.index,sa.ENSEMBL))
    check('RR8_'+name+'_47_libraries_registry_order',np.array_equal(z.columns,sorted_registry['sample']))
    if name=='donor':check('RR8_primary_full_counts_match_QC_retained_UMI_totals',np.array_equal(z.sum().to_numpy(),sorted_qc.set_index('sample').loc[z.columns,'retained_UMIs']))
sorted_panel=pd.read_csv(B/'tables/sorted_cell_panel_complete.csv')
check('RR8_3QC_2methods_21744_planned_rows',len(sorted_panel)==21744 and sorted_panel.groupby(['QC','method']).size().eq(3624).all())
for (gate,method),s in sorted_panel.groupby(['QC','method']):
    check('RR8_'+gate+'_'+method+'_independent_full_family_BH',np.allclose(s.q_full_3624_family,multipletests(s.pvalue.fillna(1),method='fdr_bh')[1]))
    check('RR8_'+gate+'_'+method+'_unique_keys',not s[['tissue','contrast','analyte']].duplicated().any())
check('RR8_unavailable_p1_complete_family',sorted_panel.loc[sorted_panel.status.eq('unavailable_or_nonunique')].pvalue.eq(1).all())
for gene in ['AIF1','PYCARD']:
    s=sorted_panel[(sorted_panel.tissue=='Myeloid_Forebrain')&(sorted_panel.contrast=='flight_balanced_age')&(sorted_panel.analyte==gene)]
    check('RR8_forebrain_'+gene+'_all_six_fits_direction_and_q',len(s)==6 and s.log2FC.lt(0).all() and s.q_full_3624_family.lt(.05).all())
sx=pd.read_csv(B/'tables/sorted_cell_exact_tests.csv')
check('RR8_exact_tests_no_cell_pseudoreplication',len(sx)==7248 and sx.n_donors.isin([7,8]).all() and sx.n_exact_labelings.isin([18,36]).all() and np.allclose(sx.exact_p*sx.n_exact_labelings,np.round(sx.exact_p*sx.n_exact_labelings)))
for gate,s in sx.groupby('QC'):check('RR8_exact_'+gate+'_3624_BH',np.allclose(s.q_full_3624_family,multipletests(s.exact_p,method='fdr_bh')[1]))
pair=pd.read_csv(B/'tables/sorted_paired_region_interactions.csv')
check('RR8_direct_region_tests_2718_rows',len(pair)==2718 and pair.n_paired_donors.isin([7,8]).all())
for (gate,error),s in pair.groupby(['QC','error_model']):check('RR8_paired_'+gate+'_'+error+'_453_BH',len(s)==453 and np.allclose(s.q_full_453_family,multipletests(s.pvalue,method='fdr_bh')[1]))
donor=pd.read_csv(B/'tables/sorted_target_donor_expression.csv')
check('RR8_donor_expression_unique_physical_sample_keys',not donor[['QC','tissue','donor','analyte']].duplicated().any())
for gene in ['AIF1','PYCARD']:
    d=donor[(donor.QC=='published_gate')&(donor.analyte==gene)].pivot(index='donor',columns='tissue',values='voom_logCPM')
    md=sorted_registry.drop_duplicates('donor').set_index('donor').loc[d.index]
    # Balanced two-by-two age/exposure design: direct coefficient equals mean of age-specific paired-region exposure differences.
    diff=d.Myeloid_Forebrain-d.Myeloid_Hippocampus
    effects=[]
    for age in [3,8]:
        f=(md.age==age)&(md.flight=='FLT');g=(md.age==age)&(md.flight=='GC');effects.append(diff[f].mean()-diff[g].mean())
    actual=pair[(pair.QC=='published_gate')&(pair.analyte==gene)&(pair.region_difference=='Myeloid_Forebrain minus Myeloid_Hippocampus')].effect.iloc[0]
    check('RR8_'+gene+'_paired_effect_recomputed_from_mouse_points',np.isclose(np.mean(effects),actual))
check('extra_biological_figure_PNG_editable_SVG',all((B/f'figures/figureS2_innate_immune_program.{ext}').exists() for ext in ['png','svg']))
rr8links=pd.read_csv(B/'tables/RR8_NASA_GEO_library_links.csv')
check('RR8_47_libraries_uniquely_match_six_NASA_accessions',len(rr8links)==47 and rr8links.unique_primary_GSM_link.all() and set(rr8links.matched_OSDR_accessions)=={f'OSD-{i}' for i in range(927,933)} and set(rr8links['sample'])==set(sorted_registry['sample']))
pc=pd.read_csv(B/'tables/OSD738_ASC_source_audit.csv')
check('RR1_ASC_20_records_match_submitted_workbook_four_scale_flags',len(pc)==20 and pc.ASC_source_matches_submitted_workbook.all() and pc.ASC_orders_of_magnitude_QC_flag.sum()==4 and set(pc.loc[pc.ASC_orders_of_magnitude_QC_flag,'Sample Name'])=={'M8','M18','M28','M38'})
scamera=pd.read_csv(B/'tables/sorted_pathway_CAMERA_complete.csv')
check('RR8_all50_Hallmark_sets_six_populations_four_contrasts_three_correlations',len(scamera)==3600 and scamera.groupby('correlation_model').size().eq(1200).all() and scamera['set'].nunique()==50)
for cor,s in scamera.groupby('correlation_model'):check('RR8_CAMERA_'+cor+'_1200_BH',np.allclose(s.q_full_1200_family,multipletests(s.PValue,method='fdr_bh')[1]))
sexact=pd.read_csv(B/'tables/sorted_pathway_exact_tests.csv')
check('RR8_Hallmark_exact_1200_whole_donor_labels',len(sexact)==1200 and sexact.n_donors.isin([7,8]).all() and sexact.n_exact_labels.isin([18,36]).all() and np.allclose(sexact.exact_p*sexact.n_exact_labels,np.round(sexact.exact_p*sexact.n_exact_labels)))
check('RR8_Hallmark_exact_independent_full_BH',np.allclose(sexact.q_full_1200_family,multipletests(sexact.exact_p,method='fdr_bh')[1]))
deletions=pd.read_csv(B/'tables/sorted_leave_one_mouse_effects.csv');deletion_models=pd.read_csv(B/'tables/sorted_leave_one_mouse_models.csv')
check('RR8_mouse_deletion_23_designs_22_estimable',len(deletion_models)==23 and deletion_models.status.eq('refitted').sum()==22 and deletion_models.status.eq('missing_age_exposure_design_cell').sum()==1)
check('RR8_mouse_deletion_all151_targets_two_engines',len(deletions)==23*151*2 and not deletions[['tissue','removed_donor','method','analyte']].duplicated().any())
for g,negative in [('AIF1',True),('PYCARD',True),('AIM2',False),('NINJ1',True)]:
    s=deletions[(deletions.tissue=='Myeloid_Forebrain')&(deletions.analyte==g)]
    check('RR8_'+g+'_forebrain_direction_all16_actual_deletion_refits',len(s)==16 and s.log2FC.notna().all() and (s.log2FC.lt(0) if negative else s.log2FC.gt(0)).all())
relative=pd.read_csv(B/'tables/sorted_marker_relative_tests.csv');deep=pd.read_csv(B/'tables/sorted_cell_depth_tests.csv')
for frame,qcol,size,label in [(relative,'q_full_1812_family',1812,'marker_relative'),(deep,'q_full_1359_family',1359,'cell_depth')]:
    check('RR8_'+label+'_complete_inference_families',len(frame)==2*size and frame.groupby('error_model').size().eq(size).all())
    for error,s in frame.groupby('error_model'):check('RR8_'+label+'_'+error+'_independent_BH',np.allclose(s[qcol],multipletests(s.pvalue,method='fdr_bh')[1]))
    cb=frame[(frame.tissue=='Myeloid_Cerebellum')&(frame.error_model=='HC3_t')]
    check('RR8_'+label+'_singleton_HC3_not_falsely_estimable',cb.pvalue.eq(1).all() and cb.ci_low.isna().all() and cb.ci_high.isna().all() and not cb.status.eq('tested').any())
dc=pd.read_csv(B/'tables/sorted_cell_depth_QC.csv');de=pd.read_csv(B/'tables/sorted_cell_depth_donor_endpoints.csv')
check('RR8_depth_summaries_23_libraries_8_shared_mice',len(dc)==23 and dc.donor.nunique()==8 and np.array_equal(dc.set_index('sample').n_primary_cells,sorted_qc.set_index('sample').loc[dc['sample'],'n_retained_cells']))
check('RR8_500UMI_threshold_retains_same_primary_cells',dc.n_primary_cells.eq(dc.n_500_UMI_cells).all())
check('RR8_depth_endpoint_all151_targets_3metrics_23samples',len(de)==151*3*23 and not de[['sample','analyte','endpoint']].duplicated().any())
dp=de.pivot(index=['sample','analyte'],columns='endpoint',values='value').dropna()
check('RR8_detection_probabilities_bounded_and_depth_expectation_below_raw',dp.raw_detection_fraction.between(0,1).all() and dp.expected_detection_at_500_UMIs.between(0,1).all() and (dp.expected_detection_at_500_UMIs<=dp.raw_detection_fraction+1e-10).all())
from scipy.special import gammaln
from scipy.stats import hypergeom
N=np.array([500,700,1000,5000,20000]);K=np.array([1,2,5,20,499]);valid=N-K>=500;analytic=np.ones(len(N));aN=N[valid];aK=K[valid];analytic[valid]=-np.expm1(gammaln(aN-aK+1)-gammaln(aN-aK-500+1)+gammaln(aN-500+1)-gammaln(aN+1))
check('depth_rarefaction_formula_matches_independent_hypergeom_distribution',np.allclose(analytic,hypergeom.sf(0,N,K,500),rtol=1e-8,atol=1e-9))
refvals=pd.read_csv(B/'tables/sorted_marker_relative_donor_values.csv');rawcounts=pd.read_csv(B/'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz',index_col=0);ann=sa.set_index('SYMBOL')
anchor_ids=[ann.loc[g,'ENSEMBL'] for g in ['Csf1r','Tmem119','P2ry12','Hexb','Sall1']]
for g in ['Aif1','Pycard']:
    sample='GSM8947901';expected=np.log2(rawcounts.loc[ann.loc[g,'ENSEMBL'],sample]+.5)-np.log2(rawcounts.loc[anchor_ids,sample]+.5).mean();actual=refvals[(refvals['sample']==sample)&refvals.analyte.eq(g.upper())].marker_relative_log2_ratio.iloc[0]
    check('RR8_'+g+'_reference_ratio_reconstructed_from_raw_counts',np.isclose(expected,actual))
pd.DataFrame(checks).to_csv(B/'tables/verification_checks.csv',index=False)
(B/'tables/verification_summary.json').write_text(json.dumps({'checks':len(checks),'passed':sum(x['passes'] for x in checks),'failed':sum(not x['passes'] for x in checks)},indent=2));print('PASSED',len(checks),'independent identity and numerical checks',flush=True)
