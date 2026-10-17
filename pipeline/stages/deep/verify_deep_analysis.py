from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Independent family, count-conservation, identity and figure-value audits."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
B=stage_root('deep');R=stage_root('extended');checks=[]
def check(name,value):
    checks.append({'check':name,'passes':bool(value)})
    assert value,name
def family(data,keys,q,size,p='pvalue'):
    for label,s in data.groupby(keys):
        check(f'{q} {label}: planned family',len(s)==size)
        check(f'{q} {label}: independent BH',np.allclose(s[q],multipletests(s[p],method='fdr_bh')[1],rtol=1e-9,atol=1e-12))
        check(f'{q} {label}: finite bounded p/q',np.isfinite(s[[p,q]]).all().all() and s[[p,q]].ge(0).all().all() and s[[p,q]].le(1).all().all())
old=pd.read_csv(R/'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz',index_col=0)
for directory in [B,B/'depth_sensitivity']:
    qc=pd.read_csv(directory/'tables/state_mouse_QC.csv')
    cells=pd.read_csv(directory/'prepared/forebrain_state_cell_registry.csv.gz')
    check(str(directory)+' eight mice/23036 cells',len(cells)==23036 and cells.donor.nunique()==8 and not cells[['sample','barcode']].duplicated().any())
    features=pd.read_csv(directory/'tables/state_nonpanel_variable_features.csv')
    panel=pd.read_csv(stage_root('core') / 'tables/verified_panel_mapping.csv')
    check(str(directory)+' 2000 non-panel features',len(features)==2000 and not features.SYMBOL.isin(panel.mouse_symbol_primary).any())
    for k in [4,6]:
        counts=pd.read_csv(directory/'prepared'/f'state_K{k}_pseudobulk.csv.gz',index_col=0)
        check(str(directory)+f' K{k} integer counts',np.isfinite(counts).all().all() and counts.ge(0).all().all() and np.equal(counts,np.floor(counts)).all().all())
        for sample in cells['sample'].unique():
            cols=qc[qc.K.eq(k)&qc['sample'].eq(sample)].key.tolist()
            check(str(directory)+f' K{k} {sample} exact count conservation',np.array_equal(counts[cols].sum(axis=1),old[sample]))
        a=pd.read_csv(directory/'tables/state_panel_complete.csv');a=a[a.K.eq(k)]
        family(a,['method'],'q_full_state_family',151*k*4)
        check(str(directory)+f' K{k} unique state/gene/contrast keys',not a[['method','state','analyte','contrast']].duplicated().any())
        check(str(directory)+f' K{k} missing state estimates retained',a.loc[a.status.ne('tested'),'pvalue'].eq(1).all() and a.loc[a.status.ne('tested'),'log2FC'].isna().all())
    comp=pd.read_csv(directory/'tables/state_composition_complete_tests.csv')
    family(comp,['K'],'HC3_q_full_1812_family',1812,'HC3_pvalue')
    family(comp,['K'],'exact_q_full_1812_family',1812,'exact_pvalue')
    check(str(directory)+' exact whole-mouse assignment resolution',comp.exact_label_assignments.eq(36).all() and np.allclose(comp.exact_pvalue*36,np.round(comp.exact_pvalue*36)))
    values=pd.read_csv(directory/'tables/state_composition_mouse_values.csv')
    raw=values[values.measure.eq('observed_composition')].pivot(index=['analyte','donor'],columns='K',values='log2_mean_cell_CP10K_plus1')
    check(str(directory)+' observed mixture invariant to state partition',np.allclose(raw[4],raw[6],equal_nan=True,rtol=1e-10,atol=1e-10))
    model=pd.read_csv(directory/'tables/state_model_QC.csv')
    for r in model.itertuples():
        ss=qc[qc.K.eq(r.K)&qc.state.eq(r.state)]
        check(str(directory)+f' K{r.K} S{r.state} all-mouse gate',bool(r.model_available)==bool(len(ss)==8 and ss.n_cells.ge(50).all()))
    for record in json.loads((directory/'tables/state_preparation_summary.json').read_text())['raw_input_hashes']:
        check('Primary source hash '+resolve_record(record['path']).name,hashlib.sha256(resolve_record(record['path']).read_bytes()).hexdigest()==record['sha256'])
first=pd.read_csv(B/'tables/additional_panel_complete.csv');family(first,['method'],'q_full_1057_family',1057)
regions=pd.read_csv(B/'tables/terminal_region_panel_complete.csv');family(regions,['method'],'q_full_755_family',755)
expanded=pd.read_csv(B/'tables/additional_panel_all_regions.csv');family(expanded,['method'],'q_expanded_1812_family',1812)
check('Expanded complete 151-panel keys',not expanded[['OSD','method','contrast','analyte']].duplicated().any())
check('Expanded unavailable estimates never zero-filled',expanded.loc[expanded.status.ne('tested'),'log2FC'].isna().all() and expanded.loc[expanded.status.ne('tested'),'pvalue'].eq(1).all())
md=pd.read_csv(B/'tables/additional_mouse_registry.csv')
check('RR23 29 libraries collapsed to26mice',len(md[md.OSD.eq('OSD-525')])==26 and md[md.OSD.eq('OSD-525')].n_technical_libraries.sum()==29)
check('MHU3 24 unique cerebrum mice',len(md[md.OSD.eq('OSD-457')])==24 and md[md.OSD.eq('OSD-457')].donor.is_unique)
for sid in [525,457]:
    c=pd.read_csv(B/'prepared'/f'OSD-{sid}_mouse_counts.csv.gz',index_col=0)
    check(f'OSD{sid} count totals match primary linked registry',np.array_equal(c.sum().to_numpy(),md.set_index('key').loc[c.columns,'raw_mouse_counts'].to_numpy()))
ind=json.loads((B/'tables/terminal_region_animal_independence.json').read_text())
check('RRRM2 18 new animals, RR10 nine shared bulk animals',ind[0]['new_mice_added']==0 and ind[0]['n_animals']==9 and ind[1]['new_mice_added']==18 and ind[1]['overlap_with_original_region_animals']==0)
missing=pd.read_csv(B/'tables/terminal_region_missing_bulk_samples.csv')
check('RR10 FL09 absent from bulk, retained as missing record',len(missing)==1 and missing.donor.iloc[0]=='RR-10_FL-09')
delete=pd.read_csv(B/'tables/additional_mouse_deletion_expanded_complete.csv.gz')
check('Every50newcohort mouse deletions x2engines',delete[['deleted_OSD','deleted_mouse']].drop_duplicates().shape[0]==50 and len(delete)==50*2*1812)
family(delete,['deleted_OSD','deleted_mouse','method'],'q_expanded_1812_family',1812)
article=pd.read_csv(B/'tables/article_5pct_panel_complete.csv');family(article,['method'],'q_full_3624_family',3624)
aq=pd.read_csv(B/'tables/article_5pct_mouse_QC.csv')
check('Article5percent gate 47libraries eightmice 62351cells',len(aq)==47 and aq.donor.nunique()==8 and aq.n_cells.sum()==62351)
ac=pd.read_csv(B/'prepared/article_5pct_donor_pseudobulk.csv.gz',index_col=0)
check('Article5percent count conservation toQC totals',np.array_equal(ac.sum().to_numpy(),aq.set_index('sample').loc[ac.columns,'retained_UMIs']))
for r in pd.read_csv(B/'tables/figure8_displayed_count_effects.csv').itertuples():
    data=pd.read_csv((B if r.analysis=='Primary' else B/'depth_sensitivity')/'tables/state_panel_complete.csv')
    s=data[data.K.eq(r.K)&data.state.eq(r.state)&data.method.eq(r.method)&data.analyte.eq(r.analyte)&data.contrast.eq('flight_balanced_age')]
    check('Figure8 displayed '+str(r),len(s)==1 and np.allclose([s.log2FC.iloc[0],s.q_full_state_family.iloc[0]],[r.log2FC,r.q]))
for r in pd.read_csv(B/'tables/figure9_displayed_contrasts.csv').itertuples():
    s=expanded[expanded.OSD.eq(r.OSD)&expanded.method.eq(r.method)&expanded.contrast.eq(r.contrast)&expanded.analyte.eq(r.analyte)]
    check('Figure9 displayed '+str(r),len(s)==1 and np.allclose([s.log2FC.iloc[0],s.q_expanded_1812_family.iloc[0]],[r.log2FC,r.q]))
for r in pd.read_csv(B/'tables/figure9_displayed_mouse_values.csv').itertuples():
    E=pd.read_csv(B/'prepared'/f'{r.OSD}_voom_logCPM.csv.gz',index_col=0);ann=pd.read_csv(B/'prepared/additional_gene_annotation.csv')
    ids=ann.loc[ann.SYMBOL.eq(r.gene)&ann.ENSEMBL.isin(E.index),'ENSEMBL'];sample=md[md.OSD.eq(r.OSD)&md.donor.eq(r.donor)].key.iloc[0]
    check('Figure9 mouse '+r.OSD+r.donor+r.gene,len(ids)==1 and np.isclose(E.loc[ids.iloc[0],sample],r.voom_logCPM))
for r in pd.read_csv(B/'tables/figure10_displayed_effects.csv').itertuples():
    if r.source=='additional_panel_all_regions':
        s=expanded[expanded.OSD.eq(r.OSD)&expanded.analyte.eq(r.analyte)&expanded.contrast.eq(r.contrast)&expanded.method.eq(r.method)];qcol='q_expanded_1812_family'
    else:
        data=pd.read_csv((R if r.source=='sorted_cell_panel_complete' else B)/'tables'/(r.source+'.csv'))
        s=data[data.tissue.eq(r.tissue)&data.analyte.eq(r.analyte)&data.contrast.eq(r.contrast)&data.method.eq(r.method)]
        if r.source=='sorted_cell_panel_complete':s=s[s.QC.eq(r.QC)]
        qcol='q_full_3624_family'
    check('Figure10 effect '+str(r),len(s)==1 and np.allclose([s.log2FC.iloc[0],s[qcol].iloc[0]],[r.log2FC,r.q]))
tr=pd.read_csv(B/'tables/terminal_region_mouse_registry.csv');E=pd.read_csv(B/'prepared/OSD-613_voom_logCPM.csv.gz',index_col=0)
ann=pd.read_csv(B/'prepared/additional_gene_annotation.csv')
for r in pd.read_csv(B/'tables/figure10_displayed_mouse_values.csv').itertuples():
    ids=ann.loc[ann.SYMBOL.eq(r.gene)&ann.ENSEMBL.isin(E.index),'ENSEMBL'];sample=tr[tr.OSD.eq(r.OSD)&tr.donor.eq(r.donor)].key.iloc[0]
    check('Figure10 mouse '+r.donor,len(ids)==1 and np.isclose(E.loc[ids.iloc[0],sample],r.voom_logCPM))
s3=pd.read_csv(B/'tables/figureS3_displayed_complete_effects.csv')
check('FigureS3 all 151x12x2 effects and q exact',len(s3)==3624 and np.allclose(s3.log2FC,expanded.log2FC,equal_nan=True) and np.allclose(s3.q_expanded_1812_family,expanded.q_expanded_1812_family))
cons=pd.read_csv(B/'tables/forebrain_all_eight_QC_method_consensus.csv')
oldsorted=pd.read_csv(R/'tables/sorted_cell_panel_complete.csv');oldsorted=oldsorted[oldsorted.tissue.eq('Myeloid_Forebrain')&oldsorted.contrast.eq('flight_balanced_age')]
newsorted=article[article.tissue.eq('Myeloid_Forebrain')&article.contrast.eq('flight_balanced_age')]
full=pd.concat([oldsorted,newsorted]).groupby('analyte').q_full_3624_family.agg(['size','max'])
check('Eight method/QC complete consensus recalculated',len(cons)==151 and all(row.all_eight_corrected==bool(full.loc[row.analyte,'size']==8 and full.loc[row.analyte,'max']<.05) for row in cons.itertuples()))
check('44 distinct mouse genes survive all eight fits',sum(cons.all_eight_corrected)==44 and cons.loc[cons.all_eight_corrected,'mouse_symbol_primary'].nunique()==44)
report={'checks':len(checks),'passed':sum(c['passes'] for c in checks),'details':checks,'scope':'Independent re-calculation, exact raw-count conservation, biological identity and no missing-effect imputation'}
(B/'tables/deep_verification.json').write_text(json.dumps(report,indent=2))
print(f"Deep validation: {report['passed']}/{report['checks']} passed")
