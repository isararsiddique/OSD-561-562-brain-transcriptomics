from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Independent source/design/statistical consistency checks for the study."""
from pathlib import Path
import json,pickle,hashlib,gzip
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
OUT=stage_root('core');ROOT=inputs_root();TAB=OUT/'tables';PRE=OUT/'prepared'
checks=[]
def check(name,condition,detail=''):
    checks.append({'check':name,'passed':bool(condition),'detail':detail})
    assert condition,name
b=pd.read_csv(TAB/'bulk_panel_complete.csv');s=pd.read_csv(TAB/'spatial_panel_complete.csv');m=pd.read_csv(TAB/'bulk_sample_registry.csv');mapp=pd.read_csv(TAB/'verified_panel_mapping.csv').fillna('')
check('corrected_panel_151_unique_human_ids',len(mapp)==151 and mapp.human_ensembl.nunique()==151)
check('primary_mapping_143_unique_mouse_assignments',mapp.mouse_symbol_primary.ne('').sum()==143 and mapp.loc[mapp.mouse_symbol_primary.ne(''),'mouse_symbol_primary'].nunique()==143)
check('bulk_full_family_4228',len(b)==4228)
check('spatial_full_family_4228',len(s)==4228)
for name,f in [('bulk',b),('spatial',s)]:
    q=multipletests(f.pvalue.fillna(1),method='fdr_bh')[1]
    check(name+'_BH_independent_library_agreement',np.allclose(q,f.q_panel_primary_full_family,atol=1e-12))
check('bulk_17_primary_two_flight',b.q_panel_primary_full_family.lt(.05).sum()==17 and ((b.contrast_kind=='flight') & b.q_panel_primary_full_family.lt(.05)).sum()==2)
check('spatial_primary_and_Welch_zero',s.q_panel_primary_full_family.lt(.05).sum()==0 and s.q_welch_panel_full_family.lt(.05).sum()==0)
check('53_libraries_29_animals_24_paired',len(m)==53 and m.animal.nunique()==29 and m.groupby('animal').study.nunique().eq(2).sum()==24)
check('source_identity_consistent_between_regions',m.groupby('animal')[['flight','age','cohort']].nunique().to_numpy().max()==1)
check('within_region_unique_animals',not m.duplicated(['study','animal']).any())
check('seven_location_conflicts',m.location_conflict.sum()==7)
for osd in [561,562]:
    with (PRE/f'OSD-{osd}_STAR_dds.pkl').open('rb') as f:dds=pickle.load(f)
    counts=pd.read_csv(PRE/f'OSD-{osd}_STAR_endogenous.csv.gz',index_col=0).T
    counts=counts.loc[:,counts.ge(10).sum(axis=0).ge(3)]
    check(f'{osd}_cache_raw_count_identity',np.array_equal(dds.X,counts.loc[dds.obs_names,dds.var_names].to_numpy()))
    check(f'{osd}_only_endogenous_genes',pd.Index(dds.var_names).str.startswith('ENSMUSG').all())
    X=np.asarray(dds.obsm['design_matrix']);lfc=np.asarray(dds.varm['LFC']);groups=dds.obs.group.astype(str).to_numpy()
    check(f'{osd}_eight_cell_design_full_rank',np.linalg.matrix_rank(X)==8)
    for tag,g in b[b.study==f'OSD-{osd}'].groupby('contrast'):
        def row(name):return X[groups==name].mean(axis=0)
        r=g.iloc[0]
        if r.contrast_kind=='interaction':
            age=tag.split('age')[1];c=row('FLT_ISS_T_'+age)-row('GC_ISS_T_'+age)-row('FLT_LAR_'+age)+row('GC_LAR_'+age)
        else:c=row(r.numerator)-row(r.denominator)
        expected=lfc@c/np.log(2)
        hit=g[g.gene_id.isin(dds.var_names)];ix=pd.Index(dds.var_names).get_indexer(hit.gene_id)
        check(f'{osd}_{tag}_effect_from_fitted_coefficients',np.allclose(expected[ix],hit.log2FoldChange,atol=1e-9))
for osd,code in [(682,'CA'),(685,'DG'),(698,'FCT'),(699,'CT')]:
    a=pd.read_csv(ROOT/'OSD-682-685-698-699-brain-spatial/data/geo'/f'GSE239336_{code}_GeneExpression_Q3norm.txt.gz',sep='\t',index_col=0)
    y=np.log2(a.loc['Gsdmd'].to_numpy());fl=np.array([c.startswith('FLT') for c in a.columns]);dr=np.array(['_TX_' in c for c in a.columns])
    direct=y[fl&dr].mean()-y[~fl&dr].mean()-y[fl&~dr].mean()+y[~fl&~dr].mean()
    actual=s[(s.study==f'OSD-{osd}') & (s.analyte=='GSDMD') & (s.contrast=='flight_by_BuOE_interaction')].iloc[0]
    check(f'{osd}_full_DiD_independent_means',np.isclose(direct,actual.log2FoldChange,atol=1e-10))
    w=ttest_ind(y[fl&~dr],y[~fl&~dr],equal_var=False)
    hit=s[(s.study==f'OSD-{osd}')&(s.analyte=='GSDMD')&(s.contrast=='flight_in_saline')].iloc[0]
    check(f'{osd}_Welch_independent_scipy_check',np.isclose(w.pvalue,hit.welch_pvalue,atol=1e-10) and np.isclose(w.df,hit.welch_df))
checksum=pd.read_csv(TAB/'primary_checksum_verification.csv')
check('both_bulk_primary_deposited_MD5_match',checksum.matches_deposited_md5.all())
ext=json.loads((TAB/'OSD-352_context_summary.json').read_text())
check('external_five_animals_no_confirmation',ext['n_flight_animals']==3 and ext['n_ground_animals']==2 and ext['n_q05']==0)
(TAB/'verification_checks.json').write_text(json.dumps({'all_passed':all(c['passed'] for c in checks),'n_checks':len(checks),'checks':checks},indent=2))
print('Passed',len(checks),'source, design and statistical consistency checks.')
