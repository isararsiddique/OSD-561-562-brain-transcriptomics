from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Disclosed post-result checks of secondary composition-adjusted target associations."""
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.api as sm
from scipy.stats import false_discovery_control
B=stage_root('extended');V1=stage_root('core')
map=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('');E=pd.read_csv(B/'tables/bulk_RNA_fraction_estimates.csv')
rows=[]
for sid in [561,562]:
    md=pd.read_csv(V1/f'tables/OSD-{sid}_sample_metadata.csv',index_col=0);ix=md.index[md.group.isin(['FLT_ISS_T_12','GC_ISS_T_12'])]
    ann=pd.read_csv(V1/f'prepared/OSD-{sid}_annotation.csv').set_index('ENSEMBL')
    raw=pd.read_csv(V1/f'prepared/OSD-{sid}_STAR_endogenous.csv.gz',index_col=0);raw.index=ann.reindex(raw.index).SYMBOL;raw=raw[raw.index.notna()].groupby(level=0).sum()
    cpm=raw.div(raw.sum(),axis=1)*1e6
    normalized=pd.read_csv(V1/f'prepared/OSD-{sid}_normalized_counts.csv.gz',index_col=0);normalized.index=ann.reindex(normalized.index).SYMBOL;normalized=normalized[normalized.index.notna()].groupby(level=0).sum()
    for method in ['NNLS','ridge_NNLS']:
        immune=E[(E.study==f'OSD-{sid}')&(E.cell_type=='Immune cells')&(E.method==method)].set_index('sample').loc[ix,'estimated_RNA_fraction']
        X=sm.add_constant(pd.DataFrame({'flight':md.loc[ix,'flight'].eq('FLT').astype(float),'immune_RNA':immune},index=ix))
        for norm,Y in [('CPM',cpm),('DESeq2_normalized_counts',normalized)]:
            for r in map.itertuples():
                valid=r.mouse_symbol_primary in Y.index
                if valid:
                    y=np.log2(Y.loc[r.mouse_symbol_primary,ix]+1);fit=sm.OLS(y,X).fit();rob=fit.get_robustcov_results(cov_type='HC3',use_t=True)
                    loo=[]
                    for j in range(len(ix)):
                        keep=np.arange(len(ix))!=j;loo.append(sm.OLS(y.iloc[keep],X.iloc[keep]).fit().params['flight'])
                for se in ['ordinary_OLS','HC3_t_residual_df']:
                    row={'study':f'OSD-{sid}','method':method,'normalization':norm,'error_model':se,'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'pvalue':1.,'effect':np.nan,'LOO_effect_min':np.nan,'LOO_effect_max':np.nan,'max_Cooks_distance':np.nan,'status':'unavailable'}
                    if valid:row.update(pvalue=fit.pvalues['flight'] if se=='ordinary_OLS' else rob.pvalues[1],effect=fit.params['flight'],LOO_effect_min=min(loo),LOO_effect_max=max(loo),max_Cooks_distance=fit.get_influence().cooks_distance[0].max(),status='secondary_conditional_association')
                    rows.append(row)
a=pd.DataFrame(rows);a['q_full_302_family']=a.groupby(['method','normalization','error_model']).pvalue.transform(lambda p:false_discovery_control(p.fillna(1)))
a.to_csv(B/'tables/adjusted_panel_normalization_HC3_LOO.csv',index=False)
print(a[a.analyte.isin(['IRF1','NOS3','SOD1','SPHK2','UCP2'])&a.study.eq('OSD-561')].to_string(index=False),flush=True)
