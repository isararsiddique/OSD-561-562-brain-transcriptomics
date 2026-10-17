from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Panel-excluded reference RNA mixtures, held-donor validation and secondary adjustment."""
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy.optimize import nnls
from scipy.stats import rankdata,false_discovery_control
import statsmodels.api as sm
B=stage_root('extended');V1=stage_root('core')
rng=np.random.default_rng(20261016)
panel=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('')
anno=pd.read_csv(B/'prepared/OSD352_RNA_annotation.csv').set_index('ENSEMBL')
raw=pd.read_csv(B/'prepared/OSD352_RNA_pseudobulk.csv.gz',index_col=0)
raw.index=anno.loc[raw.index,'SYMBOL'];raw=raw.groupby(level=0).sum()
ref=raw/raw.sum(axis=0)*1e6
registry=pd.read_csv(B/'tables/OSD352_l2_registry.csv')
types=sorted(registry.cell_type.unique());donors=sorted(registry.bio_origin.unique())
assert len(types)==5 and len(donors)==5
shared=set(ref.index)-set(panel.mouse_symbol_primary)-{''}
for sid in [561,562,352,563,564]:
    src=V1 if sid in [561,562] else B
    a=pd.read_csv(src/f'prepared/OSD-{sid}_annotation.csv');shared &= set(a.SYMBOL.dropna())
shared=sorted(g for g in shared if not g.lower().startswith(('mt-','rpl','rps')))
def reference(ds):
    return pd.DataFrame({ct:ref[[f'{d}|{ct}' for d in ds]].mean(axis=1) for ct in types})
def markers(R):
    sub=R.loc[shared];selected=set()
    for ct in types:
        specificity=(sub[ct]+.1)/(sub.drop(columns=ct).max(axis=1)+.1)
        keep=sub[ct].ge(10)&specificity.ge(2)
        selected.update(specificity[keep].sort_values(ascending=False).head(80).index)
    return sorted(selected)
def estimate(R,y,genes,method):
    A=R.loc[genes].to_numpy();v=y.loc[genes].to_numpy();w=1/np.sqrt(np.maximum(A.max(axis=1),1))
    A=A*w[:,None];v=v*w
    if method=='ridge_NNLS':
        lam=.05*np.sum(A*A)/len(types);A=np.vstack([A,np.eye(len(types))*np.sqrt(lam)]);v=np.r_[v,np.zeros(len(types))]
    x,_=nnls(A,v);return x/x.sum() if x.sum()>0 else np.full(len(types),np.nan)
ground=registry.loc[registry.flight_status.eq('ground'),'bio_origin'].unique().tolist()
cv=[]
for hold in donors:
    training=[d for d in ground if d!=hold]
    R=reference(training);genes=markers(R)
    H=reference([hold])
    for sim in range(100):
        truth=rng.dirichlet(np.ones(len(types)))
        # RNA fractions, not cell fractions; independent donor profile and multinomial counting noise.
        expectation=H.to_numpy()@truth;draw=rng.multinomial(1000000,expectation/expectation.sum())
        y=pd.Series(draw,index=H.index)/draw.sum()*1e6
        for method in ['NNLS','ridge_NNLS']:
            pred=estimate(R,y,genes,method)
            for i,ct in enumerate(types):cv.append({'held_out_donor':hold,'training_ground_donors':';'.join(training),'simulation':sim,'method':method,'cell_type':ct,'true_RNA_fraction':truth[i],'estimated_RNA_fraction':pred[i],'n_training_markers':len(genes)})
cv=pd.DataFrame(cv);cv.to_csv(B/'tables/deconvolution_held_donor_validation.csv',index=False)
cv['absolute_error']=abs(cv.true_RNA_fraction-cv.estimated_RNA_fraction)
cv.groupby(['method','cell_type']).absolute_error.agg(['mean','median','max']).to_csv(B/'tables/deconvolution_validation_summary.csv')
# Fit final reference to ground donors only; target panel never used to select signatures.
ground=registry.loc[registry.flight_status.eq('ground'),'bio_origin'].unique().tolist();R=reference(ground);genes=markers(R)
R.loc[genes].to_csv(B/'prepared/deconvolution_ground_reference.csv');pd.Series(genes,name='signature_gene').to_csv(B/'tables/deconvolution_signature_genes.csv',index=False)
assert not set(genes)&set(panel.mouse_symbol_primary)
est=[];adjust=[]
for sid in [561,562,352,563,564]:
    src=V1 if sid in [561,562] else B
    counts=pd.read_csv(src/f'prepared/OSD-{sid}_STAR_endogenous.csv.gz',index_col=0)
    a=pd.read_csv(src/f'prepared/OSD-{sid}_annotation.csv').set_index('ENSEMBL')
    counts.index=a.reindex(counts.index).SYMBOL;counts=counts.loc[counts.index.notna()].groupby(level=0).sum()
    cpm=counts/counts.sum(axis=0)*1e6
    md=pd.read_csv(src/f'tables/OSD-{sid}_{"sample_metadata" if sid in [561,562] else "metadata"}.csv',index_col=0)
    for method in ['NNLS','ridge_NNLS']:
        E=pd.DataFrame({s:estimate(R,cpm[s],genes,method) for s in cpm.columns},index=types).T
        for s,row in E.iterrows():
            for ct in types:est.append({'study':f'OSD-{sid}','sample':s,'animal':md.loc[s,'animal'],'flight':md.loc[s,'flight'],'group':md.loc[s,'group'],'method':method,'cell_type':ct,'estimated_RNA_fraction':row[ct]})
        if sid not in [561,562]:continue
        ix=md.index[md.group.isin(['FLT_ISS_T_12','GC_ISS_T_12'])];f=md.loc[ix,'flight'].eq('FLT').astype(float)
        immune=E.loc[ix,'Immune cells'];X=sm.add_constant(pd.DataFrame({'flight':f,'estimated_immune_RNA':immune},index=ix))
        assert np.linalg.matrix_rank(X)==3
        for r in panel.itertuples():
            row={'study':f'OSD-{sid}','method':method,'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'n_animals':len(ix),'pvalue':1.,'flight_effect_adjusted':np.nan,'ci_low':np.nan,'ci_high':np.nan,'flight_effect_unadjusted':np.nan,'flight_immune_correlation':float(np.corrcoef(f,immune)[0,1]),'test_status':'unavailable'}
            if r.mouse_symbol_primary in cpm.index:
                y=np.log2(cpm.loc[r.mouse_symbol_primary,ix]+1)
                fit=sm.OLS(y,X).fit();u=sm.OLS(y,sm.add_constant(f)).fit();ci=fit.conf_int().loc['flight']
                row.update(pvalue=fit.pvalues['flight'],flight_effect_adjusted=fit.params['flight'],ci_low=ci.iloc[0],ci_high=ci.iloc[1],flight_effect_unadjusted=u.params['flight'],test_status='secondary_OLS_reference_RNA_covariate')
            adjust.append(row)
pd.DataFrame(est).to_csv(B/'tables/bulk_RNA_fraction_estimates.csv',index=False)
ad=pd.DataFrame(adjust);ad['q_full_family']=ad.groupby('method').pvalue.transform(lambda p:false_discovery_control(p.fillna(1)))
ad.to_csv(B/'tables/composition_adjusted_young_terminal_panel.csv',index=False)
# Author-labelled captured nuclei proportions: sample-processing/capture composition, not absolute brain cell counts.
n=registry.pivot(index='bio_origin',columns='cell_type',values='n_nuclei');frac=n.div(n.sum(axis=1),axis=0)
f=registry.drop_duplicates('bio_origin').set_index('bio_origin').loc[frac.index].flight_status.eq('flight').to_numpy()
rows=[]
for ct in types:
    v=frac[ct].to_numpy();effect=v[f].mean()-v[~f].mean();null=[]
    for ids in itertools.combinations(range(5),3):
        m=np.isin(np.arange(5),ids);null.append(v[m].mean()-v[~m].mean())
    rows.append({'cell_type':ct,'flight_minus_ground_fraction':effect,'exact_p':np.mean(np.abs(null)>=abs(effect)-1e-12)})
t=pd.DataFrame(rows);t['q_five_types']=false_discovery_control(t.exact_p);t.to_csv(B/'tables/captured_nucleus_composition_exact.csv',index=False)
frac.to_csv(B/'tables/captured_nucleus_fractions.csv')
summary={'signature_genes':len(genes),'reference_ground_donors':ground,'CV_held_out_donors':5,'synthetic_mixtures_per_donor':100,'CV_mean_absolute_error':cv.groupby('method').absolute_error.mean().to_dict(),'all_fractions':'RNA fractions under a whole-brain reference model, not absolute cell fractions; out-of-domain cerebellum/hippocampus and capture effects remain.','adjustment':'Secondary log2(CPM+1) OLS, eight-group source reduced to young terminal subgroup; panel excluded from signatures; not mediation or handling control.'}
(B/'tables/composition_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
