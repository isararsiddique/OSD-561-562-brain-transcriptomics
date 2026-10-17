from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Exact age-stratified donor tests and direct paired regional interactions."""
from pathlib import Path
import itertools,json
import numpy as np,pandas as pd
import statsmodels.api as sm
from scipy.stats import false_discovery_control
B=stage_root('extended');V1=stage_root('core')
mapping=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('');ann=pd.read_csv(B/'prepared/GSE295428_sorted_annotation.csv')
symbols=ann.groupby('SYMBOL').ENSEMBL.agg(list).to_dict();qc=pd.read_csv(B/'tables/GSE295428_sorted_cell_QC.csv')
rows=[];pairrows=[];vals=[]
def age_masks(md):
    indices=[np.flatnonzero(md.age.eq(age).to_numpy()) for age in [3,8]]
    choices=[list(itertools.combinations(ix,int(md.iloc[ix].flight.eq('FLT').sum()))) for ix in indices]
    return [(np.isin(np.arange(len(md)),sum((list(v) for v in combination),[]))) for combination in itertools.product(*choices)]
for gate in ['published_gate','alternative_gate']:
    matrices={};metadata={}
    for tissue,s in qc.groupby('tissue'):
        p=B/f'prepared/Sorted_{tissue.replace("-","_")}_{gate}_voom_logCPM.csv.gz';Y=pd.read_csv(p,index_col=0)
        md=s.set_index('sample').loc[Y.columns].copy();md['sample']=md.index;metadata[tissue]=md;matrices[tissue]=Y
        masks=age_masks(md);f=md.flight.eq('FLT').to_numpy()
        def effects(y,mask):
            e=[]
            for age in [3,8]:
                ix=md.age.eq(age).to_numpy();e.append(y[mask&ix].mean()-y[(~mask)&ix].mean())
            return np.array([e[0],e[1],np.mean(e),e[1]-e[0]])
        for r in mapping.itertuples():
            hits=symbols.get(r.mouse_symbol_primary,[]);valid=len(hits)==1 and hits[0] in Y.index
            if valid:
                y=Y.loc[hits[0]].to_numpy();observed=effects(y,f);null=np.array([effects(y,mask) for mask in masks]);p=np.mean(np.abs(null)>=np.abs(observed)-1e-12,axis=0)
                for i,row in enumerate(md.itertuples()):vals.append({'QC':gate,'tissue':tissue,'sample':row.Index,'donor':row.donor,'flight':row.flight,'age':row.age,'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'voom_logCPM':y[i]})
            for k,tag in enumerate(['flight_age3','flight_age8','flight_balanced_age','flight_by_age']):rows.append({'QC':gate,'tissue':tissue,'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'contrast':tag,'n_donors':len(md),'n_exact_labelings':len(masks),'effect':observed[k] if valid else np.nan,'exact_p':p[k] if valid else 1.,'status':'tested' if valid else 'unavailable'})
    for region1,region0 in [('Myeloid_Cerebellum','Myeloid_Hippocampus'),('Myeloid_Forebrain','Myeloid_Hippocampus'),('Myeloid_Forebrain','Myeloid_Cerebellum')]:
        a=metadata[region1].set_index('donor');b=metadata[region0].set_index('donor');donors=sorted(set(a.index)&set(b.index));a=a.loc[donors];b=b.loc[donors];assert np.array_equal(a.age,b.age) and np.array_equal(a.flight,b.flight)
        design=sm.add_constant(pd.DataFrame({'flight':a.flight.eq('FLT').astype(float),'old':a.age.eq(8).astype(float)},index=donors));masks=age_masks(a)
        for r in mapping.itertuples():
            ids=symbols.get(r.mouse_symbol_primary,[]);valid=len(ids)==1 and ids[0] in matrices[region1].index and ids[0] in matrices[region0].index
            if valid:
                y=matrices[region1].loc[ids[0],a['sample']].to_numpy()-matrices[region0].loc[ids[0],b['sample']].to_numpy();fit=sm.OLS(y,design).fit();rob=fit.get_robustcov_results(cov_type='HC3',use_t=True)
                # Refit the same age-adjusted difference model for every whole-donor flight labeling.
                null=[]
                for mask in masks:
                    D=design.copy();D['flight']=mask.astype(float);null.append(sm.OLS(y,D).fit().params['flight'])
                pv=np.mean(abs(np.array(null))>=abs(fit.params['flight'])-1e-12)
            for model in ['ordinary_OLS','HC3_t','exact_age_stratified']:
                row={'QC':gate,'region_difference':region1+' minus '+region0,'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'n_paired_donors':len(donors),'effect':fit.params['flight'] if valid else np.nan,'pvalue':(fit.pvalues['flight'] if model=='ordinary_OLS' else rob.pvalues[1] if model=='HC3_t' else pv) if valid else 1.,'error_model':model,'status':'paired_region_by_flight_association' if valid else 'unavailable'}
                pairrows.append(row)
a=pd.DataFrame(rows);a['q_full_3624_family']=a.groupby('QC').exact_p.transform(false_discovery_control);a.to_csv(B/'tables/sorted_cell_exact_tests.csv',index=False)
d=pd.DataFrame(pairrows);d['q_full_453_family']=d.groupby(['QC','error_model']).pvalue.transform(lambda p:false_discovery_control(p.fillna(1)));d.to_csv(B/'tables/sorted_paired_region_interactions.csv',index=False)
pd.DataFrame(vals).to_csv(B/'tables/sorted_target_donor_expression.csv',index=False)
print(d[d.analyte.isin(['AIF1','PYCARD'])].to_string(index=False),flush=True)
