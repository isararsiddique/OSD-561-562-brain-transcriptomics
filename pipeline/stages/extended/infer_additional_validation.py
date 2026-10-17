from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Non-shrunk mouse-level marker-relative and detection/depth validation."""
from pathlib import Path
import json
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.stats import false_discovery_control
B=stage_root('extended');V1=stage_root('core')
panel=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('');ann=pd.read_csv(B/'prepared/GSE295428_sorted_annotation.csv');lookup=ann.groupby('SYMBOL').ENSEMBL.agg(list).to_dict()
counts=pd.read_csv(B/'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz',index_col=0);md=pd.read_csv(B/'tables/GSE295428_sorted_cell_QC.csv').set_index('sample');endpoints=pd.read_csv(B/'tables/sorted_cell_depth_donor_endpoints.csv')
anchors=['Csf1r','Tmem119','P2ry12','Hexb','Sall1'];assert not set(anchors)&set(panel.mouse_symbol_primary)
anchor_ids=[lookup[x][0] for x in anchors];assert all(len(lookup[x])==1 for x in anchors)
contrasts={'flight_age3':[-1,0,1,0],'flight_age8':[0,-1,0,1],'flight_balanced_age':[-.5,-.5,.5,.5],'flight_by_age':[1,-1,-1,1]}
ratio_rows=[];ratio_values=[];endpoint_rows=[]
def test(y,D,C,ok):
    if ok:
        fit=sm.OLS(y,D).fit();leverage=np.sum(D*(D@np.linalg.inv(D.T@D)),axis=1)
        robust_available=not np.any(np.isclose(leverage,1))
        variance_available=fit.ssr>1e-20
        robust=fit.get_robustcov_results(cov_type='HC3',use_t=True) if robust_available and variance_available else None
    for model in ['ordinary_OLS','HC3_t']:
        p=1.;effect=lo=hi=np.nan;status='unavailable_or_sparse'
        if ok:effect=float(np.asarray(C)@fit.params)
        if ok and not variance_available:status='zero_residual_variance_no_inference'
        elif ok and model=='HC3_t' and not robust_available:status='HC3_unestimable_singleton_design_cell'
        elif ok:
            t=(fit if model=='ordinary_OLS' else robust).t_test(C);effect=float(np.asarray(t.effect).item());p=float(np.asarray(t.pvalue).item());ci=np.asarray(t.conf_int()).ravel();lo,hi=map(float,ci)
            status='tested'
            if not np.isfinite(p):p=1.;status='nonfinite_inference';lo=hi=np.nan
        yield {'error_model':model,'effect':effect,'ci_low':lo,'ci_high':hi,'pvalue':p,'status':status}
for tissue in sorted(t for t in md.tissue.unique() if t.startswith('Myeloid')):
    s=md[md.tissue==tissue];group=s.flight+'_'+s.age.astype(str);D=np.column_stack([group.eq(g).astype(float) for g in ['GC_3','GC_8','FLT_3','FLT_8']]);assert np.linalg.matrix_rank(D)==4
    denom=np.log2(counts.loc[anchor_ids,s.index]+.5).mean(axis=0).to_numpy();ids=s.index
    for p in panel.itertuples():
        hit=lookup.get(p.mouse_symbol_primary,[]);ok=len(hit)==1 and (counts.loc[hit[0],ids]>=10).sum()>=3
        if ok:
            y=np.log2(counts.loc[hit[0],ids].to_numpy()+.5)-denom
            for sample,val in zip(ids,y):ratio_values.append({'sample':sample,'donor':s.loc[sample,'donor'],'tissue':tissue,'flight':s.loc[sample,'flight'],'age':s.loc[sample,'age'],'analyte':p.analyte,'marker_relative_log2_ratio':val})
        for tag,C in contrasts.items():
            for r in test(y if ok else None,D,C,ok):ratio_rows.append(dict(r,tissue=tissue,analyte=p.analyte,mouse_symbol=p.mouse_symbol_primary,contrast=tag,n_donors=len(s)))
        e=endpoints[(endpoints.tissue==tissue)&(endpoints.analyte==p.analyte)]
        for endpoint,v in e.groupby('endpoint'):
            v=v.set_index('sample').loc[ids];valid=ok and v.value.notna().all() and v.n_cells.ge(50).all()
            for r in test(v.value.to_numpy(),D,contrasts['flight_balanced_age'],valid):endpoint_rows.append(dict(r,tissue=tissue,analyte=p.analyte,mouse_symbol=p.mouse_symbol_primary,endpoint=endpoint,n_donors=len(s)))
r=pd.DataFrame(ratio_rows);r['q_full_1812_family']=r.groupby('error_model').pvalue.transform(false_discovery_control);r.to_csv(B/'tables/sorted_marker_relative_tests.csv',index=False)
pd.DataFrame(ratio_values).to_csv(B/'tables/sorted_marker_relative_donor_values.csv',index=False)
e=pd.DataFrame(endpoint_rows);e['q_full_1359_family']=e.groupby('error_model').pvalue.transform(false_discovery_control);e.to_csv(B/'tables/sorted_cell_depth_tests.csv',index=False)
(B/'tables/sorted_marker_relative_definition.json').write_text(json.dumps({'anchors':anchors,'response':'log2(target pseudobulk count+0.5) - mean(log2(anchor pseudobulk count+0.5))','scope':'Relative transcript stoichiometry to biological myeloid denominators, not absolute per-cell expression or validated housekeeping normalization','design':'four age/exposure groups; all151 targets x3 myeloid regions x4 contrasts per error model','exposure_unit':'mouse'},indent=2))
print('Marker-relative HC3 hits',int(r.query('error_model=="HC3_t"').q_full_1812_family.lt(.05).sum()))
print('Cell-depth HC3 hits',int(e.query('error_model=="HC3_t"').q_full_1359_family.lt(.05).sum()))
print(r[(r.tissue=='Myeloid_Forebrain')&(r.contrast=='flight_balanced_age')&r.analyte.isin(['AIF1','PYCARD','AIM2','NINJ1'])].to_string(index=False))
print(e[(e.tissue=='Myeloid_Forebrain')&e.analyte.isin(['AIF1','PYCARD','AIM2','NINJ1'])&(e.error_model=='HC3_t')].to_string(index=False))
