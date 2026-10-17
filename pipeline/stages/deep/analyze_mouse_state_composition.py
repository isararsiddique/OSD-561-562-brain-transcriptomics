from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Mouse-level proportions and transparent fixed-composition sensitivity."""
from pathlib import Path
from itertools import combinations, product
import json, sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
B=stage_root('deep')
root=B
if '--depth-adjusted' in sys.argv:B=root/'depth_sensitivity'
qc=pd.read_csv(B/'tables/state_mouse_QC.csv')
v=pd.read_csv(B/'tables/state_target_mouse_expression.csv')
panel=pd.read_csv(stage_root('core') / 'tables/verified_panel_mapping.csv')
tests=['flight_age3','flight_age8','flight_balanced_age','flight_by_age']
C=np.array([[-1,0,1,0],[0,-1,0,1],[-.5,-.5,.5,.5],[1,-1,-1,1]])
def design(s):
    return np.column_stack([(s.flight.eq(f)&s.age.eq(a)).to_numpy().astype(float)
                            for f,a in [('GC',3),('GC',8),('FLT',3),('FLT',8)]])
def exact_effects(y,s):
    strata=[np.flatnonzero(s.age.eq(a)) for a in [3,8]]
    out=[]
    for chosen in product(*[list(combinations(idx,2)) for idx in strata]):
        flight=np.zeros(len(s),bool);flight[np.concatenate(chosen)]=True
        sp=s.copy();sp['flight']=np.where(flight,'FLT','GC')
        out.append(C@np.linalg.lstsq(design(sp),y,rcond=None)[0])
    return np.array(out)
def fits(y,s):
    f=sm.OLS(y,design(s)).fit(cov_type='HC3',use_t=True)
    perm=exact_effects(y,s)
    ans=[]
    for i,c in enumerate(C):
        t=f.t_test(c);effect=float(np.asarray(t.effect).item())
        lo,hi=np.asarray(t.conf_int()).ravel()
        ans.append({'contrast':tests[i],'effect':effect,'ci_low':float(lo),'ci_high':float(hi),
                    'HC3_pvalue':float(np.asarray(t.pvalue).item()),
                    'exact_pvalue':float(np.mean(np.abs(perm[:,i])>=abs(effect)-1e-12)),
                    'exact_label_assignments':len(perm)})
    return ans
rows=[];values=[];props=[];coverage=[]
for k in [4,6]:
    s=qc[qc.K.eq(k)].drop_duplicates('donor').sort_values('donor')[['sample','donor','flight','age']].reset_index(drop=True)
    counts=qc[qc.K.eq(k)].pivot(index='donor',columns='state',values='n_cells').reindex(s.donor)
    proportions=counts.div(counts.sum(axis=1),axis=0)
    common=counts.columns[counts.min(axis=0)>=50].tolist()
    all_observed=bool((counts>0).all().all())
    coverage.append({'K':k,'all_states_observed_in_each_mouse':all_observed,
                     'all_states_have_50_cells_in_each_mouse':bool((counts>=50).all().all()),
                     'minimum_state_mouse_cells':int(counts.min().min()),'common_support_states':common,
                     'uniform_all_state_interpretation':'Sparse-state exploratory sensitivity; requires >=1 observed cell in every mouse/state. Not reliable subtype adjustment.' if all_observed else 'Unestimable: at least one state has no cells in a mouse; no imputation performed.'})
    for state in range(1,k+1):
        y=proportions[state].to_numpy()*100
        for result in fits(y,s):props.append(dict(K=k,state=state,units='percentage points',**result))
    for p in panel.itertuples():
        z=v[v.K.eq(k)&v.analyte.eq(p.analyte)].pivot(index='donor',columns='state',values='mean_cell_CP10K').reindex(s.donor)
        ok=np.isfinite(z.to_numpy()).all()
        # A state with no cells contributes no counts to the observed mixture.
        # This is not an estimate/imputation of that absent state's expression.
        observed_ok=np.all(np.isfinite(z.to_numpy()) | (counts.to_numpy()==0))
        common_ok=bool(common) and np.isfinite(z[common].to_numpy()).all()
        ys={'observed_composition':np.log2((z*proportions).sum(axis=1).to_numpy()+1) if observed_ok else None,
            'uniform_all_states':np.log2(z.mean(axis=1).to_numpy()+1) if ok and all_observed else None,
            'uniform_common_support':np.log2(z[common].mean(axis=1).to_numpy()+1) if common_ok else None}
        for measure,y in ys.items():
            if y is None:
                for contrast in tests:rows.append({'K':k,'analyte':p.analyte,'mouse_symbol':p.mouse_symbol_primary,'measure':measure,'contrast':contrast,
                                                   'effect':np.nan,'ci_low':np.nan,'ci_high':np.nan,'HC3_pvalue':1.,'exact_pvalue':1.,'exact_label_assignments':36,'status':'unavailable_no_complete_observations'})
            else:
                for result in fits(y,s):rows.append(dict(K=k,analyte=p.analyte,mouse_symbol=p.mouse_symbol_primary,measure=measure,status='tested',**result))
                for j,donor in enumerate(s.donor):values.append({'K':k,'analyte':p.analyte,'measure':measure,'donor':donor,'flight':s.flight.iloc[j],'age':int(s.age.iloc[j]),'log2_mean_cell_CP10K_plus1':y[j]})
a=pd.DataFrame(rows)
for k in [4,6]:
    idx=a.K.eq(k);assert idx.sum()==1812
    for p,q in [('HC3_pvalue','HC3_q_full_1812_family'),('exact_pvalue','exact_q_full_1812_family')]:
        a.loc[idx,q]=multipletests(a.loc[idx,p],method='fdr_bh')[1]
a.to_csv(B/'tables/state_composition_complete_tests.csv',index=False)
pd.DataFrame(values).to_csv(B/'tables/state_composition_mouse_values.csv',index=False)
props=pd.DataFrame(props)
for p,q in [('HC3_pvalue','HC3_q_full_40_family'),('exact_pvalue','exact_q_full_40_family')]:
    props[q]=multipletests(props[p],method='fdr_bh')[1]
assert len(props)==40
props.to_csv(B/'tables/state_proportion_mouse_tests.csv',index=False)
(B/'tables/state_composition_coverage.json').write_text(json.dumps(coverage,indent=2))
print(a[a.analyte.isin(['AIF1','PYCARD','AIM2','NINJ1'])&a.contrast.eq('flight_balanced_age')][['K','analyte','measure','effect','HC3_q_full_1812_family','exact_q_full_1812_family','status']].to_string(index=False))
