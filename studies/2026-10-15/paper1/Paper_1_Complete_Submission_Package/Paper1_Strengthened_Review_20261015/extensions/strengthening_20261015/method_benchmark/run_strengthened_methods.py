#!/usr/bin/env python3
"""Evaluate the certified constrained fit and independent-cohort prediction.

The additional cohorts test the algorithm on independently sampled data. Their
exposure effects are not reclassified as artifacts and are not removed for a
biological discovery claim. All tuning stays within each outer training set.
"""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import run_method_benchmark as b
from constrained_regression import solve_constrained

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
ALPHAS=[None,10.,1.,.1,.01]

def external_data(n):
    directory=ROOT/f'extensions/strengthening_20261015/external_data/OSD-{n}'
    counts=pd.read_csv(directory/'counts.csv.gz',index_col=0)
    meta=pd.read_csv(directory/'metadata.csv',dtype={'age':str}).set_index('sample_id')
    meta=meta.loc[meta.primary_analysis_included.astype(str).str.lower().eq('true')].copy()
    counts=counts.loc[:,meta.index]
    meta['animal']=meta.independent_unit_id
    if n=='352':label=meta.spaceflight.eq('Space Flight').astype(int)
    else:label=(~meta.radiation.str.contains('non-irradiated|sham',case=False,na=False)).astype(int)
    if label.nunique()!=2:raise ValueError(f'Invalid exposure labels for{n}')
    design=[np.ones(len(meta))]
    if n=='202':
        design += [pd.to_numeric(meta.collection_after_treatment).eq(4).to_numpy(float),
                   meta.hindlimb_unloading.eq('Hindlimb Unloaded').to_numpy(float)]
    meta['label']=label.to_numpy();meta.index=meta.animal.to_numpy()
    return dict(n=n,c=counts.to_numpy(float),genes=counts.index.to_numpy(),meta=meta,
                age_full=np.column_stack(design),canonical=counts.index.str.match(r'^ENSMUSG\d+$'),
                source_hash=hashlib.sha256(counts.to_numpy().tobytes()).hexdigest())

class ConstrainedCache:
    def __init__(self):self.cache={};self.certificates=[]
    def get(self,data,m):
        key=(data['n'],tuple(sorted(m['train'])))
        if key not in self.cache:
            fits={str(None):np.zeros((len(m['sensor']),len(m['train'])))}
            for alpha in ALPHAS[1:]:
                Q,certificate=solve_constrained(m['S'],m['Y'],m['scale'],alpha=alpha)
                self.certificates.append(dict(accession='OSD-'+data['n'],training_units='|'.join(data['meta'].iloc[m['train']].animal),**certificate))
                if not certificate['converged']:raise ValueError('Constrained fit did not meet declared tolerance')
                fits[str(alpha)]=Q
            self.cache[key]=fits
        return self.cache[key]

def constrained_tune(data,train,base_cache,cache):
    errors=np.zeros(len(ALPHAS));baseline=0.
    for left in train:
        inner=np.array([i for i in train if i!=left]);m=base_cache.get(data,inner);t=b.target(data,m,[left]);fits=cache.get(data,m);baseline+=t['b0']
        for j,alpha in enumerate(ALPHAS):errors[j]+=b.risk(m,t,t['S']@fits[str(alpha)])
    choice=int(np.argmin(errors))
    return ALPHAS[choice],errors,baseline

def main():
    protocol={
      'date':'2026-10-15','status':'exploratory prospective computational amendment after earlier benchmark; not preregistration',
      'algorithm':'Exact convex gain-constrained ridge, distinct from clipping unconstrained ridge coefficients',
      'objective':'||Y-SC||F^2/(2n)+alpha*||D C||F^2/2 subject ||C||2<=sqrt(3)',
      'hyperparameters':ALPHAS,'gain_cap':2,'optimality_gap_relative_tolerance':1e-7,
      'outer_original':'All362 previously specified ground holdouts',
      'outer_independent':'Leave one independent unit out in all design-selected OSD202,352,374,611; use NASA SampleTable retained units',
      'independent_unit':'202/352 animals;374/611 two-animal pools. No pool member is a separate replicate.',
      'nuisance':'202 intercept,time4months,unloading; other independent datasets intercept; original age specification unchanged',
      'label_for_original_score':'352 flight;202/374/611 radiation vs sham. These are exposure-associated predictors, not established nuisance signatures.',
      'nested':'Every filter, normalization reference, scaling, regression, and hyperparameter choice uses inner-training units only',
      'evaluation':'All nine original comparators plus exact constrained ridge in independent data; fixed original outputs retained',
      'inference':'Outer split and repeated-unit distributions are descriptive, not independent split-level confidence intervals',
      'biological_scope':'Independent-cohort prediction tests algorithm behavior. Predicting expression is not validation of removing genuine exposure biology.'}
    (OUT/'strengthened_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    start=time.time();base_cache=b.TrainingCache(size=128);cache=ConstrainedCache();rows=[];innerrows=[];units=[];splitroles=[]
    oldroles=pd.read_csv(ROOT/'extensions/ground_holdout/split_animal_roles.csv')
    for n in ['561','562','525','202','352','374','611']:
        original=n in ['561','562','525'];data=b.get_data(n) if original else external_data(n)
        if original:
            group=oldroles[oldroles.accession.eq('OSD-'+n)];splits=[]
            for sid,part in group.groupby('split_id'):
                train=np.sort([data['meta'].index.get_loc(a) for a in part[part.role.eq('train')].animal]);test=np.sort([data['meta'].index.get_loc(a) for a in part[part.role.eq('test')].animal]);splits.append((sid,train,test))
        else:splits=[(i,np.delete(np.arange(len(data['meta'])),i),np.array([i])) for i in range(len(data['meta']))]
        for sid,train,test in splits:
            m=base_cache.get(data,train);t=b.target(data,m,test);alpha,errors,base=constrained_tune(data,train,base_cache,cache);Q=cache.get(data,m)[str(alpha)]
            selected={'spectral_constrained_ridge_CV':dict(M=t['S']@Q,gain=b.gain_from_factors(Q,m['YY']),parameter=str(alpha))}
            if not original:selected.update(b.selected_methods(m,t,b.tune(data,train,base_cache)))
            for method,f in selected.items():
                prediction=f['M']@m['Y'];error=t['R']-prediction;sse=float(np.sum(error**2));ratio=sse/t['b0'];ratio=1. if abs(ratio-1)<1e-12 else ratio
                rows.append(dict(accession='OSD-'+n,scope='original_ground' if original else 'independent_exposure_cohort',split_id=int(sid),method=method,n_train=len(train),n_test=len(test),n_sensor=len(m['sensor']),n_genes=len(m['idx']),MSE_ratio=ratio,baseline_SSE=t['b0'],prediction_SSE=sse,operator_gain=f['gain'],parameter=f['parameter']))
                for j,i in enumerate(test):units.append(dict(accession='OSD-'+n,scope='original_ground' if original else 'independent_exposure_cohort',split_id=int(sid),method=method,unit=data['meta'].iloc[i].animal,baseline_MSE=float(np.mean(t['R'][j]**2)),prediction_MSE=float(np.mean(error[j]**2))))
            for a,e in zip(ALPHAS,errors):innerrows.append(dict(accession='OSD-'+n,split_id=int(sid),alpha=str(a),inner_SSE=e,inner_baseline_SSE=base,selected=a==alpha))
            for role,ids in [('train',train),('test',test)]:
                for i in ids:splitroles.append(dict(accession='OSD-'+n,split_id=int(sid),role=role,unit=data['meta'].iloc[i].animal))
            if sid%20==0:print(n,sid,len(splits),'elapsed',time.time()-start,flush=True)
        pd.DataFrame(rows).to_csv(OUT/'strengthened_outer_results.csv',index=False);pd.DataFrame(units).to_csv(OUT/'strengthened_unit_errors.csv',index=False);pd.DataFrame(innerrows).to_csv(OUT/'strengthened_inner_CV.csv',index=False);pd.DataFrame(cache.certificates).to_csv(OUT/'convex_optimality_certificates.csv',index=False);pd.DataFrame(splitroles).to_csv(OUT/'strengthened_split_roles.csv',index=False)
        cache.cache.clear()
    frame=pd.DataFrame(rows);summary=[]
    for (n,scope,method),part in frame.groupby(['accession','scope','method']):
        v=part.MSE_ratio.to_numpy();u=part.operator_gain.to_numpy();summary.append(dict(accession=n,scope=scope,method=method,n_splits=len(part),mean_MSE_ratio=float(v.mean()),pooled_SSE_ratio=float(part.prediction_SSE.sum()/part.baseline_SSE.sum()),median_MSE_ratio=float(np.median(v)),q025_MSE_ratio=float(np.quantile(v,.025)),q975_MSE_ratio=float(np.quantile(v,.975)),fraction_improved=float(np.mean(v<1-1e-12)),mean_gain=float(u.mean()),maximum_gain=float(u.max())))
    pd.DataFrame(summary).to_csv(OUT/'strengthened_summary.csv',index=False)
    (OUT/'strengthened_run_summary.json').write_text(json.dumps({'status':'PASS','n_outer_rows':len(rows),'n_independent_cohort_units':64,'n_certified_fits':len(cache.certificates),'all_convex_fits_converged':True,'max_relative_gap':max(x['relative_gap'] for x in cache.certificates),'elapsed_seconds':time.time()-start},indent=2)+'\n')
    print(pd.DataFrame(summary).to_string(index=False),flush=True)

if __name__=='__main__':main()
