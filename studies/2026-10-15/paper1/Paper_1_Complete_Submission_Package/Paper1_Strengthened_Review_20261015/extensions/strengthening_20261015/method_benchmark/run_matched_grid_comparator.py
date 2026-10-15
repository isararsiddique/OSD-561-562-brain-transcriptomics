#!/usr/bin/env python3
"""Nested matched-grid ridge/clipped-ridge comparison to the completed convex fit."""
import hashlib,json,time
import numpy as np
import pandas as pd
import run_method_benchmark as b
from run_strengthened_methods import external_data,ALPHAS
from constrained_regression import solve_constrained

OUT=b.OUT;ROOT=b.ROOT

def main():
    start=time.time();roles=pd.read_csv(OUT/'strengthened_split_roles.csv');cache=b.TrainingCache();rows=[];innerrows=[];fixed=[];unitrows=[]
    convex=pd.read_csv(OUT/'strengthened_outer_results.csv');convex=convex[convex.method=='spectral_constrained_ridge_CV'].set_index(['accession','split_id'])
    for n in ['561','562','525','202','352','374','611']:
        data=b.get_data(n) if n in ['561','562','525'] else external_data(n);rr=roles[roles.accession=='OSD-'+n]
        for sid,r in rr.groupby('split_id',sort=True):
            train=np.sort([data['meta'].index.get_loc(a) for a in r[r.role=='train'].unit]);test=np.sort([data['meta'].index.get_loc(a) for a in r[r.role=='test'].unit]);m=cache.get(data,train);t=b.target(data,m,test);errors={'ridge':np.zeros(len(ALPHAS)),'clipped_ridge':np.zeros(len(ALPHAS))};base=0
            for left in train:
                mi=cache.get(data,np.array([i for i in train if i!=left]));ti=b.target(data,mi,[left]);base+=ti['b0']
                for name,fits in [('ridge',mi['ridge']),('clipped_ridge',mi['bounded_ridge'])]:
                    for j,a in enumerate(ALPHAS):errors[name][j]+=b.risk(mi,ti,ti['S']@fits[str(a)])
            for key,method,fits in [('ridge','matched_grid_ridge_CV',m['ridge']),('clipped_ridge','matched_grid_gain2_clipped_ridge_CV',m['bounded_ridge'])]:
                selected=int(np.argmin(errors[key]));a=ALPHAS[selected];Q=fits[str(a)];P=t['S']@Q@m['Y'];E=t['R']-P;sse=float(np.sum(E*E));ratio=sse/t['b0'];ratio=1. if abs(ratio-1)<1e-12 else ratio
                ref=convex.loc[('OSD-'+n,sid)];rows.append(dict(accession='OSD-'+n,scope=ref.scope,split_id=int(sid),method=method,n_train=len(train),n_test=len(test),MSE_ratio=ratio,baseline_SSE=t['b0'],prediction_SSE=sse,operator_gain=b.gain_from_factors(Q,m['YY']),selected_alpha=str(a),convex_MSE_ratio=float(ref.MSE_ratio),convex_minus_comparator_MSE_ratio=float(ref.MSE_ratio-ratio)))
                for j,i in enumerate(test):unitrows.append(dict(accession='OSD-'+n,split_id=int(sid),method=method,unit=data['meta'].iloc[i].animal,baseline_MSE=float(np.mean(t['R'][j]**2)),prediction_MSE=float(np.mean(E[j]**2))))
                for j,a in enumerate(ALPHAS):innerrows.append(dict(accession='OSD-'+n,split_id=int(sid),method=method,alpha=str(a),inner_SSE=float(errors[key][j]),inner_baseline_SSE=base,selected=j==selected))
            if sid==rr.split_id.min():
                for a in ALPHAS[1:]:
                    Qc,cert=solve_constrained(m['S'],m['Y'],m['scale'],alpha=a);Qp=m['bounded_ridge'][str(a)]
                    def obj(Q):
                        C=Q@m['Y'];loss=np.sum((m['Y']-m['S']@C)**2)/(2*len(train));penalty=a*np.sum((m['scale'][:,None]*C)**2)/2;return float(loss+penalty),float(np.linalg.svd(C,compute_uv=False)[0])
                    oc,nc=obj(Qc);op,np_=obj(Qp);bound=cert['Frank_Wolfe_duality_gap']+1e-9*max(1,cert['objective_baseline'])
                    fixed.append(dict(accession='OSD-'+n,split_id=int(sid),alpha=a,constrained_objective=oc,clipped_ridge_objective=op,constrained_minus_clipped_objective=oc-op,constrained_spectral_norm=nc,clipped_spectral_norm=np_,certified_gap=cert['Frank_Wolfe_duality_gap'],converged=cert['converged'],verified_dominance_with_gap=bool(oc<=op+bound),convex_objective_reconstruction_error=abs(oc-cert['objective'])))
            if sid%40==0:print(n,sid,'elapsed',round(time.time()-start,1),flush=True)
        pd.DataFrame(rows).to_csv(OUT/'matched_grid_outer_results.csv',index=False);pd.DataFrame(innerrows).to_csv(OUT/'matched_grid_inner_CV.csv',index=False);pd.DataFrame(fixed).to_csv(OUT/'matched_grid_fixed_alpha_objectives.csv',index=False);pd.DataFrame(unitrows).to_csv(OUT/'matched_grid_unit_errors.csv',index=False)
    frame=pd.DataFrame(rows);summary=[]
    for (n,scope,method),part in frame.groupby(['accession','scope','method']):
        diff=part.convex_minus_comparator_MSE_ratio.to_numpy();v=part.MSE_ratio.to_numpy();summary.append(dict(accession=n,scope=scope,method=method,n_splits=len(part),mean_MSE_ratio=float(v.mean()),pooled_SSE_ratio=float(part.prediction_SSE.sum()/part.baseline_SSE.sum()),mean_operator_gain=float(part.operator_gain.mean()),convex_mean_MSE_ratio=float(part.convex_MSE_ratio.mean()),paired_convex_minus_comparator_mean=float(diff.mean()),paired_difference_median=float(np.median(diff)),paired_difference_q025=float(np.quantile(diff,.025)),paired_difference_q975=float(np.quantile(diff,.975)),fraction_convex_lower_error=float(np.mean(diff< -1e-12)),fraction_comparator_improved=float(np.mean(v<1-1e-12))))
    pd.DataFrame(summary).to_csv(OUT/'matched_grid_summary.csv',index=False);ff=pd.DataFrame(fixed)
    status=dict(status='PASS' if ff.verified_dominance_with_gap.all() and ff.converged.all() else 'FAIL',n_outer_rows=len(frame),n_outer_holdouts=frame[['accession','split_id']].drop_duplicates().shape[0],n_fixed_alpha_checks=len(ff),all_fixed_alpha_dominance_verified=bool(ff.verified_dominance_with_gap.all()),max_fixed_alpha_convex_minus_clipped_objective=float(ff.constrained_minus_clipped_objective.max()),max_objective_reconstruction_error=float(ff.convex_objective_reconstruction_error.max()),elapsed_seconds=time.time()-start,protocol_sha256=hashlib.sha256((OUT/'matched_grid_protocol.json').read_bytes()).hexdigest())
    (OUT/'matched_grid_run_summary.json').write_text(json.dumps(status,indent=2)+'\n');print(pd.DataFrame(summary).to_string(index=False),flush=True);print(json.dumps(status,indent=2),flush=True);assert status['status']=='PASS'

if __name__=='__main__':main()
