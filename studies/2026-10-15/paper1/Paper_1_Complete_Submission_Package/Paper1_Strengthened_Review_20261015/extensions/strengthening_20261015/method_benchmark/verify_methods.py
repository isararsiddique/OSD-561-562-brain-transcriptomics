#!/usr/bin/env python3
"""Independent dense regression, sklearn-PLS2 and gain checks."""
from pathlib import Path
import importlib.util,json
import numpy as np
import pandas as pd
import run_method_benchmark as r


def run():
    from sklearn.cross_decomposition import PLSRegression
    rng=np.random.default_rng(901531);rows=[]
    for repeat in range(20):
        n,s,g=12,7,30;X=rng.normal(size=(n,s));X-=X.mean(0);scale=np.exp(rng.normal(size=s));X=X/scale;Y=rng.normal(size=(n,g))+X@rng.normal(size=(s,g));Y-=Y.mean(0);YY=Y@Y.T
        for alpha in [.01,1.,100.]:
            Q=r.ridge_factors(X,scale,alpha);coef=np.linalg.solve(X.T@X+n*alpha*np.eye(s),X.T@Y)/scale[:,None];err=float(np.max(abs(Q@Y-coef)));rows.append(dict(check='ridge_dual_vs_dense_primal',case=repeat,parameter=alpha,error=err,passed=err<1e-10))
        pcr,rank=r.pcr_factors(X,scale);pls=r.pls2_factors(X,YY)
        for k in [1,2,3,4]:
            U,d,Vt=np.linalg.svd(X,full_matrices=False);score=X@Vt[:k].T;theta=np.linalg.lstsq(score,Y,rcond=None)[0];coef=Vt[:k].T@theta/scale[:,None];err=float(np.max(abs(pcr[k]@Y-coef)));rows.append(dict(check='PCR_vs_explicit_score_regression',case=repeat,parameter=k,error=err,passed=err<1e-10))
            model=PLSRegression(n_components=k,scale=False,tol=1e-13,max_iter=5000).fit(X,Y);test=rng.normal(size=(5,s));ours=test@pls[k]@Y;reference=model.predict(test);err=float(np.max(abs(ours-reference)));rows.append(dict(check='PLS2_covariance_deflation_vs_sklearn_NIPALS',case=repeat,parameter=k,error=err,passed=err<1e-5))
        Q=r.ridge_factors(X,scale,.01);B=Y.T@Q.T;A=np.block([[np.zeros((s,s)),np.zeros((s,g))],[-B,np.eye(g)]]);gain=float(np.linalg.svd(A,compute_uv=False)[0]);err=abs(gain-r.gain_from_factors(Q,YY));rows.append(dict(check='multi_sensor_dense_operator_gain',case=repeat,parameter=0,error=err,passed=err<1e-10))
        Qc=r.bound_factors(Q,YY);Bc=Y.T@Qc.T;Ac=np.block([[np.zeros((s,s)),np.zeros((s,g))],[-Bc,np.eye(g)]]);err=float(np.linalg.norm(Ac@Ac-Ac));b=np.r_[np.zeros(s),rng.normal(size=g)];rows.append(dict(check='bounded_projection_idempotence',case=repeat,parameter=0,error=err,passed=err<1e-10));rows.append(dict(check='bounded_sensor_null_biology_preservation',case=repeat,parameter=0,error=float(np.linalg.norm(Ac@b-b)),passed=np.allclose(Ac@b,b)))
        gainc=float(np.linalg.svd(Ac,compute_uv=False)[0]);rows.append(dict(check='bounded_dense_gain_cap2',case=repeat,parameter=0,error=max(0,gainc-2),passed=gainc<=2+1e-10))
        U,ds,Vt=np.linalg.svd(B,full_matrices=False);closest=(U*np.minimum(ds,np.sqrt(3)))@Vt;err=float(np.max(abs(Bc-closest)));rows.append(dict(check='factor_bound_vs_direct_SVD_spectral_projection',case=repeat,parameter=0,error=err,passed=err<1e-10))
        w=rng.normal(size=s+g);w/=np.linalg.norm(w);h=rng.normal(size=s+g);h-=w*(w@h);h*=8/np.linalg.norm(h);gamma=w+h;lam=rng.random();gl=w+lam*h;As=np.eye(len(w))-np.outer(gl,w);err=float(np.max(abs(np.linalg.svd(As,compute_uv=False)[0]-np.sqrt(1+lam*lam*64))));rows.append(dict(check='risk_anchored_dense_gain',case=repeat,parameter=lam,error=err,passed=err<1e-10));rows.append(dict(check='risk_anchored_exact_score_removal',case=repeat,parameter=lam,error=float(np.linalg.norm(w@As)),passed=np.linalg.norm(w@As)<1e-10))
    # Actual count preprocessing checked independently against the previous train_reference implementation.
    spec=importlib.util.spec_from_file_location('prior_ground_reference',r.ROOT/'extensions/ground_holdout/run_ground_holdout.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    for n in ['561','562','525']:
        data=r.get_data(n);roles=pd.read_csv(r.ROOT/'extensions/ground_holdout/split_animal_roles.csv');s=roles[(roles.accession==f'OSD-{n}')&(roles.split_id==0)];train=np.sort([data['meta'].index.get_loc(a) for a in s[s.role=='train'].animal]);test=np.sort([data['meta'].index.get_loc(a) for a in s[s.role=='test'].animal]);m=r.build_training(data,train);t=r.target(data,m,test);counts=pd.DataFrame(data['c'],index=data['genes'],columns=data['meta'].index);gs,yt,ye,sft,sfe,_=old.train_reference(counts,train,test)
        err=float(max(np.max(abs(m['sf']-sft)),np.max(abs(t['sf']-sfe))));rows.append(dict(check='training_only_factors_vs_independent_implementation',case=n,parameter=0,error=err,passed=err<1e-12))
        poisoned={**data,'c':data['c'].copy()};poisoned['c'][::3,test]*=7;mm=r.build_training(poisoned,train);err=float(np.max(abs(mm['logref']-m['logref'])));rows.append(dict(check='target_perturbation_training_reference_invariance',case=n,parameter=0,error=err,passed=err==0 and np.array_equal(mm['idx'],m['idx']) and np.array_equal(mm['sf'],m['sf'])))
        cache=r.TrainingCache();cached=cache.get(data,train[::-1]);direct=r.build_training(data,train);err=float(np.max(abs(cached['YY']-direct['YY'])));rows.append(dict(check='cache_sorted_training_key_algebra',case=n,parameter=0,error=err,passed=err<1e-12))
    frame=pd.DataFrame(rows);frame.to_csv(r.OUT/'independent_method_verification.csv',index=False);summary=dict(status='PASS' if frame.passed.all() else 'FAIL',n_checks=len(frame),n_failed=int((~frame.passed).sum()),max_PLS2_vs_sklearn_prediction_error=float(frame.loc[frame.check=='PLS2_covariance_deflation_vs_sklearn_NIPALS','error'].max()),seed=901531)
    (r.OUT/'independent_method_verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2));print(frame[~frame.passed].to_string(index=False));assert frame.passed.all()

if __name__=='__main__':run()
