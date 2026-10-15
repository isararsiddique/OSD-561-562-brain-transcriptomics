#!/usr/bin/env python3
"""Inductive fixed-sensor method benchmark with fully nested training estimates."""
import argparse
from collections import OrderedDict
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from run_revision import load_dataset,PANELS

RIDGE_GRID=[None,1000.,100.,10.,1.,.1,.01,.001]
COMPONENT_GRID=[0,1,2,3,4]
SCALAR_GRID=np.array([0.,.25,.5,.75,1.])
METHODS=['intercept_age_baseline','original_score','score_gain2','existing_grid_CV_scalar','risk_CV_anchored_loading','ridge_CV','PCR_CV','PLS2_CV','gain2_ridge_CV']
ANNOT=pd.read_csv(ROOT/'inputs/gene_annotation.csv.gz').set_index('ensembl_id').gene_symbol.to_dict()


def pls2_factors(X,YY,max_components=4):
    """PLS2 regression by exact covariance-maximizing sequential X/Y deflation.

    Return Q_k (sensor x training-sample) with B_scaled=Q_k @ Y.
    YY=Y @ Y.T contains training outcomes only. No target outcome enters Q.
    """
    n,s=X.shape;Xd=X.copy();F=np.eye(n);D=np.eye(s);Q=np.zeros((s,n));out={0:Q.copy()}
    for k in range(1,max_components+1):
        C=Xd.T@F
        H=C@YY@C.T
        eig,V=np.linalg.eigh((H+H.T)/2);w=V[:,-1]
        if eig[-1] <= max(1e-20,abs(eig[-1])*1e-14):
            out[k]=Q.copy();continue
        t=Xd@w;den=t@t
        if den<1e-15:out[k]=Q.copy();continue
        p=Xd.T@t/den;loading_y=t@F/den
        rotation=D@w
        Q+=np.outer(rotation,loading_y)
        Xd-=np.outer(t,p);F-=np.outer(t,loading_y)
        D=D@(np.eye(s)-np.outer(w,p))
        out[k]=Q.copy()
    return out


def gain_from_factors(Q,YY):
    """Raw-unit B=Y.T @ Q.T; exact sensor-removing projection gain."""
    H=Q@YY@Q.T;lam=max(0,float(np.linalg.eigvalsh((H+H.T)/2)[-1]));return float(np.sqrt(1+lam))


def bound_factors(Q,YY,max_gain=2.):
    """Euclidean closest B under ||B||_2 <= sqrt(L^2-1), via spectral clipping."""
    H=Q@YY@Q.T;eig,V=np.linalg.eigh((H+H.T)/2);sing=np.sqrt(np.maximum(eig,0));cap=np.sqrt(max_gain**2-1)
    scale=np.ones(len(sing));positive=sing>cap;scale[positive]=cap/sing[positive]
    return (V*scale)@V.T@Q


def ridge_factors(X,scale,alpha):
    if alpha is None:return np.zeros((X.shape[1],len(X)))
    K=X@X.T+len(X)*alpha*np.eye(len(X));return np.linalg.solve(K,X).T/scale[:,None]


def pcr_factors(X,scale,max_components=4):
    U,s,Vt=np.linalg.svd(X,full_matrices=False);out={0:np.zeros((X.shape[1],len(X)))}
    rank=int((s>max(X.shape)*np.finfo(float).eps*s[0]).sum())
    for k in range(1,max_components+1):
        kk=min(k,rank);out[k]=(Vt[:kk].T/s[:kk])@U[:,:kk].T/scale[:,None]
    return out,rank


def score_gain(norm_h,k):
    return float(np.linalg.svd(np.array([[1-k,0],[-k*norm_h,1.]]),compute_uv=False)[0])


def age_matrix(meta):
    return np.column_stack([np.ones(len(meta)),meta.age.eq('29').to_numpy(float)]) if meta.age.nunique()>1 else np.ones((len(meta),1))


def get_data(n):
    d=load_dataset(n);meta=d['meta'];mask=meta.flight.eq('GC') if n!='525' else meta.flight.isin(['HGC','VGC']);meta=meta.loc[mask].copy();c=d['counts'].loc[:,meta.index]
    meta.index=meta.animal.to_numpy()
    meta['label']=(meta.cohort.eq('ISS_T') if n!='525' else meta.flight.eq('HGC')).astype(int)
    return dict(n=n,c=c.to_numpy(float),genes=c.index.to_numpy(),meta=meta,age_full=age_matrix(meta),canonical=c.index.astype(str).str.match(r'^ENSMUSG\d+$'),source_hash=hashlib.sha256(c.to_numpy().tobytes()).hexdigest())


def build_training(data,train):
    c=data['c'];train=np.asarray(train,int);keep=data['canonical']&(c[:,train].mean(axis=1)>=10);idx=np.flatnonzero(keep)
    positive=idx[(c[idx][:,train]>0).all(axis=1)];logref=np.log(c[positive][:,train]).mean(axis=1)
    sf=np.exp(np.median(np.log(c[positive][:,train])-logref[:,None],axis=0))
    y=np.log2(c[idx][:,train]/sf+1).T;xt=data['age_full'][train];base=np.linalg.lstsq(xt,y,rcond=None)[0];res=y-xt@base
    syms=np.array([ANNOT.get(g,'') for g in data['genes'][idx]]);sensor=np.flatnonzero(np.isin(syms,PANELS['Sensor31']));outside=np.setdiff1d(np.arange(len(idx)),sensor)
    S=res[:,sensor];Y=res[:,outside];YY=Y@Y.T;YYall=YY+S@S.T;scale=np.sqrt(np.mean(S*S,axis=0));zero=int((scale<1e-10).sum());scale=np.where(scale<1e-10,1,scale);X=S/scale
    gt=data['meta'].iloc[train].label.to_numpy(float);g0=gt-xt@np.linalg.lstsq(xt,gt,rcond=None)[0];beta=S.T@g0/(g0@g0);w=beta/np.linalg.norm(beta);strain=S@w;ss=strain@strain
    gamma_norm=np.sqrt(max(1.,float(strain@YYall@strain/ss**2)));norm_h=np.sqrt(max(0,gamma_norm**2-1));cap_lambda=min(1,np.sqrt(3)/norm_h) if norm_h>0 else 1.
    ridge={str(a):ridge_factors(X,scale,a) for a in RIDGE_GRID};bounded_ridge={key:bound_factors(Q,YY) for key,Q in ridge.items()};pcr,rank=pcr_factors(X,scale);pls={k:Q/scale[:,None] for k,Q in pls2_factors(X,YY).items()}
    return dict(train=train,idx=idx,positive=positive,logref=logref,sf=sf,base=base,sensor=sensor,outside=outside,Y=Y,S=S,X=X,scale=scale,YY=YY,YYall=YYall,w=w,strain=strain,ss=ss,gamma_norm=gamma_norm,norm_h=norm_h,cap_lambda=cap_lambda,ridge=ridge,bounded_ridge=bounded_ridge,pcr=pcr,pls=pls,rank=rank,zero_variance_predictors=zero)


class TrainingCache:
    def __init__(self,size=128):self.cache=OrderedDict();self.size=size;self.hits=0;self.misses=0
    def get(self,data,train):
        key=(data['n'],tuple(sorted(map(int,train))))
        if key in self.cache:self.hits+=1;self.cache.move_to_end(key);return self.cache[key]
        self.misses+=1;m=build_training(data,key[1]);self.cache[key]=m
        if len(self.cache)>self.size:self.cache.popitem(last=False)
        return m


def target(data,m,test):
    c=data['c'];test=np.asarray(test,int);sf=[]
    for col in test:
        nz=c[m['positive'],col]>0;sf.append(np.exp(np.median(np.log(c[m['positive'][nz],col])-m['logref'][nz])))
    y=np.log2(c[m['idx']][:,test]/np.asarray(sf)+1).T;r=y-data['age_full'][test]@m['base'];S=r[:,m['sensor']];R=r[:,m['outside']];score=S@m['w']
    return dict(test=test,sf=np.asarray(sf),R=R,S=S,score=score,b0=float(np.sum(R*R)),cross=m['Y']@R.T)


def risk(m,t,M):
    # Direct SSE identity avoids materializing G outcome predictions in inner CV.
    return float(t['b0']-2*np.sum(M*t['cross'].T)+np.trace(M@m['YY']@M.T))


def tune(data,train,cache):
    errs={'ridge':np.zeros(len(RIDGE_GRID)),'bounded_ridge':np.zeros(len(RIDGE_GRID)),'PCR':np.zeros(len(COMPONENT_GRID)),'PLS2':np.zeros(len(COMPONENT_GRID)),'scalar':np.zeros(len(SCALAR_GRID))};baseline=0;cs=[];ps=[];foldrows=[]
    for left in train:
        inner=np.array([i for i in train if i!=left]);m=cache.get(data,inner);t=target(data,m,[left]);baseline+=t['b0'];M=np.outer(t['score'],m['strain']/m['ss']);C=float(np.sum(M*t['cross'].T));P=float(np.trace(M@m['YY']@M.T));cs.append(C);ps.append(P)
        for j,k in enumerate(SCALAR_GRID):errs['scalar'][j]+=risk(m,t,k*M)
        for key,fits in [('ridge',m['ridge']),('bounded_ridge',m['bounded_ridge'])]:
            for j,a in enumerate(RIDGE_GRID):errs[key][j]+=risk(m,t,t['S']@fits[str(a)])
        for key,fits in [('PCR',m['pcr']),('PLS2',m['pls'])]:
            for j,k in enumerate(COMPONENT_GRID):errs[key][j]+=risk(m,t,t['S']@fits[k])
        foldrows.append(dict(heldout_training_animal=data['meta'].iloc[left].animal,inner_train=len(inner),n_genes=len(m['idx']),n_outside=len(m['outside']),baseline_SSE=t['b0'],response_prediction_cross=C,prediction_squared_sum=P))
    choose={key:int(np.argmin(values)) for key,values in errs.items()};C=sum(cs);P=sum(ps);lam=float(np.clip(C/P,0,1)) if P>1e-20 else 0.
    omission=[float(np.clip((C-c)/(P-p),0,1)) if P-p>1e-20 else 0. for c,p in zip(cs,ps)]
    return dict(errs=errs,choose=choose,baseline=baseline,lambda_risk=lam,lambda_omission_min=min(omission),lambda_omission_max=max(omission),foldrows=foldrows,risk_C=C,risk_P=P,risk_selected_SSE=baseline-2*lam*C+lam*lam*P)


def selected_methods(m,t,cv):
    k=float(SCALAR_GRID[cv['choose']['scalar']]);lam=cv['lambda_risk'];M=np.outer(t['score'],m['strain']/m['ss'])
    out={
      'intercept_age_baseline':dict(M=np.zeros_like(M),gain=1.,parameter='baseline',sensor_score_removal=False,individual_sensor_removal=False),
      'original_score':dict(M=M,gain=m['gamma_norm'],parameter='1',sensor_score_removal=True,individual_sensor_removal=False),
      'score_gain2':dict(M=m['cap_lambda']*M,gain=np.sqrt(1+m['cap_lambda']**2*m['norm_h']**2),parameter=str(m['cap_lambda']),sensor_score_removal=True,individual_sensor_removal=False),
      'existing_grid_CV_scalar':dict(M=k*M,gain=score_gain(m['norm_h'],k),parameter=str(k),sensor_score_removal=bool(k==1),individual_sensor_removal=False),
      'risk_CV_anchored_loading':dict(M=lam*M,gain=np.sqrt(1+lam**2*m['norm_h']**2),parameter=str(lam),sensor_score_removal=True,individual_sensor_removal=False)}
    for method,key,grid,candidates in [('ridge_CV','ridge',RIDGE_GRID,m['ridge']),('gain2_ridge_CV','bounded_ridge',RIDGE_GRID,m['bounded_ridge']),('PCR_CV','PCR',COMPONENT_GRID,m['pcr']),('PLS2_CV','PLS2',COMPONENT_GRID,m['pls'])]:
        value=grid[cv['choose'][key]];Q=candidates[str(value)] if key in ['ridge','bounded_ridge'] else candidates[value]
        out[method]=dict(M=t['S']@Q,gain=gain_from_factors(Q,m['YY']),parameter='baseline' if value is None else str(value),sensor_score_removal=True,individual_sensor_removal=True,Q=Q)
    return out


def summarize(frame):
    rows=[]
    for (n,method),part in frame.groupby(['accession','method']):
        row=dict(accession=n,method=method,n_overlapping_splits=len(part),uncertainty='overlapping-split robustness range')
        for col in ['outside_prediction_MSE_ratio','operator_gain','heldout_outside_contrast_energy_removed_pct','fraction_outside_genes_prediction_worse']:
            values=part[col].to_numpy();row[col+'_mean']=float(values.mean());row[col+'_median']=float(np.median(values));row[col+'_q025']=float(np.quantile(values,.025));row[col+'_q975']=float(np.quantile(values,.975))
        row['fraction_splits_prediction_improved']=float(part.outside_prediction_MSE_ratio.lt(1-1e-12).mean());rows.append(row)
    return pd.DataFrame(rows)


def main(limit=None):
    start=time.time();roles=pd.read_csv(ROOT/'extensions/ground_holdout/split_animal_roles.csv');old=pd.read_csv(ROOT/'extensions/ground_holdout/holdout_split_results.csv');cache=TrainingCache();results=[];cvrows=[];foldrows=[];animalrows=[];checks=[];factors_saved=[]
    for n in ['561','562','525']:
        data=get_data(n);allroles=roles[roles.accession==f'OSD-{n}'];splitids=sorted(allroles.split_id.unique());splitids=splitids if limit is None else splitids[:limit]
        for sid in splitids:
            r=allroles[allroles.split_id==sid];train=np.sort([data['meta'].index.get_loc(a) for a in r[r.role=='train'].animal]);test=np.sort([data['meta'].index.get_loc(a) for a in r[r.role=='test'].animal]);m=cache.get(data,train);t=target(data,m,test);cv=tune(data,train,cache);methods=selected_methods(m,t,cv)
            g=data['meta'].iloc[test].label.to_numpy(float);xe=data['age_full'][test];gres=g-xe@np.linalg.lstsq(xe,g,rcond=None)[0];rawbeta=gres@t['R']/(gres@gres)
            for method in METHODS:
                f=methods[method];P=f['M']@m['Y'];E=t['R']-P;sse=float(np.sum(E*E));ratio=sse/t['b0'];ratio=1. if abs(ratio-1)<1e-12 else ratio;b=gres@E/(gres@gres)
                results.append(dict(accession=f'OSD-{n}',split_id=int(sid),method=method,n_training=len(train),n_test=len(test),n_genes=len(m['idx']),n_sensors=len(m['sensor']),n_outside=len(m['outside']),outside_baseline_MSE=t['b0']/t['R'].size,outside_prediction_MSE=sse/t['R'].size,outside_prediction_MSE_ratio=ratio,operator_gain=f['gain'],selected_parameter=f['parameter'],weighted_sensor_score_removed=f['sensor_score_removal'],all_individual_sensors_removed=f['individual_sensor_removal'],heldout_outside_contrast_energy_removed_pct=100*(1-b@b/(rawbeta@rawbeta)),fraction_outside_genes_prediction_worse=float((np.mean(E*E,axis=0)>np.mean(t['R']**2,axis=0)).mean()),inner_risk_lambda=cv['lambda_risk'],inner_risk_lambda_fold_omission_min=cv['lambda_omission_min'],inner_risk_lambda_fold_omission_max=cv['lambda_omission_max'],zero_variance_predictors=m['zero_variance_predictors'],training_sensor_predictor_rank=m['rank']))
                checks.append(dict(check='SSE_gram_identity',accession=f'OSD-{n}',split_id=int(sid),method=method,error=abs(sse-risk(m,t,f['M']))/max(1,sse),passed=abs(sse-risk(m,t,f['M']))/max(1,sse)<1e-10))
                for j,i in enumerate(test):animalrows.append(dict(accession=f'OSD-{n}',split_id=int(sid),method=method,animal=data['meta'].iloc[i].animal,baseline_MSE=float(np.mean(t['R'][j]**2)),prediction_MSE=float(np.mean(E[j]**2))))
                if method in ['original_score','score_gain2','existing_grid_CV_scalar']:
                    label={'original_score':'score','score_gain2':'bounded_score','existing_grid_CV_scalar':'training_cv_shrinkage'}[method];ref=old[(old.accession==f'OSD-{n}')&(old.split_id==sid)&(old.method==label)].iloc[0]
                    err=abs(ratio-ref.outside_prediction_MSE_ratio_to_nuisance_only);checks.append(dict(check='existing_result_reproduction',accession=f'OSD-{n}',split_id=int(sid),method=method,error=err,passed=err<1e-10))
            for key,values in cv['errs'].items():
                grid=SCALAR_GRID if key=='scalar' else (RIDGE_GRID if key in ['ridge','bounded_ridge'] else COMPONENT_GRID)
                for j,(param,err) in enumerate(zip(grid,values)):cvrows.append(dict(accession=f'OSD-{n}',split_id=int(sid),algorithm=key,parameter='baseline' if param is None else str(param),inner_CV_SSE=float(err),inner_CV_SSE_ratio=float(err/cv['baseline']),selected=j==cv['choose'][key],inner_target_count=len(train)))
            cvrows.append(dict(accession=f'OSD-{n}',split_id=int(sid),algorithm='continuous_risk_CV',parameter=str(cv['lambda_risk']),inner_CV_SSE=cv['risk_selected_SSE'],inner_CV_SSE_ratio=cv['risk_selected_SSE']/cv['baseline'],selected=True,inner_target_count=len(train)))
            for row in cv['foldrows']:foldrows.append(dict(accession=f'OSD-{n}',split_id=int(sid),**row))
            checks.append(dict(check='continuous_CV_not_worse_than_grid_inner_SSE',accession=f'OSD-{n}',split_id=int(sid),method='risk_CV_anchored_loading',error=max(0,cv['risk_selected_SSE']-cv['errs']['scalar'].min())/max(1,cv['baseline']),passed=cv['risk_selected_SSE']<=cv['errs']['scalar'].min()+1e-8*cv['baseline']))
            if sid==0:
                # Training-only factorization snapshots let independent reviewers certify dense operators.
                chosen={method:methods[method]['Q'] for method in ['ridge_CV','gain2_ridge_CV','PCR_CV','PLS2_CV']}
                np.savez_compressed(OUT/f'OSD-{n}_representative_factorization.npz',Y=m['Y'],S=m['S'],X=m['X'],scale=m['scale'],w=m['w'],score_train=m['strain'],**chosen);factors_saved.append(f'OSD-{n}_representative_factorization.npz')
            if sid%20==0:
                part=pd.DataFrame(results);table=part[part.accession==f'OSD-{n}'].groupby('method').outside_prediction_MSE_ratio.mean();print(f'OSD-{n} split{sid+1}/{len(splitids)} elapsed{time.time()-start:.1f}s\n{table.to_string()}',flush=True)
        # Every dataset checkpoint retains all results already computed.
        pd.DataFrame(results).to_csv(OUT/'outer_split_results.csv',index=False);pd.DataFrame(cvrows).to_csv(OUT/'nested_CV_candidates.csv',index=False);pd.DataFrame(foldrows).to_csv(OUT/'nested_CV_risk_fold_components.csv',index=False);pd.DataFrame(animalrows).to_csv(OUT/'outer_animal_prediction_errors.csv',index=False);pd.DataFrame(checks).to_csv(OUT/'verification_checks.csv',index=False)
        print('Completed',n,'cachehits',cache.hits,'misses',cache.misses,flush=True)
    frame=pd.DataFrame(results);summary=summarize(frame);summary.to_csv(OUT/'method_summary.csv',index=False)
    animal=pd.DataFrame(animalrows).groupby(['accession','method','animal'],as_index=False)[['baseline_MSE','prediction_MSE']].mean();animal['ratio_of_animal_average_errors']=animal.prediction_MSE/animal.baseline_MSE;animal.to_csv(OUT/'animal_average_prediction_summary.csv',index=False)
    check=pd.DataFrame(checks);status=dict(status='PASS' if check.passed.all() else 'FAIL',n_checks=len(check),n_failed=int((~check.passed).sum()),n_splits=frame[['accession','split_id']].drop_duplicates().shape[0],n_method_rows=len(frame),elapsed_seconds=time.time()-start,cache_hits=cache.hits,cache_misses=cache.misses,protocol_sha256=hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest(),factorization_snapshots=factors_saved,source_root=str(ROOT),max_existing_result_discrepancy=float(check.loc[check.check=='existing_result_reproduction','error'].max()),exploratory=True)
    (OUT/'run_summary.json').write_text(json.dumps(status,indent=2)+'\n');print(summary[['accession','method','outside_prediction_MSE_ratio_mean','operator_gain_mean']].to_string(index=False),flush=True);print(json.dumps(status,indent=2),flush=True)
    if not check.passed.all():raise AssertionError(check[~check.passed].to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int);main(parser.parse_args().limit)
