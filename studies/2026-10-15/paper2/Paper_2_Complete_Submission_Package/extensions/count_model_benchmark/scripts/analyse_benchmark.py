"""Exact animal-label sensitivity tests on standard count-likelihood scores.

This is NOT a LeafCutter/rMATS implementation or a DM likelihood-ratio test.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import json,hashlib,itertools,math,time,argparse
import numpy as np
import pandas as pd
from scipy.special import digamma,gammaln
from scipy.optimize import minimize
from scipy.stats import spearmanr
BASE=Path(__file__).resolve().parents[1];PACKAGE=BASE.parents[1]
CONTRASTS=[('FLT','GC'),('FLT','VIV'),('GC','VIV')]
TOL=1e-10

def adjust(p,method):
    p=np.asarray(p);n=len(p);o=np.argsort(p,kind='stable');s=p[o]
    if method=='holm':q=np.maximum.accumulate(s*(n-np.arange(n)))
    else:q=np.minimum.accumulate((s*n/(np.arange(n)+1)*(np.sum(1/np.arange(1,n+1)) if method=='by' else 1))[::-1])[::-1]
    ans=np.empty(n);ans[o]=np.minimum(1,q);return ans

def dm_objective(loga,x):
    alpha=np.exp(loga);A=alpha.sum();N=x.sum(axis=1)
    ll=np.sum(gammaln(A)-gammaln(N+A))+np.sum(gammaln(x+alpha)-gammaln(alpha))
    grad=alpha*(np.sum(digamma(A)-digamma(N+A))+np.sum(digamma(x+alpha)-digamma(alpha),axis=0))
    return -ll,-grad

def fit_dm(x):
    pi=x.sum(axis=0)/x.sum();fits=[]
    for initial_A in (10.,1000.):
        r=minimize(dm_objective,np.log(np.maximum(pi*initial_A,1e-6)),args=(x,),jac=True,method='L-BFGS-B',bounds=[(np.log(1e-6),np.log(1e5))]*x.shape[1],options=dict(ftol=1e-11,gtol=1e-6,maxiter=1000,maxls=40))
        fits.append((r,initial_A))
    good=[(r,st) for r,st in fits if r.success and np.isfinite(r.fun) and np.isfinite(r.x).all()]
    r,st=min(good or fits,key=lambda rr:rr[0].fun)
    alpha=np.exp(r.x);A=alpha.sum();q=digamma(x+alpha)-digamma(alpha)
    S=alpha*(q-q@(alpha/A)[:,None])
    # q@(alpha/A)[:,None] yields one common row mean.
    diag=dict(fit_valid=bool(good),concentration=A,rho=1/(A+1),selected_start=st,log_likelihood_without_multinomial_constant=-float(r.fun),n_iterations=int(r.nit),gradient_max_abs=float(np.abs(r.jac).max()),alpha_bound_hit=bool(((alpha<=1e-6*(1+1e-5))|(alpha>=1e5*(1-1e-5))).any()),start10_converged=bool(fits[0][0].success),start1000_converged=bool(fits[1][0].success),start10_loglik=-float(fits[0][0].fun),start1000_loglik=-float(fits[1][0].fun),best_solution_loglik_gap=abs(float(fits[0][0].fun-fits[1][0].fun)),message=str(r.message))
    return S,alpha,diag

def designs(groups):
    out=[]
    for a,b in CONTRASTS:
        ix=np.flatnonzero(np.isin(groups,[a,b]));nA=np.count_nonzero(groups[ix]==a);nB=len(ix)-nA
        comb=np.array(list(itertools.combinations(range(len(ix)),nA)),dtype=np.int16)
        D=np.full((len(comb),len(ix)),-1/nB);D[np.arange(len(comb))[:,None],comb]=1/nA
        d=np.where(groups[ix]==a,1/nA,-1/nB);out.append((a+'_vs_'+b,ix,D,d))
    return out

def whiten(s):
    s=s-s.mean(axis=0);C=s.T@s/(len(s)-1);ev,U=np.linalg.eigh(C);keep=ev>max(0,ev.max())*1e-10
    w=s@U[:,keep]/np.sqrt(ev[keep]);rank=int(keep.sum());condition=float(ev[keep].max()/ev[keep].min()) if rank else np.inf
    return w,rank,condition

def family(scores,valid,design,ids,model,rep):
    rows=[];globalrows=[];nulllist=[]
    for contrast,ix,D,d in design:
        ws=[];off=[0];ranks=[];conditions=[]
        for e,s in enumerate(scores):
            w,rank,cond=whiten(s[ix]);ws.append(w);off.append(off[-1]+rank);ranks.append(rank);conditions.append(cond)
        ranks=np.array(ranks);off=np.array(off);h=np.concatenate(ws,axis=1);obs=np.add.reduceat((d@h)**2,off[:-1]);ge=np.zeros(len(ids),dtype=np.int64)
        valid2=valid&(ranks>0);sobs=float(np.mean(obs[valid2]/ranks[valid2]));gnull=np.zeros(len(D))
        for first in range(0,len(D),1024):
            proj=D[first:first+1024]@h;stats=np.add.reduceat(proj*proj,off[:-1],axis=1)
            ge+=(stats>=obs[None,:]-TOL*np.maximum(1,obs)[None,:]).sum(axis=0)
            gnull[first:first+1024]=(stats[:,valid2]/ranks[None,valid2]).mean(axis=1)
        ng=int(np.count_nonzero(gnull>=sobs-TOL*max(1,sobs)))
        for e,cid in enumerate(ids):rows.append(dict(replicate=rep,model=model,contrast=contrast,cluster_id=cid,valid=bool(valid2[e]),rank=int(ranks[e]),covariance_condition=conditions[e],statistic=obs[e],n_ge_observed=int(ge[e]) if valid2[e] else len(D),n_assignments=len(D),exact_p=ge[e]/len(D) if valid2[e] else 1.))
        globalrows.append(dict(replicate=rep,model=model,contrast=contrast,n_valid_clusters=int(valid2.sum()),statistic=sobs,n_ge_observed=ng,n_assignments=len(D),exact_p=ng/len(D)))
        if rep=='full':nulllist.append(pd.DataFrame(dict(model=model,contrast=contrast,assignment_index=np.arange(len(D)),statistic=gnull)))
    local=pd.DataFrame(rows);local['BH_q_full_2769']=adjust(local.exact_p,'bh');local['BY_q_full_2769']=adjust(local.exact_p,'by')
    glob=pd.DataFrame(globalrows);glob['Holm_p_3']=adjust(glob.exact_p,'holm')
    return local,glob,nulllist

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--skip-thinning',action='store_true');args=ap.parse_args()
    t=time.monotonic();out=BASE/'results';out.mkdir(exist_ok=True)
    z=np.load(PACKAGE/'extensions/multiway_splicing/tables/frozen_features.npz');samples=z['samples'];ids=z['cluster_ids'];off=z['offsets'];E=len(ids)
    raw=pd.read_csv(PACKAGE/'extensions/multiway_transport/tables/RR23_all_frozen_alternative_counts.csv.gz');x=raw[list(samples)].to_numpy(dtype=float).T
    assert E==923 and x.shape==(26,3376) and (x>=0).all() and np.array_equal(np.repeat(ids,np.diff(off)),raw.cluster_id)
    P=x/np.repeat(np.add.reduceat(x,off[:-1],axis=1),np.diff(off),axis=1);assert np.max(np.abs(np.sqrt(P)-z['H']))<1e-13
    groups=pd.read_csv(PACKAGE/'inputs/animal_metadata.csv').set_index('sample').loc[samples,'group'].to_numpy();design=designs(groups)
    defs=pd.read_csv(PACKAGE/'extensions/multiway_splicing/tables/selected_cluster_definitions.csv').set_index('cluster_id').loc[ids]
    original=pd.read_csv(PACKAGE/'extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz');diags=[];alpharows=[];dm=[];mn=[];valid=[]
    for e,cid in enumerate(ids):
        xe=x[:,off[e]:off[e+1]];s,alpha,diag=fit_dm(xe);dm.append(s);valid.append(diag['fit_valid']);pi=xe.sum(axis=0)/xe.sum();mn.append(xe-xe.sum(axis=1)[:,None]*pi)
        diag.update(replicate='full',cluster_id=cid,n_alternatives=xe.shape[1],min_selected_reads=int(xe.sum(axis=1).min()),median_selected_reads=float(np.median(xe.sum(axis=1))),max_selected_reads=int(xe.sum(axis=1).max()));diags.append(diag)
        for j,a in enumerate(alpha):alpharows.append(dict(replicate='full',cluster_id=cid,alternative_index=j,alpha=a,model_mean_pi=a/alpha.sum()))
        if (e+1)%100==0:print('full fits',e+1,'elapsed',round(time.monotonic()-t,1),flush=True)
    np.savez_compressed(out/'full_fitted_score_vectors.npz',dm=np.concatenate(dm,axis=1),multinomial=np.concatenate(mn,axis=1),offsets=off,cluster_ids=ids,samples=samples)
    allglobals=[];alllocals=[];nulls=[]
    for name,scores,v in [('DM_score',dm,np.array(valid)),('multinomial_score',mn,np.ones(E,dtype=bool))]:
        loc,glob,nul=family(scores,v,design,ids,name,'full');alllocals.append(loc);allglobals.append(glob);nulls+=nul
        print(name,loc.groupby('contrast').apply(lambda d:int((d.BY_q_full_2769<=.05).sum()),include_groups=False).to_dict(),flush=True)
    complete=pd.concat(alllocals,ignore_index=True).merge(original[['cluster_id','contrast','exact_p','BY_q_full_3E','BH_q_full_3E','native_mean_total_variation','max_abs_proportion_difference']],on=['cluster_id','contrast'],validate='many_to_one',suffixes=('','_original'))
    complete.to_csv(out/'full_model_comparison_complete.csv.gz',index=False);pd.concat(nulls,ignore_index=True).to_csv(out/'full_global_exact_null.csv.gz',index=False)
    # Stricter gate keeps original family; no refitting or outcome-driven selection.
    strict=alllocals[0].copy();gate=dict(zip(ids,defs.min_dominance>=.95-1e-12));strict['dominance95_valid']=strict.cluster_id.map(gate);strict.loc[~strict.dominance95_valid,'exact_p']=1
    strict['BH_q_full_2769']=adjust(strict.exact_p,'bh');strict['BY_q_full_2769']=adjust(strict.exact_p,'by');strict.to_csv(out/'DM_dominance95_complete.csv.gz',index=False)
    thinrows=[]
    if not args.skip_thinning:
        idx=np.load(PACKAGE/'extensions/multiway_robustness/results/frozen_count_index.npz');gates=pd.read_csv(PACKAGE/'extensions/multiway_robustness/results/thinning_measurement_gates.csv.gz')
        for rep in range(1,11):
            thin=np.load(PACKAGE/f'extensions/multiway_robustness/thinned_counts/replicate_{rep:02d}_seed_{2026101500+rep}.npz')['counts'];xt=thin[idx['selected_rows']].T.astype(float);del thin
            gate=gates[gates.replicate==rep].set_index('cluster_id').loc[ids,'measurement_valid'].to_numpy(bool);scores=[];v=[]
            for e,cid in enumerate(ids):
                xe=xt[:,off[e]:off[e+1]];s,alpha,diag=fit_dm(xe);scores.append(s);v.append(diag['fit_valid']);diag.update(replicate=str(rep),cluster_id=cid,n_alternatives=xe.shape[1],min_selected_reads=int(xe.sum(axis=1).min()),median_selected_reads=float(np.median(xe.sum(axis=1))),max_selected_reads=int(xe.sum(axis=1).max()));diags.append(diag)
                for j,a in enumerate(alpha):alpharows.append(dict(replicate=str(rep),cluster_id=cid,alternative_index=j,alpha=a,model_mean_pi=a/alpha.sum()))
            loc,glob,_=family(scores,np.array(v)&gate,design,ids,'DM_score',str(rep));thinrows.append(loc);allglobals.append(glob)
            print('thinning',rep,'fits_valid',sum(v),'measurement_valid',sum(gate),'BY',sum(loc.BY_q_full_2769<=.05),'elapsed',round(time.monotonic()-t,1),flush=True)
            pd.concat(thinrows,ignore_index=True).to_csv(out/'DM_thinning_complete_frozen_families.csv.gz',index=False)
            pd.DataFrame(diags).to_csv(out/'all_model_fit_diagnostics.csv.gz',index=False)
    pd.DataFrame(diags).to_csv(out/'all_model_fit_diagnostics.csv.gz',index=False);pd.DataFrame(alpharows).to_csv(out/'all_fitted_alpha.csv.gz',index=False)
    pd.concat(allglobals,ignore_index=True).to_csv(out/'all_global_results.csv',index=False)
    hitrows=[]
    for (mod,c),d in complete.groupby(['model','contrast']):
        old=d.BY_q_full_3E<=.05;new=d.BY_q_full_2769<=.05
        hitrows.append(dict(model=mod,contrast=c,n_tests=len(d),n_valid=int(d.valid.sum()),n_original_BY=int(old.sum()),n_model_BY=int(new.sum()),n_shared_BY=int((new&old).sum()),spearman_negative_log_p=float(spearmanr(-np.log10(d.exact_p),-np.log10(d.exact_p_original)).statistic),median_covariance_condition=float(d.covariance_condition.median())))
    pd.DataFrame(hitrows).to_csv(out/'model_concordance_summary.csv',index=False)
    tcaf='6:2:42663893:start';tc=complete[complete.cluster_id==tcaf];tc.to_csv(out/'Tcaf1_full_model_comparison.csv',index=False)
    if thinrows:
        thin=pd.concat(thinrows,ignore_index=True);thin[thin.cluster_id==tcaf].to_csv(out/'Tcaf1_all10_thinning_comparison.csv',index=False)
        orig=original[['cluster_id','contrast','BY_q_full_3E']];track=thin.merge(orig,on=['cluster_id','contrast'],validate='many_to_one');track[track.BY_q_full_3E<=.05].to_csv(out/'original43_DM_thinning_tracking.csv',index=False)
        s=thin.groupby(['replicate','contrast']).agg(n_valid=('valid','sum'),n_BY=('BY_q_full_2769',lambda s:int((s<=.05).sum())),n_BH=('BH_q_full_2769',lambda s:int((s<=.05).sum()))).reset_index();s.to_csv(out/'thinning_family_summary.csv',index=False)
    summary=dict(n_clusters=E,n_animals=len(samples),n_alternatives=x.shape[1],n_tests_per_family=3*E,full_DM_BY=int((alllocals[0].BY_q_full_2769<=.05).sum()),full_multinomial_BY=int((alllocals[1].BY_q_full_2769<=.05).sum()),full_DM_fit_failures=int(E-sum(valid)),strict95_DM_BY=int((strict.BY_q_full_2769<=.05).sum()),thinning_replicates=len(thinrows),elapsed_seconds=time.monotonic()-t,plan_sha256=hashlib.sha256((BASE/'RETROSPECTIVE_ANALYSIS_PLAN.md').read_bytes()).hexdigest(),fit_mean_is_pooled=True,thinning_scope='ten original A-half draws; B complements not analyzed in this benchmark',exact_calibration='fixed-size animal-label enumeration; conditional on exchangeability',scientific_status='retrospective model sensitivity; not independent biological replication or a causal adjustment')
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
