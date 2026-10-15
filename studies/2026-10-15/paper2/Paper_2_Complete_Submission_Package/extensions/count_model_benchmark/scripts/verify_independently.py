"""Independent checks; no import from the analysis implementation."""
from pathlib import Path
import json,hashlib,itertools,math,gzip
import numpy as np
import pandas as pd
from scipy.special import gammaln,softmax
from scipy.stats import dirichlet_multinomial,rankdata
BASE=Path(__file__).resolve().parents[1];P=BASE.parents[1];R=BASE/'results'
checks=[]
def check(name,ok,detail=''):
    checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))

def by_independent(p):
    n=len(p);harmonic=sum(1/k for k in range(1,n+1));order=sorted(range(n),key=lambda i:p[i]);ans=np.ones(n);running=1.
    for rank in range(n,0,-1):
        i=order[rank-1];running=min(running,p[i]*n*harmonic/rank);ans[i]=running
    return ans

def scalar_ll(x,alpha):
    A=sum(alpha);return sum(gammaln(A)-gammaln(sum(row)+A)+sum(gammaln(v+a)-gammaln(a) for v,a in zip(row,alpha)) for row in x)

def main():
    manifests=json.loads((BASE/'inputs/frozen_manifest.json').read_text())
    payload=json.loads((BASE/'inputs/frozen_payload_sha256.json').read_text())
    for r in manifests:
        raw=(P/r['path']).read_bytes();same_bytes=hashlib.sha256(raw).hexdigest()==r['sha256']
        same_payload=False
        if r['path'] in payload:
            frozen=payload[r['path']]
            same_payload=(frozen['original_compressed_sha256']==r['sha256'] and hashlib.sha256(gzip.decompress(raw)).hexdigest()==frozen['decompressed_sha256'])
        check('input SHA256 or identical frozen gzip payload '+r['path'],same_bytes or same_payload,'byte identical' if same_bytes else 'decompressed bytes identical; gzip header differs')
    z=np.load(P/'extensions/multiway_splicing/tables/frozen_features.npz');ids=z['cluster_ids'];off=z['offsets'];samples=z['samples']
    scores=np.load(R/'full_fitted_score_vectors.npz');counts=pd.read_csv(P/'extensions/multiway_transport/tables/RR23_all_frozen_alternative_counts.csv.gz');X=counts[list(samples)].to_numpy(float).T
    groups=pd.read_csv(P/'inputs/animal_metadata.csv').set_index('sample').loc[samples,'group'].to_numpy()
    fits=pd.read_csv(R/'all_model_fit_diagnostics.csv.gz');alpha=pd.read_csv(R/'all_fitted_alpha.csv.gz');local=pd.read_csv(R/'full_model_comparison_complete.csv.gz');glob=pd.read_csv(R/'all_global_results.csv');null=pd.read_csv(R/'full_global_exact_null.csv.gz')
    af=alpha[alpha.replicate=='full'];ff=fits[fits.replicate=='full'];tfits=fits[fits.replicate!='full']
    check('923 full and 9230 thinning fitted clusters',len(ff)==923 and len(tfits)==9230)
    check('Every alpha is finite and within prespecified bounds',np.isfinite(alpha.alpha).all() and alpha.alpha.between(1e-6*(1-1e-9),1e5*(1+1e-9)).all())
    check('all 10153 fits converged',fits.fit_valid.all())
    check('Complete 2769-test families for both models',local.groupby('model').size().eq(2769).all())
    check('ranks equal K minus 1',all(int(r.rank)==off[np.flatnonzero(ids==r.cluster_id)[0]+1]-off[np.flatnonzero(ids==r.cluster_id)[0]]-1 for r in local.itertuples()))
    check('counts reconstruct frozen composition',np.allclose(X/np.repeat(np.add.reduceat(X,off[:-1],axis=1),np.diff(off),axis=1),z['H']**2,atol=1e-14))
    # Every original native effect is unchanged; no alternative effect was substituted.
    original=pd.read_csv(P/'extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz').set_index(['cluster_id','contrast'])
    check('all native composition effects unchanged',all(abs(r.native_mean_total_variation-original.loc[(r.cluster_id,r.contrast),'native_mean_total_variation'])<1e-13 for r in local.itertuples()))
    for model,df in local.groupby('model'):
        check(model+' independent BY adjustment',np.allclose(by_independent(df.exact_p.to_numpy()),df.BY_q_full_2769,atol=1e-14))
        check(model+' exact fractions',np.allclose(df.exact_p,df.n_ge_observed/df.n_assignments,atol=1e-14))
    selected=sorted(set([0,17,100,461,922,int(np.flatnonzero(ids=='6:2:42663893:start')[0])]))
    verification_stats=[]
    for e in selected:
        cid=ids[e];x=X[:,off[e]:off[e+1]];a=af[af.cluster_id==cid].sort_values('alternative_index').alpha.to_numpy();S=scores['dm'][:,off[e]:off[e+1]]
        ll=scalar_ll(x,a);lib=dirichlet_multinomial.logpmf(x,a,x.sum(axis=1)).sum();comb=np.sum(gammaln(x.sum(axis=1)+1)-gammaln(x+1).sum(axis=1))
        check(cid+' scalar likelihood versus SciPy distribution',abs(ll+comb-lib)<1e-6,abs(ll+comb-lib))
        check(cid+' stored likelihood',abs(ll-ff.set_index('cluster_id').loc[cid,'log_likelihood_without_multinomial_constant'])<1e-6)
        # Perturb each composition logit while holding concentration fixed.
        A=sum(a);eta=np.log(a/A);step=1e-5;num=np.zeros_like(S)
        for j in range(len(a)):
            plus=eta.copy();minus=eta.copy();plus[j]+=step;minus[j]-=step
            ap=A*softmax(plus);am=A*softmax(minus)
            num[:,j]=(dirichlet_multinomial.logpmf(x,ap,x.sum(axis=1))-dirichlet_multinomial.logpmf(x,am,x.sum(axis=1)))/(2*step)
        check(cid+' logit score finite differences',np.allclose(S,num,rtol=3e-4,atol=3e-4),float(abs(S-num).max()))
        check(cid+' redundant score dimension sums to zero',np.max(abs(S.sum(axis=1)))<1e-10)
        for groupA,groupB in [('FLT','GC'),('FLT','VIV'),('GC','VIV')]:
            c=groupA+'_vs_'+groupB;ix=np.flatnonzero(np.isin(groups,[groupA,groupB]));n=len(ix);nA=np.count_nonzero(groups[ix]==groupA);nB=n-nA
            for model,key in [('DM_score','dm'),('multinomial_score','multinomial')]:
                s=scores[key][ix,off[e]:off[e+1]];C=np.cov(s,rowvar=False,ddof=1);inverse=np.linalg.pinv(C,rcond=1e-10,hermitian=True)
                delta=s[groups[ix]==groupA].mean(axis=0)-s[groups[ix]==groupB].mean(axis=0);obs=float(delta@inverse@delta)
                r=local[(local.cluster_id==cid)&(local.contrast==c)&(local.model==model)].iloc[0]
                check(cid+' '+c+' '+model+' direct covariance statistic',abs(obs-r.statistic)<1e-8)
                combination=np.array(list(itertools.combinations(range(n),nA)));SA=s[combination].sum(axis=1);difference=SA/nA-(s.sum(axis=0)[None,:]-SA)/nB;stats=np.einsum('ni,ij,nj->n',difference,inverse,difference)
                nge=int(np.count_nonzero(stats>=obs-1e-10*max(1,obs)));check(cid+' '+c+' '+model+' direct exact exceedance count',nge==r.n_ge_observed,(nge,int(r.n_ge_observed)))
                # Exhaustive reference labels must give a super-uniform p distribution,
                # using all ties inclusively (no randomization within ties).
                p=rankdata(-np.round(stats,10),method='max')/len(stats)
                rates=[float(np.mean(p<=level+1e-15)) for level in [.01,.05,.1]]
                check(cid+' '+c+' '+model+' exhaustive null calibration',all(rate<=level+1e-12 for rate,level in zip(rates,[.01,.05,.1])),rates)
                verification_stats.append(dict(cluster_id=cid,contrast=c,model=model,direct_statistic=obs,direct_n_ge=nge,n_assignments=len(stats),reference_type1_at01=rates[0],reference_type1_at05=rates[1],reference_type1_at10=rates[2]))
    # Independently check every global p against saved full null arrays.
    for r in glob[glob.replicate=='full'].itertuples():
        vv=null[(null.model==r.model)&(null.contrast==r.contrast)].statistic.to_numpy();nge=np.count_nonzero(vv>=r.statistic-1e-10*max(1,r.statistic));check('global null count '+r.model+' '+r.contrast,nge==r.n_ge_observed)
    thin=pd.read_csv(R/'DM_thinning_complete_frozen_families.csv.gz');gates=pd.read_csv(P/'extensions/multiway_robustness/results/thinning_measurement_gates.csv.gz')
    check('10 complete frozen thinning families',len(thin)==27690 and thin.groupby('replicate').size().eq(2769).all())
    for rep,d in thin.groupby('replicate'):
        check('thinning '+str(rep)+' independent BY',np.allclose(by_independent(d.exact_p.to_numpy()),d.BY_q_full_2769,atol=1e-13))
        gate=gates[gates.replicate==rep].set_index('cluster_id').measurement_valid
        bad=~d.cluster_id.map(gate).to_numpy();check('thinning '+str(rep)+' failed gates remain p1',np.all(d.exact_p.to_numpy()[bad]==1))
    strict=pd.read_csv(R/'DM_dominance95_complete.csv.gz');check('95% gate retains 2769 denominator',len(strict)==2769 and np.all(strict.loc[~strict.dominance95_valid,'exact_p']==1));check('95% gate independent BY',np.allclose(by_independent(strict.exact_p.to_numpy()),strict.BY_q_full_2769,atol=1e-13))
    pd.DataFrame(verification_stats).to_csv(R/'independent_selected_exact_checks.csv',index=False);pd.DataFrame(checks).to_csv(R/'independent_verification.csv',index=False)
    summary=dict(n_checks=len(checks),n_passed=sum(r['passed'] for r in checks),all_passed=all(r['passed'] for r in checks),scope='Independent scalar likelihood, SciPy distribution, finite difference scores, direct covariance inverse and complete selected permutation counts; all output family corrections and input hashes.')
    (R/'independent_verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));assert summary['all_passed'],pd.DataFrame(checks).query('passed == False').to_string(index=False)
if __name__=='__main__':main()
