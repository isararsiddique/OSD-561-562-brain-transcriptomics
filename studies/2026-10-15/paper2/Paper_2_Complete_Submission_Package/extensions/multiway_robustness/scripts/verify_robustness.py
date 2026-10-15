"""Independent count reconstruction, closed-form exact tests, SciPy correction."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import argparse,itertools,json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
BASE=Path(__file__).resolve().parents[1];PACKAGE=BASE.parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=PACKAGE);args=p.parse_args();old=args.base/'extensions/multiway_splicing';R=BASE/'results'
    checks=[]
    def check(name,ok,detail=''):
        assert bool(ok),name+' '+str(detail)
        checks.append(dict(check=name,passed=True,detail=str(detail)))
    df=pd.read_csv(args.base/'inputs/all_animal_junction_counts.csv.gz',index_col=0);x=df.to_numpy(np.int64);z=np.load(R/'frozen_count_index.npz');sr=z['selected_rows'];off=z['offsets'];br=z['boundary_rows'];bo=z['boundary_offsets'];ids=z['cluster_ids'];samples=z['samples'];E=len(ids)
    check('original count dimensions and frozen family',x.shape==(805562,26) and E==923 and np.array_equal(samples,df.columns))
    orig=np.load(old/'tables/frozen_features.npz');H=orig['H'];meta=pd.read_csv(args.base/'inputs/animal_metadata.csv').set_index('sample').loc[samples];groups=meta.group.to_numpy();alts=pd.read_csv(old/'tables/animal_alternative_compositions.csv.gz')[['cluster_id','alternative_index','junction_id']].drop_duplicates().set_index(['cluster_id','alternative_index'])
    for e,cid in enumerate(ids):
        actual=list(df.index[sr[off[e]:off[e+1]]]);expected=[alts.loc[(cid,j),'junction_id'] for j in range(off[e+1]-off[e])]
        assert actual==expected,cid
    check('all frozen alternatives independently matched to original rows',True,f'{len(sr)} alternative occurrences')
    parts=df.index.str.split(':',expand=True)
    for boundary,j in [('start',2),('end',3)]:
        key=np.array([str(a)+':'+str(b)+':'+str(c)+':'+boundary for a,b,c in zip(parts.get_level_values(0),parts.get_level_values(1),parts.get_level_values(j))]);need=set(c for c in ids if c.endswith(':'+boundary));mask=np.array([k in need for k in key]);independent=pd.Series(np.flatnonzero(mask)).groupby(key[mask],sort=False).apply(list)
        for e,cid in enumerate(ids):
            if cid in need:assert np.array_equal(np.array(independent[cid]),br[bo[e]:bo[e+1]]),cid
    check('all canonical boundary denominators independently complete',True,f'{len(br)} canonical junction occurrences')
    spaces=[]
    for a,b in [('FLT','GC'),('FLT','VIV'),('GC','VIV')]:
        ix=np.flatnonzero(np.isin(groups,[a,b]));nA=np.count_nonzero(groups[ix]==a);nB=len(ix)-nA;comb=np.array(list(itertools.combinations(range(len(ix)),nA)));d=np.where(groups[ix]==a,1/nA,-1/nB);spaces.append((a+'_vs_'+b,ix,nA,nB,comb,d))
        check('full label space '+a+'_vs_'+b,len(comb)==(48620 if nB==9 else 24310))
    def reconstruct(counts):
        hs=[];tot=[];dom=[]
        for e,cid in enumerate(ids):
            c=counts[sr[off[e]:off[e+1]]];n=c.sum(axis=0);alln=counts[br[bo[e]:bo[e+1]]].sum(axis=0);P=np.divide(c,n[None,:],out=np.zeros_like(c,dtype=float),where=n[None,:]>0);hs.append(np.sqrt(P).T);tot.append(n);dom.append(np.divide(n,alln,out=np.zeros(26),where=alln>0))
        return np.concatenate(hs,axis=1),np.array(tot),np.array(dom)
    hr,tr,dr=reconstruct(x);check('original sqrt features reconstructed',np.max(np.abs(hr-H))<1e-14)
    thin_table=pd.read_csv(R/'local_thinning_complete_fixed_families.csv.gz');thin_globals=pd.read_csv(R/'global_thinning_all_replicates.csv');l95=pd.read_csv(R/'local_dominance95_complete_fixed_family.csv.gz');g95=pd.read_csv(R/'global_dominance95.csv');manifest=pd.read_csv(R/'thinned_count_manifest.csv');nverified=0
    def verify_family(h,valid,local,glob,name):
        nonlocal nverified
        check('fixed family size '+name,len(local)==2769 and not local.duplicated(['cluster_id','contrast']).any())
        check('invalid measurements assigned p=1 '+name,(local.loc[~local.measurement_valid,'exact_p']==1).all())
        for method,col in [('bh','BH_q_full_2769'),('by','BY_q_full_2769')]:check('SciPy '+method+' fixed family '+name,np.max(np.abs(false_discovery_control(local.exact_p.to_numpy(),method=method)-local[col]))<1e-12)
        ksum=np.zeros((26,26))
        for e in np.flatnonzero(valid):
            q=h[:,off[e]:off[e+1]];ksum+=q@q.T
        kg=ksum/valid.sum();check('PSD '+name,np.linalg.eigvalsh(kg).min()>-1e-12)
        for c,ix,nA,nB,comb,d in spaces:
            lc=local[local.contrast==c].set_index('cluster_id').loc[ids];check('independent measurement gates '+name+' '+c,np.array_equal(lc.measurement_valid.to_numpy(),valid))
            coef=1/nA+1/nB;membership=np.zeros((len(comb),len(ix)));membership[np.arange(len(comb))[:,None],comb]=1
            for e,cid in enumerate(ids):
                q=h[ix,off[e]:off[e+1]];obs=np.sum((q[d>0].mean(axis=0)-q[d<0].mean(axis=0))**2);saved=lc.iloc[e]
                assert abs(obs-saved.statistic)<1e-13,(name,c,cid,'statistic')
                if valid[e]:
                    summed_A=membership@q;delta=summed_A/nA-(q.sum(axis=0)[None,:]-summed_A)/nB;ps=np.square(delta).sum(axis=1);ge=np.count_nonzero(ps>=obs-1e-12)
                    assert ge==saved.n_ge_observed,(name,c,cid,ge,saved.n_ge_observed)
                else:assert saved.n_ge_observed==len(comb)
                nverified+=1
            g=glob[glob.contrast==c].iloc[0];k=kg[np.ix_(ix,ix)];aa=k[comb[:,:,None],comb[:,None,:]].sum(axis=(1,2));ac=k.sum(axis=1)[comb].sum(axis=1);ps=coef*coef*aa-2*coef*ac/nB+k.sum()/nB**2;t=d@k@d
            check('independent global exact '+name+' '+c,np.count_nonzero(ps>=t-1e-12)==g.n_ge_observed and abs(t-g.statistic)<1e-13)
        pg=glob.exact_p;holm=np.array([min(1,max((3-np.count_nonzero(pg<q))*q for q in pg if q<=v)) for v in pg]);check('independent Holm '+name,np.max(np.abs(holm-glob.Holm_p_3))<1e-12)
    verify_family(H,(dr>=.95-1e-12).all(axis=1),l95,g95,'dominance95')
    for row in manifest.itertuples():
        path=BASE/row.file;check('thinned file hash '+str(row.replicate),hashlib.sha256(path.read_bytes()).hexdigest()==row.sha256);thin=np.load(path)['counts'];check('full matrix thinning bounds '+str(row.replicate),thin.shape==x.shape and np.all(thin<=x))
        rng=np.random.Generator(np.random.PCG64(row.seed));again=rng.binomial(x,.5).astype(np.uint32);check('full seeded matrix reproducible '+str(row.replicate),np.array_equal(thin,again));del again
        h,t,dom=reconstruct(thin);valid=(t>=20).all(axis=1)&(dom>=.9-1e-12).all(axis=1)
        # Every duplicate junction row is identically reused in all cluster counts.
        vals=thin[sr];unique,counts=np.unique(sr,return_counts=True)
        for j in unique[counts>1]:assert np.array_equal(vals[sr==j],np.repeat(thin[j][None,:],np.count_nonzero(sr==j),axis=0))
        check('shared-junction thinning reused exactly '+str(row.replicate),True,f'{np.count_nonzero(counts>1)} reused junctions')
        verify_family(h,valid,thin_table[thin_table.replicate==row.replicate],thin_globals[thin_globals.replicate==row.replicate],f'thinning{row.replicate}')
        print('independently verified thinning',row.replicate,flush=True)
    check('all new local exact tests independently verified',nverified==11*2769,f'{nverified} tests; invalid entries independently checked as p=1')
    boot=np.load(R/'animal_bootstrap_replicates.npz');W=boot['weights'];check('whole-animal bootstrap fixed group sizes',all(np.all(W[:,groups==g].sum(axis=1)==np.count_nonzero(groups==g)) for g in ['FLT','GC','VIV']) and W.shape==(2000,26))
    rng=np.random.Generator(np.random.PCG64(2026101520));expected=np.zeros_like(W)
    for g in ['FLT','GC','VIV']:
        ix=np.flatnonzero(groups==g);draw=rng.choice(ix,size=(2000,len(ix)),replace=True)
        for r in range(2000):np.add.at(expected[r],draw[r],1)
    check('bootstrap seed reproducible',np.array_equal(expected,W))
    ci=pd.read_csv(R/'all_cluster_animal_bootstrap_intervals.csv.gz');signed=pd.read_csv(R/'all_alternative_signed_bootstrap_intervals.csv.gz');globalci=pd.read_csv(R/'global_animal_bootstrap_intervals.csv');K=H@H.T/E
    for j,(c,ix,nA,nB,comb,d) in enumerate(spaces):
        a,b=c.split('_vs_');ia=np.flatnonzero(groups==a);ib=np.flatnonzero(groups==b);delta=W[:,ia]@H[ia]/nA-W[:,ib]@H[ib]/nB;observed=H[ia].mean(axis=0)-H[ib].mean(axis=0);distance=[];error=[]
        for e,cid in enumerate(ids):
            sl=slice(off[e],off[e+1]);distance.append(np.sqrt(np.square(delta[:,sl]).sum(axis=1)/2));error.append(np.sqrt(np.square(delta[:,sl]-observed[sl]).sum(axis=1)/2))
        distance=np.array(distance).T;error=np.array(error).T;check('all bootstrap distances/errors '+c,np.max(np.abs(distance-boot['distances'][:,j]))<1e-13 and np.max(np.abs(error-boot['centered_error_distances'][:,j]))<1e-13)
        sc=ci[ci.contrast==c].set_index('cluster_id').loc[ids];q=np.quantile(distance,[.025,.975],axis=0);rad=np.quantile(error,.95,axis=0);check('all cluster bootstrap intervals '+c,np.max(np.abs(q[0]-sc.percentile_low))<1e-13 and np.max(np.abs(q[1]-sc.percentile_high))<1e-13 and np.max(np.abs(rad-sc.centered_error_radius95))<1e-13 and np.max(np.abs(np.maximum(0,sc.observed_distance-rad)-sc.confidence_ball_distance_low))<1e-13 and np.max(np.abs(np.minimum(1,sc.observed_distance+rad)-sc.confidence_ball_distance_high))<1e-13)
        nd=W[:,ia]@(H[ia]**2)/nA-W[:,ib]@(H[ib]**2)/nB;sq=np.quantile(nd,[.025,.975],axis=0);ss=signed[signed.contrast==c].set_index(['cluster_id','alternative_index']).loc[[(cid,k) for cid,n in zip(ids,np.diff(off)) for k in range(n)]];check('all signed alternative bootstrap intervals '+c,np.max(np.abs(sq[0]-ss.percentile_low))<1e-13 and np.max(np.abs(sq[1]-ss.percentile_high))<1e-13)
    diag=pd.read_csv(R/'dispersion_and_studentization_diagnostics.csv');disp=pd.read_csv(R/'pooled_group_dispersion.csv').set_index('group');localdisp=pd.read_csv(R/'all_cluster_group_dispersion.csv.gz');pg=[];pt=[]
    for c,ix,nA,nB,comb,d in spaces:
        a,b=c.split('_vs_');h=H[ix];k=K[np.ix_(ix,ix)];sel=np.flatnonzero(d>0);rem=np.flatnonzero(d<0)
        Va=np.square(h[sel]-h[sel].mean(axis=0)).sum()/(nA-1)/E;Vb=np.square(h[rem]-h[rem].mean(axis=0)).sum()/(nB-1)/E;t=np.square(h[sel].mean(axis=0)-h[rem].mean(axis=0)).sum()/E;obsd=abs(Va-Vb)/2;obst=t/(Va/nA+Vb/nB);row=diag[diag.contrast==c].iloc[0]
        check('direct dispersion and studentized statistic '+c,abs(Va/2-row.variance_A_squared_hellinger)<1e-13 and abs(Vb/2-row.variance_B_squared_hellinger)<1e-13 and abs(obst-row.studentized_centroid_statistic)<1e-12)
        nd=nt=0
        for sel in comb:
            rem=np.setdiff1d(np.arange(len(ix)),sel);A=k[np.ix_(sel,sel)];B=k[np.ix_(rem,rem)];C=k[np.ix_(sel,rem)];VA=(np.trace(A)-A.sum()/nA)/(nA-1);VB=(np.trace(B)-B.sum()/nB)/(nB-1);T=A.sum()/nA**2+B.sum()/nB**2-2*C.sum()/(nA*nB)
            nd+=abs(VA-VB)/2>=obsd-1e-12;nt+=T/(VA/nA+VB/nB)>=obst-1e-12
        pg.append(nd/len(comb));pt.append(nt/len(comb));check('separate direct dispersion/studentized full spaces '+c,abs(nd/len(comb)-row.dispersion_exact_reference_p)<1e-14 and abs(nt/len(comb)-row.studentized_exact_reference_p)<1e-14)
    for probs,col in [(np.array(pg),'dispersion_Holm_3'),(np.array(pt),'studentized_Holm_3')]:
        adj=np.array([min(1,max((3-np.count_nonzero(probs<q))*q for q in probs if q<=v)) for v in probs]);check('independent Holm '+col,np.max(np.abs(adj-diag[col]))<1e-12)
    loo=pd.read_csv(R/'LOO_crosscheck.csv');prior=pd.read_csv(old/'tables/leave_one_animal_out_effects.csv');merged=loo.merge(prior[['contrast','omitted_sample','hellinger_centroid_distance']],on=['contrast','omitted_sample'],suffixes=['_new','_old'],validate='one_to_one');check('inherited leave-one-animal-out effects crosschecked',len(merged)==52 and np.max(np.abs(merged.hellinger_centroid_distance_new-merged.hellinger_centroid_distance_old))<1e-13)
    result=dict(all_passed=True,n_checks=len(checks),n_local_tests_verified=nverified,checks=checks);(R/'independent_verification.json').write_text(json.dumps(result,indent=2));pd.DataFrame(checks).to_csv(R/'independent_verification.csv',index=False);print('ALL CHECKS PASS',len(checks),flush=True)

if __name__=='__main__':main()
