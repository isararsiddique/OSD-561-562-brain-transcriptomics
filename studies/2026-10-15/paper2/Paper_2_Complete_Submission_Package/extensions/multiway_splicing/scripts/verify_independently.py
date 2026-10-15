"""Independent checks: original counts, closed-form kernel sums, SciPy FDR."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
from collections import defaultdict
import argparse, itertools, json, math, ast
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
BASE=Path(__file__).resolve().parents[1];ROOT=BASE.parents[2]

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=ROOT/'research/paper2_beta_20261014/inputs');args=p.parse_args()
    checks=[]
    def record(name,value,detail=''):
        assert bool(value),f'{name}: {detail}'
        checks.append(dict(check=name,passed=True,detail=detail))
    df=pd.read_csv(args.base/'all_animal_junction_counts.csv.gz',index_col=0);x=df.to_numpy()
    defs=pd.read_csv(BASE/'tables/selected_cluster_definitions.csv').set_index('cluster_id')
    candidate=pd.read_csv(BASE/'tables/all_candidate_clusters.csv.gz').set_index('cluster_id')
    lookup=pd.read_csv(BASE/'annotation/junction_gene_lookup.csv.gz').fillna('').set_index('junction_id')
    qualifying=np.count_nonzero(x>=5,axis=1)>=3
    parts=[j.split(':') for j in df.index]; retained={}; independent_candidates=set()
    for b,bcol in [('start',2),('end',3)]:
        q=defaultdict(list);keys=[':'.join([t[0],t[1],t[bcol]]) for t in parts]
        for i,key in enumerate(keys):
            if qualifying[i]:q[key].append(i)
        q={k:v for k,v in q.items() if len(v)>=3};totals={k:np.zeros(26,dtype=np.int64) for k in q}
        for i,key in enumerate(keys):
            if key in totals:totals[key]+=x[i]
        for key,ii in q.items():
            cid=key+':'+b;independent_candidates.add(cid);y=x[ii];s=y.sum(axis=0)
            dominance=np.divide(s,totals[key],out=np.zeros(26),where=totals[key]>0);proportions=np.divide(y,s[None,:],out=np.zeros_like(y,dtype=float),where=s[None,:]>0)
            means=proportions.mean(axis=1)
            keep=np.min(s)>=20 and np.min(dominance)>=.9-1e-12 and np.min(means)>=.01-1e-12 and np.count_nonzero(means>=.05-1e-12)>=3 and np.max(means)<=.9+1e-12
            # Reverse animal columns and apply all pooled filters afresh.
            yr=y[:,::-1];sr=yr.sum(axis=0);dr=np.divide(sr,totals[key][::-1],out=np.zeros(26),where=totals[key][::-1]>0);mr=np.divide(yr,sr[None,:],out=np.zeros_like(yr,dtype=float),where=sr[None,:]>0).mean(axis=1)
            reverse_keep=np.min(sr)>=20 and np.min(dr)>=.9-1e-12 and np.min(mr)>=.01-1e-12 and np.count_nonzero(mr>=.05-1e-12)>=3 and np.max(mr)<=.9+1e-12
            record_flag=keep==reverse_keep
            if not record_flag:raise AssertionError('animal order changed selection '+cid)
            if keep:retained[cid]=ii
        print('independent selection',b,len(retained),flush=True)
    record('candidate boundary completeness',independent_candidates==set(candidate.index),f'{len(independent_candidates)} candidates')
    record('independent frozen selection equals all definitions',set(retained)==set(defs.index),f'{len(retained)} retained clusters')
    record('all candidate decisions invariant to reversed animal order',True,f'{len(independent_candidates)} evaluated; labels never enter these functions')
    tree=ast.parse((BASE/'scripts/select_clusters.py').read_text());read_args=[ast.unparse(n.args[0]) for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='read_csv']
    record('selector does not read animal group metadata',set(read_args)=={'counts_path','lookup_path'},str(read_args))
    z=np.load(BASE/'tables/frozen_features.npz');H=z['H'];off=z['offsets'];ids=z['cluster_ids'];samples=z['samples'];E=len(ids)
    record('frozen sample order matches original matrix',np.array_equal(samples,df.columns))
    feature_errors=[];mapping_errors=[]
    for e,cid in enumerate(ids):
        ii=retained[cid];native=x[ii]/x[ii].sum(axis=0)[None,:]
        feature_errors.append(np.max(np.abs(H[:,off[e]:off[e+1]]**2-native.T)))
        sets=[set(filter(None,lookup.loc[df.index[i],'matching_gene_ids'].split(';'))) for i in ii];common=set.intersection(*sets)
        row=defs.loc[cid];actual_gene=next(iter(common)) if len(common)==1 else ''
        mapping_errors.append(len(common)!=row.gene_mapping_count or actual_gene!=(row.gene_id if pd.notna(row.gene_id) else ''))
    record('all local sqrt features reconstruct original counts',max(feature_errors)<1e-14,f'max error {max(feature_errors):.3g}')
    record('all-alternative intersection annotation',not any(mapping_errors),f'{E} cluster maps checked against bundled lookup; raw GTF provenance is hashed')
    K=H@H.T/E;eig=np.linalg.eigvalsh(K)
    record('primary animal kernel PSD',eig.min()>=-1e-12,f'min eigenvalue {eig.min():.12g}')
    record('equal-cluster kernel unit diagonal',np.max(np.abs(np.diag(K)-1))<1e-12)
    meta=pd.read_csv(args.base/'animal_metadata.csv').set_index('sample').loc[samples];groups=meta.group.to_numpy()
    glob=pd.read_csv(BASE/'tables/global_exact_results.csv');local=pd.read_csv(BASE/'tables/local_exact_results_complete.csv.gz');sens=pd.read_csv(BASE/'tables/global_sensitivity_results.csv')
    local_kernels=np.stack([H[:,off[e]:off[e+1]]@H[:,off[e]:off[e+1]].T for e in range(E)])
    ordered=defs.loc[ids];mapped=(ordered.gene_mapping_count==1).to_numpy();gm=ordered.loc[mapped,'gene_id'];ng=gm.nunique();w=np.array([1/(ng*np.count_nonzero(gm==g)) for g in gm])
    secondary={'mapped_subset_equal_cluster':local_kernels[mapped].mean(axis=0),'equal_gene_weight':np.sum(local_kernels[mapped]*w[:,None,None],axis=0),'dominance_95':local_kernels[(ordered.min_dominance>=.95-1e-12).to_numpy()].mean(axis=0)}
    record('equal-gene weights give each mapped gene equal total weight',all(abs(w[(gm==g).to_numpy()].sum()-1/ng)<1e-14 for g in gm.unique()))
    for name,ks in secondary.items():record('secondary PSD '+name,np.linalg.eigvalsh(ks).min()>=-1e-12)
    distro=pd.read_csv(BASE/'tables/global_exact_null_distributions.csv.gz')
    local_checked=0
    for row in glob.itertuples():
        a,b=row.contrast.split('_vs_');ix=np.flatnonzero(np.isin(groups,[a,b]));nA=(groups[ix]==a).sum();nB=len(ix)-nA;k=K[np.ix_(ix,ix)]
        sel=np.flatnonzero(groups[ix]==a);rem=np.flatnonzero(groups[ix]==b)
        direct=((H[ix[sel]].mean(axis=0)-H[ix[rem]].mean(axis=0))**2).sum()/E
        record('direct feature equals kernel statistic '+row.contrast,abs(direct-row.statistic)<1e-13)
        # Independent fixed-size enumeration from closed-form within/across sums.
        total=k.sum();rowsums=k.sum(axis=1);coef=(1/nA+1/nB)
        brute=[];nobserved=0
        for c in itertools.combinations(range(len(ix)),int(nA)):
            c=np.array(c);aa=k[np.ix_(c,c)].sum();ac=rowsums[c].sum()
            brute.append(coef*coef*aa-2*coef*ac/nB+total/(nB*nB))
            nobserved+=np.array_equal(c,sel)
        brute=np.array(brute);saved=distro.loc[distro.contrast==row.contrast,'statistic'].to_numpy()
        record('closed-form full exact distribution '+row.contrast,len(brute)==math.comb(len(ix),int(nA)) and np.max(np.abs(brute-saved))<2e-13,f'{len(brute)} assignments, maximum difference {np.max(np.abs(brute-saved)):.3g}')
        record('observed assignment included once '+row.contrast,nobserved==1)
        record('independent exact global p '+row.contrast,np.count_nonzero(brute>=direct-1e-12)==row.n_ge_observed)
        # Independently check each local event with its own small Gram matrix.
        comb=np.array(list(itertools.combinations(range(len(ix)),int(nA))))
        for name,ks in secondary.items():
            sk=ks[np.ix_(ix,ix)];sa=sk[comb[:,:,None],comb[:,None,:]].sum(axis=(1,2));sc=sk.sum(axis=1)[comb].sum(axis=1)
            sp=coef*coef*sa-2*coef*sc/nB+sk.sum()/(nB*nB)
            sd=np.where(groups[ix]==a,1/nA,-1/nB);st=sd@sk@sd
            sr=sens[(sens.sensitivity==name)&(sens.contrast==row.contrast)].iloc[0]
            record('secondary exact p '+row.contrast+' '+name,np.count_nonzero(sp>=st-1e-12)==sr.n_ge_observed and abs(st-sr.statistic)<1e-13)
        for e,cid in enumerate(ids):
            h=H[ix,off[e]:off[e+1]];kl=h@h.T;aa=kl[comb[:,:,None],comb[:,None,:]].sum(axis=(1,2));ac=kl.sum(axis=1)[comb].sum(axis=1)
            pl=coef*coef*aa-2*coef*ac/nB+kl.sum()/(nB*nB)
            tl=np.sum((h[sel].mean(axis=0)-h[rem].mean(axis=0))**2)
            stored=local[(local.contrast==row.contrast)&(local.cluster_id==cid)].iloc[0]
            assert abs(tl-stored.statistic)<1e-13
            assert np.count_nonzero(pl>=tl-1e-12)==stored.n_ge_observed,(cid,row.contrast)
            local_checked+=1
        print('independently checked',row.contrast,flush=True)
    record('every local exact p verified separately',local_checked==3*E,f'{local_checked} tests, complete assignment spaces, independent closed-form sums')
    for name,col in [('bh','BH_q_full_3E'),('by','BY_q_full_3E')]:
        expected=false_discovery_control(local.exact_p.to_numpy(),method=name)
        record('SciPy '+name+' full family',np.max(np.abs(expected-local[col]))<1e-12,f'family size {len(local)}')
    expected_holm=np.array([min(1,max((3-sum(glob.exact_p<q))*q for q in glob.exact_p if q<=v)) for v in glob.exact_p])
    record('independent Holm primary three',np.max(np.abs(expected_holm-glob.holm_p_3))<1e-12)
    for name,sub in sens.groupby('sensitivity'):
        sh=np.array([min(1,max((3-sum(sub.exact_p<q))*q for q in sub.exact_p if q<=v)) for v in sub.exact_p])
        record('independent Holm secondary '+name,np.max(np.abs(sh-sub.holm_p_within_sensitivity_3))<1e-12)
    for c,sub in local.groupby('contrast'):
        for name,col in [('bh','BH_q_within_contrast_E'),('by','BY_q_within_contrast_E')]:record('SciPy '+name+' '+c,np.max(np.abs(false_discovery_control(sub.exact_p.to_numpy(),method=name)-sub[col]))<1e-12)
    result=dict(all_passed=True,n_checks=len(checks),checks=checks)
    (BASE/'tables/independent_verification.json').write_text(json.dumps(result,indent=2))
    pd.DataFrame(checks).to_csv(BASE/'tables/independent_verification.csv',index=False)
    print('ALL CHECKS PASS',len(checks),flush=True)

if __name__=='__main__':main()
