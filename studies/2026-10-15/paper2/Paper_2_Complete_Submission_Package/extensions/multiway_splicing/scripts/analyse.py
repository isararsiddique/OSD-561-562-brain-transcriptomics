"""Exact fixed-size permutation inference on frozen local sqrt compositions."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
from pathlib import Path
import itertools, math, json, argparse
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1]; ROOT=BASE.parents[2]
TOL=1e-12
CONTRASTS=[('FLT','GC'),('FLT','VIV'),('GC','VIV')]

def adjust(p,method):
    p=np.asarray(p); order=np.argsort(p,kind='stable'); ranked=p[order]; m=len(p)
    if method=='holm': q=np.maximum.accumulate(ranked*(m-np.arange(m)))
    else:
        fac=sum(1/i for i in range(1,m+1)) if method=='by' else 1
        q=np.minimum.accumulate((ranked*m*fac/(np.arange(m)+1))[::-1])[::-1]
    out=np.empty(m);out[order]=np.minimum(q,1);return out

def assignments(n,a):
    comb=np.array(list(itertools.combinations(range(n),a)),dtype=np.int16)
    d=np.full((len(comb),n),-1/(n-a),dtype=float)
    d[np.arange(len(comb))[:,None],comb]=1/a
    assert len(d)==math.comb(n,a)
    return d

def exact(K,ix,groups,a,b,D):
    lab=groups[ix]; nA=(lab==a).sum(); nB=(lab==b).sum()
    d=np.where(lab==a,1/nA,-1/nB)
    k=K[np.ix_(ix,ix)]; obs=float(d@k@d)
    perm=np.einsum('ij,ij->i',D@k,D)
    ge=int((perm>=obs-TOL).sum())
    return dict(contrast=a+'_vs_'+b,n_A=int(nA),n_B=int(nB),statistic=obs,hellinger_centroid_distance=float(np.sqrt(max(obs,0)/2)),exact_p=ge/len(perm),n_assignments=len(perm),n_ge_observed=ge),perm,d

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--base',type=Path,default=ROOT/'research/paper2_beta_20261014/inputs');args=parser.parse_args()
    z=np.load(BASE/'tables/frozen_features.npz'); H=z['H']; offsets=z['offsets']; samples=z['samples']; ids=z['cluster_ids']; E=len(ids)
    defs=pd.read_csv(BASE/'tables/selected_cluster_definitions.csv').set_index('cluster_id').loc[ids]
    meta=pd.read_csv(args.base/'animal_metadata.csv').set_index('sample').loc[samples]
    meta[['Source Name','group','Comment[Euthanasia Date]','Comment[BSP Dissection Date]']].to_csv(BASE/'tables/analysis_animals.csv')
    groups=meta.group.to_numpy()
    K=(H@H.T)/E
    local_K=np.stack([H[:,offsets[e]:offsets[e+1]]@H[:,offsets[e]:offsets[e+1]].T for e in range(E)])
    assert np.allclose(K,local_K.mean(axis=0),atol=1e-14)
    pd.DataFrame(K,index=samples,columns=samples).to_csv(BASE/'tables/equal_cluster_animal_kernel.csv')
    mapped=(defs.gene_mapping_count==1).to_numpy(); dom95=(defs.min_dominance>=.95-TOL).to_numpy()
    gm=defs.loc[mapped,'gene_id']; freq=gm.value_counts(); weights=gm.map(lambda g:1/freq[g]).to_numpy(float);weights/=weights.sum()
    kernels={'equal_cluster_primary':K,'mapped_subset_equal_cluster':local_K[mapped].mean(axis=0),'equal_gene_weight':np.einsum('e,eij->ij',weights,local_K[mapped]),'dominance_95':local_K[dom95].mean(axis=0)}
    primary=[]; sensitivities=[]; locals_out=[]; distributions=[]; compmeans=[]; loo=[]
    for a,b in CONTRASTS:
        ix=np.flatnonzero(np.isin(groups,[a,b])); nA=int((groups[ix]==a).sum()); D=assignments(len(ix),nA)
        row,perm,d=exact(K,ix,groups,a,b,D);row.update(n_clusters=E);primary.append(row)
        distributions.append(pd.DataFrame(dict(contrast=a+'_vs_'+b,assignment_index=np.arange(len(perm)),statistic=perm,hellinger_centroid_distance=np.sqrt(np.maximum(perm,0)/2))))
        for name,k in kernels.items():
            if name=='equal_cluster_primary':continue
            r,_,_=exact(k,ix,groups,a,b,D)
            r.update(sensitivity=name,n_clusters=int(mapped.sum()) if name in ['mapped_subset_equal_cluster','equal_gene_weight'] else int(dom95.sum()),n_genes=int(gm.nunique()) if name in ['mapped_subset_equal_cluster','equal_gene_weight'] else np.nan)
            sensitivities.append(r)
        h=H[ix];obs=np.add.reduceat((d@h)**2,offsets[:-1])
        ge=np.zeros(E,dtype=np.int64)
        for start in range(0,len(D),1024):
            proj=D[start:start+1024]@h
            stats=np.add.reduceat(proj*proj,offsets[:-1],axis=1)
            ge+=(stats>=obs[None,:]-TOL).sum(axis=0)
        native=H*H
        for e,cid in enumerate(ids):
            sl=slice(offsets[e],offsets[e+1]); ma=native[groups==a,sl].mean(axis=0);mb=native[groups==b,sl].mean(axis=0);delta=ma-mb
            locals_out.append(dict(cluster_id=cid,contrast=a+'_vs_'+b,statistic=float(obs[e]),hellinger_centroid_distance=float(np.sqrt(max(obs[e],0)/2)),native_mean_total_variation=float(np.abs(delta).sum()/2),max_abs_proportion_difference=float(np.abs(delta).max()),exact_p=ge[e]/len(D),n_ge_observed=int(ge[e]),n_assignments=len(D)))
            for j,(av,bv,dv) in enumerate(zip(ma,mb,delta)):compmeans.append(dict(cluster_id=cid,contrast=a+'_vs_'+b,alternative_index=j,mean_proportion_A=av,mean_proportion_B=bv,delta_A_minus_B=dv))
        for omit in ix:
            yi=ix[ix!=omit];lab=groups[yi];dy=np.where(lab==a,1/(lab==a).sum(),-1/(lab==b).sum()); t=float(dy@K[np.ix_(yi,yi)]@dy)
            loo.append(dict(contrast=a+'_vs_'+b,omitted_sample=samples[omit],omitted_group=groups[omit],n_A=int((lab==a).sum()),n_B=int((lab==b).sum()),statistic=t,hellinger_centroid_distance=np.sqrt(max(t,0)/2)))
        print(row,flush=True)
    primary=pd.DataFrame(primary);primary['holm_p_3']=adjust(primary.exact_p,'holm');primary['holm_significant_05']=primary.holm_p_3<=.05
    sens=pd.DataFrame(sensitivities)
    for name,idx in sens.groupby('sensitivity').groups.items():sens.loc[idx,'holm_p_within_sensitivity_3']=adjust(sens.loc[idx,'exact_p'],'holm')
    local=pd.DataFrame(locals_out)
    local['BH_q_full_3E']=adjust(local.exact_p,'bh');local['BY_q_full_3E']=adjust(local.exact_p,'by')
    for c,idx in local.groupby('contrast').groups.items():
        local.loc[idx,'BH_q_within_contrast_E']=adjust(local.loc[idx,'exact_p'],'bh');local.loc[idx,'BY_q_within_contrast_E']=adjust(local.loc[idx,'exact_p'],'by')
    local['BY_significant_full_family_05']=local.BY_q_full_3E<=.05
    local=local.merge(defs.reset_index(),on='cluster_id',validate='many_to_one')
    primary.to_csv(BASE/'tables/global_exact_results.csv',index=False)
    sens.to_csv(BASE/'tables/global_sensitivity_results.csv',index=False)
    local.to_csv(BASE/'tables/local_exact_results_complete.csv.gz',index=False)
    pd.concat(distributions,ignore_index=True).to_csv(BASE/'tables/global_exact_null_distributions.csv.gz',index=False)
    pd.DataFrame(compmeans).to_csv(BASE/'tables/alternative_group_mean_differences.csv.gz',index=False)
    pd.DataFrame(loo).to_csv(BASE/'tables/leave_one_animal_out_effects.csv',index=False)
    eig=np.linalg.eigvalsh(K); cen=np.eye(26)-np.ones((26,26))/26; w,v=np.linalg.eigh(cen@K@cen);order=np.argsort(w)[::-1]; w=w[order];v=v[:,order];scores=v[:,:2]*np.sqrt(np.maximum(w[:2],0))
    pd.DataFrame(dict(sample=samples,group=groups,PCo1=scores[:,0],PCo2=scores[:,1])).to_csv(BASE/'tables/kernel_principal_coordinates.csv',index=False)
    summary=dict(n_clusters=E,n_local_tests=len(local),n_local_BH_full=int((local.BH_q_full_3E<=.05).sum()),n_local_BY_full=int(local.BY_significant_full_family_05.sum()),n_unique_BY_clusters=int(local.loc[local.BY_significant_full_family_05,'cluster_id'].nunique()),n_primary_Holm_positive=int(primary.holm_significant_05.sum()),primary_kernel_min_eigenvalue=float(eig.min()),primary_kernel_max_diagonal_error=float(np.max(np.abs(np.diag(K)-1))),pcoa_variance_fraction=(w[:2]/w[w>1e-12].sum()).tolist(),tie_tolerance=TOL,n_dom95_clusters=int(dom95.sum()),n_mapped_clusters=int(mapped.sum()),n_mapped_genes=int(gm.nunique()))
    (BASE/'tables/inference_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
