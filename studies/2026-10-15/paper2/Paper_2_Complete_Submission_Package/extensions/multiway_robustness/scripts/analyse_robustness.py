"""Frozen-family local robustness; full-matrix binomial thinning; animal bootstrap."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import argparse,itertools,json,hashlib,math
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1];PACKAGE=BASE.parents[1]
CONTRASTS=[('FLT','GC'),('FLT','VIV'),('GC','VIV')];TOL=1e-12

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()

def adjust(p,method):
    p=np.asarray(p);m=len(p);order=np.argsort(p,kind='stable');r=p[order]
    if method=='holm':q=np.maximum.accumulate(r*(m-np.arange(m)))
    else:q=np.minimum.accumulate((r*m/(np.arange(m)+1)*(sum(1/i for i in range(1,m+1)) if method=='by' else 1))[::-1])[::-1]
    out=np.empty(m);out[order]=np.minimum(q,1);return out

def labels(groups):
    out=[]
    for a,b in CONTRASTS:
        ix=np.flatnonzero(np.isin(groups,[a,b]));nA=np.count_nonzero(groups[ix]==a);nB=len(ix)-nA;comb=np.array(list(itertools.combinations(range(len(ix)),nA)),dtype=np.int16)
        D=np.full((len(comb),len(ix)),-1/nB);D[np.arange(len(comb))[:,None],comb]=1/nA
        d=np.where(groups[ix]==a,1/nA,-1/nB)
        out.append((a+'_vs_'+b,ix,nA,nB,comb,D,d))
    return out

def exact_family(H,off,valid,spaces,ids,rep):
    E=len(ids);locals=[];globals=[];K=H@H.T/E
    if not valid.all():K=sum(H[:,off[e]:off[e+1]]@H[:,off[e]:off[e+1]].T for e in np.flatnonzero(valid))/valid.sum()
    for contrast,ix,nA,nB,comb,D,d in spaces:
        h=H[ix];obs=np.add.reduceat((d@h)**2,off[:-1]);ge=np.zeros(E,dtype=np.int64)
        for start in range(0,len(D),1024):
            q=D[start:start+1024]@h;s=np.add.reduceat(q*q,off[:-1],axis=1);ge+=(s>=obs[None,:]-TOL).sum(axis=0)
        for e,cid in enumerate(ids):locals.append(dict(replicate=rep,contrast=contrast,cluster_id=cid,measurement_valid=bool(valid[e]),statistic=obs[e],hellinger_centroid_distance=np.sqrt(max(obs[e],0)/2),n_ge_observed=int(ge[e]) if valid[e] else len(D),n_assignments=len(D),exact_p=ge[e]/len(D) if valid[e] else 1.))
        k=K[np.ix_(ix,ix)];t=float(d@k@d);per=np.einsum('ij,ij->i',D@k,D);nge=np.count_nonzero(per>=t-TOL)
        globals.append(dict(replicate=rep,contrast=contrast,n_clusters=int(valid.sum()),statistic=t,hellinger_centroid_distance=np.sqrt(max(t,0)/2),n_ge_observed=int(nge),n_assignments=len(D),exact_p=nge/len(D)))
    local=pd.DataFrame(locals);local['BH_q_full_2769']=adjust(local.exact_p,'bh');local['BY_q_full_2769']=adjust(local.exact_p,'by');local['BY_positive']=local.BY_q_full_2769<=.05
    glob=pd.DataFrame(globals);glob['Holm_p_3']=adjust(glob.exact_p,'holm')
    return local,glob

def features(counts,selected_rows,off,boundary_rows,boundary_off):
    y=counts[selected_rows];total=np.add.reduceat(y,off[:-1],axis=0);canonical=np.add.reduceat(counts[boundary_rows],boundary_off[:-1],axis=0)
    dom=np.divide(total,canonical,out=np.zeros_like(total,dtype=float),where=canonical>0)
    repeat=np.repeat(total,np.diff(off),axis=0);P=np.divide(y,repeat,out=np.zeros_like(y,dtype=float),where=repeat>0)
    return np.sqrt(P).T,total,dom

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=PACKAGE);args=p.parse_args();inp=args.base/'inputs';old=args.base/'extensions/multiway_splicing'
    assert (BASE/'RETROSPECTIVE_PLAN_2026-10-15.md').exists()
    df=pd.read_csv(inp/'all_animal_junction_counts.csv.gz',index_col=0);x=df.to_numpy(dtype=np.int64);assert x.shape==(805562,26)
    z=np.load(old/'tables/frozen_features.npz');H=z['H'];off=z['offsets'];ids=z['cluster_ids'];samples=z['samples'];E=len(ids);assert E==923 and np.array_equal(samples,df.columns)
    defs=pd.read_csv(old/'tables/selected_cluster_definitions.csv').set_index('cluster_id').loc[ids];meta=pd.read_csv(inp/'animal_metadata.csv').set_index('sample').loc[samples];groups=meta.group.to_numpy();spaces=labels(groups)
    a=pd.read_csv(old/'tables/animal_alternative_compositions.csv.gz')[['cluster_id','alternative_index','junction_id']].drop_duplicates().set_index(['cluster_id','alternative_index'])
    lookup={j:i for i,j in enumerate(df.index)};selected_rows=np.array([lookup[a.loc[(cid,j),'junction_id']] for cid,n in zip(ids,np.diff(off)) for j in range(n)],dtype=np.int64)
    boundarymap={cid:[] for cid in ids};keys=set(ids)
    for i,j in enumerate(df.index):
        ch,st,start,end=j.split(':')
        k=f'{ch}:{st}:{start}:start'
        if k in keys:boundarymap[k].append(i)
        k=f'{ch}:{st}:{end}:end'
        if k in keys:boundarymap[k].append(i)
    boundary_rows=np.array([i for cid in ids for i in boundarymap[cid]],dtype=np.int64);boundary_off=np.r_[0,np.cumsum([len(boundarymap[cid]) for cid in ids])]
    np.savez_compressed(BASE/'results/frozen_count_index.npz',selected_rows=selected_rows,offsets=off,boundary_rows=boundary_rows,boundary_offsets=boundary_off,cluster_ids=ids,samples=samples)
    h,tot,dom=features(x,selected_rows,off,boundary_rows,boundary_off);assert np.max(np.abs(h-H))<1e-14
    original=pd.read_csv(old/'tables/local_exact_results_complete.csv.gz');hits=original[original.BY_q_full_3E<=.05].copy();assert len(hits)==43 and hits.cluster_id.nunique()==41
    hits.to_csv(BASE/'results/original_43_BY_tests.csv',index=False)
    passing95=(dom>=.95-TOL).all(axis=1);l95,g95=exact_family(H,off,passing95,spaces,ids,'dominance95')
    for r in l95.itertuples():
        if r.measurement_valid:assert abs(r.exact_p-original[(original.cluster_id==r.cluster_id)&(original.contrast==r.contrast)].exact_p.iloc[0])<1e-14
    l95.to_csv(BASE/'results/local_dominance95_complete_fixed_family.csv.gz',index=False);g95.to_csv(BASE/'results/global_dominance95.csv',index=False)
    counts_hash=sha(inp/'all_animal_junction_counts.csv.gz');repglobals=[];repslocals=[];measurement=[];thin_manifest=[]
    for rep in range(1,11):
        seed=2026101500+rep;rng=np.random.Generator(np.random.PCG64(seed));thin=rng.binomial(x,.5).astype(np.uint32)
        path=BASE/'thinned_counts'/f'replicate_{rep:02d}_seed_{seed}.npz';np.savez_compressed(path,counts=thin,samples=samples)
        hr,tr,dr=features(thin,selected_rows,off,boundary_rows,boundary_off);valid=(tr>=20).all(axis=1)&(dr>=.90-TOL).all(axis=1)
        lr,gr=exact_family(hr,off,valid,spaces,ids,rep);repslocals.append(lr);repglobals.append(gr)
        for e,cid in enumerate(ids):measurement.append(dict(replicate=rep,cluster_id=cid,min_selected_reads=int(tr[e].min()),min_dominance=float(dr[e].min()),coverage_pass=bool((tr[e]>=20).all()),dominance_pass=bool((dr[e]>=.90-TOL).all()),measurement_valid=bool(valid[e])))
        thin_manifest.append(dict(replicate=rep,seed=seed,probability=.5,file=str(path.relative_to(BASE)),sha256=sha(path),counts_input_sha256=counts_hash,n_clusters_valid=int(valid.sum()),original_total_counts=int(x.sum()),thinned_total_counts=int(thin.sum())))
        print('thinning',rep,'valid',valid.sum(),'BY',lr.BY_positive.sum(),flush=True)
    allthin=pd.concat(repslocals,ignore_index=True);allthin.to_csv(BASE/'results/local_thinning_complete_fixed_families.csv.gz',index=False)
    pd.concat(repglobals,ignore_index=True).to_csv(BASE/'results/global_thinning_all_replicates.csv',index=False);pd.DataFrame(measurement).to_csv(BASE/'results/thinning_measurement_gates.csv.gz',index=False);pd.DataFrame(thin_manifest).to_csv(BASE/'results/thinned_count_manifest.csv',index=False)
    track=allthin.merge(hits[['cluster_id','contrast','BY_q_full_3E']],on=['cluster_id','contrast'],how='inner',validate='many_to_one');track.rename(columns={'BY_q_full_3E':'original_BY_q'},inplace=True);track.to_csv(BASE/'results/original_hit_thinning_tracking_complete.csv',index=False)
    st=track.groupby(['cluster_id','contrast']).agg(n_replicates=('replicate','size'),n_measurement_valid=('measurement_valid','sum'),n_BY_positive=('BY_positive','sum'),min_BY_q=('BY_q_full_2769','min'),max_BY_q=('BY_q_full_2769','max')).reset_index();st=st.merge(l95[['cluster_id','contrast','measurement_valid','BY_q_full_2769','BY_positive']],on=['cluster_id','contrast'],validate='one_to_one');st.rename(columns={'measurement_valid':'valid_dominance95','BY_q_full_2769':'BY_q_dominance95','BY_positive':'BY_positive_dominance95'},inplace=True);st.to_csv(BASE/'results/original_hit_stability_summary.csv',index=False)
    # Joint animal bootstrap: count weights are common to every cluster/contrast.
    rng=np.random.Generator(np.random.PCG64(2026101520));B=2000;W=np.zeros((B,26),dtype=np.int16)
    for g in ['FLT','GC','VIV']:
        ix=np.flatnonzero(groups==g);draw=rng.choice(ix,size=(B,len(ix)),replace=True)
        for r in range(B):np.add.at(W[r],draw[r],1)
    bdist=np.empty((B,3,E));berr=np.empty_like(bdist);global_dist=np.empty((B,3));global_err=np.empty_like(global_dist);cis=[];signed=[]
    for ci,(c,ix,nA,nB,comb,D,d) in enumerate(spaces):
        ga,gb=c.split('_vs_');ia=np.flatnonzero(groups==ga);ib=np.flatnonzero(groups==gb);Db=W[:,ia]/nA;Db2=W[:,ib]/nB;delta=Db@H[ia]-Db2@H[ib];observed=d@H[ix]
        t=np.add.reduceat(delta*delta,off[:-1],axis=1);error=np.add.reduceat((delta-observed[None,:])**2,off[:-1],axis=1)
        bdist[:,ci,:]=np.sqrt(t/2);berr[:,ci,:]=np.sqrt(error/2);global_dist[:,ci]=np.sqrt(t.mean(axis=1)/2);global_err[:,ci]=np.sqrt(error.mean(axis=1)/2)
        native=H*H;native_d=Db@native[ia]-Db2@native[ib];native_obs=d@native[ix]
        for e,cid in enumerate(ids):
            obs=np.sqrt(np.sum(observed[off[e]:off[e+1]]**2)/2);lo,hi=np.quantile(bdist[:,ci,e],[.025,.975]);radius=np.quantile(berr[:,ci,e],.95)
            cis.append(dict(cluster_id=cid,contrast=c,observed_distance=obs,percentile_low=lo,percentile_high=hi,centered_error_radius95=radius,confidence_ball_distance_low=max(0,obs-radius),confidence_ball_distance_high=min(1,obs+radius),bootstrap_replicates=B))
            for j,k in enumerate(range(off[e],off[e+1])):
                q=np.quantile(native_d[:,k],[.025,.975]);signed.append(dict(cluster_id=cid,alternative_index=j,contrast=c,observed_delta_A_minus_B=native_obs[k],percentile_low=q[0],percentile_high=q[1],bootstrap_replicates=B))
    pd.DataFrame(cis).to_csv(BASE/'results/all_cluster_animal_bootstrap_intervals.csv.gz',index=False);pd.DataFrame(signed).to_csv(BASE/'results/all_alternative_signed_bootstrap_intervals.csv.gz',index=False)
    np.savez_compressed(BASE/'results/animal_bootstrap_replicates.npz',weights=W,distances=bdist,centered_error_distances=berr,global_distances=global_dist,global_error_distances=global_err,cluster_ids=ids,contrasts=np.array([s[0] for s in spaces]),seed=np.array(2026101520))
    # Group dispersion is in squared Hellinger units, bias adjusted n/(n-1).
    K=H@H.T/E;localdisp=[];globaldisp=[]
    for g in ['FLT','GC','VIV']:
        ii=np.flatnonzero(groups==g);n=len(ii);variance=(np.trace(K[np.ix_(ii,ii)])-K[np.ix_(ii,ii)].sum()/n)/(n-1)/2
        boot= W[:,ii]@K[np.ix_(ii,ii)];quad=np.einsum('ij,ij->i',boot,W[:,ii]);vb=(n-quad/n)/(n-1)/2
        qs=np.quantile(vb,[.025,.975]);globaldisp.append(dict(group=g,n_animals=n,unbiased_mean_squared_hellinger_radius=variance,bootstrap_percentile_low=qs[0],bootstrap_percentile_high=qs[1]))
        for e,cid in enumerate(ids):
            he=H[ii,off[e]:off[e+1]];vg=((he*he).sum()-np.square(he.sum(axis=0)).sum()/n)/(n-1)/2;localdisp.append(dict(cluster_id=cid,group=g,n_animals=n,unbiased_mean_squared_hellinger_radius=vg))
    pd.DataFrame(globaldisp).to_csv(BASE/'results/pooled_group_dispersion.csv',index=False);pd.DataFrame(localdisp).to_csv(BASE/'results/all_cluster_group_dispersion.csv.gz',index=False)
    dg=[];loo=[]
    for c,ix,nA,nB,comb,D,d in spaces:
        k=K[np.ix_(ix,ix)];coef=1/nA+1/nB;aa=k[comb[:,:,None],comb[:,None,:]].sum(axis=(1,2));ac=k.sum(axis=1)[comb].sum(axis=1);bb=k.sum()-2*ac+aa
        va=(nA-aa/nA)/(nA-1);vb=(nB-bb/nB)/(nB-1);Tperm=coef*coef*aa-2*coef*ac/nB+k.sum()/nB**2
        sel=np.flatnonzero(d>0);rem=np.flatnonzero(d<0);A=k[np.ix_(sel,sel)].sum();B2=k[np.ix_(rem,rem)].sum();VA=(nA-A/nA)/(nA-1);VB=(nB-B2/nB)/(nB-1);Tobs=d@k@d;disp=abs(VA-VB)/2;student=Tobs/(VA/nA+VB/nB)
        dg.append(dict(contrast=c,variance_A_squared_hellinger=VA/2,variance_B_squared_hellinger=VB/2,variance_ratio_A_to_B=VA/VB,absolute_variance_difference=disp,dispersion_exact_reference_p=np.count_nonzero(abs(va-vb)/2>=disp-TOL)/len(D),studentized_centroid_statistic=student,studentized_exact_reference_p=np.count_nonzero(Tperm/(va/nA+vb/nB)>=student-TOL)/len(D),n_assignments=len(D)))
        for omit in ix:
            yi=ix[ix!=omit];ga,gb=c.split('_vs_');lab=groups[yi];dy=np.where(lab==ga,1/np.count_nonzero(lab==ga),-1/np.count_nonzero(lab==gb));effect=np.sqrt(max(dy@K[np.ix_(yi,yi)]@dy,0)/2)
            loo.append(dict(contrast=c,omitted_sample=samples[omit],hellinger_centroid_distance=effect))
    dg=pd.DataFrame(dg);dg['dispersion_Holm_3']=adjust(dg.dispersion_exact_reference_p,'holm');dg['studentized_Holm_3']=adjust(dg.studentized_exact_reference_p,'holm');dg.to_csv(BASE/'results/dispersion_and_studentization_diagnostics.csv',index=False)
    pd.DataFrame(loo).to_csv(BASE/'results/LOO_crosscheck.csv',index=False)
    globalci=[]
    for ci,space in enumerate(spaces):
        c,ix,_,_,_,_,d=space;obs=np.sqrt(d@K[np.ix_(ix,ix)]@d/2);q=np.quantile(global_dist[:,ci],[.025,.975]);rad=np.quantile(global_err[:,ci],.95);globalci.append(dict(contrast=c,observed_distance=obs,percentile_low=q[0],percentile_high=q[1],centered_error_radius95=rad,confidence_ball_distance_low=max(0,obs-rad),confidence_ball_distance_high=min(1,obs+rad)))
    pd.DataFrame(globalci).to_csv(BASE/'results/global_animal_bootstrap_intervals.csv',index=False)
    refs=[inp/'all_animal_junction_counts.csv.gz',inp/'animal_metadata.csv',old/'tables/frozen_features.npz',old/'tables/selected_cluster_definitions.csv',old/'tables/local_exact_results_complete.csv.gz',old/'tables/animal_alternative_compositions.csv.gz',old/'tables/leave_one_animal_out_effects.csv']
    pd.DataFrame([dict(package_relative_path=str(f.relative_to(args.base)),sha256=sha(f),bytes=f.stat().st_size) for f in refs]).to_csv(BASE/'results/frozen_input_links.csv',index=False)
    summary=dict(n_clusters=E,n_local_tests_per_family=3*E,n_original_positive_tests=len(hits),n_original_unique_clusters=int(hits.cluster_id.nunique()),n_valid_dominance95=int(passing95.sum()),n_BY_dominance95=int(l95.BY_positive.sum()),thinning_replicates=10,thinning_probability=.5,bootstrap_replicates=2000,n_original_hits_recovered95=int(st.BY_positive_dominance95.sum()),n_original_hits_recovered_all10=int((st.n_BY_positive==10).sum()),n_original_hits_recovered_any10=int((st.n_BY_positive>0).sum()),kernel_min_eigenvalue=float(np.linalg.eigvalsh(K).min()),plan_sha256=sha(BASE/'RETROSPECTIVE_PLAN_2026-10-15.md'))
    (BASE/'results/summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
