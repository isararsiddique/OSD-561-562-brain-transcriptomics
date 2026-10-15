"""Full-family read-split repeatability, animal influence and BY resolution."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[3];E=B/'extensions/effect_stability';R=E/'results'
spec=importlib.util.spec_from_file_location('robust',B/'extensions/multiway_robustness/scripts/analyse_robustness.py');rob=importlib.util.module_from_spec(spec);spec.loader.exec_module(rob)
CONTRASTS=[('FLT','GC'),('FLT','VIV'),('GC','VIV')]
def native(counts,index):
 y=counts[index['selected_rows']].astype(float);off=index['offsets'];tot=np.add.reduceat(y,off[:-1],axis=0)
 den=np.repeat(tot,np.diff(off),axis=0);p=np.divide(y,den,out=np.zeros_like(y),where=den>0).T
 bt=np.add.reduceat(counts[index['boundary_rows']],index['boundary_offsets'][:-1],axis=0)
 dom=np.divide(tot,bt,out=np.zeros_like(tot),where=bt>0)
 valid=(tot>=20).all(1)&(dom>=.9-1e-12).all(1)
 return p,valid,tot,dom

def contrasts(p,groups):return [p[groups==a].mean(0)-p[groups==b].mean(0) for a,b in CONTRASTS]
def metrics(ref,other,off):
 dot=np.add.reduceat(ref*other,off[:-1]);aa=np.add.reduceat(ref*ref,off[:-1]);bb=np.add.reduceat(other*other,off[:-1])
 norm=np.sqrt(aa*bb);cos=np.divide(dot,norm,out=np.full_like(dot,np.nan),where=norm>1e-20)
 proj=np.divide(dot,aa,out=np.full_like(dot,np.nan),where=aa>1e-20)
 return dict(cosine=cos,relative_projection=proj,native_total_variation_pp=50*np.add.reduceat(abs(other),off[:-1]),reference_total_variation_pp=50*np.add.reduceat(abs(ref),off[:-1]),orthogonal_residual_norm=np.sqrt(np.maximum(bb-np.divide(dot*dot,aa,out=np.zeros_like(dot),where=aa>1e-20),0)))

def envelope(frame,family):
 p=np.sort(frame.exact_p.to_numpy());m=len(p);hm=np.sum(1/np.arange(1,m+1));rat=np.arange(1,m+1)*.05/(m*hm*p);j=int(np.argmax(rat));q=rob.adjust(frame.exact_p,'by')
 return dict(family=family,n_tests=m,harmonic_penalty=hm,minimum_exact_p=p[0],minimum_BY_q=q.min(),maximum_BY_support_ratio=rat[j],support_rank=j+1,support_raw_p=p[j],BY_positive=int((q<=.05).sum()),n_at_minimum_p=int(np.sum(abs(p-p[0])<1e-14)),minimum_tests_at_observed_minimum_for_BY=int(np.ceil(m*hm*p[0]/.05)),valid_tests=int(frame.measurement_valid.sum()) if 'measurement_valid' in frame else m)

def main():
 assert (E/'RETROSPECTIVE_PLAN_2026-10-15.md').exists()
 idx=np.load(B/'extensions/multiway_robustness/results/frozen_count_index.npz');off=idx['offsets'];ids=idx['cluster_ids'];samples=idx['samples'];groups=pd.read_csv(B/'inputs/animal_metadata.csv').set_index('sample').loc[samples,'group'].to_numpy()
 counts=pd.read_csv(B/'inputs/all_animal_junction_counts.csv.gz',index_col=0);assert np.array_equal(counts.columns,samples);x=counts.to_numpy(np.int64)
 p,valid,tot,dom=native(x,idx);assert valid.all();orig_delta=contrasts(p,groups);old=pd.read_csv(B/'extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz');hits=set(map(tuple,old.loc[old.BY_q_full_3E<=.05,['contrast','cluster_id']].to_numpy()))
 names=[a+'_vs_'+b for a,b in CONTRASTS];spaces=rob.labels(groups);effect=[];paired=[];complocal=[];compglobal=[];gates=[];envelopes=[envelope(old,'original')]
 oldthin=pd.read_csv(B/'extensions/multiway_robustness/results/local_thinning_complete_fixed_families.csv.gz')
 for rep in range(1,11):
  seed=2026101500+rep;draw=np.load(B/f'extensions/multiway_robustness/thinned_counts/replicate_{rep:02d}_seed_{seed}.npz')['counts'].astype(np.int64)
  complement=x-draw;assert complement.min()>=0
  pa,va,ta,da=native(draw,idx);pb,vb,tb,db=native(complement,idx);d_a=contrasts(pa,groups);d_b=contrasts(pb,groups)
  for half,ps,vs,ts,ds,dh in [('A',pa,va,ta,da,d_a),('B',pb,vb,tb,db,d_b)]:
   for ci,c in enumerate(names):
    m=metrics(orig_delta[ci],dh[ci],off)
    for k,cid in enumerate(ids):effect.append(dict(replicate=rep,half=half,contrast=c,cluster_id=cid,measurement_valid=bool(vs[k]),primary_BY_test=(c,cid) in hits,**{key:float(v[k]) for key,v in m.items()}))
   for k,cid in enumerate(ids):gates.append(dict(replicate=rep,half=half,cluster_id=cid,valid=bool(vs[k]),minimum_selected_reads=int(ts[k].min()),minimum_dominance=float(ds[k].min())))
  for ci,c in enumerate(names):
   m=metrics(d_a[ci],d_b[ci],off)
   for k,cid in enumerate(ids):paired.append(dict(replicate=rep,contrast=c,cluster_id=cid,both_halves_measured=bool(va[k]&vb[k]),primary_BY_test=(c,cid) in hits,**{key:float(v[k]) for key,v in m.items()}))
  lb,gb=rob.exact_family(np.sqrt(pb),off,vb,spaces,ids,rep);complocal.append(lb);compglobal.append(gb)
  envelopes.extend([envelope(oldthin[oldthin.replicate==rep],f'split_{rep:02d}_A'),envelope(lb,f'split_{rep:02d}_B')])
  print('split',rep,'complement valid',vb.sum(),'BY',lb.BY_positive.sum(),flush=True)
 ef=pd.DataFrame(effect);pair=pd.DataFrame(paired);ef.to_csv(R/'all_frozen_effect_vectors_split_metrics.csv.gz',index=False);pair.to_csv(R/'all_complementary_half_concordance.csv.gz',index=False);pd.DataFrame(gates).to_csv(R/'all_read_split_measurement_gates.csv.gz',index=False)
 pd.concat(complocal,ignore_index=True).to_csv(R/'complement_B_local_exact_full_families.csv.gz',index=False);pd.concat(compglobal,ignore_index=True).to_csv(R/'complement_B_global_exact.csv',index=False);pd.DataFrame(envelopes).to_csv(R/'BY_resolution_envelopes.csv',index=False)
 influence=[]
 for ci,(a,b) in enumerate(CONTRASTS):
  c=names[ci]
  for animal in np.flatnonzero(np.isin(groups,[a,b])):
   ix=np.arange(len(groups))!=animal;delta=p[ix&(groups==a)].mean(0)-p[ix&(groups==b)].mean(0);m=metrics(orig_delta[ci],delta,off)
   for k,cid in enumerate(ids):influence.append(dict(contrast=c,omitted_sample=samples[animal],omitted_group=groups[animal],cluster_id=cid,primary_BY_test=(c,cid) in hits,**{key:float(v[k]) for key,v in m.items()}))
 influence=pd.DataFrame(influence);influence.to_csv(R/'all_animal_leave_one_out_native_effects.csv.gz',index=False)
 summaries=[]
 for label,df,validcol in [('original_to_half',ef,'measurement_valid'),('half_to_half',pair,'both_halves_measured'),('leave_one_animal_out',influence,None)]:
  for c in names:
   for hit in [False,True]:
    f=df[(df.contrast==c)&(df.primary_BY_test if hit else np.ones(len(df),bool))].copy()
    if f.empty:continue
    n=len(f)
    if validcol:f=f[f[validcol]]
    if f.empty:continue
    summaries.append(dict(analysis=label,contrast=c,subset='original_BY_tests' if hit else 'all923',instances_total=n,instances_measured=len(f),positive_alignment_fraction=(f.cosine>0).mean(),median_cosine=f.cosine.median(),cosine_q05=f.cosine.quantile(.05),cosine_q95=f.cosine.quantile(.95),median_relative_projection=f.relative_projection.median(),projection_q05=f.relative_projection.quantile(.05),projection_q95=f.relative_projection.quantile(.95),median_TV_pp=f.native_total_variation_pp.median()))
 sm=pd.DataFrame(summaries);sm.to_csv(R/'effect_stability_summary.csv',index=False)
 tc='6:2:42663893:start';ef[ef.cluster_id==tc].to_csv(R/'Tcaf1_all20_halves_native_effects.csv',index=False);pair[pair.cluster_id==tc].to_csv(R/'Tcaf1_all10_complementary_pairs.csv',index=False);influence[influence.cluster_id==tc].to_csv(R/'Tcaf1_all_animal_omissions.csv',index=False)
 # Full-family ranked curves, including p=1 gates, are diagnostic only.
 rank=[]
 families=[('original',old)]+[(f'split_{r:02d}_A',oldthin[oldthin.replicate==r]) for r in range(1,11)]+[(f'split_{r:02d}_B',d) for r,d in enumerate(complocal,1)]
 for name,df in families:
  ps=np.sort(df.exact_p);m=len(ps);critical=np.arange(1,m+1)*.05/(m*np.sum(1/np.arange(1,m+1)))
  rank.extend(dict(family=name,rank=i+1,p_value=float(pv),BY_critical=float(cv),passes_step_up_line=bool(pv<=cv)) for i,(pv,cv) in enumerate(zip(ps,critical)))
 pd.DataFrame(rank).to_csv(R/'complete_ordered_pvalue_curves.csv.gz',index=False)
 summary={'n_clusters':len(ids),'full_local_family':len(old),'n_original_BY_tests':len(hits),'paired_read_splits':10,'complement_B_BY_calls':int(sum(d.BY_positive.sum() for d in complocal)),'complement_B_Holm_global_positive':int(sum((d.Holm_p_3<=.05).sum() for d in compglobal)),'effect_metric_rows':len(ef),'paired_metric_rows':len(pair),'animal_influence_rows':len(influence),'counts_are_disjoint_halves_of_same_animals':True,'new_biological_replication':False,'plan_sha256':hashlib.sha256((E/'RETROSPECTIVE_PLAN_2026-10-15.md').read_bytes()).hexdigest()}
 (R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(sm.to_string(index=False))
if __name__=='__main__':main()
