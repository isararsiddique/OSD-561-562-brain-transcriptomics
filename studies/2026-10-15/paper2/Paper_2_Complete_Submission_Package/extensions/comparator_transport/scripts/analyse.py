"""Frozen comparator-specific transport, including distinct-animal RR10 cohort."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import importlib.util,itertools,json,hashlib
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[3];E=B/'extensions/comparator_transport';R=E/'results';OLD=B/'extensions/multiway_transport'
spec=importlib.util.spec_from_file_location('oldtransport',OLD/'scripts/analyse.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old);old.NBOOT=4000
TOL=1e-12

def weights(meta):
 observed=np.zeros(len(meta));allocations=np.zeros((1,len(meta)));design=[]
 for st in sorted(meta.stratum.unique()):
  ix=np.flatnonzero(meta.stratum.to_numpy()==st);nf=int(meta.iloc[ix].exposure.eq('FLT').sum());ng=len(ix)-nf;w=len(ix)/len(meta);assert nf>=2 and ng>=2
  combos=list(itertools.combinations(ix,nf));D=np.zeros((len(combos),len(meta)));D[:,ix]=-w/ng
  for j,co in enumerate(combos):D[j,list(co)]=w/nf
  observed[ix]=np.where(meta.iloc[ix].exposure.eq('FLT'),w/nf,-w/ng)
  allocations=(allocations[:,None,:]+D[None,:,:]).reshape(-1,len(meta));design.append(dict(stratum=st,n_FLT=nf,n_GC=ng,weight=w,n_assignments=len(combos)))
 return observed,allocations,pd.DataFrame(design)

def main():
 assert (E/'RETROSPECTIVE_PLAN_2026-10-15.md').exists()
 defs=pd.read_csv(OLD/'inputs/RR23_frozen_cluster_definitions.csv',dtype={'chrom':str});defs['boundary_key']=defs.chrom+':'+defs.strand.astype(str)+':'+defs.boundary.astype(str)
 comp=pd.read_csv(OLD/'inputs/RR23_frozen_alternative_compositions.csv.gz');members=comp[['cluster_id','alternative_index','junction_id','matching_gene_ids']].drop_duplicates().sort_values(['cluster_id','alternative_index']).reset_index(drop=True)
 assert list(members.cluster_id.drop_duplicates())==list(defs.cluster_id)
 off=np.r_[0,defs.n_alternatives.cumsum()];sizes=np.diff(off);ids=defs.cluster_id.to_numpy();rr=old.reconstruct(OLD/'inputs/RR23_animal_junction_counts.csv.gz',defs,members,805562)
 rmeta=pd.read_csv(OLD/'inputs/RR23_animal_metadata.csv').set_index('sample',drop=False).loc[rr['samples']]
 cohorts={};inputfiles=[OLD/'inputs/RR23_animal_junction_counts.csv.gz',OLD/'inputs/RR23_animal_metadata.csv',OLD/'inputs/RR23_frozen_cluster_definitions.csv',OLD/'inputs/RR23_frozen_alternative_compositions.csv.gz']
 for tag,path,rows,mp in [('RRRM2',OLD/'inputs/RRRM2_animal_junction_counts.csv.gz',668237,OLD/'inputs/RRRM2_animal_metadata.csv'),('RR10',B/'extensions/cohort_search/inputs/RR10/RR10_animal_junction_counts.csv.gz',590325,B/'extensions/cohort_search/inputs/RR10/RR10_animal_metadata.csv')]:
  val=old.reconstruct(path,defs,members,rows);meta=pd.read_csv(mp).rename(columns={'index':'sample'}).set_index('sample',drop=False).loc[val['samples']].copy()
  if tag=='RR10':meta['exposure']=meta.group;meta['stratum']='all'
  gate=(val['selected_counts'].min(1)>=20)&(val['dominance'].min(1)>=.9-TOL);assert gate.any()
  observed,alloc,design=weights(meta);assert len(alloc)==(252 if tag=='RR10' else 11200)
  design.to_csv(R/f'{tag}_allocation_design.csv',index=False);meta.to_csv(R/f'{tag}_included_animals.csv',index=False)
  gt=defs[['cluster_id']].copy();gt['minimum_selected_reads']=val['selected_counts'].min(1);gt['minimum_dominance']=val['dominance'].min(1);gt['common_measurable']=gate;gt.to_csv(R/f'{tag}_all923_measurement_gates.csv',index=False)
  np.savez_compressed(R/f'{tag}_frozen_measurements.npz',counts=val['counts'],proportions=val['fractions'],selected=val['selected_counts'],boundary=val['all_boundary_counts'],gate=gate,samples=np.array(val['samples']),offsets=off,cluster_ids=ids)
  cohorts[tag]=dict(val=val,meta=meta,gate=gate,observed=observed,alloc=alloc)
  inputfiles.extend([path,mp]);print(tag,'measurable',gate.sum(),'animals',len(meta),flush=True)
 intersection=cohorts['RRRM2']['gate']&cohorts['RR10']['gate'];pd.DataFrame(dict(cluster_id=ids,in_both_validation_cohorts=intersection)).to_csv(R/'common_both_validation_panel.csv',index=False)
 globals=[];locals=[];animal_scores=[];boots=[];nulls=[];train=[]
 for ci,(tag,co) in enumerate(cohorts.items()):
  val,meta,gate,observed,alloc=[co[k] for k in ['val','meta','gate','observed','alloc']];rng=np.random.default_rng(2026101541+ci);W=old.bootstrap_weights(meta,rng)
  for compgroup in ['GC','VIV']:
   rflt=rmeta.group.eq('FLT').to_numpy();rcontrol=rmeta.group.eq(compgroup).to_numpy();base=rr['fractions'][:,rcontrol].mean(1);d=rr['fractions'][:,rflt].mean(1)-base;norm=np.sqrt(np.add.reduceat(d*d,off[:-1]));u=np.divide(d,np.repeat(norm,sizes),out=np.zeros_like(d),where=np.repeat(norm,sizes)>TOL)
   direction='FLT_minus_'+compgroup
   if ci==0:
    tr=members.copy();tr['direction']=direction;tr['RR23_native_delta']=d;tr['RR23_unit_vector']=u;tr['RR23_comparator_mean']=base;train.append(tr)
   scores=np.add.reduceat(u[:,None]*(val['fractions']-base[:,None]),off[:-1],axis=0);scores[~gate]=np.nan
   local_effect=np.full(923,np.nan);local_effect[gate]=scores[gate]@observed
   lp=np.ones(923);lp2=np.ones(923)
   for j in range(0,923,40):
    k=np.arange(j,min(j+40,923));validk=k[gate[k]]
    if not len(validk):continue
    ns=scores[validk]@alloc.T;obs=local_effect[validk,None];lp[validk]=(ns>=obs-TOL).mean(1);lp2[validk]=(abs(ns)>=abs(obs)-TOL).mean(1)
   lc=pd.DataFrame(dict(cohort=tag,direction=direction,cluster_id=ids,measurement_valid=gate,effect=local_effect,exact_p=lp,exact_two_sided_p=lp2,n_assignments=len(alloc)));locals.append(lc)
   # Resample full source-animal vectors; discovery directions never use validation outcomes.
   trmeta=rmeta.loc[rmeta.group.isin(['FLT',compgroup])].copy();trix=rmeta.index.get_indexer(trmeta.index);trmeta['exposure']=np.where(trmeta.group.eq('FLT'),'FLT','GC');trmeta['stratum']='all';WT=old.bootstrap_weights(trmeta,rng)
   np.savez_compressed(R/f'{tag}_{direction}_bootstrap_weights.npz',validation_weights=W,training_weights=WT,training_indices=trix,validation_samples=np.array(val['samples']),training_samples=trmeta.index.to_numpy(dtype=str),seed=np.array(2026101541+ci))
   for panel,chosen in [('cohort_measured',gate),('intersection_both_cohorts',intersection)]:
    if not chosen.any():continue
    sc=scores[chosen].mean(0);effect=float(sc@observed);nul=alloc@sc;fixed=W@sc
    aa=np.repeat(chosen,sizes);td=WT@rr['fractions'][aa][:,trix].T;vd=W@val['fractions'][aa].T;osub=np.r_[0,sizes[chosen].cumsum()]
    tn=np.sqrt(np.add.reduceat(td*td,osub[:-1],axis=1));joint=np.divide(np.add.reduceat(td*vd,osub[:-1],axis=1),tn,out=np.zeros_like(tn),where=tn>TOL).mean(1)
    gl=dict(cohort=tag,direction=direction,panel=panel,n_clusters=int(chosen.sum()),n_animals=len(meta),n_FLT=int(meta.exposure.eq('FLT').sum()),n_GC=int(meta.exposure.eq('GC').sum()),n_assignments=len(alloc),effect_pp=100*effect,conditional_CI95_low_pp=100*np.quantile(fixed,.025),conditional_CI95_high_pp=100*np.quantile(fixed,.975),two_cohort_CI95_low_pp=100*np.quantile(joint,.025),two_cohort_CI95_high_pp=100*np.quantile(joint,.975),exact_p=float((nul>=effect-TOL).mean()),exact_two_sided_p=float((abs(nul)>=abs(effect)-TOL).mean()))
    globals.append(gl);print(gl,flush=True)
    for j,sample in enumerate(val['samples']):animal_scores.append(dict(cohort=tag,direction=direction,panel=panel,sample=sample,exposure=meta.iloc[j].exposure,stratum=meta.iloc[j].stratum,score=float(sc[j]),contrast_weight=float(observed[j])))
    boots.extend(dict(cohort=tag,direction=direction,panel=panel,draw=j,conditional_effect_pp=100*f,two_cohort_effect_pp=100*t) for j,(f,t) in enumerate(zip(fixed,joint)))
    nulls.extend(dict(cohort=tag,direction=direction,panel=panel,assignment=j,effect_pp=100*v) for j,v in enumerate(nul))
 globals=pd.DataFrame(globals)
 for panel in globals.panel.unique():
  mask=globals.panel.eq(panel);globals.loc[mask,'Holm_four_global_p']=old.adjust(globals.loc[mask,'exact_p'],'holm')
 local=pd.concat(locals,ignore_index=True)
 for meth in ['bh','by','holm']:local[meth.upper()+'_joint3692']=old.adjust(local.exact_p,meth)
 local.to_csv(R/'all3692_local_directional_tests.csv.gz',index=False);globals.to_csv(R/'all_comparator_global_transport.csv',index=False);pd.DataFrame(animal_scores).to_csv(R/'all_animal_scores.csv',index=False);pd.DataFrame(boots).to_csv(R/'all_bootstrap_draws.csv.gz',index=False);pd.DataFrame(nulls).to_csv(R/'all_exact_global_nulls.csv.gz',index=False);pd.concat(train).to_csv(R/'RR23_only_frozen_directions.csv.gz',index=False)
 local[local.cluster_id.eq('6:2:42663893:start')].to_csv(R/'Tcaf1_all_comparators_and_cohorts.csv',index=False)
 pd.DataFrame([dict(file=str(p.relative_to(B)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputfiles]).to_csv(R/'frozen_input_hashes.csv',index=False)
 summary=dict(RRRM2_animals=19,RR10_animals=10,RRRM2_common_clusters=int(cohorts['RRRM2']['gate'].sum()),RR10_common_clusters=int(cohorts['RR10']['gate'].sum()),common_both_validation=int(intersection.sum()),local_family=3692,primary_global_Holm_calls=int(((globals.panel=='cohort_measured')&(globals.Holm_four_global_p<=.05)).sum()),local_BH_calls=int((local.BH_joint3692<=.05).sum()),local_BY_calls=int((local.BY_joint3692<=.05).sum()),RR10_shares_discovery_launch=True,RR10_distinct_source_animals=True,RR10_exposure_date_alias=True,retrospective_plan_sha256=hashlib.sha256((E/'RETROSPECTIVE_PLAN_2026-10-15.md').read_bytes()).hexdigest())
 (R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
