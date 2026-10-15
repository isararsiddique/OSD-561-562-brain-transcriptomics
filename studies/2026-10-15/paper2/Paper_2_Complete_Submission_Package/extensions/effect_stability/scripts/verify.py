"""Independent count, metric, correction and family checks (no analysis imports)."""
from pathlib import Path
import json,itertools
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[3];E=B/'extensions/effect_stability';R=E/'results'
checks=[]
def check(name,c):
 assert c,name
 checks.append(name)
def by(p):
 p=np.asarray(p);m=len(p);o=np.argsort(p);v=p[o]*m*sum(1/j for j in range(1,m+1))/np.arange(1,m+1);v=np.minimum.accumulate(v[::-1])[::-1];z=np.empty(m);z[o]=np.minimum(v,1);return z

def main():
 ef=pd.read_csv(R/'all_frozen_effect_vectors_split_metrics.csv.gz');pr=pd.read_csv(R/'all_complementary_half_concordance.csv.gz');lo=pd.read_csv(R/'all_animal_leave_one_out_native_effects.csv.gz');loc=pd.read_csv(R/'complement_B_local_exact_full_families.csv.gz');env=pd.read_csv(R/'BY_resolution_envelopes.csv')
 check('all effect, pair, influence and complementary families complete',len(ef)==55380 and len(pr)==27690 and len(lo)==47996 and len(loc)==27690 and len(env)==21)
 check('metrics finite and bounded cosines',np.isfinite(ef.cosine).all() and ef.cosine.between(-1-1e-12,1+1e-12).all())
 old=pd.read_csv(B/'extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz');hit=set(map(tuple,old.loc[old.BY_q_full_3E<=.05,['contrast','cluster_id']].to_numpy()))
 check('posthoc subset exactlyoriginal43 tests',set(map(tuple,ef.loc[ef.primary_BY_test,['contrast','cluster_id']].drop_duplicates().to_numpy()))==hit)
 for rep,df in loc.groupby('replicate'):
  check(f'complement{rep} full BY independently reconstructed',len(df)==2769 and np.allclose(by(df.exact_p),df.BY_q_full_2769,atol=1e-12))
  check(f'complement{rep} gate failures retained',np.all(df.loc[~df.measurement_valid,'exact_p']==1))
  e=env[env.family==f'split_{rep:02d}_B'].iloc[0];ps=np.sort(df.exact_p);harm=sum(1/j for j in range(1,len(ps)+1));ratio=np.max(np.arange(1,len(ps)+1)*.05/(len(ps)*harm*ps))
  check(f'complement{rep} BY crossing equivalence',np.isclose(ratio,e.maximum_BY_support_ratio) and bool(ratio>=1)==bool((by(df.exact_p)<=.05).any()))
 # Direct reconstruction by coordinate membership, independent of frozen_count_index.
 x=pd.read_csv(B/'inputs/all_animal_junction_counts.csv.gz',index_col=0);m=pd.read_csv(B/'inputs/animal_metadata.csv').set_index('sample').loc[x.columns];groups=m.group.to_numpy()
 long=pd.read_csv(B/'extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz');member=long[['cluster_id','alternative_index','junction_id']].drop_duplicates().sort_values(['cluster_id','alternative_index'])
 cids=sorted(member.cluster_id.unique());testids=[cids[0],cids[len(cids)//2],cids[-1],'6:2:42663893:start'];lookup={v:i for i,v in enumerate(x.index)}
 for rep in range(1,11):
  draw=np.load(B/f'extensions/multiway_robustness/thinned_counts/replicate_{rep:02d}_seed_{2026101500+rep}.npz')['counts'].astype(np.int64);com=x.to_numpy()-draw
  check(f'split{rep} all counts partition exactly',np.all(com>=0) and np.array_equal(com+draw,x.to_numpy()))
  for cid in testids:
   rows=[lookup[j] for j in member[member.cluster_id==cid].junction_id];raw=x.to_numpy()[rows].astype(float);aa=draw[rows].astype(float);bb=com[rows].astype(float);ps=[a/a.sum(axis=0) for a in [raw,aa,bb]]
   for g1,g2 in [('FLT','GC'),('FLT','VIV'),('GC','VIV')]:
    c=g1+'_vs_'+g2;ds=[p[:,groups==g1].mean(1)-p[:,groups==g2].mean(1) for p in ps];target=pr[(pr.replicate==rep)&(pr.cluster_id==cid)&(pr.contrast==c)].iloc[0]
    cosine=float(np.dot(ds[1],ds[2])/(np.linalg.norm(ds[1])*np.linalg.norm(ds[2])))
    check(f'split{rep} {cid} {c} paired cosine',abs(cosine-target.cosine)<1e-10)
    for half,d in zip(['A','B'],ds[1:]):
     rec=ef[(ef.replicate==rep)&(ef.cluster_id==cid)&(ef.contrast==c)&(ef.half==half)].iloc[0];proj=float(np.dot(ds[0],d)/np.dot(ds[0],ds[0]));check(f'split{rep} {cid} {c}{half} projection',abs(proj-rec.relative_projection)<1e-10)
 # Exact full-vector permutation result for two complementary TcAF1 contrasts.
 rows=[lookup[j] for j in member[member.cluster_id=='6:2:42663893:start'].junction_id];p=com[rows].astype(float);p/=p.sum(axis=0)
 for a,b in [('FLT','GC'),('FLT','VIV')]:
  ix=np.flatnonzero(np.isin(groups,[a,b]));n=(groups[ix]==a).sum();h=np.sqrt(p[:,ix]).T;obs=np.square(h[groups[ix]==a].mean(0)-h[groups[ix]==b].mean(0)).sum();ge=0;alloc=0
  for combo in itertools.combinations(range(len(ix)),n):
   sel=np.zeros(len(ix),bool);sel[list(combo)]=True;t=np.square(h[sel].mean(0)-h[~sel].mean(0)).sum();ge+=t>=obs-1e-12;alloc+=1
  rr=loc[(loc.replicate==10)&(loc.cluster_id=='6:2:42663893:start')&(loc.contrast==a+'_vs_'+b)].iloc[0]
  check('independent exact Tcaf1 complement10 '+a+b,abs(ge/alloc-rr.exact_p)<1e-12)
 result={'status':'PASS','checks_passed':len(checks),'checks':checks,'scope':'Independent numerical checks, not biological replication.'};(R/'verification.json').write_text(json.dumps(result,indent=2));print('PASS',len(checks),'independent checks')
if __name__=='__main__':main()
