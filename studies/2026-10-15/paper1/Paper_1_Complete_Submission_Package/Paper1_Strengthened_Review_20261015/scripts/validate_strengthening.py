#!/usr/bin/env python3
"""Recompute full extension families, identity accounting and comparator summaries."""
import json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

def validate(root=None):
 root=Path(root or Path(__file__).resolve().parents[1]);x=root/'extensions/strengthening_20261015';checks=[]
 def ck(name,value):checks.append(dict(name=name,passed=bool(value)))
 def read(p):return pd.read_csv(x/p)
 def js(p):return json.loads((x/p).read_text())
 def equal(a,b):return np.allclose(a,b,rtol=1e-9,atol=1e-11,equal_nan=True)
 md=read('external_data/harmonized_metadata.csv');primary=md[md.primary_analysis_included]
 ck('66 acquired units,64 primary units',len(md)==66 and len(primary)==64)
 ck('Primary accession units40,6,8,10',primary.groupby('accession').size().to_dict()=={'OSD-202':40,'OSD-352':6,'OSD-374':8,'OSD-611':10})
 ck('Primary units46animals18pools',primary.replicate_unit.eq('animal').sum()==46 and primary.replicate_unit.ne('animal').sum()==18)
 ck('No duplicate acquired expression unit',not md.duplicated(['accession','independent_unit_id']).any())
 for a,q in md.groupby('accession'):
  c=pd.read_csv(x/'external_data'/a/'counts.csv.gz',index_col=0)
  arr=c.to_numpy();ck(a+' count columns match metadata',set(c.columns)==set(q.sample_id));ck(a+' integer nonnegative counts',np.isfinite(arr).all() and (arr>=0).all() and np.equal(arr,np.floor(arr)).all())
 for r in js('external_data/output_manifest.json'):
  # Original acquisition paths are historical provenance; resolve archived copies.
  rel=r['path'].split('/external_data/',1)[1];p=x/'external_data'/rel
  ck('Acquisition source hash '+rel,p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'])
 d=read('method_benchmark/strengthened_outer_results.csv');s=read('method_benchmark/strengthened_summary.csv');roles=read('method_benchmark/strengthened_split_roles.csv');inner=read('method_benchmark/strengthened_inner_CV.csv')
 ck('1002new outer rows retained',len(d)==1002 and not d.duplicated(['accession','split_id','method']).any())
 for (a,sp),q in roles.groupby(['accession','split_id']):
  train=set(q[q.role.eq('train')].unit);test=set(q[q.role.eq('test')].unit);z=d[(d.accession==a)&(d.split_id==sp)]
  ck(a+f' split{sp} disjoint complete roles',not train&test and len(train)>0 and len(test)>0 and z.n_train.eq(len(train)).all() and z.n_test.eq(len(test)).all())
  f=inner[(inner.accession==a)&(inner.split_id==sp)]
  ck(a+f' split{sp} one minimum inner choice',len(f)==5 and f.selected.sum()==1 and equal(f[f.selected].inner_SSE.iloc[0],f.inner_SSE.min()))
 ck('64external uniqueheldout units',d[d.scope.eq('independent_exposure_cohort')][['accession','split_id']].drop_duplicates().shape[0]==64)
 ck('Every MSEratio directly agreeswith SSE',equal(d.MSE_ratio,d.prediction_SSE/d.baseline_SSE))
 for row in s.itertuples():
  q=d[(d.accession==row.accession)&(d.method==row.method)]
  ck(row.accession+'/'+row.method+' fullsummarymean',len(q)==row.n_splits and equal(q.MSE_ratio.mean(),row.mean_MSE_ratio) and equal(q.prediction_SSE.sum()/q.baseline_SSE.sum(),row.pooled_SSE_ratio))
 cert=read('method_benchmark/convex_optimality_certificates.csv')
 ck('16084certified finite convex fits',len(cert)==16084 and cert.converged.all() and np.isfinite(cert.relative_gap).all())
 ck('All convex optimality tolerances satisfied',cert.relative_gap.le(1.00001e-7).all() and cert.coefficient_spectral_norm.le(np.sqrt(3)+1e-9).all())
 ck('All new constrained outer gains satisfy bound',d[d.method.eq('spectral_constrained_ridge_CV')].operator_gain.le(2+1e-9).all())
 verify=js('method_benchmark/convex_independent_verification.json');ck('96independentfullcoefficient checks pass',verify['status']=='PASS' and verify['n_checks']==96 and all(q['passed'] for q in verify['checks']))
 matched=read('method_benchmark/matched_grid_outer_results.csv');summary=read('method_benchmark/matched_grid_summary.csv');obj=read('method_benchmark/matched_grid_fixed_alpha_objectives.csv')
 ck('852matchedgridrows426splitspreserved',len(matched)==852 and matched[['accession','split_id']].drop_duplicates().shape[0]==426)
 for r in summary.itertuples():
  q=matched[(matched.accession==r.accession)&(matched.method==r.method)]
  ck(r.accession+'/'+r.method+' matchedsummary',equal(q.MSE_ratio.mean(),r.mean_MSE_ratio) and equal(q.prediction_SSE.sum()/q.baseline_SSE.sum(),r.pooled_SSE_ratio) and equal(q.convex_minus_comparator_MSE_ratio.mean(),r.paired_convex_minus_comparator_mean))
 ck('28fixedalphafullobjectivechecks',len(obj)==28 and obj.verified_dominance_with_gap.all() and obj.converged.all() and obj.convex_objective_reconstruction_error.max()<1e-9)
 B=read('correction_benchmark/all_repetitions.csv');bs=read('correction_benchmark/benchmark_summary.csv');pr=js('correction_benchmark/protocol.json')
 ck('1400independentcountdatasets14methodseach',len(B)==19600 and B.groupby(['scenario','n_train','repetition']).size().eq(14).all() and B[['scenario','n_train','repetition']].drop_duplicates().shape[0]==1400 and set(B.method)==set(pr['methods']) and set(B.scenario)==set(pr['scenarios']))
 ck('All19600methodstatusesexplicitandPASS',B.status.eq('PASS').all())
 ck('600realcountcalibrationgenesretained',len(read('correction_benchmark/real_count_calibration.csv'))==600)
 ck('Simulatedconstrainedoperatorgainsbounded',B[B.method.eq('ground_constrained_ridge_CV')].operator_gain.le(2+1e-8).all())
 metrics=['effect_RMSE','outside_effect_RMSE','biological_FDP','power','null_gene_nominal_rejection']
 for r in bs.to_dict('records'):
  q=B[(B.scenario==r['scenario'])&(B.n_train==r['n_train'])&(B.method==r['method'])]
  ok=len(q)==100
  for m in metrics:
   v=q[m].dropna();ok=ok and equal(v.mean(),r[m+'_mean']) and equal(v.std(ddof=1)/np.sqrt(len(v)) if len(v) else np.nan,r[m+'_MCse'])
  ck(r['scenario']+'/'+str(r['n_train'])+'/'+r['method']+' completeMCmeans',ok)
 provenance=js('correction_benchmark/comparator_provenance.json')
 for f,key in [('neuroCombat_0_2_12_source.py','implementation_sha256'),('RUVg_primary_source.R','RUVg_source_sha256')]:ck('Comparator sourcehash '+f,hashlib.sha256((x/'correction_benchmark'/f).read_bytes()).hexdigest()==provenance[key])
 b=read('biological_validation/all_biological_effects.csv');f=b.p_wild.notna()
 ck('1488biologicalrows1047finiteprobabilities',len(b)==1488 and f.sum()==1047)
 for m,col in [('bh','q_BH_complete_new_family'),('by','q_BY_complete_new_family')]:ck('Independentrecomputecompletebiological '+m,equal(stats.false_discovery_control(b.loc[f,'p_wild'],method=m),b.loc[f,col]))
 bv=js('biological_validation/independent_verification.json');ck('1337independentbiologicalcheckspass',bv['passed'] and bv['n_checks']==1337 and not bv['failed'])
 for rel,dig in js('biological_validation/source_hashes.json').items():ck('Biologicalsourcehash '+rel,hashlib.sha256((root/rel).read_bytes()).hexdigest()==dig)
 (root/'results/strengthening_validation.json').write_text(json.dumps(dict(passed=all(c['passed'] for c in checks),n_checks=len(checks),failed=[c for c in checks if not c['passed']],checks=checks),indent=2)+'\n')
 return checks
if __name__=='__main__':
 import sys
 c=validate();failed=[x for x in c if not x['passed']];print(json.dumps(dict(n_checks=len(c),failed=failed),indent=2));sys.exit(bool(failed))
