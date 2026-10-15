#!/usr/bin/env python3
"""Independent explicit regression verification of the efficient exact-null code."""
from pathlib import Path
import itertools, json
import numpy as np
import pandas as pd
import run_independent_audit as r


def main():
    rows=[]
    for n in ('561','562'):
        _,c,m,_=r.read_source(n);m=m.loc[m.flight.eq('GC')];c=c.loc[:,m.index];age=m.age.to_numpy()
        null=pd.read_csv(r.OUT/f'OSD-{n}_exact_label_null.csv',dtype={'labels':str});names=list(r.PANELS)+['WholeGenome']
        for row_id in [0,123,399]:
            labels=np.array(list(null.loc[row_id,'labels']),int)
            cells=[np.flatnonzero((labels==g)&(age==a)) for g,a in itertools.product([0,1],sorted(set(age)))];sums=np.zeros(len(names));ns=0
            for chosen in itertools.product(*cells):
                test=np.array(chosen);train=np.setdiff1d(np.arange(len(m)),test)
                gs,yt,ye,_,_=r.expression_reference(c.to_numpy(float),c.index.to_numpy(),train,test);panels=r.panel_indices(gs)
                x=np.column_stack([np.ones(len(train)),age[train]=='29',labels[train]])
                betas=np.linalg.lstsq(x,yt.T,rcond=None)[0][-1]
                for j,ix in enumerate(list(panels.values())+[np.arange(len(gs))]):
                    weight=betas[ix]/np.linalg.norm(betas[ix]);scores=weight@ye[ix]
                    for a in sorted(set(age)):
                        pos=np.flatnonzero((age[test]==a)&(labels[test]==1))[0];neg=np.flatnonzero((age[test]==a)&(labels[test]==0))[0]
                        sums[j]+=(scores[pos]>scores[neg])+0.5*(scores[pos]==scores[neg])
                ns+=1
            assert ns==81
            maxerr=float(np.max(abs(sums/(2*ns)-null.loc[row_id,names].to_numpy(float))))
            rows.append(dict(accession=f'OSD-{n}',assignment_id=row_id,check='full_regression_81_splits_vs_cached_Gram_null',max_absolute_difference=maxerr,passed=maxerr<1e-12))
        vals=null[names].to_numpy();lookup={label:i for i,label in enumerate(null.labels)}
        maxsym=max(float(np.max(abs(vals[i]-vals[lookup[''.join('1' if g=='0' else '0' for g in label)]]))) for i,label in enumerate(null.labels))
        rows.append(dict(accession=f'OSD-{n}',assignment_id=-1,check='all_400_label_complement_invariance',max_absolute_difference=maxsym,passed=maxsym<1e-12))
        rows.append(dict(accession=f'OSD-{n}',assignment_id=-1,check='conditional_null_mean_AUC_half_all_9_scores',max_absolute_difference=float(np.max(abs(vals.mean(axis=0)-.5))),passed=bool(np.allclose(vals.mean(axis=0),.5,atol=1e-12))))
    rows=pd.DataFrame(rows);rows.to_csv(r.OUT/'exact_ground_independent_regression_checks.csv',index=False)
    (r.OUT/'exact_ground_verification_summary.json').write_text(json.dumps(dict(status='PASS' if rows.passed.all() else 'FAIL',n_checks=len(rows),n_explicit_refitted_splits=486,maximum_error=float(rows.max_absolute_difference.max())),indent=2)+'\n')
    print(rows.to_string(index=False))
    assert rows.passed.all()

if __name__=='__main__':main()
