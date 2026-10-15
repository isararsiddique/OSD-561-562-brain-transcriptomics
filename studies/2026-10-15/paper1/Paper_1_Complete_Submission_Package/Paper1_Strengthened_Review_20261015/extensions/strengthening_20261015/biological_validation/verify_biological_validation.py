#!/usr/bin/env python3
"""Independent source-count, dense-sandwich and restricted-null verification."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy import stats
OUT=Path(__file__).resolve().parent;BASE=OUT.parents[2]
def source_path(record,key):
    p=Path(record[key]);return p if p.is_absolute()else BASE/p
checks=[]
def check(name,ok,error=None):checks.append(dict(name=name,passed=bool(ok),error=error))
results=pd.read_csv(OUT/'all_biological_effects.csv')
scores=pd.read_csv(OUT/'all_animal_endpoint_scores.csv')
availability=pd.read_csv(OUT/'endpoint_availability.csv')
configuration=json.loads((OUT/'run_configuration.json').read_text())
finite=results.p_wild.notna()
for method,column in [('bh','q_BH_complete_new_family'),('by','q_BY_complete_new_family')]:
    q=stats.false_discovery_control(results.loc[finite,'p_wild'].to_numpy(),method=method)
    check('Independent SciPy '+method+' full family',np.allclose(q,results.loc[finite,column],atol=1e-12))
check('Full wild and bootstrap configuration',configuration['wild']==99999 and configuration['bootstrap']==2000)
check('Unavailable rows retain unavailable probabilities',results.loc[~results.status.isin(['estimable','derived IUT']),'p_wild'].isna().all())
check('Every pooled observation labelled as pool, no independent-animal count',results.loc[results.experimental_unit.eq('pool'),'n_animals'].isna().all())
for rel,digest in json.loads((OUT/'source_hashes.json').read_text()).items():
    check('Source hash '+rel,hashlib.sha256((BASE/rel).read_bytes()).hexdigest()==digest)

for n in ['561','562']:
    count=pd.read_csv(BASE/f'inputs/OSD-{n}_counts.csv.gz',index_col=0)
    md=pd.read_csv(BASE/f'inputs/OSD-{n}_metadata.csv.gz',dtype=str).set_index('sample').loc[count.columns]
    count=count.T.groupby(md['Source Name'],sort=True).sum().T
    meta=md.reset_index().drop_duplicates('Source Name').set_index('Source Name').loc[count.columns]
    keep=count.index.astype(str).str.fullmatch(r'ENSMUSG\d+')&(count.mean(axis=1)>=10)
    count=count.loc[keep];refcount=count.loc[(count>0).all(axis=1)].to_numpy(float)
    geom=np.exp(np.mean(np.log(refcount),axis=1))
    # Sorted geometric-ratio median is independent of the production log-median expression.
    ratios=np.sort(refcount/geom[:,None],axis=0);g=len(ratios)
    sf=ratios[g//2]if g%2 else np.sqrt(ratios[g//2-1]*ratios[g//2])
    y=np.log2(count.to_numpy()/sf+1)
    available=availability[availability.accession=='OSD-'+n]
    outcome={}
    source=scores[scores.accession=='OSD-'+n].set_index('animal').loc[count.columns]
    for a in available.itertuples():
        ids=str(a.ensembl_ids).split('|');ix=count.index.get_indexer(ids)
        if (ix<0).any():continue
        values=y[ix].mean(axis=0);outcome[a.endpoint]=values
        error=float(np.max(abs(values-source[a.endpoint].to_numpy())))
        check('Source-count reconstructed '+n+'/'+a.endpoint,error<1e-11,error)
    age=meta['Factor Value[Age]'].astype(str).to_numpy()
    f=np.array(['_FL_'in a for a in count.columns],float)
    c=np.array(['_ISS_'in a for a in count.columns],float)
    xadd=np.column_stack([np.ones(len(meta)),age=='29',f,c,f*c]).astype(float)
    cells=[(a,cc,ff)for a in ['12','29']for cc in [0,1]for ff in [0,1]]
    xcell=np.column_stack([((age==a)&(c==cc)&(f==ff)).astype(float)for a,cc,ff in cells])
    for endpoint,values in outcome.items():
        part=results[(results.accession=='OSD-'+n)&(results.endpoint==endpoint)
             &(results.sensitivity=='all_animals')&results.kind.eq('direct_flight_collection_interaction')]
        for row in part.itertuples():
            if row.status!='estimable':continue
            subset=np.ones(len(meta),bool)
            if row.estimand=='additive_age_interaction':x=xadd.copy();ell=np.array([0,0,0,0,1.])
            elif row.estimand in ['age12_interaction','age29_interaction']:
                selected='12'if row.estimand.startswith('age12')else'29';subset=age==selected
                x=xcell[subset,:4]if selected=='12'else xcell[subset,4:];ell=np.array([1.,-1,-1,1])
            else:
                x=xcell.copy();ell=np.array([1.,-1,-1,1,1,-1,-1,1])/2 if row.estimand=='equal_age_interaction'else np.array([-1.,1,1,-1,1,-1,-1,1])
            value=values[subset]
            if row.composition_adjusted:
                if endpoint=='Erythroid7':continue
                markers=['Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Alas2','Slc4a1','Gypa']
                other=[k for k in markers if k!=endpoint]
                cov=np.vstack([outcome[k][subset]for k in other]).mean(axis=0)
                x=np.column_stack([x,cov]);ell=np.r_[ell,0.]
            coef=np.linalg.lstsq(x,value,rcond=None)[0];effect=float(ell@coef)
            bread=np.linalg.inv(x.T@x);res=value-x@coef;hat=np.diag(x@bread@x.T)
            meat=x.T@np.diag((res/(1-hat))**2)@x
            se=float(np.sqrt(ell@bread@meat@bread@ell))
            check('Dense source effect '+n+'/'+endpoint+'/'+row.estimand+'/'+str(row.composition_adjusted),abs(effect-row.effect_log2)<1e-10,float(abs(effect-row.effect_log2)))
            check('Independent matrix HC3 '+n+'/'+endpoint+'/'+row.estimand+'/'+str(row.composition_adjusted),abs(se-row.se)<1e-10,float(abs(se-row.se)))
    for estimand,plan_index in [('additive_age_interaction',2),('equal_age_interaction',8)]:
        row=results[(results.accession=='OSD-'+n)&results.endpoint.eq('Circadian')&results.estimand.eq(estimand)
                    &results.sensitivity.eq('all_animals')&~results.composition_adjusted].iloc[0]
        x=xadd if estimand.startswith('additive')else xcell
        ell=np.array([0.,0,0,0,1])if estimand.startswith('additive')else np.array([1.,-1,-1,1,1,-1,-1,1])/2
        value=outcome['Circadian'];pinv=np.linalg.pinv(x);bread=np.linalg.inv(x.T@x)
        h=np.diag(x@pinv);a=x@bread@ell;den=float(ell@bread@ell)
        beta=float(ell@pinv@value);fitted=x@pinv@value
        fit0=fitted-a*beta/den;h0=h-a*a/den;nullres=(value-fit0)/np.sqrt(1-h0)
        observed=beta/float(row.se);weightsvec=a/(1-h);residual_operator=np.eye(len(x))-x@pinv
        rng=np.random.default_rng(configuration['seed']+int(n)*1000+plan_index);exceed=0
        for start in range(0,configuration['wild'],1000):
            batch=min(1000,configuration['wild']-start)
            weights=rng.choice([-1.,1.],size=(batch,len(x)))
            fake=fit0+weights*nullres;eff=fake@(pinv.T@ell);res=fake@residual_operator.T
            se=np.sqrt(res**2@weightsvec**2);t=eff/se
            exceed+=int((abs(t)>=abs(observed)-1e-12).sum())
        probability=(exceed+1)/(configuration['wild']+1)
        check('Independent constrained-null wild p '+n+'/'+estimand,abs(probability-row.p_wild)<1e-12,abs(probability-row.p_wild))

for record in json.loads((OUT/'independent_datasets.json').read_text()):
    acc=record['accession'];count=pd.read_csv(source_path(record,'counts'),index_col=0)
    md=pd.read_csv(source_path(record,'metadata'),dtype=str).set_index('sample').loc[count.columns]
    count=count.T.groupby(md.animal,sort=True).sum().T
    meta=md.reset_index().drop_duplicates('animal').set_index('animal').loc[count.columns]
    keep=count.index.astype(str).str.fullmatch(r'ENSMUSG\d+')&(count.mean(axis=1)>=10);count=count.loc[keep]
    positive=count.loc[(count>0).all(axis=1)].to_numpy(float)
    geom=np.exp(np.log(positive).mean(axis=1));ratios=np.sort(positive/geom[:,None],axis=0);g=len(ratios)
    sf=ratios[g//2]if g%2 else np.sqrt(ratios[g//2-1]*ratios[g//2]);y=np.log2(count.to_numpy()/sf+1)
    source=scores[scores.accession==acc].set_index('animal').loc[count.columns];outcome={}
    for a in availability[availability.accession==acc].itertuples():
        ix=count.index.get_indexer(str(a.ensembl_ids).split('|'))
        if (ix<0).any():continue
        val=y[ix].mean(axis=0);outcome[a.endpoint]=val;error=float(abs(val-source[a.endpoint].to_numpy()).max())
        check('External source count score '+acc+'/'+a.endpoint,error<1e-11,error)
    part=results[(results.accession==acc)&results.status.eq('estimable')]
    for row in part.itertuples():
        if row.endpoint not in outcome:continue
        mask=np.ones(len(meta),bool);value=outcome[row.endpoint]
        if acc=='OSD-352':
            flight=meta.flight.eq('FLT').to_numpy(float);x=np.column_stack([np.ones(len(meta)),flight]);ell=np.array([0.,1])
            if row.sensitivity.startswith('RIN'):
                rin=meta.rin.astype(float).to_numpy();x=np.column_stack([x,rin-rin.mean()]);ell=np.r_[ell,0.]
        elif acc in ['OSD-374','OSD-611']:
            if row.sensitivity=='NASA_sampletable_primary':mask=meta.primary_analysis_included.str.lower().eq('true').to_numpy()
            dose=meta.radiation_dose_gy.astype(float).to_numpy();levels=sorted(set(dose))
            high,low=[float(s)for s in row.estimand.replace('dose','').split('_minus_')]
            x=np.column_stack([(dose[mask]==a).astype(float)for a in levels]);ell=np.array([float(a==high)-float(a==low)for a in levels])
        else:
            time=meta.collection_after_treatment.to_numpy();r=meta.radiation.ne('non-irradiated').to_numpy(int)
            u=meta.hindlimb_unloading.eq('Hindlimb Unloaded').to_numpy(int)
            months=['1','4']
            if row.estimand.endswith('month1')and'_minus_'not in row.estimand:months=['1'];mask=time=='1'
            if row.estimand.endswith('month4')and'_minus_'not in row.estimand:months=['4'];mask=time=='4'
            keys=[(t,hu,rad)for t in months for hu in [0,1]for rad in [0,1]]
            x=np.column_stack([((time==t)&(u==hu)&(r==rad))[mask].astype(float)for t,hu,rad in keys])
            v=np.array([-1.,1,0,0])if row.estimand.startswith('radiation_at')else(np.array([-1.,0,1,0])if row.estimand.startswith('unloading')else np.array([1.,-1,-1,1]))
            ell=v if len(months)==1 else(np.r_[-v,v]if row.estimand.endswith('month4_minus_month1')else np.r_[v,v]/2)
            if row.sensitivity=='lane_adjusted':
                lanes=pd.get_dummies(meta.loc[mask,'sequencing_lanes'],drop_first=True).to_numpy(float);x=np.column_stack([x,lanes]);ell=np.r_[ell,np.zeros(lanes.shape[1])]
        value=value[mask];coef=np.linalg.lstsq(x,value,rcond=None)[0];effect=float(ell@coef)
        bread=np.linalg.inv(x.T@x);res=value-x@coef;hat=np.diag(x@bread@x.T)
        se=float(np.sqrt(ell@bread@(x.T@np.diag((res/(1-hat))**2)@x)@bread@ell))
        check('External dense source effect '+acc+'/'+row.endpoint+'/'+row.estimand+'/'+row.sensitivity,abs(effect-row.effect_log2)<1e-10,abs(effect-row.effect_log2))
        check('External independent sandwich '+acc+'/'+row.endpoint+'/'+row.estimand+'/'+row.sensitivity,abs(se-row.se)<1e-10,abs(se-row.se))

summary=dict(passed=all(x['passed']for x in checks),n_checks=len(checks),failed=[x for x in checks if not x['passed']],checks=checks,
             scope='Source counts, original panel means, dense HC3 sandwich, independent constrained-null wild draws and SciPy full-family correction')
(OUT/'independent_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(dict(passed=summary['passed'],n_checks=summary['n_checks'],n_failed=len(summary['failed']),first_failures=summary['failed'][:4]),indent=2))
