#!/usr/bin/env python3
"""Real-count-calibrated known-truth benchmark with explicit information tiers.

This does not label PCA as SVA. Published RUVg is reimplemented from its
matrix/log-input source; neuroCombat is the official Python package.
"""
import contextlib
import hashlib
import io
import json
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from statistics_core import learn_anchor, bounded_gamma, bh

SCENARIOS={
 'balanced_stable':dict(correlation=0.,flip=False,overlap=False,heteroskedastic=False),
 'correlated_stable':dict(correlation=.5,flip=False,overlap=False,heteroskedastic=False),
 'strongly_correlated_stable':dict(correlation=.83,flip=False,overlap=False,heteroskedastic=False),
 'correlated_changed_loading':dict(correlation=.5,flip=True,overlap=False,heteroskedastic=False),
 'correlated_sensor_biology':dict(correlation=.5,flip=False,overlap=True,heteroskedastic=False),
 'correlated_batch_dispersion':dict(correlation=.5,flip=False,overlap=False,heteroskedastic=True),
 'biological_null':dict(correlation=.5,flip=False,overlap=False,heteroskedastic=False,null=True),
}
METHODS=['raw','ground_score','ground_cap2','ground_cv_score',
         'known_batch_OLS','neuroCombat_EB','neuroCombat_noEB',
         'RUVg_matrix_k1_oracle_controls','RUVg_covariate_k1_oracle_controls',
         'PCA_covariate_k1']

def fit(y,x,z=None):
    design=np.column_stack([np.ones(len(x)),x]) if z is None else np.column_stack([np.ones(len(x)),x,z])
    if np.linalg.matrix_rank(design)<design.shape[1]:
        raise ValueError('Unidentifiable biological contrast')
    b=np.linalg.lstsq(design,y.T,rcond=None)[0]
    e=y.T-design@b;df=len(x)-design.shape[1]
    var=np.sum(e*e,axis=0)/df*np.linalg.inv(design.T@design)[1,1]
    t=np.divide(b[1],np.sqrt(var),out=np.zeros(y.shape[0]),where=var>0)
    return b[1],2*stats.t.sf(abs(t),df)

def norm_reference(training,target):
    # Every comparator sees the same reference; no target-label selection.
    positive=(training>0).all(axis=1)
    reference=np.log(training[positive]).mean(axis=1)
    def transform(c):
        logs=[]
        for col in c[positive].T:
            nz=col>0;logs.append(np.median(np.log(col[nz])-reference[nz]))
        sf=np.exp(logs)
        return np.log2(c/sf+1)
    return transform(training),transform(target)

def ruvg_log_matrix(y,controls,k=1):
    # Same algebra as RUVSeq matrix method with isLog=TRUE, center=TRUE,
    # tolerance=1e-8, drop=0. Scaling log base changes no subspace.
    center=y.T-y.T.mean(axis=0)
    u,s,_=np.linalg.svd(center[:,controls],full_matrices=False)
    k=min(k,int(np.sum(s>1e-8)))
    if not k:raise ValueError('RUVg control rank zero')
    w=u[:,:k];alpha=np.linalg.solve(w.T@w,w.T@y.T)
    return (y.T-w@alpha).T,w

def pca_factor(y,k=1):
    center=y.T-y.T.mean(axis=0)
    u,_,_=np.linalg.svd(center,full_matrices=False)
    return u[:,:k]

def draw_counts(rng,base,alpha,biobeta,nuisance,x,z,shift=False,hetero=False,sensor=None):
    load=nuisance.copy()
    if shift:
        outside=np.ones(len(base),bool);outside[sensor]=False;load[outside]*=-1
    sf=np.exp(rng.normal(0,.25,len(x)))
    # Low-rank latent biological variation remains independent of labels.
    mu=base[:,None]*sf*np.exp2(biobeta[:,None]*x+load[:,None]*z)
    dispersion=alpha[:,None]*np.ones((1,len(x)))
    if hetero:dispersion[:,z>0]*=2
    shape=1/dispersion
    return rng.negative_binomial(shape,shape/(shape+mu))

def label_design(n,correlation):
    half=n//2
    # Symmetric contingency table; actual sample correlation is recorded.
    matching=int(round(n*(1+correlation)/4))
    matching=min(half-1,max(1,matching));mismatching=half-matching
    x=np.repeat([-.5,.5],half)
    z=np.array([-.5]*matching+[.5]*mismatching+[-.5]*mismatching+[.5]*matching)
    return x,z

def score_cv(yt,z,sensor):
    grid=np.linspace(0,1,11);risk=np.zeros(len(grid))
    outside=np.setdiff1d(np.arange(len(yt)),sensor)
    for held in range(yt.shape[1]):
        keep=np.arange(yt.shape[1])!=held
        _,w,gamma=learn_anchor(yt[:,keep],z[keep],None,sensor)
        resid=yt[:,held]-yt[:,keep].mean(axis=1)
        pred=gamma[outside]*(w@resid)
        risk+=np.sum((resid[outside,None]-pred[:,None]*grid)**2,axis=0)
    return float(grid[np.argmin(risk)])

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--repetitions',type=int,default=100);args=ap.parse_args()
    protocol={'scope':'Exploratory benchmark fixed before outcomes; no acceptance prediction',
              'seed':15102026561,'scenarios':SCENARIOS,'methods':METHODS,
              'repetitions_per_setting':args.repetitions,'n_train':[8,24],'n_target':24,
              'n_genes':600,'n_sensor':30,'n_oracle_negative_controls':120,
              'calibration':'OSD-561 actual count means and within-age/flight/cohort moment dispersions; genes sampled without outcome selection',
              'ground_information':'Control-only source cohort; no target batch or biological labels to fit correction',
              'target_information':'Known-batch OLS and official ComBat use target batch labels and protect the biological design; RUVg uses true synthetic negative-control genes',
              'PCA':'Target-expression PCA factor covariate k=1; not SVA',
              'RUVg':'Python algebra reimplementation of RUVSeq matrix/log-input source; k=1; matrix correction and covariate adjustment both shown',
              'truth':'Known log2 mean-count biological coefficient; linear log-expression estimators are deliberately model-misspecified for NB observations',
              'FDR':'Within each dataset and method all600 gene probabilities form BH family; FDP compares with true coefficient support',
              'uncertainty':'Independent repetition distributions and Monte Carlo intervals; simulated samples are not real animals',
              'missing':'Algorithm failure retained, never replaced by another method',
              'limitations':'Generic NB mechanism calibrated to empirical marginal count statistics; not full RNA-seq process or proof of external biological recovery'}
    (HERE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    from neuroCombat import neuroCombat
    import neuroCombat as module
    modulepath=Path(module.__file__).parent/'neuroCombat.py'
    # Save verified comparator source with hash and package identity.
    (HERE/'neuroCombat_0_2_12_source.py').write_bytes(modulepath.read_bytes())
    (HERE/'comparator_provenance.json').write_text(json.dumps({
        'neuroCombat_package':'0.2.12','package_url':'https://pypi.org/project/neuroCombat/0.2.12/',
        'implementation_sha256':hashlib.sha256(modulepath.read_bytes()).hexdigest(),
        'RUVg_source_url':'https://raw.githubusercontent.com/drisso/RUVSeq/master/R/RUVg-methods.R',
        'RUVg_source_sha256':hashlib.sha256((HERE/'RUVg_primary_source.R').read_bytes()).hexdigest()},indent=2)+'\n')
    counts=pd.read_csv(ROOT/'inputs/OSD-561_counts.csv.gz',index_col=0)
    md=pd.read_csv(ROOT/'results/sample_manifest.csv');md=md[md.accession=='OSD-561'].set_index('sample').loc[counts.columns]
    groups=md.groupby(['flight','cohort','age']).indices
    c=counts.to_numpy(float);mean=c.mean(axis=1);within=np.zeros(len(c));dof=0
    for idx in groups.values():
        within+=np.sum((c[:,idx]-c[:,idx].mean(axis=1,keepdims=True))**2,axis=1);dof+=len(idx)-1
    alpha=np.clip((within/max(dof,1)-mean)/(mean*mean+1e-12),.02,.5)
    eligible=np.flatnonzero(counts.index.str.match(r'^ENSMUSG\d+$')&(mean>=200))
    rng=np.random.default_rng(protocol['seed']);chosen=rng.choice(eligible,600,replace=False)
    base=mean[chosen];dispersion=alpha[chosen]
    pd.DataFrame({'gene':counts.index[chosen],'count_mean':base,'NB_dispersion':dispersion}).to_csv(HERE/'real_count_calibration.csv',index=False)
    rows=[]
    for scenario,config in SCENARIOS.items():
      for ntrain in protocol['n_train']:
       for rep in range(args.repetitions):
        rg=np.random.default_rng(protocol['seed']+1000000*list(SCENARIOS).index(scenario)+1000*ntrain+rep)
        sensor=np.arange(30);controls=np.arange(30,150);outside=np.arange(30,600)
        beta=np.zeros(600);bio=np.arange(150,210);beta[bio]=rg.choice([-1,1],len(bio))*.8
        if config['overlap']:beta[sensor[:8]]=rg.choice([-1,1],8)*.8
        if config.get('null'):beta[:]=0
        load=rg.normal(.2,.25,600);load[sensor]=rg.normal(.8,.1,30)
        ztrain=np.repeat([-.5,.5],ntrain//2);xtrain=np.zeros(ntrain)
        xt,zt=label_design(24,config['correlation'])
        train=draw_counts(rg,base,dispersion,beta,load,xtrain,ztrain,sensor=sensor)
        target=draw_counts(rg,base,dispersion,beta,load,xt,zt,config['flip'],config['heteroskedastic'],sensor)
        yt,ye=norm_reference(train,target)
        _,w,gamma=learn_anchor(yt,ztrain,None,sensor);shrink=score_cv(yt,ztrain,sensor)
        covariance=pd.DataFrame({'batch':(zt>0).astype(int),'biology':(xt>0).astype(int)})
        for method in METHODS:
         record={'scenario':scenario,'n_train':ntrain,'n_target':24,'repetition':rep,
                 'method':method,'actual_batch_biology_correlation':float(np.corrcoef(xt,zt)[0,1]),
                 'original_loading_gain':float(np.linalg.norm(gamma)),'cv_strength':shrink,'status':'PASS'}
         try:
          if method=='raw':effect,p=fit(ye,xt)
          elif method in ('ground_score','ground_cap2','ground_cv_score'):
            g=gamma if method=='ground_score' else bounded_gamma(w,gamma,2) if method=='ground_cap2' else w+shrink*(gamma-w)
            yc=ye-np.outer(g,w@ye);effect,p=fit(yc,xt);record['operator_gain']=float(np.linalg.norm(g))
          elif method=='known_batch_OLS':effect,p=fit(ye,xt,zt)
          elif method.startswith('neuroCombat'):
            with contextlib.redirect_stdout(io.StringIO()),warnings.catch_warnings():
              warnings.simplefilter('ignore')
              yc=neuroCombat(ye,covariance,'batch',categorical_cols=['biology'],eb=method=='neuroCombat_EB')['data']
            effect,p=fit(yc,xt)
          elif method.startswith('RUVg'):
            yc,factor=ruvg_log_matrix(ye,controls)
            effect,p=fit(yc,xt) if method.startswith('RUVg_matrix') else fit(ye,xt,factor)
          else:effect,p=fit(ye,xt,pca_factor(ye))
          null=beta==0;nonnull=~null;discover=bh(p)<=.05
          record.update(effect_RMSE=float(np.sqrt(np.mean((effect-beta)**2))),
                        outside_effect_RMSE=float(np.sqrt(np.mean((effect[outside]-beta[outside])**2))),
                        biological_FDP=float(np.sum(discover&null)/max(1,np.sum(discover))),
                        power=float(np.sum(discover&nonnull)/np.sum(nonnull)) if nonnull.any() else np.nan,
                        n_discoveries=int(discover.sum()),null_gene_nominal_rejection=float(np.mean(p[null]<.05)),
                        signed_biological_bias=float(np.mean(effect[nonnull]*np.sign(beta[nonnull])-abs(beta[nonnull]))) if nonnull.any() else np.nan)
         except Exception as exc:record.update(status='FAIL',failure=str(exc))
         rows.append(record)
        if rep%25==0:print(scenario,ntrain,rep,'complete',flush=True)
    results=pd.DataFrame(rows);results.to_csv(HERE/'all_repetitions.csv',index=False)
    metrics=['effect_RMSE','outside_effect_RMSE','biological_FDP','power','null_gene_nominal_rejection']
    summaries=[]
    for (scenario,ntrain,method),part in results.groupby(['scenario','n_train','method']):
        row=dict(scenario=scenario,n_train=ntrain,method=method,n_repetitions=len(part),n_failures=int(part.status.eq('FAIL').sum()))
        for metric in metrics:
            values=part[metric].dropna().to_numpy();row[metric+'_mean']=float(values.mean()) if len(values) else np.nan
            row[metric+'_MCse']=float(values.std(ddof=1)/np.sqrt(len(values))) if len(values)>1 else np.nan
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(HERE/'benchmark_summary.csv',index=False)
    # Paired differences use independently generated datasets as their units.
    paired=[]
    pivot=results.pivot(index=['scenario','n_train','repetition'],columns='method',values='effect_RMSE')
    for (scenario,ntrain),part in pivot.groupby(level=[0,1]):
      for method in METHODS[1:]:
        diff=(part[method]-part['raw']).dropna();paired.append(dict(scenario=scenario,n_train=ntrain,method=method,
          mean_RMSE_difference_to_raw=float(diff.mean()),MCse=float(diff.std(ddof=1)/np.sqrt(len(diff))),
          n_independent_pairs=len(diff),fraction_repetitions_better_than_raw=float((diff<0).mean())))
    pd.DataFrame(paired).to_csv(HERE/'paired_method_differences.csv',index=False)
    print('Complete',len(results),'method rows',int(results.status.eq('FAIL').sum()),'failures',flush=True)

if __name__=='__main__':main()
