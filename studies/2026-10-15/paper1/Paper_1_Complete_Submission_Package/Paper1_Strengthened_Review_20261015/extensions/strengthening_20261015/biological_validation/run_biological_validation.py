#!/usr/bin/env python3
"""Complete animal-level biological sensitivity; no manuscript/package edits."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import linalg, stats

OUT=Path(__file__).resolve().parent
BASE=OUT.parents[2]
sys.path.insert(0,str(BASE/'scripts'))
from run_revision import load_dataset, PANELS, indices

MARKERS=['Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Alas2','Slc4a1','Gypa']
SEED=20261015561562
def source_path(record,key):
    p=Path(record[key]);return p if p.is_absolute()else BASE/p

def bh(values, by=False):
    p=np.asarray(values,float);finite=np.isfinite(p);q=np.full(len(p),np.nan)
    order=np.argsort(p[finite],kind='stable');s=p[finite][order];n=len(s)
    ranks=np.arange(1,n+1)
    v=np.minimum(1,np.minimum.accumulate((s*n/ranks)[::-1])[::-1])
    if by:v=np.minimum(1,v*np.sum(1/ranks))
    z=np.empty(n);z[order]=v;q[finite]=z;return q

def independent_normalize(c):
    keep=c.index.astype(str).str.fullmatch(r'ENSMUSG\d+')&(c.mean(axis=1)>=10)
    c=c.loc[keep].astype(float);positive=(c>0).all(axis=1)
    logs=np.log(c.loc[positive].to_numpy());ref=logs.mean(axis=1)
    factors=np.exp(np.median(logs-ref[:,None],axis=0))
    return np.log2(c.to_numpy()/factors+1),factors,c.index

def summarize(d, annotation):
    symbols=annotation.set_index('ensembl_id').gene_symbol.to_dict()
    endpoints={p:idx for p,idx in indices(d['y'].index,annotation).items()}
    syms=np.array([symbols.get(g,'') for g in d['y'].index])
    for marker in MARKERS:endpoints[marker]=np.flatnonzero(syms==marker)
    endpoints['Erythroid7']=np.flatnonzero(np.isin(syms,MARKERS))
    names=list(endpoints);scores=np.vstack([d['y'].iloc[ix].mean(axis=0).to_numpy()
        if len(ix) else np.full(len(d['meta']),np.nan) for ix in endpoints.values()])
    return names,scores,endpoints,syms

def matrix_fit(y,x,contrast):
    """Outcome rows, animal columns. General contrast, independent HC3 form."""
    if np.linalg.matrix_rank(x)<x.shape[1]:return None
    inv=np.linalg.inv(x.T@x);map_=inv@x.T;coefs=y@map_.T
    fitted=coefs@x.T;res=y-fitted;h=np.einsum('ij,ji->i',x,map_)
    effect=coefs@contrast
    if np.any(h>=1-1e-9):return dict(effect=effect,status='unit leverage; uncertainty unavailable')
    influence=(x@inv@contrast)/(1-h)
    se=np.sqrt(res**2@influence**2)
    statistic=np.divide(effect,se,out=np.full(len(y),np.nan),where=se>0)
    k=len(contrast);basis=linalg.null_space(contrast.reshape(1,-1));x0=x@basis
    map0=np.linalg.pinv(x0);fit0=(y@map0.T)@x0.T
    h0=np.einsum('ij,ji->i',x0,map0)
    null_res=(y-fit0)/np.sqrt(1-h0)
    rss=np.sum(res**2,axis=1);rss0=np.sum((y-fit0)**2,axis=1)
    r2=np.divide(rss0-rss,rss0,out=np.full(len(y),np.nan),where=rss0>0)
    return dict(effect=effect,se=se,statistic=statistic,coefs=coefs,res=res,h=h,
                map=map_,influence=influence,fit0=fit0,null_res=null_res,
                partial_R2=r2,status='estimable',df=len(x)-k)

def inference(y,x,contrast,strata,nwild,nboot,seed):
    fit=matrix_fit(y,x,contrast)
    if fit is None:return dict(status='rank deficient; effect and uncertainty unavailable')
    if fit['status']!='estimable':return fit
    rng=np.random.default_rng(seed);exceed=np.zeros(len(y),int)
    for start in range(0,nwild,1000):
        size=min(1000,nwild-start);weights=rng.choice([-1.,1.],size=(size,len(x)))
        ys=fit['fit0'][None,:,:]+fit['null_res'][None,:,:]*weights[:,None,:]
        coef=ys@fit['map'].T;res=ys-coef@x.T
        se=np.sqrt(res**2@fit['influence']**2);eff=coef@contrast
        ts=np.divide(eff,se,out=np.zeros_like(eff),where=se>0)
        exceed+=(abs(ts)>=abs(fit['statistic'])[None,:]-1e-12).sum(axis=0)
    boot=[];failed=0
    for _ in range(nboot):
        chosen=np.concatenate([rng.choice(s,len(s),replace=True) for s in strata])
        xb=x[chosen]
        if np.linalg.matrix_rank(xb)<x.shape[1]:failed+=1;continue
        boot.append(np.linalg.lstsq(xb,y[:,chosen].T,rcond=None)[0].T@contrast)
    boot=np.array(boot)
    fit.update(p_wild=(exceed+1)/(nwild+1),
        p_wild_MC_lower95=np.where(exceed==0,0,stats.beta.ppf(.025,exceed,nwild-exceed+1)),
        p_wild_MC_upper95=np.where(exceed==nwild,1,stats.beta.ppf(.975,exceed+1,nwild-exceed)),
        bootstrap_lower95=np.quantile(boot,.025,axis=0),
        bootstrap_upper95=np.quantile(boot,.975,axis=0),
        n_successful_bootstraps=len(boot),n_failed_bootstraps=failed,
        p_HC3_t_approximate=2*stats.t.sf(abs(fit['statistic']),fit['df']))
    return fit

def cells(meta,ages):
    keys=[(a,f,c) for a in ages for c in ['LAR','ISS_T'] for f in ['GC','FLT']]
    x=np.column_stack([((meta.age==a)&(meta.flight==f)&(meta.cohort==c)).to_numpy(float)
                       for a,f,c in keys])
    return x,keys

def plans(d,core=False):
    m=d['meta'];acc=m.accession.iloc[0];result=[]
    if acc=='OSD-202':
        time=m.collection_after_treatment.astype(str)
        rad=(m.radiation!='non-irradiated').to_numpy(int)
        unload=(m.hindlimb_unloading=='Hindlimb Unloaded').to_numpy(int)
        unit=np.array([-1.,1,0,0]);unloadunit=np.array([-1.,0,1,0]);inter=np.array([1.,-1,-1,1])
        base=[]
        for tt in ['1','4']:
            mask=(time==tt).to_numpy();keys=[(tt,u,r)for u in [0,1]for r in [0,1]]
            x=np.column_stack([((time==t).to_numpy()&(unload==u)&(rad==r))[mask].astype(float)for t,u,r in keys])
            support={str(key):int(x[:,j].sum())for j,key in enumerate(keys)}
            for label,c in [('radiation_at_normal_loading',unit),('unloading_without_radiation',unloadunit),('radiation_by_unloading',inter)]:
                base.append(dict(name=label+'_month'+tt,mask=mask,x=x,c=c,cell_support=support))
        keys=[(tt,u,r)for tt in ['1','4']for u in [0,1]for r in [0,1]]
        x=np.column_stack([((time==t).to_numpy()&(unload==u)&(rad==r)).astype(float)for t,u,r in keys]);mask=np.ones(len(m),bool)
        support={str(key):int(x[:,j].sum())for j,key in enumerate(keys)}
        for label,v in [('radiation_at_normal_loading',unit),('unloading_without_radiation',unloadunit),('radiation_by_unloading',inter)]:
            for suffix,c in [('equal_time',np.r_[v,v]/2),('month4_minus_month1',np.r_[-v,v])]:
                base.append(dict(name=label+'_'+suffix,mask=mask,x=x,c=c,cell_support=support))
        for p in base:
            for lane in [False,True]:
                x=p['x'];c=p['c']
                if lane:
                    dummy=pd.get_dummies(m.loc[p['mask'],'sequencing_lanes'],drop_first=True).to_numpy(float)
                    x=np.column_stack([x,dummy]);c=np.r_[c,np.zeros(dummy.shape[1])]
                result.append(dict(p,x=x,c=c,sensitivity='lane_adjusted'if lane else'all_animals',
                    composition=False,kind='Earth_radiation_unloading_factorial',
                    strata=['radiation','hindlimb_unloading','collection_after_treatment']))
        return result
    if acc in ['OSD-374','OSD-611']:
        dose=m.radiation_dose_gy.astype(float)
        subsets=['NASA_sampletable_primary','all_deposited_pools_sensitivity']if acc=='OSD-611'else['all_deposited_pools']
        levels=sorted(dose.unique())
        for sensitivity in subsets:
            mask=np.ones(len(m),bool)if sensitivity!='NASA_sampletable_primary'else m.primary_analysis_included.str.lower().eq('true').to_numpy()
            mm=m.loc[mask];x=np.column_stack([(dose.loc[mask]==a).to_numpy(float)for a in levels]);support={str(a):int(x[:,j].sum())for j,a in enumerate(levels)}
            for positive,negative in itertools.combinations(levels,2):
                c=np.array([float(a==negative)-float(a==positive)for a in levels])
                result.append(dict(name='dose'+str(negative)+'_minus_dose'+str(positive),sensitivity=sensitivity,
                    mask=mask,x=x,c=c,cell_support=support,kind='Earth_radiation_pool_context',composition=False,
                    strata=['radiation_dose_gy']))
        return result
    conflict=set(json.loads((BASE/'extensions/flight_collection_interaction/analysis_summary.json').read_text())[
        'conflicting_animals_either_region'])
    sensitivities=['all_animals','exclude_any_region_metadata_conflict'] if core else ['all_animals']
    for sensitivity in sensitivities:
        keep=np.ones(len(m),bool) if sensitivity=='all_animals' else ~m.animal.isin(conflict).to_numpy()
        if core:
            # Mission-matched flight contrasts, conditional on collection cohort.
            for cohort in ['ISS_T','LAR']:
                mask=keep&(m.cohort==cohort).to_numpy();mm=m.loc[mask]
                f=(mm.flight=='FLT').to_numpy(float)
                cols=[np.ones(len(mm))]
                if mm.age.nunique()>1:cols.append((mm.age=='29').to_numpy(float))
                cols.append(f);x=np.column_stack(cols);c=np.zeros(x.shape[1]);c[-1]=1
                result.append(dict(name='FLT_minus_GC_'+cohort,sensitivity=sensitivity,
                    mask=mask,x=x,c=c,kind='flight_minus_matched_ground',composition=False,
                    strata=['flight','age'],cell_support=mm.groupby(['flight','age']).size().to_dict()))
            # Original additive-age estimand, retained as sensitivity.
            mm=m.loc[keep];f=(mm.flight=='FLT').to_numpy(float);c0=(mm.cohort=='ISS_T').to_numpy(float)
            x=np.column_stack([np.ones(len(mm)),(mm.age=='29').to_numpy(float),f,c0,f*c0]);c=np.array([0,0,0,0,1.])
            complete_support={(a,f,c):int(((mm.age==a)&(mm.flight==f)&(mm.cohort==c)).sum())
                              for a in ['12','29']for f in ['GC','FLT']for c in ['LAR','ISS_T']}
            base=[dict(name='additive_age_interaction',mask=keep,x=x,c=c,cell_support=complete_support)]
            for age in ['12','29']:
                mask=keep&(m.age==age).to_numpy();mm=m.loc[mask];x,keys=cells(mm,[age]);c=np.array([1,-1,-1,1.])
                support={(a,f,c):int(((mm.age==a)&(mm.flight==f)&(mm.cohort==c)).sum())for a,f,c in keys}
                base.append(dict(name='age'+age+'_interaction',mask=mask,x=x,c=c,cell_support=support))
            mm=m.loc[keep];x,keys=cells(mm,['12','29']);c=np.array([1,-1,-1,1,1,-1,-1,1.])/2
            support={(a,f,c):int(((mm.age==a)&(mm.flight==f)&(mm.cohort==c)).sum())for a,f,c in keys}
            base.append(dict(name='equal_age_interaction',mask=keep,x=x,c=c,cell_support=support))
            c=np.array([-1,1,1,-1,1,-1,-1,1.])
            base.append(dict(name='age29_minus_age12_interaction',mask=keep,x=x,c=c,cell_support=support))
            for row in base:
                for composition in [False,True]:
                    result.append(dict(row,sensitivity=sensitivity,kind='direct_flight_collection_interaction',
                        composition=composition,strata=['flight','cohort','age']))
        else:
            comparisons=[('FLT','HGC'),('FLT','VGC')]if acc=='OSD-525'else[('FLT','GC')]
            for positive,negative in comparisons:
                mask=keep&m.flight.isin([positive,negative]).to_numpy();mm=m.loc[mask]
                cols=[np.ones(len(mm))]
                if mm.age.nunique()>1:cols.append((mm.age=='29').to_numpy(float))
                cols.append((mm.flight==positive).to_numpy(float));x=np.column_stack(cols);c=np.zeros(x.shape[1]);c[-1]=1
                result.append(dict(name=positive+'_minus_'+negative,sensitivity=sensitivity,mask=mask,x=x,c=c,
                    kind='flight_minus_matched_ground' if negative!='VGC'else'flight_minus_vivarium_control',
                    composition=False,strata=['flight','age'],cell_support=mm.groupby(['flight','age']).size().to_dict()))
                if acc=='OSD-352'and 'rin'in mm and mm.rin.notna().all()and(mm.rin!='').all():
                    rin=mm.rin.astype(float).to_numpy();rin=rin-rin.mean();xr=np.column_stack([x,rin]);cr=np.r_[c,0.]
                    result.append(dict(name=positive+'_minus_'+negative,sensitivity='RIN_adjusted_no_exposure_overlap',
                        mask=mask,x=xr,c=cr,kind='flight_minus_matched_ground',composition=False,
                        strata=['flight','age'],cell_support=mm.groupby(['flight','age']).size().to_dict()))
    return result

def execute(args):
    ann=pd.read_csv(BASE/'inputs/gene_annotation.csv.gz')
    datasets={n:load_dataset(n) for n in ['561','562','563','564','612','613','525']}
    # An independently audited acquisition manifest can add standardized deposits.
    extra=OUT/'independent_datasets.json'
    if extra.exists():
        for record in json.loads(extra.read_text()):
            c=pd.read_csv(source_path(record,'counts'),index_col=0)
            if not c.index.is_unique or not c.columns.is_unique or not np.isfinite(c.to_numpy()).all()or(c<0).any().any()or not np.array_equal(c.to_numpy(),np.round(c.to_numpy())):
                raise ValueError('Independent count source is not unique finite nonnegative integer data')
            m=pd.read_csv(source_path(record,'metadata'),dtype=str).set_index('sample')
            if set(c.columns)-set(m.index):raise ValueError('Independent counts lack metadata')
            m=m.loc[c.columns];collapsed=c.T.groupby(m.animal,sort=True).sum().T
            meta=m.reset_index().drop_duplicates('animal').set_index('animal',drop=False).loc[collapsed.columns]
            collapsed.columns=[record['accession']+':'+a for a in collapsed.columns];meta.index=collapsed.columns
            y,sf,genes=independent_normalize(collapsed)
            datasets[record['accession'].split('-')[-1]]=dict(y=pd.DataFrame(y,index=genes,columns=collapsed.columns),
                counts=collapsed,meta=meta,sf=sf,record=record)
    effects=[];omitanimal=[];omitgene=[];score_rows=[];avail=[];support=[];checks=[];hashes={}
    for n,d in datasets.items():
        m=d['meta'];acc=m.accession.iloc[0];names,scores,endpoints,symbols=summarize(d,ann)
        raw,sf,genes=independent_normalize(d['counts'])
        if not np.allclose(raw,d['y'].loc[genes].to_numpy(),atol=1e-12):raise ValueError('Independent normalization mismatch')
        checks.append(dict(check=acc+' independent normalization',passed=True,max_error=float(abs(raw-d['y'].to_numpy()).max())))
        for name,ix in endpoints.items():
            expected=PANELS.get(name,[name] if name!='Erythroid7'else MARKERS)
            present=sorted(set(symbols[ix]));avail.append(dict(accession=acc,endpoint=name,n_genes=len(ix),
                expected_symbols='|'.join(expected),present_symbols='|'.join(present),
                missing_symbols='|'.join(s for s in expected if s not in present),
                ensembl_ids='|'.join(d['y'].index[ix])))
        for j,row in enumerate(m.itertuples()):
            q=row._asdict();q['analysis_sample']=m.index[j];q.update({name:float(scores[i,j])for i,name in enumerate(names)})
            score_rows.append(q)
        for plan_index,p in enumerate(plans(d,core=n in ['561','562'])):
            mask=p['mask'];mm=m.loc[mask];base_x=p['x'];base_c=p['c'];yy=scores[:,mask]
            for cell,count in p['cell_support'].items():support.append(dict(accession=acc,estimand=p['name'],
                sensitivity=p['sensitivity'],cell=str(cell),n_animals=count))
            strata=[group.index.to_numpy()for _,group in mm.reset_index(drop=True).groupby(p['strata'],sort=True)]
            if p['composition']:
                tasks=[(list(range(len(PANELS))),yy[:len(PANELS)],yy[names.index('Erythroid7')])]
                for marker in MARKERS:
                    i=names.index(marker);other=[names.index(s)for s in MARKERS if s!=marker]
                    tasks.append(([i],yy[i:i+1],np.nanmean(yy[other],axis=0)))
                # Adjusting the marker mean for itself has no meaningful estimand.
                i=names.index('Erythroid7');tasks.append(([i],yy[i:i+1],None))
            else:tasks=[(list(range(len(names))),yy,None)]
            print(acc,p['name'],p['sensitivity'],'marker adjustment'if p['composition']else'raw',flush=True)
            independent_tasks=[]
            for positions,outcomes,covariate in tasks:
                valid=np.isfinite(outcomes).all(axis=1)
                if covariate is not None and not np.isfinite(covariate).all():valid[:]=False
                if valid.any():independent_tasks.append(([positions[j]for j in np.flatnonzero(valid)],outcomes[valid],covariate))
                for j in np.flatnonzero(~valid):independent_tasks.append(([positions[j]],outcomes[j:j+1],covariate))
            for positions,outcomes,covariate in independent_tasks:
                x=base_x if covariate is None else np.column_stack([base_x,covariate])
                c=base_c if covariate is None else np.append(base_c,0.)
                if p['composition'] and covariate is None:fit=dict(status='marker score cannot adjust for itself')
                elif not np.isfinite(outcomes).all()or not np.isfinite(x).all():fit=dict(status='endpoint or covariate unavailable')
                else:fit=inference(outcomes,x,c,strata,args.wild,args.bootstrap,SEED+int(n)*1000+plan_index)
                for pos,endpoint_index in enumerate(positions):
                    name=names[endpoint_index]
                    common=dict(accession=acc,mission=mm.mission.iloc[0],tissue=mm.tissue.iloc[0],
                        endpoint=name,estimand=p['name'],sensitivity=p['sensitivity'],
                        composition_adjusted=p['composition'],kind=p['kind'],n_independent_units=len(mm),
                        n_animals=np.nan if 'replicate_unit'in mm and mm.replicate_unit.str.contains('pool').any()else len(mm),
                        experimental_unit='pool' if 'replicate_unit'in mm and mm.replicate_unit.str.contains('pool').any()else'animal',
                        n_genes=len(endpoints[name]),status=fit['status'],n_wild_draws=args.wild,
                        animals='|'.join(mm.animal),source_overlap_note='one region per fit; mission contexts remain distinct',
                        source_design_warning=d.get('record',{}).get('design_warning',''),
                        minimum_design_cell_animals=min(p['cell_support'].values()),
                        age_cell_warning='empty or singleton design cells; additive-age estimands extrapolate or pool age'
                        if min(p['cell_support'].values())<2 else '')
                    if 'effect'in fit:common['effect_log2']=float(fit['effect'][pos])
                    if fit['status']=='estimable':
                        for field in ['se','partial_R2','p_wild','p_wild_MC_lower95','p_wild_MC_upper95',
                                      'bootstrap_lower95','bootstrap_upper95','p_HC3_t_approximate']:
                            common[field]=float(fit[field][pos])
                        common['max_animal_leverage']=float(fit['h'].max())
                        common['n_successful_bootstraps']=fit['n_successful_bootstraps']
                        common['n_failed_bootstraps']=fit['n_failed_bootstraps']
                        if covariate is not None:
                            basefit=matrix_fit(outcomes,base_x,base_c)
                            basecoef=np.linalg.lstsq(base_x,outcomes.T,rcond=None)[0].T
                            rssbase=np.sum((outcomes-basecoef@base_x.T)**2,axis=1)
                            rss=np.sum(fit['res']**2,axis=1)
                            common['marker_partial_R2_given_full_base_design']=float((rssbase[pos]-rss[pos])/rssbase[pos])
                            if basefit is not None and 'effect'in basefit:
                                common['unadjusted_same_design_effect_log2']=float(basefit['effect'][pos])
                                common['marker_adjustment_absolute_change_log2']=float(fit['effect'][pos]-basefit['effect'][pos])
                    effects.append(common)
                    if 'effect'not in fit:continue
                    for j,animal in enumerate(mm.animal):
                        retain=np.arange(len(mm))!=j;xr=x[retain]
                        effect=np.linalg.lstsq(xr,outcomes[pos,retain],rcond=None)[0]@c if np.linalg.matrix_rank(xr)==x.shape[1]else np.nan
                        omitanimal.append(dict(common,omitted_independent_unit=animal,omitted_effect_log2=float(effect),
                            omission_status='estimable'if np.isfinite(effect)else'rank deficient'))
                    ix=endpoints[name]
                    if len(ix)>1:
                        constituent=d['y'].iloc[ix, np.flatnonzero(mask)].to_numpy()
                        for j,gene in enumerate(ix):
                            leave=(constituent.sum(axis=0)-constituent[j])/(len(ix)-1)
                            effect=np.linalg.lstsq(x,leave,rcond=None)[0]@c
                            omitgene.append(dict(common,omitted_ensembl_id=d['y'].index[gene],
                                omitted_symbol=symbols[gene],omitted_effect_log2=float(effect)))
    # An additional identical gene definition within each anatomical group.
    flight_datasets={n:d for n,d in datasets.items()if 'record'not in d or d['record']['exposure_kind']=='spaceflight'}
    by_tissue={}
    for n,d in flight_datasets.items():by_tissue.setdefault(d['meta'].tissue.iloc[0],[]).append(n)
    harmonized=[];screens=[];differences=[]
    for tissue,cohorts in by_tissue.items():
        if len(cohorts)<2:continue
        common=flight_datasets[cohorts[0]]['y'].index
        for n in cohorts[1:]:common=common.intersection(flight_datasets[n]['y'].index)
        for n in cohorts:
            d=dict(flight_datasets[n]);d['y']=d['y'].loc[common]
            endpoints=indices(common,ann);names=list(PANELS)
            yy=np.vstack([d['y'].iloc[endpoints[name]].mean(axis=0).to_numpy()for name in names])
            for plan_index,p in enumerate(plans(d,core=n in ['561','562'])):
                if p['kind']!='flight_minus_matched_ground'or p['sensitivity']!='all_animals'or p['composition']:continue
                mm=d['meta'].loc[p['mask']];strata=[s.index.to_numpy()for _,s in mm.reset_index(drop=True).groupby(p['strata'])]
                fit=inference(yy[:,p['mask']],p['x'],p['c'],strata,args.wild,args.bootstrap,SEED+int(n)*1000+plan_index)
                for j,name in enumerate(names):
                    row=dict(accession=mm.accession.iloc[0],mission=mm.mission.iloc[0],tissue=tissue,
                        endpoint=name,estimand=p['name']+'_common_gene_definition',sensitivity='all_animals',
                        composition_adjusted=False,kind='harmonized_flight_context',n_independent_units=len(mm),
                        n_animals=len(mm),experimental_unit='animal',n_genes=len(endpoints[name]),
                        n_genome_aligned_genes=len(common),status=fit['status'],n_wild_draws=args.wild,
                        animals='|'.join(mm.animal),endpoint_gene_ids='|'.join(common[endpoints[name]]))
                    if 'effect'in fit:row['effect_log2']=float(fit['effect'][j])
                    if fit['status']=='estimable':
                        for field in ['se','partial_R2','p_wild','p_wild_MC_lower95','p_wild_MC_upper95','bootstrap_lower95','bootstrap_upper95','p_HC3_t_approximate']:
                            row[field]=float(fit[field][j])
                    harmonized.append(row);effects.append(row)
        h=pd.DataFrame(harmonized);h=h[h.tissue==tissue]
        for name,part in h.groupby('endpoint'):
            valid=part.status.eq('estimable').all();same=part.effect_log2.gt(0).all()or part.effect_log2.lt(0).all()
            probability=float(part.p_wild.max())if valid and same else 1.
            screen=dict(accession='cross_context_screen',mission='multiple observational contexts',tissue=tissue,
                endpoint=name,estimand='all_available_flight_contexts_direction_IUT',sensitivity='all_animals',
                composition_adjusted=False,kind='harmonized_panel_direction_screen',status='derived IUT',
                p_wild=probability,n_components=len(part),n_independent_missions=part.mission.nunique(),
                direction_concordant=bool(same),minimum_component_effect=float(part.effect_log2.min()),
                maximum_component_effect=float(part.effect_log2.max()),endpoint_gene_ids=part.endpoint_gene_ids.iloc[0],
                interpretation='conservative common-direction screen; no pooled causal exposure estimand')
            screens.append(screen);effects.append(screen)
            for (_,a),(_,b) in itertools.combinations(part.iterrows(),2):
                shared=set(a.animals.split('|'))&set(b.animals.split('|'))
                delta=float(a.effect_log2-b.effect_log2);se=np.sqrt(a.se*a.se+b.se*b.se)if not shared else np.nan
                differences.append(dict(tissue=tissue,endpoint=name,context_a=a.accession+' '+a.estimand,
                    context_b=b.accession+' '+b.estimand,effect_difference=delta,
                    descriptive_Wald_lower95=delta-1.96*se,descriptive_Wald_upper95=delta+1.96*se,
                    n_shared_source_units=len(shared),different_missions=a.mission!=b.mission,
                    interpretation='context heterogeneity; approximate descriptive uncertainty, not a meta-analysis'))
    frame=pd.DataFrame(effects)
    frame['q_BH_complete_new_family']=bh(frame.get('p_wild',pd.Series(np.nan,index=frame.index)))
    frame['q_BY_complete_new_family']=bh(frame.get('p_wild',pd.Series(np.nan,index=frame.index)),by=True)
    frame['new_family_finite_tests']=frame.p_wild.notna().sum()
    frame.to_csv(OUT/'all_biological_effects.csv',index=False)
    frame[frame.kind=='harmonized_flight_context'].to_csv(OUT/'harmonized_flight_panel_effects.csv',index=False)
    frame[frame.kind=='harmonized_panel_direction_screen'].to_csv(OUT/'harmonized_panel_direction_screen.csv',index=False)
    pd.DataFrame(differences).to_csv(OUT/'all_cross_context_effect_differences.csv',index=False)
    pd.DataFrame(omitanimal).to_csv(OUT/'all_animal_omissions.csv',index=False)
    pd.DataFrame(omitgene).to_csv(OUT/'all_gene_omissions.csv',index=False)
    pd.DataFrame(score_rows).to_csv(OUT/'all_animal_endpoint_scores.csv',index=False)
    pd.DataFrame(avail).to_csv(OUT/'endpoint_availability.csv',index=False)
    pd.DataFrame(support).drop_duplicates().to_csv(OUT/'design_cell_support.csv',index=False)
    for path in [BASE/'inputs/gene_annotation.csv.gz',OUT/'EXPLORATORY_PROTOCOL.md',
                 *[BASE/'inputs'/f'OSD-{n}_{kind}.csv.gz'for n in datasets if n in ['561','562','563','564','612','613','525']for kind in ['counts','metadata']]]:
        hashes[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
    for path in [OUT/'EXTERNAL_DESIGN_AMENDMENT.md',Path(__file__),OUT/'prepare_external_metadata.py']:
        hashes[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
    if extra.exists():
        for record in json.loads(extra.read_text()):
            for name in ['counts','metadata','source_metadata']:
                path=source_path(record,name);hashes[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT/'source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    (OUT/'run_configuration.json').write_text(json.dumps(dict(seed=SEED,wild=args.wild,
        bootstrap=args.bootstrap,n_effect_rows=len(frame),n_finite_p=int(frame.p_wild.notna().sum()),
        checks=checks,warning='Exploratory observational inference; no biological cell fractions or causal mediation'),indent=2)+'\n')
    print('Finished',len(frame),'effect rows',frame.p_wild.notna().sum(),'finite new tests',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--wild',type=int,default=99999)
    parser.add_argument('--bootstrap',type=int,default=2000);execute(parser.parse_args())
