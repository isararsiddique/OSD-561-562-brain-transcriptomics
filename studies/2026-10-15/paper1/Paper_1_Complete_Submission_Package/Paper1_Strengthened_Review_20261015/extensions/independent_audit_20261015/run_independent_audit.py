#!/usr/bin/env python3
"""Read-only scientific audit. Only outputs alongside this script are written."""
from pathlib import Path
import hashlib, itertools, json, sys
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from run_revision import PANELS, load_dataset   # used only for independent comparisons


def adjust(p, by=False):
    p=np.asarray(p,float);ans=np.full(len(p),np.nan);ids=np.flatnonzero(np.isfinite(p));order=ids[np.argsort(p[ids])]
    mult=sum(1/np.arange(1,len(ids)+1)) if by else 1.
    ans[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)*mult/np.arange(1,len(order)+1))[::-1])[::-1])
    return ans


def read_source(n):
    raw=pd.read_csv(ROOT/f'inputs/OSD-{n}_counts.csv.gz',index_col=0)
    md=pd.read_csv(ROOT/f'inputs/OSD-{n}_metadata.csv.gz',dtype=str).set_index('sample').loc[raw.columns]
    animals=sorted(md['Source Name'].unique())
    collapsed=pd.DataFrame({a:raw.loc[:,md.index[md['Source Name'].eq(a)]].sum(axis=1) for a in animals})
    records=[]
    for a in animals:
        row=md.loc[md['Source Name'].eq(a)].iloc[0]
        if n in ('561','562','613'):
            flight='FLT' if '_FL_' in a else 'GC';cohort='ISS_T' if '_ISS_' in a else 'LAR';age=str(row['Factor Value[Age]'])
        elif n=='525':
            flight={'F':'FLT','G':'HGC','V':'VGC'}[a.split('_')[-1][0]];cohort=flight;age='single age'
        else:flight='FLT' if '_FL-' in a else 'GC';cohort='matched collection';age='single age'
        location=row.get('Factor Value[Euthanasia Location]','');expected='On ISS' if flight=='FLT' and cohort=='ISS_T' else 'On Earth'
        records.append(dict(animal=a,flight=flight,cohort=cohort,age=age,location_conflict=bool(pd.notna(location) and n in ('561','562','613') and location!=expected),n_libraries=int(md['Source Name'].eq(a).sum())))
    meta=pd.DataFrame(records).set_index('animal',drop=False)
    return raw,collapsed,meta,md


ANNOT=pd.read_csv(ROOT/'inputs/gene_annotation.csv.gz').set_index('ensembl_id').gene_symbol.to_dict()

def panel_indices(genes):
    symbols=np.array([ANNOT.get(g,'') for g in genes])
    return {name:np.flatnonzero(np.isin(symbols,syms)) for name,syms in PANELS.items()}


def expression_reference(c,genes,train,test):
    # Explicit sorted medians; training data determine both filter and reference.
    keep=np.array([str(g).startswith('ENSMUSG') and str(g)[7:].isdigit() for g in genes]) & (c[:,train].mean(axis=1)>=10)
    filtered=c[keep];genes=np.asarray(genes)[keep]
    positive=(filtered[:,train]>0).all(axis=1)
    lr=np.log(filtered[positive][:,train]).mean(axis=1)
    cols=np.r_[train,test];factors=[]
    for col in cols:
        nz=filtered[positive,col]>0
        ordered=np.sort(np.log(filtered[positive,col][nz])-lr[nz]);m=len(ordered)
        med=ordered[m//2] if m%2 else (ordered[m//2-1]+ordered[m//2])/2
        factors.append(np.exp(med))
    y=np.log2(filtered[:,cols]/np.asarray(factors)+1)
    return genes,y[:,:len(train)],y[:,len(train):],filtered[:,train].mean(axis=1),np.asarray(factors)


def provenance():
    checks=[];account=[];allanimals=set();core=[]
    for file in ['provenance.json','upstream_count_provenance.json']:
        for rec in json.loads((ROOT/'inputs'/file).read_text()):
            path=ROOT/'inputs'/rec.get('bundled_path',rec.get('bundled_upstream_path'))
            expected=rec.get('bundled_sha256',rec.get('bundled_upstream_sha256'))
            checks.append(dict(check='hash',file=str(path.relative_to(ROOT)),passed=hashlib.sha256(path.read_bytes()).hexdigest()==expected))
    data={}
    for n in ('525','561','562','613','563','564','612'):
        raw,c,m,md=read_source(n);data[n]=(raw,c,m,md)
        prod=load_dataset(n)
        a=prod['counts'].copy();a.columns=[x.split(':',1)[1] for x in a.columns]
        diff=float(np.max(abs(a.loc[c.index,c.columns].to_numpy()-c.to_numpy())))
        valid=bool(np.isfinite(raw.to_numpy()).all() and (raw.to_numpy()>=0).all() and np.equal(raw.to_numpy(),np.round(raw.to_numpy())).all())
        checks.append(dict(check='integer_nonnegative_and_exact_collapse',file=f'OSD-{n}',passed=valid and diff==0,max_difference=diff))
        if n in ('561','562'):
            up=pd.read_csv(ROOT/f'inputs/upstream_OSD-{n}_STAR.csv.gz',index_col=0)
            checks.append(dict(check='official_STAR_equality',file=f'OSD-{n}',passed=bool(np.array_equal(up.loc[raw.index,raw.columns].to_numpy(),raw.to_numpy()))))
            core.append(m)
        account.append(dict(accession=f'OSD-{n}',libraries=raw.shape[1],animals=len(m),raw_genes=len(raw),technical_extra=raw.shape[1]-len(m)))
        allanimals.update(m.animal)
    pd.DataFrame(account).to_csv(OUT/'independent_dataset_accounting.csv',index=False)
    # Recompute complete gene-family adjustments from component p values.
    for path,pcol,qcol in [('results/regional_gene_replication.csv.gz','replication_max_p','q_BH_global_replication_family'),('extensions/flight_collection_interaction/gene_flight_collection_interactions_exploratory.csv.gz','pvalue','q_BH_complete_gene_interaction_family')]:
        tab=pd.read_csv(ROOT/path);delta=float(np.nanmax(abs(adjust(tab[pcol])-tab[qcol])))
        checks.append(dict(check='independent_BH',file=path,passed=delta<1e-12,max_difference=delta,n_evaluable=int(tab[pcol].notna().sum()),n_q05=int(tab[qcol].le(.05).sum())))
    print('Source audit:',len(allanimals),'unique animals;',len(checks),'checks',flush=True)
    return data,checks,len(allanimals)


def exact_ground(data,checks):
    summaries=[]
    for n in ('561','562'):
        _,counts,meta,_=data[n];meta=meta.loc[meta.flight.eq('GC')].copy();counts=counts.loc[:,meta.index]
        labels=meta.cohort.eq('ISS_T').to_numpy(int);age=meta.age.to_numpy();ageids=[np.flatnonzero(age==a) for a in sorted(set(age))]
        assert all(len(i)==6 and labels[i].sum()==3 for i in ageids)
        assignments=[]
        for a,b in itertools.product(itertools.combinations(ageids[0],3),itertools.combinations(ageids[1],3)):
            g=np.zeros(len(meta),int);g[list(a)+list(b)]=1;assignments.append(g)
        assignments=np.asarray(assignments);names=list(PANELS)+['WholeGenome']
        values=np.zeros((len(assignments),len(names)));nsplit=np.zeros(len(assignments),int)
        c=counts.to_numpy(float);genes=counts.index.to_numpy()
        for sid,(a,b) in enumerate(itertools.product(itertools.combinations(ageids[0],2),itertools.combinations(ageids[1],2))):
            test=np.array(a+b);train=np.setdiff1d(np.arange(len(meta)),test)
            gs,yt,ye,_,_=expression_reference(c,genes,train,test);panels=panel_indices(gs)
            gram=np.stack([yt[ix].T@ye[ix] for ix in list(panels.values())+[np.arange(len(gs))]])
            valid=(assignments[:,test[:2]].sum(axis=1)==1)&(assignments[:,test[2:]].sum(axis=1)==1);ids=np.flatnonzero(valid)
            orientation=2*assignments[ids][:,train]-1
            score=np.einsum('nt,ptu->npu',orientation,gram)/4
            # Train8 is balanced within each age; nuisance age contributions cancel.
            for pair in (np.array([0,1]),np.array([2,3])):
                pos=np.argmax(assignments[ids][:,test[pair]],axis=1);neg=1-pos
                delta=score[np.arange(len(ids))[:,None],np.arange(len(names))[None,:],pair[pos][:,None]]-score[np.arange(len(ids))[:,None],np.arange(len(names))[None,:],pair[neg][:,None]]
                values[ids]+=(delta>1e-12)+.5*(abs(delta)<=1e-12)
            nsplit[ids]+=1
            if sid%75==0:print(f'OSD-{n} exact null preprocessing {sid+1}/225',flush=True)
        assert np.all(nsplit==81);values/=2*nsplit[:,None]
        observed=np.flatnonzero((assignments==labels).all(axis=1));assert len(observed)==1
        obs=values[observed[0]];p=(values>=obs[None,:]-1e-12).mean(axis=0)
        table=pd.DataFrame(values,columns=names);table.insert(0,'assignment_id',np.arange(len(table)));table['is_observed']=np.arange(len(table))==observed[0];table['labels']=''
        table['labels']=[''.join(map(str,g)) for g in assignments];table.to_csv(OUT/f'OSD-{n}_exact_label_null.csv',index=False)
        old=pd.read_csv(ROOT/'extensions/ground_holdout/holdout_summary.csv').query('method == "raw" and accession == @("OSD-"+n)') if False else pd.read_csv(ROOT/'extensions/ground_holdout/holdout_summary.csv')
        old=old[(old.method=='raw')&(old.accession==f'OSD-{n}')].iloc[0]
        for name,metric in [('Sensor31','sensor_age_stratified_AUC_mean'),('WholeGenome','whole_genome_age_stratified_AUC_mean')]:
            delta=float(obs[names.index(name)]-old[metric]);checks.append(dict(check='ground_holdout_AUC_independent',file=f'OSD-{n}:{name}',passed=abs(delta)<1e-12,max_difference=abs(delta)))
        for i,name in enumerate(names):summaries.append(dict(accession=f'OSD-{n}',panel_or_comparator=name,mean_crossfitted_age_stratified_AUC=obs[i],p_exact_complete_label_enumeration=p[i],n_label_assignments=len(assignments),n_overlapping_holdouts_per_assignment=81,n_ground_animals=len(meta),null_AUC_mean=values[:,i].mean(),null_AUC_q975=np.quantile(values[:,i],.975),n_null_at_least_observed=int((values[:,i]>=obs[i]-1e-12).sum())))
        print('Exact ground',n,'observed:',dict(zip(names,obs)),'p:',dict(zip(names,p)),flush=True)
    out=pd.DataFrame(summaries);out['q_BH_complete_audit_18_tests']=adjust(out.p_exact_complete_label_enumeration);out['q_BY_complete_audit_18_tests']=adjust(out.p_exact_complete_label_enumeration,True);out.to_csv(OUT/'exact_ground_sensor_association.csv',index=False)
    return out


def age_effects(data,checks):
    rows=[];means=[];full=[]
    conflict=set(pd.concat([data[n][2] for n in ('561','562')]).query('location_conflict').animal)
    for n in ('561','562'):
        _,c,m,_=data[n]
        # Match the fixed full-deposit normalization used in original interaction.
        gs,y,_,_,_=expression_reference(c.to_numpy(float),c.index.to_numpy(),np.arange(len(m)),np.array([],int));panels=panel_indices(gs)
        yp=np.stack([y[ix].mean(axis=0) for ix in panels.values()])
        for sensitivity in ['all_animals','exclude_any_region_metadata_conflict']:
            take=np.ones(len(m),bool) if sensitivity=='all_animals' else ~m.animal.isin(conflict).to_numpy();mm=m.iloc[np.flatnonzero(take)];yy=yp[:,take]
            f=mm.flight.eq('FLT').to_numpy(float);co=mm.cohort.eq('ISS_T').to_numpy(float);a=mm.age.eq('29').to_numpy(float)
            additive=np.column_stack([np.ones(len(mm)),a,f,co,f*co]);addbeta=np.linalg.lstsq(additive,yy.T,rcond=None)[0][-1]
            old=pd.read_csv(ROOT/'extensions/flight_collection_interaction/panel_flight_collection_interactions.csv')
            old=old[(old.accession==f'OSD-{n}')&(old.sensitivity==sensitivity)].set_index('panel')
            delta=float(np.max(abs(addbeta-old.loc[list(PANELS)].interaction_log2.to_numpy())))
            checks.append(dict(check='interaction_full_design_independent',file=f'OSD-{n}:{sensitivity}',passed=delta<1e-10,max_difference=delta))
            effects={}
            for ageval in ['12','29']:
                cellmeans={};sizes={}
                for fv,cv in itertools.product(['GC','FLT'],['LAR','ISS_T']):
                    cell=(mm.age==ageval)&(mm.flight==fv)&(mm.cohort==cv);sizes[(fv,cv)]=int(cell.sum())
                    cellmeans[(fv,cv)]=yy[:,cell].mean(axis=1) if cell.any() else np.full(len(PANELS),np.nan)
                    for j,p in enumerate(PANELS):means.append(dict(accession=f'OSD-{n}',sensitivity=sensitivity,panel=p,age=ageval,flight=fv,cohort=cv,n_animals=sizes[(fv,cv)],mean_panel_log_expression=cellmeans[(fv,cv)][j]))
                effect=cellmeans[('FLT','ISS_T')]-cellmeans[('GC','ISS_T')]-cellmeans[('FLT','LAR')]+cellmeans[('GC','LAR')];effects[ageval]=effect
                for j,p in enumerate(PANELS):rows.append(dict(accession=f'OSD-{n}',sensitivity=sensitivity,panel=p,age=ageval,direct_age_specific_interaction=effect[j],all_four_cells_represented=all(v>0 for v in sizes.values()),min_cell_n=min(sizes.values()),original_additive_age_interaction=addbeta[j]))
            for j,p in enumerate(PANELS):rows.append(dict(accession=f'OSD-{n}',sensitivity=sensitivity,panel=p,age='equal_age_standardized',direct_age_specific_interaction=(effects['12'][j]+effects['29'][j])/2,all_four_cells_represented=bool(np.isfinite(effects['12'][j]+effects['29'][j])),min_cell_n=np.nan,original_additive_age_interaction=addbeta[j]))
            if sensitivity=='all_animals':
                x=np.column_stack([np.ones(len(mm)),a,f,co,a*f,a*co,f*co,a*f*co]);bread=np.linalg.inv(x.T@x);fit=np.linalg.lstsq(x,yy.T,rcond=None)[0];r=yy.T-x@fit;h=np.einsum('ij,jk,ik->i',x,bread,x)
                for j,p in enumerate(PANELS):
                    sandwich=bread@x.T@np.diag(r[:,j]**2/(1-h)**2)@x@bread
                    contrasts={'young':np.array([0,0,0,0,0,0,1,0.]),'old':np.array([0,0,0,0,0,0,1,1.]),'equal_age_mean':np.array([0,0,0,0,0,0,1,.5]),'age_difference':np.array([0,0,0,0,0,0,0,1.])}
                    for target,v in contrasts.items():
                        b=float(v@fit[:,j]);se=float(np.sqrt(v@sandwich@v));full.append(dict(accession=f'OSD-{n}',panel=p,estimand=target,effect_log2=b,HC3_se=se,HC3_normal_lower95=b-1.96*se,HC3_normal_upper95=b+1.96*se,p_HC3_normal_approximate=2*stats.norm.sf(abs(b/se)),n_animals=len(mm),residual_df=len(mm)-8))
    pd.DataFrame(rows).to_csv(OUT/'age_specific_panel_interactions.csv',index=False);pd.DataFrame(means).to_csv(OUT/'age_design_cell_panel_means.csv',index=False)
    full=pd.DataFrame(full);full['q_BH_all_64_saturated_estimands']=adjust(full.p_HC3_normal_approximate);full.to_csv(OUT/'saturated_age_panel_interactions.csv',index=False)
    return pd.DataFrame(rows),full


def risk_decomposition(data,checks):
    splits=pd.read_csv(ROOT/'extensions/ground_holdout/split_animal_roles.csv');original=pd.read_csv(ROOT/'extensions/ground_holdout/holdout_split_results.csv');rows=[];animalrows=[];tertiles=[]
    for n in ('561','562','525'):
        _,c,m,_=data[n];mask=m.flight.eq('GC') if n!='525' else m.flight.isin(['HGC','VGC']);m=m.loc[mask];c=c.loc[:,m.index]
        labels=m.cohort.eq('ISS_T').to_numpy(float) if n!='525' else m.flight.eq('HGC').to_numpy(float)
        roles=splits[splits.accession==f'OSD-{n}']
        for sid,s in roles.groupby('split_id',sort=True):
            train=np.array([m.index.get_loc(a) for a in s[s.role=='train'].animal]);test=np.array([m.index.get_loc(a) for a in s[s.role=='test'].animal])
            # Sorting isn't statistically needed; this follows saved source order.
            gs,yt,ye,mean_c,sf=expression_reference(c.to_numpy(float),c.index.to_numpy(),train,test);panels=panel_indices(gs);ix=panels['Sensor31'];outside=np.setdiff1d(np.arange(len(gs)),ix)
            xt=np.column_stack([np.ones(len(train)),m.iloc[train].age.eq('29').to_numpy(float)]) if n!='525' else np.ones((len(train),1));xe=np.column_stack([np.ones(len(test)),m.iloc[test].age.eq('29').to_numpy(float)]) if n!='525' else np.ones((len(test),1))
            fit=np.linalg.lstsq(xt,yt.T,rcond=None)[0];y0=yt-(xt@fit).T;g0=labels[train]-xt@np.linalg.lstsq(xt,labels[train],rcond=None)[0];beta=y0@g0/(g0@g0)
            w=np.zeros(len(gs));w[ix]=beta[ix];w/=np.linalg.norm(w);scoretrain=w@y0;gamma=y0@scoretrain/(scoretrain@scoretrain);resid=ye-(xe@fit).T;score=w@resid;prediction=np.outer(gamma,score)
            r=resid[outside];p=prediction[outside];base=float(np.mean(r*r));predsq=float(np.mean(p*p));cross=float(np.mean(r*p));mse=float(np.mean((r-p)**2));ratio=mse/base
            old=original[(original.accession==f'OSD-{n}')&(original.split_id==sid)&(original.method=='score')].iloc[0]
            err=abs(ratio-old.outside_prediction_MSE_ratio_to_nuisance_only);checks.append(dict(check='holdout_risk_independent_reconstruction',file=f'OSD-{n}:split{sid}',passed=err<1e-10,max_difference=err))
            gene_base=np.mean(r*r,axis=1);gene_model=np.mean((r-p)**2,axis=1)
            rows.append(dict(accession=f'OSD-{n}',split_id=sid,n_training=len(train),n_test=len(test),n_outside_genes=len(outside),baseline_MSE=base,prediction_squared_MSE_component=predsq,twice_response_prediction_cross_component=2*cross,score_MSE=mse,score_MSE_ratio=ratio,frac_outside_genes_prediction_worse=float((gene_model>gene_base).mean()),decomposition_error=abs(mse-(base+predsq-2*cross))))
            # Tertiles are formed from training mean counts only.
            order=np.argsort(mean_c[outside],kind='stable')
            for j,idx in enumerate(np.array_split(order,3)):
                rj=r[idx];pj=p[idx];tertiles.append(dict(accession=f'OSD-{n}',split_id=sid,training_expression_tertile=j+1,n_genes=len(idx),mean_count_min=float(mean_c[outside][idx].min()),mean_count_max=float(mean_c[outside][idx].max()),baseline_MSE=float(np.mean(rj*rj)),score_MSE=float(np.mean((rj-pj)**2)),score_MSE_ratio=float(np.mean((rj-pj)**2)/np.mean(rj*rj)),frac_prediction_worse=float((np.mean((rj-pj)**2,axis=1)>np.mean(rj*rj,axis=1)).mean())))
            for j,i in enumerate(test):animalrows.append(dict(accession=f'OSD-{n}',split_id=sid,animal=m.iloc[i].animal,age=m.iloc[i].age,cohort=m.iloc[i].cohort,baseline_MSE=float(np.mean(r[:,j]**2)),score_MSE=float(np.mean((r[:,j]-p[:,j])**2))))
        print('Risk decomposition complete',n,flush=True)
    rows=pd.DataFrame(rows);rows.to_csv(OUT/'holdout_prediction_risk_decomposition.csv',index=False)
    a=pd.DataFrame(animalrows);a.to_csv(OUT/'heldout_animal_prediction_errors.csv',index=False);summ=[]
    for (acc,animal),part in a.groupby(['accession','animal']):summ.append(dict(accession=acc,animal=animal,age=part.age.iloc[0],cohort=part.cohort.iloc[0],n_repeated_holdouts=len(part),mean_baseline_MSE=part.baseline_MSE.mean(),mean_score_MSE=part.score_MSE.mean(),ratio_of_animal_mean_errors=part.score_MSE.mean()/part.baseline_MSE.mean()))
    pd.DataFrame(summ).to_csv(OUT/'heldout_animal_average_prediction_errors.csv',index=False);tertiles=pd.DataFrame(tertiles);tertiles.to_csv(OUT/'training_expression_tertile_prediction_errors.csv',index=False)
    return rows,pd.DataFrame(summ),tertiles


def main():
    data,checks,total=provenance();ground=exact_ground(data,checks);age,full=age_effects(data,checks);risk,animal,tertile=risk_decomposition(data,checks)
    checkframe=pd.DataFrame(checks);checkframe.to_csv(OUT/'independent_verification_checks.csv',index=False)
    failed=checkframe[~checkframe.passed]
    result=dict(status='PASS' if failed.empty else 'FAIL',n_checks=len(checkframe),n_failed=len(failed),n_unique_animals=total,protocol_sha256=hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest(),source_root=str(ROOT),versions={'numpy':np.__version__,'pandas':pd.__version__},max_holdout_ratio_discrepancy=float(checkframe.loc[checkframe.check=='holdout_risk_independent_reconstruction','max_difference'].max()))
    (OUT/'audit_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    print(ground.to_string(index=False),flush=True)
    print(age[(age.panel=='Circadian')].to_string(index=False),flush=True)
    print(risk.groupby('accession')[['score_MSE_ratio','frac_outside_genes_prediction_worse']].mean(),flush=True)
    print(animal.groupby('accession').ratio_of_animal_mean_errors.agg(['mean','median','min','max']),flush=True)
    if not failed.empty:raise AssertionError(failed.to_string(index=False))

if __name__=='__main__':main()
