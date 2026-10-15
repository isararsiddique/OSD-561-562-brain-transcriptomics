"""Independent comparator-transport reconstruction; imports no analysis routines."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import pandas as pd,numpy as np,json,itertools,math,hashlib,zipfile
B=Path(__file__).resolve().parents[1];P=B.parents[1];R=B/'results';checks=[]
def ck(name,ok,detail=''):
    checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name+': '+str(detail))
def read_counts(path,ids,junction_order,sizes):
    selected={};all_boundary=pd.DataFrame();sample_names=None;rows=0
    wanted=set(junction_order);frozen=set(ids)
    for c in pd.read_csv(path,index_col=0,chunksize=100000):
        if sample_names is None:sample_names=c.columns.tolist();all_boundary=pd.DataFrame(0,index=ids,columns=sample_names,dtype=np.int64)
        rows+=len(c);present=c.index.intersection(wanted)
        for j in present:selected[j]=c.loc[j].to_numpy(dtype=np.int64)
        bits=c.index.to_series().str.split(':',expand=True)
        for pos,suffix in [(2,':start'),(3,':end')]:
            key=bits[0]+':'+bits[1]+':'+bits[pos]+suffix;keep=key.isin(frozen)
            if keep.any():
                agg=c.loc[keep].groupby(key.loc[keep]).sum();all_boundary.loc[agg.index]+=agg
    alt=np.stack([selected.get(j,np.zeros(len(sample_names),dtype=np.int64)) for j in junction_order]);off=np.r_[0,np.cumsum(sizes)];tot=np.add.reduceat(alt,off[:-1],axis=0);den=np.repeat(tot,sizes,axis=0);fraction=np.divide(alt,den,out=np.zeros_like(alt,dtype=float),where=den>0);dom=np.divide(tot,all_boundary.to_numpy(),out=np.zeros_like(tot,dtype=float),where=all_boundary.to_numpy()>0)
    ck(path.stem+' selected <= boundary',np.all(tot<=all_boundary.to_numpy()))
    return dict(samples=sample_names,counts=alt,selected=tot,canonical=all_boundary.to_numpy(),P=fraction,dominance=dom,gate=(tot>=20).all(1)&(dom>=.9-1e-12).all(1),rows=rows)
def assignments(meta):
    obs=np.zeros(len(meta));choices=[];design=[]
    for st in sorted(meta.stratum.unique()):
        ix=np.flatnonzero(meta.stratum.values==st);nf=int((meta.iloc[ix].exposure=='FLT').sum());ng=len(ix)-nf;weight=len(ix)/len(meta)
        obs[ix]=np.where(meta.iloc[ix].exposure=='FLT',weight/nf,-weight/ng)
        choices.append([(tuple(sel),ix,nf,ng,weight) for sel in itertools.combinations(ix,nf)]);design.append((st,nf,ng,weight,math.comb(len(ix),nf)))
    D=[]
    for combo in itertools.product(*choices):
        v=np.zeros(len(meta))
        for sel,ix,nf,ng,weight in combo:v[ix]=-weight/ng;v[list(sel)]=weight/nf
        D.append(v)
    D=np.array(D);ck('allocation weights sum zero '+str(len(meta)),np.max(abs(D.sum(1)))<1e-14);ck('observed assignment included '+str(len(meta)),np.any(np.max(abs(D-obs),axis=1)<1e-14));return obs,D,design
def p_adjust(p,method):
    p=np.asarray(p);n=len(p);order=np.argsort(p,kind='stable');out=np.empty(n);q=[]
    if method=='holm':
        current=0
        for rank,i in enumerate(order):current=max(current,(n-rank)*p[i]);q.append(min(current,1))
    else:
        harmonic=sum(1/k for k in range(1,n+1)) if method=='by' else 1;current=1
        for rank in range(n-1,-1,-1):current=min(current,p[order[rank]]*n*harmonic/(rank+1));q.append(current)
        q=q[::-1]
    out[order]=q;return out

def reconstruct_everything():
    z=np.load(P/'extensions/multiway_splicing/tables/frozen_features.npz');ids=z['cluster_ids'];off=z['offsets'];sizes=np.diff(off)
    lookup=pd.read_csv(P/'extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz')[['cluster_id','alternative_index','junction_id']].drop_duplicates().set_index(['cluster_id','alternative_index'])
    order=[(cid,j) for cid,n in zip(ids,sizes) for j in range(n)];junctions=lookup.loc[order,'junction_id'].to_numpy();ck('frozen 923 clusters / 3376 occurrences',len(ids)==923 and len(junctions)==3376)
    data={};data['RR23']=read_counts(P/'inputs/all_animal_junction_counts.csv.gz',ids,junctions,sizes);data['RRRM2']=read_counts(P/'extensions/multiway_transport/inputs/RRRM2_animal_junction_counts.csv.gz',ids,junctions,sizes);data['RR10']=read_counts(P/'extensions/cohort_search/inputs/RR10/RR10_animal_junction_counts.csv.gz',ids,junctions,sizes)
    ck('raw source row counts',[data[c]['rows'] for c in ['RR23','RRRM2','RR10']]==[805562,668237,590325]);ck('RR23 all gates',data['RR23']['gate'].all());ck('frozen RR23 H squared reproduced',np.allclose(data['RR23']['P'].T,z['H']**2,atol=1e-14,rtol=0))
    rr=pd.read_csv(P/'inputs/animal_metadata.csv').set_index('sample').loc[data['RR23']['samples']];v0=pd.read_csv(P/'extensions/multiway_transport/inputs/RRRM2_animal_metadata.csv').set_index('sample');v=v0.loc[data['RRRM2']['samples']]
    r10=pd.read_csv(P/'extensions/cohort_search/inputs/RR10/RR10_animal_metadata.csv').rename(columns={'index':'sample'}).set_index('sample').loc[data['RR10']['samples']];r10['exposure']=r10.group;r10['stratum']='all'
    ck('RR23 groups',rr.group.value_counts().to_dict()==dict(FLT=9,GC=9,VIV=8));ck('RRRM2 groups',v.exposure.value_counts().to_dict()==dict(FLT=10,GC=9));ck('RR10 groups',r10.exposure.value_counts().to_dict()==dict(FLT=5,GC=5));ck('RRRM2 metadata exclusion',set(v.index)==set(v0.index[~((v0.cohort=='ISS_T')&(v0.age==29))]));ck('all distinct source animals',len(set(rr['Source Name'])|set(v['Source Name'])|set(r10['Source Name']))==55)
    ck('RR10 source animal and RFID no RR23 overlap',not set(r10['Source Name'])&set(rr['Source Name']) and not set(r10['Comment[RFID]'].astype(str))&set(rr['Comment[RFID]'].astype(str)))
    zmeta=zipfile.ZipFile(P/'extensions/cohort_search/inputs/RR10/OSD-563_metadata_OSD-563-ISA.zip');raw=pd.read_csv(zmeta.open('s_OSD-563.txt'),sep='\t').set_index('Sample Name').loc[r10.index]
    for col in ['Source Name','Comment[Euthanasia Date]','Comment[BSP Dissection Date]','Characteristics[Sex]','Characteristics[Strain]','Factor Value[Spaceflight]','Comment[RFID]']:ck('RR10 original ISA '+col,np.array_equal(raw[col].astype(str),r10[col].astype(str)))
    ck('RR10 euthanasia date alias',np.array_equal(r10.exposure=='FLT',r10['Comment[Euthanasia Date]']=='02-Jan-2021'));ck('RR10 dissection date alias',np.array_equal(r10.exposure=='FLT',r10['Comment[BSP Dissection Date]']=='12-May-2021'))
    for cohort,meta,expected in [('RRRM2',v,11200),('RR10',r10,252)]:
        data[cohort]['meta']=meta;obs,D,design=assignments(meta);data[cohort].update(obs=obs,D=D,design=design);ck(cohort+' exact assignment count',len(D)==expected)
    ck('RRRM2 three block allocations',[x[-1] for x in data['RRRM2']['design']]==[56,10,20]);ck('RRRM2 previous common count',data['RRRM2']['gate'].sum()==281)
    for comp in ['GC','VIV']:
        meanf=data['RR23']['P'][:,rr.group=='FLT'].mean(1);meanref=data['RR23']['P'][:,rr.group==comp].mean(1);delta=meanf-meanref;norm=np.sqrt(np.add.reduceat(delta*delta,off[:-1]));den=np.repeat(norm,sizes);u=np.divide(delta,den,out=np.zeros_like(delta),where=den>1e-12)
        data['RR23'][comp]=dict(delta=delta,norm=norm,reference=meanref,unit=u)
        ck('training unit norm '+comp,np.allclose(np.add.reduceat(u*u,off[:-1]),norm>1e-12,atol=1e-12));ck('training zero-sum directions '+comp,np.max(abs(np.add.reduceat(u,off[:-1])))<1e-10)
        for cohort in ['RRRM2','RR10']:
            d=data[cohort];scores=np.add.reduceat(u[:,None]*(d['P']-meanref[:,None]),off[:-1],axis=0);uncentered=np.add.reduceat(u[:,None]*d['P'],off[:-1],axis=0);ck(cohort+comp+' centering cancels all local effects',np.allclose(scores@d['obs'],uncentered@d['obs'],atol=1e-14))
            d[comp]=dict(scores=scores)
            for panel,keep in [('primary',d['gate']),('intersection',data['RRRM2']['gate']&data['RR10']['gate'])]:
                avg=scores[keep].mean(0);observed=avg@d['obs'];null=d['D']@avg;p=np.mean(null>=observed-1e-12);p2=np.mean(abs(null)>=abs(observed)-1e-12)
                d[comp][panel]=dict(score=avg,effect=observed,null=null,p=p,p2=p2,keep=keep)
            effects=scores@d['obs'];lp=np.ones(len(ids));lp2=np.ones(len(ids))
            for j in np.flatnonzero(d['gate']):
                null=d['D']@scores[j];lp[j]=np.mean(null>=effects[j]-1e-12);lp2[j]=np.mean(abs(null)>=abs(effects[j])-1e-12)
            d[comp].update(local_effect=effects,local_p=lp,local_p2=lp2)
    ck('previous RRRM2 FLT-GC one-sided p reproduced',data['RRRM2']['GC']['primary']['p']==.73625)
    return data,ids,off,junctions,rr

def bootweights(meta,rng):
    W=np.zeros((4000,len(meta)))
    for st in sorted(meta.stratum.unique()):
        ix=np.flatnonzero(meta.stratum.values==st);prop=len(ix)/len(meta)
        for label,sign in [('FLT',1),('GC',-1)]:
            sub=ix[meta.iloc[ix].exposure.values==label];counts=rng.multinomial(len(sub),np.repeat(1/len(sub),len(sub)),size=4000);W[:,sub]=sign*prop*counts/len(sub)
    return W

def verify_outputs(data,ids,off,junctions,rr):
    gs=pd.read_csv(R/'all_comparator_global_transport.csv');ls=pd.read_csv(R/'all3692_local_directional_tests.csv.gz');ns=pd.read_csv(R/'all_exact_global_nulls.csv.gz');sc=pd.read_csv(R/'all_animal_scores.csv');bs=pd.read_csv(R/'all_bootstrap_draws.csv.gz');tr=pd.read_csv(R/'RR23_only_frozen_directions.csv.gz');summary=json.loads((R/'summary.json').read_text());sizes=np.diff(off)
    ck('complete eight global rows',len(gs)==8 and not gs.duplicated(['cohort','direction','panel']).any());ck('complete 3692 local slots',len(ls)==3692 and not ls.duplicated(['cohort','direction','cluster_id']).any());ck('all 923 clusters in each local family',ls.groupby(['cohort','direction']).size().eq(923).all())
    for meth in ['holm','by','bh']:
        vals=p_adjust(ls.exact_p,meth);ck('joint 3692 '+meth,np.allclose(vals,ls[meth.upper()+'_joint3692'],atol=1e-13))
    ck('joint local Holm unattainable',3692/11200>.05 and (ls.HOLM_joint3692>.05).all())
    for panel,tab in gs.groupby('panel'):ck('four-test global Holm '+panel,np.allclose(p_adjust(tab.exact_p,'holm'),tab.Holm_four_global_p,atol=1e-13))
    common=data['RRRM2']['gate']&data['RR10']['gate'];table=pd.read_csv(R/'common_both_validation_panel.csv').set_index('cluster_id').loc[ids];ck('intersection panel exactly label-free overlap',np.array_equal(common,table.in_both_validation_cohorts))
    for comp in ['GC','VIV']:
        expected=data['RR23'][comp];tab=tr[tr.direction=='FLT_minus_'+comp].set_index(['cluster_id','alternative_index']).loc[[(c,j) for c,n in zip(ids,sizes) for j in range(n)]]
        for col,k in [('RR23_native_delta','delta'),('RR23_unit_vector','unit'),('RR23_comparator_mean','reference')]:ck('training '+comp+' '+k,np.allclose(expected[k],tab[col],atol=1e-13,rtol=0))
        ck('training alternative membership '+comp,np.array_equal(junctions,tab.junction_id))
    for ci,cohort in enumerate(['RRRM2','RR10']):
        d=data[cohort];meta=d['meta'];gates=pd.read_csv(R/f'{cohort}_all923_measurement_gates.csv').set_index('cluster_id').loc[ids];frozen=np.load(R/f'{cohort}_frozen_measurements.npz');ck(cohort+' raw sample order',np.array_equal(frozen['samples'],d['samples']))
        for col,k in [('counts','counts'),('selected','selected'),('boundary','canonical')]:ck(cohort+' raw measurement '+col,np.array_equal(frozen[col],d[k]))
        ck(cohort+' positive-total proportions',np.allclose(frozen['proportions'][np.isfinite(frozen['proportions'])],d['P'][np.isfinite(frozen['proportions'])],atol=1e-14));ck(cohort+' only zero totals give undefined proportions',np.array_equal(~np.isfinite(frozen['proportions']),np.repeat(d['selected']==0,sizes,axis=0)))
        ck(cohort+' blind gates',np.array_equal(gates.common_measurable,d['gate']));ck(cohort+' gate reads',np.array_equal(gates.minimum_selected_reads,d['selected'].min(1)));ck(cohort+' gate dominance',np.allclose(gates.minimum_dominance,d['dominance'].min(1)))
        actualmeta=pd.read_csv(R/f'{cohort}_included_animals.csv').set_index('sample').loc[d['samples']];ck(cohort+' frozen source identities',np.array_equal(actualmeta['Source Name'],meta['Source Name']));ck(cohort+' frozen exposure and blocks',np.array_equal(actualmeta[['exposure','stratum']],meta[['exposure','stratum']]))
        ds=pd.read_csv(R/f'{cohort}_allocation_design.csv');ck(cohort+' exported assignments',ds.n_assignments.prod()==len(d['D']))
        rng=np.random.default_rng(2026101541+ci);W=bootweights(meta,rng)
        for comp in ['GC','VIV']:
            direction='FLT_minus_'+comp;lc=ls[(ls.cohort==cohort)&(ls.direction==direction)].set_index('cluster_id').loc[ids];expected=d[comp];ck(cohort+comp+' all local one-sided exact p',np.array_equal(np.rint(lc.exact_p*len(d['D'])).astype(int),np.rint(expected['local_p']*len(d['D'])).astype(int)) and np.allclose(lc.exact_p,expected['local_p'],atol=2e-16,rtol=0));ck(cohort+comp+' all local two-sided exact p',np.allclose(lc.exact_two_sided_p,expected['local_p2'],atol=1e-16));ck(cohort+comp+' local measurable effects',np.allclose(lc.loc[d['gate'],'effect'],expected['local_effect'][d['gate']],atol=1e-13));ck(cohort+comp+' local unmeasurable bookkeeping',lc.loc[~d['gate'],'effect'].isna().all() and lc.loc[~d['gate'],'exact_p'].eq(1).all())
            train_ix=np.flatnonzero(rr.group.isin(['FLT',comp]));tm=rr.iloc[train_ix].copy();tm['exposure']=np.where(tm.group=='FLT','FLT','GC');tm['stratum']='all';WT=bootweights(tm,rng)
            zw=np.load(R/f'{cohort}_{direction}_bootstrap_weights.npz')
            ck(cohort+comp+' saved validation bootstrap weights',np.allclose(zw['validation_weights'],W,atol=1e-15,rtol=0));ck(cohort+comp+' saved training bootstrap weights',np.allclose(zw['training_weights'],WT,atol=1e-15,rtol=0));ck(cohort+comp+' training indices and source samples',np.array_equal(zw['training_indices'],train_ix) and np.array_equal(zw['training_samples'],rr.index[train_ix]));ck(cohort+comp+' validation sample identities',np.array_equal(zw['validation_samples'],d['samples']));ck(cohort+comp+' seed',int(zw['seed'])==2026101541+ci)
            train_delta=WT@data['RR23']['P'][:,train_ix].T;validation_delta=W@d['P'].T;norm=np.sqrt(np.add.reduceat(train_delta**2,off[:-1],axis=1));numer=np.add.reduceat(train_delta*validation_delta,off[:-1],axis=1);jointlocal=np.divide(numer,norm,out=np.zeros_like(numer),where=norm>1e-12)
            for pname,key in [('cohort_measured','primary'),('intersection_both_cohorts','intersection')]:
                e=expected[key];row=gs[(gs.cohort==cohort)&(gs.direction==direction)&(gs.panel==pname)].iloc[0];label=cohort+comp+pname
                ck(label+' panel size',row.n_clusters==e['keep'].sum());ck(label+' effect',abs(row.effect_pp-100*e['effect'])<1e-11);ck(label+' exact p',abs(row.exact_p-e['p'])<1e-15 and abs(row.exact_two_sided_p-e['p2'])<1e-15);ck(label+' animal and assignment counts',row.n_animals==len(meta) and row.n_assignments==len(d['D']))
                nn=ns[(ns.cohort==cohort)&(ns.direction==direction)&(ns.panel==pname)].sort_values('assignment');ck(label+' complete exact null',len(nn)==len(d['D']) and np.allclose(nn.effect_pp,100*e['null'],atol=1e-12))
                ss=sc[(sc.cohort==cohort)&(sc.direction==direction)&(sc.panel==pname)].set_index('sample').loc[d['samples']];ck(label+' all frozen animal scores',np.allclose(ss.score,e['score'],atol=1e-13));ck(label+' contrast weights',np.allclose(ss.contrast_weight,d['obs'],atol=1e-15))
                bb=bs[(bs.cohort==cohort)&(bs.direction==direction)&(bs.panel==pname)].sort_values('draw');fixed=W@e['score'];joint=jointlocal[:,e['keep']].mean(1);ck(label+' all 4000 fixed bootstrap draws',len(bb)==4000 and np.allclose(bb.conditional_effect_pp,100*fixed,atol=1e-11));ck(label+' all 4000 retrained bootstrap draws',np.allclose(bb.two_cohort_effect_pp,100*joint,atol=1e-11));ck(label+' fixed interval quantiles',np.allclose(np.quantile(fixed,[.025,.975])*100,[row.conditional_CI95_low_pp,row.conditional_CI95_high_pp],atol=1e-11));ck(label+' retrained interval quantiles',np.allclose(np.quantile(joint,[.025,.975])*100,[row.two_cohort_CI95_low_pp,row.two_cohort_CI95_high_pp],atol=1e-11))
    tc=pd.read_csv(R/'Tcaf1_all_comparators_and_cohorts.csv');expected=ls[ls.cluster_id=='6:2:42663893:start'];ck('Tcaf1 all four slots retained',len(tc)==4 and tc.equals(expected.reset_index(drop=True)))
    for row in pd.read_csv(R/'frozen_input_hashes.csv').itertuples():ck('source hash '+row.file,hashlib.sha256((P/row.file).read_bytes()).hexdigest()==row.sha256 and (P/row.file).stat().st_size==row.bytes)
    ck('plan hash',hashlib.sha256((B/'RETROSPECTIVE_PLAN_2026-10-15.md').read_bytes()).hexdigest()==summary['retrospective_plan_sha256']);ck('summary overlap counts',[summary['RRRM2_common_clusters'],summary['RR10_common_clusters'],summary['common_both_validation']]==[int(data['RRRM2']['gate'].sum()),int(data['RR10']['gate'].sum()),int(common.sum())]);ck('summary flags preserve interpretation',summary['RR10_shares_discovery_launch'] and summary['RR10_distinct_source_animals'] and summary['RR10_exposure_date_alias']);ck('summary discovery counts',summary['primary_global_Holm_calls']==int(((gs.panel=='cohort_measured')&(gs.Holm_four_global_p<=.05)).sum()) and summary['local_BH_calls']==int((ls.BH_joint3692<=.05).sum()) and summary['local_BY_calls']==int((ls.BY_joint3692<=.05).sum()))
    return dict(status='passed',independent_checks=len(checks),failed=0,checks=checks,statistical_limits=['No exposure-isolated causal inference: RR10 collection/dissection dates are aliased with exposure.','RR10 uses distinct source animals but shares the discovery launch and is not independent-launch replication.','All local Holm discovery is unattainable in the 3692-slot family: 3692/11200 = '+str(3692/11200)+'.','Retrospectively added comparator and cohort families have no project-wide confirmatory error guarantee.','Independent implementation recreated every primary/intersection exact allocation, local p value, animal score and all fixed/retrained bootstrap draws.'])

if __name__=='__main__':
    try:
        data,ids,off,junctions,rr=reconstruct_everything();result=verify_outputs(data,ids,off,junctions,rr)
    except Exception as exc:
        result=dict(status='failed',independent_checks=len(checks),failed=sum(not c['passed'] for c in checks),error=repr(exc),checks=checks);(R/'independent_verification.json').write_text(json.dumps(result,indent=2));raise
    (R/'independent_verification.json').write_text(json.dumps(result,indent=2));print('PASSED',len(checks),'independent transport checks')
