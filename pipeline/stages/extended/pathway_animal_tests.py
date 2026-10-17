from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Positive-set singscore and exact whole-animal flight label tests."""
from pathlib import Path
import itertools,json
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
B=stage_root('extended');V1=stage_root('core')
sets={a[0]:a[2:] for line in (B/'sources/mh.all.v2025.1.Mm.symbols.gmt').read_text().splitlines() if (a:=line.split('\t'))}
rows=[];values=[]
def run(sid,X,ann,md,core=False,cell=False):
    # Identical >=10 in >=3 samples count gate to the count-model branch.
    X=X.loc[X.ge(10).sum(axis=1).ge(3)];N=len(X);symbols=ann.reindex(X.index).SYMBOL
    ranks=rankdata(X.to_numpy(),axis=0,method='average')
    for name,genes in sets.items():
        sel=symbols.isin(genes).to_numpy();m=int(sel.sum());score=(ranks[sel].mean(axis=0)-(N+1)/2)/(N-m) if m>=10 and m<N else np.full(X.shape[1],np.nan)
        for i,s in enumerate(X.columns):values.append({'study':sid,'set':name,'sample':s,'flight':md.loc[s,'flight'],'group':md.loc[s,'group'],'positive_set_singscore':score[i],'n_covered_genes':m,'n_background_genes':N})
        strata=[]
        if core:
            for cohort in ['ISS_T','LAR']:
                for age in ['12','29']:strata.append((f'flight_{cohort}_age{age}',md.group.isin([f'FLT_{cohort}_{age}',f'GC_{cohort}_{age}']).to_numpy()))
        else:strata=[('flight_vs_ground',np.ones(len(md),bool))]
        for tag,ix in strata:
            y=score[ix];f=md.iloc[np.flatnonzero(ix)].flight.eq('FLT').to_numpy();effect=float(y[f].mean()-y[~f].mean()) if m>=10 else np.nan
            null=[]
            if m>=10:
                for ids in itertools.combinations(range(len(y)),int(f.sum())):
                    mask=np.isin(np.arange(len(y)),ids);null.append(y[mask].mean()-y[~mask].mean())
            rows.append({'study':sid,'set':name,'contrast':tag,'n_flight':int(f.sum()),'n_ground':int((~f).sum()),'n_covered_genes':m,'singscore_effect':effect,'exact_p':float(np.mean(np.abs(null)>=abs(effect)-1e-12)) if null else 1.,'n_exact_labels':len(null),'family':'cells_250' if cell else 'core_flight_400' if core else 'independent_150'})
for sid in [561,562,352,563,564]:
    src=V1 if sid in [561,562] else B
    X=pd.read_csv(src/f'prepared/OSD-{sid}_STAR_endogenous.csv.gz',index_col=0)
    md=pd.read_csv(src/f'tables/OSD-{sid}_{"sample_metadata" if sid in [561,562] else "metadata"}.csv',index_col=0).loc[X.columns]
    ann=pd.read_csv(src/f'prepared/OSD-{sid}_annotation.csv').set_index('ENSEMBL')
    run(f'OSD-{sid}',X,ann,md,sid in [561,562])
ann=pd.read_csv(B/'prepared/OSD352_RNA_annotation.csv').set_index('ENSEMBL')
for level,cts in [('l2',['CNS neurons','CNS glia','Immune cells']),('l3',['Astroependymal cells','Oligodendrocytes'])]:
    X=pd.read_csv(B/f'prepared/OSD352_RNA_{"subtype_" if level=="l3" else ""}pseudobulk.csv.gz',index_col=0)
    md=pd.read_csv(B/f'tables/OSD352_{level}_registry.csv').set_index('column')
    for ct in cts:
        s=md[md.cell_type.eq(ct)].copy();s['flight']=s.flight_status.map({'flight':'FLT','ground':'GC'});s['group']=s.flight
        name='OSD352_'+level+'_'+''.join(c if c.isalnum() else '_' for c in ct)
        run(name,X.loc[:,s.index],ann,s,cell=True)
a=pd.DataFrame(rows);a['q_full_family']=a.groupby('family').exact_p.transform(false_discovery_control)
a.to_csv(B/'tables/pathway_animal_exact_tests.csv',index=False)
pd.DataFrame(values).to_csv(B/'tables/pathway_singscore_values.csv',index=False)
print(a.groupby('family').q_full_family.apply(lambda q:int(q.lt(.05).sum())).to_dict(),flush=True)
