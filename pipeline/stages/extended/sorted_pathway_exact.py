from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""RR8 complete Hallmark scores with age-stratified whole-donor inference."""
from pathlib import Path
import itertools
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
B=stage_root('extended')
sets={a[0]:a[2:] for l in (B/'sources/mh.all.v2025.1.Mm.symbols.gmt').read_text().splitlines() if (a:=l.split('\t'))}
X=pd.read_csv(B/'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz',index_col=0);ann=pd.read_csv(B/'prepared/GSE295428_sorted_annotation.csv').set_index('ENSEMBL');md=pd.read_csv(B/'tables/GSE295428_sorted_cell_QC.csv').set_index('sample')
rows=[];values=[]
for tissue in sorted(md.tissue.unique()):
    s=md[(md.tissue==tissue)&md.n_retained_cells.ge(50)];C=X.loc[:,s.index];C=C.loc[C.ge(10).sum(axis=1).ge(3)];N=len(C);ranks=rankdata(C.to_numpy(),axis=0,method='average');symbols=ann.loc[C.index].SYMBOL
    f=s.flight.eq('FLT').to_numpy();ages=[s.age.eq(a).to_numpy() for a in [3,8]]
    choices=[list(itertools.combinations(np.flatnonzero(ix),int(f[ix].sum()))) for ix in ages];masks=[np.isin(np.arange(len(s)),sum((list(v) for v in combo),[])) for combo in itertools.product(*choices)]
    def effects(y,mask):
        a,b=[y[mask&ix].mean()-y[(~mask)&ix].mean() for ix in ages];return np.array([a,b,(a+b)/2,b-a])
    for name,genes in sets.items():
        hit=symbols.isin(genes).to_numpy();m=int(hit.sum());ok=10<=m<N;y=(ranks[hit].mean(axis=0)-(N+1)/2)/(N-m) if ok else np.full(len(s),np.nan)
        observed=effects(y,f) if ok else np.full(4,np.nan);pv=np.mean(np.abs(np.array([effects(y,mask) for mask in masks]))>=np.abs(observed)-1e-12,axis=0) if ok else np.ones(4)
        for sample,score in zip(s.index,y):values.append({'tissue':tissue,'sample':sample,'donor':s.loc[sample,'donor'],'flight':s.loc[sample,'flight'],'age':s.loc[sample,'age'],'set':name,'singscore':score})
        for k,tag in enumerate(['flight_age3','flight_age8','flight_balanced_age','flight_by_age']):rows.append({'tissue':tissue,'set':name,'contrast':tag,'n_donors':len(s),'n_exact_labels':len(masks),'n_covered_genes':m,'score_effect':observed[k],'exact_p':pv[k],'status':'tested' if ok else 'insufficient_coverage'})
a=pd.DataFrame(rows);a['q_full_1200_family']=false_discovery_control(a.exact_p);a.to_csv(B/'tables/sorted_pathway_exact_tests.csv',index=False);pd.DataFrame(values).to_csv(B/'tables/sorted_pathway_donor_scores.csv',index=False)
print('RR8 exact pathway corrected hits',int(a.q_full_1200_family.lt(.05).sum()))
