from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Post-result distribution-free gene sensitivity, keeping animals as n."""
from pathlib import Path
import itertools
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
OUT=stage_root('core');TAB=OUT/'tables';PRE=OUT/'prepared'
bulk=pd.read_csv(TAB/'bulk_panel_complete.csv');rows=[]
for osd in [561,562]:
    norm=pd.read_csv(PRE/f'OSD-{osd}_normalized_counts.csv.gz',index_col=0)
    meta=pd.read_csv(TAB/f'OSD-{osd}_sample_metadata.csv',index_col=0)
    for tag,g in bulk[(bulk.study==f'OSD-{osd}') & (bulk.contrast_kind=='flight')].groupby('contrast'):
        num,den=g.iloc[0][['numerator','denominator']];sub=meta[meta.group.isin([num,den])]
        f=sub.group.eq(num).to_numpy();assign=list(itertools.combinations(range(len(sub)),int(f.sum())))
        masks=np.array([np.isin(np.arange(len(sub)),x) for x in assign])
        w=masks/f.sum()-(~masks)/(~f).sum()
        for r in g.itertuples():
            row={'study':r.study,'contrast':tag,'analyte':r.analyte,'mouse_symbol':r.mouse_symbol,
                 'n_flight':int(f.sum()),'n_ground':int((~f).sum()),'n_assignments':len(assign),
                 'log_normalized_effect':np.nan,'exact_p':1.,'status':'unavailable'}
            if r.gene_id in norm.index:
                y=np.log2(norm.loc[r.gene_id,sub.index].to_numpy()+1)
                effect=y[f].mean()-y[~f].mean();null=w@y
                row.update(log_normalized_effect=effect,exact_p=float(np.mean(abs(null)>=abs(effect)-1e-12)),status='exact_conditional_animal_sensitivity')
            rows.append(row)
out=pd.DataFrame(rows);assert len(out)==1208
out['q_exact_flight_family']=false_discovery_control(out.exact_p)
out.to_csv(TAB/'bulk_gene_exact_sensitivity.csv',index=False)
print('Exact gene sensitivity q<.05',out.q_exact_flight_family.lt(.05).sum())
