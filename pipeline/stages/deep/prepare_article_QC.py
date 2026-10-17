from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Resolve article <=5% versus GEO <=50% mitochondrial threshold via sensitivity."""
from pathlib import Path
import gzip,json
import numpy as np,pandas as pd
from scipy.io import mmread
B=stage_root('deep');R=stage_root('extended')
md=pd.read_csv(R/'tables/GSE295428_current_sample_registry.csv');ann=pd.read_csv(R/'prepared/GSE295428_sorted_annotation.csv')
pb={};qc=[]
for s in md.itertuples():
    p=next((R/'sources/GSE295428_sorted_cells').glob(s.sample+'*_matrix.mtx.gz'))
    with gzip.open(p,'rb') as f:x=mmread(f).tocsc()
    lib=np.asarray(x.sum(axis=0)).ravel();nfeat=np.diff(x.indptr)
    mt=np.asarray(x[ann.SYMBOL.str.lower().str.startswith('mt-').to_numpy()].sum(axis=0)).ravel()/np.maximum(lib,1)
    keep=(nfeat>=200)&(nfeat<=2500)&(mt<=.05)
    pb[s.sample]=np.asarray(x[:,keep].sum(axis=1)).ravel()
    qc.append({'sample':s.sample,'donor':s.donor,'tissue':s.tissue,'flight':s.flight,'age':s.age,'n_cells':int(keep.sum()),'retained_UMIs':int(lib[keep].sum())})
    print(s.sample,s.tissue,int(keep.sum()),flush=True)
pd.DataFrame(pb,index=ann.ENSEMBL).to_csv(B/'prepared/article_5pct_donor_pseudobulk.csv.gz')
pd.DataFrame(qc).to_csv(B/'tables/article_5pct_mouse_QC.csv',index=False)
(B/'tables/article_GEO_QC_discrepancy.json').write_text(json.dumps({'article_source':'https://doi.org/10.1038/s41467-026-68737-1','article_methods':'Cells removed at <200 or >2500 detected features, or >5% mitochondrial RNA','GEO_source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295428','GEO_processing':'Threshold says >50% mitochondrial RNA','primary_analysis':'Completed primary gate follows deposited GEO processing statement; it is not an exact reproduction of article Methods','new_sensitivity':'200-2500 features, <=5% mitochondrial UMI; current attached Cell Ranger matrices, no claim to reproduce SoupX or DoubletFinder','n_libraries':len(qc),'n_mice':len(set(q['donor'] for q in qc)),'n_retained_cells':sum(q['n_cells'] for q in qc)},indent=2))
