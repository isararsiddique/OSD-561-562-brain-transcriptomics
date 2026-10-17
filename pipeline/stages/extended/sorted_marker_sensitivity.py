from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Target-excluded marker gate within enriched myeloid samples; no pure-cell claim."""
from pathlib import Path
import gzip,json
import numpy as np,pandas as pd
from scipy.io import mmread
B=stage_root('extended');V1=stage_root('core')
md=pd.read_csv(B/'tables/GSE295428_sorted_cell_QC.csv');ann=pd.read_csv(B/'prepared/GSE295428_sorted_annotation.csv')
panel=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('')
gate=['Csf1r','Tmem119','P2ry12','Hexb','Sall1'];assert not set(gate)&set(panel.mouse_symbol_primary)
marker_context=gate+['Lyz2','Mrc1','Lyve1','Pf4','Ccr2','Ly6c2','S100a8','S100a9','Snap25','Rbfox3','Syt1','Aqp4','Gfap','Pecam1','Cldn5','Mbp','Plp1']
pb={s:np.zeros(len(ann),np.int64) for s in md['sample']};nmarker={};context=[]
for r in md.itertuples():
    if not r.tissue.startswith('Myeloid_'):nmarker[r.sample]=0;continue
    p=next((B/'sources/GSE295428_sorted_cells').glob(r.sample+'*_matrix.mtx.gz'))
    with gzip.open(p,'rb') as f:X=mmread(f).tocsc()
    nf=np.diff(X.indptr);lib=np.asarray(X.sum(axis=0)).ravel();mt=np.asarray(X[ann.SYMBOL.str.lower().str.startswith('mt-').to_numpy()].sum(axis=0)).ravel()/np.maximum(lib,1)
    qc=(nf>=200)&(nf<=2500)&(mt<=.5);ids=ann.index[ann.SYMBOL.isin(gate)].to_numpy();assert len(ids)==5
    positive=np.asarray((X[ids]>0).sum(axis=0)).ravel()>=2;keep=qc&positive
    pb[r.sample]=np.asarray(X[:,keep].sum(axis=1)).ravel();nmarker[r.sample]=int(keep.sum())
    total=np.asarray(X[:,qc].sum(axis=1)).ravel();cpm=total/total.sum()*1e6
    for g in marker_context:
        ids=ann.index[ann.SYMBOL.eq(g)]
        if len(ids)==1:context.append({'sample':r.sample,'donor':r.donor,'tissue':r.tissue,'flight':r.flight,'age':r.age,'gene':g,'log2CPM_plus1':float(np.log2(cpm[ids[0]]+1)),'marker_positive_fraction':keep.sum()/qc.sum()})
    print(r.sample,r.tissue,r.donor,int(qc.sum()),int(keep.sum()),flush=True)
md['n_marker_positive_cells']=md['sample'].map(nmarker);md.to_csv(B/'tables/GSE295428_sorted_cell_QC.csv',index=False)
pd.DataFrame(pb,index=ann.ENSEMBL).to_csv(B/'prepared/GSE295428_sorted_marker_positive_pseudobulk.csv.gz')
pd.DataFrame(context).to_csv(B/'tables/sorted_myeloid_marker_context.csv',index=False)
(B/'tables/sorted_marker_gate_definition.json').write_text(json.dumps({'marker_genes':gate,'all_excluded_from_151_CNS_panel':True,'gate':'At least two of five marker genes detected within the published QC gate; restrict original FACS Myeloid libraries.','role':'Post-result sensitivity analysis, not a validated microglial annotation. Marker regulation or capture can affect selection. Does not prove pure cells or causal within-cell regulation.'},indent=2))
