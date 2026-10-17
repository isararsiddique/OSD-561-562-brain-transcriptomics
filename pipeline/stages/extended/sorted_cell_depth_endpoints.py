from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Mouse summaries of detection and depth-standardized expression; no cell-level exposure tests."""
from pathlib import Path
import gzip,json
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy.special import gammaln
B=stage_root('extended');V1=stage_root('core')
panel=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('');md=pd.read_csv(B/'tables/GSE295428_sorted_cell_QC.csv');md=md[md.tissue.str.startswith('Myeloid')]
rows=[];quality=[]
for sample in md.itertuples():
    folder=B/'sources/GSE295428_sorted_cells';features=next(folder.glob(sample.sample+'_*_features.tsv.gz'));matrix=next(folder.glob(sample.sample+'_*_matrix.mtx.gz'))
    ann=pd.read_csv(features,sep='\t',header=None,names=['ENSEMBL','SYMBOL','type'])
    with gzip.open(matrix,'rb') as f:X=mmread(f,spmatrix=True).tocsc()
    lib=np.asarray(X.sum(axis=0)).ravel();nf=np.diff(X.indptr);mt=np.asarray(X[ann.SYMBOL.str.lower().str.startswith('mt-').to_numpy()].sum(axis=0)).ravel()/np.maximum(lib,1)
    keep=(nf>=200)&(nf<=2500)&(mt<=.5);depth=keep&(lib>=500);assert keep.sum()==sample.n_retained_cells
    quality.append({'sample':sample.sample,'donor':sample.donor,'tissue':sample.tissue,'flight':sample.flight,'age':sample.age,'n_primary_cells':int(keep.sum()),'n_500_UMI_cells':int(depth.sum()),'median_primary_UMIs':float(np.median(lib[keep])),'median_primary_features':float(np.median(nf[keep])),'median_primary_mito_fraction':float(np.median(mt[keep]))})
    lookup=ann.reset_index().groupby('SYMBOL')['index'].agg(list).to_dict()
    for p in panel.itertuples():
        hit=lookup.get(p.mouse_symbol_primary,[]);ok=len(hit)==1
        metrics={'raw_detection_fraction':np.nan,'mean_cell_log2_CP10k_plus1':np.nan,'expected_detection_at_500_UMIs':np.nan}
        if ok:
            k=X[hit[0],:].toarray().ravel();metrics['raw_detection_fraction']=float((k[keep]>0).mean());metrics['mean_cell_log2_CP10k_plus1']=float(np.log2(1+1e4*k[keep]/lib[keep]).mean())
            n=lib[depth];z=k[depth];prob=np.zeros(len(n));pos=z>0;certain=pos&(n-z<500);calc=pos&~certain;prob[certain]=1
            logp=gammaln(n[calc]-z[calc]+1)-gammaln(n[calc]-z[calc]-500+1)+gammaln(n[calc]-500+1)-gammaln(n[calc]+1)
            prob[calc]=-np.expm1(np.minimum(logp,0));assert np.all((prob>=0)&(prob<=1))
            metrics['expected_detection_at_500_UMIs']=float(prob.mean()) if len(prob)>=50 else np.nan
        for metric,value in metrics.items():rows.append({'sample':sample.sample,'donor':sample.donor,'tissue':sample.tissue,'flight':sample.flight,'age':sample.age,'analyte':p.analyte,'mouse_symbol':p.mouse_symbol_primary,'endpoint':metric,'value':value,'n_cells':int(depth.sum() if metric.startswith('expected') else keep.sum())})
    print(sample.sample,sample.tissue,int(keep.sum()),int(depth.sum()),flush=True)
pd.DataFrame(rows).to_csv(B/'tables/sorted_cell_depth_donor_endpoints.csv',index=False);pd.DataFrame(quality).to_csv(B/'tables/sorted_cell_depth_QC.csv',index=False)
(B/'tables/sorted_cell_depth_definition.json').write_text(json.dumps({'rarefaction_UMIs':500,'method':'Analytical expectation: 1 - choose(N-k,500)/choose(N,500); primary-QC cells with N>=500','exposure_unit':'Mouse, never individual cell','endpoints':['raw_detection_fraction','mean_cell_log2_CP10k_plus1','expected_detection_at_500_UMIs'],'families':'All 151 targets x 3 myeloid populations x 3 endpoints = 1359 per donor error specification','scope':'Depth sensitivity on the same biological samples; does not remove cell-state mixtures, ambient RNA, dissociation or day confounding.'},indent=2))
