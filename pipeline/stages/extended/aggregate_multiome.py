from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""All RNA and ATAC donor pseudobulk, bounded sparse H5 blocks."""
from pathlib import Path
from collections import Counter
import hashlib,json
import numpy as np
import pandas as pd
import h5py

B=stage_root('extended'); ROOT=inputs_root()
SRC=ROOT/'research/brain_integrative_2026/data/OSD-352'
meta=pd.read_csv(SRC/'multiomics_brain_metadata_all_cells.csv',index_col=0)
assert meta.index.is_unique
groups=sorted(set(zip(meta.bio_origin,meta['predicted.celltype.l2'])))
fine=sorted(set(zip(meta.bio_origin,meta['predicted.celltype.l3'])))
gi={g:i for i,g in enumerate(groups)}; fi={g:i for i,g in enumerate(fine)}
donor=meta.groupby('bio_origin').cell_flight_status.agg(lambda a:a.unique()[0])
assert meta.groupby('bio_origin').cell_flight_status.nunique().eq(1).all()
codes={s:gi[(r.bio_origin,r['predicted.celltype.l2'])] for s,r in meta.iterrows()}
fcodes={s:fi[(r.bio_origin,r['predicted.celltype.l3'])] for s,r in meta.iterrows()}
path=SRC/'GLDS-352_snRNA-Seq_filtered_feature_bc_matrix.h5'
with h5py.File(path,'r') as h:
    m=h['matrix']; names=m['features/name'][:].astype(str); ids=m['features/id'][:].astype(str)
    types=m['features/feature_type'][:].astype(str); bar=m['barcodes'][:].astype(str); ptr=m['indptr'][:]
    assert set(meta.index).issubset(bar)
    rnaidx=np.flatnonzero(types=='Gene Expression'); peakidx=np.flatnonzero(types=='Peaks')
    rna_map=np.full(len(names),-1,int);rna_map[rnaidx]=np.arange(len(rnaidx))
    broad=np.zeros((len(names),len(groups)),np.int64)
    finer=np.zeros((len(rnaidx),len(fine)),np.int64)
    detected=np.zeros((len(rnaidx),len(groups)),np.int64)
    for start in range(0,len(bar),512):
        end=min(start+512,len(bar));lo,hi=ptr[start],ptr[end]
        ix=m['indices'][lo:hi];v=m['data'][lo:hi]
        col=np.repeat(np.arange(end-start),np.diff(ptr[start:end+1]))
        g=np.array([codes.get(b,-1) for b in bar[start:end]])[col]
        fg=np.array([fcodes.get(b,-1) for b in bar[start:end]])[col]
        ok=g>=0;flat=ix[ok]*len(groups)+g[ok]
        broad+=np.bincount(flat,weights=v[ok],minlength=broad.size).reshape(broad.shape).astype(np.int64)
        ri=rna_map[ix];ok=(fg>=0)&(ri>=0)
        finer+=np.bincount(ri[ok]*len(fine)+fg[ok],weights=v[ok],minlength=finer.size).reshape(finer.shape).astype(np.int64)
        ok=(g>=0)&(ri>=0)
        detected+=np.bincount(ri[ok]*len(groups)+g[ok],weights=(v[ok]>0),minlength=detected.size).reshape(detected.shape).astype(np.int64)
cols=[f'{d}|{ct}' for d,ct in groups];finecols=[f'{d}|{ct}' for d,ct in fine]
pd.DataFrame(broad[rnaidx],index=ids[rnaidx],columns=cols).to_csv(B/'prepared/OSD352_RNA_pseudobulk.csv.gz')
pd.DataFrame(broad[peakidx],index=names[peakidx],columns=cols).to_csv(B/'prepared/OSD352_ATAC_pseudobulk.csv.gz')
pd.DataFrame(finer,index=ids[rnaidx],columns=finecols).to_csv(B/'prepared/OSD352_RNA_subtype_pseudobulk.csv.gz')
pd.DataFrame(detected,index=ids[rnaidx],columns=cols).to_csv(B/'prepared/OSD352_RNA_detected_nuclei.csv.gz')
pd.DataFrame({'ENSEMBL':ids[rnaidx],'SYMBOL':names[rnaidx]}).to_csv(B/'prepared/OSD352_RNA_annotation.csv',index=False)
pd.DataFrame({'feature_id':ids,'symbol_or_peak':names,'feature_type':types}).to_csv(B/'tables/OSD352_feature_registry.csv.gz',index=False)
for level,gr,c in [('l2',groups,cols),('l3',fine,finecols)]:
    n=Counter(zip(meta.bio_origin,meta[f'predicted.celltype.{level}']))
    rows=[]
    for i,((d,ct),column) in enumerate(zip(gr,c)):
        rows.append({'column':column,'bio_origin':d,'cell_type':ct,'flight_status':donor[d],
                     'n_nuclei':n[(d,ct)],'RNA_total_UMIs':int((broad[rnaidx,i] if level=='l2' else finer[:,i]).sum()),
                     'ATAC_total_counts':int(broad[peakidx,i].sum()) if level=='l2' else None})
    pd.DataFrame(rows).to_csv(B/f'tables/OSD352_{level}_registry.csv',index=False)
summary={'H5_SHA256':hashlib.sha256(path.read_bytes()).hexdigest(),'matched_nuclei':len(meta),'all_barcodes':len(bar),
         'RNA_features':len(rnaidx),'ATAC_peaks':len(peakidx),'broad_donor_celltype_groups':len(groups),'subtype_donor_groups':len(fine),
         'RNA_total':int(broad[rnaidx].sum()),'ATAC_total':int(broad[peakidx].sum()),
         'provenance':'Author bio_origin and predicted.celltype labels. No invented publication animal aliases. All assays share five donor groups.'}
(B/'tables/multiome_aggregation_summary.json').write_text(json.dumps(summary,indent=2))
print(summary,flush=True)
