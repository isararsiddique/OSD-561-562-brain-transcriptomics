from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Current GEO sorted-brain-cell records, independent cell QC and donor pseudobulk."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import re,json,gzip,io,hashlib,requests
import numpy as np,pandas as pd
from scipy.io import mmread
B=stage_root('extended')
text=(B/'sources/GSE295428_current_sample_SOFT.txt').read_text()
samples=[]
for block in re.split(r'(?=\^SAMPLE = )',text):
    m=re.match(r'\^SAMPLE = (GSM\d+)',block)
    if not m:continue
    def val(key):
        s=re.search(r'^!'+re.escape(key)+r' = (.*)$',block,re.M);return s[1].strip() if s else ''
    tissue=val('Sample_source_name_ch1')
    if not any(t in tissue.lower() for t in ['myeloid_cerebellum','myeloid_forebrain','myeloid_hippocampus']):continue
    title=val('Sample_title');gsm=m[1];mouse=re.search(r'Mouse (\d+)',title)[1]
    fields=dict(re.findall(r'^!Sample_characteristics_ch1 = ([^:]+): (.+)\s*$',block,re.M));fields={k:v.strip() for k,v in fields.items()}
    treatment=fields['treatment'];group='FLT' if treatment=='spaceflight' else 'GC';assert treatment in ['spaceflight','control']
    urls=re.findall(r'^!Sample_supplementary_file_\d+ = (.+)\s*$',block,re.M);urls=[u.strip().replace('ftp://','https://') for u in urls]
    assert len(urls)==3
    samples.append({'sample':gsm,'donor':group+'_Mouse'+mouse,'tissue':tissue,'title':title,'flight':group,'age':int(fields['age']),'batch':fields.get('batch',''),'sex':fields.get('Sex',''),'strain':fields.get('strain',''),'files':urls})
assert len(samples)==47,(len(samples),[(s['sample'],s['tissue']) for s in samples])
folder=B/'sources/GSE295428_sorted_cells';folder.mkdir(exist_ok=True)
def download(u):
    p=folder/u.rsplit('/',1)[-1]
    if p.exists():
        try:
            with gzip.open(p,'rb') as f:
                while f.read(4*1024*1024):pass
            return {'URL':u,'file':record_path(p),'bytes':p.stat().st_size,'SHA256':hashlib.sha256(p.read_bytes()).hexdigest(),'gzip_CRC_checked':True}
        except Exception:pass
    for attempt in range(3):
        try:
            with request_get(u,stream=True,timeout=(25,90)) as r:
                r.raise_for_status()
                with p.open('wb') as f:
                    for c in r.iter_content(1024*1024):f.write(c)
                if r.headers.get('Content-Length'):assert p.stat().st_size==int(r.headers['Content-Length'])
            with gzip.open(p,'rb') as f:
                while f.read(4*1024*1024):pass
            return {'URL':u,'file':record_path(p),'bytes':p.stat().st_size,'SHA256':hashlib.sha256(p.read_bytes()).hexdigest(),'gzip_CRC_checked':True}
        except Exception:
            if attempt==2:raise
with ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(download,[u for s in samples for u in s['files']]))
(B/'sources/sorted_cell_download_provenance.json').write_text(json.dumps(records,indent=2))
pd.DataFrame([{k:v for k,v in s.items() if k!='files'} for s in samples]).to_csv(B/'tables/GSE295428_current_sample_registry.csv',index=False)
qc=[];pb={};relaxed={};annotation=None
for s in samples:
    paths={u.rsplit('_',1)[-1]:folder/u.rsplit('/',1)[-1] for u in s['files']}
    feat=pd.read_csv(paths['features.tsv.gz'],sep='\t',header=None,names=['ENSEMBL','SYMBOL','type'])
    bar=pd.read_csv(paths['barcodes.tsv.gz'],sep='\t',header=None)[0]
    with gzip.open(paths['matrix.mtx.gz'],'rb') as f:X=mmread(f).tocsc()
    assert X.shape==(len(feat),len(bar)) and feat.ENSEMBL.is_unique and bar.is_unique and (X.data>=0).all()
    if annotation is None:annotation=feat
    else:assert np.array_equal(feat.ENSEMBL,annotation.ENSEMBL)
    nfeature=np.diff(X.indptr);library=np.asarray(X.sum(axis=0)).ravel();mt=np.asarray(X[feat.SYMBOL.str.lower().str.startswith('mt-').to_numpy()].sum(axis=0)).ravel()/np.maximum(library,1)
    keep=(nfeature>=200)&(nfeature<=2500)&(mt<=.5)
    broad=(nfeature>=200)&(nfeature<=6000)&(mt<=.2)
    pb[s['sample']]=np.asarray(X[:,keep].sum(axis=1)).ravel();relaxed[s['sample']]=np.asarray(X[:,broad].sum(axis=1)).ravel()
    qc.append(dict(s,files=None,n_input_barcodes=len(bar),n_retained_cells=int(keep.sum()),n_sensitivity_cells=int(broad.sum()),raw_UMIs=int(library.sum()),retained_UMIs=int(library[keep].sum()),median_retained_mt=float(np.median(mt[keep])) if keep.any() else None))
    print(s['sample'],s['tissue'],s['donor'],len(bar),int(keep.sum()),flush=True)
pd.DataFrame(pb,index=annotation.ENSEMBL).to_csv(B/'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz')
pd.DataFrame(relaxed,index=annotation.ENSEMBL).to_csv(B/'prepared/GSE295428_sorted_sensitivity_pseudobulk.csv.gz')
annotation.to_csv(B/'prepared/GSE295428_sorted_annotation.csv',index=False)
qc=pd.DataFrame(qc).drop(columns='files');qc.to_csv(B/'tables/GSE295428_sorted_cell_QC.csv',index=False)
assert qc.donor.nunique()==8
(B/'tables/sorted_cell_preparation_summary.json').write_text(json.dumps({'GEO':'GSE295428','sample_libraries':len(qc),'unique_donor_labels':8,'brain_tissue_populations':qc.tissue.value_counts().to_dict(),'input_barcodes':int(qc.n_input_barcodes.sum()),'retained_cells':int(qc.n_retained_cells.sum()),'main_QC':'GEO-record feature/mitochondrial gate: 200-2500 features, <=50% mitochondrial counts. Samples require >=50 retained cells for exposure fitting.','sensitivity_QC':'200-6000 features, <=20% mitochondrial counts. Predeclared alternative, not selected by result.','processing_scope':'Independent QC and pseudobulk from current attached Cell Ranger gene matrices. No claim of exact author SoupX/DoubletFinder/Seurat reproduction; postprocessing status of attachments not established. FACS labels represent enriched populations, not a proven pure microglial subtype.','sample_identity':'Current GEO title, characteristics and sample supplementary filenames; Use verified deposited sample identities. Shared mouse labels across tissues remain eight exposure donors.'},indent=2))
