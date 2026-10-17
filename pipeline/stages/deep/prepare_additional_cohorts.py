from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Prospectively selected additional cohorts; primary ISA sample/animal links."""
from pathlib import Path
import gzip, re, json, hashlib
import numpy as np
import pandas as pd
B=stage_root('deep')
gtf=stage_root('extended') / 'sources/Mus_musculus.GRCm38.102.gtf.gz'
genes=[]
with gzip.open(gtf,'rt') as f:
    for line in f:
        if line.startswith('#'):continue
        fields=line.rstrip().split('\t')
        if fields[2]!='gene':continue
        a=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
        genes.append({'ENSEMBL':a['gene_id'],'SYMBOL':a.get('gene_name',''),'gene_biotype':a.get('gene_biotype','')})
ann=pd.DataFrame(genes);assert ann.ENSEMBL.is_unique
ann.to_csv(B/'prepared/additional_gene_annotation.csv',index=False)
reg=[];audits=[]
for sid in [525,457]:
    fname=B/'sources'/f'GLDS-{sid}_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv'
    x=pd.read_csv(fname,index_col=0);assert x.index.is_unique and x.columns.is_unique
    assert np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all() and np.equal(x.to_numpy(),np.floor(x.to_numpy())).all()
    isa=pd.read_csv(B/'sources'/f's_OSD-{sid}.txt',sep='\t').set_index('Sample Name',drop=False)
    runs=pd.read_csv(B/'sources'/f'GLDS-{sid}_rna_seq_bulkRNASeq_v2_runsheet.csv').set_index('Sample Name',drop=False)
    if sid==457:
        isa=isa[isa['Factor Value[Tissue]'].eq('Cerebrum')]
    assert isa.index.is_unique and isa.index.isin(x.columns).all() and isa.index.isin(runs.index).all()
    selected=x.loc[:,isa.index]
    mouse=selected.index.str.match('^ENSMUSG')
    selected=selected.loc[mouse]
    pb={}
    for donor,s in isa.groupby('Source Name',sort=True):
        assert s['Factor Value[Spaceflight]'].nunique()==1
        exposure=s['Factor Value[Spaceflight]'].iloc[0]
        flight={'Space Flight':'FLT','Ground Control':'GC','Ground control':'GC','Vivarium Control':'VIV'}[exposure]
        genotype='WT' if sid==525 or s['Factor Value[Genotype]'].iloc[0]=='Wild Type' else 'Nrf2KO'
        key=f'OSD{sid}_{donor}'
        assert (runs.loc[s.index,'Factor Value[Spaceflight]'].str.lower()==exposure.lower()).all()
        pb[key]=selected.loc[:,s.index].sum(axis=1)
        reg.append({'OSD':f'OSD-{sid}','key':key,'donor':donor,'flight':flight,'genotype':genotype,
                    'region':'Cerebellum' if sid==525 else 'Cerebrum','sample_libraries':';'.join(s.index),
                    'n_technical_libraries':len(s),'strain':s['Characteristics[Strain]'].iloc[0],
                    'sex':s['Characteristics[Sex]'].iloc[0],
                    'preservation':s['Parameter Value[Sample Preservation Method]'].iloc[0],
                    'euthanasia_date':s['Comment[Euthanasia Date]'].iloc[0] if sid==525 else '',
                    'dissection_date':s['Comment[BSP Dissection Date]'].iloc[0] if sid==525 else '',
                    'raw_mouse_counts':int(pb[key].sum())})
    out=pd.DataFrame(pb);assert out.sum().sum()==selected.sum().sum()
    out.to_csv(B/'prepared'/f'OSD-{sid}_mouse_counts.csv.gz')
    audits.append({'OSD':f'OSD-{sid}','input_rows':len(x),'input_libraries':x.shape[1],
                   'selected_tissue_libraries':selected.shape[1],'independent_animals':out.shape[1],
                   'retained_mouse_gene_rows':len(out),'excluded_nonmouse_ERCC_or_other_rows':int((~mouse).sum()),
                   'integer_counts_verified':True,'technical_replicates_summed':True,
                   'all_retained_counts_conserved':True,'unique_primary_ISA_sample_links':True,
                   'source_SHA256':hashlib.sha256(fname.read_bytes()).hexdigest()})
pd.DataFrame(reg).to_csv(B/'tables/additional_mouse_registry.csv',index=False)
(B/'tables/additional_count_preparation_audit.json').write_text(json.dumps(audits,indent=2))
(B/'sources/additional_annotation_provenance.json').write_text(json.dumps({'path':record_path(gtf),'SHA256':hashlib.sha256(gtf.read_bytes()).hexdigest(),'use':'Existing official Ensembl GRCm38.102 annotation, keyed by stable gene ID; nonunique symbols remain unavailable for panel inference'},indent=2))
print(pd.DataFrame(reg).groupby(['OSD','flight','genotype']).size().to_string(),flush=True)
print(json.dumps(audits,indent=2),flush=True)
