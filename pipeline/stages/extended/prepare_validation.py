from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Fresh ISA identity and deposited checksums for independent mission RNA."""
from pathlib import Path
import pandas as pd,numpy as np,hashlib,json
B=stage_root('extended');ROOT=inputs_root()
checks=[];allmeta=[]
for sid,glds in [(352,352),(563,568),(564,569)]:
    src=B/'sources'/f'OSD-{sid}'
    p=src/f'GLDS-{glds}_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv'
    side=src/f'GLDS-{glds}_rna_seq_processed_md5sum_GLbulkRNAseq.tsv'
    records=dict(line.split('\t')[:2] for line in side.read_text().splitlines() if '\t' in line)
    expected=records[p.name.split('_rna_seq_',1)[1]];actual=hashlib.md5(p.read_bytes()).hexdigest();assert actual==expected
    checks.append({'study':sid,'file':record_path(p),'deposited_MD5':expected,'actual_MD5':actual,'matches':expected==actual})
    count=pd.read_csv(p,index_col=0);count=count.loc[count.index.str.match(r'^ENSMUSG\d+(?:\.\d+)?$')]
    assert count.index.is_unique and count.columns.is_unique and (count>=0).all().all()
    isa=pd.read_csv(B/f'sources/OSD-{sid}_sample_ISA.csv').set_index('Sample Name')
    md=isa.loc[count.columns].copy();assert md['Source Name'].is_unique
    md['animal']=md['Source Name'];md['study']=f'OSD-{sid}';md['flight']=md['Factor Value[Spaceflight]'].map({'Space Flight':'FLT','Ground Control':'GC'});assert md.flight.notna().all()
    md['group']=md.flight;md.index.name='sample'
    md.to_csv(B/f'tables/OSD-{sid}_metadata.csv');allmeta.append(md.reset_index())
    count.to_csv(B/f'prepared/OSD-{sid}_STAR_endogenous.csv.gz')
    # Ensembl/SYMBOL annotation from the deposited NASA differential expression table.
    ann=pd.read_csv(ROOT/f'research/brain_integrative_2026/data/OSD-{sid}/GLDS-{glds}_rna_seq_differential_expression_GLbulkRNAseq.csv',usecols=['ENSEMBL','SYMBOL','GENENAME']).drop_duplicates('ENSEMBL')
    ann.to_csv(B/f'prepared/OSD-{sid}_annotation.csv',index=False)
    print(sid,count.shape,md.flight.value_counts().to_dict(),flush=True)
pd.DataFrame(checks).to_csv(B/'tables/validation_checksum_verification.csv',index=False)
pd.concat(allmeta,ignore_index=True).to_csv(B/'tables/validation_sample_registry.csv',index=False)
