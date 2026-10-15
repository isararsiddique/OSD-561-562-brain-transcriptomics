"""Independent handoff checks and an extension file manifest."""
from pathlib import Path
import json,hashlib,csv,zipfile,re
import pandas as pd
B=Path(__file__).resolve().parents[1];D=B/'inputs/RR10';T=B/'tables'
plan_sha=hashlib.sha256((B/'RETROSPECTIVE_SEARCH_PLAN.md').read_bytes()).hexdigest();assert plan_sha=='a103dce4736d4d5399262945c63bb25b9777dfa57adbefcf7692423345559d90'
prov=pd.read_csv(T/'RR10_download_provenance.csv');sj=prov[prov.file.str.endswith('SJ.out.tab')];assert len(sj)==10;assert sj.publisher_md5_match.astype(str).str.lower().eq('true').all()
m=pd.read_csv(D/'RR10_animal_junction_counts.csv.gz',index_col=0);a=pd.read_csv(D/'RR10_animal_metadata.csv');assert m.columns.tolist()==a['sample'].tolist();assert a.group.value_counts().to_dict()=={'FLT':5,'GC':5};assert a['Comment[RFID]'].nunique()==10
# Verify the independently constructed checks persisted against current matrix.
c=pd.read_csv(T/'RR10_matrix_checks.csv').set_index('sample');assert m.sum().to_dict()==c.retained_read_sum.to_dict()
summary=json.loads((T/'RR10_freeze_summary.json').read_text());assert hashlib.sha256((D/'RR10_animal_junction_counts.csv.gz').read_bytes()).hexdigest()==summary['matrix_sha256'];assert hashlib.sha256((D/'RR10_animal_metadata.csv').read_bytes()).hexdigest()==summary['metadata_sha256']
md5=pd.read_csv(D/'GLDS-568_rna_seq_processed_md5sum_GLbulkRNAseq.tsv',sep='\t',header=None,names=['file','md5']).set_index('file').md5.to_dict();zip_path=D/'GLDS-568_rna_seq_processing_info_GLbulkRNAseq.zip';assert hashlib.md5(zip_path.read_bytes()).hexdigest()==md5['processing_info_GLbulkRNAseq.zip']
z=zipfile.ZipFile(zip_path);s=z.read('nextflow_processing_info_GLbulkRNAseq.txt').decode();star=[x for x in s.splitlines() if x.strip().startswith('STAR --') and '--readFilesIn' in x];assert len(star)==10;assert all('--twopassMode Basic' in x and 'GRCm39.dna.primary_assembly_and_ERCC92' in x for x in star)
(D/'RR10_STAR_command_evidence.txt').write_text('\n'.join(star)+'\n')
meta=json.loads((B/'sources/OSD-563_meta.json').read_text())['study']['OSD-563']['studies'][0];protocol=next(x['description'] for x in meta['protocols'] if x['name']=='GeneLab RNAseq data processing protocol');(D/'RR10_publisher_processing_protocol.txt').write_text(protocol+'\n');assert '2.7.11b' in protocol and 'Ensembl release 112' in protocol
ex=json.loads((T/'RR10_publisher_checksum_exception.json').read_text());assert ex['actual_md5']!=ex['publisher_sidecar_md5'];assert hashlib.md5((D/ex['file']).read_bytes()).hexdigest()==ex['actual_md5']
cat=pd.read_csv(T/'complete_accession_decisions.csv');assert len(cat)==73 and cat.accession.nunique()==73 and cat.decision.notna().all();lit=pd.read_csv(T/'complete_literature_title_abstract_screen.csv');assert len(lit)==60 and lit.id.nunique()==60
results=dict(status='PASS_WITH_DOCUMENTED_PUBLISHER_METADATA_CHECKSUM_EXCEPTION',frozen_plan_hash_unchanged=True,star_files_md5_matched=10,count_matrix_shape=list(m.shape),sample_order_identity_and_5_plus_5_design_verified=True,stored_independent_reconstruction_totals_match=True,matrix_and_metadata_hashes_match=True,processing_zip_md5_matches=True,ten_STAR_commands_two_pass_and_GRCm39_verified=True,protocol_STAR2711b_Ensembl112_verified=True,software_versions_checksum_exception_retained=True,complete_accession_decisions=73,complete_date_bounded_title_abstract_decisions=60)
(T/'handoff_verification.json').write_text(json.dumps(results,indent=2)+'\n')
rows=[]
for p in sorted(B.rglob('*')):
 if p.is_file() and p.name!='FILE_MANIFEST.csv' and '__pycache__' not in str(p):rows.append(dict(path=str(p.relative_to(B)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(rows).to_csv(T/'FILE_MANIFEST.csv',index=False);print(json.dumps(results,indent=2));print('manifested',len(rows),'files')
