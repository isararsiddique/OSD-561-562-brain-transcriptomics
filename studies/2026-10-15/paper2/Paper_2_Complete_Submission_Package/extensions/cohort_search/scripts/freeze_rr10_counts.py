"""Freeze outcome-independent RR10 full junction matrix and official metadata."""
from pathlib import Path
import json,csv,hashlib,zipfile
import pandas as pd,numpy as np
B=Path(__file__).resolve().parents[1];D=B/'inputs'/'RR10';T=B/'tables'
md5=pd.read_csv(D/'GLDS-568_rna_seq_processed_md5sum_GLbulkRNAseq.tsv',sep='\t',header=None,names=['file','md5']).set_index('file').md5.to_dict()
prov=pd.read_csv(T/'RR10_download_provenance.csv').fillna('')
cols=[];checks=[]
for p in sorted(D.glob('*SJ.out.tab')):
 raw=p.read_bytes();short=p.name.removeprefix('GLDS-568_rna_seq_'); expected=md5[short];actual=hashlib.md5(raw).hexdigest();assert expected==actual
 d=pd.read_csv(p,sep='\t',header=None,names=['chr','start','end','strand','motif','annotated','unique','multi','overhang'],dtype={'chr':str})
 retain=d.strand.isin([1,2]) & d.motif.between(1,6) & d.annotated.eq(1) & d.overhang.ge(8) & d.unique.gt(0) & d.chr.isin([str(x) for x in range(1,20)]+['X','Y','MT'])
 keep=d.loc[retain];key=keep.chr+':'+keep.strand.astype(str)+':'+keep.start.astype(str)+':'+keep.end.astype(str)
 sample=short.removesuffix('_SJ.out.tab');ser=pd.Series(keep.unique.to_numpy(),index=key,name=sample).groupby(level=0).sum();cols.append(ser)
 checks.append(dict(sample=sample,publisher_md5=expected,actual_md5=actual,publisher_md5_match=True,total_deposited_rows=len(d),retained_rows=len(keep),unique_retained_junctions=len(ser),retained_read_sum=int(ser.sum()),unique_source_file=True))
 for i in prov.index[prov.file.eq(p.name)]:prov.loc[i,['publisher_md5','publisher_md5_match']]=[expected,True]
m=pd.concat(cols,axis=1).fillna(0).astype('int64').sort_index();m.to_csv(D/'RR10_animal_junction_counts.csv.gz',compression={'method':'gzip','mtime':0})
z=zipfile.ZipFile(D/'OSD-563_metadata_OSD-563-ISA.zip');meta=pd.read_csv(z.open('s_OSD-563.txt'),sep='\t').rename(columns={'Sample Name':'sample'});meta=meta.set_index('sample').loc[m.columns].rename_axis('sample').reset_index();meta['group']=meta['Factor Value[Spaceflight]'].map({'Space Flight':'FLT','Ground Control':'GC'});assert meta.group.notna().all();assert meta['Comment[RFID]'].nunique()==10;assert meta['Source Name'].nunique()==10
meta.to_csv(D/'RR10_animal_metadata.csv',index=False)
# Archive the source assay table and exact processing reference text.
a=next(n for n in z.namelist() if 'rna-sequencing-' in n);(D/'RR10_original_assay.tsv').write_bytes(z.read(a))
rr23=pd.read_csv(B/'../../inputs/animal_metadata.csv');assert not set(meta['Source Name'])&set(rr23['Source Name']);assert not set(meta['Comment[RFID]'].astype(str))&set(rr23['Comment[RFID]'].astype(str))
for s in checks:np.testing.assert_array_equal(m[s['sample']].sum(),s['retained_read_sum'])
# Independent per-line reconstruction of all retained columns (no pandas filter reuse).
for p in sorted(D.glob('*SJ.out.tab')):
 recon={}
 for line in p.open():
  c,st,en,strand,motif,annot,unique,multi,over=line.rstrip().split('\t');
  if strand in ('1','2') and 1<=int(motif)<=6 and annot=='1' and int(over)>=8 and int(unique)>0 and c in [str(x) for x in range(1,20)]+['X','Y','MT']:
   k=f'{c}:{strand}:{st}:{en}';recon[k]=recon.get(k,0)+int(unique)
 sample=p.name.removeprefix('GLDS-568_rna_seq_').removesuffix('_SJ.out.tab')
 observed=m[sample];assert {k:int(v) for k,v in observed.items() if v}==recon
pd.DataFrame(checks).to_csv(T/'RR10_matrix_checks.csv',index=False);prov.to_csv(T/'RR10_download_provenance.csv',index=False)
summary=dict(accession='OSD-563',label='RR10 distinct-animal same-launch transport cohort',mission='RR-10',launch='SpaceX-21 shared with RR23',shape=list(m.shape),n_flt=5,n_gc=5,publisher_md5_verified_star_files=10,all_columns_independently_reconstructed=True,source_name_and_RFID_overlap_with_RR23=0,sex='Female',strain='B6129SF2/J',genotype='Wild Type',age_at_launch_weeks='14-15',shared_marker_injections=['Alizarin Red','Calcein','BrdU'],causal_limit='Exposure perfectly aliased with euthanasia date (FLT 2021-01-02; GC 2021-01-05) and dissection date (FLT 2021-05-12; GC 2021-05-14).',gene_family_outcomes_examined=False,filters='STAR unique reads; annotated=1; motif1..6; strand1/2; overhang>=8; chr1..19,X,Y,MT',matrix_sha256=hashlib.sha256((D/'RR10_animal_junction_counts.csv.gz').read_bytes()).hexdigest(),metadata_sha256=hashlib.sha256((D/'RR10_animal_metadata.csv').read_bytes()).hexdigest(),statistical_resolution='5+5 animals: 252 unrestricted assignments, one-sided min p 1/252; two-sided absolute statistic min p at least 2/252; no 923-test Holm discovery attainable.')
(T/'RR10_freeze_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
