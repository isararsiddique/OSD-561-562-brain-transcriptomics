from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Recheck retained local assay inputs against retained primary-server downloads."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,gzip,requests,pandas as pd
OUT=stage_root('core');ROOT=inputs_root();SRC=OUT/'sources';TAB=OUT/'tables'
jobs=[]
for n in [561,562]:
    obj=json.loads((SRC/f'OSD-{n}_files.json').read_text())
    files=obj['studies'][f'OSD-{n}']['study_files']
    f=next(f for f in files if f['file_name'].endswith('ISA.zip'))
    url='https://osdr.nasa.gov'+f['remote_url']
    local=ROOT/'data'/f'OSD-{n}'/f'OSD-{n}_metadata_OSD-{n}-ISA.zip'
    jobs.append((f'OSD-{n}_ISA',url,local,False))
for code in ['CA','DG','FCT','CT']:
    name=f'GSE239336_{code}_GeneExpression_Q3norm.txt.gz'
    url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239336/suppl/'+name
    jobs.append((name,url,ROOT/'OSD-682-685-698-699-brain-spatial/data/geo'/name,True))
def verify(job):
    label,url,local,decompress=job
    cached=SRC/(label+'.download')
    if network_allowed():
        resp=request_get(url,timeout=45);resp.raise_for_status();raw=resp.content
        cached.write_bytes(raw)
    else:
        raw=cached.read_bytes()
    old=local.read_bytes()
    a=gzip.decompress(raw) if decompress else raw;b=gzip.decompress(old) if decompress else old
    result={'input':label,'url':url,'local_path':str(local),'fresh_bytes':len(raw),'raw_local_sha256':hashlib.sha256(old).hexdigest(),
            'raw_fresh_sha256':hashlib.sha256(raw).hexdigest(),'content_matches_fresh_primary':a==b,'decompressed_comparison':decompress}
    assert a==b,label+' differs from primary deposit';return result
with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(verify,jobs))
pd.DataFrame(rows).to_csv(TAB/'fresh_primary_input_verification.csv',index=False)
# Verify the retained cell annotation against its originally verified ZIP CRC.
record=json.loads((ROOT/'research/brain_integrative_2026/tables/rr3_metadata_download.json').read_text())
local=ROOT/'research/brain_integrative_2026/data/OSD-352/multiomics_brain_metadata_all_cells.csv'
raw=local.read_bytes();import zlib
assert len(raw)==record['member_bytes'] and zlib.crc32(raw)&0xffffffff==record['member_crc32']
assert hashlib.sha256(raw).hexdigest()==record['sha256']
(SRC/'OSD-352_annotation_zip_provenance.json').write_text(json.dumps(record,indent=2))
print('Verified six snapshot assay inputs and retained OSD-352 annotation ZIP CRC/SHA identity.')
