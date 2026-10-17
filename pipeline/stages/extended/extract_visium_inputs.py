from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Retrieve only counts and spot coordinates with ZIP CRC and byte-range validation."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,requests,struct,zlib,hashlib
B=stage_root('extended');inv=json.loads((B/'sources/visium_archive_inventory.json').read_text());url=inv['url'];size=inv['bytes']
def chunk(start,end):
    for attempt in range(3):
        try:
            r=request_get(url,headers={'Range':f'bytes={start}-{end}'},timeout=90)
            assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {start}-{end}/{size}'
            assert len(r.content)==end-start+1
            return r.content
        except Exception:
            if attempt==2:raise
def get(m):
    name=m['file'];p=B/'sources'/name
    if p.exists() and p.stat().st_size==m['size'] and zlib.crc32(p.read_bytes())&0xffffffff==m['CRC']:return {'member':name,'file':record_path(p),'CRC_verified':True,'SHA256':hashlib.sha256(p.read_bytes()).hexdigest()}
    head=chunk(m['offset'],m['offset']+29);a=struct.unpack('<IHHHHHIIIHH',head);assert a[0]==0x04034b50
    start=m['offset']+30+a[-2]+a[-1];end=start+m['compressed_size']
    data=b''.join(chunk(k,min(k+4*1024*1024,end)-1) for k in range(start,end,4*1024*1024))
    raw=zlib.decompress(data,-15) if m['compression']==8 else data
    assert len(raw)==m['size'] and zlib.crc32(raw)&0xffffffff==m['CRC']
    p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
    return {'member':name,'file':record_path(p),'CRC_verified':True,'SHA256':hashlib.sha256(raw).hexdigest()}
selected=[m for m in inv['members'] if m['file'].startswith('visium_data/') and m['file'].endswith(('filtered_feature_bc_matrix.h5','tissue_positions_list.csv','scalefactors_json.json'))]
with ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(get,selected))
(B/'sources/visium_extraction_provenance.json').write_text(json.dumps({'archive_url':url,'records':records,'extraction':'Individual validated ZIP members, no image download'},indent=2));print('Verified Visium members',len(records),flush=True)
