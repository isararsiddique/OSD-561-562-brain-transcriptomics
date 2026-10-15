"""Acquire complete RR10 cerebellum STAR junction files, irrespective of outcomes."""
from pathlib import Path
import json,csv,urllib.request,hashlib,datetime,concurrent.futures
B=Path(__file__).resolve().parents[1];D=B/'inputs'/'RR10';D.mkdir(parents=True,exist_ok=True)
inv=json.loads((B/'sources/OSD-563_files.json').read_text())['studies']['OSD-563']['study_files']
files=[f for f in inv if any(x in f['file_name'] for x in ['SJ.out.tab','processed_md5sum','runsheet','software_versions','ISA.zip','processing_metadata','processing_info'])]
def fetch(f):
 p=D/f['file_name'];url='https://osdr.nasa.gov'+f['remote_url'];row=dict(file=f['file_name'],url=url,expected_bytes=f['file_size'],status='',retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  b=p.read_bytes() if p.exists() else urllib.request.urlopen(url,timeout=60).read();p.write_bytes(b);row.update(status='downloaded',bytes=len(b),md5=hashlib.md5(b).hexdigest(),sha256=hashlib.sha256(b).hexdigest())
 except Exception as e:row.update(status='error',error=str(e))
 print(row['file'],row['status'],flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:rows=list(e.map(fetch,files))
md5p=next(D.glob('*processed_md5sum*'));md5={}
for l in md5p.read_text().splitlines():
 t=l.split()
 if len(t)>=2:md5[t[0].lstrip('*')]=t[1]
for r in rows:
 fname=r['file'];exp=md5.get(fname) or md5.get(fname.removeprefix('GLDS-568_rna_seq_'));r['publisher_md5']=exp or '';r['publisher_md5_match']=(r.get('md5')==exp) if exp else ''
fields=list(dict.fromkeys(k for r in rows for k in r))
with (B/'tables/RR10_download_provenance.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
