"""Reproduce the bounded NASA OSDR / Europe PMC search without outcome selection."""
from pathlib import Path
import csv, json, hashlib, urllib.request, urllib.parse, datetime, concurrent.futures
B=Path(__file__).resolve().parents[1]; S=B/'sources'; T=B/'tables'; S.mkdir(exist_ok=True); T.mkdir(exist_ok=True)
JOBS=[]
for term in ['cerebellum','cerebellar','brain']:
 JOBS.append((f'OSDR_search_{term}.json','https://osdr.nasa.gov/osdr/data/search?'+urllib.parse.urlencode({'term':term,'size':10000}),'NASA keyword '+term))
for i,q in enumerate(['(cerebellum OR cerebellar) AND (spaceflight OR "space flight") AND (mouse OR mice OR murine)','(RR-10 OR RR-23 OR RRRM-2 OR "SpaceX-16") AND (cerebellum OR cerebellar)']):
 JOBS.append((f'EuropePMC_search_{i+1}.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode(dict(query=q,format='json',resultType='core',pageSize=1000)),q))
# Metadata-confirmed candidate classes; fetch all before eligibility judgement.
for n in [33,222,352,525,536,561,562,563,564,612,613,682,685,698,699,916,927,928,929,930,931,932]:
 for k in ['meta','files']:JOBS.append((f'OSD-{n}_{k}.json',f'https://osdr.nasa.gov/osdr/data/osd/{k}/{n}',f'Candidate metadata/inventory OSD-{n}'))
def fetch(job):
 name,url,query=job; p=S/name; row=dict(file=name,url=url,query=query,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='',bytes=0,sha256='',error='')
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 scientific metadata audit'})
  with urllib.request.urlopen(req,timeout=55) as r:b=r.read();row['status']=r.status
  json.loads(b);p.write_bytes(b);row.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
 except Exception as e:row['status']='ERROR';row['error']=str(e)
 print(name,row['status'],row['bytes'],flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as e:rows=list(e.map(fetch,JOBS))
with (T/'retrieval_log.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Finished',len(rows))
