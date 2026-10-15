"""Correct a broad query and record date-bounded supplementary literature discovery."""
from pathlib import Path
import json,csv,urllib.parse,urllib.request,hashlib,datetime,concurrent.futures
B=Path(__file__).resolve().parents[1]
qs=['(TITLE_ABS:spaceflight OR TITLE_ABS:"space flight") AND (TITLE_ABS:cerebellum OR TITLE_ABS:cerebellar OR TITLE_ABS:brain) AND (TITLE_ABS:mouse OR TITLE_ABS:mice OR TITLE_ABS:murine) AND FIRST_PDATE:[1900-01-01 TO 2026-10-15]', '("RR-10" OR "RR-23" OR "RRRM-2" OR "SpaceX-16") AND (spaceflight OR "space flight") AND (cerebellum OR cerebellar) AND FIRST_PDATE:[1900-01-01 TO 2026-10-15]']
def fetch(x):
 i,q=x;url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode(dict(query=q,format='json',resultType='core',pageSize=1000));b=urllib.request.urlopen(url,timeout=55).read();p=B/'sources'/f'EuropePMC_supplement_{i+1}.json';p.write_bytes(b);d=json.loads(b);return dict(file=p.name,url=url,query=q,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status=200,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),error='',hitCount=d['hitCount'],returned=len(d['resultList']['result']))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:rows=list(e.map(fetch,enumerate(qs)))
with (B/'tables/literature_supplement_retrieval_log.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print([(x['hitCount'],x['returned']) for x in rows])
