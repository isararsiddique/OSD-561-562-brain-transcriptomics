from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Resolve RR8 GEO brain libraries to primary NASA OSDR records, without new replication."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,re,datetime
import requests,pandas as pd
B=stage_root('extended')
def retrieve(sid):
    out={}
    for kind in ['meta','files']:
        url=f'https://osdr.nasa.gov/osdr/data/osd/{kind}/{sid}';p=B/f'sources/OSD-{sid}_{kind}.json'
        if not p.exists():
            r=request_get(url,timeout=50);r.raise_for_status();p.write_bytes(r.content)
        out[kind]=json.loads(p.read_text())
    s=out['meta']['study'][f'OSD-{sid}']['studies'][0]
    gsms=sorted(set(re.findall(r'GSM\d+',json.dumps(out))))
    return {'OSD':f'OSD-{sid}','title':s['title'],'assay_platforms':'; '.join(str(a['technologyPlatform']) for a in s['assays']),'measurement_types':'; '.join(str(a['measurementType']) for a in s['assays']),'linked_GSMs':gsms,'meta_SHA256':hashlib.sha256((B/f'sources/OSD-{sid}_meta.json').read_bytes()).hexdigest()}
with ThreadPoolExecutor(max_workers=3) as e:records=list(e.map(retrieve,range(927,933)))
reg=pd.read_csv(B/'tables/GSE295428_current_sample_registry.csv');rows=[]
for r in reg.itertuples():
    hits=[x for x in records if r.sample in x['linked_GSMs']]
    rows.append({'sample':r.sample,'tissue':r.tissue,'donor':r.donor,'flight':r.flight,'age':r.age,'matched_OSDR_accessions':';'.join(x['OSD'] for x in hits),'unique_primary_GSM_link':len(hits)==1})
pd.DataFrame(rows).to_csv(B/'tables/RR8_NASA_GEO_library_links.csv',index=False)
(B/'sources/RR8_NASA_GEO_accession_audit.json').write_text(json.dumps(records,indent=2))
print(pd.DataFrame(rows).groupby(['tissue','matched_OSDR_accessions']).size().to_string());print('Unique primary GSM links',sum(r['unique_primary_GSM_link'] for r in rows),'of',len(rows))
