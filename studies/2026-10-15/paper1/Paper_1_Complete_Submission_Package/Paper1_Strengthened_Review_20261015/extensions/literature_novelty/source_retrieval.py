import concurrent.futures
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT=Path(__file__).parent
SOURCES={
 'masarapu_publisher': 'https://www.nature.com/articles/s41467-024-48916-8',
 'collection_fulltext': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9701605/fullTextXML',
 'holley_fulltext': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9365836/fullTextXML',
 'marsh_fulltext': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11645269/fullTextXML',
 'dachet_fulltext': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7988150/fullTextXML',
 'nygaard_fulltext': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4679072/fullTextXML',
 'batch_causal_fulltext': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12319767/fullTextXML',
}
DOIS=['10.1038/s41467-024-48916-8','10.1016/j.crmeth.2022.100325','10.1038/s41526-022-00217-4','10.1038/s41593-022-01022-8','10.1038/s41598-021-85801-6','10.1093/biostatistics/kxv027','10.1093/biostatistics/kxr034','10.1038/nbt.2931','10.1371/journal.pgen.0030161','10.1016/j.neuron.2018.04.001','10.1038/nrg.2016.150','10.1016/j.jgg.2019.08.002','10.1162/imag_a_00458','10.3390/cells14060417','10.1038/s41467-026-68737-1']
for doi in DOIS:
 SOURCES['crossref_'+doi.replace('/','_')]='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='')
def get(item):
 name,url=item
 try:
  raw=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'OSDR-literature-review/1.0'}),timeout=90).read()
  ext='.json' if name.startswith('crossref') else ('.html' if name.endswith('publisher') else '.xml')
  (ROOT/'cache'/(name+ext)).write_bytes(raw)
  print(name,len(raw))
  return {'name':name,'url':url,'retrieved_date':'2026-10-15','sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'ok':True}
 except Exception as e:
  print(name,str(e)); return {'name':name,'url':url,'ok':False,'error':str(e)}
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: results=list(ex.map(get,SOURCES.items()))
 (ROOT/'source_provenance.json').write_text(json.dumps(results,indent=2))
