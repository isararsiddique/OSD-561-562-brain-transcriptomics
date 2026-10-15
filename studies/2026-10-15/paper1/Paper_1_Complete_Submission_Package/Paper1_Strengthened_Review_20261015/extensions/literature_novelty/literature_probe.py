"""Read-only public literature search; cache exact responses for traceability."""
import concurrent.futures
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).parent
CACHE = ROOT / "cache"
QUERIES = {
    "masarapu_2024": "DOI:10.1038/s41467-024-48916-8",
    "collection_2022": "DOI:10.1016/j.crmeth.2022.100325",
    "holley_2022": "DOI:10.1038/s41526-022-00217-4",
    "marsh_2022": "DOI:10.1038/s41593-022-01022-8",
    "dachet_2021": "DOI:10.1038/s41598-021-85801-6",
    "space_brain_circadian": '(spaceflight OR "space flight") AND (brain OR hippocamp*) AND circadian',
    "mouse_space_brain": '(spaceflight OR "space flight") AND (mouse OR mice) AND (brain OR hippocamp*) FIRST_PDATE:[2020-01-01 TO 2026-10-15]',
    "batch_confidence": 'TITLE_ABS:("batch effects") AND TITLE_ABS:("exaggerated confidence")',
    "batch_confounding": 'TITLE_ABS:("batch effect") AND TITLE_ABS:(confounding)',
    "ground_control_space": '(spaceflight OR "space flight") AND ("ground control" OR "habitat control") AND (transcriptomic OR transcriptome)',
}

def get(name, query):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": query, "format": "json", "resultType": "core", "pageSize": 100})
    req = urllib.request.Request(url, headers={"User-Agent": "OSDR-literature-review/1.0 (public secondary analysis)"})
    raw = urllib.request.urlopen(req, timeout=90).read()
    (CACHE / f"{name}.json").write_bytes(raw)
    obj = json.loads(raw)
    print(name, "hits", obj["hitCount"])
    records = obj["resultList"]["result"]
    for r in records[:25]:
        print(r.get("id"), r.get("pmcid"), r.get("doi"), r.get("title"))
    return {"name": name, "query": query, "url": url, "retrieved_date": "2026-10-15", "sha256": hashlib.sha256(raw).hexdigest(), "hit_count": obj["hitCount"], "records_returned": len(records)}

if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        provenance = list(ex.map(lambda item: get(*item), QUERIES.items()))
    (ROOT / "search_provenance.json").write_text(json.dumps(provenance, indent=2))
