#!/usr/bin/env python3
"""Worker C discovery helper: Europe PMC / Crossref keyword searches with on-disk cache.

Used only to FIND candidate DOIs. Every candidate is then verified by c_verify.py
(Crossref by DOI + Europe PMC by DOI + Crossref 'updates' filter + doi.org handle).

Usage: python c_search.py epmc "<query>" [n]
       python c_search.py crossref "<bibliographic query>" [n]
"""
import hashlib
import json
import pathlib
import sys
import time
import urllib.parse

import requests

HERE = pathlib.Path(__file__).resolve().parent.parent
CACHE = HERE / "cache" / "search"
CACHE.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Paper1-bibliography-verification/1.0 (read-only; python-requests)"}
_last = [0.0]


def _get(url, params):
    key = hashlib.sha1((url + json.dumps(params, sort_keys=True)).encode()).hexdigest()
    f = CACHE / f"{key}.json"
    if f.exists():
        return json.loads(f.read_text())
    for attempt in range(6):
        wait = 1.2 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        r = requests.get(url, params=params, headers=UA, timeout=60)
        _last[0] = time.time()
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(min(60, 5 * 2 ** attempt))
            continue
        break
    r.raise_for_status()
    data = r.json()
    f.write_text(json.dumps({"url": url, "params": params, "response": data}))
    return {"url": url, "params": params, "response": data}


def epmc(query, n=10):
    d = _get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
             {"query": query, "format": "json", "pageSize": n, "resultType": "lite"})
    out = []
    for r in d["response"].get("resultList", {}).get("result", []):
        out.append(f"{r.get('pubYear')} | {r.get('authorString','')[:60]} | {r.get('title','')[:110]} | "
                   f"{r.get('journalTitle','')} | doi={r.get('doi')} | pmid={r.get('pmid')} | src={r.get('source')} | cited={r.get('citedByCount')}")
    return d["response"].get("hitCount"), out


def crossref(query, n=8):
    d = _get("https://api.crossref.org/works",
             {"query.bibliographic": query, "rows": n,
              "select": "DOI,title,container-title,issued,author,type"})
    out = []
    for it in d["response"]["message"]["items"]:
        au = it.get("author", [{}])
        first = au[0].get("family", "") if au else ""
        yr = it.get("issued", {}).get("date-parts", [[None]])[0][0]
        out.append(f"{yr} | {first} | {(it.get('title') or [''])[0][:110]} | "
                   f"{(it.get('container-title') or [''])[0]} | {it.get('type')} | doi={it.get('DOI')}")
    return None, out


if __name__ == "__main__":
    mode, q = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    hits, rows = (epmc if mode == "epmc" else crossref)(q, n)
    print(f"## {mode}: {q}  (hits={hits})")
    for row in rows:
        print("  " + row)
