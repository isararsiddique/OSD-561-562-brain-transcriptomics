"""Read-only literature/metadata query helper for Worker B (novelty / prior art).

Every call is appended to ../B_search_log.json with endpoint, query, URL, UTC time,
HTTP status, hit count and the top hits, and the raw response is cached under ../cache/.
No project data is ever sent; only public search strings.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import time
import urllib.parse
from datetime import datetime, timezone

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LOG = os.path.join(ROOT, "B_search_log.json")
CACHE = os.path.join(ROOT, "cache")
os.makedirs(CACHE, exist_ok=True)
UA = {"User-Agent": "Paper1-novelty-check/1.0 (read-only literature search)"}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_last = [0.0]


def _sleep():
    dt = time.time() - _last[0]
    if dt < 0.4:
        time.sleep(0.4 - dt)
    _last[0] = time.time()


def _load_log():
    if os.path.exists(LOG):
        with open(LOG) as fh:
            return json.load(fh)
    return []


def _append(entry):
    log = _load_log()
    entry["n"] = len(log) + 1
    log.append(entry)
    with open(LOG, "w") as fh:
        json.dump(log, fh, indent=1, ensure_ascii=False)


def _get(url, params=None, accept_json=True, tries=4):
    err = None
    for k in range(tries):
        _sleep()
        try:
            r = requests.get(url, params=params, headers=UA, timeout=60)
            if r.status_code in (429, 500, 502, 503, 504):
                err = f"HTTP {r.status_code}"
                time.sleep(2 * (k + 1))
                continue
            return r
        except requests.RequestException as e:  # network error
            err = repr(e)
            time.sleep(2 * (k + 1))
    raise RuntimeError(err)


def _cache_name(tag, text):
    h = hashlib.sha256(text.encode()).hexdigest()[:12]
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", tag)[:60]
    return os.path.join(CACHE, f"{safe}_{h}")


def epmc(query, tag, purpose="", page_size=25, synonym=False):
    """Europe PMC search. Full-text is searched by default by Europe PMC for OA records;
    use field prefixes (e.g. TITLE_ABS:) to restrict."""
    params = {"query": query, "format": "json", "pageSize": page_size,
              "resultType": "core", "synonym": str(synonym).lower()}
    # NB: adding sort=RELEVANCE makes this endpoint return HTTP 503 (tested 15 Oct 2026);
    # the default order is already relevance-ranked.
    url = EPMC + "/search?" + urllib.parse.urlencode(params)
    entry = {"date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "endpoint": "EuropePMC REST /search", "tag": tag, "query": query,
             "url": url, "purpose": purpose}
    try:
        r = _get(EPMC + "/search", params)
        entry["http_status"] = r.status_code
        d = r.json()
        fn = _cache_name("epmc_" + tag, url) + ".json"
        with open(fn, "w") as fh:
            json.dump(d, fh)
        entry["cache"] = os.path.relpath(fn, ROOT)
        entry["hit_count"] = int(d.get("hitCount", 0))
        hits = []
        for x in d.get("resultList", {}).get("result", []):
            hits.append({"pmid": x.get("pmid"), "pmcid": x.get("pmcid"), "doi": x.get("doi"),
                         "year": x.get("pubYear"), "journal": x.get("journalTitle") or (x.get("journalInfo") or {}).get("journal", {}).get("title"),
                         "title": x.get("title"), "first_author": (x.get("authorString") or "").split(",")[0]})
        entry["top_hits"] = hits[:page_size]
    except Exception as e:  # record failures plainly
        entry["http_status"] = None
        entry["error"] = str(e)
        d = None
    _append(entry)
    return entry, d


def pubmed(term, tag, purpose="", retmax=20):
    params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": retmax}
    url = EUTILS + "/esearch.fcgi?" + urllib.parse.urlencode(params)
    entry = {"date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "endpoint": "NCBI E-utilities esearch (pubmed)", "tag": tag, "query": term,
             "url": url, "purpose": purpose}
    try:
        r = _get(EUTILS + "/esearch.fcgi", params)
        entry["http_status"] = r.status_code
        d = r.json()
        es = d.get("esearchresult", {})
        entry["hit_count"] = int(es.get("count", 0))
        ids = es.get("idlist", [])
        entry["pmids"] = ids
        if ids:
            r2 = _get(EUTILS + "/esummary.fcgi", {"db": "pubmed", "id": ",".join(ids), "retmode": "json"})
            s = r2.json().get("result", {})
            entry["top_hits"] = [{"pmid": i, "title": s.get(i, {}).get("title"),
                                  "year": (s.get(i, {}).get("pubdate") or "")[:4],
                                  "journal": s.get(i, {}).get("source"),
                                  "doi": next((a.get("value") for a in s.get(i, {}).get("articleids", []) if a.get("idtype") == "doi"), None)}
                                 for i in ids]
    except Exception as e:
        entry["http_status"] = None
        entry["error"] = str(e)
        d = None
    _append(entry)
    return entry, d


def fetch(url, tag, purpose="", ext="txt", params=None):
    """Fetch a single record / full text / metadata page and cache it."""
    entry = {"date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "endpoint": "GET", "tag": tag, "url": url + (("?" + urllib.parse.urlencode(params)) if params else ""),
             "purpose": purpose}
    text = None
    try:
        r = _get(url, params)
        entry["http_status"] = r.status_code
        text = r.text
        fn = _cache_name("get_" + tag, entry["url"]) + "." + ext
        with open(fn, "w") as fh:
            fh.write(text)
        entry["cache"] = os.path.relpath(fn, ROOT)
        entry["bytes"] = len(r.content)
        entry["sha256"] = hashlib.sha256(r.content).hexdigest()
    except Exception as e:
        entry["http_status"] = None
        entry["error"] = str(e)
    _append(entry)
    return entry, text


def epmc_fulltext(pmcid, tag, purpose=""):
    return fetch(f"{EPMC}/{pmcid}/fullTextXML", tag, purpose, ext="xml")


def ncbi_pmc(pmc_num, tag, purpose=""):
    return fetch(EUTILS + "/efetch.fcgi", tag, purpose, ext="xml",
                 params={"db": "pmc", "id": str(pmc_num).replace("PMC", ""), "rettype": "xml"})


def pubmed_abstract(pmid, tag, purpose=""):
    return fetch(EUTILS + "/efetch.fcgi", tag, purpose, ext="xml",
                 params={"db": "pubmed", "id": str(pmid), "rettype": "abstract", "retmode": "xml"})


def strip_xml(x):
    x = re.sub(r"<xref[^>]*>.*?</xref>", "", x, flags=re.S)
    x = re.sub(r"<[^>]+>", " ", x)
    return re.sub(r"\s+", " ", x)


def sentences_with(text, patterns, width=400, maxn=8):
    out = []
    for p in patterns:
        for m in re.finditer(p, text, flags=re.I):
            s = max(0, m.start() - width // 2)
            out.append(text[s:m.end() + width // 2])
            if len(out) >= maxn:
                return out
    return out


def epmc_record(pmid=None, doi=None, tag="rec", purpose="read abstract at source"):
    q = f"EXT_ID:{pmid} AND SRC:MED" if pmid else f'DOI:"{doi}"'
    e, d = epmc(q, tag, purpose, page_size=3)
    rec = None
    if d:
        rs = d.get("resultList", {}).get("result", [])
        rec = rs[0] if rs else None
    return rec


def crossref(biblio, tag, purpose="", rows=3):
    params = {"query.bibliographic": biblio, "rows": rows,
              "select": "DOI,title,author,issued,container-title,type"}
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
    entry = {"date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "endpoint": "Crossref /works query.bibliographic", "tag": tag, "query": biblio,
             "url": url, "purpose": purpose}
    try:
        r = _get("https://api.crossref.org/works", params)
        entry["http_status"] = r.status_code
        d = r.json()
        items = d.get("message", {}).get("items", [])
        entry["hit_count"] = d.get("message", {}).get("total-results")
        entry["top_hits"] = [{"doi": i.get("DOI"), "title": (i.get("title") or [None])[0],
                              "year": ((i.get("issued") or {}).get("date-parts") or [[None]])[0][0],
                              "journal": (i.get("container-title") or [None])[0],
                              "first_author": ((i.get("author") or [{}])[0]).get("family")} for i in items]
    except Exception as e:
        entry["http_status"] = None
        entry["error"] = str(e)
    _append(entry)
    return entry
