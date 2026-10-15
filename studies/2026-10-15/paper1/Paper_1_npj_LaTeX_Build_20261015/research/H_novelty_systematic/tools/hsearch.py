"""H novelty search: rate-limited, cached, logged, capped read-only API client.

Every HTTP GET counts against a hard cap of 600 calls (state in cache/call_count.txt).
Records are streamed to cache/records.jsonl; one log line per query to cache/search_log.jsonl.
"""
import hashlib, json, os, time, urllib.parse, urllib.request, datetime

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(OUT, 'cache')
os.makedirs(CACHE, exist_ok=True)
CNT = os.path.join(CACHE, 'call_count.txt')
REC = os.path.join(CACHE, 'records.jsonl')
LOG = os.path.join(CACHE, 'search_log.jsonl')
CAP = 600
_last = [0.0]
DATE_EPMC = 'FIRST_PDATE:[1900-01-01 TO 2026-10-15]'


def _count():
    return int(open(CNT).read()) if os.path.exists(CNT) else 0


def get(url, params=None):
    """GET with cache; returns (parsed_json_or_text, from_cache)."""
    full = url + ('?' + urllib.parse.urlencode(params) if params else '')
    key = hashlib.sha256(full.encode()).hexdigest()[:24]
    path = os.path.join(CACHE, 'raw', key)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        txt = open(path, encoding='utf-8').read()
        return (json.loads(txt) if txt[:1] in '{[' else txt), True
    n = _count()
    if n >= CAP:
        raise RuntimeError('API call cap of %d reached' % CAP)
    wait = 0.35 - (time.time() - _last[0])
    if wait > 0:
        time.sleep(wait)
    req = urllib.request.Request(full, headers={'User-Agent': 'H-novelty-search/1.0 (read-only literature search)'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                txt = r.read().decode('utf-8', 'replace')
            break
        except Exception as e:  # noqa
            if attempt == 2:
                open(CNT, 'w').write(str(n + 1))
                raise
            time.sleep(2 + 3 * attempt)
    _last[0] = time.time()
    open(CNT, 'w').write(str(n + 1))
    with open(path, 'w', encoding='utf-8') as f:
        f.write(txt)
    with open(os.path.join(CACHE, 'raw_index.tsv'), 'a') as f:
        f.write(key + '\t' + full + '\n')
    return (json.loads(txt) if txt[:1] in '{[' else txt), False


def _log(entry):
    entry['run_time_utc'] = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')


def _emit(rec):
    with open(REC, 'a', encoding='utf-8') as f:
        f.write(json.dumps(rec, ensure_ascii=False) + '\n')


def epmc(qid, claim, query, cap=2000, date=True, note=''):
    q = '(' + query + ')' + (' AND ' + DATE_EPMC if date else '')
    base = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
    cursor, got, hit, ids, calls = '*', 0, None, [], 0
    while True:
        d, _ = get(base, {'query': q, 'format': 'json', 'resultType': 'core', 'pageSize': 1000,
                          'cursorMark': cursor, 'synonym': 'false'})
        calls += 1
        hit = d.get('hitCount', 0)
        res = d.get('resultList', {}).get('result', [])
        for r in res:
            rec = {'db': 'EuropePMC', 'qid': qid, 'claim': claim, 'id': r.get('id'), 'src': r.get('source'),
                   'pmid': r.get('pmid'), 'pmcid': r.get('pmcid'), 'doi': (r.get('doi') or '').lower(),
                   'title': r.get('title'), 'year': r.get('pubYear'),
                   'journal': (r.get('journalInfo') or {}).get('journal', {}).get('title') or r.get('bookOrReportDetails', {}).get('publisher'),
                   'pubtype': ';'.join((r.get('pubTypeList') or {}).get('pubType', [])),
                   'abstract': (r.get('abstractText') or '')[:4000]}
            _emit(rec); ids.append(rec['doi'] or rec['id']); got += 1
        nxt = d.get('nextCursorMark')
        if not res or not nxt or nxt == cursor or got >= cap or got >= hit:
            break
        cursor = nxt
    _log({'qid': qid, 'claim': claim, 'database': 'Europe PMC', 'string': query,
          'filters': DATE_EPMC if date else 'none', 'hit_count': hit, 'retrieved': got,
          'truncated': got < hit, 'api_calls': calls, 'note': note, 'records': ids})
    return hit, got


def epmc_count(query, date=True):
    q = '(' + query + ')' + (' AND ' + DATE_EPMC if date else '')
    d, _ = get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',
               {'query': q, 'format': 'json', 'resultType': 'lite', 'pageSize': 25, 'synonym': 'false'})
    return d.get('hitCount', 0), d.get('resultList', {}).get('result', [])


def _inv(ii):
    if not ii:
        return ''
    pos = {}
    for w, ps in ii.items():
        for p in ps:
            pos[p] = w
    return ' '.join(pos[k] for k in sorted(pos))[:4000]


def openalex(qid, claim, filt=None, search=None, cap=200, note=''):
    params = {'per-page': 200, 'select': 'id,doi,title,publication_year,abstract_inverted_index,primary_location,type,ids'}
    f = 'to_publication_date:2026-10-15'
    if filt:
        f = filt + ',' + f
    params['filter'] = f
    if search:
        params['search'] = search
    d, _ = get('https://api.openalex.org/works', params)
    hit = d.get('meta', {}).get('count', 0)
    ids = []
    for w in d.get('results', [])[:cap]:
        src = ((w.get('primary_location') or {}).get('source') or {}) or {}
        rec = {'db': 'OpenAlex', 'qid': qid, 'claim': claim, 'id': w.get('id'),
               'pmid': ((w.get('ids') or {}).get('pmid') or '').replace('https://pubmed.ncbi.nlm.nih.gov/', ''),
               'pmcid': '', 'doi': (w.get('doi') or '').replace('https://doi.org/', '').lower(),
               'title': w.get('title'), 'year': w.get('publication_year'), 'journal': src.get('display_name'),
               'pubtype': w.get('type'), 'abstract': _inv(w.get('abstract_inverted_index'))}
        _emit(rec); ids.append(rec['doi'] or rec['id'])
    _log({'qid': qid, 'claim': claim, 'database': 'OpenAlex', 'string': search or filt,
          'filters': f, 'hit_count': hit, 'retrieved': len(ids), 'truncated': len(ids) < hit,
          'api_calls': 1, 'note': note, 'records': ids})
    return hit, len(ids)


def crossref(qid, claim, q, rows=100, note=''):
    d, _ = get('https://api.crossref.org/works', {'query.bibliographic': q, 'rows': rows,
               'filter': 'until-pub-date:2026-10-15',
               'select': 'DOI,title,issued,container-title,type,abstract'})
    m = d.get('message', {})
    hit = m.get('total-results', 0)
    ids = []
    for w in m.get('items', []):
        y = ((w.get('issued') or {}).get('date-parts') or [[None]])[0][0]
        rec = {'db': 'Crossref', 'qid': qid, 'claim': claim, 'id': w.get('DOI'), 'pmid': '', 'pmcid': '',
               'doi': (w.get('DOI') or '').lower(), 'title': (w.get('title') or [''])[0], 'year': y,
               'journal': (w.get('container-title') or [''])[0], 'pubtype': w.get('type'),
               'abstract': (w.get('abstract') or '')[:4000]}
        _emit(rec); ids.append(rec['doi'])
    _log({'qid': qid, 'claim': claim, 'database': 'Crossref', 'string': q,
          'filters': 'until-pub-date:2026-10-15; relevance top %d' % rows, 'hit_count': hit,
          'retrieved': len(ids), 'truncated': True, 'api_calls': 1, 'note': note, 'records': ids})
    return hit, len(ids)


def log_misc(entry):
    _log(entry)
