"""Stage 2: deduplicate records by DOI, then PMID, then PMCID, then normalised title+year."""
import json, re, collections
recs = [json.loads(l) for l in open('../cache/records.jsonl')]
def nt(t): return re.sub(r'[^a-z0-9]', '', (t or '').lower())[:120]
groups, key2g = [], {}
for r in recs:
    keys = [k for k in ['doi:' + r['doi'] if r.get('doi') else None,
                        'pmid:' + str(r['pmid']) if r.get('pmid') else None,
                        'pmcid:' + r['pmcid'] if r.get('pmcid') else None,
                        'ty:%s:%s' % (nt(r['title']), r.get('year')) if nt(r['title']) else None] if k]
    g = next((key2g[k] for k in keys if k in key2g), None)
    if g is None:
        g = len(groups); groups.append([])
    groups[g].append(r)
    for k in keys: key2g[k] = g
out = []
for i, g in enumerate(groups):
    best = max(g, key=lambda r: len(r.get('abstract') or ''))
    out.append({'uid': 'H%05d' % i, 'doi': next((r['doi'] for r in g if r.get('doi')), ''),
                'pmid': next((str(r['pmid']) for r in g if r.get('pmid')), ''),
                'pmcid': next((r['pmcid'] for r in g if r.get('pmcid')), ''),
                'title': best['title'], 'year': best['year'], 'journal': best.get('journal'),
                'pubtype': best.get('pubtype'), 'abstract': best.get('abstract') or '',
                'sources': sorted({r['db'] for r in g}), 'qids': sorted({r['qid'] for r in g}),
                'claims': sorted({r['claim'] for r in g}), 'epmc_src': next((r['src'] for r in g if r.get('src')), '')})
with open('../cache/unique.jsonl', 'w') as f:
    for o in out: f.write(json.dumps(o, ensure_ascii=False) + '\n')
c = collections.Counter(cl for o in out for cl in o['claims'])
print('identified', len(recs), 'unique', len(out), 'per-claim unique', dict(c),
      'by db', collections.Counter(r['db'] for r in recs))
