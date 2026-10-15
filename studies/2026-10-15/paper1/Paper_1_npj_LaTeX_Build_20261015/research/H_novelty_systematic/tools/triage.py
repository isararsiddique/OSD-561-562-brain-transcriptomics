"""Stage 3a: registered rule-based title/abstract pass (>=2 concept blocks of the claim, or any C4 accession/mission term)."""
import json, re, collections
R = lambda p: re.compile(p, re.I)
B = {
 'SPACE': R(r'space ?flight|microgravity|international space station|rodent research|\bRRRM|rodent habitat|space station|\bISS\b|on-?orbit'),
 'COLLECT': R(r'euthan|sample collection|collection (technique|protocol|method)|dissect|live animal return|tissue preservation|carcass|ketamine|xylazine|isoflurane|an(a)?esthe|pentobarbital|decapitat|post-?mortem|sacrific|harvest|tissue collection|confound|artifact|artefact'),
 'BRAIN': R(r'brain|hippocamp|cerebell|cortex|cortical|nervous system|neural|neuron'),
 'IEG': R(r'immediate[- ]early|primary response gene|c-?fos|\bfos\b|egr-?1|npas4|junb|fosb|nr4a1|\barc\b'),
 'OMICS': R(r'transcriptom|rna-?seq|rna sequencing|gene expression|mrna|microarray|genomic|omics|methylation|proteom|single-cell'),
 'CORR': R(r'batch|unwanted variation|\bRUV|surrogate variable|combat|confounder|latent (factor|confound)|hidden (factor|confound)|harmoni[sz]'),
 'GAIN': R(r'operator norm|spectral norm|spectral gain|amplif|oblique projection|lipschitz|singular value|variance inflation|condition number|inflat'),
 'LOAD': R(r'loading|precision|imprecis|estimation error|uncertaint|sample size|detectab|removab|removal|false positive|exaggerat|over-?correct|pitfall'),
 'FW': R(r'frank-? ?wolfe|conditional gradient|duality gap|convex|constrain'),
 'HOLD': R(r'hold-?out|held-?out|cross-?validation|leakage|negative control|benchmark|validation|noise ceiling'),
}
ACCM = R(r'\b(OSD|GLDS)-?(561|562|563|564|612|613|525|556|557|568|569|588|589)\b|RRRM-?2|rodent research reference mission|\bRR-?(10|23)\b|rodent research[- ](10|23)\b')
CLB = {'C1': ['SPACE', 'COLLECT', 'BRAIN', 'IEG', 'OMICS'], 'C2': ['CORR', 'GAIN', 'OMICS', 'LOAD'],
       'C3': ['FW', 'CORR', 'SPACE', 'HOLD', 'GAIN'], 'C4': ['SPACE', 'BRAIN', 'OMICS']}
out = []
for l in open('../cache/unique.jsonl'):
    u = json.loads(l); t = (u['title'] or '') + ' ' + (u['abstract'] or '')
    hitb = {b for b in B if B[b].search(t)}
    passed, why = [], []
    for c in u['claims']:
        n = [b for b in CLB[c] if b in hitb]
        if len(n) >= 2 or (c == 'C4' and ACCM.search(t)):
            passed.append(c); why.append('%s:%s' % (c, '+'.join(n) + ('+ACC' if ACCM.search(t) else '')))
    u['blocks'] = sorted(hitb); u['pass3a'] = passed; u['why3a'] = ';'.join(why)
    out.append(u)
with open('../cache/triaged.jsonl', 'w') as f:
    for u in out: f.write(json.dumps(u, ensure_ascii=False) + '\n')
print('pass', sum(1 for u in out if u['pass3a']), 'of', len(out),
      collections.Counter(c for u in out for c in u['pass3a']))
