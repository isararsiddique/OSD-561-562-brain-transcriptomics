"""Stage 4 helper: fetch full text (Europe PMC fullTextXML or a URL) via the capped client and print
only sentences matching a regex, so the screener reads the decisive passages."""
import re, sys, html
sys.path.insert(0, __import__('os').path.dirname(__file__))
import hsearch as h


def text_of(src):
    if src.startswith('PMC'):
        t, _ = h.get('https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML' % src)
    else:
        t, _ = h.get(src)
    if not isinstance(t, str):
        t = str(t)
    t = re.sub(r'<[^>]+>', ' ', t)
    return re.sub(r'\s+', ' ', html.unescape(t))


def grep(src, pat, width=0, maxn=40):
    t = text_of(src)
    sents = re.split(r'(?<=[.;!?])\s+(?=[A-Z(\[])', t)
    rx = re.compile(pat, re.I)
    n = 0
    for i, s in enumerate(sents):
        if rx.search(s):
            ctx = ' '.join(sents[max(0, i - width):i + width + 1])
            print('-', ctx[:700]); n += 1
            if n >= maxn:
                break
    print('[%d chars, %d matches shown]' % (len(t), n))


if __name__ == '__main__':
    grep(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 0,
         int(sys.argv[4]) if len(sys.argv) > 4 else 40)
