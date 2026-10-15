#!/usr/bin/env python3
"""Build ../cache/manual_texts.json: supporting text for references whose abstract is not carried by
Crossref / Europe PMC / OpenAlex, or whose manuscript claim rests on full text.

Every text is extracted programmatically from a file cached under ../cache/ (publisher HTML, Semantic
Scholar JSON), except Szyld 2006, whose Springer landing page serves a bot wall to scripted clients; that
abstract was read with the interactive web_fetch tool on 15 Oct 2026 from
https://link.springer.com/article/10.1007/s11075-006-9046-2 and is transcribed verbatim below (flagged).
"""
import html
import json
import pathlib
import re

C = pathlib.Path(__file__).resolve().parent.parent / "cache"


def plain(t):
    t = re.sub(r"<script.*?</script>", " ", t, flags=re.S)
    t = re.sub(r"<style.*?</style>", " ", t, flags=re.S)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", html.unescape(t)).strip()


out = {}

# van den Brink 2017 (Nat. Methods Brief Communication; PubMed/Crossref carry no abstract): opening paragraph
t = (C / "publisher" / "10.1038_nmeth.4437.html").read_text()
paras = [plain(p) for p in re.findall(r"<p>(.*?)</p>", t, flags=re.S)]
lead = next(p for p in paras if p.startswith("In many gene expression studies"))
out["vandenBrink2017dissociation"] = [{"kind": "publisher page, opening paragraph (no indexed abstract)",
                                       "source": "https://www.nature.com/articles/nmeth.4437", "text": lead}]

# Masarapu 2024: manuscript claim about bulk concordance rests on full text, not the abstract
t = plain((C / "fulltext" / "Masarapu2024_nature.html").read_text())
i = t.find("Comparison of 825 spaceflight multiomics DEGs")
j = t.find("cell-specific analysis of central nervous system responses to spaceflight", i)
out["Masarapu2024spatial"] = [{"kind": "publisher full text (Results)",
                               "source": "https://www.nature.com/articles/s41467-024-48916-8",
                               "text": t[i:j + len("cell-specific analysis of central nervous system responses to spaceflight")] + "."}]

# O'Donoghue & Candes 2015: abstract from Semantic Scholar Graph API
s2 = json.loads((C / "semanticscholar" / "10.1007_s10208-013-9150-3.json").read_text())
out["ODonoghue2015restart"] = [{"kind": "abstract (Semantic Scholar Graph API)", "source": s2["url"],
                                "text": re.sub(r"\s+", " ", s2["json"]["abstract"]).strip()}]

# Jaggi 2013 (PMLR, no DOI): abstract from the PMLR proceedings page
t = (C / "nondoi" / "pmlr_jaggi13.html").read_text()
m = re.search(r'<div id="abstract"[^>]*>(.*?)</div>', t, flags=re.S)
out["Jaggi2013frankwolfe"] = [{"kind": "abstract (PMLR proceedings page)",
                               "source": "https://proceedings.mlr.press/v28/jaggi13.html", "text": plain(m.group(1))}]

# Szyld 2006: transcribed verbatim from the web_fetch tool output (see module docstring)
out["Szyld2006oblique"] = [{"kind": "abstract (Springer landing page via web_fetch tool; transcribed verbatim, "
                                    "tool tokenisation artefacts such as 'i. e. ,' left as received)",
                            "source": "https://link.springer.com/article/10.1007/s11075-006-9046-2",
                            "text": "Given an oblique projector P on a Hilbert space, i. e. , an operator satisfying P ²=P, "
                                    "which is neither null nor the identity, it holds that ||P|| = ||I –P||. This useful "
                                    "equality, while not widely-known, has been proven repeatedly in the literature. Many "
                                    "published proofs are reviewed, and simpler ones are presented."}]

(C / "manual_texts.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
for k, v in out.items():
    print(k, "|", v[0]["text"][:160])
