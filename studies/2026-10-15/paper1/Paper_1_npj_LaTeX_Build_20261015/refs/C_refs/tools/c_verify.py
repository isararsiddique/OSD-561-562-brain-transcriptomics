#!/usr/bin/env python3
"""Worker C: fetch and verify every Paper 1 reference at source, then emit the bibliography files.

Sources (read-only, cached as JSON under ../cache/, polite rate, no e-mail address sent):
  * Crossref works API by DOI           -> bibliographic fields, type, update-to / relation
  * Crossref works?filter=updates:<DOI> -> correction / retraction / expression-of-concern notices
  * doi.org handle API                  -> DOI resolves (responseCode 1)
  * Europe PMC REST (resultType=core)   -> PMID, PMCID, abstract, pubType, commentCorrectionList
  * OpenAlex by DOI (abstract fallback only, when neither Crossref nor Europe PMC carries one)
  * DataCite + OSDR API                 -> NASA OSDR dataset DOIs (datasets kept in a separate list)

Usage:
  python c_verify.py fetch      # fetch + dump abstracts to ../C_abstracts_review.md
  python c_verify.py build      # verify quotes, write C_references.bib, C_datasets.bib, C_reference_audit.json,
                                # C_references_nature.txt, C_REPORT.md (report tables)
"""
import hashlib
import html
import json
import pathlib
import re
import sys
import time
import urllib.parse

import requests

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from c_entries import ENTRIES  # noqa: E402
from c_quotes import QUOTES, ABBR, NOTES, DROPPED  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent.parent
CACHE = OUT / "cache"
MANUSCRIPT = pathlib.Path(""
                          "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015/MANUSCRIPT_COMPLETE.md")
UA = {"User-Agent": "Paper1-bibliography-verification/1.0 (read-only; python-requests)"}
ACCESSED = "15 October 2026"
_last = {}


def _safe(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)[:150] + "_" + hashlib.sha1(s.encode()).hexdigest()[:8]


def get_json(kind, ident, url, params=None, spacing=1.2, allow_404=True):
    d = CACHE / kind
    d.mkdir(parents=True, exist_ok=True)
    f = d / (_safe(ident) + ".json")
    if f.exists():
        return json.loads(f.read_text())
    host = urllib.parse.urlparse(url).netloc
    for attempt in range(7):
        wait = spacing - (time.time() - _last.get(host, 0))
        if wait > 0:
            time.sleep(wait)
        try:
            r = requests.get(url, params=params, headers=UA, timeout=60)
        except requests.RequestException as e:  # transient network error
            _last[host] = time.time()
            time.sleep(min(60, 5 * 2 ** attempt))
            err = str(e)
            continue
        _last[host] = time.time()
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(min(90, 5 * 2 ** attempt))
            err = f"HTTP {r.status_code}"
            continue
        break
    else:
        raise RuntimeError(f"{kind} {ident}: gave up ({err})")
    rec = {"url": r.url, "status": r.status_code, "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if r.status_code == 404 and allow_404:
        rec["json"] = None
    else:
        r.raise_for_status()
        rec["json"] = r.json()
    f.write_text(json.dumps(rec))
    return rec


def clean(txt):
    if not txt:
        return ""
    txt = re.sub(r"<jats:title>.*?</jats:title>", " ", txt, flags=re.S)
    txt = re.sub(r"<[^>]+>", " ", txt)
    txt = html.unescape(txt)
    return re.sub(r"\s+", " ", txt).strip()


def fetch(doi):
    q = urllib.parse.quote(doi, safe="")
    cr = get_json("crossref", doi, f"https://api.crossref.org/works/{q}")
    upd = get_json("crossref_updates", doi, "https://api.crossref.org/works",
                   params={"filter": f"updates:{doi}", "rows": 50})
    hd = get_json("doi_handle", doi, f"https://doi.org/api/handles/{q}", spacing=0.5)
    ep = get_json("europepmc", doi, "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                  params={"query": f'DOI:"{doi}"', "resultType": "core", "format": "json", "pageSize": 5},
                  spacing=0.4)
    return cr, upd, hd, ep


def openalex(doi):
    q = urllib.parse.quote(f"https://doi.org/{doi}", safe=":/")
    return get_json("openalex", doi, f"https://api.openalex.org/works/{q}", spacing=0.5)


def oa_abstract(rec):
    j = rec.get("json") or {}
    inv = j.get("abstract_inverted_index")
    if not inv:
        return ""
    pos = {}
    for w, idxs in inv.items():
        for i in idxs:
            pos[i] = w
    return clean(" ".join(pos[i] for i in sorted(pos)))


def pick_epmc(ep):
    res = ((ep.get("json") or {}).get("resultList") or {}).get("result") or []
    med = [r for r in res if r.get("source") == "MED"]
    return (med or [r for r in res if r.get("source") != "PPR"] or [None])[0]


def year_of(m):
    for k in ("published-print", "journal-issue", "published-online", "issued"):
        v = m.get(k)
        if k == "journal-issue" and v:
            v = v.get("published-print")
        if v and v.get("date-parts") and v["date-parts"][0] and v["date-parts"][0][0]:
            return v["date-parts"][0][0], k
    return None, None


def first_year(m, k):
    v = m.get(k)
    try:
        return v["date-parts"][0][0]
    except Exception:
        return None


MANUAL = json.loads((CACHE / "manual_texts.json").read_text()) if (CACHE / "manual_texts.json").exists() else {}


def record_nondoi(e):
    """Archival proceedings without a DOI (PMLR): fields read from the cached proceedings page."""
    t = (CACHE / "nondoi" / "pmlr_jaggi13.html").read_text()
    meta = lambda k: re.findall(rf'<meta name="{k}" content="([^"]*)"', t)
    i = t.find("@InProceedings")
    bib = html.unescape(t[i:i + 900])
    gv = lambda k: (re.search(rf"{k}\s*=\s*\{{(.*?)\}},", bib, flags=re.S) or [None, None])[1]
    auth = []
    for a in meta("citation_author"):
        g, f = a.rsplit(" ", 1)
        auth.append({"family": f, "given": g, "is_group": False})
    return {
        "key": e["key"], "doi": None, "nondoi_id": e["nondoi"], "origin": e["origin"], "old_no": None,
        "topic": e.get("topic"), "section": e["section"], "supported_claim": e["claim"],
        "crossref_status": None, "doi_resolves": None, "crossref_type": "proceedings-article",
        "title": clean(meta("citation_title")[0]), "journal": gv("booktitle"), "publisher": gv("publisher"),
        "series": gv("series"), "volume": gv("volume"), "issue": gv("number"),
        "pages": (gv("pages") or "").replace("--", "-"), "pages_source": "PMLR page", "article_number": None,
        "year": int(gv("year")), "year_source": "PMLR page", "year_online": None, "year_print": None,
        "authors": auth, "n_authors_crossref": len(auth), "crossref_update_to_field": None,
        "crossref_relation": {}, "crossref_updated_by_field": None, "crossref_update_notices": [],
        "pmid": None, "pmcid": None, "epmc_source": None, "epmc_title": "", "epmc_journal": None,
        "epmc_volume": None, "epmc_issue": None, "epmc_year": None, "epmc_pages": None, "epmc_n_authors": None,
        "epmc_pubtypes": [], "epmc_comment_corrections": [], "abstract_epmc": "", "abstract_crossref": "",
        "abstract_openalex": "", "issn": None, "url": "https://proceedings.mlr.press/v28/jaggi13.html",
        "manual_texts": MANUAL.get(e["key"], []),
    }


def record(e):
    if e.get("nondoi"):
        return record_nondoi(e)
    doi = e["doi"]
    cr, upd, hd, ep = fetch(doi)
    m = (cr.get("json") or {}).get("message") or {}
    ep1 = pick_epmc(ep)
    abs_ep = clean((ep1 or {}).get("abstractText"))
    abs_cr = clean(m.get("abstract"))
    abs_oa = ""
    if not abs_ep and not abs_cr:
        abs_oa = oa_abstract(openalex(doi))
    notices = []
    for it in ((upd.get("json") or {}).get("message") or {}).get("items", []):
        for u in it.get("update-to", []):
            if u.get("DOI", "").lower() == doi.lower():
                notices.append({"notice_doi": it.get("DOI"), "type": u.get("type"),
                                "label": u.get("label"), "date": (u.get("updated") or {}).get("date-time")})
    ep_cc = []
    if ep1:
        for c in ((ep1.get("commentCorrectionList") or {}).get("commentCorrection") or []):
            ep_cc.append({"type": c.get("type"), "id": c.get("id"), "source": c.get("source"),
                          "reference": c.get("reference")})
    ep_types = (ep1 or {}).get("pubTypeList", {}).get("pubType", []) if ep1 else []
    if isinstance(ep_types, str):
        ep_types = [ep_types]
    yr, yr_src = year_of(m)
    pages, pages_src = m.get("page"), "Crossref page"
    if not pages and m.get("article-number"):
        pages, pages_src = m.get("article-number"), "Crossref article-number"
    if not pages and ep1 and ep1.get("pageInfo"):
        pages, pages_src = re.sub(r"^Article\s*", "", ep1["pageInfo"]), "Europe PMC pageInfo"
    s2f = CACHE / "semanticscholar" / (doi.replace("/", "_") + "_journal.json")
    if not pages and s2f.exists():
        pages, pages_src = (json.loads(s2f.read_text())["json"].get("journal") or {}).get("pages"), "Semantic Scholar journal.pages"
    if not pages:
        pages_src = "none available"
    authors = m.get("author", [])
    ji = (ep1 or {}).get("journalInfo") or {}
    rec = {
        "key": e["key"], "doi": doi, "origin": e["origin"], "old_no": e.get("old_no"),
        "topic": e.get("topic"), "section": e["section"], "supported_claim": e["claim"],
        "crossref_status": cr.get("status"), "doi_resolves": bool(hd.get("json") and hd["json"].get("responseCode") == 1),
        "crossref_type": m.get("type"),
        "title": clean((m.get("title") or [""])[0]),
        "journal": clean((m.get("container-title") or [""])[0]),
        "publisher": m.get("publisher"),
        "volume": m.get("volume"), "issue": m.get("issue"),
        "pages": pages, "pages_source": pages_src, "article_number": m.get("article-number"),
        "year": yr, "year_source": yr_src,
        "year_online": first_year(m, "published-online"), "year_print": first_year(m, "published-print"),
        "authors": [{"family": a.get("family") or a.get("name"), "given": a.get("given", ""),
                     "is_group": "family" not in a} for a in authors],
        "n_authors_crossref": len(authors),
        "crossref_update_to_field": m.get("update-to"),
        "crossref_relation": {k: v for k, v in (m.get("relation") or {}).items()},
        "crossref_updated_by_field": m.get("updated-by"),
        "crossref_update_notices": notices,
        "pmid": (ep1 or {}).get("pmid"), "pmcid": (ep1 or {}).get("pmcid"),
        "epmc_source": (ep1 or {}).get("source"),
        "epmc_title": clean((ep1 or {}).get("title")),
        "epmc_journal": (ji.get("journal") or {}).get("title"),
        "epmc_volume": ji.get("volume"), "epmc_issue": ji.get("issue"), "epmc_year": ji.get("yearOfPublication"),
        "epmc_pages": (ep1 or {}).get("pageInfo"),
        "epmc_n_authors": len(((ep1 or {}).get("authorList") or {}).get("author") or []) if ep1 else None,
        "epmc_pubtypes": ep_types,
        "epmc_comment_corrections": ep_cc,
        "abstract_epmc": abs_ep, "abstract_crossref": abs_cr, "abstract_openalex": abs_oa,
        "issn": m.get("ISSN"),
        "manual_texts": MANUAL.get(e["key"], []),
    }
    return rec


def retraction_status(r):
    bad = []
    for n in r["crossref_update_notices"]:
        t = (n.get("type") or "").lower()
        if "retract" in t or "withdraw" in t or "concern" in t or "removal" in t:
            bad.append(f"Crossref {n['type']} {n['notice_doi']}")
    for c in r["epmc_comment_corrections"]:
        t = (c.get("type") or "").lower()
        if "retract" in t or "concern" in t or "withdrawn" in t:
            bad.append(f"EuropePMC '{c['type']}' {c.get('id')}")
    if any("retracted" in p.lower() for p in r["epmc_pubtypes"]):
        bad.append("EuropePMC pubType Retracted Publication")
    corr = [f"Crossref {n['type']} {n['notice_doi']}" for n in r["crossref_update_notices"]
            if (n.get("type") or "").lower() in ("correction", "erratum", "corrigendum", "addendum")]
    corr += [f"EuropePMC '{c['type']}' {c.get('id')}" for c in r["epmc_comment_corrections"]
             if (c.get("type") or "").lower().startswith(("erratum", "correction", "corrected", "update"))]
    if bad:
        status = "RETRACTED/EOC: " + "; ".join(bad)
    elif corr:
        status = "not retracted; correction(s) recorded: " + "; ".join(corr)
    else:
        status = "no retraction, correction or expression of concern found (Crossref updates filter + Europe PMC)"
    return status, bool(bad)


# ------------------------------------------------------------------------- formatting
LATEX_MAP = {"&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#", "$": r"\$"}


def tex(s):
    s = s or ""
    return "".join(LATEX_MAP.get(ch, ch) for ch in s)


def initials(given):
    parts = re.split(r"([ -])", given.strip())
    out = ""
    for p in parts:
        if p == " ":
            out += " "
        elif p == "-":
            out = out.rstrip() + "-"
        elif p:
            if re.fullmatch(r"[A-Z]\.?", p):
                out += p[0] + "."
            elif p.endswith(".") and len(p) <= 3:
                out += p
            else:
                out += p[0].upper() + "."
    return re.sub(r"\s+", " ", out).strip()


def nature_authors(auth):
    def one(a):
        if a["is_group"]:
            return a["family"]
        return f"{a['family']}, {initials(a['given'])}" if a["given"] else a["family"]
    if len(auth) > 5:
        return one(auth[0]) + " et al."
    names = [one(a) for a in auth]
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + " & " + names[-1]


def pages_of(r):
    p = r.get("pages") or r.get("article_number") or ""
    parts = re.split(r"-+|\u2013", p) if p else []
    if len(parts) == 2 and parts[0] == parts[1]:
        p = parts[0]
    return p.replace("--", "\u2013").replace("-", "\u2013") if p else ""


def fix_chars(s):
    return (s or "").replace("\u2010", "-").replace("\u2011", "-")


def nature_line(r):
    title = fix_chars(r["title"]).rstrip(".")
    end = "" if title.endswith(("?", "!")) else "."
    abbr = ABBR.get(r["key"]) or ABBR.get(r["journal"]) or r["journal"]
    vol = r.get("volume") or ""
    pg = pages_of(r)
    core = f"{abbr} {vol}, {pg}" if vol else f"{abbr} {pg}"
    if r["crossref_type"] in ("book", "monograph"):
        return f"{nature_authors(r['authors'])} {title}{end} ({r['publisher']}, {r['year']}). https://doi.org/{r['doi']}"
    if not r["doi"]:
        return f"{nature_authors(r['authors'])} {title}{end} Proc. Mach. Learn. Res. {vol}, {pg} ({r['year']}). {r['url']}"
    return f"{nature_authors(r['authors'])} {title}{end} {core.strip().rstrip(',')} ({r['year']}). https://doi.org/{r['doi']}"


def bib_entry(r):
    au = " and ".join((("{" + tex(a["family"]) + "}") if a["is_group"] else f"{tex(a['family'])}, {tex(a['given'])}")
                      for a in r["authors"])
    pg = pages_of(r).replace("\u2013", "--")
    if not r["doi"]:
        f = [("author", au), ("title", "{" + tex(r["title"]) + "}"), ("booktitle", tex(r["journal"])),
             ("series", tex(r["series"])), ("volume", r["volume"]), ("pages", pg), ("year", str(r["year"])),
             ("publisher", r["publisher"]), ("url", r["url"])]
        body = ",\n".join(f"  {k} = {{{v}}}" for k, v in f if v)
        return f"@inproceedings{{{r['key']},\n{body}\n}}\n"
    if r["crossref_type"] in ("book", "monograph"):
        f = [("author", au), ("title", "{" + tex(fix_chars(r["title"])) + "}"), ("publisher", tex(r["publisher"])),
             ("year", str(r["year"])), ("doi", r["doi"])]
        typ = "book"
    else:
        f = [("author", au), ("title", "{" + tex(fix_chars(r["title"])) + "}"), ("journal", tex(r["journal"])),
             ("year", str(r["year"])), ("volume", r.get("volume") or ""), ("number", r.get("issue") or ""),
             ("pages", pg), ("doi", r["doi"])]
        if r.get("pmid"):
            f.append(("pmid", r["pmid"]))
        typ = "article"
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in f if v)
    return f"@{typ}{{{r['key']},\n{body}\n}}\n"


# ------------------------------------------------------------------------- old-reference comparison
def old_reference_texts():
    txt = MANUSCRIPT.read_text()
    sec = txt.split("## References", 1)[1].split("## Figures", 1)[0]
    out = {}
    for m in re.finditer(r"^(\d+)\. (.+)$", sec, flags=re.M):
        out[int(m.group(1))] = m.group(2).strip()
    return out


def compare_old(r, old):
    issues = []
    if r["volume"] and not re.search(rf"\b{re.escape(r['volume'])}\b", old):
        issues.append(f"volume {r['volume']} not in manuscript string")
    pg = r.get("pages") or r.get("article_number")
    if pg:
        first = re.split(r"[-\u2013]", pg)[0]
        if first not in old:
            issues.append(f"page/article {pg} not in manuscript string")
    if r["year"] and f"({r['year']})" not in old:
        issues.append(f"year {r['year']} (Crossref {r['year_source']}) differs from manuscript string")
    fam = (r["authors"][0]["family"] if r["authors"] else "") or ""
    if fam and fam.split()[0].lower() not in old.lower():
        issues.append(f"first author {fam} not in manuscript string")
    if r["doi"].lower() not in old.lower():
        issues.append("DOI string differs")
    return issues


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


def main(mode):
    recs = []
    for e in ENTRIES:
        if e["key"] in DROPPED:
            continue
        print("..", e["key"], file=sys.stderr)
        recs.append(record(e))
    if mode == "fetch":
        lines = ["# Abstracts retrieved at source (Worker C review file)\n"]
        for r in recs:
            ab = r["abstract_epmc"] or r["abstract_crossref"] or r["abstract_openalex"]
            src = "EuropePMC" if r["abstract_epmc"] else "Crossref" if r["abstract_crossref"] else "OpenAlex" if r["abstract_openalex"] else "NONE"
            lines.append(f"## {r['key']}  ({r['crossref_type']}; {r['journal']} {r['year']}; PMID {r['pmid']}; abstract: {src})\n")
            lines.append(f"Title: {r['title']}\n\n{ab}\n")
            for mt in r.get("manual_texts", []):
                lines.append(f"[{mt['kind']}; {mt['source']}]\n\n{mt['text']}\n")
        (OUT / "C_abstracts_review.md").write_text("\n".join(lines))
        print(f"fetched {len(recs)}", file=sys.stderr)
        return
    build(recs)


def build(recs):
    old_txt = old_reference_texts()
    audit, bib, nat, problems = [], [], [], []
    for r in recs:
        qs = QUOTES.get(r["key"]) or []
        qs = [qs] if isinstance(qs, str) else list(qs)
        texts = {"abstract (Europe PMC)": r["abstract_epmc"], "abstract (Crossref)": r["abstract_crossref"],
                 "abstract (OpenAlex)": r["abstract_openalex"]}
        for mt in r.get("manual_texts", []):
            texts[f"{mt['kind']} <{mt['source']}>"] = mt["text"]
        found = []
        for q in qs:
            hit = next((src for src, ab in texts.items() if ab and norm(q) in norm(ab)), None)
            found.append(hit)
        q = qs[0] if qs else None
        q_src = found[0] if found else None
        nwords = max((len(x.split()) for x in qs), default=0)
        all_found = bool(qs) and all(found)
        status, retracted = retraction_status(r)
        checks = {
            "doi_resolves_doi_org": r["doi_resolves"] if r["doi"] else None,
            "crossref_record_found": (r["crossref_status"] == 200) if r["doi"] else None,
            "crossref_type_peer_reviewed_journal": r["crossref_type"] in ("journal-article", "proceedings-article"),
            "title_matches_europepmc": (None if not r["epmc_title"] else
                                        norm(re.sub(r"\W", "", r["title"].lower())) == norm(re.sub(r"\W", "", r["epmc_title"].lower()))),
            "volume_matches_europepmc": None if not r["epmc_volume"] else str(r["epmc_volume"]) == str(r["volume"]),
            "year_matches_europepmc": None if not r["epmc_year"] else str(r["epmc_year"]) == str(r["year"]),
            "author_count_matches_europepmc": None if not r["epmc_n_authors"] else r["epmc_n_authors"] == r["n_authors_crossref"],
            "quote_found_in_abstract": all_found,
            "quote_le_40_words": 0 < nwords <= 40,
            "not_retracted": not retracted,
        }
        old_issues = compare_old(r, old_txt[r["old_no"]]) if r["origin"] == "old" else []
        a = {
            "key": r["key"], "origin": r["origin"], "manuscript_ref_no": r["old_no"],
            "doi": r["doi"], "non_doi_identifier": r.get("nondoi_id"), "url": r.get("url"), "pmid": r["pmid"], "pmcid": r["pmcid"],
            "verified_fields": {k: r[k] for k in ("title", "journal", "year", "year_source", "year_online", "year_print",
                                                  "volume", "issue", "pages", "pages_source", "article_number", "publisher", "crossref_type")}
                               | {"authors": [f"{x['family']}, {x['given']}".strip(", ") for x in r["authors"]],
                                  "n_authors": r["n_authors_crossref"]},
            "europepmc_fields": {"title": r["epmc_title"], "journal": r["epmc_journal"], "year": r["epmc_year"],
                                 "volume": r["epmc_volume"], "issue": r["epmc_issue"], "pages": r["epmc_pages"],
                                 "n_authors": r["epmc_n_authors"], "pubtypes": r["epmc_pubtypes"], "source": r["epmc_source"]},
            "crossref_match": checks,
            "retraction_status": status,
            "crossref_update_notices": r["crossref_update_notices"],
            "europepmc_comment_corrections": r["epmc_comment_corrections"],
            "supporting_quote": q, "quote_source": q_src, "quote_words": nwords,
            "additional_quotes": [{"quote": x, "source": f} for x, f in zip(qs[1:], found[1:])],
            "supported_claim": r["supported_claim"], "section": r["section"], "topic_group": r["topic"],
            "manuscript_string_issues": old_issues,
            "manuscript_string": old_txt.get(r["old_no"]) if r["origin"] == "old" else None,
            "note": NOTES.get(r["key"], ""),
            "nature_reference": nature_line(r),
        }
        peer = checks["crossref_type_peer_reviewed_journal"]
        hard = [k for k in ("doi_resolves_doi_org", "crossref_record_found", "quote_found_in_abstract",
                            "quote_le_40_words", "not_retracted") if checks[k] is False]
        soft = [k for k in ("title_matches_europepmc", "volume_matches_europepmc", "year_matches_europepmc",
                            "author_count_matches_europepmc") if checks[k] is False]
        biblio_hard = [k for k in hard if k not in ("quote_found_in_abstract", "quote_le_40_words")]
        if not hard:
            st = "VERIFIED"
        elif not biblio_hard and not qs:
            st = "BIBLIOGRAPHIC FIELDS VERIFIED; supporting text NOT VERIFIED (no abstract retrievable at source)"
        else:
            st = "FAILED: " + ", ".join(hard)
        a["verification_status"] = st + ("" if peer else " (NOT a peer-reviewed journal article: " + str(r["crossref_type"]) + ")")
        a["discrepancies_vs_europepmc"] = soft
        a["counts_as_peer_reviewed"] = peer and not hard
        a["bibliographic_only"] = peer and not biblio_hard and bool(hard)
        if hard or soft or old_issues:
            problems.append((r["key"], hard, soft, old_issues))
        audit.append(a)
        bib.append(bib_entry(r))
        nat.append(a["nature_reference"])
    (OUT / "C_reference_audit.json").write_text(json.dumps(audit, indent=1, ensure_ascii=False))
    (OUT / "C_references.bib").write_text(
        "% Paper 1 verified bibliography (Worker C, 15 October 2026). Fields from Crossref, cross-checked against\n"
        "% Europe PMC; see C_reference_audit.json. Keys AuthorYEARkeyword. NASA OSDR datasets are in C_datasets.bib.\n\n"
        + "\n".join(bib))
    (OUT / "C_references_nature.txt").write_text(
        "\n".join(f"{i}. {l}" for i, l in enumerate(nat, 1)) + "\n")
    json.dump([{"key": k, "hard": h, "soft": s, "old": o} for k, h, s, o in problems],
              open(OUT / "cache" / "problems.json", "w"), indent=1)
    for k, h, s, o in problems:
        print(f"PROBLEM {k}: hard={h} soft={s} old={o}")
    print(f"entries={len(audit)} peer_reviewed_verified={sum(a['counts_as_peer_reviewed'] for a in audit)}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fetch")
