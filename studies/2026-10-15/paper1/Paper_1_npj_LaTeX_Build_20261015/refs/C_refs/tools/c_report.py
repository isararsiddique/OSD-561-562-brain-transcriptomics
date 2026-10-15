#!/usr/bin/env python3
"""Worker C: write ../C_REPORT.md from C_reference_audit.json and cache/datasets_audit.json."""
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from c_quotes import CANDIDATES_NOT_USED  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent.parent
A = json.loads((OUT / "C_reference_audit.json").read_text())
D = json.loads((OUT / "cache" / "datasets_audit.json").read_text())


def short(a):
    au = a["verified_fields"]["authors"]
    fam = au[0].split(",")[0] if au else "?"
    return f"{fam}{' et al.' if len(au) > 2 else (' & ' + au[1].split(',')[0] if len(au) == 2 else '')} {a['verified_fields']['year']}"


def status(a):
    s = a["verification_status"]
    if s.startswith("VERIFIED") and a["counts_as_peer_reviewed"]:
        tag = "VERIFIED"
    elif a["bibliographic_only"]:
        tag = "BIBLIO ONLY (no abstract at source)"
    elif "monograph" in s or "book" in s:
        tag = "VERIFIED (book, not journal)"
    else:
        tag = s
    if "correction(s) recorded" in a["retraction_status"]:
        tag += "; correction exists"
    if a["discrepancies_vs_europepmc"]:
        tag += "; EPMC differs: " + ",".join(x.replace("_matches_europepmc", "") for x in a["discrepancies_vs_europepmc"])
    return tag


old = [a for a in A if a["origin"] == "old"]
new = [a for a in A if a["origin"] == "new"]
full = [a for a in A if a["counts_as_peer_reviewed"]]
journal_full = [a for a in full if a["verified_fields"]["crossref_type"] == "journal-article"]
procs = [a for a in full if a["verified_fields"]["crossref_type"] == "proceedings-article"]
biblio = [a for a in A if a["bibliographic_only"]]
books = [a for a in A if a["verified_fields"]["crossref_type"] in ("monograph", "book")]
corr = [a for a in A if "correction(s) recorded" in a["retraction_status"]]
retr = [a for a in A if a["retraction_status"].startswith("RETRACTED")]
old_issues = [a for a in old if a["manuscript_string_issues"]]
multi_q = sum(1 for a in A for _ in a.get("additional_quotes", []))

L = []
w = L.append
w("# Worker C - verified bibliography for Paper 1 (15 October 2026)\n")
w("## Headline verdicts\n")
w(f"- **{len(full)} references fully verified and peer-reviewed** ({len(journal_full)} journal articles + {len(procs)} archival "
  f"proceedings paper without a DOI). Target of >= 48 met. Each has a resolving DOI (or publisher page), a Crossref record, "
  f"no retraction/expression of concern, and a <= 40-word quote that the script found verbatim in the abstract or cached full text.")
w(f"- Re-verified from MANUSCRIPT_COMPLETE.md: {len(old)} of 27 non-dataset references. "
  f"{sum(1 for a in old if a['counts_as_peer_reviewed'])} fully verified; Wold 2001 (ref 19) bibliographic fields only "
  f"(no abstract retrievable at source); Boyd & Vandenberghe 2004 (ref 20) is a Cambridge University Press book, verified "
  f"but not a peer-reviewed journal article. All 27 manuscript strings agree with Crossref on volume, first page/article "
  f"number, year, first author and DOI ({len(old_issues)} mismatches).")
w(f"- New: {len(new)} references added ({sum(1 for a in new if a['counts_as_peer_reviewed'])} fully verified; "
  f"Frank & Wolfe 1956 and Izenman 1975 bibliographic fields only).")
w(f"- Retractions / expressions of concern: **{len(retr)}**. Published corrections exist for {len(corr)} entries "
  f"({', '.join(a['key'] for a in corr)}); none changes the cited finding as far as the correction titles show (correction "
  f"texts not read).")
w("- **Datasets (refs 23-33): all 11 OSDR studies have DataCite DOIs that resolve and land on the correct OSD page. The "
  "manuscript gives a DOI for only 2 (OSD-525, OSD-612).** Five DOIs (OSD-561, 562, 563, 564, 613) are absent from the "
  "OSDR study metadata itself and were found in DataCite by exact title, then accepted only because DataCite's landing URL "
  "is that accession's OSDR page. npj Microgravity guidelines say a dataset citation should include a DOI or accession code "
  "where available, so add all 11.")
w("- **Manuscript sentences that the cited source does not fully support** (fix before submission):")
w("  1. Ref 1 (Afshinnekoo 2020): the abstract's six features are oxidative stress, DNA damage, mitochondrial dysregulation, "
  "epigenetic changes, telomere length and microbiome shifts. It does not list 'immune responses' or 'neuronal adaptation'. "
  "Reword, or co-cite da Silveira 2020 (innate immunity, chronic inflammation, circadian rhythm).")
w("  2. Ref 5 (Tullai 2007) is growth-factor stimulation of cultured cells, not neuronal activity. It supports the IEG / "
  "delayed-primary / secondary classification; keep it for that and let Tyssowski 2018 carry 'neuronal activity'.")
w("  3. Ref 10 (Masarapu 2024): 'limited concordance with bulk expression from the same mice' is in the full text, not the "
  "abstract (11 shared DEGs out of 825 and 629, only 2 with the same direction). Supported, but cite the full-text result.")
w("  4. Ref 22 (Johnson 2007) is cited for 'the official neuroCombat 0.2.12 implementation'. neuroCombat is the Fortin et "
  "al. implementation: add Fortin 2018 (NeuroImage) at that sentence.")
w("  5. Ref 34 (Winkler 2014) is cited for 'the Monte Carlo plus-one alternative'. The plus-one rule is Phipson & Smyth "
  "2010: add it there.")
w("- **Evidence bearing on the protocol-confound reframe (state it as found):**")
w("  - Lai Polo et al. 2020 (iScience) already showed that rodent spaceflight RNA-seq is confounded by sample collection: "
  "'Changes due to preservation method eclipsed those between flight and ground samples'. The general claim 'collection "
  "protocol confounds rodent spaceflight transcriptomes' is therefore **not new**; Paper 1's contribution has to be the "
  "brain-specific IEG signature, the exact aliasing in RRRM-2 (euthanasia agent and freezing method with cohort), and the "
  "post hoc correction-gain result. Cite Lai Polo 2020 prominently and position against it (relevant to Worker B).")
w("  - Staib-Lasarzik 2014: a few minutes of anaesthesia before sampling change hippocampal mRNA, mainly early-regulated "
  "genes; Hamaya 2000: anaesthetics transiently raise brain c-fos mRNA; Bunting 2016: high-dose isoflurane blocks "
  "initiation of Arc/Egr1 transcription. **None of these tested ketamine/xylazine**, so they support the principle, not the "
  "specific agent.")
w("  - Counter-evidence that must be cited with any ketamine argument: Nakao 1993 found **no** c-Fos protein in rat "
  "hippocampus after 100 mg/kg ketamine (c-Fos was present in neocortex and thalamus). Davoudian 2023 found c-Fos "
  "elevation (incl. CA1) after **subanaesthetic** ketamine. Both are protein read-outs on a timescale longer than an "
  "overdose-to-dissection interval; neither settles mRNA within minutes.")
w("  - Chakkingal Bhaskaran 2023 (pigs): euthanasia method did **not** change protein-coding expression; storage condition "
  "(RNAlater vs liquid nitrogen) did. Dachet 2021: activity-dependent transcripts fall within hours ex vivo. Together these "
  "point at freezing method / time-to-preservation as at least as plausible as the anaesthetic. Neither RRRM-2 factor can "
  "be separated from the other by the deposited design.")
w("- npj Microgravity submission guidelines (read 15 Oct 2026, cached) state no reference-count cap; they allow datasets in "
  "the numbered list. With 73 literature references plus 11 datasets the list would be 84 entries; the parent should cut "
  "to what the rewritten text actually cites.")
w("")
w("## Counts\n")
w("| quantity | n |\n|---|---|")
for k, v in [("non-dataset references in audit", len(A)), ("  re-verified from manuscript (refs 1-22, 34-38)", len(old)),
             ("  new", len(new)), ("fully verified, peer-reviewed (journal + proceedings)", len(full)),
             ("  journal articles", len(journal_full)), ("  archival proceedings, no DOI (PMLR)", len(procs)),
             ("bibliographic fields verified, supporting text NOT verified", len(biblio)),
             ("book (not counted as peer-reviewed)", len(books)),
             ("retracted / expression of concern", len(retr)), ("with a published correction", len(corr)),
             ("additional quotes checked beyond the primary quote", multi_q),
             ("indexed in Europe PMC with a PMID", sum(1 for a in A if a["pmid"])),
             ("  of which with a PMCID (free full text)", sum(1 for a in A if a["pmcid"])),
             ("Crossref type posted-content (preprint)", sum(1 for a in A if a["verified_fields"]["crossref_type"] == "posted-content")),
             ("OSDR dataset entries (separate list)", len(D)), ("  with a resolving dataset DOI", sum(1 for d in D if d["doi_resolves"])),
             ("  DOI already in manuscript", sum(1 for d in D if d["doi_in_manuscript"])),
             ("dropped from the manuscript's 27", 0)]:
    w(f"| {k} | {v} |")
w("\nNothing from the manuscript's 27 was dropped: all exist, resolve and match. Boyd 2004 stays as a book reference; Wold "
  "2001 stays (standard PLS reference) with its supporting text marked NOT VERIFIED.\n")
w("## Reference table\n")
w("Status VERIFIED = DOI resolves (doi.org handle API), Crossref record matches, Europe PMC cross-check (title, volume, year, "
  "author count) where indexed, no retraction/EoC in Crossref `updates` filter or Europe PMC, and the quote in "
  "`C_reference_audit.json` was found verbatim at source. Topic groups: a spaceflight; b collection/anaesthesia; c activity-"
  "regulated transcription; d unwanted variation; e validation/inference; f linear algebra/optimisation; g circadian/erythroid.\n")
w("| # | key | short citation | journal | year | origin | grp | status | supported claim (sentence / topic) |")
w("|---|---|---|---|---|---|---|---|---|")
for i, a in enumerate(A, 1):
    o = f"old ref {a['manuscript_ref_no']}" if a["origin"] == "old" else "new"
    claim = a["supported_claim"].replace("|", "/")
    w(f"| {i} | `{a['key']}` | {short(a)} | {a['verified_fields']['journal']} | {a['verified_fields']['year']} | {o} | "
      f"{a['topic_group']} | {status(a)} | {claim} |")
w("\n## Datasets (kept separately, `C_datasets.bib`)\n")
w("| ms ref | accession | DOI | DOI source | resolves | in manuscript | title matches manuscript |")
w("|---|---|---|---|---|---|---|")
for d in D:
    src = "OSDR study metadata" if d["doi_source"] == "OSDR live metadata" else ("DataCite (title + landing URL)" if d["doi_source"] else "none")
    w(f"| {d['manuscript_ref_no']} | {d['accession']} | {d['doi']} | {src} | {d['doi_resolves']} | {d['doi_in_manuscript']} | "
      f"{d['title_matches_manuscript']} |")
w("\nNature-style dataset lines (for the reference list):\n")
for d in D:
    w(f"- NASA Open Science Data Repository. {d['title_osdr_live']}. NASA OSDR {d['accession']} "
      f"https://doi.org/{d['doi']} ({d['datacite_year']}).")
w("\n## Candidates from the brief that were replaced, re-sourced or not used\n")
w("| candidate | outcome |\n|---|---|")
for c, why in CANDIDATES_NOT_USED:
    w(f"| {c} | {why} |")
w("\n## Per-entry caveats\n")
for a in A:
    if a["note"]:
        w(f"- `{a['key']}`: {a['note']}")
w("\n## Method\n")
w("- One script, `tools/c_verify.py` (fetch, then build), with all API responses cached as JSON under `cache/` "
  "(crossref, crossref_updates, doi_handle, europepmc, openalex, semanticscholar, datacite, osdr_meta). Rate: >= 1.2 s "
  "between Crossref calls, back-off on HTTP 429/5xx; plain User-Agent, no e-mail address sent.")
w("- Bibliographic fields are taken from Crossref, not typed; Europe PMC is the independent cross-check. Year = print-issue "
  "year where Crossref has one, else online year (Molania 2023, O'Donoghue 2015 have earlier online dates; Bunting 2016 is "
  "in the 2015 volume of Front. Behav. Neurosci.).")
w("- Quotes are in `tools/c_quotes.py`; the build asserts each is an exact substring of an abstract (Europe PMC, Crossref, "
  "OpenAlex) or of a cached publisher text built by `tools/c_texts.py` (van den Brink 2017 opening paragraph, Masarapu "
  "2024 full-text result, O'Donoghue 2015 Semantic Scholar abstract, Jaggi 2013 PMLR abstract). Szyld 2006's abstract was "
  "read with the web_fetch tool because Springer serves a bot wall to scripts; it is transcribed verbatim and flagged.")
w("- Retraction/correction: Crossref `works?filter=updates:<DOI>` (notices that update the work) plus Europe PMC "
  "`commentCorrectionList` and `pubTypeList`. Correction notices were identified, not read.")
w("- Datasets: `tools/c_datasets.py` reads live OSDR metadata, the deposited snapshot in the source package, DataCite and "
  "the doi.org handle API.")
w("- `C_references.bib` + `C_datasets.bib` compile with tectonic/BibTeX (unsrt, TeX Gyre Heros) with no missing glyphs "
  "(test in scratch, not kept).")
w("- Not done: the correction notices themselves were not read; abstracts were read but full texts were not (except "
  "Masarapu 2024); Wold 2001, Frank & Wolfe 1956 and Izenman 1975 supporting text could not be retrieved (ScienceDirect and "
  "Wiley return 403 to scripted and web_fetch requests). Unavailable MCP servers in this session: fusion (failed to "
  "start), cloudflare-api and github (awaiting authorization); none was needed.")
w("\n## Files (all in `refs/C_refs/`)\n")
for f, d in [("C_references.bib", "73 literature entries, keys AuthorYEARkeyword"),
             ("C_datasets.bib", "11 @misc OSDR entries with DOI, URL, howpublished, note 'Accessed 15 October 2026'"),
             ("C_reference_audit.json", "per entry: key, DOI, PMID, PMCID, Crossref fields, Europe PMC fields, match checks, "
                                        "retraction status, notices, quote(s) and where found, supported claim, section, notes"),
             ("C_references_nature.txt", "Nature-style list in the same order as the table"),
             ("C_abstracts_review.md", "every abstract retrieved, for re-reading"),
             ("tools/c_entries.py", "keys, DOIs, sections, claims"), ("tools/c_quotes.py", "quotes, abbreviations, notes, unused candidates"),
             ("tools/c_verify.py", "fetch + verify + emit"), ("tools/c_texts.py", "builds cache/manual_texts.json"),
             ("tools/c_datasets.py", "dataset verification"), ("tools/c_report.py", "this report"),
             ("tools/c_search.py", "discovery searches (Europe PMC / Crossref), cached under cache/search"),
             ("cache/", "raw API responses and publisher pages")]:
    w(f"- `{f}`: {d}")
w("\nRebuild: `cd refs/C_refs/tools && P=python && $P c_texts.py && "
  "$P c_verify.py build && $P c_datasets.py && $P c_report.py` (offline from cache).")
(OUT / "C_REPORT.md").write_text("\n".join(L) + "\n")
print("wrote", OUT / "C_REPORT.md", len(L), "lines")
