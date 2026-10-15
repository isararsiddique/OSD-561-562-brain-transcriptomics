# H_REPORT — Systematic novelty search for Paper 1 (searched up to 15 October 2026)

## Verdicts

- **C1 Empirical: NOVEL.** No record measures how large a collection- or euthanasia-protocol effect is on brain
  immediate-early or rapid primary-response mRNA in spaceflight rodents or their never-flown controls. No record
  compares such an effect with the matched flight effect either. The closest record is Grandke 2026
  (Nat Commun, RRRM-1). It used small RNA from 13 tissues and both of its arms used ketamine/xylazine; it reported a
  pooled TERM-versus-LAR variance share but no size specific to brain: "However, the extraction-tissue interaction
  and extraction itself contributed 2.5% and 0.8% of the variance. This underscores the impact of the extraction
  timepoint (TERM: day 21, LAR: day 43) and the preservation method (see Methods) on miRNA expression."
  Keep the claim worded as brain **mRNA** / rapid primary-response genes. Do not say "brain" without that
  qualifier, because Grandke's pooled small-RNA figure includes brain.
- **C2 Conceptual: NOVEL.** No omics correction paper reports or proposes the spectral norm, operator gain or
  amplification factor of the correction operator as a routine diagnostic. No such paper separates near-perfect
  detectability from loading imprecision either. The closest records are adjacent, and the paper should cite the
  first two:
  - Li et al., two-step batch correction (bioRxiv 10.1101/2021.01.24.428009). It works at the level of the
    operator, but the operator is orthogonal-type with no gain above 1, and the paper proposes no diagnostic:
    "the samples in adjusted data are correlated, with the correlation matrix defined by M = (I − H12)(I − H12)T".
  - Batch Probing Score (bioRxiv 10.1101/2025.05.12.653389). It treats detectability as a sign that correction
    failed: "the most popular batch correctors fail to fully remove batch signal. In most cases, the batch of origin
    remains clearly identifiable after correction, even though standard evaluation metrics cannot detect…".
  - RiTex 2026 (radiomics harmonisation, PMC13302101). It uses an operator norm only inside a proof:
    "‖dT‖op the operator norm of the Fréchet derivative of T".
- **C3 Methodological: NOVEL**, for the combination as claimed. No record uses a mission's own flight-naive animals,
  collected under both protocols, as ground-only holdouts. No record fits an omics correction under a hard
  spectral-norm cap certified by a Frank–Wolfe duality gap. The closest spaceflight record is the GeneLab
  batch-correction evaluation (Front Astron Space Sci 2023, 10.3389/fspas.2023.1200132). It tabulates Animal Return
  (LAR versus ISS-T) and preservation as batch variables across seven liver datasets, but it has no holdout and no
  gain bound. The components remain established prior art, as the paper already states.
- **C4 Sample: NARROWED (provisional).** All seven deposits have an empty OSDR publication field. No Europe PMC,
  OpenAlex or Crossref record cites any of the 14 OSD or GLDS accessions. One closed-access record could not be
  ruled out: Piekarski, Guarnieri, Arif, Aslinger & McNerney 2026, "Spaceflight modulates
  neurodegeneration-associated transcriptomic pathways across brain cell types", Life Sci Space Res,
  doi:10.1016/j.lssr.2026.08.007.
  - Its indexed keywords are "single-nucleus transcriptomics", "brain cell types" and "age effects". Its reference
    list includes Masarapu 2024.
  - It may therefore use the single-nucleus component of OSD-612 (RR-10) or OSD-613 (RRRM-2), or it may use RR-3
    data. Neither the abstract nor the full text is available from any source searched.
  - No quote is possible.
  - **Surviving wording:** "The bulk RNA-seq libraries of the seven core deposits (OSD-561, -562, -563, -564, -612,
    -613, -525; GLDS-556, -557, -568, -569, -588, -589, -525) had no linked publication in OSDR, and no paper or
    preprint citing any of these accessions was found in Europe PMC, OpenAlex or Crossref on 15 October 2026."
  - Before submission, obtain the Piekarski full text. If it analyses OSD-612 or OSD-613, either remove those two
    from the claim or restrict the claim to their bulk libraries and cite the paper.

No record ANTICIPATES any claim.

## PRISMA-style counts

| Stage | n |
|---|---|
| Records identified by the registered strings | 7,557 (Europe PMC 2,498; OpenAlex 2,379; Crossref 2,680) |
| NASA OSDR study records checked | 7 (0 linked publications) |
| Records identified by stage-4 follow-up searches | 3 (one of them new) |
| Unique records after deduplication (DOI, then PMID, then PMCID, then title+year) | 6,249 (+1 follow-up = 6,250) |
| Screened on title/abstract | 6,250 |
| – excluded by the registered stage-3a rule | 3,711 |
| – stage-3c audit of 3a exclusions | 120 checked by hand, 0 false exclusions |
| – stage-3b manual title screen, abstract read where possibly relevant | 2,411 passing 3a, plus 7 rescued by title term |
| Full text assessed | 14 (13 read in full, 1 closed-access record assessed on metadata only) |
| Relevant, not anticipating | 66 |
| Unresolved (C4) | 1 |
| Decisive (ANTICIPATES) | 0 |

API calls: 160 of the 600 cap; 157 succeeded and are cached, and 3 Europe PMC fullTextXML requests returned
HTTP 500. Calls were at least 0.35 s apart, and all were read-only GETs.

Full texts read: Grandke 2026, Masarapu 2024 (RR-3), Kremsky 2023 (RR-18), Lai Polo 2020, Overbey 2022,
Earth-reference-controls 2020, Ganti 2026 preprint (OSDR OXPHOS census; liver, muscle and retina only; none of the
seven), GeneLab batch correction 2023, two-step correction 2021, Batch Probing Score 2025, RiTex 2026, deconfounded
regression 2025, and the ferroptosis paper in npj Microgravity 2026 (OSD-682/685/698/699; not one of the seven).
For the Piekarski 2026 paper only the metadata and reference list could be read.

## Positive-control recall

Twelve of twelve controls were retrieved by at least one string of their claim set, so no string revision was
needed.

| Control | Retrieved by (database:query) |
|---|---|
| Lai Polo 2020 | EPMC C1-2, C3-3; OpenAlex OA-C1-2; Crossref CR-C1-2 |
| Grandke 2026 | EPMC C1-3, C4-2; OpenAlex OA-C1-3 |
| Overbey 2022 | EPMC C1-2, C4-2; OpenAlex OA-C1-2 |
| Choi 2020 | EPMC C1-8 only (methods-restricted string). This is the weakest recall. |
| Kremsky 2023 | EPMC C1-3, C1-4, C1-8; OpenAlex; Crossref |
| Masarapu 2024 | EPMC C1-3, C3-3; OpenAlex |
| Holley 2022 | EPMC C1-3, C3-3; OpenAlex; Crossref |
| Staib-Lasarzik 2014 | EPMC C1-5, C1-6; OpenAlex; Crossref |
| Nygaard 2016 | EPMC C2-4; OpenAlex; Crossref |
| Zindler 2020 | EPMC C2-4; OpenAlex |
| Gerard 2021 | OpenAlex OA-C2-5; Crossref (not indexed in Europe PMC) |
| Wang 2017 | EPMC C2-5; OpenAlex; Crossref |

Recall was measured as DOI membership in the retrieved sets. No Europe PMC query was truncated, so this equals the
registered DOI-ANDed test for Europe PMC. For the truncated OpenAlex and Crossref sets it is a lower bound.

**External recall check** against B_closest_prior_art.json, opened only after screening was complete: 22 of its 32
items were retrieved. These 10 were not:
- the oblique-projection mathematics papers (10.1007/s11075-006-9046-2, 10.1109/78.286957);
- the Varma cross-validation bias paper and the ML leakage paper;
- Zames 1966;
- Miao 2018 on proxy confounders;
- dissociation and activity-gene artefact papers (10.1038/nmeth.4437, 10.1038/s41593-022-01022-8,
  10.1016/j.neuron.2018.04.001);
- the ketamine/isoflurane c-Fos mapping paper (eLife 88420);
- the cortex/hippocampus rebound-pattern paper (10.3389/fnmol.2022.892345);
- the pig euthanasia study (10.1186/s12864-023-09794-4).

All 10 are component or non-spaceflight context: mathematics, machine learning, dissociation, pig. The strings did
not target them, and none bears on a verdict.

This search found five records that B lacked:
- Piekarski 2026, the unresolved C4 record;
- the Batch Probing Score and two-step correction papers (C2);
- the GeneLab batch-correction evaluation (C3);
- the Ganti 2026 OSDR census.

## Deviations from the protocol

1. The brief's Python environment (runtime-735c0e1b/p1env) disappeared mid-run, so system python3 3.14 was used.
   Nothing was installed.
2. The Kremsky and Gerard control DOIs were taken from src/paper1.bib, because they are absent from
   refs/C_refs/C_references.bib.
3. Recall was checked by DOI membership instead of re-querying with DOI-ANDed strings (see above).
4. OpenAlex parsed some boolean strings loosely: OA-C2-5 returned 95,891 hits and OA-C1-5 returned 3,809. In those
   cases only the registered top 200 by relevance were screened. Crossref has no boolean operators and was always
   relevance-truncated (top 100, or top 20 for accessions).
5. Stage 3b screened titles for every record and read abstracts only where a title was possibly relevant. Not every
   abstract was read.
6. Stage-4 follow-ups outside the registered strings, all logged, were:
   - a Europe PMC author search (Piekarski / Guarnieri);
   - the OpenAlex reference list of the LSSR paper;
   - one Semantic Scholar record lookup (abstract null). Semantic Scholar is not a protocol database and was used
     for this single record only;
   - an OpenAlex search for spaceflight single-nucleus brain work;
   - a title lookup for the npj Microgravity ferroptosis paper.

## Limits

- Not searched: Scopus, Web of Science, Google Scholar, Embase and Dimensions. They require a login, and scraping was
  not allowed.
- OpenAlex full-text coverage is partial, and Crossref matching is relevance-only.
- Indexing lags for 2026 work. The Aug 2026 ferroptosis paper and the LSSR paper are absent from Europe PMC, and the
  LSSR paper is closed access. Papers from the last 6–8 weeks may be missing.
- Accession matching works only when an accession appears in indexed text. A paper that reuses the data without
  naming OSD or GLDS numbers would be missed. The C4-2 to C4-4 mission strings partly offset this.
- Data observation, not a novelty issue: OSDR metadata for OSD-561, -562 and -613 lists Mission SpaceX-18
  (07/25/2019–08/27/2019) and Project RR-17. This looks inconsistent with an RRRM-2 flight, so check it before
  quoting mission dates.
- Side note for C1 interpretation: in Grandke (RRRM-1), both LAR and TERM used ketamine/xylazine, with acepromazine
  added for TERM.

## Files

- `H_PROTOCOL.md` and `H_PROTOCOL.sha256`: hashed at 2026-10-15T05:00:07Z, before any search; the hash still
  verifies.
- `H_search_log.json`: 122 logged queries with string, filters, time, hit count and record IDs. Every HTTP request
  is also listed in `cache/raw_index.tsv`.
- `H_screening.csv`: 6,250 rows.
- `cache/`: raw responses and the records, unique and triaged sets.
- `tools/`: client, queries, deduplication, triage, full-text grep and output builder.
