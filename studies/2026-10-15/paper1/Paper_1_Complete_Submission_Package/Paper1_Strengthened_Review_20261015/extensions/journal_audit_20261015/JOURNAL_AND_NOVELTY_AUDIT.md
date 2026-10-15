# npj Microgravity format and novelty audit

Verified 15 October 2026, Asia/Kolkata. This audit reads the existing manuscript and source evidence; it does not independently rerun biological analyses, certify acceptance, submit the paper, or contact anyone. Publisher HTML, extracted text and hashes accompany this report. Earlier manuscript counts are a timestamped snapshot because the parent agent is editing concurrently.

## Exact Article format

The authoritative [Content types](https://www.nature.com/npjmgrav/content-types) page states: “We do not impose strict limits on word count or page numbers.” There is **no fixed Article main-text word target**. A main narrative of approximately 4,400 words does not breach a publisher limit. Concision is recommended.

| Element | Exact publisher instruction |
| --- | --- |
| Title | Up to 15 words; free of punctuation, idioms and puns. |
| Abstract | Unstructured; no subheadings; up to 150 words. |
| Introduction | No subheadings. |
| Results | Use subheadings. |
| Discussion | No subheadings; no separate limitations or conclusions sections. Limitations may be explained in connected Discussion prose. |
| Methods | Use subheadings. **All Methods must be in the main manuscript file; Supplementary Methods are not permitted.** Put mathematical derivations needed to reproduce the method into a main Methods subsection. |
| Data availability | Mandatory. |
| Code availability | Required where applicable. The current general reporting policy specifies a separate section after Data availability and before References; the Article content-type table also lists it separately. |
| Acknowledgments | Declare funding here. **A separate Funding statement is not permitted.** |
| Author contributions | Mandatory; identify all authors individually by initials. |
| Competing interests | Mandatory; disclose financial and non-financial interests for each author, or state none. |
| References | Guide of 60, not strictly enforced; numbered in order of first citation; Nature style. |
| Figure legends | Up to 350 words for each complete figure. |

The Article table lists Introduction–Results–Discussion–Methods. The [journal homepage](https://www.nature.com/npjmgrav/) explicitly permits either that order **or IMRaD**, with Methods before Results. There is no need to force a different Methods position when the chosen order is consistent.

The [Submission guidelines](https://www.nature.com/npjmgrav/for-authors-and-referees/submission-guidelines) allow flexible initial formatting, but the user's request is to prepare the paper properly; apply the final requirements now. A combined editable Word file and PDF are suitable for peer review. Main figure legends go after References; put each legend on the same page as its figure in the combined review file. Tables belong at the end of the text document. A cover letter is required. Related overlapping manuscripts under consideration or in press must be disclosed by the authors.

## Artwork and statistical reporting

The official guidance requires RGB images at **at least 300 dpi**, consistent **Arial or Helvetica**, an optimum **8 pt** font at intended output size, and thinnest final lines **at least 1 pt**. Use a white background, accessible colors, comparable color visibility, and avoid red–green contrasts or rainbow scales. Use lowercase bold panel labels. Submit each complete multipanel figure in **one file on one page**, with logically connected panels. Editable vector files are preferred; PDF/SVG and high-resolution PNG are accepted. The page says to prepare for A4 page width rather than imposing a universal numerical figure width. It does not specify a fixed maximum number of main figures; unnecessary display items should be avoided.

The statistical instructions specifically require plots to show the full data distribution. Bar charts need overlaid individual points or conversion to a dot plot or box-and-whisker plot. For simulation figures, show repetition distributions or point/distribution summaries, and define Monte Carlo mean uncertainty separately from real-animal uncertainty. Split distributions must never imply independent validation samples.

The main Methods must include a **Statistics and reproducibility** subsection. Report each test's comparison, biological n, one- or two-sidedness, alpha, exact probability, multiplicity family and rationale; justify approximations for small samples. Define every error bar and center/spread statistic. The study's library, animal, animal–region, paired-animal, overlapping-split and independent-simulation counts have different meanings and must remain distinct in captions.

Supplementary Information should be one merged PDF without Methods. Large spreadsheet tables should be separately named **Supplementary Data XX**, with titles and legends in the PDF, rather than labeled Supplementary Tables. All cited supplementary content must be available for reviewers. The Reporting Summary is required during revision and published with accepted research; the editorial policy checklist is required before review. Both are encouraged at initial submission. ARRIVE Essential 10 is strongly encouraged for preclinical animal work; secondary public-data analyses should report source ethics faithfully and clearly identify which original animal procedures were outside the reanalysis team's control. Corresponding-author ORCID must be linked before final submission.

**Graphical abstract:** no graphical abstract is required in the inspected Article/content-type/submission guidance. A scientific overview schematic may be retained when it helps communicate the actual analysis. A suggested homepage featured image is a distinct, optional final-stage item, for which the Image Licence to Publish form is requested. Avoid portraying a hypothesized causal pathway as experimentally established. Current Nature Portfolio [AI policy](https://www.nature.com/nature-portfolio/editorial-policies/ai) permits disclosed, verified assistance under human oversight and does not permit fabricated results or undisclosed AI replacement of scholarly judgment. The existing substantive Codex disclosure is appropriate; authors still need to check the analyses and interpretation.

## Strongest defensible novelty

The strongest framing is the **empirical separation of repeatable collection-cohort scores from transferable genome-wide correction**, supported by animal-aware inference, training-only preprocessing, correction geometry, uncertainty propagation and explicit failure benchmarks. It is stronger and better differentiated than a universal neuroinflammation or molecular collection-clock claim.

The reported result is concrete: core ground collection discrimination remains perfect over the prespecified balanced holdouts, while original outside-sensor predictions worsen relative to the nuisance-only baseline in cerebellum and RR-23. The three settings use only 30 unique source ground animals and overlapping splits. A repeatable signature is therefore insufficient evidence that its fitted genome-wide loading identifies transferable nuisance. This is a within-resource counterexample, not external mission validation.

The spectral-gain result gives a mathematical explanation and a transparent guardrail. The Euclidean closest-loading construction is standard constrained projection geometry. Present its application and failure analysis as the contribution, rather than claiming invention of oblique projection, universal superiority over SVA/RUV, biological FDR control, or recovery of the true flight response. The selected cap trades reduced instability for incomplete nuisance removal; the retained high-gain failures support the manuscript's credibility.

| Closest source | What it already establishes and how the present paper differs |
| --- | --- |
| [Holley et al. 2022](https://doi.org/10.1038/s41526-022-00217-4) | RR-9 male C57BL/6 mice, 35-day flight, right caudal half hemisphere, NanoString 757-gene neuroinflammation panel, n=5–6 per group. CO2 euthanasia was about 38 ± 4 h after splashdown. The Discussion explicitly recognizes acute landing/collection effects and proposes orbit euthanasia with immediate preservation. It neither analyses these RRRM-2/RR-10/RR-23 counts nor tests ground-loading transfer. Do not claim to be first to recognize recovery/collection ambiguity or to refute this different assay. |
| [Dachet et al. 2021](https://doi.org/10.1038/s41598-021-85801-6) | Human surgically removed cortex followed over 0–24 h already showed preferential neuronal activity-transcript loss and glial-expression increases despite relatively stable housekeeping genes/RIN. This establishes tissue vulnerability; it does not calibrate a mouse spaceflight collection clock. |
| [Marsh et al. 2022](https://doi.org/10.1038/s41593-022-01022-8) | Ex vivo brain dissociation creates microglial signatures; a preservation protocol was experimentally validated. Those cell-resolved experiments do not establish that every bulk glial signal here is an artifact. |
| [Overbey et al. 2022](https://doi.org/10.1016/j.crmeth.2022.100325) | Spaceflight single-cell/spatial collection constraints were already identified. The present contribution quantifies correction stability in an observational bulk resource. |
| [Nygaard et al. 2016](https://doi.org/10.1093/biostatistics/kxv027) and existing SVA/RUV references | Batch adjustment bias and exaggerated downstream confidence, and control-based nuisance estimation, are established. This paper offers a spaceflight application with a measurable loading-gain safeguard, not a first general batch-confounding theorem. |
| [Gerstner et al. 2016](https://doi.org/10.1186/s12864-016-3065-8) | Previously overlooked close context: RUV already integrated mouse cortex/hippocampus/whole-brain sleep datasets, with time-of-day-matched controls and sleep-recovery dynamics. Clock-related expression can reflect sleep, stress and sampling time; clock-gene association is not uniquely a microgravity mechanism. This paper does not validate the present bounded operator. |
| [Masarapu et al. 2024](https://doi.org/10.1038/s41467-024-48916-8) | RR-3 spatial/single-nucleus multiomics and bulk comparisons share source animals. They provide anatomical context, not independent validation of this correction or of these cohorts. |
| [Jimenez-Chavez et al. 2026](https://doi.org/10.1038/s41526-026-00639-4) | Newly published 12 October 2026: PBMCs from ten astronauts, with age/mission-duration associations and more pronounced responses in younger astronauts. Relevant current context only; the publicly inspected page gives the abstract, not full main text. It provides no brain or collection-correction validation. |

The prior extension's verified Grandke 2026 small-RNA and Roy 2025 protein results remain useful context. Roy reports decreased inflammasome activation, so it must not be cited as demonstrating increased neuroinflammation. These assays and cohorts should not be pooled as equivalent validation.

## Claims and remaining scientific limits

Retain the rapid-response ground association as observational. Keep the circadian flight-by-collection interaction as a candidate, with its full family correction, probability/Monte Carlo uncertainty, count-model comparison and altered age support after exclusions. Strong hemoglobin interactions are compatible with composition or collection-linked biology, but do not identify blood fraction, contamination or a cell source. Paired tests found no corrected regional panel effect, and the cross-cohort screen found no corrected reproducible gene; a nominal ranked list cannot become a validated signature.

No search establishes exhaustive priority. The targeted Europe PMC endpoint returned HTTP 503, and a Crossref search was broad, partially rate-limited and poorly specific. Exact primary-source comparisons, the earlier cached literature review and the current publisher check support differentiated framing, but not “first ever” absence claims. The full source records preserve these failures.

Before actual submission, authors still need to confirm their contributions and coauthor consent, final journal fit, source ethics, complete review access to the code/data, and a real versioned archive identifier. Do not fabricate a repository or DOI. Journal ranking or more significant results cannot guarantee Q1 acceptance.

## Saved evidence

- `source_provenance.json`: official/current sources, URLs, retrieval status and hashes.
- `content_types.txt`, `submission_guidelines.txt`, `journal_home.txt`, `reporting_policy.txt`: exact publisher requirements.
- `holley2022.txt`, `sleep_ruv2016.txt`, `jimenez2026.txt`: newly retrieved primary publisher evidence.
- `prior_primary_text_evidence.json`: hashes, titles, DOI and abstracts from cached Dachet/Marsh/Overbey/Nygaard full text.
- `manuscript_format_snapshot.json`: timestamped counts and former Supplementary Methods references; not a final-format certification.
- `suggested_abstract_150_words.txt`: 138-word optional abstract draft, to be updated if the ongoing statistical audit changes a reported result.
- `targeted_search_provenance.json`: exact additional searches, including failure statuses.

