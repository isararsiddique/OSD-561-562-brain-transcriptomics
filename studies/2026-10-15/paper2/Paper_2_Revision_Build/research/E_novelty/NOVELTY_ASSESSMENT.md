# Paper 2: novelty assessment (15 Oct 2026)

Evidence files in this folder: `epmc_novelty_search2.json` (title/abstract queries A–N), `epmc_fulltext_accessions.json`,
`epmc_fulltext_cerebellum_flight_rnaseq.json`, `epmc_method_prior_art.json`, `osdr_meta_{525,561,563}.json`,
`crossref_new_refs.json`. These were targeted searches of Europe PMC (title/abstract plus open full text) and the OSDR
publication links. They are not a systematic review, and unpublished or non-indexed work may exist.

| Axis | Verdict | Evidence | Wording placed in the paper |
|---|---|---|---|
| Population / sample | **ESTABLISHED (bounded)** | OSDR metadata lists 0 publications for OSD-525, -561 and -563. A Europe PMC full-text search for the OSD-/GLDS- accession strings returned 0 records. The searches found no cerebellar RNA-seq or junction study in spaceflight mice. The closest brain precedent (Holley 2022, npj Microgravity) used NanoString on whole brain from male C57BL/6 mice. Two of our cohorts are female, and RR10 is B6129SF2/J. | "no linked publication ... on 15 October 2026"; "extending beyond the male C57BL/6 animals of the closest brain-expression precedent" |
| Empirical | **ESTABLISHED (bounded)** | Prior spaceflight and unloading splicing work covers muscle (Henrich 2022; Sun 2021, hindlimb unloading) and Arabidopsis (Beisel 2019). Query A returned no brain splicing hit. Housing differences between AEM/habitat and vivarium were already reported (Beheshti 2018, transcriptome in other tissues; Hong 2021, ovary), so the GC–VIV difference is an **extension to brain RNA processing**, not a new phenomenon. Its size is new: binary GC–VIV 4.0746 pp = 73% of FLT–GC (5.5442) and above FLT–VIV (3.8144); multialternative GC–VIV 0.06551 vs FLT–VIV 0.06615. | "habitat and vivarium controls differ in cerebellar RNA processing by a margin comparable to flight, extending housing-associated differences reported in other tissues" |
| Methodological | **PARTIAL: design contribution, not new statistical theory** | Discrete-p FDR is an active field (Gilbert 2005; Heyse-type and Chen 2018/2020 procedures; Döhler 2018; Jiang 2017 DS-FDR). Targeted searches found no published closed-form tie-block attainability criterion k* for BY under exact permutation. However, k* follows in two lines from known facts, so it is framed as a design and reporting check. It is combined with complementary read partitions and simulation (Fig. 9). | "builds on, rather than replaces, discrete FDR procedures" |
| Theoretical | **REFINEMENT only** | The paper refines one interpretive assumption: that FLT−GC is a clean flight effect with GC as a neutral reference. It also separates effect repeatability from decision repeatability. No new theory of spaceflight biology is claimed. | Folded into the conceptual item: "refines the common practice of treating the habitat control as a neutral reference" |
| Conceptual | **ESTABLISHED (bounded)** | Query E (transportability/external validation × spaceflight) returned no space-omics paper using frozen-direction transport. The paper treats identifiability (exact rank deficiency), attainability (k*) and transport (comparator-symmetric frozen test) as three separate evidence layers. | "integrated into one evidence standard for space-omics reanalysis" |

## Not claimed, and why
- "First" anything: the search is bounded and the deposits are public.
- A causal microgravity effect: flight and collection date are aliased (rank 3 of 4).
- A new algorithm class: the DM score and the composition tests are established mathematics.

## Where it appears in the manuscript
Abstract (sample and housing sentences); Introduction (four-item novelty list); Results (73% comparison, with a
dependence caveat); Discussion (paragraph comparing against prior work, with the limits of the search); Methods (the
15 Oct novelty-search record). Six new references were verified against Crossref (`crossref_new_refs.json`).
