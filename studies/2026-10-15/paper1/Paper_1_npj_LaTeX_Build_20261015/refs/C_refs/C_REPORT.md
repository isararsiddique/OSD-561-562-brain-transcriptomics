# Worker C - verified bibliography for Paper 1 (15 October 2026)

## Headline verdicts

- **69 references fully verified and peer-reviewed** (68 journal articles + 1 archival proceedings paper without a DOI). Target of >= 48 met. Each has a resolving DOI (or publisher page), a Crossref record, no retraction/expression of concern, and a <= 40-word quote that the script found verbatim in the abstract or cached full text.
- Re-verified from MANUSCRIPT_COMPLETE.md: 27 of 27 non-dataset references. 25 fully verified; Wold 2001 (ref 19) bibliographic fields only (no abstract retrievable at source); Boyd & Vandenberghe 2004 (ref 20) is a Cambridge University Press book, verified but not a peer-reviewed journal article. All 27 manuscript strings agree with Crossref on volume, first page/article number, year, first author and DOI (0 mismatches).
- New: 46 references added (44 fully verified; Frank & Wolfe 1956 and Izenman 1975 bibliographic fields only).
- Retractions / expressions of concern: **0**. Published corrections exist for 6 entries (Afshinnekoo2020spaceflight, Marsh2022artifactual, Ronca2019behavior, Lacar2016nuclear, Hrvatin2018visual, Button2013power); none changes the cited finding as far as the correction titles show (correction texts not read).
- **Datasets (refs 23-33): all 11 OSDR studies have DataCite DOIs that resolve and land on the correct OSD page. The manuscript gives a DOI for only 2 (OSD-525, OSD-612).** Five DOIs (OSD-561, 562, 563, 564, 613) are absent from the OSDR study metadata itself and were found in DataCite by exact title, then accepted only because DataCite's landing URL is that accession's OSDR page. npj Microgravity guidelines say a dataset citation should include a DOI or accession code where available, so add all 11.
- **Manuscript sentences that the cited source does not fully support** (fix before submission):
  1. Ref 1 (Afshinnekoo 2020): the abstract's six features are oxidative stress, DNA damage, mitochondrial dysregulation, epigenetic changes, telomere length and microbiome shifts. It does not list 'immune responses' or 'neuronal adaptation'. Reword, or co-cite da Silveira 2020 (innate immunity, chronic inflammation, circadian rhythm).
  2. Ref 5 (Tullai 2007) is growth-factor stimulation of cultured cells, not neuronal activity. It supports the IEG / delayed-primary / secondary classification; keep it for that and let Tyssowski 2018 carry 'neuronal activity'.
  3. Ref 10 (Masarapu 2024): 'limited concordance with bulk expression from the same mice' is in the full text, not the abstract (11 shared DEGs out of 825 and 629, only 2 with the same direction). Supported, but cite the full-text result.
  4. Ref 22 (Johnson 2007) is cited for 'the official neuroCombat 0.2.12 implementation'. neuroCombat is the Fortin et al. implementation: add Fortin 2018 (NeuroImage) at that sentence.
  5. Ref 34 (Winkler 2014) is cited for 'the Monte Carlo plus-one alternative'. The plus-one rule is Phipson & Smyth 2010: add it there.
- **Evidence bearing on the protocol-confound reframe (state it as found):**
  - Lai Polo et al. 2020 (iScience) already showed that rodent spaceflight RNA-seq is confounded by sample collection: 'Changes due to preservation method eclipsed those between flight and ground samples'. The general claim 'collection protocol confounds rodent spaceflight transcriptomes' is therefore **not new**; Paper 1's contribution has to be the brain-specific IEG signature, the exact aliasing in RRRM-2 (euthanasia agent and freezing method with cohort), and the post hoc correction-gain result. Cite Lai Polo 2020 prominently and position against it (relevant to Worker B).
  - Staib-Lasarzik 2014: a few minutes of anaesthesia before sampling change hippocampal mRNA, mainly early-regulated genes; Hamaya 2000: anaesthetics transiently raise brain c-fos mRNA; Bunting 2016: high-dose isoflurane blocks initiation of Arc/Egr1 transcription. **None of these tested ketamine/xylazine**, so they support the principle, not the specific agent.
  - Counter-evidence that must be cited with any ketamine argument: Nakao 1993 found **no** c-Fos protein in rat hippocampus after 100 mg/kg ketamine (c-Fos was present in neocortex and thalamus). Davoudian 2023 found c-Fos elevation (incl. CA1) after **subanaesthetic** ketamine. Both are protein read-outs on a timescale longer than an overdose-to-dissection interval; neither settles mRNA within minutes.
  - Chakkingal Bhaskaran 2023 (pigs): euthanasia method did **not** change protein-coding expression; storage condition (RNAlater vs liquid nitrogen) did. Dachet 2021: activity-dependent transcripts fall within hours ex vivo. Together these point at freezing method / time-to-preservation as at least as plausible as the anaesthetic. Neither RRRM-2 factor can be separated from the other by the deposited design.
- npj Microgravity submission guidelines (read 15 Oct 2026, cached) state no reference-count cap; they allow datasets in the numbered list. With 73 literature references plus 11 datasets the list would be 84 entries; the parent should cut to what the rewritten text actually cites.

## Counts

| quantity | n |
|---|---|
| non-dataset references in audit | 73 |
|   re-verified from manuscript (refs 1-22, 34-38) | 27 |
|   new | 46 |
| fully verified, peer-reviewed (journal + proceedings) | 69 |
|   journal articles | 68 |
|   archival proceedings, no DOI (PMLR) | 1 |
| bibliographic fields verified, supporting text NOT verified | 3 |
| book (not counted as peer-reviewed) | 1 |
| retracted / expression of concern | 0 |
| with a published correction | 6 |
| additional quotes checked beyond the primary quote | 19 |
| indexed in Europe PMC with a PMID | 61 |
|   of which with a PMCID (free full text) | 51 |
| Crossref type posted-content (preprint) | 0 |
| OSDR dataset entries (separate list) | 11 |
|   with a resolving dataset DOI | 11 |
|   DOI already in manuscript | 2 |
| dropped from the manuscript's 27 | 0 |

Nothing from the manuscript's 27 was dropped: all exist, resolve and match. Boyd 2004 stays as a book reference; Wold 2001 stays (standard PLS reference) with its supporting text marked NOT VERIFIED.

## Reference table

Status VERIFIED = DOI resolves (doi.org handle API), Crossref record matches, Europe PMC cross-check (title, volume, year, author count) where indexed, no retraction/EoC in Crossref `updates` filter or Europe PMC, and the quote in `C_reference_audit.json` was found verbatim at source. Topic groups: a spaceflight; b collection/anaesthesia; c activity-regulated transcription; d unwanted variation; e validation/inference; f linear algebra/optimisation; g circadian/erythroid.

| # | key | short citation | journal | year | origin | grp | status | supported claim (sentence / topic) |
|---|---|---|---|---|---|---|---|---|
| 1 | `Afshinnekoo2020spaceflight` | Afshinnekoo et al. 2020 | Cell | 2020 | old ref 1 | a | VERIFIED; correction exists | Reviews of spaceflight biology identify oxidative stress, mitochondrial regulation, immune responses and neuronal adaptation as recurring themes [1]. |
| 2 | `Gebre2025osdr` | Gebre et al. 2025 | Nucleic Acids Research | 2025 | old ref 2 | a | VERIFIED | NASA's Open Science Data Repository (OSDR) and consensus RNA-seq processing make comparative analyses possible [2,3]. |
| 3 | `Overbey2021genelab` | Overbey et al. 2021 | iScience | 2021 | old ref 3 | a | VERIFIED | NASA's Open Science Data Repository (OSDR) and consensus RNA-seq processing make comparative analyses possible [2,3]. |
| 4 | `Tyssowski2018activity` | Tyssowski et al. 2018 | Neuron | 2018 | old ref 4 | c | VERIFIED | Neuronal activity can induce distinct rapid and secondary gene-expression programs [4,5]. |
| 5 | `Tullai2007primary` | Tullai et al. 2007 | Journal of Biological Chemistry | 2007 | old ref 5 | c | VERIFIED | Neuronal activity can induce distinct rapid and secondary gene-expression programs [4,5]. |
| 6 | `Marsh2022artifactual` | Marsh et al. 2022 | Nature Neuroscience | 2022 | old ref 6 | b | VERIFIED; correction exists | Tissue processing can generate artifactual or confounding glial signatures [6]. |
| 7 | `Overbey2022collection` | Overbey et al. 2022 | Cell Reports Methods | 2022 | old ref 7 | a/b | VERIFIED | Sample collection during spaceflight introduces additional logistical constraints that affect the suitability of material for single-cell and spatial assays [7]. |
| 8 | `Holley2022brain` | Holley et al. 2022 | npj Microgravity | 2022 | old ref 8 | a | VERIFIED | Previous mouse studies reported spaceflight-associated expression changes and oxidative injury in brain tissue [8,9]. |
| 9 | `Mao2020bbb` | Mao et al. 2020 | The FASEB Journal | 2020 | old ref 9 | a | VERIFIED | Previous mouse studies reported spaceflight-associated expression changes and oxidative injury in brain tissue [8,9]. |
| 10 | `Masarapu2024spatial` | Masarapu et al. 2024 | Nature Communications | 2024 | old ref 10 | a | VERIFIED | Spatial and single-nucleus multiomics in RR-3 localized responses to regions and cell populations, and showed limited concordance with bulk expression from the same mice [10]. |
| 11 | `Grandke2026mirna` | Grandke et al. 2026 | Nature Communications | 2026 | old ref 11 | a | VERIFIED | Recent small-RNA analysis across 686 mouse samples from 13 organs reported predominantly age-independent spaceflight responses [11]. |
| 12 | `Roy2025inflammasome` | Roy et al. 2025 | Cells | 2025 | old ref 12 | a | VERIFIED | A protein-level RR-1 study reported decreased brain inflammasome activation, with no significant flight-ground differences in the assayed cytokines [12]. |
| 13 | `Leek2010batch` | Leek et al. 2010 | Nature Reviews Genetics | 2010 | old ref 13 | d | VERIFIED | Unwanted variation can dominate high-dimensional expression [13,14]. |
| 14 | `Nygaard2016batch` | Nygaard et al. 2016 | Biostatistics | 2016 | old ref 14 | d | VERIFIED | Adjustments that preserve a chosen group contrast can exaggerate downstream confidence [13,14]. |
| 15 | `Leek2007sva` | Leek & Storey 2007 | PLoS Genetics | 2007 | old ref 15 | d | VERIFIED | Surrogate-variable and control-gene methods illustrate how auxiliary information can help model nuisance structure [15,16,17]. |
| 16 | `Risso2014ruvseq` | Risso et al. 2014 | Nature Biotechnology | 2014 | old ref 16 | d | VERIFIED | Control-gene methods model nuisance structure [15-17]; RUVg follows the published RUVSeq algorithm [16]. |
| 17 | `GagnonBartsch2012ruv` | Gagnon-Bartsch & Speed 2012 | Biostatistics | 2012 | old ref 17 | d | VERIFIED | Control-gene methods model nuisance structure [15,16,17]. |
| 18 | `Hoerl1970ridge` | Hoerl & Kennard 1970 | Technometrics | 1970 | old ref 18 | f | VERIFIED | Ridge regression [18] is compared under identical held-out animals; ridge references motivate the penalty. |
| 19 | `Wold2001pls` | Wold et al. 2001 | Chemometrics and Intelligent Laboratory Systems | 2001 | old ref 19 | f | BIBLIO ONLY (no abstract at source) | Principal-component regression and partial least squares [19] are compared under identical held-out animals. |
| 20 | `Boyd2004convex` | Boyd & Vandenberghe 2004 |  | 2004 | old ref 20 | f | VERIFIED (book, not journal) | A multi-sensor convex formulation places a spectral gain constraint inside the prediction objective [20,21]; the spectral ball is a standard convex constraint [20]. |
| 21 | `Beck2009fista` | Beck & Teboulle 2009 | SIAM Journal on Imaging Sciences | 2009 | old ref 21 | f | VERIFIED | Accelerated steps use a monotone objective restart and step size 1/lambda_max(H), following established proximal-gradient principles [21]. |
| 22 | `Johnson2007combat` | Johnson et al. 2007 | Biostatistics | 2007 | old ref 22 | d | VERIFIED | Known-truth count simulations compare the algorithm with established batch and control-gene methods [16,17,22]. |
| 23 | `Winkler2014permutation` | Winkler et al. 2014 | NeuroImage | 2014 | old ref 34 | e | VERIFIED | Exact whole-animal permutation inference for panel contrasts [34]. |
| 24 | `Benjamini1995fdr` | Benjamini & Hochberg 1995 | Journal of the Royal Statistical Society Series B: Statistical Methodology | 1995 | old ref 35 | e | VERIFIED | All 232 panel probabilities across primary contrasts and sensitivities were corrected together using BH [35]. |
| 25 | `Love2014deseq2` | Love et al. 2014 | Genome Biology | 2014 | old ref 36 | e | VERIFIED | PyDESeq2 0.5.4 fitted seven per-deposit count models [36,37]. |
| 26 | `Muzellec2023pydeseq2` | Muzellec et al. 2023 | Bioinformatics | 2023 | old ref 37 | e | VERIFIED | PyDESeq2 0.5.4 fitted seven per-deposit count models [36,37]. |
| 27 | `Morris2019simulation` | Morris et al. 2019 | Statistics in Medicine | 2019 | old ref 38 | e | VERIFIED | Simulations follow explicit aims, data-generating mechanisms, estimands, methods and performance measures [38]. |
| 28 | `GarrettBakelman2019twins` | Garrett-Bakelman et al. 2019 | Science | 2019 | new | a | VERIFIED | A year-long mission produced transcriptional, metabolic and other molecular changes in a human astronaut compared with his twin (spaceflight context for the opening paragraph). |
| 29 | `Overbey2024soma` | Overbey et al. 2024 | Nature | 2024 | new | a | VERIFIED | Large public spaceflight omics resources now support comparative re-analysis (with [2,3]). |
| 30 | `daSilveira2020mitochondria` | da Silveira et al. 2020 | Cell | 2020 | new | a | VERIFIED | Multi-omic re-analysis of GeneLab rodent and astronaut data identifies mitochondrial stress as a recurring spaceflight theme (with [1]). |
| 31 | `Choi2020habitat` | Choi et al. 2020 | Scientific Reports | 2020 | new | a/b | VERIFIED | NASA's RR-1 validation established on-orbit euthanasia with tissues frozen in situ and of acceptable RNA quality, the collection mode used for ISS-terminal cohorts. |
| 32 | `LaiPolo2020collection` | Lai Polo et al. 2020 | iScience | 2020 | new | a/b | VERIFIED | Sample collection technique has already been shown to confound RNA-seq of rodent spaceflight experiments; the ISS-T versus LAR contrast in RRRM-2 is aliased with euthanasia agent and freezing method. |
| 33 | `Ronca2019behavior` | Ronca et al. 2019 | Scientific Reports | 2019 | new | a | VERIFIED; correction exists | Mice in the ISS Rodent Habitat were more active than identically housed ground controls and developed spaceflight-specific circling behaviour, a behavioural difference between flight and ground animals that can itself drive activity-regulated genes. |
| 34 | `Roberts2017mri` | Roberts et al. 2017 | New England Journal of Medicine | 2017 | new | a | VERIFIED | Long-duration spaceflight is associated with structural brain changes in astronauts, motivating rodent brain transcriptomics. |
| 35 | `StaibLasarzik2014anesthesia` | Staib-Lasarzik et al. 2014 | Journal of Neurotrauma | 2014 | new | b | VERIFIED | The anaesthetic used for euthanasia alone changes brain mRNA expression in healthy mice, so ketamine/xylazine (ISS-T) versus gas anaesthesia (LAR) is a candidate cause of the cohort effect. |
| 36 | `Hamaya2000anesthetics` | Hamaya et al. 2000 | Anesthesia & Analgesia | 2000 | new | b | VERIFIED | General anaesthetics alter immediate-early gene expression in vital organs, so the IEG signal is plausibly agent dependent. |
| 37 | `Bunting2016isoflurane` | Bunting et al. 2016 | Frontiers in Behavioral Neuroscience | 2016 | new | b | VERIFIED; EPMC differs: year | High-dose isoflurane prevents initiation of Arc and Zif268 (Egr1) transcription in rat hippocampus and cortex, so a gas-anaesthesia cohort (LAR) is expected to carry less newly induced IEG mRNA than an injectable-overdose cohort (ISS-T). |
| 38 | `Davoudian2023ketamine` | Davoudian et al. 2023 | ACS Chemical Neuroscience | 2023 | new | b | VERIFIED | Subanaesthetic ketamine elevates c-Fos in many mouse brain regions, including hippocampal CA1 (evidence that ketamine is an IEG inducer; dose and read-out differ from a ketamine/xylazine overdose). |
| 39 | `Overmyer2015euthanasia` | Overmyer et al. 2015 | PLOS ONE | 2015 | new | b | VERIFIED | Anaesthesia and euthanasia method alter the molecular (metabolite) profile of mouse tissues, with ketamine the most variable pre-euthanasia method; brain was not assayed. |
| 40 | `ChakkingalBhaskaran2023euthanasia` | Chakkingal Bhaskaran et al. 2023 | BMC Genomics | 2023 | new | b | VERIFIED | Storage condition (RNAlater vs liquid-nitrogen snap freezing) changed differential expression in several pig tissues, whereas euthanasia method (nitrogen anoxia vs T-61) changed only snRNAs; RRRM-2 aliases euthanasia agent and freezing method with collection cohort. |
| 41 | `vandenBrink2017dissociation` | van den Brink et al. 2017 | Nature Methods | 2017 | new | b | VERIFIED | Enzymatic tissue dissociation strongly alters the transcriptome of a cell subpopulation (muscle stem cells), a technical artefact that can be mistaken for biology. |
| 42 | `Wu2017actseq` | Wu et al. 2017 | Neuron | 2017 | new | b/c | VERIFIED | Conventional whole-cell dissociation artificially induces transcriptional perturbations that confound detection of activity-induced changes in brain cells. |
| 43 | `Lacar2016nuclear` | Lacar et al. 2016 | Nature Communications | 2016 | new | b/c | VERIFIED; correction exists | Brief novel-environment exposure induces IEGs (Fos, Arc, Egr1) and large transcriptome changes in mouse dentate granule cells, so short pre-collection experience is itself a hippocampal IEG stimulus. |
| 44 | `Hrvatin2018visual` | Hrvatin et al. 2018 | Nature Neuroscience | 2018 | new | c | VERIFIED; correction exists | Sensory stimulation evokes divergent transcriptional responses in every cortical cell type, including non-neuronal cells, so an activity signature in bulk tissue mixes cell types. |
| 45 | `Ferreira2018ischemia` | Ferreira et al. 2018 | Nature Communications | 2018 | new | b | VERIFIED | Death and post-mortem ischaemic time change tissue transcriptomes, so time from death to preservation is a nuisance variable. |
| 46 | `Dachet2021postmortem` | Dachet et al. 2021 | Scientific Reports | 2021 | new | b | VERIFIED | In fresh human neocortex held ex vivo to mimic the post-mortem interval, neuronal activity-dependent transcripts fell selectively within a few hours; time to preservation is therefore a nuisance variable for an IEG panel. |
| 47 | `Yap2018activity` | Yap & Greenberg 2018 | Neuron | 2018 | new | c | VERIFIED | Review of activity-dependent transcription in the vertebrate brain (background for the rPRG/SRG panels, with [4,5]). |
| 48 | `Saha2011poised` | Saha et al. 2011 | Nature Neuroscience | 2011 | new | c | VERIFIED | Rapid IEGs such as Arc are poised for near-instantaneous transcription by paused RNA polymerase II, so they can respond within a short peri-collection window. |
| 49 | `Zhang2020combatseq` | Zhang et al. 2020 | NAR Genomics and Bioinformatics | 2020 | new | d | VERIFIED | Count-based batch adjustment for RNA-seq (ComBat-seq) is an established comparator. |
| 50 | `Molania2023prps` | Molania et al. 2023 | Nature Biotechnology | 2023 | new | d | VERIFIED | Unwanted variation compromises downstream analyses of large RNA-seq studies; RUV-III with pseudo-replicates of pseudo-samples (PRPS) removes it using samples of shared biology. |
| 51 | `Wang2017confounder` | Wang et al. 2017 | The Annals of Statistics | 2017 | new | d | VERIFIED | Latent confounders correlated with the primary variable bias multiple testing; adjustment is identified only under stated conditions (negative controls, or sparse non-nulls). |
| 52 | `Gerstner2016sleep` | Gerstner et al. 2016 | BMC Genomics | 2016 | new | d/c | VERIFIED | Removal of unwanted variation (RUV) improved detection of differential expression after sleep deprivation in mouse cortex. |
| 53 | `Zindler2020combat` | Zindler et al. 2020 | BMC Bioinformatics | 2020 | new | d | VERIFIED | Batch correction with ComBat can systematically introduce false-positive results, in unbalanced and even balanced designs (with [14]). |
| 54 | `Peixoto2015analysis` | Peixoto et al. 2015 | Nucleic Acids Research | 2015 | new | d | VERIFIED | In experimental-neuroscience RNA-seq, batch effects can dominate the signal of interest and the normalisation choice changes power and reproducibility. |
| 55 | `Leek2014svaseq` | Leek 2014 | Nucleic Acids Research | 2014 | new | d | VERIFIED | Surrogate-variable estimation for sequencing data (svaseq) is the count-data analogue of [15]. |
| 56 | `Fortin2018neurocombat` | Fortin et al. 2018 | NeuroImage | 2018 | new | d | VERIFIED | The neuroCombat implementation used as comparator originates in the cortical-thickness harmonisation work (cite with [22]). |
| 57 | `Varma2006bias` | Varma & Simon 2006 | BMC Bioinformatics | 2006 | new | e | VERIFIED | Tuning inside the validation loop biases error estimates; hence nested selection with training-only preprocessing. |
| 58 | `Ambroise2002selection` | Ambroise & McLachlan 2002 | Proceedings of the National Academy of Sciences | 2002 | new | e | VERIFIED | Gene selection outside cross-validation produces selection bias; all filtering and normalisation references are fitted within training animals. |
| 59 | `Benjamini2001dependency` | Benjamini & Yekutieli 2001 | The Annals of Statistics | 2001 | new | e | VERIFIED | BY correction controls FDR under arbitrary dependence; no result in the new family passed BY correction. |
| 60 | `Lazic2018unit` | Lazic et al. 2018 | PLOS Biology | 2018 | new | e | VERIFIED | Pseudoreplication inflates N and produces false positives; the animal, not the library or tissue region, is the experimental unit here. |
| 61 | `PercieduSert2020arrive` | Percie du Sert et al. 2020 | PLOS Biology | 2020 | new | e | VERIFIED | ARRIVE 2.0 is the reporting checklist for in vivo experiments (Essential 10: minimum requirement); cite for reporting of experimental units and design recommendations. |
| 62 | `Button2013power` | Button et al. 2013 | Nature Reviews Neuroscience | 2013 | new | e | VERIFIED; correction exists | Small samples reduce power and inflate effect sizes of nominally significant findings in neuroscience. |
| 63 | `Phipson2010permutation` | Phipson & Smyth 2010 | Statistical Applications in Genetics and Molecular Biology | 2010 | new | e | VERIFIED | Permutation p-values should never be zero; the Monte Carlo plus-one alternative [34] follows this rule. |
| 64 | `Efron1979bootstrap` | Efron 1979 | The Annals of Statistics | 1979 | new | e | VERIFIED | Percentile intervals use the bootstrap resampling of animals. |
| 65 | `Szyld2006oblique` | Szyld 2006 | Numerical Algorithms | 2006 | new | f | VERIFIED | For a nontrivial oblique projector P, norm(P) = norm(I - P); this identity links the norm of the score-removal operator to the norm of the removed score projector. |
| 66 | `Frank1956quadratic` | Frank & Wolfe 1956 | Naval Research Logistics Quarterly | 1956 | new | f | BIBLIO ONLY (no abstract at source) | Original Frank-Wolfe (conditional-gradient) algorithm, cited with Jaggi 2013 for the duality-gap certificate. |
| 67 | `Izenman1975rrr` | Izenman 1975 | Journal of Multivariate Analysis | 1975 | new | f | BIBLIO ONLY (no abstract at source) | Constrained multi-sensor regression is related to reduced-rank multivariate linear regression, which constrains the coefficient matrix spectrum by rank rather than by gain. |
| 68 | `ODonoghue2015restart` | O’Donoghue & Candès 2015 | Foundations of Computational Mathematics | 2015 | new | f | VERIFIED | Accelerated steps use an objective-based restart (with [21]). |
| 69 | `Zhang2014circadian` | Zhang et al. 2014 | Proceedings of the National Academy of Sciences | 2014 | new | g | VERIFIED | Forty-three per cent of mouse protein-coding genes oscillate with time of day in some organ, largely organ-specifically, so dissection time can confound a circadian panel contrast. |
| 70 | `Takahashi2017clock` | Takahashi 2017 | Nature Reviews Genetics | 2017 | new | g | VERIFIED | The mammalian circadian clock is a cell-autonomous transcription-translation autoregulatory feedback loop (background for the circadian panel). |
| 71 | `Biagioli2009globin` | Biagioli et al. 2009 | Proceedings of the National Academy of Sciences | 2009 | new | g | VERIFIED; EPMC differs: title | Haemoglobin transcripts are not exclusively erythroid: alpha- and beta-globin are expressed by dopaminergic neurons and some astrocytes and oligodendrocytes, so an Hb-based erythroid score is an imperfect blood-content proxy. |
| 72 | `Nakao1993ketamine` | Nakao et al. 1993 | Neuroscience Letters | 1993 | new | b | VERIFIED | Counter-evidence to cite with [Davoudian]: high-dose ketamine did not induce c-Fos protein in rat hippocampus, so the ketamine/xylazine explanation of the hippocampal rPRG excess is not established by prior data. |
| 73 | `Jaggi2013frankwolfe` | Jaggi 2013 | Proceedings of the 30th International Conference on Machine Learning | 2013 | new | f | VERIFIED | The Frank-Wolfe (conditional-gradient) duality gap is a certificate of suboptimality for constrained convex problems, including spectral/nuclear-norm balls. |

## Datasets (kept separately, `C_datasets.bib`)

| ms ref | accession | DOI | DOI source | resolves | in manuscript | title matches manuscript |
|---|---|---|---|---|---|---|
| 23 | OSD-202 | 10.26030/ewfb-7g23 | OSDR study metadata | True | False | True |
| 24 | OSD-352 | 10.26030/jm59-zy54 | OSDR study metadata | True | False | True |
| 25 | OSD-374 | 10.26030/cvn0-sp98 | OSDR study metadata | True | False | True |
| 26 | OSD-611 | 10.26030/t0z5-5g08 | OSDR study metadata | True | False | True |
| 27 | OSD-525 | 10.26030/3r5e-6b84 | OSDR study metadata | True | True | True |
| 28 | OSD-561 | 10.26030/yj0g-1273 | DataCite (title + landing URL) | True | False | True |
| 29 | OSD-562 | 10.26030/ay53-g970 | DataCite (title + landing URL) | True | False | True |
| 30 | OSD-563 | 10.26030/ycn8-v994 | DataCite (title + landing URL) | True | False | True |
| 31 | OSD-564 | 10.26030/48p5-yy05 | DataCite (title + landing URL) | True | False | True |
| 32 | OSD-612 | 10.26030/66jw-j748 | OSDR study metadata | True | True | True |
| 33 | OSD-613 | 10.26030/cevt-ja62 | DataCite (title + landing URL) | True | False | True |

Nature-style dataset lines (for the reference list):

- NASA Open Science Data Repository. Low dose (0.04 Gy) irradiation (LDR) and hindlimb unloading (HLU) microgravity in mice: brain transcriptomic and epigenomic data. NASA OSDR OSD-202 https://doi.org/10.26030/ewfb-7g23 (2018).
- NASA Open Science Data Repository. Multi-modal transcriptional and chromatin accessibility analysis of brains from mice flown on the RR-3 mission. NASA OSDR OSD-352 https://doi.org/10.26030/jm59-zy54 (2021).
- NASA Open Science Data Repository. Bi-directional and shared epigenomic signatures following proton irradiation. NASA OSDR OSD-374 https://doi.org/10.26030/cvn0-sp98 (2022).
- NASA Open Science Data Repository. Postsynaptic density radiation signature following space irradiation. NASA OSDR OSD-611 https://doi.org/10.26030/t0z5-5g08 (2025).
- NASA Open Science Data Repository. Transcriptional profiling of cerebellums from mice flown on the RR-23 mission. NASA OSDR OSD-525 https://doi.org/10.26030/3r5e-6b84 (2022).
- NASA Open Science Data Repository. Bulk RNA sequencing and spatially resolved transcriptional profiling of cerebellums from mice flown on the RRRM-2 mission. NASA OSDR OSD-561 https://doi.org/10.26030/yj0g-1273 (2022).
- NASA Open Science Data Repository. Bulk RNA sequencing and spatially resolved transcriptional profiling of hippocampi from mice flown on the RRRM-2 mission. NASA OSDR OSD-562 https://doi.org/10.26030/ay53-g970 (2022).
- NASA Open Science Data Repository. Bulk RNA sequencing and spatially resolved transcriptional profiling of cerebellums from mice flown on the RR-10 mission. NASA OSDR OSD-563 https://doi.org/10.26030/ycn8-v994 (2022).
- NASA Open Science Data Repository. Bulk RNA sequencing and spatially resolved transcriptional profiling of hippocampi from mice flown on the RR-10 mission. NASA OSDR OSD-564 https://doi.org/10.26030/48p5-yy05 (2022).
- NASA Open Science Data Repository. Bulk RNA sequencing and single nuclei transcriptomics and epigenomics of brain tissue from mice flown on the RR-10 mission. NASA OSDR OSD-612 https://doi.org/10.26030/66jw-j748 (2023).
- NASA Open Science Data Repository. Bulk RNA sequencing and single nuclei transcriptomics and epigenomics of brain tissue from mice flown on the RRRM-2 mission. NASA OSDR OSD-613 https://doi.org/10.26030/cevt-ja62 (2023).

## Candidates from the brief that were replaced, re-sourced or not used

| candidate | outcome |
|---|---|
| West & Greenberg 2011, Cold Spring Harb. Perspect. Biol. | Commissioned review; Yap & Greenberg 2018 Neuron covers the same ground in a preferred journal. |
| Molania et al. 2019, Nucleic Acids Res. (RUV-III for NanoString) | NanoString-specific; replaced by Molania et al. 2023 Nat. Biotechnol. (RUV-III/PRPS for RNA-seq). |
| Jaggi 2013 ICML | Included, but as archival proceedings without a DOI (PMLR v28); verified on the PMLR page rather than Crossref. |
| Zindler et al. 2020 'PeerJ' | Venue in the brief was wrong; the verified paper is BMC Bioinformatics 21, 271. |
| Choi et al. 2020 'RR-1 rodent habitat validation' | Included; venue is Scientific Reports 10, 2336. |
| Ketamine-induced IEG expression in rodent brain | No transcriptomic study of ketamine/xylazine overdose found in Europe PMC title searches; used Davoudian 2023 (subanaesthetic ketamine c-Fos) plus counter-evidence Nakao 1993 (high-dose ketamine, no hippocampal c-Fos). |
| Kremsky et al. 2023 Int. J. Mol. Sci. | Spaceflight brain antioxidant study; not needed for any Paper 1 sentence. |
| Erythroid/blood-contamination in brain RNA-seq (perfusion studies) | Europe PMC title search found no peer-reviewed study quantifying blood-transcript contamination in perfused vs non-perfused rodent brain bulk RNA-seq; Biagioli 2009 used for the narrower point that Hb transcripts are also neural/glial. |
| Recht, Fazel & Parrilo 2010 SIAM Rev. | Nuclear/spectral norm duality is covered by Boyd & Vandenberghe and Jaggi 2013; not added. |

## Per-entry caveats

- `Afshinnekoo2020spaceflight`: MANUSCRIPT CLAIM ONLY PARTLY SUPPORTED: the abstract's six features are oxidative stress, DNA damage, mitochondrial dysregulation, epigenetic changes, telomere length and microbiome shifts; it does not list immune responses or neuronal adaptation. Reword to the six features, or co-cite daSilveira2020mitochondria (innate immunity, inflammation, circadian rhythm).
- `Tullai2007primary`: Growth-factor stimulation of cultured cells, not neuronal activity; supports the IEG / delayed-primary / secondary classification, not 'neuronal activity' per se (Tyssowski 2018 covers the neuronal part).
- `Holley2022brain`: Male C57BL/6, 35-day SpaceX-12 mission, live return with tissue collected 38 +/- 4 h after splashdown (a LAR-type design), NanoString panel.
- `Mao2020bbb`: Same SpaceX-12 cohort as Holley 2022 (live return, 38 +/- 4 h).
- `Masarapu2024spatial`: The bulk-concordance part of the manuscript sentence is supported by the full text (11 shared DEGs, 2 with the same direction), not by the abstract. Female mice, RR-3.
- `Grandke2026mirna`: 'Predominantly age-independent' is a fair paraphrase of 'Age-dependent changes, smaller in magnitude compared to age-independent changes'. Samples are female mice.
- `Roy2025inflammasome`: Claim supported. Note the comparison is also against a launch-day baseline group.
- `Nygaard2016batch`: Microarray setting; the argument is design-generic.
- `Wold2001pls`: No abstract carried by Crossref, Europe PMC, OpenAlex or Semantic Scholar; ScienceDirect returns HTTP 400/403 to scripted and web_fetch requests. Bibliographic fields verified; supporting text NOT VERIFIED. Standard PLS tutorial; acceptable for 'partial least squares' but the parent should confirm the PLS variant used.
- `Boyd2004convex`: Cambridge University Press monograph (Crossref type 'monograph'). Retain as a book reference; NOT counted as a peer-reviewed journal article.
- `Johnson2007combat`: Manuscript cites [22] for 'the official neuroCombat 0.2.12 implementation'. neuroCombat is the Fortin et al. implementation; add Fortin2018neurocombat beside [22] there.
- `Winkler2014permutation`: Manuscript cites [34] for 'the Monte Carlo plus-one alternative'; the plus-one rule is Phipson2010permutation. Add it there.
- `Benjamini1995fdr`: Quote reproduces the source spelling 'Bonferronitype'.
- `Choi2020habitat`: RR-1 (female C57BL/6NTac 32 wk and C57BL/6J 16 wk): on-orbit euthanasia with tissues frozen in situ. Journal is Scientific Reports.
- `LaiPolo2020collection`: Central prior evidence for the reframe: preservation/collection protocol effects eclipsed flight-ground differences in rodent spaceflight RNA-seq. Peer-reviewed version (iScience) used; a bioRxiv preprint (10.1101/2020.07.18.209775) also exists and is NOT cited.
- `Ronca2019behavior`: An Author Correction exists (Sci. Rep. 9, 10.1038/s41598-019-45362-1); recorded, not a retraction.
- `StaibLasarzik2014anesthesia`: Hippocampus, qPCR of a 10-transcript panel in mice; anaesthetics tested were isoflurane, midazolam/fentanyl/medetomidine and chloral hydrate - NOT ketamine/xylazine. Supports the general principle only.
- `Hamaya2000anesthetics`: Rat; isoflurane, pentobarbital and propofol (no ketamine/xylazine). c-fos/c-jun RT-PCR.
- `Bunting2016isoflurane`: Rat catFISH; Frontiers in Behavioral Neuroscience volume 9 (2015 volume), article 363, published online January 2016 (Crossref year 2016). Key uses 2016.
- `Davoudian2023ketamine`: Subanaesthetic ketamine, c-Fos protein (whole-brain imaging), not a ketamine/xylazine overdose and not mRNA. Pair with Nakao1993ketamine (counter-evidence).
- `Overmyer2015euthanasia`: Metabolomics, not transcriptomics; tissues were muscle, liver, heart, adipose and serum (no brain).
- `ChakkingalBhaskaran2023euthanasia`: Pigs (piglets), 3'mRNA-seq of pituitary, hypothalamus, liver, lung. Euthanasia method did NOT change protein-coding expression; storage condition (RNAlater vs LN2) did. Cite for storage/preservation, and state the null euthanasia result honestly.
- `vandenBrink2017dissociation`: Muscle stem cells, not brain. Nature Methods Brief Communication without an indexed abstract; the supporting text is the article's opening paragraph on nature.com.
- `Wu2017actseq`: Mouse amygdala/brain, single-cell dissociation artefact.
- `Ferreira2018ischemia`: Human GTEx post-mortem tissues; supports time-to-preservation as a covariate, by analogy only for rodent collection.
- `Dachet2021postmortem`: Human neocortex held ex vivo after surgical resection; activity-dependent transcripts fall within hours. Relevant to time-to-freezing (cryochiller vs LN2) by analogy.
- `Molania2023prps`: Crossref published-online 2022, print issue 2023 (Nat. Biotechnol. 41, 82-95); year 2023 used.
- `Wang2017confounder`: Identification requires negative controls or sparse non-nulls; a nuisance perfectly aliased with the primary variable is not covered by either condition (inference by the parent, not a quote).
- `Zindler2020combat`: Brief candidate listed 'PeerJ'; the paper is in BMC Bioinformatics (Crossref verified). DNA methylation arrays.
- `Fortin2018neurocombat`: Source of the neuroCombat implementation used in the comparator benchmark.
- `Szyld2006oblique`: Springer landing page serves a bot wall to scripted clients; abstract read with the web_fetch tool and transcribed verbatim (tool tokenisation artefacts kept). Bibliographic fields from Crossref.
- `Frank1956quadratic`: 1956 article without an abstract in any index; Wiley returns HTTP 403. Bibliographic fields verified; supporting text NOT VERIFIED. Jaggi2013frankwolfe (abstract verified) supports the duality-gap statement; cite both or Jaggi alone.
- `Izenman1975rrr`: No abstract in Crossref/Europe PMC/OpenAlex/Semantic Scholar; ScienceDirect open-archive PDF returns HTTP 403 to scripted requests. Bibliographic fields verified; supporting text NOT VERIFIED. Optional reference; drop if the Discussion does not mention reduced-rank regression.
- `ODonoghue2015restart`: Crossref published-online 2013, print issue 2015 (Found. Comput. Math. 15, 715-732); year 2015 used. Abstract from Semantic Scholar.
- `Zhang2014circadian`: 12 mouse organs including brainstem, cerebellum and hypothalamus (organ list from the paper, not the abstract).
- `Biagioli2009globin`: Supports caution in calling Hba/Hbb an erythroid (blood) score. Alas2, Slc4a1 and Gypa are not covered by this source.
- `Nakao1993ketamine`: COUNTER-EVIDENCE. Rat, c-Fos protein immunohistochemistry. Include so the discussion does not overstate the ketamine explanation; mRNA within minutes of an overdose is a different read-out.
- `Jaggi2013frankwolfe`: PMLR volume 28 (ICML 2013), no DOI exists. Verified against the PMLR proceedings page (title, author, pages 427-435, abstract). Archival peer-reviewed proceedings; counted separately from journal articles.

## Method

- One script, `tools/c_verify.py` (fetch, then build), with all API responses cached as JSON under `cache/` (crossref, crossref_updates, doi_handle, europepmc, openalex, semanticscholar, datacite, osdr_meta). Rate: >= 1.2 s between Crossref calls, back-off on HTTP 429/5xx; plain User-Agent, no e-mail address sent.
- Bibliographic fields are taken from Crossref, not typed; Europe PMC is the independent cross-check. Year = print-issue year where Crossref has one, else online year (Molania 2023, O'Donoghue 2015 have earlier online dates; Bunting 2016 is in the 2015 volume of Front. Behav. Neurosci.).
- Quotes are in `tools/c_quotes.py`; the build asserts each is an exact substring of an abstract (Europe PMC, Crossref, OpenAlex) or of a cached publisher text built by `tools/c_texts.py` (van den Brink 2017 opening paragraph, Masarapu 2024 full-text result, O'Donoghue 2015 Semantic Scholar abstract, Jaggi 2013 PMLR abstract). Szyld 2006's abstract was read with the web_fetch tool because Springer serves a bot wall to scripts; it is transcribed verbatim and flagged.
- Retraction/correction: Crossref `works?filter=updates:<DOI>` (notices that update the work) plus Europe PMC `commentCorrectionList` and `pubTypeList`. Correction notices were identified, not read.
- Datasets: `tools/c_datasets.py` reads live OSDR metadata, the deposited snapshot in the source package, DataCite and the doi.org handle API.
- `C_references.bib` + `C_datasets.bib` compile with tectonic/BibTeX (unsrt, TeX Gyre Heros) with no missing glyphs (test in scratch, not kept).
- Not done: the correction notices themselves were not read; abstracts were read but full texts were not (except Masarapu 2024); Wold 2001, Frank & Wolfe 1956 and Izenman 1975 supporting text could not be retrieved (ScienceDirect and Wiley return 403 to scripted and web_fetch requests). Unavailable MCP servers in this session: fusion (failed to start), cloudflare-api and github (awaiting authorization); none was needed.

## Files (all in `refs/C_refs/`)

- `C_references.bib`: 73 literature entries, keys AuthorYEARkeyword
- `C_datasets.bib`: 11 @misc OSDR entries with DOI, URL, howpublished, note 'Accessed 15 October 2026'
- `C_reference_audit.json`: per entry: key, DOI, PMID, PMCID, Crossref fields, Europe PMC fields, match checks, retraction status, notices, quote(s) and where found, supported claim, section, notes
- `C_references_nature.txt`: Nature-style list in the same order as the table
- `C_abstracts_review.md`: every abstract retrieved, for re-reading
- `tools/c_entries.py`: keys, DOIs, sections, claims
- `tools/c_quotes.py`: quotes, abbreviations, notes, unused candidates
- `tools/c_verify.py`: fetch + verify + emit
- `tools/c_texts.py`: builds cache/manual_texts.json
- `tools/c_datasets.py`: dataset verification
- `tools/c_report.py`: this report
- `tools/c_search.py`: discovery searches (Europe PMC / Crossref), cached under cache/search
- `cache/`: raw API responses and publisher pages

Rebuild: `cd refs/C_refs/tools && P=/Users/zain/.kiro/crew/scratch/runtime-735c0e1b/p1env/bin/python && $P c_texts.py && $P c_verify.py build && $P c_datasets.py && $P c_report.py` (offline from cache).
