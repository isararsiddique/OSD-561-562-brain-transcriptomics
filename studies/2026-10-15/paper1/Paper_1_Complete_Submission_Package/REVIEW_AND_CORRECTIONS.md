# Review of Antigravity's refreshed manuscripts

The 49-page Paper 1 and 46-page Paper 2 are preserved in a hashed source snapshot. All 221 reported sync checks were reproduced; the following issues remain outside that audit's effective coverage. The author supplied by the user is Isarar Siddique, with affiliation and email blank.

## 1. Author and affiliation

Both frozen frontmatters name Zainulabedin Khan and the Mumbai laboratory. The user supplied Isarar Siddique and requested blank affiliation and email.

Required correction: Replace author details consistently across manuscript sources, running headers, PDFs and refreshed packages. Do not retain the unconfirmed collaborative group or invent an affiliation.

## 2. Paper 2 rendered abstract is empty

PAPER_2_MANUSCRIPT.tex contains \abstract{} and \keyword{} even though the Markdown has an abstract. The 46-page PDF starts with an empty abstract box.

Required correction: Generate both manuscript formats from one evidence-backed abstract. This is repaired in the complete three-arm Paper 2 manuscript.

## 3. Draft displays publication identifiers

The template contains a generated 10.1038/s41526-026-0042-x string and received/accepted/published dates, plus journal publication branding. These do not establish an actual publication.

Required correction: Use neutral manuscript styling; remove assigned-looking DOI, volume/article identity, journal logo and acceptance metadata from draft outputs.

## 4. Paper 1 vivarium date contradicts raw metadata

Paper 1 introduction claims VIV was collected on 20 January 2021. All eight vivarium animals in the deposited ISA table have euthanasia and BSP dissection on 17 January, as does the updated sample manifest. Flight is 14 January.

Required correction: Correct both ground dates to 17 January. Preserve the flight/date alias; do not call the ground-ground contrast a three-day collection difference.

## 5. Paper 2 network headline does not follow current table

Frozen table: 288 rows, 7 local q<0.05 rows, 0 pooled q<0.05 rows, minimum pooled q=0.058775510. Abstract/graphical text says eight local rows and minimum pooled q=0.0792.

Required correction: Update all counts from the same frozen table and treat local families as exploratory. The completed paper uses the independently executed three-arm junction analysis instead of this unsupported headline.

## 6. Claimed significant RR23 uncoupling is absent in the frozen rerun

LIONESS_SingleSample_Edge: p=0.01439737, local q=0.05758947; Standardized_CoFluctuation_Product: p=0.02690251, local q=0.10761004

Required correction: Report the estimates and their corrected uncertainty. Do not describe a nominal/local result as a definitive selective flight mechanism.

## 7. Paper 1 mixes incompatible regression summaries

Abstract/conclusion: 22 contrasts, eight predictors, R2=0.6528. Frozen current table: n=25, two reported predictors, R2=0.653520112. For the displayed F(8,13)=3.056, the upper-tail p is 0.035937120, not 3.49e-7.

Required correction: Choose one reproducible model, derive every summary from it, and distinguish point-biserial association from a fitted-model F test. Shared animals across contrasts also require dependence-aware inference.

## 8. Residual correlation inference breaks design strata

p2_algo2a_within_animal_synchrony.py demeans by eight flight/route/age cells, then calls an unrestricted permutation of all residuals. q_bh_8cell_resid is calculated from ordinary Pearson p values rather than the saved permutation p values.

Required correction: Permute pairings within complete design cells or use a justified residual model with appropriate degrees of freedom. Recompute the specified family rather than retain the current residual p values.

## 9. Disjoint-arm GeoMx contrast does not remove slide aliasing

The spatial code permutes 12 sections across four labels and calls sections independent across slides. Each group occupies one slide. The main text calls the replacement contrast unconfounded and interprets DG attenuation as significant countermeasure damping.

Required correction: Disjoint arms remove shared-baseline covariance, but not the group/slide alias. Report spatial contrast ratios descriptively and remove causal rescue p-value claims unless a valid independent replication design is available.

## 10. Graphical abstract retains incorrect cohort identities

Paper 2 graphical text reverses RR23/OSD-525 and RR10/OSD-563/564/612, calls GSE295428 optical nerve/retina, and gives OSD-562 n=44. Current census supports OSD-525 RR23, OSD-563/564/612 RR10 and OSD-562 27 bulk assay rows.

Required correction: Rebuild graphical figures from the authoritative sample manifest. Do not retain invented modalities or placeholder profile-bar totals.

## 11. Kinetic class enrichments do not measure a terminal time window

Paper 1 graphical text calls expression-class gradients diagnostic of a 15-60 minute perimortem window. No observed per-animal transcription-onset timeline or calibrated terminal-delay experiment supports that temporal attribution.

Required correction: Describe enrichment of acute-response classes. Distinguish a plausible handling hypothesis from measured elapsed time and causal artifact identification.

## 12. 221 passing checks have limited coverage

The saved audit was independently rerun on the frozen snapshot and again passed all 221 checks. Its suites check selected existence, disjointness, identity regexes and numerical patterns in MD/TeX. They do not check rendered abstracts, user-approved authors, publication status, every narrative value or causal identification.

Required correction: Add independent scientific, raw-metadata and rendered-artifact checks. A passing count is not evidence of a validated biological contribution.

## Completed Paper 2

The replacement scientific manuscript is ../paper2_beta_20261014/manuscript/PAPER_2_MANUSCRIPT.md, with its matching PDF and editable LaTeX. It reports the independently executed three-arm junction analysis, uses the requested author, and omits publisher-assigned draft metadata. Its entire quantitative claim is generated from the verified all-animal analysis, not from length or the 221-check count.

Focused review of high-impact claims and artifact consistency, not an exhaustive revalidation of every 95-page argument. Source snapshots and hashes preserve the reviewed versions during concurrent updates.