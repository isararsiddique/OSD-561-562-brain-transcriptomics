# I — Background-corrected biological interpretation layer for Paper 2

Everything here tests the Paper 2 junction gene sets against the **measured gene
universe** (1,980 symbols from the 1,983 genes in the gene-level omnibus), never
against the mouse genome. A genome background would have manufactured enrichment
for every neuronal and RNA-handling term, because the junction pipeline only ever
measured genes that are expressed and multi-junction in cerebellum.

**Headline: nothing is significant.** 0 of 6,682 functional term tests and 0 of
102 cerebellar cell-type tests pass correction, for the flight set or the housing
set. The informative part is panel (b) of the figure: every observed odds ratio
sits below the smallest odds ratio its own test could have detected, so this is a
**power-limited null, not evidence of no biology**. The manuscript's existing
sentence — that the candidate identities "motivate RNA-handling, neuronal and
mitochondrial follow-up; they do not constitute tested pathway enrichment or
cell-type localization" — is now a tested statement rather than a disclaimer.

---

## 1. Gene sets (definitions, not new discoveries)

Built by `scripts/01_build_gene_sets.py` from
`extensions/gene_omnibus/results/gene_omnibus_complete.csv.gz`
(sha256 `05e245d0…38de9e85`), read-only.

| Set | n | Definition |
|---|---|---|
| `FLT_GC_BY` | 47 | joint BY q ≤ 0.05 in FLT−GC |
| `FLT_GC_top100` | **106** | strongest FLT−GC exact evidence; p ≤ 1.6454e-4 |
| `GC_VIV_top100` | 100 | strongest GC−VIV exact evidence; p ≤ 1.0366e-2 |
| background | 1,980 | every measured gene carrying a symbol |

Tie handling: `p_exact` is discrete (permutation-based), and 12 genes share the
p value of the 100th-ranked FLT−GC gene. All of them are kept, so the set is 106
rather than 100; breaking the tie on an effect-size statistic would have mixed
two kinds of evidence. The GC−VIV cut had 1 gene at the boundary, so that set is
exactly 100. Three measured genes lack a symbol and are dropped from the query
sets *and* the background, keeping the universe identical on both sides.

## 2. Functional enrichment — null under g:SCS and under BH

`scripts/03_gprofiler_enrichment.py`, `scripts/04_gprofiler_sensitivity_BH.py`.

- g:Profiler REST `gost/profile`, organism `mmusculus`, `domain_scope='custom'`
  with the 1,980-gene background, sources GO:BP/CC/MF, KEGG, Reactome,
  `significance_threshold_method='g_SCS'`, threshold 0.05.
- **API version `e114_eg62_p19_9dc92ac8`**, Ensembl 114 / GRCm39; GO classes
  `releases/2026-01-23`, KEGG FTP release 2026-03-15, Reactome 2026-3-20.
  Queried **2026-10-15 UTC**. Effective domain size 1,992 (symbol-to-Ensembl
  expansion of the 1,980 submitted). Raw JSON in `raw/`.
- Only public gene symbols were transmitted.

| Set | terms returned | g:SCS significant | strongest term | overlap | OR | uncorrected p | BH q |
|---|---|---|---|---|---|---|---|
| FLT−GC BY (47) | 2,151 | **0** | translational initiation (GO:0006413) | 4/10 | 30.9 | 5.2e-05 | 0.090 |
| FLT−GC top (106) | 3,505 | **0** | mitochondrion (GO:0005739) | 24/214 | 2.7 | 1.8e-04 | 0.205 |
| GC−VIV top (100) | 2,948 | **0** | neg. regulation of cell communication (GO:0010648) | 17/151 | 2.8 | 8.7e-04 | 0.444 |

Because `all_results=True` returns the *corrected* p capped at 1.0, the
uncorrected hypergeometric p and the Haldane-corrected odds ratio were recomputed
exactly from each term's 2×2 table, and BH applied within each set over terms of
5–500 genes (1,750 / 2,646 / 2,286 tests = 6,682). BH is strictly more permissive
than g:SCS, so this is a stronger null: minimum q = 0.090, 0.205, 0.444.

Genes behind the strongest term (`tables/headline_term_driving_genes.csv`):
*Eif2a, Eif4a2, Eif4h, Rps6kb1* for translational initiation; 13 of the 47 BY
genes annotate to mitochondrion (*Ank2, Csde1, Dap3, Huwe1, Nbr1, Ndufa5, Opa1,
Osgepl1, Ptcd3, Rps6kb1, Stxbp1, Tomm6, Trim27*). Both directions point the same
way as the manuscript's existing RNA-handling and mitochondrial reading, and
neither survives correction.

## 3. Cerebellar cell types — null, and markedly underpowered

`scripts/02_extract_kozareva_markers.py`, `scripts/05_celltype_enrichment.py`.

Markers obtained **at source**: Kozareva et al. 2021, *Nature* 598:214–219,
Supplementary Table 2 ("Highly differentially expressed genes across all
clusters"), sheet `Supplementary_Table2` of `41586_2021_3220_MOESM4_ESM.xlsx`
(sha256 recorded in `markers/kozareva_extraction_manifest.json`). The published
workbook is strict OOXML, which openpyxl silently reports as having zero sheets,
so it is parsed from the package XML; 36,945 marker rows extracted. No marker
was supplied from memory.

Marker rule: `padj ≤ 0.05`, `logFC > 0`, top 200 per subtype by logFC, 46 subtypes
collapsed to 17 classes, then intersected with the measured background. A fixed
per-subtype cap is used because the published table is itself truncated at 500
rows for many subtypes, so an "all significant markers" rule would hand larger
marker sets to whichever subtypes happen not to be truncated.

One-sided Fisher exact per class, BH across classes within each set; run twice,
once with all markers and once with class-unique markers only (102 tests total).

**0 significant at BH q ≤ 0.05.** Strongest signals:

| Set | class | overlap | expected | OR | p | BH q |
|---|---|---|---|---|---|---|
| GC−VIV top (class-unique) | Granule | 4/28 | 1.41 | 3.53 | 0.049 | 0.839 |
| GC−VIV top (all markers) | Granule | 7/79 | 3.99 | 2.00 | 0.100 | 1.00 |
| FLT−GC BY (all markers) | Bergmann glia | 2/44 | 1.04 | 2.45 | 0.281 | 1.00 |

The granule-cell lean in the **housing** set (*Cadm3, Ndrg3, Ppp3ca, Rif1*,
plus *Plcb4, Ptprd, Stxbp1*) is the only thing that even approaches nominal
significance, and it does not survive correction across 17 classes.

Marker sets surviving the background restriction are small — 2 (Ependymal) to 93
(OPC) genes — which is what destroys the power.

## 4. Why the null is uninformative about biology but informative about design

`scripts/06_power_and_overlap.py`.

For every class × set, the smallest overlap reaching uncorrected p ≤ 0.05 implies
a **minimum detectable odds ratio** of 2.2–95 (median 4.1 for the 47-gene BY set,
2.8 and 3.0 for the 106- and 100-gene sets). All 51 observed odds ratios fall
below their own floor. For functional terms the same floor runs from ~8.5-fold
for a 10-gene term in the 47-gene set to ~1.5-fold for a 500-gene term — so only
a very strong, very concentrated signal was ever visible.

**A real limit, honestly stated:** these sets cannot distinguish "no enrichment"
from "enrichment this design cannot see". Detecting a 1.5–2-fold cell-type lean
at 17 classes would need a query set in the several-hundred-gene range, which the
BY-controlled junction analysis does not produce at n = 9 per arm.

## 5. Flight and housing sets overlap heavily — and that is partly structural

40 of the 106 strongest FLT−GC genes are also among the 100 strongest GC−VIV
genes (5.4 expected if independent, OR 18.3); 21 of the 47 BY genes are
(2.4 expected, OR 19.0).

**This is not a clean test and is reported descriptively only.** FLT−GC and
GC−VIV share the ground-control arm, so a gene whose GC estimate is noisy enters
both rankings under the null; the two rankings are positively dependent by
construction. The p value is stored in
`tables/power_and_overlap_summary.json` flagged
`p_one_sided_NOT_A_VALID_TEST` with that caveat attached, and must not be read as
evidence of a shared biological programme. What it does support is the
manuscript's existing caution that housing hardware is not a passive comparator:
the housing contrast ranks many of the same genes.

---

## Draft Results paragraph

> Neither the flight-associated nor the housing-associated gene set is
> functionally or cell-type enriched once the measured gene universe is used as
> the background. Testing the 47 joint-BY FLT–GC genes, the 106 genes with the
> strongest FLT–GC exact evidence and the 100 strongest GC–VIV genes against all
> 1,980 measured genes, no Gene Ontology, KEGG or Reactome term reached
> significance under g:Profiler's g:SCS correction, and none reached
> Benjamini–Hochberg q ≤ 0.05 across the 6,682 terms of 5–500 genes
> (minimum q = 0.090, 0.205 and 0.444 respectively). The strongest uncorrected
> signals were translational initiation in the BY set (GO:0006413; Eif2a, Eif4a2,
> Eif4h and Rps6kb1 of 10 measured members; odds ratio 30.9, p = 5.2 × 10⁻⁵) and
> mitochondrion in the wider FLT–GC set (GO:0005739; 24 of 214, odds ratio 2.7,
> p = 1.8 × 10⁻⁴), consistent in direction with the candidate identities but not
> surviving correction. Over-representation of the markers of 17 cerebellar
> classes from a published single-cell atlas was likewise null: no class passed
> BH correction in either marker definition (102 tests), the strongest being
> granule-cell markers in the housing set (4 of 28 class-unique markers,
> odds ratio 3.5, uncorrected p = 0.049, q = 0.839). These nulls are
> power-limited rather than evidential. The smallest detectable odds ratio ranged
> from 2.2 to 95 across class tests, and every observed odds ratio fell below its
> own detectability floor (Supplementary Fig. 6b), so the analysis excludes only
> concentrated, large-magnitude enrichment. Separately, 40 of the 106 strongest
> FLT–GC genes were also among the 100 strongest GC–VIV genes; because the two
> contrasts share the ground-control arm these rankings are positively dependent
> by construction, so the overlap is reported as a description of the housing
> comparator rather than as a test of shared biology. Junction-level change is
> not expression change, and these tests ask only whether the gene identities are
> unusual among the genes this assay measured.

## Draft Methods paragraph

> **Functional and cell-type over-representation.** Gene sets were taken from the
> existing gene-level omnibus results without recomputation: genes passing joint
> Benjamini–Yekutieli correction in FLT–GC (n = 47), and the genes with the
> strongest exact evidence in FLT–GC and in GC–VIV. Because the exact p values are
> discrete, all genes tied with the 100th-ranked gene were retained, giving 106
> and 100 genes respectively. The background for every test was the measured gene
> universe — all 1,980 of the 1,983 analysed genes carrying an Ensembl symbol —
> never the mouse genome, since the junction pipeline measures only expressed,
> multi-junction genes and a genome background would manufacture neuronal and
> RNA-handling enrichment. Functional enrichment used the g:Profiler g:GOSt REST
> API (version e114_eg62_p19_9dc92ac8; Ensembl 114, GRCm39; GO classes
> 2026-01-23, KEGG 2026-03-15, Reactome 2026-03-20; queried 15 October 2026) with
> organism mmusculus, domain_scope 'custom' supplying the background explicitly,
> sources GO:BP, GO:CC, GO:MF, KEGG and Reactome, and the g:SCS multiple-testing
> method at threshold 0.05.<sup>kolberg2023,raudvere2019,reimand2007</sup>
> Because g:Profiler returns corrected p values capped at 1.0 when all results are
> requested, the uncorrected one-sided hypergeometric p and the
> Haldane–Anscombe-corrected odds ratio were recomputed from each term's 2 × 2
> contingency table, and Benjamini–Hochberg correction applied within each gene
> set over terms containing 5–500 genes, as a deliberately more permissive
> sensitivity analysis.<sup>benjamini1995</sup> Cerebellar cell-type markers were
> taken from Supplementary Table 2 of Kozareva et al., retaining genes with
> adjusted p ≤ 0.05 and positive log fold change and, within each of the 46
> annotated subtypes, the 200 with the largest log fold change; subtypes were
> collapsed to 17 classes.<sup>kozareva2021</sup> A fixed per-subtype cap was used
> because the published table is truncated at 500 rows for many subtypes. Marker
> sets and query sets were both intersected with the measured background before
> testing. Over-representation was assessed by one-sided Fisher exact test with
> Benjamini–Hochberg correction across classes within each gene set, repeated
> using only markers assigned to a single class. To make the null interpretable,
> the smallest overlap reaching an uncorrected p ≤ 0.05 and its implied odds ratio
> were computed for every class and term size. The overlap between the FLT–GC and
> GC–VIV gene sets is reported descriptively only, because the two contrasts share
> the ground-control arm and the rankings are therefore positively dependent under
> the null. Scripts, raw API responses and result tables are deposited with the
> analysis package.

## New references (all DOIs verified via Crossref, 2026-10-15; `raw/crossref_doi_verification.json`)

```bibtex
@article{kozareva2021transcriptomic,
  author  = {Kozareva, Velina and Martin, Caroline and Osorno, Tomas and
             Rudolph, Stephanie and Guo, Chong and Vanderburg, Charles and
             Nadaf, Naeem and Regev, Aviv and Regehr, Wade G. and Macosko, Evan},
  title   = {A transcriptomic atlas of mouse cerebellar cortex comprehensively
             defines cell types},
  journal = {Nature},
  volume  = {598},
  number  = {7879},
  pages   = {214--219},
  year    = {2021},
  doi     = {10.1038/s41586-021-03220-z}
}

@article{kolberg2023gprofiler,
  author  = {Kolberg, Liis and Raudvere, Uku and Kuzmin, Ivan and Adler, Priit
             and Vilo, Jaak and Peterson, Hedi},
  title   = {g:Profiler---interoperable web service for functional enrichment
             analysis and gene identifier mapping (2023 update)},
  journal = {Nucleic Acids Research},
  volume  = {51},
  number  = {W1},
  pages   = {W207--W212},
  year    = {2023},
  doi     = {10.1093/nar/gkad347}
}

@article{raudvere2019gprofiler,
  author  = {Raudvere, Uku and Kolberg, Liis and Kuzmin, Ivan and Arak, Tambet
             and Adler, Priit and Peterson, Hedi and Vilo, Jaak},
  title   = {g:Profiler: a web server for functional enrichment analysis and
             conversions of gene lists (2019 update)},
  journal = {Nucleic Acids Research},
  volume  = {47},
  number  = {W1},
  pages   = {W191--W198},
  year    = {2019},
  doi     = {10.1093/nar/gkz369}
}

@article{reimand2007gprofiler,
  author  = {Reimand, J{\"u}ri and Kull, Meelis and Peterson, Hedi and
             Hansen, Jaanus and Vilo, Jaak},
  title   = {g:Profiler---a web-based toolset for functional profiling of gene
             lists from large-scale experiments},
  journal = {Nucleic Acids Research},
  volume  = {35},
  number  = {suppl\_2},
  pages   = {W193--W200},
  year    = {2007},
  doi     = {10.1093/nar/gkm226}
}

@article{ashburner2000gene,
  author  = {Ashburner, Michael and Ball, Catherine A. and Blake, Judith A. and
             Botstein, David and Butler, Heather and Cherry, J. Michael and
             Davis, Allan P. and Dolinski, Kara and Dwight, Selina S. and
             Eppig, Janan T. and others},
  title   = {Gene Ontology: tool for the unification of biology},
  journal = {Nature Genetics},
  volume  = {25},
  number  = {1},
  pages   = {25--29},
  year    = {2000},
  doi     = {10.1038/75556}
}

@article{kanehisa2000kegg,
  author  = {Kanehisa, Minoru and Goto, Susumu},
  title   = {KEGG: Kyoto Encyclopedia of Genes and Genomes},
  journal = {Nucleic Acids Research},
  volume  = {28},
  number  = {1},
  pages   = {27--30},
  year    = {2000},
  doi     = {10.1093/nar/28.1.27}
}

@article{milacic2024reactome,
  author  = {Milacic, Marija and Beavers, Deidre and Conley, Patrick and
             Gong, Chuqiao and Gillespie, Marc and Griss, Johannes and
             Haw, Robin and Jassal, Bijay and Matthews, Lisa and
             May, Bruce and others},
  title   = {The Reactome Pathway Knowledgebase 2024},
  journal = {Nucleic Acids Research},
  volume  = {52},
  number  = {D1},
  pages   = {D672--D678},
  year    = {2024},
  doi     = {10.1093/nar/gkad1025}
}
```

Author lists above were taken from the Crossref records, not from recall, and
match them in order. One discrepancy to be aware of: the Crossref deposit for
10.1093/nar/28.1.27 lists only `Kanehisa, M.`, while the article itself is by
Kanehisa and Goto; the entry follows the article. Journal, volume, pages and
year for all seven entries are as returned by Crossref.

`benjamini1995controlling` (10.1111/j.2517-6161.1995.tb02031.x) and
`benjamini2001control` (10.1214/aos/1013699998) are already cited in the
manuscript; both were re-verified via Crossref in this pass.

---

## Files

```
research/I_enrichment/
  scripts/01_build_gene_sets.py              gene sets + measured background
  scripts/02_extract_kozareva_markers.py     strict-OOXML parse of the atlas tables
  scripts/03_gprofiler_enrichment.py         g:Profiler REST, custom background
  scripts/04_gprofiler_sensitivity_BH.py     exact recomputation + BH sensitivity
  scripts/05_celltype_enrichment.py          Fisher/BH cerebellar class tests
  scripts/06_power_and_overlap.py            minimum detectable effects, set overlap
  scripts/07_verify_dois_crossref.py         Crossref DOI verification
  scripts/08_figure_S6.py                    Figure S6 (183 mm, Arial >= 7 pt)
  raw/   gprofiler_*.json, gprofiler_run_metadata.json,
         gprofiler_data_versions_mmusculus.json, crossref_doi_verification.json
  markers/ Kozareva2021_SupplementaryTables_1-4.xlsx (+ SI PDF),
           Kozareva2021_Supplementary_Table1-4.csv, kozareva_extraction_manifest.json
  tables/  gene_sets_*, background_measured_genes.txt, set_*.txt,
           gprofiler_all_terms.csv.gz, gprofiler_recomputed_with_BH.csv.gz,
           gprofiler_top20_per_set.csv, gprofiler_significant_terms.csv (empty),
           celltype_enrichment_results.csv, celltype_enrichment_meta.json,
           celltype_minimum_detectable_effect.csv, cerebellar_marker_sets_used.csv,
           term_minimum_detectable_effect.csv, power_and_overlap_summary.json,
           headline_term_driving_genes.csv
  figures/ Figure_S6_enrichment.pdf/.png, Figure_S6_caption.txt
```

Nothing outside `research/I_enrichment/` was written; the submission package and
`src/PAPER_2_npj.tex` are untouched. Adding Supplementary Fig. 6 to the
manuscript is a separate decision, since the honest framing is a negative result
with a stated detectability floor.

## Unavailable in this session

MCP servers `fusion` (failed to start), `cloudflare-api` and `github` (awaiting
authorization) were not mounted. None was needed: g:Profiler, Springer Nature and
Crossref were reached over HTTPS directly.
