# Spatial Transcriptomics of the Mouse Brain under Spaceflight and Antioxidant Treatment

Re-analysis of NASA OSDR accessions **OSD-682, OSD-685, OSD-698 and OSD-699**,
centred on a pre-specified CNS-injury / extracellular-vesicle (EV) mRNA target panel.

> **The analysis lives in
> [`OSD-682-685-698-699-brain-spatial/`](OSD-682-685-698-699-brain-spatial) —
> start with [its README](OSD-682-685-698-699-brain-spatial/README.md).**

The four accessions are not four experiments. They are four brain regions of one
NanoString GeoMx Digital Spatial Profiling run (GEO
[GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336)): CA1,
dentate gyrus, frontal cortex and cerebral cortex, in a 2 × 2 factorial of
spaceflight and the antioxidant BuOE — 48 ROIs, 15,782 targets, n = 3 per cell.

## Headline findings

**The differential-expression statistics deposited with the study cannot be used
as they stand.** The `Adjusted pvalue` column is smaller than its own raw
`Pvalue` for 77% of targets, which no multiple-testing correction can produce,
and it is rounded to two decimals. Read as deposited, 1,317–2,361 targets look
significant per contrast; correct Benjamini-Hochberg control yields **zero**
across the entire experiment.

**Transcriptome-wide, the design has almost no power — and that is quantifiable.**
At n = 3 a target would need a **6–8-fold change** to survive correction across
15,782 targets. The empty genome-wide result reflects the design, not the biology.

**Restricting to the pre-specified target panel is what makes the data informative.**
Testing the 111 measurable panel members instead of 15,782 targets cuts
multiplicity ~140-fold, and eight results then reach BH < 0.05 — **JUNB**,
**MMP12** and **neurogranin**. All have genome-wide adjusted p-values of
0.32–0.87, so none is visible without the restriction.

**Six Ensembl IDs in the source target spreadsheet point to unrelated genes** and
were corrected against Ensembl (Iba1 → NMRAL1, Neurogranin → DIDO1, SBDP → DCX,
SAA → FLCN, plus two apparent transposed-digit typos). One of the three
significant genes is neurogranin, which the supplied ID would have missed.

Full detail, figures, caveats and limitations:
[`OSD-682-685-698-699-brain-spatial/README.md`](OSD-682-685-698-699-brain-spatial/README.md).

## Layout

```
.
├── CNS mRNA targets & Others in EV ... .csv   # source target panel (150 analytes)
└── OSD-682-685-698-699-brain-spatial/         # the analysis
    ├── README.md                              # methods, results, limitations
    ├── Makefile                               # make setup && make all
    ├── scripts/                                # 13 pipeline stages
    ├── config/                                 # target mapping, sign convention
    └── results/                                # figures (PDF + PNG) and tables
```

## Quick start

```bash
cd OSD-682-685-698-699-brain-spatial
make setup     # uv venv on Python 3.12 + pinned requirements
make all       # data -> conventions -> analyze -> panels -> tables -> summary
```

Input data (~12 MB) is downloaded on demand and is not versioned here.

## The OSD-561/562 bulk RNA-seq project

That project previously occupied this repository root. It is preserved on the
remote and is not lost:

| Ref | Contents |
|---|---|
| `backup/osd-561-562-main` branch | full project as it stood, 112 files |
| `osd-561-562-main-backup` tag | same commit, `81f83ef` |

To recover it:

```bash
git checkout backup/osd-561-562-main
```

## Data availability

Inputs are openly available from the NASA Open Science Data Repository
([OSD-682](https://osdr.nasa.gov/bio/repo/data/studies/OSD-682),
[OSD-685](https://osdr.nasa.gov/bio/repo/data/studies/OSD-685),
[OSD-698](https://osdr.nasa.gov/bio/repo/data/studies/OSD-698),
[OSD-699](https://osdr.nasa.gov/bio/repo/data/studies/OSD-699)) and from NCBI GEO
([GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336)).
Ortholog mapping uses the [Ensembl REST API](https://rest.ensembl.org). No raw
data is redistributed; the download scripts fetch it on demand. Please cite OSDR,
GEO and the original data contributors when reusing.

## Authors

**Isarar Siddique** — [@isararsiddique](https://github.com/isararsiddique) ·
<siddisrar786@gmail.com>

## License

Not yet set. The analysis code may be released under a permissive license of your
choice; the underlying OSDR and GEO data remain subject to their own terms.
