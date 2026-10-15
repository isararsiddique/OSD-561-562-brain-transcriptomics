# Design-selected external mouse-brain data protocol

15 October 2026, Asia/Kolkata. This acquisition defines eligible external cohorts from organism, tissue, assay, exposure/control design, source identities and availability of unnormalized counts. It does not select cohorts, genes or samples by the strength or direction of expression results. Subsequent expression analyses remain exploratory; this document is not a preregistration.

## Selection and retained evidence

Four available adult mouse-brain bulk RNA-seq datasets were selected: OSD-202, OSD-352, OSD-374 and OSD-611. Original RRRM-2/RR-10/RR-23 deposits are not counted as external cohorts. Single-cell/spatial observations from the same animals are not additional animal replicates. The acquired matrices contain **66 columns**, of which **64 are primary eligible analysis units**: 40 + 6 + 8 + 10. These units combine individual animals and pooled animals and must not be reported as 64 independent mice.

| Accession | Design and primary units | Relevant interpretation limits |
| --- | --- | --- |
| [OSD-202](https://osdr.nasa.gov/bio/repo/data/studies/OSD-202) | Earth brain; 40 female C57BL/6J animals, recorded age six months. Radiation 0/0.04 Gy cobalt-57 gamma × normally loaded/hindlimb unloaded × one/four months after treatment. Twenty animals at each collection time; eight cells have 2–6 animals. | This is an independent factorial Earth experiment, not a flight replication. Radiation/unloading lasted 21 days according to the study description. Collection time, elapsed ageing and collection context cannot be separated. Three sequencing lanes partly align with design cells. |
| [OSD-352](https://osdr.nasa.gov/bio/repo/data/studies/OSD-352) | RR-3 Brain; six female BALB/c animals, recorded age at launch 12 weeks: three flight and three habitat ground controls. Duration records are 39–42 days. STAR counts already have six animal-level columns. | Independent of the original missions, but the same RR-3 animals supplied published spatial/multiome material. The 18 ISA sample records include spatial technical sections; do **not** sum already collapsed bulk counts again. Brain is its own anatomical context, not interchangeable with hippocampus/cerebellum. |
| [OSD-374](https://osdr.nasa.gov/bio/repo/data/studies/OSD-374) | Earth hippocampus; male C57BL/6J, six months; four sham and four 1 Gy proton pools, sampled two weeks after irradiation. Each of eight libraries pools two animals; 16 source-animal labels are mapped to eight disjoint libraries. | Inference uses eight pools. Assay records repeat each library for its two member samples; they are not 16 expression replicates. This is radiation context, not independent spaceflight confirmation. |
| [OSD-611](https://osdr.nasa.gov/bio/repo/data/studies/OSD-611) | Earth hippocampus; male C57BL/6J, six months; two weeks after silicon-28 irradiation. Twelve STAR columns were acquired. The NASA SampleTable contains ten: three sham, five 0.3 Gy, two 0.6 Gy pools. | The source collection protocol states two hippocampi per RNA-seq library, but ISA supplies one generic label per pool and does not identify both members. Ten primary libraries are ten pooled experimental units, not ten individually identified mice or 20 replicates. The two extra STAR columns remain available and flagged, without inventing a reason for their SampleTable absence. |

For OSD-611, `primary_analysis_included` and `present_in_nasa_sampletable` identify the ten primary columns. `counts.csv.gz` deliberately retains all twelve acquired columns. Analysis code must apply the flag or explicitly label a complete twelve-column sensitivity. The absent SampleTable columns are Mouse_hippocampus_2 (sham) and Mouse_hippocampus_5 (0.6 Gy). Their exclusion is an availability/design rule, not a demonstrated QC failure or a result-based exclusion.

OSD-687 was checked directly through NASA APIs and rejected: it profiles histone modifications in Inspiration4 human blood T cells, rather than mouse-brain RNA. Other locally catalogued deposits with kidney tissue (OSD-771), embryonic hypothalamic cell culture (OSD-792), human neural organoids (OSD-863/871), microarrays, miRNA, protein, behaviour or spatial-only assays do not meet the selected adult mouse-brain bulk-mRNA count criterion. This is a bounded audit of available candidates, not a claim that every possible public brain study has been exhaustively screened.

## Identity, metadata and numerical processing

`acquire_external_data.py` copies exact NASA processed STAR **unnormalized integer gene counts**, ISA ZIPs, NASA SampleTables and study/file metadata. It uses no normalized counts, VST matrices or prior differential-expression discoveries. The clean CSV.gz changes storage/index labeling only: no gene filtering, normalization, imputation, technical summation or coefficient estimation occurs during acquisition. It checks unique row/column identifiers and finite, nonnegative integer values. “Unnormalized counts” are processed alignment/count outputs; raw FASTQ alignment was not repeated.

ISA Sample Name maps to count columns by whitespace canonicalization and deposited merged-file naming. Source Name determines individual animals when available. Repeated ISA `Unit` headers are read by their position immediately after the associated field, preserving age, collection-time and dose units. Dose source values/units and their Gy conversion remain separate columns. Empty metadata fields remain missing rather than being inferred. Names such as Mouse 1 are local study/pool labels, not cross-study RFID identities; accession-qualified identifiers prevent accidental merging of unrelated labels.

Each cohort has `counts.csv.gz`, `metadata.csv`, preserved sample/assay records and full `study_protocols.json`. Metadata include source identities, pool status, age/units, sex, strain, anatomy, exposure/control labels, collection-time units, radiation dose units, preservation, available euthanasia information, library batch, sequencing lane/instrument, raw-file names and available RIN. `sample_to_source_membership.csv` records mapped members; OSD-611 members remain explicitly unknown.

## Collection and quality limitations

OSD-352 flight RIN values are 6.3–6.8, whereas ground values are 7.5–8.1, with no overlap. All six bulk libraries share LIB-BATCH-021, but source sacrifice dates differ between flight and ground groups. Frozen carcasses were thawed at room temperature for 15–20 minutes before dissection; the deposited collection protocol also contains a generic heart-tissue sentence. Preserve the original wording and use the actual assay/sample Brain designation. A same-library-batch result does not eliminate preservation or collection confounding. RIN adjustment here requires extrapolation and cannot establish an unconfounded flight effect; the RIN separation does not itself prove every observed effect is artifactual.

OSD-202 available RIN values range 7.0–8.8. Sequencing-lane sensitivities must retain cells with sparse support and report unavailable uncertainty when the adjusted design cannot identify a contrast. OSD-374/611 report liquid-nitrogen preservation and cervical dislocation followed by decapitation; RIN is unavailable in the harmonized records. Differences in age, sex, strain, assay selection and brain material prevent simple pooling with original flight cohorts.

External exposure-associated scores may test prediction/transport behavior, but they are **not established nuisance-only scores**. Radiation and flight alter biology. Improved held-out expression prediction does not show that subtracting the predicted component recovers the true biological exposure effect.

## Provenance and reproduction

`acquisition_provenance.json` records each local source path, NASA URL, expected NASA byte size, copied-file byte size, SHA256 and acquisition mode. `output_manifest.json` hashes derived acquisition outputs. `cohort_inventory.csv` describes acquired matrices; primary flags determine the smaller analysis set. Local source files were read and copied without modifying the original repository/package. OSD-687 API responses preserve the rejection evidence. Dataset citations should credit the original investigators/NASA OSDR rather than imply new animal experiments.

To reproduce acquisition from the same local cache:

```sh
PYTHONPATH=/tmp/p1_completion_py /Users/zain/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 acquire_external_data.py
```

Run from this directory, or pass the absolute script path. `--source-root` selects an equivalent NASA cache, `--output` selects an acquisition output folder, and `--online` permits retrieval when a source file is missing. Scientific inference, selection of genes and figure generation occur in the separate benchmark/biological-analysis extensions.

Final manuscript formatting must follow the verified npj Microgravity requirements cached in the earlier `extensions/journal_audit_20261015/JOURNAL_AND_NOVELTY_AUDIT.md`: title ≤15 words, unstructured abstract ≤150 words, all Methods in the main file, no separate Funding statement, readable RGB figures ≥300 dpi, defined experimental n and complete distributions/uncertainty. No graphical abstract is mandated by the inspected guidance. Public data/code archive identifiers must be real; none is invented here.
