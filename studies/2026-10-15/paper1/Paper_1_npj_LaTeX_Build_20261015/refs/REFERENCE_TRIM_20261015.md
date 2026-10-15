# Reference trim (15 Oct 2026): 89 → 53 references

The authors asked for about 50 references, keeping the top and relevant ones. The result has 53:
- **42 peer-reviewed sources**, still above the earlier minimum of 40;
- **11 NASA OSDR dataset citations**, which the journal requires for every deposit used.

The pre-trim sources and PDF are in `archive/v3_before_reference_trim/`. Every removed key was checked to be no
longer cited. Where a sentence would have lost all of its support, it was reworded or re-supported as listed
below. The unused entries remain in `paper1.bib`; BibTeX prints only cited entries.

## Removed (36)

| Reason | References |
|---|---|
| Redundant background (other cited sources carry the claim) | Afshinnekoo 2020 (Cell), Overbey 2024 SOMA (Nature), Yap 2018 (Neuron), Leek 2014 svaseq, Zhang 2020 ComBat-seq, Molania 2023 RUV-III, Gerard 2021 RUV*, Gerstner 2016, Zindler 2020, Varma 2006, Kapoor 2023, Tullai 2007 |
| Lower-tier or less relevant to brain transcriptomes | Mao 2020 (FASEB J; same cohort as Holley 2022), Overmyer 2015 (PLOS ONE; metabolomics, no brain), Chakkingal Bhaskaran 2023 (BMC Genomics; pig), Hamaya 2000, Davoudian 2023 (ACS Chem Neurosci), Nakao 1993, Dachet 2021, Ferreira 2018, Lacar 2016, Overbey 2022 (scRNA collection perspective), Choi 2020 (replaced by Lai Polo 2020 and Grandke 2026, which describe the same two collection routes) |
| Standard methods that need no citation | Hoerl 1970 (ridge), Wold 2001 (PLS), Efron 1979 (bootstrap), Winkler 2014 and Phipson 2010 (permutation), Beck 2009 and O'Donoghue 2015 (accelerated gradient), Frank 1956 (Frank–Wolfe; Jaggi 2013 kept), Morris 2019 (simulation reporting), Button 2013 (low power) |
| Tangential "related construction" or recommendation support | Izenman 1975, Gouk 2021 (that Methods clause was removed), Ye 2023 (remeasurement; the recommendation is kept without a citation) |

## Wording changes forced by the trim

- **Introduction:**
  - The list of brain findings now reads "oxidative-stress pathways"; "blood–brain barrier proteins" was dropped
    with Mao 2020.
  - The IEG sentence now reads "respond to sensory experience … tissue dissociation and ex vivo handling";
    "novelty" and "post-mortem delay" were dropped with Lacar, Dachet and Ferreira.
  - The sentence on collection strategy (Overmyer, Chakkingal Bhaskaran) was deleted; Lai Polo 2020 already carries
    the point.
- **Discussion:**
  - The pig-tissue clause was deleted.
  - The ketamine sentence now cites Hu 2024 only.
- **Methods:** the clause naming reduced-rank and Lipschitz-constrained regression was deleted.

## Kept despite a mid-tier journal (each is the only source for its claim)

- Kremsky 2023 (IJMS): IEGs among spaceflight-associated brain genes.
- Staib-Lasarzik 2014 (J Neurotrauma): anaesthesia for euthanasia changes brain mRNA.
- Bunting 2016 (Front Behav Neurosci): isoflurane blocks IEG initiation, an alternative reading of the
  protocol effect.
- Lai Polo 2020 and Overbey 2021 (iScience): the closest prior art, and the GeneLab pipeline that produced the
  counts.
- Roy 2025 (Cells): spaceflight brain inflammasome study by the corresponding author.

The five signature papers that the pending published-signature paragraph needs (Hrvatin, Marsh, van den Brink, Wu,
Tyssowski) are all kept.
