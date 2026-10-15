# D — Collection protocol versus spaceflight (run inline by the parent, 15 Oct 2026)

Worker D never started (host memory below the gateway's spawn floor), so the parent ran this analysis
directly. Script: `D_analysis.py` (seed 15102026; reuses the package's validated loaders and estimators
read-only). Runtime about 25 s.

## Verdict

The working hypothesis holds for the rapid primary-response (rPRG) panel, with one qualification.
In RRRM-2 ground mice that never flew, the ISS-T protocol exceeded the LAR protocol by 0.704 log2
(cerebellum; bootstrap 0.450 to 0.962; exact p = 0.005) and 1.118 (hippocampus; 0.666 to 1.539;
p = 0.010). Flight effects with collection matched (within LAR) were 0.188 (-0.146 to 0.519; p = 0.38)
and 0.281 (0.084 to 0.486; p = 0.060). The direct dominance contrast (ground ISS-T minus flight LAR) was
positive for rPRG in cerebellum (0.503; p = 0.010; q = 0.034), hippocampus (0.837; p = 0.015;
q = 0.042) and pooled (0.662; p = 0.0043; q = 0.034); 9 of 24 dominance tests passed BH.

Qualification: the IEG signal is not only a protocol signal. RR-23 habitat ground controls exceeded
vivarium controls by 0.641 log2 (q = 0.012) with agent, preservation and collection date matched, and the
RR-23 flight effect on rPRG flips sign with the control group (-0.399 vs HGC, q = 0.017; +0.242 vs VGC,
q = 0.21). The defensible claim is that IEG/HSP panels track collection protocol and housing at least as
strongly as spaceflight, not that spaceflight has no IEG effect. The hippocampal matched flight effect
(0.28, bootstrap interval excluding zero) may be a small real flight or return effect.

## What was done

1. Reproduction gate: all 80 package panel rows for the four RRRM-2 contrasts (two regions plus pooled)
   reproduced to 2.2e-16 (effects) and 1.1e-16 (exact p). `D_tables/D1_reproduction_vs_package.csv`.
2. Dominance test over the complete family of 8 panels x 3 analyses (24 tests, all reported; not preregistered). `D_tables/D2_*.csv`.
3. 2x2 decomposition with 2,000 whole-animal bootstraps stratified by flight x cohort x age; difference
   P_G - F_L and frequency P_G > F_L (1.000 for rPRG in every analysis). `D_tables/D3_*.csv`.
4. Gene-level effects for rPRG, dPRG, HSP and SRG genes: |P_G| > |F_L| in 18/19 cerebellar and 19/20
   hippocampal rPRG genes (42/56 and 53/57 across all four panels). `D_tables/D4_*.csv`.
5. Cross-mission classification of all 17 primary contrasts by deposited protocol factors (package panel
   family, BH over 232 tests). `D_tables/D5_*.csv`, `D5b_class_summary.csv`.

## Not done (from the original worker brief)

- Independent published gene signatures (Marsh 2022, van den Brink 2017, Hrvatin 2018) were not
  downloaded or tested. The panels used are literature-defined (Tyssowski 2018), so the gene-set
  circularity concern is limited, but an orthogonal-signature test remains open.
- OSD-352 (RR-3) was not added to the cross-mission forest: its flight and ground sacrifice dates differ
  (17/19 vs 22/23 May 2016) and RIN is fully separated, so it is not protocol-matched.

## Figure

`../../figures/Figure_2.{pdf,png,svg}` (built by `../../figures/build_figures.py fig2`).
