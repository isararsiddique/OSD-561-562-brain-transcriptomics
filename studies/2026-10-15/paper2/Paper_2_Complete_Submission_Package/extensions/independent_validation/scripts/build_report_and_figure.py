"""Render the four-panel scientific figure and report from checked tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = Path(__file__).resolve().parents[1]
T, F = B / "tables", B / "figures"


def main():
    F.mkdir(exist_ok=True)
    s = json.loads((T / "summary.json").read_text())
    checks = json.loads((T / "independent_checks.json").read_text())
    assert checks["all_checks_passed"]
    event = pd.read_csv(T / "complete_common_panel_effects.csv.gz")
    boot = pd.read_csv(T / "animal_bootstrap_draws.csv.gz")
    animal = pd.read_csv(T / "RRRM2_heldout_animal_scores.csv")
    strata = pd.read_csv(T / "descriptive_stratum_score_effects.csv")
    sens = pd.read_csv(T / "measurement_sensitivities.csv")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 12, "axes.labelsize": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "savefig.facecolor": "white"})
    teal, orange, ink, grey = "#007E87", "#CC6D29", "#243441", "#889BA8"
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 10.5))
    fig.subplots_adjust(top=.87, bottom=.18, left=.16, right=.97, hspace=.62, wspace=.38)
    fig.suptitle("The complete RR23 junction profile does not show aggregate transfer to RRRM2", x=.53, y=.976,
                 fontsize=15, weight="bold", color=ink)
    fig.text(.53, .939, "Independent launch · cerebellum · fixed RR23 directions · 19 held-out animals in three route–age strata",
             ha="center", fontsize=10.5, color=ink)
    ax = axes[0, 0]
    counts = [s["frozen_definitions"], s["RR23_measured_all26"], s["common_measured_events"]]
    bars = ax.barh([2, 1, 0], counts, color=[grey, teal, orange], height=.58)
    ax.set_yticks([2, 1, 0], ["RR23 all26\ndefinitions", "Measured in all\n26 RR23 animals", "Measured in both\nRR23 and RRRM2"])
    ax.set_xlim(0, 6600)
    ax.set_xlabel("Junction-pair events")
    ax.set_title("A  Measurement transport limits the panel", loc="left", weight="bold", color=ink, pad=14)
    for bar, value in zip(bars, counts):
        ax.text(value + 95, bar.get_y() + bar.get_height() / 2, f"{value:,}", va="center", weight="bold", color=ink)
    ax.text(.0, -.24, f"{100*s['common_transport_fraction']:.1f}% of the RR23 measured panel is measurable in both cohorts.\nRRRM2 failures: 2,436 coverage only; 75 dominance only; 1,226 both.",
            transform=ax.transAxes, fontsize=9.0, color=ink, va="top")
    ax = axes[0, 1]
    a = 100 * event.RR23_delta_FLT_minus_GC
    b = 100 * event.RRRM2_delta_stratified_FLT_minus_GC
    ax.scatter(a, b, s=10, color=teal, alpha=.40, edgecolors="none")
    limit = max(abs(a).max(), abs(b).max()) * 1.07
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.axhline(0, color=grey, lw=.8)
    ax.axvline(0, color=grey, lw=.8)
    ax.plot([-limit, limit], [-limit, limit], "--", color=grey, lw=.8)
    ax.set_xlabel("RR23 FLT − GC fraction (percentage points)")
    ax.set_ylabel("RRRM2 stratified FLT − GC (percentage points)")
    ax.set_title(f"B  Complete common panel: {len(event):,} events", loc="left", weight="bold", color=ink, pad=14)
    ax.text(.04, .95, f"Spearman ρ = {s['spearman_rho']:.3f}\nPearson r = {s['pearson_r']:.3f}\nSame direction = {100*s['same_direction_fraction']:.1f}%", transform=ax.transAxes,
            va="top", fontsize=10, bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
    ax = axes[1, 0]
    ci = s["spearman_animal_bootstrap_CI95"]
    ax.hist(boot.spearman_rho, bins=37, density=True, color=teal, alpha=.85, edgecolor="white", linewidth=.45)
    ax.axvline(0, color=ink, lw=1, linestyle="--")
    ax.axvline(s["spearman_rho"], color=orange, lw=2)
    for q in ci:
        ax.axvline(q, color=orange, linestyle=":", lw=1.5)
    ax.set_xlabel("Common-panel Spearman correlation")
    ax.set_ylabel("Animal-bootstrap density")
    ax.set_title("C  Animal resampling leaves correlation near zero", loc="left", weight="bold", color=ink, pad=14)
    ax.text(.04, .95, f"2,000 draws within exposure × route–age cells\nPercentile 95% CI [{ci[0]:.3f}, {ci[1]:.3f}]", transform=ax.transAxes,
            va="top", fontsize=9.5, bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
    ax = axes[1, 1]
    order = ["ISS_T_12", "LAR_12", "LAR_29"]
    rng = np.random.default_rng(2015)
    for pos, st in enumerate(order):
        sub = animal.loc[animal.stratum.eq(st)]
        # Common within-stratum center uses no exposure label and changes no contrast.
        center = sub.learned_direction_score.mean()
        for exp, color, offset in [("GC", grey, -.17), ("FLT", orange, .17)]:
            group = sub.loc[sub.exposure.eq(exp)]
            y = 100 * (group.learned_direction_score - center)
            x = pos + offset + rng.uniform(-.035, .035, len(group))
            ax.scatter(x, y, s=45, color=color, edgecolors="white", lw=.5, label=exp if pos == 0 else None, zorder=3)
            ax.plot([pos + offset - .07, pos + offset + .07], [y.mean(), y.mean()], color=ink, lw=1.6, zorder=4)
    ax.axhline(0, color=grey, lw=.8)
    ax.set_xticks(range(3), ["ISS-T · 12 wk\n5 FLT / 3 GC", "LAR · 12 wk\n2 FLT / 3 GC", "LAR · 29 wk\n3 FLT / 3 GC"])
    ax.set_ylabel("Fixed direction score (pp)\ncentered within stratum")
    ax.set_title("D  Held-out animal score does not validate", loc="left", weight="bold", color=ink, pad=14)
    ax.legend(frameon=False, loc="upper right", ncol=2, fontsize=9)
    ax.text(.01, -.30, f"Stratified score Δ = {s['learned_direction_score_delta_percentage_points']:.3f} pp; exact one-sided p = {s['learned_direction_score_exact_one_sided_p']:.3f}\n11,200 labels preserving FLT counts within every stratum", transform=ax.transAxes,
            va="top", fontsize=9, color=ink)
    fig.text(.16, .029, "All eligible events used; no RR23 significance filter. Measurement: ≥20 pair reads and ≥90% boundary dominance in every animal.\nRR23 flight and collection day are aliased. Results describe cohort associations and do not establish causal flight effects or functional injury.",
             fontsize=9.0, color=ink, va="bottom")
    fig.savefig(F / "independent_launch_validation_four_panel.png", dpi=300)
    fig.savefig(F / "independent_launch_validation_four_panel.svg")
    fig.savefig(F / "Figure_Independent_Validation.png", dpi=300)
    fig.savefig(F / "Figure_Independent_Validation.svg")
    fig.savefig(F / "Figure_Independent_Validation.pdf")
    plt.close(fig)
    caption = ("Independent-launch validation of the complete RR23 all26-defined junction panel in RRRM2/OSD-561 cerebellum. "
               "(A) Measurement transport: 5,724 frozen definitions, 5,249 measured across all 26 RR23 animals, and 1,512 meeting the same "
               "all-animal coverage and dominance gates in both studies (28.8% of the RR23 measured panel). RRRM2 failure counts are disjoint. "
               "(B) Every common-panel event effect; RR23 FLT−GC means and the RRRM2 animal-count-weighted average of FLT−GC means within "
               "route×age strata. Spearman ρ=0.010 and Pearson r=−0.011; the dashed diagonal marks equal effects. "
               "(C) Spearman correlation from 2,000 independent whole-animal bootstrap draws within exposure groups in RR23 and "
               "exposure×route×age cells in RRRM2; orange lines show the observed correlation and percentile 95% interval [−0.069, 0.074], "
               "and the black dashed line marks zero. Events are not treated as independent biological replicates. "
               "(D) RR23-only learned directions projected with equal event weights into every held-out RRRM2 animal. Points are animals "
               "and bars are group means; a common stratum mean is subtracted only for display. The weighted score contrast is 0.128 "
               "percentage points (animal-bootstrap 95% interval −0.137 to 0.422), exact one-sided p=0.2525 over all 11,200 "
               "assignments preserving FLT counts within each route×age stratum (two-sided p=0.5022). "
               "No RR23 significance, large-effect or direction-positive subset was selected. The measured common panel shows no "
               "aggregate transfer. Measurement requires at least 20 pair reads and at least 90% dominance among all retained "
               "annotated canonical alternatives sharing the boundary in every animal. RR23 flight and collection day are aliased, "
               "so this comparison cannot establish causal flight effects or functional injury.\n")
    (F / "Figure_Independent_Validation_caption.txt").write_text(caption)
    score_ci = np.array(s["learned_direction_score_animal_bootstrap_CI95"]) * 100
    stratum_lines = "\n".join(f"| {row.stratum} | {row.n_FLT}/{row.n_GC} | {100*row.score_delta:.3f} | [{100*row.CI95_low:.3f}, {100*row.CI95_high:.3f}] |" for row in strata.itertuples())
    sens_lines = "\n".join(f"| {row.panel} | {row.n_events:,} | {row.spearman_rho:.4f} | {100*row.score_delta:.3f} | {row.exact_one_sided_p:.4f} |" for row in sens.itertuples())
    report = f"""# Independent-launch validation of the RR23 all26 junction panel

Author: Isarar Siddique  
Affiliation:  
Email:  
Analysis date: 2026-10-15

## Conclusion

The complete RR23 FLT–GC junction-fraction profile did not show aggregate transfer to independently launched RRRM2/OSD-561 cerebellum under the dated retrospective plan. Among {s['common_measured_events']:,} events measurable in both cohorts, Spearman ρ was {s['spearman_rho']:.4f} (2,000 animal-bootstrap percentile 95% CI {ci[0]:.4f} to {ci[1]:.4f}); Pearson r was {s['pearson_r']:.4f}. The prespecified held-out learned-direction score gave a stratified FLT–GC contrast of {s['learned_direction_score_delta_percentage_points']:.3f} percentage points (95% CI {score_ci[0]:.3f} to {score_ci[1]:.3f}), exact one-sided p = {s['learned_direction_score_exact_one_sided_p']:.4f} and two-sided p = {s['learned_direction_score_exact_two_sided_p']:.4f}. This supports a clear limit on cohort-level transport, rather than a generalized conserved flight signature.

This is a new calculation for the all26-defined panel. The historical different-panel RRRM2 result was not imported or substituted. No significant, large-effect, direction-positive, or annotated-gene subset of RR23 was selected. A null result does not prove absence of a biological response or absence of all event-specific replication. RR23 flight status is aliased with collection day; no causal flight or functional injury claim is identified by these data.

## Design and measurement

The RR23 feature definition used all 26 biological animals (9 FLT, 9 GC, 8 VIV). Its learned contrast used the 9 FLT and 9 GC animals. The frozen 5,724 event definitions are ordered pairs of annotated canonical splice junctions sharing a boundary. Fractions are junction_a / (junction_a + junction_b); they are local junction fractions, not complete isoform abundances. The source STAR junction counts have technical libraries summed before analysis, and every retained column maps to one distinct animal.

RRRM2 contains 26 animal metadata records. The entire ISS_T × 29-week stratum was excluded under the existing rule because its four FLT records have conflicting deposited versus source-inferred euthanasia location; the three corresponding GC animals were excluded with that unusable stratum. The included 19 animals comprise 10 FLT and 9 GC across ISS_T × 12-week (5/3), LAR × 12-week (2/3), and LAR × 29-week (3/3) strata. The effect averages within-stratum FLT–GC contrasts with fixed weights 8/19, 5/19, and 6/19. These labels retain their original cohort identifiers; they do not resolve all collection and processing differences.

Measurement required ≥20 reads across the fixed pair and ≥0.90 pair dominance among all retained annotated canonical junctions sharing its boundary, in every included animal. RR23 measurements include all 26 animals. RRRM2 boundary alternatives were independently reconstructed from its complete animal-level junction-count file.

- Frozen RR23 definitions: {s['frozen_definitions']:,}.
- Measured across all RR23 animals: {s['RR23_measured_all26']:,}.
- Measured in RRRM2 among the frozen definitions: {s['RRRM2_measured_definitions']:,}.
- Common primary panel: {s['common_measured_events']:,}, {100*s['common_transport_fraction']:.2f}% of the RR23 measured panel.
- RR23 measured events failing RRRM2 coverage: {s['RR23_measured_RRRM2_coverage_fail']:,}; failing dominance: {s['RR23_measured_RRRM2_dominance_fail']:,}; failing both: {s['RR23_measured_RRRM2_both_fail']:,}. These failure counts overlap. The disjoint failure counts are 2,436 coverage only, 75 dominance only, and 1,226 both.

The shared panel is determined by outcome-independent measurement thresholds, not inferred significance. Poorly measured features were excluded from transfer estimation rather than interpreted as biologically unchanged. Requiring all animals to pass restricts conclusions to the measured panel and can disproportionately retain abundant or stable features.

## Held-out score and inferential unit

For each common event, the sign of the RR23 FLT–GC mean difference was learned once. Each held-out animal score is the mean signed difference between its fraction and the RR23 GC mean across every common event. Directions and centering use RR23 only; RRRM2 exposure labels do not train the score. Events are equally weighted and their correlation is preserved within each animal. Subtracting a common reference does not change a FLT–GC contrast.

Every one of the 11,200 assignments preserving FLT counts within route–age strata was enumerated. The observed assignment is included. {checks['exact_null_one_sided_extreme_count']:,} assignments were as or more positive than the observed contrast; {checks['exact_null_two_sided_extreme_count']:,} were as or more extreme in absolute value. The primary expected-direction test was one-sided at alpha 0.05. This is an animal-level permutation inference under conditional exchangeability within strata; labels were not physically randomized by this retrospective analysis, and residual collection/processing confounding remains possible.

Profile intervals resampled whole animals within RR23 exposure groups and RRRM2 exposure × route–age cells, keeping panel and stratum weights fixed (2,000 draws, seed {s['bootstrap_seed']}). RR23 training effects were re-estimated in the profile bootstrap. The held-out score interval keeps the primary RR23 directions fixed and resamples RRRM2 animals only. No event-independence correlation p values, gene-level independent-replicate p values, read-level binomial replication, or technical-replicate inference were used. Same-direction concordance was {100*s['same_direction_fraction']:.2f}% (animal-bootstrap 95% CI {100*s['same_direction_animal_bootstrap_CI95'][0]:.2f}% to {100*s['same_direction_animal_bootstrap_CI95'][1]:.2f}%).

## Secondary results

The entire common panel also received held-out exact directional event p values with directions frozen from RR23 and Holm correction across all {s['common_measured_events']:,} common events. Zero events passed Holm-adjusted p < 0.05. The smallest unadjusted directional event p was {s['minimum_event_directional_p']:.6f}; it is not a corrected replication claim.

All route–age score contrasts are shown descriptively; none was selected to rescue the aggregate result. Small cell sizes make these percentile intervals unstable and they should not be interpreted as independently powered stratum validation.

| Route–age stratum | FLT/GC | Score difference (pp) | Animal-bootstrap 95% CI (pp) |
|---|---:|---:|---:|
{stratum_lines}

The two planned label-invariant measurement sensitivities also do not replace the primary result. Their one-sided p values are additionally Holm-adjusted together in the saved sensitivity table.

| Panel | Common events | Spearman ρ | Score difference (pp) | Exact one-sided p |
|---|---:|---:|---:|---:|
{sens_lines}

## Interpretation and limitations

Independent launch helps test generality across missions but is not biological matching. RR23 metadata describe male C57BL/6J mice aged 16–17 weeks at launch, while RRRM2 metadata describe female C57BL/6NTac mice at 12 or 29 weeks, with distinct recovery routes, tissue sampling and preservation. Any loss of transfer can combine mission, sex, strain, age, collection, processing, and measurement effects; this analysis cannot assign the cause. Conversely, a positive score could not have established a flight-specific causal effect.

Only 28.8% of the RR23 measured panel can be evaluated under the same all-animal thresholds in RRRM2. The near-zero profile correlation has an interval spanning modest positive and negative associations. The primary null therefore establishes neither exact equivalence to zero nor a universal absence of response. The scientifically supportable contribution is a rigorously bounded, retrospectively planned test of junction-profile transportability, including its measurement limitations. Any paper should present this null validation prominently and should not claim a conserved cerebellar injury mechanism.

## Verification and outputs

Independent verification passed all checks. A separate script did not import the analysis code: it read raw gzip CSV counts, summed boundary alternatives with dictionaries, reproduced all 5,724 RR23 fractions and all 5,724 RRRM2 fractions/pair counts/dominance rows, recalculated every held-out animal score and common event contrast, independently enumerated the 11,200 score assignments, verified 13 deterministic event nulls and all Holm values, and checked saved bootstrap interval quantiles. SHA256 hashes of all source inputs and the dated plan are unchanged.

A standalone copied `beta/extensions/independent_validation` layout reproduced the complete summary exactly with `--base ../..`, using only the five required beta files and frozen extension inputs. All independent checks passed in that copied layout; the result is saved in `tables/portable_layout_check.json`. The temporary test copy was removed.

The figure has four panels: measurement transport; all common-panel event effects; correlation bootstrap uncertainty; and individual held-out scores grouped by route–age. Panel D subtracts a common mean within each stratum for display only; all inferential calculations use the original fixed score. The mean bars are descriptive. The complete feature and animal tables preserve every eligible feature and animal.

Main artifacts are `figures/independent_launch_validation_four_panel.png` (300 dpi) and `.svg`; `tables/summary.json`; `tables/complete_feature_mapping_and_gates.csv.gz`; `tables/complete_common_panel_effects.csv.gz`; `tables/RRRM2_heldout_animal_scores.csv`; `tables/RRRM2_exact_score_null.csv.gz`; `tables/animal_bootstrap_draws.csv.gz`; `tables/measurement_sensitivities.csv`; `tables/input_manifest_sha256.csv`; and `tables/independent_checks.json`.

## Reproduce

Run the three scripts in this order with Python 3.11+ and numpy, pandas, scipy, and matplotlib, and `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`. Paths below are relative to this extension. For the first two scripts, supply `--base /path/to/beta` pointing to the beta analysis/package root containing `inputs/` and `results/`:

1. `scripts/analyze_validation.py --base /path/to/beta`
2. `scripts/verify_independently.py --base /path/to/beta`
3. `scripts/build_report_and_figure.py`

Input paths and byte-level SHA256 checksums are recorded with portable beta/extension scopes in `tables/input_manifest_sha256.csv`. Frozen RRRM2 count and metadata inputs, publisher inventory and MD5 checksum file, original file URLs and hash verification, and independent retained-STAR-to-collapsed-matrix verification are saved inside `inputs/`. The raw STAR files are not needed for reproduction from the checked collapsed matrix. The analysis reads frozen source files without modifying them. Its frozen plan is `RETROSPECTIVE_ANALYSIS_PLAN.md`; prospective preregistration is not claimed.
"""
    (B / "INDEPENDENT_VALIDATION_REPORT.md").write_text(report)
    print(json.dumps({"figure": str(F / "independent_launch_validation_four_panel.png"),
                      "report": str(B / "INDEPENDENT_VALIDATION_REPORT.md"),
                      "independent_checks_passed": True}, indent=2))


if __name__ == "__main__":
    main()
