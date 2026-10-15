"""Journal-sized four-panel figure, complete alternative summaries and concise report."""
import os
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import gaussian_kde

B = Path(__file__).resolve().parents[1]
T, F = B / "tables", B / "figures"


def main():
    summary = json.loads((T / "summary.json").read_text())
    checks = json.loads((T / "independent_numerical_checks.json").read_text())
    assert checks["all_checks_passed"]
    results = pd.read_csv(T / "all_923_cluster_directional_results.csv")
    common = results.loc[results.common_measurable].copy()
    animals = pd.read_csv(T / "RRRM2_equal_cluster_animal_scores.csv")
    bootstrap = pd.read_csv(T / "animal_bootstrap_global_draws.csv.gz")
    sensitivity = pd.read_csv(T / "measurement_sensitivities.csv")
    strata = pd.read_csv(T / "descriptive_stratum_effects.csv")
    comparison = pd.read_csv(T / "binary_multiway_measurement_transport_comparison.csv")
    # Complete descriptive component means; no component-wise hypothesis tests.
    alt = pd.read_csv(T / "all_3376_alternative_group_composition_summaries.csv.gz")
    tcaf = results.loc[results.gene_symbol.eq("Tcaf1")].copy()
    tcaf_alt = alt.loc[alt.cluster_id.isin(tcaf.loc[tcaf.common_measurable, "cluster_id"])]
    tcaf_alt.to_csv(T / "Tcaf1_all_frozen_alternative_mean_contrasts.csv", index=False)
    # Arial/Helvetica at final physical size; fail explicitly rather than silently changing typography.
    chosen = None
    for family in ["Arial", "Helvetica"]:
        try:
            font_manager.findfont(family, fallback_to_default=False)
            chosen = family
            break
        except ValueError:
            continue
    if chosen is None:
        raise RuntimeError("Publication figure requires an installed Arial or Helvetica font.")
    plt.rcParams.update({"font.family": chosen, "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
                         "axes.linewidth": 1, "xtick.major.width": 1, "ytick.major.width": 1,
                         "lines.linewidth": 1, "patch.linewidth": 1, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    blue, orange, slate = "#247887", "#C75B36", "#788B98"
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 163 / 25.4))
    fig.subplots_adjust(left=.115, right=.985, top=.94, bottom=.12, wspace=.46, hspace=.54)
    ax = axes[0, 0]
    coverage_count = 923 - summary["RRRM2_coverage_fail"]
    values = [923, coverage_count, len(common)]
    bars = ax.barh([2, 1, 0], values, height=.55, color=[slate, slate, blue], linewidth=1)
    ax.set_yticks([2, 1, 0], ["Frozen RR23\nclusters", "RRRM2\ncoverage gate", "Measured in\nboth cohorts"])
    ax.set_xlim(0, 1080)
    ax.set_xlabel("Number of clusters")
    ax.set_title("a) Measurement transport", loc="left", weight="bold", pad=10)
    for bar, value in zip(bars, values):
        ax.text(value + 25, bar.get_y() + bar.get_height() / 2, f"{value:,}", va="center", fontsize=8)
    ax = axes[0, 1]
    ax.scatter(100 * common.RR23_training_native_L2_norm, 100 * common.RRRM2_directional_effect,
               s=14, color=blue, alpha=.65, edgecolors="none", linewidths=1)
    if tcaf.common_measurable.any():
        t = tcaf.loc[tcaf.common_measurable].iloc[0]
        ax.scatter([100 * t.RR23_training_native_L2_norm], [100 * t.RRRM2_directional_effect],
                   s=48, marker="D", facecolors="none", edgecolors=orange, linewidths=1.1, label="Tcaf1 illustration")
        ax.legend(loc="lower right", frameon=False, handletextpad=.3)
    ax.axhline(0, color=slate, linewidth=1)
    ax.set_xlabel("RR23 native contrast norm (pp)")
    ax.set_ylabel("RRRM2 directional contrast (pp)")
    ax.set_title("b) All measurable clusters", loc="left", weight="bold", pad=10)
    ax.text(.04, .96, f"BH: {summary['local_BH_full923_05']}; BY: {summary['local_BY_full923_05']}", transform=ax.transAxes,
            va="top", fontsize=8, bbox=dict(facecolor="white", edgecolor="none", alpha=.95))
    ax = axes[1, 0]
    conditional = 100 * bootstrap.conditional_fixed_direction_global_effect.to_numpy()
    joint = 100 * bootstrap.joint_retrained_direction_global_effect.to_numpy()
    extent = [min(np.quantile(conditional, .002), np.quantile(joint, .002)),
              max(np.quantile(conditional, .998), np.quantile(joint, .998))]
    grid = np.linspace(extent[0], extent[1], 400)
    ax.plot(grid, gaussian_kde(conditional)(grid), color=blue, linewidth=1.2, label="Fixed directions")
    ax.plot(grid, gaussian_kde(joint)(grid), color=slate, linestyle="--", linewidth=1.2, label="Both cohorts")
    ax.axvline(0, color=slate, linestyle=":", linewidth=1)
    ax.axvline(summary["global_equal_cluster_effect_pp"], color=orange, linewidth=1.2)
    ax.set_xlabel("Equal-cluster directional contrast (pp)")
    ax.set_ylabel("Animal-bootstrap density")
    ax.set_title("c) Global score uncertainty", loc="left", weight="bold", pad=10)
    ax.legend(loc="upper right", frameon=False)
    ax.text(.03, .04, f"Exact one-sided p = {summary['global_exact_one_sided_p']:.3f}", transform=ax.transAxes,
            fontsize=8, bbox=dict(facecolor="white", edgecolor="none", alpha=.95))
    ax = axes[1, 1]
    order = ["ISS_T_12", "LAR_12", "LAR_29"]
    rng = np.random.default_rng(3115)
    for pos, st in enumerate(order):
        subset = animals.loc[animals.stratum.eq(st)]
        center = subset.equal_cluster_directional_score.mean()
        for exposure, color, shift in [("GC", "#D55E00", -.16), ("FLT", "#0072B2", .16)]:
            group = subset.loc[subset.exposure.eq(exposure)]
            y = 100 * (group.equal_cluster_directional_score - center)
            ax.scatter(pos + shift + rng.uniform(-.025, .025, len(group)), y, s=26,
                       color=color, edgecolors="white", linewidths=1, label=exposure if pos == 0 else None, zorder=3)
            ax.plot([pos + shift - .06, pos + shift + .06], [y.mean(), y.mean()], color="black", linewidth=1.2, zorder=4)
    ax.axhline(0, color=slate, linewidth=1)
    ax.set_xticks(range(3), ["ISS-T\n12 wk", "LAR\n12 wk", "LAR\n29 wk"])
    ax.set_ylabel("Centered animal score (pp)")
    ax.set_title("d) Held-out biological animals", loc="left", weight="bold", pad=10)
    ax.legend(loc="upper right", frameon=False, ncol=2, columnspacing=.5, handletextpad=.1)
    # Text and line audit in the physical plotting canvas.
    text_sizes = [text.get_fontsize() for text in fig.findobj(matplotlib.text.Text) if text.get_text()]
    assert min(text_sizes) >= 8
    F.mkdir(exist_ok=True)
    fig.savefig(F / "Figure_Multiway_Transport.pdf")
    fig.savefig(F / "Figure_Multiway_Transport.svg")
    fig.savefig(F / "Figure_Multiway_Transport.png", dpi=600)
    plt.close(fig)
    figure_metadata = dict(width_mm=183, height_mm=163, dpi_png=600, font_family=chosen,
                           minimum_text_font_points=min(text_sizes), configured_minimum_line_points=1,
                           lowercase_bold_panel_labels=["a)", "b)", "c)", "d)"],
                           figure_title_inside_canvas=False, footer_inside_canvas=False,
                           PDF_SVG_vector=True, PNG_expected_pixels=[round(183 / 25.4 * 600), round(163 / 25.4 * 600)])
    (T / "figure_format_metadata.json").write_text(json.dumps(figure_metadata, indent=2) + "\n")
    conditional_ci = np.array(summary["conditional_fixed_direction_bootstrap_CI95"]) * 100
    joint_ci = np.array(summary["joint_retrained_direction_bootstrap_CI95"]) * 100
    caption = f"""Independent-launch transport of frozen RR23 multiway junction compositions into RRRM2 cerebellum. (a) The complete family comprises 923 clusters; {coverage_count} pass RRRM2 coverage and {len(common)} (30.44%) meet both cohort measurement gates. Each animal requires ≥20 reads across the fixed alternatives and ≥90% dominance among all retained annotated canonical alternatives sharing the boundary. (b) Every common-cluster RR23 FLT–GC native contrast norm and its held-out directional effect. Unit Euclidean directions were learned solely from nine RR23 FLT and nine GC animals, without significance selection. No cluster passes BH or BY across all 923; failed measurements retain p=1. Tcaf1 is a historical illustration, not a selected validation family. (c) Kernel density displays of 2,000 whole-animal bootstrap draws: solid, held-out animals with directions fixed; dashed, both cohorts with RR23 directions re-estimated. Orange marks the observed equal-cluster contrast ({summary['global_equal_cluster_effect_pp']:.3f} pp); dotted line marks zero. Conditional percentile 95% interval [{conditional_ci[0]:.3f}, {conditional_ci[1]:.3f}] pp; two-cohort interval [{joint_ci[0]:.3f}, {joint_ci[1]:.3f}] pp. Exact positive-direction p={summary['global_exact_one_sided_p']:.5f} across all 11,200 stratified assignments (absolute two-sided p={summary['global_exact_absolute_two_sided_p']:.6f}). (d) All 19 held-out animals, centered by a common stratum mean only for display; points are animals and ticks are group means. Strata contain 5/3, 2/3 and 3/3 FLT/GC animals with weights 8/19, 5/19 and 6/19. Blue denotes FLT and orange GC. Positive projections mean alignment with RR23; percentage points scale a linear composition projection, not percentages of genes altered. Clusters and bootstrap draws are not independent biological replicates. RR23 flight and collection day are aliased; this comparison identifies neither a causal flight effect nor functional injury.\n"""
    assert len(caption.split()) <= 350
    (F / "Figure_Multiway_Transport_caption.txt").write_text(caption)
    tcaf_text = "Tcaf1 was not measurable under the common gate and no held-out composition was forced."
    if tcaf.common_measurable.any():
        t = tcaf.loc[tcaf.common_measurable].iloc[0]
        tcaf_text = (f"The historically selected Tcaf1 cluster `{t.cluster_id}` is measurable (three fixed alternatives; RRRM2 minimum selected reads {int(t.RRRM2_minimum_selected_reads)}, minimum dominance {100*t.RRRM2_minimum_dominance:.2f}%). Its held-out directional effect is {100*t.RRRM2_directional_effect:.3f} pp (conditional 95% CI {100*t.conditional_bootstrap_CI95_low:.3f} to {100*t.conditional_bootstrap_CI95_high:.3f}), one-sided p={t.exact_one_sided_directional_p:.6f}, two-sided p={t.exact_absolute_two_sided_p:.6f}; BH/BY/Holm q values are all 1. Its native held-out contrast has cosine alignment {t.RRRM2_direction_cosine:.3f} with RR23, a descriptive opposite alignment with wide uncertainty, not evidence of a reversed causal mechanism. All three alternatives and all 45 RR23/RRRM2 animal compositions are preserved; no alternative was chosen to make the comparison favorable.")
    sensitivity_lines = "\n".join(f"| {row.panel} | {row.common_clusters} | {100*row.global_effect:.3f} | {row.exact_one_sided_p:.6f} |" for row in sensitivity.itertuples())
    stratum_lines = "\n".join(f"| {row.stratum} | {row.n_FLT}/{row.n_GC} | {100*row.effect:.3f} | [{100*row.CI95_low:.3f}, {100*row.CI95_high:.3f}] |" for row in strata.itertuples())
    report = f"""# Frozen multiway composition transport: scientific report

Author: Isarar Siddique  
Affiliation:  
Email:  
Analysis date: 2026-10-15

## Main result

The measured multiway RR23 FLT–GC composition profile does not show aggregate positive-direction transfer to independently launched RRRM2 cerebellum. The equal-cluster held-out effect is **{summary['global_equal_cluster_effect_pp']:.3f} projected percentage points**, with exact one-sided **p={summary['global_exact_one_sided_p']:.5f}** and absolute two-sided p={summary['global_exact_absolute_two_sided_p']:.6f} across all 11,200 stratified assignments. The 2,000-draw fixed-direction animal-bootstrap 95% interval is [{conditional_ci[0]:.3f}, {conditional_ci[1]:.3f}] pp. Independently resampling both cohorts and re-estimating RR23 directions gives a wider interval [{joint_ci[0]:.3f}, {joint_ci[1]:.3f}] pp. Neither interval excludes zero; the negative point estimate does not establish a reversed response.

Only **{len(common)}/923 clusters ({100*summary['common_measurement_fraction']:.2f}%)** are measurable in both cohorts. Of these, {100*summary['positive_directional_cluster_fraction']:.2f}% have positive learned-direction effects (conditional animal-bootstrap 95% interval {100*summary['positive_directional_cluster_fraction_conditional_bootstrap_CI95'][0]:.2f}% to {100*summary['positive_directional_cluster_fraction_conditional_bootstrap_CI95'][1]:.2f}%). Mean native cosine alignment is {summary['mean_native_cosine_alignment']:.4f}, reported descriptively without cluster-independence p values.

## Design and estimand

The retrospective plan was recorded before inspecting these RRRM2 outcomes; previous binary transport and RR23 multiway results were already known. All 923 all26-defined clusters and all 3,376 alternative occurrences were retained. The analysis did not use RR23 significance, favorable effect size, positive held-out direction, gene annotation or Tcaf1 status to select the testing family.

Native composition divides each frozen alternative count by the sum over that cluster's fixed alternatives. Coverage and dominance were evaluated in every animal: ≥20 selected-alternative reads and ≥90% dominance among all retained annotated canonical alternatives sharing the boundary. RR23 gates use all 26 animals; RRRM2 uses all 19 included animals. Coverage fails in 604 clusters and dominance in 312, with 274 failing both. The disjoint failures are 330 coverage only, 38 dominance only and 274 both. Inadequately measured clusters receive p=1 and no interpretable held-out effect; an absent single alternative can have a valid zero fraction when its complete cluster passes.

For each cluster, the RR23 9-FLT minus 9-GC native mean vector is normalized by its Euclidean norm. No training direction is estimated from RRRM2 labels. The projected held-out score subtracts the RR23 GC mean composition and takes a dot product with that frozen unit direction. The global score averages all common clusters equally. None of the 923 training directions met the predeclared near-zero threshold. Unit normalization avoids weighting a cluster by its RR23 effect magnitude but can leave uncertain directions for small training effects; the separate two-cohort bootstrap exposes this uncertainty rather than silently selecting stable or larger effects.

The prior RRRM2 inclusion rule was preserved: the entire ISS_T×29-week stratum, including its corresponding controls, remains excluded because its four FLT records have conflicting location metadata. Included strata are ISS_T_12 (5 FLT/3 GC), LAR_12 (2/3), and LAR_29 (3/3). Within-stratum contrasts receive fixed animal-count weights 8/19, 5/19 and 6/19. Whole-animal allocation preserves every within-animal alternative/cluster dependency. Its interpretation requires conditional exchangeability within these strata; retrospective labels were not physically randomized.

## Complete local correction and exact-test resolution

All **923** directional local p values were corrected together, with failed measurements retained at p=1. There are **zero BH and zero BY discoveries** at 0.05; the smallest raw directional p is {summary['smallest_local_one_sided_p']:.6f}. BH is an exploratory dependence-assumption sensitivity, while BY is the conservative step-up map under valid p values. Zero local calls do not establish that every true effect is zero.

Zero Holm calls are structurally guaranteed at the declared family size and label space: the smallest attainable raw one-sided p is 1/11,200=0.000089286, whereas the first Holm step requires p≤0.05/923=0.000054171. Its minimum attainable first-step adjusted p is **923/11,200=0.0824107**. This first-step bound applies regardless of the observed effect and must not be described as biological evidence of absence. BY can still reject through its step-up rule if sufficiently many small p values occur; none did here.

The earlier full 1,512-common-feature binary directional family has the same issue: its minimum attainable first-step Holm-adjusted p is **1,512/11,200=0.135**, exceeding 0.05. This statement concerns the first sorted p value and the consequent inability to initiate any Holm rejection, not a universal minimum of every unadjusted or BH/BY p value. It does not affect the separate one-score global tests. No correction family was reduced to rescue significance.

## Measurement comparison and candidate illustration

The binary panel measured 1,512 of 5,249 RR23 measurable features (28.81%); this multiway panel measures 281 of 923 (30.44%). These are distinct observables with different selection and geometry, so the proportions are descriptive and do not show that one assay has statistically superior transport or that their effect magnitudes are directly comparable.

{tcaf_text}

## Sensitivities and strata

Two planned stricter gates retain the primary direction definition and use the complete stricter common panel. Neither rescues aggregate positive-direction transport; their two one-sided sensitivity p values have joint Holm values of 1.

| Measurement panel | Common clusters | Directional contrast (pp) | Exact one-sided p |
|---|---:|---:|---:|
{sensitivity_lines}

Stratum effects are descriptive, based on very small cells, and are not independent validation datasets:

| Stratum | FLT/GC | Directional contrast (pp) | Conditional bootstrap 95% interval (pp) |
|---|---:|---:|---:|
{stratum_lines}

## Interpretation and reproducibility

The additional real-data contribution is a frozen multialternative composition-direction test and an explicit audit of transportability, measurement failure and finite exact-test resolution. It provides an independent-launch boundary on the generality of RR23 observations, not a new conserved flight-injury mechanism. Mission, sex, strain, age, recovery/collection route, tissue sampling, preservation and sequencing coverage differ between RR23 and RRRM2 and cannot be separately identified as causes of failed transfer. RR23 flight and collection date are aliased. No protein, functional injury, cell-intrinsic or isolated causal flight claim follows. RR10 shares SpaceX-21 with RR23 and is not the independent-launch validation cohort.

All {checks['independent_check_count']} independent numerical checks passed. A separate implementation read all 805,562 RR23 and 668,237 RRRM2 count rows in chunks and independently recomputed every frozen alternative count, composition, full-boundary total and dominance gate. It checked every RR23-only unit vector and held-out projection, directly enumerated all 11,200 global allocations, reproduced 17 row-order-selected local exact p values and all 923 BH/BY/Holm values, and verified both whole-animal bootstrap series and intervals. Maximum independent unit-direction arithmetic difference was {checks['maximum_independent_unit_direction_absolute_difference']:.2e}, below the declared 1e−12 tolerance. All frozen-source hashes and the dated plan remain unchanged.

The extension is self-contained: 15 frozen input files include actual count matrices, full animal metadata, all cluster memberships/features, annotation lookup and NASA publisher provenance. The inherited RRRM2 provenance documents 26 raw junction files matching publisher MD5 and complete reconstruction of the 19 retained columns. Reproduction requires no path to the source repository. Run `scripts/analyse.py`, `scripts/verify_independently.py`, then `scripts/build_outputs.py` from any directory. The figure uses 183 mm width, Arial/Helvetica at ≥8 pt, ≥1 pt lines, lowercase bold panel labels, editable vector PDF/SVG and 600 dpi PNG, with title/caption outside the canvas.

Complete outputs include `tables/all_923_cluster_directional_results.csv`, `all_923_cluster_measurement_gates.csv`, `all_3376_alternative_group_composition_summaries.csv.gz`, all animal/alternative counts and fractions, frozen training vectors, all held-out animal scores, the full 11,200-assignment null, saved 2,000-draw conditional and two-cohort bootstraps, Tcaf1's complete observable case, and independent check records. `figures/Figure_Multiway_Transport_caption.txt` is the separate manuscript caption.
"""
    (B / "SCIENTIFIC_REPORT.md").write_text(report)
    print(json.dumps(dict(figure=str(F / "Figure_Multiway_Transport.pdf"), report=str(B / "SCIENTIFIC_REPORT.md"), figure_format=figure_metadata), indent=2))


if __name__ == "__main__":
    main()
