#!/usr/bin/env python3
"""Build a review draft and original scientific figures from revised CSVs.

Numerical assertions are looked up from results, not copied from legacy prose.
The output is an exploratory author-review draft, not a submission-ready paper.
"""
import html
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate,
                              Spacer, Table, TableStyle)
import reportlab

HERE = Path(__file__).resolve().parents[1]
R = HERE / "results"
CLAIMS = []


def table_lookup(frame, filename, filters, field):
    select = np.ones(len(frame), dtype=bool)
    for name, value in filters.items():
        select &= frame[name].to_numpy() == value
    row = frame.loc[select]
    if len(row) != 1:
        raise ValueError(f"Ambiguous manuscript lookup {filters}: {len(row)} rows")
    number = float(row.iloc[0][field])
    CLAIMS.append({"source_file": filename, "filters": json.dumps(filters, sort_keys=True),
                   "field": field, "value": number})
    return number


def main():
    upstream = json.loads((HERE / "inputs/upstream_count_provenance.json").read_text())
    if len(upstream) != 2 or any(not r["equal_gene_and_sample_ids"] or r["max_count_difference"] != 0 for r in upstream):
        raise ValueError("Core NASA count origin has not been verified")
    panels = pd.read_csv(R / "panel_results.csv")
    gradients = pd.read_csv(R / "gradient_results.csv")
    anchors = pd.read_csv(R / "ground_anchor_sensitivity.csv")
    splits = pd.read_csv(R / "OSD525_disjoint_anchor_partitions.csv")
    accounting = pd.read_csv(R / "dataset_accounting.csv")
    estimator = pd.read_csv(R / "estimator_sensitivity.csv")
    run = json.loads((R / "run_summary.json").read_text())
    authors = json.loads((HERE / "authors.json").read_text())
    figures = HERE / "figures"
    figures.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 13,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "savefig.dpi": 600})
    colorset = {"rPRG": "#1b7898", "HSP": "#d88824", "SRG": "#a44e85"}
    captions = []

    def save(fig, name, caption):
        fig.savefig(figures / f"{name}.png", dpi=600, bbox_inches="tight", facecolor="white")
        fig.savefig(figures / f"{name}.svg", bbox_inches="tight", facecolor="white")
        plt.close(fig)
        captions.append((name, caption))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), gridspec_kw={"width_ratios": [1, 1.5]})
    core = [53, 29]
    axes[0].bar(["Region libraries", "Source animals"], core, color=["#657787", "#1b7898"], width=.65)
    for i, v in enumerate(core):
        axes[0].text(i, v + 1, str(v), ha="center", weight="bold")
    axes[0].set_ylim(0, 65)
    axes[0].set_ylabel("Count in OSD-561 + OSD-562")
    axes[0].set_title("a  Shared animals across regions", loc="left", fontsize=14)
    labels = ["Pooled ground", "Pooled flight", "Hippocampus flight"]
    contrasts = ["RRRM2_pooled_GC_ISST_vs_LAR", "RRRM2_pooled_FLT_ISST_vs_LAR", "OSD562_hippocampus_FLT_ISST_vs_LAR"]
    for j, panel in enumerate(["rPRG", "HSP", "SRG"]):
        sub = panels[(panels.sensitivity == "design_cohort") & (panels.panel == panel)].set_index("contrast").loc[contrasts]
        ypos = np.arange(3) + (j - 1) * .18
        axes[1].errorbar(sub.mean_log2_effect, ypos,
                         xerr=np.vstack([sub.mean_log2_effect - sub.bootstrap_lower95,
                                         sub.bootstrap_upper95 - sub.mean_log2_effect]),
                         fmt="o", capsize=3, color=colorset[panel], label=panel, markersize=5)
    axes[1].axvline(0, color="#888888", lw=.8)
    axes[1].set_yticks(range(3), labels)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Mean panel log2 effect (ISS-T − LAR)")
    axes[1].legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(.5, -.18))
    axes[1].set_title("b  Response-panel effects", loc="left", fontsize=14)
    fig.tight_layout(w_pad=2)
    save(fig, "figure1_sample_accounting_and_panel_effects",
         "Figure 1. Sample accounting and response-panel effects. (a) The cerebellum and hippocampus matrices contain 53 region libraries from 29 Source Name animals. (b) Points are mean gene-level log2 effects; bars are 95% percentile intervals from 1,000 animal bootstraps within label and age. All regions of a sampled animal remain together. Panel labels are inherited exploratory gene sets, not measured transcriptional time windows. Intervals are conditional on normalization and cohort exchangeability.")

    primary = panels[panels.sensitivity == "design_cohort"]
    matrix = primary.pivot(index="contrast", columns="panel", values="mean_log2_effect")
    qmatrix = primary.pivot(index="contrast", columns="panel", values="q_BH_all_tests_in_family")
    chosen = ["rPRG", "HSP", "SRG", "Glia_Complement", "BBB_Endothelium", "Circadian"]
    labels = [name.replace("RRRM2_pooled", "RRRM-2 pooled").replace("OSD561_cerebellum", "OSD-561 Cb")
              .replace("OSD562_hippocampus", "OSD-562 Hip").replace("_ISST_vs_LAR", " ISS-T / LAR")
              .replace("_FLT_vs_GC", " FLT / GC").replace("_vs_", " / ").replace("_", " ") for name in matrix.index]
    fig, ax = plt.subplots(figsize=(9.4, 6.5))
    vmax = np.abs(matrix[chosen].to_numpy()).max()
    im = ax.imshow(matrix[chosen], cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(range(len(labels)), labels, fontsize=13)
    ax.set_xticks(range(len(chosen)), ["rPRG", "HSP", "SRG", "Glia / complement", "BBB / endothelium", "Circadian"], rotation=25, ha="right")
    for i, name in enumerate(matrix.index):
        for j, p in enumerate(chosen):
            value = matrix.loc[name, p]
            star = "*" if qmatrix.loc[name, p] < .05 else ""
            ax.text(j, i, f"{value:.2f}{star}", ha="center", va="center", fontsize=12,
                    color="white" if abs(value) > .6 * vmax else "black")
    fig.colorbar(im, ax=ax, label="Mean panel log2 effect", shrink=.8, pad=.03)
    fig.tight_layout()
    save(fig, "figure2_cross_cohort_panel_effects",
         "Figure 2. Unadjusted panel effects across 17 primary design contrasts. Numbers show mean log2 effects; an asterisk denotes BH q < 0.05 across all 232 panel tests, including metadata sensitivities. Positive effects follow the numerator in each row. Cb, cerebellum; Hip, hippocampus; HGC, habitat ground control; VGC, vivarium control. OSD-563, OSD-564 and OSD-612 are three regions of one RR-10 animal cohort, not independent missions. Glial and endothelial transcript abundance does not directly measure inflammation or barrier permeability.")

    sub = anchors[anchors["mode"].isin(["projection", "score"])].copy()
    names = sub.comparison.drop_duplicates().tolist()
    fig, ax = plt.subplots(figsize=(10, 8.3))
    for j, (mode, color) in enumerate([("projection", "#1b7898"), ("score", "#d88824")]):
        part = sub[sub["mode"] == mode].set_index("comparison").loc[names]
        ypos = np.arange(len(names)) + (j - .5) * .22
        ax.scatter(part.outside_sensor_effect_energy_removed_pct, ypos, label=mode.capitalize(), color=color, s=24)
        for y0, (_, row) in zip(ypos, part.iterrows()):
            ax.plot([row.anchor_bootstrap_energy_lower95, row.anchor_bootstrap_energy_upper95], [y0, y0], color=color, lw=1.2)
    ax.axvline(0, color="#555555", lw=.8)
    short = [n.replace("anchor_OSD", "OSD-").replace("_to_OSD", " → OSD-")
             .replace("_FLT_recovery", " recovery contrast").replace("_FLT_vs_GC", " FLT / GC")
             .replace("OSD525_shared_HGC_anchor", "OSD-525 shared HGC (9 animals)") for n in names]
    ax.set_yticks(range(len(names)), short, fontsize=13)
    ax.invert_yaxis()
    ax.set_xlabel("Reduction in squared log2-effect energy outside Sensor31 (%)\nNegative values indicate increased energy; neither sign establishes biological accuracy")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    save(fig, "figure3_ground_anchor_transfer",
         "Figure 3. Ground-anchor sensitivity and transfer. Dots show the observed change in squared contrast coefficients outside the Sensor31 genes after genome-expression projection or ground-calibrated score subtraction. Lines show 95% percentile ranges from 200 ground-animal bootstraps within age and cohort, conditional on normalized counts and a fixed target. These ranges describe anchor uncertainty, not complete sampling intervals or calibrated gene-level inference. The within-OSD-525 evaluation shares nine HGC animals with its anchor and is displayed as such. Negative reductions are retained. Removing expression variation is not proof of removing confounding.")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for mode, color in [("projection", "#1b7898"), ("score", "#d88824")]:
        s = splits[splits["mode"] == mode]
        axes[0].scatter(s.sensor_energy_removed_pct, s.outside_sensor_energy_removed_pct, s=13, alpha=.6, color=color, label=mode.capitalize())
        axes[1].hist(s.lambda_conditional_diagnostic, bins=np.linspace(1, 6.5, 17), color=color, alpha=.55, label=mode.capitalize())
    axes[0].axhline(0, color="#777777", lw=.8)
    axes[0].set_xlabel("Sensor31 effect-energy reduction (%)")
    axes[0].set_ylabel("Outside-sensor effect-energy reduction (%)")
    axes[1].set_xlabel("Conditional λ diagnostic outside Sensor31")
    axes[1].set_ylabel("Overlapping partition count")
    axes[0].legend(frameon=False)
    axes[0].set_title("a  Sensor removal and external effects", loc="left", fontsize=14)
    axes[1].set_title("b  Conditional dispersion", loc="left", fontsize=14)
    fig.tight_layout()
    save(fig, "figure4_disjoint_anchor_partitions",
         "Figure 4. RR-23 disjoint-anchor partitions. Each of 126 partitions trains on four of nine HGC animals plus the same eight VGC animals, then evaluates the nine FLT animals against the remaining five HGC animals. Training and evaluation animals are disjoint within a partition. Partitions reuse the same cohort and do not constitute 126 independent validations. (a) Removing sensor effect energy does not ensure a reduction outside the sensor. (b) λ is derived from conditional OLS t probabilities and is a dispersion diagnostic, not a measure of biological truth or empirical false-positive rate.")

    if "--figures-only" in sys.argv:
        (R / "figure_manifest.json").write_text(json.dumps([
            {"filename": name, "caption": caption, "png_dpi": 600, "vector_format": "SVG"}
            for name, caption in captions], indent=2) + "\n")
        return

    def panel(c, p, field="mean_log2_effect", sensitivity="design_cohort"):
        return table_lookup(panels, "panel_results.csv", {"contrast": c, "panel": p, "sensitivity": sensitivity}, field)

    def gradient(c, field="rPRG_minus_SRG_mean_log2_effect", sensitivity="design_cohort"):
        return table_lookup(gradients, "gradient_results.csv", {"contrast": c, "sensitivity": sensitivity}, field)

    def anchor(c, mode, field):
        return table_lookup(anchors, "ground_anchor_sensitivity.csv", {"comparison": c, "mode": mode}, field)

    pgc, pfl, hip = "RRRM2_pooled_GC_ISST_vs_LAR", "RRRM2_pooled_FLT_ISST_vs_LAR", "OSD562_hippocampus_FLT_ISST_vs_LAR"
    table_rows = "\n".join(f"| {r.accession} | {r.mission} | {r.tissue} | {r.n_deposited_libraries_in_count_matrix} | {r.n_animals_in_region} |"
                           for r in accounting.itertuples())
    names = ", ".join(a["name"] for a in authors["authors"])
    contribution = "\n\n".join(f"**{a['name']}:** " + "; ".join(a["credit"]) + "." for a in authors["authors"])
    glia_min_q = float(panels.loc[(panels.panel == "Glia_Complement") & (panels.sensitivity == "design_cohort"), "q_BH_all_tests_in_family"].min())
    split_stats = splits.groupby("mode")[["sensor_energy_removed_pct", "outside_sensor_energy_removed_pct"]].agg(["mean", "std"])
    draft = f"""# Cohort-associated transcriptional variation and limits of ground-anchored adjustment in spaceflight mouse brain

{names}

{authors['affiliation']}

Correspondence: {authors['corresponding_author']} ({authors['corresponding_email']}). ORCID: 0009-0005-2305-9579.

**Author-review revision — 14 October 2026. Working journal target: npj Microgravity. Metadata adjudication, final validation and author declarations remain unresolved.**

## Abstract

Spaceflight brain transcriptomes can combine exposure-related biology with collection, handling and animal-cohort differences. We reanalysed seven NASA Open Science Data Repository bulk RNA-sequencing deposits with explicit animal identities, consistent technical-replicate handling and animal-restricted inference. The primary cerebellum and hippocampus deposits comprised 53 region libraries from 29 animals; seven deposited location fields conflicted with design identifiers. Ground-cohort differences included rapid-response and heat-shock expression, but secondary-response expression also changed in flight hippocampus. The pooled flight rapid-minus-secondary effect was {gradient(pfl):.3f} log2 units (two-sided animal-permutation p = {gradient(pfl, 'p_animal_permutation'):.3f}), which does not establish a selective rapid-response program. Two ground-trained expression transformations attenuated sensor effects in some settings but showed variable transfer and could increase effect energy elsewhere. In 126 overlapping RR-23 partitions with disjoint training and evaluation animals, score subtraction reduced sensor energy by {split_stats.loc['score', ('sensor_energy_removed_pct', 'mean')]:.1f}% on average while increasing outside-sensor energy by {-split_stats.loc['score', ('outside_sensor_energy_removed_pct', 'mean')]:.1f}%. These results support careful cohort auditing and transparent sensitivity analysis; they do not identify a precise peri-mortem clock or validate universal removal of technical variation. Independent collection-controlled replication is needed before assigning corrected expression changes to spaceflight mechanisms.

Keywords: spaceflight; mouse brain; bulk RNA sequencing; cohort confounding; sample provenance; ground controls.

## Introduction

Brain expression studies from spaceflight offer opportunities to examine responses to an unusual environmental exposure, but the measured contrast depends on how animals were housed, collected, preserved and assigned to cohorts. A flight-versus-ground contrast may contain several effects that cannot be distinguished with an endpoint expression matrix alone. More generally, batch effects can dominate high-throughput measurements and complicate biological interpretation [1]. An adjustment that removes variation may also remove a real response if it aligns with the learned direction.

Immediate and delayed transcriptional responses motivate inspection of response-gene panels. Neuronal activity patterns can induce different expression programs [2], and immediate-early and delayed primary response genes differ in their properties [3]. These observations do not assign a unique elapsed time, euthanasia mechanism or activity state to every gene in a heterogeneous bulk brain sample. Gene-level counts cannot resolve transcript isoforms or directly demonstrate promoter pause release. Here, previously assembled panels are treated as exploratory summaries rather than validated clocks.

The NASA Open Science Data Repository (OSDR) provides processed counts and sample metadata for multiple mouse brain regions. We examined whether explicit animal accounting and ground-cohort contrasts alter interpretation of a previous integrative analysis. We then evaluated two ground-trained linear transformations as sensitivity analyses. The objective was to characterize cohort-associated variation and the limits of transfer, rather than assume that a ground contrast uniquely measures artifact or that its subtraction recovers a true flight response.

## Results

### Correct dataset identities and animal accounting

The revision includes seven bulk RNA-sequencing deposits, totaling {run['n_libraries']} count-matrix libraries and {run['n_unique_source_animals']} distinct deposited Source Name identifiers. These are three rodent cohorts with repeated regional sampling, not seven independent missions. Cerebellum and hippocampus from OSD-561 and OSD-562 share animals: 53 libraries correspond to 29 animals. OSD-613 comprises a separate 18-animal cerebral-hemisphere subset of RRRM-2. OSD-563, OSD-564 and OSD-612 share ten RR-10 animals, with one flight animal absent from the hemisphere matrix. OSD-525 has 29 libraries from 26 RR-23 animals: nine FLT, nine HGC and eight VGC; three animals have two technical libraries.

| Accession | Cohort | Tissue | Count-matrix libraries | Animals within region |
|---|---|---|---:|---:|
{table_rows}

Table 1. Included processed bulk datasets. Regional animal counts are not additive across a cohort. Library totals refer to included count matrices, not every sample or assay listed in the public deposit.

OSD-525 is cerebellum, OSD-564 is hippocampus, and OSD-612 and OSD-613 are left cerebral hemisphere. The seven included deposits exclude unrelated kidney and human histone-assay datasets previously assigned brain identities. RRRM-2 and RR-10 have different strains and study designs. RR-10 and RR-23 are separate animal cohorts but share a SpaceX CRS-21 launch; launch independence is therefore not asserted.

Seven RRRM-2 region records from six animals have deposited euthanasia-location fields that disagree with Source Name cohort identifiers. All seven records and both fields are retained in the conflict table. Removing every region of a conflicting animal and, separately, using the location field within a region provide sensitivity analyses; neither resolves the deposited disagreement. A ground ISS-T designation names a matched ground cohort and does not mean that the ground animal was on the ISS.

### Ground and flight response-panel effects differ by region

In pooled RRRM-2 ground animals, the rapid-response panel mean effect was {panel(pgc, 'rPRG'):.3f} log2 units (animal-permutation p = {panel(pgc, 'rPRG', 'p_animal_permutation'):.4f}; BH q = {panel(pgc, 'rPRG', 'q_BH_all_tests_in_family'):.3f}). The secondary-response mean effect was {panel(pgc, 'SRG'):.3f}. Their difference was {gradient(pgc):.3f} log2 units (p = {gradient(pgc, 'p_animal_permutation'):.4f}; gradient-family q = {gradient(pgc, 'q_BH_all_tests_in_family'):.3f}). The ground rapid-panel effect after exclusion of all regions from conflicting animals remained {panel(pgc, 'rPRG', sensitivity='exclude_all_regions_of_conflicting_animals'):.3f}.

In the pooled flight recovery-exposure contrast, rapid and secondary panels had effects of {panel(pfl, 'rPRG'):.3f} and {panel(pfl, 'SRG'):.3f} log2 units. The rapid-minus-secondary difference was {gradient(pfl):.3f}, with p = {gradient(pfl, 'p_animal_permutation'):.3f} and gradient-family q = {gradient(pfl, 'q_BH_all_tests_in_family'):.3f}. Flight hippocampus showed a secondary-panel effect of {panel(hip, 'SRG'):.3f} (p = {panel(hip, 'SRG', 'p_animal_permutation'):.4f}; panel-family q = {panel(hip, 'SRG', 'q_BH_all_tests_in_family'):.3f}). Thus, a universally flat secondary response and a rapid-only flight program are not supported.

RR-23 HGC-versus-VGC rapid and secondary effects were {panel('OSD525_HGC_vs_VGC', 'rPRG'):.3f} and {panel('OSD525_HGC_vs_VGC', 'SRG'):.3f}; the rapid-minus-secondary comparison had p = {gradient('OSD525_HGC_vs_VGC', 'p_animal_permutation'):.3f}. RR-23 FLT-versus-HGC rapid and secondary effects were {panel('OSD525_FLT_vs_HGC', 'rPRG'):.3f} and {panel('OSD525_FLT_vs_HGC', 'SRG'):.3f}. These direction changes underscore dependence on the chosen comparator. The comparisons represent endpoint cohort differences, not a measured temporal progression.

### Ground-anchor attenuation does not establish accurate biological correction

For the matched hippocampus ground anchor applied to the flight recovery-exposure contrast, score subtraction removed {anchor('anchor_OSD562_to_OSD562_FLT_recovery', 'score', 'sensor_effect_energy_removed_pct'):.1f}% of squared sensor contrast energy. Projection removed {anchor('anchor_OSD562_to_OSD562_FLT_recovery', 'projection', 'sensor_effect_energy_removed_pct'):.1f}%. These quantities describe coefficient attenuation, not the percentage of technical artifact removed. A sensor participates in learning the score, so attenuation in that sensor is not independent validation.

Transfer was variable. Applying the cerebellum ground score to the hippocampus recovery-exposure contrast removed {anchor('anchor_OSD561_to_OSD562_FLT_recovery', 'score', 'sensor_effect_energy_removed_pct'):.1f}% of sensor energy, but the outside-sensor energy reduction was {anchor('anchor_OSD561_to_OSD562_FLT_recovery', 'score', 'outside_sensor_effect_energy_removed_pct'):.1f}%. Applying the hippocampus anchor to RR-23 reduced sensor energy by {anchor('anchor_OSD562_to_OSD525', 'score', 'sensor_effect_energy_removed_pct'):.1f}% while the outside-sensor reduction was {anchor('anchor_OSD562_to_OSD525', 'score', 'outside_sensor_effect_energy_removed_pct'):.1f}%. Increased energy after transfer is retained and displayed. Anchor bootstrap ranges are broad in several settings and sometimes include increased energy despite a positive point estimate.

The within-RR-23 anchor and FLT-versus-HGC target share nine HGC animals. It is not an independent validation. We therefore also partitioned HGC animals into four for anchor learning and five for evaluation, obtaining all 126 partitions. Across partitions, projection reduced sensor energy by {split_stats.loc['projection', ('sensor_energy_removed_pct', 'mean')]:.1f}% (SD {split_stats.loc['projection', ('sensor_energy_removed_pct', 'std')]:.1f}), and score subtraction by {split_stats.loc['score', ('sensor_energy_removed_pct', 'mean')]:.1f}% (SD {split_stats.loc['score', ('sensor_energy_removed_pct', 'std')]:.1f}). Outside-sensor reductions averaged {split_stats.loc['projection', ('outside_sensor_energy_removed_pct', 'mean')]:.1f}% for projection and {split_stats.loc['score', ('outside_sensor_energy_removed_pct', 'mean')]:.1f}% for score subtraction. The partitions overlap extensively and are a sensitivity distribution from one cohort, not 126 independent replications.

### Estimator consistency and limited inflammatory interpretation

After applying identical animal-level technical-library summation to STAR and RSEM matrices, RR-23 HGC-versus-VGC beta correlations were {float(estimator.iloc[0].sensor_beta_pearson_r):.3f} in the sensor and {float(estimator.iloc[0].genome_beta_pearson_r):.3f} genome-wide. This is consistency of estimators using the same libraries, not independent biological validation.

The inherited glia/complement panel did not meet q < 0.05 in any of the 13 primary unadjusted contrasts; its lowest panel-family q was {glia_min_q:.3f}. This does not exclude inflammatory biology, because bulk tissue composition, small samples, gene-set choice and directional averaging affect sensitivity. It does prevent a strong claim of a consistent neuroinflammatory signature based on this analysis alone. Corrected glial and endothelial panel effects remain exploratory fixed-anchor sensitivities, and expression changes do not establish inflammation, immune-cell activation or blood-brain barrier permeability.

## Discussion

The clearest result is that sample accounting and comparator choice materially constrain interpretation. Shared animals must remain the units of resampling when their regions appear together. Technical libraries must not count as additional animals. Study titles, tissues, assays and location fields must be resolved before combining an apparently broad catalog into evidence for a brain mechanism.

The ground cohorts show a rapid-response difference, but this does not isolate its cause. Handling, anesthetic method, environment, age, tissue composition and collection timing can contribute jointly. In RR-23, FLT, HGC and VGC euthanasia/dissection dates differ, so FLT-versus-VGC cannot be assumed to provide unconfounded biological truth. In RRRM-2, ISS-T and live-animal-return flight cohorts differ in orbital exposure and recovery history: the deposited study describes approximately 55–58 days in orbit for ISS-T and 32 days in orbit followed by 24 days of Earth recovery for live return. Their contrast is not a clean flight indicator or a minutes-to-hours terminal collection clock.

Ground-trained adjustment is attractive because the anchor excludes flight samples, but ground-only estimation does not make a direction artifact-specific. Biological responses that overlap the anchor can be subtracted. Score subtraction is an oblique linear operator; it is not equivalent to orthogonal projection and need not reduce variation outside its sensor. Furthermore, projection of an expression matrix produces projected OLS coefficients, but gene-specific changes in standard errors prevent equivalence to direct t-vector projection. These distinctions matter when reporting statistical significance after adjustment.

Our uncertainty summaries have practical limits. Panel permutation tests rely on within-age cohort-label exchangeability in observational data. Bootstrap intervals condition on normalization; anchor intervals keep the target fixed and do not propagate the entire target, normalization or gene-selection pipeline. Conditional corrected t values do not incorporate anchor estimation uncertainty. We accordingly do not claim calibrated corrected gene-level FDR, verified recovery of true biology or superiority over standard surrogate-variable or control-gene methods [5,6].

Before a strong methods claim, independent collection-controlled samples, repeated-seed null and signal simulations, standard package baselines, and whole-pipeline resampling are needed. For a specific inflammatory mechanism, additional data should distinguish cellular composition from regulation and connect expression to an independent inflammatory phenotype. Direct time-course and isoform evidence would be required for kinetic and transcript-specific interpretations.

The contribution of this revision is a reproducible, animal-aware description of cohort effects and a transparent demonstration of the limits of ground-anchor transfer. It supports testable follow-up hypotheses, while narrowing mechanistic conclusions to the measurements actually available.

## Methods

### Deposits, provenance and inclusion

We used the processed unnormalized count matrices from OSD-525, OSD-561, OSD-562, OSD-563, OSD-564, OSD-612 and OSD-613. Raw reads were not realigned in this revision. Local metadata for all deposits except OSD-525 includes earlier derived cohort fields; the present analysis preserves these records and derives design groups from Source Name identifiers. OSD-525 metadata was extracted from its deposited ISA sample table. Public metadata snapshots were retrieved on 14 October 2026. Original paths, bundled paths, SHA-256 checksums and extraction transformations are recorded in the input provenance manifest. No count-matrix library was removed without metadata. The analysis rejects duplicate identifiers, invalid counts and absent sample metadata.

### Preprocessing and gene panels

All OSD-525 technical libraries were summed by deposited Source Name before normalization, consistently for STAR and RSEM. Gene rows with canonical ENSMUSG identifiers were retained, excluding ERCC and other non-mouse controls. Within each deposit, genes with mean raw count at least ten were retained. Size factors were calculated by the median ratio to the geometric mean among genes positive in every included sample, following the median-ratio normalization principle [4]. The analysis used log2(count/size factor + 1), not a DESeq2 negative-binomial differential-expression fit. Size factors and retained gene counts are exported. Cerebellum and hippocampus were normalized separately and intersected when pooled. Annotation uses mouse GRCm39 Ensembl release 112; missing symbols are documented without silent replacement.

The legacy exploratory panels contain 20 rapid, 14 delayed primary, 14 secondary and 11 heat-shock genes. The sensor is the 31-gene rapid-plus-heat-shock union. Additional inherited glial/complement, endothelial/BBB and circadian panels are explicitly defined in the output tables. These sets were not preregistered for this revision and are not claimed to reproduce the full gene memberships or kinetics of cited experiments. Gene-level Homer1 counts do not resolve its isoforms. A gene can appear in more than one pathway panel; within-set symbols are unique.

### Contrasts, uncertainty and multiplicity

For each gene, ordinary least squares estimated the specified binary contrast with an intercept, an age indicator when applicable and tissue indicators in the pooled model. An SVD basis determines the nuisance rank; unidentifiable contrasts raise errors. Whole-animal labels were permuted within age, with every region of an animal assigned the same label. All possible assignments were enumerated for every reported contrast, including the observed assignment. Two-sided p values equal the fraction whose absolute mean panel coefficient is at least the observed absolute value, retaining ties. The rapid-minus-secondary statistic was tested using the same joint animal permutations; genes were never treated as independent replicate animals.

We evaluated 13 primary contrasts and eight metadata sensitivities. BH adjustment [7] was applied to the family of 168 panel tests and separately to 21 rapid-minus-secondary tests; overlapping sensitivities are not independent replications. Animal-cluster CR1 standard errors with a cluster-count-minus-one t approximation and within-contrast gene BH values are exported for exploratory inspection, not as primary evidence. Panel and gradient intervals use 1,000 percentile bootstrap resamples of animals within label and age, keeping all available regions. Both inference procedures assume cohort exchangeability and do not establish causal identification.

### Ground-trained transformations

Let Y be a genes-by-samples log-expression matrix. Ground samples define a nuisance-residualized matrix Y0 and an OLS cohort coefficient vector bGC. The unit genome vector is v = bGC / norm(bGC). Sensor weights w equal the ground coefficients on Sensor31 and zero elsewhere, normalized to unit length. A ground score is sGC = transpose(w) Y0. Ground loadings are gamma = Y0 sGC / transpose(sGC) sGC.

Genome-expression projection uses Yproj = Y − v transpose(v) Y. Ground-calibrated score subtraction uses Yscore = Y − gamma transpose(w) Y. OLS is refitted with the original target design after each transformation. These operators do not fit a new target sample-score covariate or justify a target N−3 degrees-of-freedom formula. On the aligned gene universe, projection is orthogonal to v, whereas score subtraction removes the w-score through an oblique operator. All target genes are aligned before fitting the anchor, and its genome and sensor normalization are calculated on that aligned universe.

Evaluation reports changes in squared target coefficient energy, separately inside and outside Sensor31. A negative reduction is an increase. These are descriptive diagnostics with no assumed artifact ground truth. The reported lambda diagnostic converts two-sided conditional OLS t probabilities to chi-squared-one quantiles and divides their median by the null median, outside Sensor31. It does not prove calibration. Corrected panel permutation results are explicitly conditional on a fixed anchor and adjusted within their own exploratory family. Shared anchor/target animals are exported. Two hundred ground-animal bootstrap resamples within cohort and age quantify conditional anchor sensitivity, while holding normalized expression and the target fixed.

### Reproducibility and AI assistance

The analysis runs from bundled processed counts, metadata and gene annotation without the original absolute download paths or Paper 2 code. Runtime versions were Python {run['python']}, NumPy {run['numpy']}, pandas {run['pandas']} and SciPy {run['scipy']}; stochastic resampling uses seed 561562. Regression checks compare the implementation with full-design least squares and independently assembled cluster-sandwich covariance, verify block-preserving animal permutations and finite p-value resolution, and distinguish expression projection from t-vector projection. Numerical assertions in this draft are generated from revised CSVs with a machine-readable lookup registry.

An OpenAI Codex assistant was used to audit the prior package, develop this revision's computational scripts, produce figures and assist with drafting. It is not an author. The human authors must validate the code, results, citations and final interpretation and accept responsibility before submission. This disclosure is included because the target journal asks that substantive LLM use be documented in Methods.

## Data availability

Source data are publicly deposited in NASA OSDR under OSD-525, OSD-561, OSD-562, OSD-563, OSD-564, OSD-612 and OSD-613. Accession pages are reachable through the OSDR portal (https://osdr.nasa.gov/bio/repo/). The local revision bundle includes the exact processed input snapshots, original-library and animal-region manifests, metadata disagreements, normalization factors, panel definitions, contrast specifications and generated result tables. A stable public archive identifier for the revision bundle must be assigned before final publication.

## Code availability

Custom Python scripts, dependency requirements, statistical checks and draft-generation code accompany this local revision package. The original-input snapshot step is separate from the portable reanalysis. A public versioned repository or archive DOI remains to be provided before final publication.

## Author contributions

{contribution}

Atish Raza's final supplied role was truncated as “Data”; Data curation has not been inferred without confirmation. The statements reflect the submitting author's supplied information. All-author agreement to author order, contributions and the final manuscript remains to be documented.

## Acknowledgements and funding

We acknowledge the original investigators and NASA OSDR for making the deposited data available. Funding source, grant identifiers and the funder's role in this reanalysis have not yet been supplied and must be confirmed by the authors.

## Ethics and competing interests

This is a secondary analysis of deposited mouse data; no new animal experiments were performed in this revision. Exact source animal-ethics approvals must be extracted from the original study records. No exemption determination or approval identifier is invented. Financial and non-financial competing-interest statements remain to be supplied by all authors.

## References

1. Leek, J. T. et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nature Reviews Genetics 11, 733–739 (2010). https://doi.org/10.1038/nrg2825
2. Tyssowski, K. M. et al. Different Neuronal Activity Patterns Induce Different Gene Expression Programs. Neuron 98, 530–546.e11 (2018). https://doi.org/10.1016/j.neuron.2018.04.001
3. Tullai, J. W. et al. Immediate-Early and Delayed Primary Response Genes Are Distinct in Function and Genomic Architecture. Journal of Biological Chemistry 282, 23981–23995 (2007). https://doi.org/10.1074/jbc.M702044200
4. Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
5. Leek, J. T. & Storey, J. D. Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis. PLoS Genetics 3, e161 (2007). https://doi.org/10.1371/journal.pgen.0030161
6. Risso, D., Ngai, J., Speed, T. P. & Dudoit, S. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology 32, 896–902 (2014). https://doi.org/10.1038/nbt.2931
7. Benjamini, Y. & Hochberg, Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society Series B 57, 289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

## Figure legends

"""
    draft = draft.replace("No count-matrix library was removed without metadata.",
                          "The core local count matrices were independently compared with freshly downloaded NASA STAR files GLDS-556 (OSD-561) and GLDS-557 (OSD-562): all gene IDs, sample IDs and numeric counts were identical. The upstream files, file-listing snapshots, URLs and checksums are included. No count-matrix library was removed without metadata.")
    draft += "\n\n".join(caption for _, caption in captions) + "\n"
    (HERE / "MANUSCRIPT_REVISED.md").write_text(draft)
    pd.DataFrame(CLAIMS).drop_duplicates(subset=["source_file", "filters", "field"]).to_csv(R / "manuscript_numeric_lookups.csv", index=False)

    fontdir = Path(reportlab.__file__).parent / "fonts"
    for name, file in [("Vera", "Vera.ttf"), ("VeraBold", "VeraBd.ttf"), ("VeraItalic", "VeraIt.ttf")]:
        pdfmetrics.registerFont(TTFont(name, str(fontdir / file)))
    pdfmetrics.registerFontFamily("Vera", normal="Vera", bold="VeraBold", italic="VeraItalic", boldItalic="VeraBold")
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="BodyCustom", fontName="Vera", fontSize=9, leading=13,
                              spaceAfter=8, textColor=colors.HexColor("#23323f")))
    styles.add(ParagraphStyle(name="TitleCustom", fontName="VeraBold", fontSize=18, leading=23, spaceAfter=16))
    styles.add(ParagraphStyle(name="H1Custom", fontName="VeraBold", fontSize=13, leading=17, spaceBefore=12, spaceAfter=7, keepWithNext=True))
    styles.add(ParagraphStyle(name="H2Custom", fontName="VeraBold", fontSize=10, leading=14, spaceBefore=9, spaceAfter=6, keepWithNext=True))
    styles.add(ParagraphStyle(name="CaptionCustom", fontName="Vera", fontSize=8, leading=11, spaceAfter=8))

    def para(text, style="BodyCustom"):
        text = html.escape(text)
        text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
        return Paragraph(text, styles[style])

    body = draft.split("## Figure legends")[0]
    story = []
    chunks = body.split("\n\n")
    for chunk in chunks:
        chunk = chunk.strip()
        if not chunk:
            continue
        if chunk.startswith("| Accession"):
            lines = [line for line in chunk.splitlines() if not line.startswith("|---")]
            rows = [[para(cell.strip(), "CaptionCustom") for cell in line.strip("|").split("|")] for line in lines]
            table = Table(rows, colWidths=[58, 57, 129, 107, 100], repeatRows=1)
            table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e6eef3")),
                                       ("VALIGN", (0, 0), (-1, -1), "TOP"),
                                       ("LINEBELOW", (0, 0), (-1, 0), .6, colors.HexColor("#657787")),
                                       ("LINEBELOW", (0, 1), (-1, -1), .25, colors.HexColor("#d7dfe5")),
                                       ("LEFTPADDING", (0, 0), (-1, -1), 5),
                                       ("RIGHTPADDING", (0, 0), (-1, -1), 5)]))
            story.extend([table, Spacer(1, 7)])
        elif chunk.startswith("# "):
            story.append(para(chunk[2:], "TitleCustom"))
        elif chunk.startswith("### "):
            story.append(para(chunk[4:], "H2Custom"))
        elif chunk.startswith("## "):
            story.append(para(chunk[3:], "H1Custom"))
        elif re.match(r"^1\. Leek", chunk):
            story.extend(para(line) for line in chunk.splitlines() if line.strip())
        else:
            story.append(para(chunk))
    for name, caption in captions:
        story.append(PageBreak())
        image_path = figures / f"{name}.png"
        from PIL import Image as PILImage
        with PILImage.open(image_path) as im:
            w, h = im.size
        width = 451
        height = width * h / w
        if height > 570:
            width *= 570 / height
            height = 570
        story.append(KeepTogether([Image(str(image_path), width=width, height=height),
                                   Spacer(1, 12), para(caption, "CaptionCustom")]))

    def footer(canvas, document):
        canvas.saveState()
        canvas.setFont("Vera", 7)
        canvas.setFillColor(colors.HexColor("#657787"))
        canvas.drawString(48, 24, "Author-review revision • 14 October 2026 • unresolved items remain")
        canvas.drawRightString(A4[0] - 48, 24, str(document.page))
        canvas.restoreState()

    SimpleDocTemplate(str(HERE / "MANUSCRIPT_REVISED.pdf"), pagesize=A4, rightMargin=48, leftMargin=48,
                      topMargin=44, bottomMargin=42, title="Cohort-associated transcriptional variation in spaceflight mouse brain",
                      author=names).build(story, onFirstPage=footer, onLaterPages=footer)
    (R / "figure_manifest.json").write_text(json.dumps([
        {"filename": name, "caption": caption, "png_dpi": 600, "vector_format": "SVG"}
        for name, caption in captions], indent=2) + "\n")
    print(f"Built revised manuscript, {len(CLAIMS)} numerical lookups and {len(captions)} figures")


if __name__ == "__main__":
    main()
