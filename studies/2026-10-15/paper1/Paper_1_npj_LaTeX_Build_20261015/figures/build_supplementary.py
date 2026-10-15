"""Generate Supplementary Information LaTeX for Paper 1 from result tables (no typed numbers)."""
import json
import re
from pathlib import Path

import pandas as pd

BUILD = Path(__file__).resolve().parents[1]
PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
FIGSRC = PKG / "publication_figures"
OUT = BUILD / "src/supplementary.tex"

main_md = (PKG / "MANUSCRIPT_COMPLETE.md").read_text().split("## Figures", 1)[1]
si_md = (PKG / "SUPPLEMENTARY_INFORMATION.md").read_text()
leg = {f"F{m.group(1)}": m.group(2).strip() for m in re.finditer(r"^Figure (\d+)\. (.+)$", main_md, flags=re.M)}
leg.update({f"S{m.group(1)}": m.group(2).strip()
            for m in re.finditer(r"^Supplementary Figure (\d+)\. (.+)$", si_md, flags=re.M)})

FIGS = [  # (file stem, legend key, short title used in the new paper)
    ("figure1_sample_accounting_and_panel_effects", "F1"),
    ("figure2_cross_cohort_panel_effects", "F2"),
    ("figureS2_all_deposit_structure", "S2"),
    ("figure3_count_model_and_biology", "F3"),
    ("figureS3_flight_collection_interaction", "S3"),
    ("figureS4_fixed_blood_marker_context", "S4"),
    ("figure4_paired_regional_interactions", "F4"),
    ("figureS1_cross_cohort_gene_screen", "S1"),
    ("figureS6_biological_robustness", "S6"),
    ("figureS7_external_biological_context", "S7"),
    ("figure6_ground_holdout", "F6"),
    ("figureS5_independent_ground_and_age_validation", "S5"),
    ("figure5_bounded_operator_and_joint_uncertainty", "F5"),
    ("figure9_ground_anchor_transfer", "F9"),
    ("figure10_disjoint_anchor_partitions", "F10"),
    ("figure11_known_truth_validation", "F11"),
    ("figure12_high_gain_domain_shift", "F12"),
]


def tex(s):
    s = s.replace("\\", "\\textbackslash{}")
    for a, b in (("&", "\\&"), ("%", "\\%"), ("_", "\\_"), ("#", "\\#"), ("$", "\\$")):
        s = s.replace(a, b)
    rep = {"≤": "$\\le$", "≥": "$\\ge$", "×": "$\\times$", "−": "$-$", "–": "--", "—": "---", "±": "$\\pm$",
           "²": "$^2$", "₂": "$_2$", "⁻": "$^-$", "¹": "$^1$", "⁰": "$^0$", "⁵": "$^5$", "⁷": "$^7$", "α": "$\\alpha$",
           "β": "$\\beta$", "γ": "$\\gamma$", "λ": "$\\lambda$", "χ": "$\\chi$", "‖": "$\\Vert$", "’": "'", "‘": "'",
           "“": "``", "”": "''", "→": "$\\rightarrow$", "·": "$\\cdot$", "µ": "$\\mu$", "Σ": "$\\Sigma$", "ᵀ": "$^{\\top}$"}
    for a, b in rep.items():
        s = s.replace(a, b)
    # format panel letters (a), (b,c) in bold as the journal requires
    s = re.sub(r"\(([a-h](?:[,–-][a-h])*)\)", lambda m: "\\textbf{" + m.group(1) + "}", s)
    return s


def short(status):
    return {"ALIASED (disjoint values)": "A", "partially overlapping (unbalanced)": "P", "matched": "M",
            "not deposited (cannot assess)": "--", "same levels in both groups (balanced/blocked or adjusted)": "B"}[status]


lines = [r"""\documentclass[10pt]{article}
\usepackage[T1]{fontenc}\usepackage{tgtermes}\usepackage{tgheros}\usepackage{newtxmath}
\usepackage[paperwidth=210mm,paperheight=279mm,left=14mm,right=14mm,top=18mm,bottom=18mm]{geometry}
\usepackage{graphicx}\usepackage{booktabs}\usepackage{tabularx}\usepackage{caption}\usepackage{xcolor}
\usepackage{fancyhdr}\usepackage[colorlinks=true,allcolors=blue!50!black]{hyperref}
\graphicspath{{../figures/supplementary/}}
\DeclareCaptionLabelSeparator{npjbar}{\space\textbar\space}
\captionsetup[figure]{name=Supplementary Fig.,labelfont=bf,labelsep=npjbar,font={small},justification=justified,singlelinecheck=false}
\captionsetup[table]{name=Supplementary Table,labelfont=bf,labelsep=npjbar,font={small,sf},justification=raggedright,singlelinecheck=false}
\pagestyle{fancy}\fancyhf{}\fancyhead[L]{\sffamily\small Supplementary Information \textendash{} author review version}
\fancyhead[R]{\sffamily\small\bfseries npj Microgravity style}\fancyfoot[R]{\sffamily\small\bfseries\thepage}
\renewcommand{\headrulewidth}{0.34pt}
\begin{document}
\setlength{\parindent}{0pt}
{\sffamily\bfseries\LARGE Supplementary Information\par}\vskip 6pt
{\sffamily\large Collection protocol outweighs spaceflight effects on brain immediate-early genes and resists post hoc correction\par}\vskip 6pt
{\sffamily\small Isarar Siddique, Raiyan Siddique, Arshi Niyaz, Jamyang Ugyen Tshomo \& Atish Raza. Convolity AI Lab, Lucknow, India.\par}
\vskip 10pt
{\small This file contains Supplementary Tables 1\textendash3, built directly from the result tables of the protocol decomposition and the design audit, and Supplementary Figs. 1\textendash17, which reproduce the validated figures and legends of the full analysis. Numbers in the legends refer to the complete results package; overlapping holdouts, bootstraps and simulations are not independent animals.\par}
"""]

# Supplementary Table 1: protocol factors per contrast
a = pd.read_csv(BUILD / "research/A_audit/A_contrast_aliasing.csv")
factors = [("euthanasia_method", "Agent"), ("design_euthanasia_location", "Kill location"),
           ("preservation_class", "Bulk freezing"), ("carcass_preservation", "Carcass"),
           ("euthanasia_date", "Date"), ("recovery_on_earth_days", "Recovery"), ("habitat", "Habitat"),
           ("age", "Age")]
piv = a.pivot_table(index="contrast", columns="factor", values="status", aggfunc="first")
order = list(dict.fromkeys(a.contrast))
lines.append(r"\begin{table}[h]\centering\caption{Protocol factors that each analysed contrast matches (M), aliases (A), only partly overlaps (P), balances or adjusts (B), or that were not deposited (--)}\label{stab:alias}")
lines.append(r"\sffamily\scriptsize\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}X" + "c" * len(factors) + "@{}}\\toprule")
lines.append("Contrast & " + " & ".join(f for _, f in factors) + r" \\\midrule")
for c in order:
    row = [tex(c)] + [short(piv.loc[c, f]) if f in piv.columns and pd.notna(piv.loc[c, f]) else "" for f, _ in factors]
    lines.append(" & ".join(row) + r" \\")
lines.append(r"\bottomrule\end{tabularx}\par\medskip\raggedright\scriptsize Derived factor by factor from the deposited ISA fields; the complete table with the values in each group is \texttt{A\_contrast\_aliasing.csv}. Cb, cerebellum; Hip, hippocampus.\end{table}")

# Supplementary Table 2: dominance tests
d = pd.read_csv(BUILD / "research/D_protocol/D_tables/D2_dominance_GC_ISST_vs_FLT_LAR.csv")
lines.append(r"\begin{table}[h]\centering\caption{Dominance test: ground ISS-T minus flight LAR, by panel and region}\label{stab:dom}")
lines.append(r"\sffamily\scriptsize\begin{tabular}{@{}l l r r r r@{}}\toprule Region & Panel & Effect (log$_2$) & Exact $p$ & BH $q$ (24 tests) & Assignments \\\midrule")
for _, r in d.iterrows():
    lines.append(f"{r.region} & {tex(r.panel)} & {r.GC_ISST_minus_FLT_LAR:.3f} & {r.p_exact:.4f} & {r.q_BH_family24:.4f} & {int(r.n_assignments)} \\\\")
lines.append(r"\bottomrule\end{tabular}\par\medskip\raggedright\scriptsize Positive values mean that the protocol effect in ground mice exceeds the matched flight effect, both measured against the shared ground LAR reference. Exact enumeration of whole-animal labels within age blocks.\end{table}")

# Supplementary Table 3: 2x2 decomposition
s = pd.read_csv(BUILD / "research/D_protocol/D_tables/D3_two_by_two_bootstrap.csv")
lines.append(r"\clearpage\begin{table}[h]\centering\caption{RRRM-2 flight-by-protocol decomposition for all eight panels}\label{stab:2x2}")
lines.append(r"\sffamily\scriptsize\begin{tabular}{@{}l l r r r r r@{}}\toprule Region & Panel & $P_G$ [95\% CI] & $F_L$ [95\% CI] & $F_I$ & $P_F$ & $P_G-F_L$ [95\% CI] \\\midrule")
for _, r in s.iterrows():
    lines.append(f"{r.region} & {tex(r.panel)} & {r.P_G:.3f} [{r.P_G_lo95:.2f}, {r.P_G_hi95:.2f}] & {r.F_L:.3f} [{r.F_L_lo95:.2f}, {r.F_L_hi95:.2f}] & {r.F_I:.3f} & {r.P_F:.3f} & {r.delta_PG_minus_FL:.3f} [{r.delta_lo95:.2f}, {r.delta_hi95:.2f}] \\\\")
lines.append(r"\bottomrule\end{tabular}\par\medskip\raggedright\scriptsize $P_G$, ground ISS-T minus ground LAR; $F_L$, flight minus ground within LAR; $F_I$, flight minus ground within ISS-T; $P_F$, flight ISS-T minus flight LAR. Intervals from 2,000 whole-animal bootstraps stratified by flight, cohort and age (seed 15102026).\end{table}")

# Supplementary figures
for k, (stem, key) in enumerate(FIGS, 1):
    legend = leg[key]
    title, _, body = legend.partition(". ")
    lines.append(r"\clearpage\begin{figure}[p]\centering")
    lines.append(r"\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{" + stem + ".png}")
    lines.append(r"\caption{\textbf{" + tex(title) + ".} " + tex(body) + "}\\label{sfig:" + str(k) + "}\\end{figure}")
lines.append(r"\end{document}")
OUT.write_text("\n".join(lines))
json.dump({str(k): {"file": stem, "legend_source": key} for k, (stem, key) in enumerate(FIGS, 1)},
          open(BUILD / "src/supplementary_figure_map.json", "w"), indent=1)
print("wrote", OUT, "with", len(FIGS), "figures")
