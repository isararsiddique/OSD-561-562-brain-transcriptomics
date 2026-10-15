"""Build Paper 1 main figures at npj Microgravity print size (180 mm double column).

Style: Arial, >= 7 pt at print size, Okabe-Ito colours with redundant marker shape / fill / line
style, rounded rectangles for diagrams. Every number is read from a results table; nothing is
typed in. Usage: python build_figures.py fig2 [fig1 fig3 ...]
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

BUILD = Path(__file__).resolve().parents[1]
PKG = Path(""
           "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
D = BUILD / "research/D_protocol"
B = BUILD / "research/B_novelty/analysis"
OUT = BUILD / "figures"
MM = 1 / 25.4

OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
      "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7",
      "grey": "#7F7F7F"}
PALETTE = set(OI.values()) | {"#FFFFFF", "#F2F2F2", "#E6E6E6", "#BFBFBF", "#595959", "#404040",
                              "#FBE6C2", "#D9EEFB", "#CCEBE1", "#F6D9C7", "#EEDCE8", "#FFFFFF00"}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 7.5, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7.5, "legend.fontsize": 7, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "lines.linewidth": 1.0, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 600,
})


def letter(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    print("wrote", name, flush=True)


def rbox(ax, x, y, w, h, text, fc, ec="#404040", fs=6.6, weight="normal", lw=0.7, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
                                fc=fc, ec=ec, lw=lw, ls=ls, transform=ax.transAxes, clip_on=False))
    ax.text(x + w / 2, y + h / 2, text, transform=ax.transAxes, ha="center", va="center",
            fontsize=fs, fontweight=weight, linespacing=1.25)


def arrow(ax, p0, p1, color, ls="-", lw=1.1, label=None, lx=0, ly=0, ha="center", fs=6.6):
    a = FancyArrowPatch(p0, p1, transform=ax.transAxes, arrowstyle="-|>", mutation_scale=7,
                        color=color, lw=lw, ls=ls, shrinkA=0, shrinkB=0, clip_on=False)
    ax.add_patch(a)
    if label:
        mx, my = (p0[0] + p1[0]) / 2 + lx, (p0[1] + p1[1]) / 2 + ly
        ax.text(mx, my, label, transform=ax.transAxes, ha=ha, va="center", fontsize=fs, color="#000000")


# ----------------------------------------------------------------------------------- Figure 2
def fig2():
    tbl = pd.read_csv(D / "D_tables/D3_two_by_two_bootstrap.csv")
    dom = pd.read_csv(D / "D_tables/D2_dominance_GC_ISST_vs_FLT_LAR.csv")
    genes = pd.read_csv(D / "D_tables/D4_gene_level_effects.csv")
    cm = pd.read_csv(D / "D_tables/D5_cross_mission_by_protocol_class.csv")

    from matplotlib.lines import Line2D
    fig = plt.figure(figsize=(180 * MM, 196 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.92, 1.08], height_ratios=[0.86, 1.14],
                          left=0.085, right=0.985, top=0.975, bottom=0.055, wspace=0.42, hspace=0.24)

    # (a) design schematic --------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0]); ax.set_axis_off()
    ax.text(-0.09, 1.0, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    colx = {"ISS_T": 0.08, "LAR": 0.58}
    rowy = {"FLT": 0.56, "GC": 0.15}
    w, h = 0.40, 0.29
    fills = {"ISS_T": "#FBE6C2", "LAR": "#D9EEFB"}
    txt = {("FLT", "ISS_T"): "Killed on orbit\nketamine/xylazine\ncryochiller\nn = 9 Cb, 9 Hc",
           ("GC", "ISS_T"): "Killed on Earth\nketamine/xylazine\nliquid nitrogen\nn = 6 Cb, 6 Hc",
           ("FLT", "LAR"): "24 d Earth recovery\ngas anaesthesia\nliquid nitrogen\nn = 5 Cb, 6 Hc",
           ("GC", "LAR"): "Killed on Earth\ngas anaesthesia\nliquid nitrogen\nn = 6 Cb, 6 Hc"}
    for (fl, co), t in txt.items():
        rbox(ax, colx[co], rowy[fl], w, h, t, fills[co], fs=7)
    for co, t in (("ISS_T", "ISS-terminal (ISS-T)"), ("LAR", "Live animal return (LAR)")):
        ax.text(colx[co] + w / 2, rowy["FLT"] + h + 0.035, t, transform=ax.transAxes, ha="center",
                fontsize=7.5, fontweight="bold")
    for fl, t in (("FLT", "Flight"), ("GC", "Ground")):
        ax.text(colx["ISS_T"] - 0.03, rowy[fl] + h / 2, t, transform=ax.transAxes, ha="right",
                va="center", fontsize=7.5, fontweight="bold", rotation=90)
    gx0, gx1 = colx["ISS_T"] + w, colx["LAR"]
    arrow(ax, (gx1, rowy["GC"] + h / 2), (gx0, rowy["GC"] + h / 2), OI["vermillion"], lw=1.2,
          label=r"$P_G$", ly=0.045, fs=7.5)
    arrow(ax, (gx1, rowy["FLT"] + h / 2), (gx0, rowy["FLT"] + h / 2), OI["vermillion"], lw=1.2,
          ls=(0, (1.5, 1.2)), label=r"$P_F$", ly=0.045, fs=7.5)
    gy0, gy1 = rowy["GC"] + h, rowy["FLT"]
    arrow(ax, (colx["LAR"] + w / 2, gy0), (colx["LAR"] + w / 2, gy1), OI["blue"], lw=1.2,
          ls=(0, (3.5, 1.5)), label=r"$F_L$", lx=0.055, fs=7.5)
    arrow(ax, (colx["ISS_T"] + w / 2, gy0), (colx["ISS_T"] + w / 2, gy1), OI["blue"], lw=1.2,
          ls=(0, (1.5, 1.2)), label=r"$F_I$", lx=0.055, fs=7.5)
    keys = [("o", True, OI["vermillion"], r"$P_G$  protocol, no flight"),
            ("o", False, OI["vermillion"], r"$P_F$  protocol, all flown"),
            ("s", False, OI["blue"], r"$F_I$  flight, freezing differs"),
            ("s", True, OI["blue"], r"$F_L$  flight, protocol matched")]
    h_keys = [Line2D([], [], ls="-", lw=0.9, color=c, marker=m, ms=4, mfc=c if f else "white", mec=c,
                     label=l) for m, f, c, l in keys]
    ax.legend(handles=h_keys, loc="upper left", bbox_to_anchor=(-0.02, 0.12), ncol=2, frameon=False,
              fontsize=7, handlelength=1.3, columnspacing=0.8, handletextpad=0.4, borderaxespad=0)

    # (b) forest of panel effects ---------------------------------------------------------
    sub = gs[0, 1].subgridspec(1, 2, wspace=0.10)
    panels = ["rPRG", "dPRG", "HSP", "SRG", "Circadian"]
    spec = [("P_G", OI["vermillion"], "o", True), ("P_F", OI["vermillion"], "o", False),
            ("F_I", OI["blue"], "s", False), ("F_L", OI["blue"], "s", True)]
    offs = np.linspace(0.28, -0.28, 4)
    for k, region in enumerate(["cerebellum", "hippocampus"]):
        axb = fig.add_subplot(sub[0, k])
        t = tbl[tbl.region == region].set_index("panel")
        dm = dom[dom.region == region].set_index("panel")
        for i, p in enumerate(panels):
            yb = len(panels) - 1 - i
            if i % 2 == 0:
                axb.axhspan(yb - 0.45, yb + 0.45, color="#F2F2F2", lw=0, zorder=0)
            for (key, col, mk, filled), off in zip(spec, offs):
                v, lo, hi = t.loc[p, key], t.loc[p, f"{key}_lo95"], t.loc[p, f"{key}_hi95"]
                axb.plot([lo, hi], [yb + off] * 2, color=col, lw=0.9, zorder=2)
                axb.plot(v, yb + off, mk, ms=3.8, mfc=col if filled else "white", mec=col, mew=0.9, zorder=3)
            if dm.loc[p, "q_BH_family24"] <= 0.05:
                axb.text(1.70, yb, r"$\ddagger$", ha="right", va="center", fontsize=8)
        axb.axvline(0, color="#595959", lw=0.6, ls=(0, (2, 2)), zorder=1)
        axb.set_ylim(-0.55, len(panels) - 0.45)
        axb.set_xlim(-0.95, 1.75)
        axb.set_xticks([-0.5, 0, 0.5, 1, 1.5])
        axb.set_yticks(range(len(panels)))
        axb.set_yticklabels(panels[::-1] if k == 0 else [])
        axb.set_xlabel("Panel effect (log$_2$)")
        axb.set_title(region.capitalize(), fontsize=7.5, fontweight="bold", pad=3)
        if k == 0:
            axb.text(-0.42, 1.0, "b", transform=axb.transAxes, fontsize=9, fontweight="bold", va="bottom")

    # (c) gene-level scatter --------------------------------------------------------------
    axc = fig.add_subplot(gs[1, 0])
    axc.text(-0.2, 1.0, "c", transform=axc.transAxes, fontsize=9, fontweight="bold", va="bottom")
    pcol = {"rPRG": OI["vermillion"], "dPRG": OI["orange"], "HSP": OI["green"], "SRG": OI["blue"]}
    rmk = {"cerebellum": "o", "hippocampus": "^"}
    for region in ("cerebellum", "hippocampus"):
        s = genes[genes.region == region]
        for p, col in pcol.items():
            ss = s[s.panels.str.split("|").str[0] == p]
            axc.scatter(ss.P_G, ss.F_L, s=13, marker=rmk[region],
                        facecolor=col if region == "hippocampus" else "white", edgecolor=col, lw=0.8, zorder=3)
    lim = (-1.6, 2.7)
    axc.plot(lim, lim, color="#595959", lw=0.7, ls=(0, (3, 2)))
    axc.axhline(0, color="#BFBFBF", lw=0.6); axc.axvline(0, color="#BFBFBF", lw=0.6)
    axc.set_xlim(lim); axc.set_ylim(-1.6, 2.1)
    axc.set_xlabel(r"Protocol effect in ground mice, $P_G$ (log$_2$)")
    axc.set_ylabel(r"Matched flight effect, $F_L$ (log$_2$)")
    for region, sym, dx, dy in (("cerebellum", "Npas4", -0.62, 0.25), ("cerebellum", "Fos", -0.62, 0.36),
                                ("hippocampus", "Fos", -0.05, -0.62), ("hippocampus", "Egr1", 0.05, 0.55)):
        r = genes[(genes.region == region) & (genes.symbol == sym)].iloc[0]
        axc.annotate(sym, (r.P_G, r.F_L), xytext=(r.P_G + dx, r.F_L + dy), fontsize=7, style="italic",
                     arrowprops=dict(arrowstyle="-", lw=0.5, color="#595959"))
    n_cb = genes[genes.region == "cerebellum"]; n_hc = genes[genes.region == "hippocampus"]
    axc.text(0.03, 0.985, f"|$P_G$| > |$F_L$| in {int((n_cb.P_G.abs() > n_cb.F_L.abs()).sum())}/{len(n_cb)}"
             f" cerebellar and\n{int((n_hc.P_G.abs() > n_hc.F_L.abs()).sum())}/{len(n_hc)} hippocampal panel genes",
             transform=axc.transAxes, va="top", fontsize=7)
    hh = [Line2D([], [], ls="", marker="s", ms=4, mfc=c, mec=c, label=p) for p, c in pcol.items()]
    hh += [Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec="#404040", label="Cerebellum"),
           Line2D([], [], ls="", marker="^", ms=4, mfc="#404040", mec="#404040", label="Hippocampus")]
    axc.legend(handles=hh, loc="lower right", frameon=False, fontsize=7, ncol=2, handletextpad=0.2,
               columnspacing=0.6, borderaxespad=0.1)

    # (d) cross-mission forest by protocol class ----------------------------------------------
    axd = fig.add_subplot(gs[1, 1])
    axd.text(-0.47, 1.0, "d", transform=axd.transAxes, fontsize=9, fontweight="bold", va="bottom")
    order = ["Protocol differs, no flight", "Protocol differs, all flown",
             "Flight, euthanasia matched, freezing differs", "Flight, protocol matched",
             "Housing only, protocol matched"]
    ccol = {order[0]: OI["vermillion"], order[1]: OI["orange"], order[2]: OI["purple"],
            order[3]: OI["blue"], order[4]: OI["green"]}
    rows, y = [], 0
    ylabels, ypos, headers = [], [], []
    for c in order:
        sub_c = cm[(cm["class"] == c) & (cm.panel == "rPRG")]
        headers.append((y, c)); y += 1
        for _, r in sub_c.iterrows():
            rows.append((y, r)); ylabels.append(r.label); ypos.append(y); y += 1
        y += 0.3
    hsp = cm[cm.panel == "HSP"].set_index("contrast")
    for yy, r in rows:
        col = ccol[r["class"]]
        sig = r.q_BH_232 <= 0.05
        axd.plot([r.lo95, r.hi95], [yy - 0.15] * 2, color=col, lw=1.0)
        axd.plot(r.effect, yy - 0.15, "o", ms=4.0, mfc=col if sig else "white", mec=col, mew=0.9)
        hr = hsp.loc[r.contrast]
        axd.plot([hr.lo95, hr.hi95], [yy + 0.2] * 2, color=col, lw=0.8, ls=(0, (1.5, 1)))
        axd.plot(hr.effect, yy + 0.2, "D", ms=3.0, mfc=col if hr.q_BH_232 <= 0.05 else "white", mec=col, mew=0.8)
    for yy, c in headers:
        short = {order[0]: "Protocol differs, no flight", order[1]: "Protocol differs, all flown",
                 order[2]: "Flight; freezing differs", order[3]: "Flight; protocol matched",
                 order[4]: "Housing only; matched"}[c]
        axd.text(-0.02, yy, short, transform=axd.get_yaxis_transform(), fontsize=7, fontweight="bold",
                 color=ccol[c], va="center", ha="right")
    axd.axvline(0, color="#595959", lw=0.6, ls=(0, (2, 2)))
    axd.set_yticks(ypos); axd.set_yticklabels(ylabels, fontsize=7)
    axd.tick_params(axis="y", length=0)
    axd.set_ylim(y - 0.1, -0.7)
    axd.set_xlim(-1.35, 2.25)
    axd.set_xticks([-1, -0.5, 0, 0.5, 1, 1.5])
    axd.set_xlabel("rPRG / HSP effect (log$_2$), 95% animal bootstrap")
    hd = [Line2D([], [], ls="-", lw=1.0, marker="o", ms=4, mfc="#404040", mec="#404040", color="#404040",
                 label="rPRG"),
          Line2D([], [], ls=(0, (1.5, 1)), lw=0.8, marker="D", ms=3, mfc="#404040", mec="#404040",
                 color="#404040", label="HSP"),
          Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec="#404040", label="open: BH q > 0.05")]
    axd.legend(handles=hd, loc="center right", bbox_to_anchor=(1.0, 0.565), frameon=False, fontsize=7,
               handlelength=1.5, borderaxespad=0)
    save(fig, "Figure_2")


# ----------------------------------------------------------------------------------- Figure 5
def fig5():
    from matplotlib.lines import Line2D
    hold = pd.read_csv(PKG / "extensions/ground_holdout/holdout_split_results.csv")
    geo = pd.read_csv(PKG / "results/operator_geometry.csv")
    nc = json.loads((B / "noise_ceiling_results.json").read_text())
    nc525 = json.loads((B / "noise_ceiling_OSD525.json").read_text())
    settings = [("OSD-561", "RRRM-2\ncerebellum"), ("OSD-562", "RRRM-2\nhippocampus"), ("OSD-525", "RR-23\ncerebellum")]
    scol = {"OSD-561": OI["vermillion"], "OSD-562": OI["blue"], "OSD-525": OI["green"]}
    smk = {"OSD-561": "o", "OSD-562": "^", "OSD-525": "s"}
    rng = np.random.default_rng(5)

    fig = plt.figure(figsize=(180 * MM, 140 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.85, 1.15], left=0.075, right=0.985, top=0.93, bottom=0.08,
                          wspace=0.30, hspace=0.55)

    # (a) oblique projection geometry -----------------------------------------------------
    ax = fig.add_subplot(gs[0, 0]); ax.set_aspect("equal")
    ax.text(-0.12, 1.0, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    w = np.array([1.0, 0.0]); hvec = np.array([0.0, 1.6]); gam = w + hvec; ng = np.linalg.norm(gam)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th), np.sin(th), color="#7F7F7F", lw=0.9, ls=(0, (2, 1.5)))
    ax.plot([0, 0], [-2.25, 2.25], color="#BFBFBF", lw=0.6)
    ax.plot([0.0, 0.0], [-ng, ng], color=OI["vermillion"], lw=3.0, solid_capstyle="butt", alpha=0.9)
    ax.plot([0.0, 0.0], [-1, 1], color="#404040", lw=1.6, solid_capstyle="butt")
    for vec, col, lab, off in ((w, OI["blue"], r"$w$", (0.06, -0.25)),
                               (gam, OI["black"], r"$\gamma=w+h$", (0.07, 0.02))):
        ax.annotate("", xy=vec, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", lw=1.1, color=col,
                                                               mutation_scale=8))
        ax.text(vec[0] + off[0], vec[1] + off[1], lab, fontsize=7.5, color=col)
    xs = np.array([-1.6, 1.0]) / ng
    ax.plot(*xs, "o", ms=3.5, color="#000000")
    ax.annotate("", xy=(0, ng), xytext=xs, arrowprops=dict(arrowstyle="-|>", lw=0.9, ls=(0, (2, 1.2)),
                                                           color=OI["vermillion"], mutation_scale=7))
    ax.text(xs[0] - 0.08, xs[1] + 0.05, "unit vector", fontsize=7, ha="right", va="bottom")
    ax.text(0.1, ng + 0.05, r"$\|A_s x\| = \|\gamma\|$", fontsize=7.5, color=OI["vermillion"], va="bottom")
    ax.set_xlim(-2.25, 2.25); ax.set_ylim(-2.25, 2.45)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0.15, -1.08, "orthogonal projection:\nhalf-length 1", fontsize=7, color="#404040", ha="left", va="top")
    ax.text(-0.15, -1.3, "score removal\n" + r"$A_s=I-\gamma w^{\top}$:" + "\nhalf-length " + r"$\|\gamma\|_2$",
            fontsize=7, color=OI["vermillion"], va="top", ha="right")

    # (b) train-test loading cosine versus the identical-loading noise ceiling ------------------
    axb = fig.add_subplot(gs[0, 1])
    axb.text(-0.12, 1.0, "b", transform=axb.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for i, (acc, name) in enumerate(settings):
        raw = hold[(hold.accession == acc) & (hold.method == "raw")]
        for j, (col_name, kind) in enumerate((("train_test_contrast_cosine", "whole_genome"),
                                              ("train_test_sensor_contrast_cosine", "Sensor31"))):
            x0 = i + (-0.17 if j == 0 else 0.17)
            v = raw[col_name].to_numpy()
            axb.scatter(x0 + rng.uniform(-0.06, 0.06, len(v)), v, s=4, color=scol[acc], alpha=0.45, lw=0,
                        marker=smk[acc], zorder=2)
            axb.plot([x0 - 0.09, x0 + 0.09], [v.mean()] * 2, color="#000000", lw=1.3, zorder=4)
            if acc == "OSD-525":
                exp_mean, exp_lo, exp_hi = (nc525["expected_cosine_if_identical_loading_analytic"], None, None) \
                    if kind == "whole_genome" else (None, None, None)
            else:
                e = nc[acc][kind]
                exp_mean = e["expected_cosine_if_identical_loading_analytic"]
                exp_lo, exp_hi = e["expected_cosine_sim_q025_q975"]
            if exp_mean is not None:
                if exp_lo is not None:
                    axb.add_patch(plt.Rectangle((x0 - 0.11, exp_lo), 0.22, exp_hi - exp_lo, fc="#E6E6E6",
                                                ec="#7F7F7F", lw=0.6, zorder=1))
                axb.plot([x0 - 0.11, x0 + 0.11], [exp_mean] * 2, color="#7F7F7F", lw=1.0, ls=(0, (2, 1)), zorder=3)
            auc = raw.sensor_age_stratified_AUC.mean()
        axb.text(i, 1.03, name.replace("\n", " ") + f"\nsensor AUC {auc:.2f}", ha="center", va="bottom", fontsize=7)
    axb.set_xticks([i + d for i in range(3) for d in (-0.17, 0.17)])
    axb.set_xticklabels(["genome", "sensor"] * 3, fontsize=7)
    axb.set_ylim(-0.25, 1.08); axb.set_xlim(-0.5, 2.5)
    axb.set_ylabel("Train-test loading cosine")
    hb = [Line2D([], [], ls="-", lw=1.3, color="#000000", label="observed mean"),
          Line2D([], [], ls=(0, (2, 1)), lw=1.0, color="#7F7F7F", label="expected if loading identical"),
          plt.Rectangle((0, 0), 1, 1, fc="#E6E6E6", ec="#7F7F7F", lw=0.6, label="95% simulated range")]
    axb.legend(handles=hb, loc="lower center", bbox_to_anchor=(0.44, -0.33), ncol=3, frameon=False, fontsize=7,
               handlelength=1.3, borderaxespad=0.0, columnspacing=0.7, handletextpad=0.4)

    # (c) operator gain -------------------------------------------------------------------
    axc = fig.add_subplot(gs[1, 0])
    axc.text(-0.2, 1.0, "c", transform=axc.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for i, (acc, name) in enumerate(settings):
        sc = hold[(hold.accession == acc) & (hold.method == "score")]
        v = sc.operator_gain.to_numpy()
        axc.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=4, color=scol[acc], alpha=0.45, lw=0,
                    marker=smk[acc])
        axc.plot([i - 0.22, i + 0.22], [v.mean()] * 2, color="#000000", lw=1.3)
        g = geo[geo.anchor_accession == acc].score_spectral_gain.to_numpy()
        axc.scatter(np.full(len(g), i + 0.33), g, s=12, marker="D", facecolor="white", edgecolor="#000000",
                    lw=0.7, zorder=4)
    axc.axhline(2, color=OI["orange"], lw=1.2, ls=(0, (4, 1.5)))
    axc.axhline(1, color="#595959", lw=0.8, ls=(0, (1, 1)))
    axc.set_xticks(range(3)); axc.set_xticklabels([n for _, n in settings], fontsize=7)
    axc.set_ylim(0, 30); axc.set_xlim(-0.5, 2.6)
    axc.set_ylabel(r"Spectral gain $\|A\|_2$ of score removal")
    hc = [Line2D([], [], ls="", marker="o", ms=3, color="#7F7F7F", label="held-out split"),
          Line2D([], [], ls="", marker="D", ms=4, mfc="white", mec="#000000", label="full-data geometry"),
          Line2D([], [], ls=(0, (4, 1.5)), lw=1.2, color=OI["orange"], label="gain cap L = 2"),
          Line2D([], [], ls=(0, (1, 1)), lw=0.8, color="#595959", label="orthogonal (gain 1)")]
    axc.legend(handles=hc, loc="upper right", frameon=False, fontsize=7, handletextpad=0.3, borderaxespad=0.2)

    # (d) held-out prediction versus genome-wide loading SNR ------------------------------------
    axd = fig.add_subplot(gs[1, 1])
    axd.text(-0.12, 1.0, "d", transform=axd.transAxes, fontsize=9, fontweight="bold", va="bottom")
    snr = {"OSD-561": nc["OSD-561"]["whole_genome"]["signal_energy_S"] / nc["OSD-561"]["whole_genome"]["noise_energy_train"],
           "OSD-562": nc["OSD-562"]["whole_genome"]["signal_energy_S"] / nc["OSD-562"]["whole_genome"]["noise_energy_train"],
           "OSD-525": nc525["signal_energy_S"] / nc525["noise_energy_train"]}
    for acc, name in settings:
        for meth, mk, filled, dx in (("score", smk[acc], True, 0.0), ("bounded_score", smk[acc], False, 0.0)):
            v = hold[(hold.accession == acc) & (hold.method == meth)].outside_prediction_MSE_ratio_to_nuisance_only
            x = snr[acc] * np.exp(rng.uniform(-0.07, 0.07, len(v)) + (0.16 if meth == "bounded_score" else -0.16))
            axd.scatter(x, v, s=4, color=scol[acc], alpha=0.35 if filled else 0.25, lw=0, marker=mk)
            xm = snr[acc] * np.exp(0.16 if meth == "bounded_score" else -0.16)
            axd.plot(xm, v.mean(), marker=mk, ms=5.5, mfc=scol[acc] if filled else "white", mec="#000000",
                     mew=0.8, ls="", zorder=5)
    axd.set_xscale("log")
    axd.axhline(1, color="#595959", lw=0.8, ls=(0, (3, 2)))
    axd.axvline(1, color="#7F7F7F", lw=0.8, ls=(0, (1, 1.5)))
    axd.set_xlim(0.12, 3.5); axd.set_ylim(0.8, 1.38)
    from matplotlib.ticker import NullFormatter, NullLocator
    axd.xaxis.set_minor_locator(NullLocator()); axd.xaxis.set_minor_formatter(NullFormatter())
    axd.set_xticks([0.15, 0.3, 1, 2]); axd.set_xticklabels(["0.15", "0.3", "1", "2"])
    axd.set_xlabel("Genome-wide loading SNR in training animals")
    axd.set_ylabel("Held-out outside-sensor MSE / baseline")
    hd = [Line2D([], [], ls="", marker="o", ms=5, mfc="#7F7F7F", mec="#000000", label="score removal (mean)"),
          Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec="#000000", label="gain-bounded score (mean)")]
    hd += [Line2D([], [], ls="", marker=smk[a], ms=4, mfc=scol[a], mec=scol[a], label=n.replace("\n", " "))
           for a, n in settings]
    axd.legend(handles=hd, loc="upper right", frameon=False, fontsize=7, borderaxespad=0.2, ncol=1)
    axd.text(0.13, 0.82, "three settings; descriptive", fontsize=7, color="#595959")
    save(fig, "Figure_4")


# ----------------------------------------------------------------------------------- Figure 1
def fig1():
    acc = pd.read_csv(PKG / "results/dataset_accounting.csv").set_index("accession")
    n = lambda a: int(acc.loc[a, "n_animals_in_region"])  # noqa: E731
    fig = plt.figure(figsize=(180 * MM, 118 * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
    ax.text(0.012, 0.975, "a", fontsize=9, fontweight="bold", va="top", transform=ax.transAxes)
    ax.text(0.035, 0.975, "Seven core NASA brain deposits (three missions) and four external datasets",
            fontsize=7.5, fontweight="bold", va="top", transform=ax.transAxes)
    groups = [
        ("RRRM-2", "female C57BL/6NTac; 12 and 29 weeks;\nISS-terminal and live-return cohorts",
         [f"OSD-561 cerebellum, {n('OSD-561')}", f"OSD-562 hippocampus, {n('OSD-562')}",
          f"OSD-613 left hemisphere, {n('OSD-613')} (ISS-T only)"], 0.012, 0.30),
        ("RR-10", "female B6129SF2/J;\ncollection matched within deposit",
         [f"OSD-563 cerebellum, {n('OSD-563')}", f"OSD-564 hippocampus, {n('OSD-564')}",
          f"OSD-612 left hemisphere, {n('OSD-612')}"], 0.322, 0.22),
        ("RR-23", "male C57BL/6J; flight, habitat\nand vivarium ground controls",
         [f"OSD-525 cerebellum, {n('OSD-525')}"], 0.552, 0.19),
        ("External (refitted only)", "Earth irradiation/unloading and RR-3",
         ["OSD-202 brain, 40 mice", "OSD-352 RR-3 brain, 6 mice", "OSD-374 hippocampus, 8 pools",
          "OSD-611 hippocampus, 10 pools"], 0.752, 0.236),
    ]
    top, hgt = 0.905, 0.30
    for name, sub, chips, x, w in groups:
        rbox(ax, x, top - hgt, w, hgt, "", "#F2F2F2", ec="#7F7F7F")
        ax.text(x + 0.008, top - 0.012, name, fontsize=7.5, fontweight="bold", va="top", transform=ax.transAxes)
        ax.text(x + 0.008, top - 0.055, sub, fontsize=7, va="top", transform=ax.transAxes, linespacing=1.2)
        for k, c in enumerate(chips):
            rbox(ax, x + 0.008, top - 0.145 - k * 0.046, w - 0.016, 0.038, c, "#FFFFFF", ec="#BFBFBF", fs=7, lw=0.6)
    ax.text(0.5, top - hgt - 0.03, "129 libraries  \u2192  126 animal\u2013region profiles  \u2192  83 distinct mice "
            "(numbers are animals per deposit; regions of one animal are never split)",
            ha="center", fontsize=7, va="top", transform=ax.transAxes)

    ax.text(0.012, 0.49, "b", fontsize=9, fontweight="bold", va="top", transform=ax.transAxes)
    ax.text(0.035, 0.49, "Four questions, designs and answers", fontsize=7.5, fontweight="bold", va="top",
            transform=ax.transAxes)
    cols = [
        ("Collection protocol or\nspaceflight?",
         "2\u00d72 flight \u00d7 protocol\ndecomposition; exact whole-\nanimal permutation; 17\ncontrasts in 3 missions",
         "Protocol effect 0.70\u20131.12 log$_2$;\nmatched flight effects\n\u22120.40 to +0.28 (Fig. 2)"),
        ("What biology survives\nprotocol control?",
         "Count models; flight \u00d7\ncollection interactions;\n54,989-test cross-cohort\nscreen",
         "Circadian and erythroid\ncandidates; no gene at FDR\n< 0.05 in every contrast\n(Fig. 3)"),
        ("Can a detectable signature\nbe removed post hoc?",
         "362 flight-naive holdouts;\nidentical-loading noise\nceiling; operator gain",
         "Detected (AUC 0.96\u20131.00);\nloading imprecise; removal\ngain 8\u201315 (Fig. 4)"),
        ("What can a bounded\ncorrection achieve?",
         "Gain-constrained convex\nregression (L = 2); 4 external\ndatasets; 1,400 count\nsimulations",
         "1.2\u20137.2% lower held-out\nerror at gain \u2264 2; no\nbiological recovery\nguarantee (Figs. 5\u20137)"),
    ]
    cw, gap, x0 = 0.226, 0.018, 0.017
    rows = [(0.425, 0.075, "#FBE6C2"), (0.300, 0.105, "#D9EEFB"), (0.150, 0.105, "#CCEBE1")]
    for i, texts in enumerate(cols):
        x = x0 + i * (cw + gap)
        for (y, hh, fc), t in zip(rows, texts):
            rbox(ax, x, y - hh, cw, hh, t, fc, fs=7, weight="bold" if fc == "#FBE6C2" else "normal")
        arrow(ax, (x + cw / 2, rows[0][0] - rows[0][1]), (x + cw / 2, rows[1][0]), "#404040", lw=0.9)
        arrow(ax, (x + cw / 2, rows[1][0] - rows[1][1]), (x + cw / 2, rows[2][0]), "#404040", lw=0.9)
    keyy = 0.012
    for k, (fc, lab) in enumerate((("#FBE6C2", "question"), ("#D9EEFB", "design and test"), ("#CCEBE1", "answer"))):
        rbox(ax, 0.30 + k * 0.15, keyy, 0.025, 0.022, "", fc)
        ax.text(0.33 + k * 0.15, keyy + 0.011, lab, fontsize=7, va="center", transform=ax.transAxes)
    save(fig, "Figure_1")


# ----------------------------------------------------------------------------------- Figure 3
def fig3():
    from matplotlib.lines import Line2D
    cmp_ = pd.read_csv(PKG / "results/count_model_panel_comparison.csv")
    fci = pd.read_csv(PKG / "extensions/flight_collection_interaction/panel_flight_collection_interactions.csv")
    blood = pd.read_csv(PKG / "extensions/flight_collection_interaction/fixed_blood_marker_animal_scores.csv")
    rep = pd.read_csv(PKG / "results/regional_gene_replication_summary.csv")
    panels = ["rPRG", "dPRG", "SRG", "HSP", "Circadian", "Glia_Complement", "BBB_Endothelium", "Sensor31"]
    pc = dict(zip(panels, [OI["vermillion"], OI["orange"], OI["blue"], OI["green"], OI["purple"],
                           OI["sky"], OI["grey"], OI["black"]]))
    pm = dict(zip(panels, ["o", "s", "^", "D", "v", "P", "X", "*"]))
    fig = plt.figure(figsize=(180 * MM, 132 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.9, 1.1], left=0.075, right=0.985, top=0.96, bottom=0.075,
                          wspace=0.34, hspace=0.42)

    # (a) count model versus linear panel effects
    ax = fig.add_subplot(gs[0, 0])
    ax.text(-0.2, 1.0, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for p in panels:
        s = cmp_[cmp_.panel == p]
        ax.scatter(s.mean_log2_effect, s.mean_NB_log2_effect, s=12, marker=pm[p], facecolor=pc[p],
                   edgecolor=pc[p], lw=0.5, alpha=0.85, label=p.replace("_", " "))
    lim = (-1.4, 1.6)
    ax.plot(lim, lim, color="#595959", lw=0.7, ls=(0, (3, 2)))
    ax.axhline(0, color="#BFBFBF", lw=0.6); ax.axvline(0, color="#BFBFBF", lw=0.6)
    ax.set_xlim(lim); ax.set_ylim(lim)
    r = np.corrcoef(cmp_.mean_log2_effect, cmp_.mean_NB_log2_effect)[0, 1]
    agree = int(cmp_.same_effect_direction.sum())
    ax.text(0.03, 0.97, f"{agree}/{len(cmp_)} same direction\nPearson r = {r:.3f}", transform=ax.transAxes,
            va="top", fontsize=7)
    ax.set_xlabel("Linear panel effect (log$_2$)"); ax.set_ylabel("Negative-binomial panel effect (log$_2$)")
    ax.legend(loc="lower right", frameon=False, fontsize=6.6, ncol=2, handletextpad=0.1, columnspacing=0.4,
              borderaxespad=0.1, markerscale=0.9)

    # (b) flight x collection interaction per panel
    axb = fig.add_subplot(gs[0, 1])
    axb.text(-0.36, 1.0, "b", transform=axb.transAxes, fontsize=9, fontweight="bold", va="bottom")
    order = ["rPRG", "dPRG", "SRG", "HSP", "Circadian", "Glia_Complement", "BBB_Endothelium", "Sensor31"]
    for i, p in enumerate(order):
        for region, dy, mk in (("cerebellum", 0.17, "o"), ("hippocampus", -0.17, "^")):
            for sens, filled, ddy in (("all_animals", True, 0.0), ("exclude_any_region_metadata_conflict", False, -0.08)):
                r_ = fci[(fci.region == region) & (fci.panel == p) & (fci.sensitivity == sens)]
                if len(r_) != 1:
                    continue
                r_ = r_.iloc[0]
                y = i + dy + ddy
                col = OI["vermillion"] if region == "cerebellum" else OI["blue"]
                axb.plot([r_.bootstrap_lower95, r_.bootstrap_upper95], [y, y], color=col, lw=0.8 if filled else 0.6,
                         ls="-" if filled else (0, (1.5, 1)))
                axb.plot(r_.interaction_log2, y, mk, ms=4 if filled else 3.4, mfc=col if filled else "white",
                         mec=col, mew=0.8)
                if filled and r_.q_BH_complete_panel_family <= 0.05:
                    axb.text(r_.bootstrap_upper95 + 0.05, y, f"q = {r_.q_BH_complete_panel_family:.3f}",
                             fontsize=7, va="center")
    axb.axvline(0, color="#595959", lw=0.6, ls=(0, (2, 2)))
    axb.set_yticks(range(len(order))); axb.set_yticklabels([o.replace("_", " ") for o in order])
    axb.set_ylim(len(order) - 0.5, -0.6)
    axb.set_xlabel("Flight \u00d7 collection interaction, (FLT\u2212GC)$_{ISS-T}$ \u2212 (FLT\u2212GC)$_{LAR}$ (log$_2$)")
    hb = [Line2D([], [], ls="-", lw=0.8, marker="o", ms=4, color=OI["vermillion"], label="cerebellum"),
          Line2D([], [], ls="-", lw=0.8, marker="^", ms=4, color=OI["blue"], label="hippocampus"),
          Line2D([], [], ls=(0, (1.5, 1)), lw=0.6, marker="o", ms=3.4, mfc="white", mec="#404040", color="#404040",
                 label="conflict animals\nexcluded")]
    axb.legend(handles=hb, loc="upper right", frameon=False, fontsize=6.8, borderaxespad=0.1)
    axb.set_xlim(-1.25, 2.05)
    axb.set_xticks([-1, -0.5, 0, 0.5, 1])

    # (c) erythroid indicator by group
    axc = fig.add_subplot(gs[1, 0])
    axc.text(-0.2, 1.0, "c", transform=axc.transAxes, fontsize=9, fontweight="bold", va="bottom")
    groups = [("FLT", "ISS_T", "FLT\nISS-T"), ("GC", "ISS_T", "GC\nISS-T"), ("FLT", "LAR", "FLT\nLAR"),
              ("GC", "LAR", "GC\nLAR")]
    rng = np.random.default_rng(3)
    for j, (acc_, col, mk) in enumerate((("OSD-561", OI["vermillion"], "o"), ("OSD-562", OI["blue"], "^"))):
        for i, (fl, co, lab) in enumerate(groups):
            s = blood[(blood.accession == acc_) & (blood.flight == fl) & (blood.cohort == co)]
            x = i + (-0.17 if j == 0 else 0.17)
            axc.scatter(x + rng.uniform(-0.05, 0.05, len(s)), s.blood_indicator_mean_log2, s=14, marker=mk,
                        facecolor="white" if s.age.astype(str).iloc[0] == "x" else col, edgecolor=col, lw=0.7,
                        alpha=0.9)
            axc.plot([x - 0.1, x + 0.1], [s.blood_indicator_mean_log2.mean()] * 2, color="#000000", lw=1.2)
    axc.set_xticks(range(4)); axc.set_xticklabels([g[2] for g in groups])
    axc.set_ylabel("Erythroid marker mean (log$_2$, 7 genes)")
    hc = [Line2D([], [], ls="", marker="o", ms=4, color=OI["vermillion"], label="cerebellum"),
          Line2D([], [], ls="", marker="^", ms=4, color=OI["blue"], label="hippocampus"),
          Line2D([], [], ls="-", lw=1.2, color="#000000", label="group mean")]
    axc.legend(handles=hc, loc="upper left", frameon=False, fontsize=7, borderaxespad=0.1)

    # (d) cross-cohort gene screen
    axd = fig.add_subplot(gs[1, 1])
    axd.text(-0.16, 1.0, "d", transform=axd.transAxes, fontsize=9, fontweight="bold", va="bottom")
    names = {"cerebellum": "Cerebellum\n(3 contrasts)", "hippocampus": "Hippocampus\n(3 contrasts)",
             "left_cerebral_hemisphere": "Left hemisphere\n(2 contrasts)"}
    bars = [("n_genes", "genes tested", "#E6E6E6", ""), ("direction_concordant", "direction concordant", OI["sky"], ""),
            ("nominal_max_p05", "nominal max-p < 0.05", OI["orange"], "//"), ("BH_FDR05", "BH FDR < 0.05", OI["vermillion"], "xx")]
    wbar = 0.2
    for i, (_, r_) in enumerate(rep.iterrows()):
        for k, (col, lab, fc, hatch) in enumerate(bars):
            v = r_[col]
            x = i + (k - 1.5) * wbar
            if v > 0:
                axd.bar(x, v, width=wbar * 0.92, color=fc, edgecolor="#404040", lw=0.6, hatch=hatch,
                        label=lab if i == 0 else None, log=True)
                axd.text(x, v * 1.25, f"{int(v):,}", ha="center", va="bottom", fontsize=6.5, rotation=90)
            else:
                axd.bar(x, 0.0, width=wbar * 0.92, color=fc, edgecolor="#404040", lw=0.6, hatch=hatch,
                        label=lab if i == 0 else None)
                axd.text(x, 0.75, "0", ha="center", va="bottom", fontsize=7, fontweight="bold")
    axd.set_yscale("log"); axd.set_ylim(0.5, 3e7)
    axd.set_xticks(range(len(rep))); axd.set_xticklabels([names[t] for t in rep.tissue])
    axd.set_ylabel("Genes (log scale)")
    axd.legend(loc="upper center", frameon=False, fontsize=6.8, ncol=2, borderaxespad=0.1, columnspacing=0.6)
    save(fig, "Figure_3")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["fig2"]:
        globals()[name]()
