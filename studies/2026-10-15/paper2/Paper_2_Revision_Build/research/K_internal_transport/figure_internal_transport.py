#!/usr/bin/env python3
"""Figure 11: internal (held-out RR23) transport benchmark against external transport."""
import json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
D = np.load(HERE / "internal_transport_draws.npz"); S = json.loads((HERE / "summary.json").read_text())
BLUE, GREY, VERM, BLACK, SKY, GREEN = "#0072B2", "#9A9A9A", "#D55E00", "#000000", "#56B4E9", "#009E73"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
fig, axs = plt.subplots(1, 2, figsize=(183 / 25.4, 80 / 25.4), sharex=True)
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.17, top=0.80, wspace=0.18)
bins = np.linspace(-3, 9, 61)
for ax, (L, k, name) in zip(axs, [("a", "RR10", "RR10 panel (702 clusters)"), ("b", "RRRM2", "RRRM2 panel (281 clusters)")]):
    s = S[k]
    ax.hist(D[f"perm_{k}"], bins=bins, density=True, color=GREY, alpha=0.6, label="label-permuted reference")
    ax.hist(D[f"gcviv_{k}"], bins=bins, density=True, color=GREEN, alpha=0.75, label="same-date habitat–vivarium axis, held out")
    ax.hist(D[f"split_{k}"], bins=bins, density=True, color=BLUE, alpha=0.85, label="flight–habitat axis, held out (15,876 splits)")
    lo, med, hi = s["lopo_q025_median_q975"]
    top = ax.get_ylim()[1]; ax.set_ylim(0, top * 1.55); y0 = top * 1.38
    ax.errorbar([med], [y0], xerr=[[med - lo], [hi - med]], fmt="s", color=SKY, ms=4, capsize=2, lw=1, label="held-out pair, 8 v 8 training (2.5–97.5%)")
    e, (cl, ch) = s["external_effect_pp"], s["external_CI95"]
    ax.errorbar([e], [top * 1.18], xerr=[[e - cl], [ch - e]], fmt="o", color=VERM, ms=5, capsize=2, lw=1.2, label="external cohort estimate (95% CI)")
    ax.axvline(0, c=BLACK, lw=0.6, ls=":")
    r = 100 * e / s["split_half_q025_median_q975"][1]
    ax.text(ch + 0.25, top * 1.18, f"point estimate = {r:.1f}% of internal median".replace("-", "\u2212"), ha="left", va="center", color=VERM)
    ax.set_title(f"{L})  {name}", loc="left", fontweight="bold", pad=4)
    ax.set_xlabel("Held-out transported score (pp)")
axs[0].set_ylabel("Density")
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.53, 1.0), handletextpad=0.3, columnspacing=1.0)
for ext in ["pdf", "png"]: fig.savefig(HERE / f"Figure_11_internal_transport.{ext}", dpi=600)
texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible()]
assert min(t.get_fontsize() for t in texts) >= 7
print({k: round(100 * S[k]["external_effect_pp"] / S[k]["split_half_q025_median_q975"][1], 1) for k in ["RR10", "RRRM2"]})
