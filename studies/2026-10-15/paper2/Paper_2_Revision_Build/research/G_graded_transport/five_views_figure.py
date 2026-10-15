#!/usr/bin/env python3
"""Figure 10: five views of the same per-cluster data (exploratory; see PLAN_2026-10-15.md).

a) scatter, b) hexbin, c) density contours, d) marginal distributions, e) 2D density heatmap.
183 mm wide, Arial >= 7 pt, Okabe-Ito colours. Every number printed on the figure is also written to
five_views_numbers.json so the manuscript can quote it exactly.
"""
import json, pathlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde, spearmanr

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parents[2] / "Paper_2_Complete_Submission_Package"
loc = pd.read_csv(PKG / "extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz")
per = pd.read_csv(HERE / "per_cluster_external_effects.csv", index_col=0)
res = pd.read_csv(HERE / "graded_transport_results.csv")
BLUE, ORANGE, GREEN, SKY, VERM, GREY, BLACK = "#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#8C8C8C", "#000000"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})

w = {c: loc[loc.contrast == c].set_index("cluster_id") for c in ["FLT_vs_GC", "FLT_vs_VIV", "GC_vs_VIV"]}
ids = w["FLT_vs_GC"].index
ev_gc = -np.log10(w["FLT_vs_GC"].exact_p); ev_cc = -np.log10(w["GC_vs_VIV"].loc[ids].exact_p)
tv_gc = 100 * w["FLT_vs_GC"].native_mean_total_variation; tv_viv = 100 * w["FLT_vs_VIV"].loc[ids].native_mean_total_variation
depth = w["FLT_vs_GC"].min_selected_reads
by_gc = w["FLT_vs_GC"].BY_q_full_3E <= 0.05
N = {}

def panel(ax, letter, title): ax.set_title(f"{letter})  {title}", loc="left", fontweight="bold", pad=4)

fig = plt.figure(figsize=(183 / 25.4, 150 / 25.4))
gs = GridSpec(2, 3, figure=fig, left=0.07, right=0.925, bottom=0.075, top=0.95, wspace=0.62, hspace=0.5)
DJ = json.loads((HERE / "disjoint_animal_checks.json").read_text())

# a) scatter: every cluster, RR23 evidence vs RR10 external effect
ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", "Scatter: every cluster")
e10 = 100 * per["RR10_FLT_minus_GC"].reindex(ids); m = e10.notna()
ax.scatter(ev_gc[m & ~by_gc], e10[m & ~by_gc], s=4, c=GREY, alpha=0.55, lw=0, label="other clusters")
ax.scatter(ev_gc[m & by_gc], e10[m & by_gc], s=9, c=VERM, lw=0, label="RR23 BY-positive")
ax.axhline(0, c=BLACK, lw=0.6, ls=":")
r = res[(res.cohort == "RR10") & (res.direction == "FLT_minus_GC")].iloc[0]
ax.text(0.03, 0.97, f"n = {int(r.n_clusters)} clusters\nSpearman ρ = {r.spearman_rho:.3f}\nexact p = {r.exact_p_one_sided:.3f}", transform=ax.transAxes, va="top")
ax.set_xlabel("RR23 flight–ground evidence, −log$_{10}$ p"); ax.set_ylabel("RR10 projected effect (pp)")
ax.legend(frameon=False, loc="lower right", handletextpad=0.2, borderaxespad=0.1)
N["a"] = dict(n=int(r.n_clusters), rho=float(r.spearman_rho), p=float(r.exact_p_one_sided), n_BY_in_RR10=int((m & by_gc).sum()))

# b) hexbin: where clusters accumulate in effect size x read depth
ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", "Hexbin: accumulation")
hb = ax.hexbin(np.log10(depth), tv_gc, gridsize=22, cmap="Blues", mincnt=1, linewidths=0.2)
ax.scatter(np.log10(depth[by_gc]), tv_gc[by_gc], s=8, facecolors="none", edgecolors=VERM, lw=0.7, label="BY-positive")
cb = fig.colorbar(hb, ax=ax, pad=0.02, fraction=0.05); cb.set_label("clusters per bin"); cb.ax.tick_params(labelsize=7)
ax.set_xticks([1.5, 2, 2.5, 3], ["32", "100", "316", "1,000"])
ax.set_xlabel("Minimum selected reads per animal"); ax.set_ylabel("Flight–ground effect, TV (pp)")
rho_b = spearmanr(depth, tv_gc).statistic
ax.text(0.03, 0.97, f"Spearman ρ = {rho_b:.2f}", transform=ax.transAxes, ha="left", va="top")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, 0.93), handletextpad=0.2, borderaxespad=0.1)
N["b"] = dict(rho_depth_tv=float(rho_b), median_depth_BY=float(depth[by_gc].median()), median_depth_other=float(depth[~by_gc].median()),
              median_tv_BY=float(tv_gc[by_gc].median()), median_tv_other=float(tv_gc[~by_gc].median()))

# c) density contours: comparator dependence of the same cluster effects
ax = fig.add_subplot(gs[0, 2]); panel(ax, "c", "Density contours")
X, Y = tv_gc.to_numpy(), tv_viv.to_numpy(); k = gaussian_kde(np.vstack([X, Y]))
gx, gy = np.mgrid[0:max(X.max(), Y.max()) * 1.02:120j, 0:max(X.max(), Y.max()) * 1.02:120j]
zz = k(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
ax.scatter(X, Y, s=2, c=GREY, alpha=0.35, lw=0)
cs = ax.contour(gx, gy, zz, levels=6, colors=[BLUE], linewidths=0.8)
lim = max(X.max(), Y.max()) * 1.02; ax.plot([0, lim], [0, lim], c=BLACK, lw=0.6, ls="--"); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
rho_c = spearmanr(X, Y).statistic; above = float((Y > X).mean())
ax.text(0.03, 0.97, f"shared FLT: ρ = {rho_c:.2f}\ndisjoint, depth-adjusted:\nρ = {DJ['c_disjoint_partial_median']:.2f}", transform=ax.transAxes, va="top")
ax.set_xlabel("Flight–habitat effect, TV (pp)"); ax.set_ylabel("Flight–vivarium effect, TV (pp)")
N["c"] = dict(rho=float(rho_c), frac_above_diag=above, n=len(X))
N["disjoint"] = DJ

# d) marginal distributions: the two external cohorts against each other
sub = gs[1, 0:2].subgridspec(4, 7, wspace=0.05, hspace=0.05)
axm = fig.add_subplot(sub[1:, 0:4]); axt = fig.add_subplot(sub[0, 0:4], sharex=axm); axr = fig.add_subplot(sub[1:, 4], sharey=axm)
a10, a2 = 100 * per["RR10_FLT_minus_GC"], 100 * per["RRRM2_FLT_minus_GC"]; both = a10.notna() & a2.notna()
axm.scatter(a10[both], a2[both], s=5, c=GREEN, alpha=0.7, lw=0)
axm.axhline(0, c=BLACK, lw=0.6, ls=":"); axm.axvline(0, c=BLACK, lw=0.6, ls=":")
bins = np.linspace(min(a10[both].min(), a2[both].min()), max(a10[both].max(), a2[both].max()), 30)
axt.hist(a10[both], bins=bins, color=SKY, edgecolor="white", lw=0.3); axr.hist(a2[both], bins=bins, color=ORANGE, edgecolor="white", lw=0.3, orientation="horizontal")
for a_ in (axt, axr): a_.axis("off")
axt.set_title("d)  Marginal distributions: two external cohorts", loc="left", fontweight="bold", pad=4)
rho_d = spearmanr(a10[both], a2[both]).statistic
axm.text(0.03, 0.97, f"n = {int(both.sum())} common clusters\nSpearman ρ = {rho_d:.3f} (descriptive)", transform=axm.transAxes, va="top")
axm.set_xlabel("RR10 projected flight–ground effect (pp)"); axm.set_ylabel("RRRM2 projected effect (pp)")
N["d"] = dict(n=int(both.sum()), rho=float(rho_d), RR10_median=float(a10[both].median()), RRRM2_median=float(a2[both].median()))

# e) 2D density heatmap: flight evidence vs control-control evidence
ax = fig.add_subplot(gs[1, 2]); panel(ax, "e", "2D density heatmap")
tv_cc = 100 * w["GC_vs_VIV"].loc[ids].native_mean_total_variation
h = ax.hist2d(tv_gc, tv_cc, bins=26, range=[[0, 30], [0, 30]], cmap="YlOrBr", cmin=1)
cb = fig.colorbar(h[3], ax=ax, pad=0.02, fraction=0.05); cb.set_label("clusters per bin"); cb.ax.tick_params(labelsize=7)
rho_e = spearmanr(tv_gc, tv_cc).statistic
ax.text(0.97, 0.97, f"shared GC: ρ = {rho_e:.2f}\ndisjoint, depth-adjusted:\nρ = {DJ['e_disjoint_partial_median']:.2f}\nnoise-only reference:\nρ = {DJ['noise_only_GCsplit_vs_FLTVIV_partial_median']:.2f}", transform=ax.transAxes, ha="right", va="top")
ax.set_xlabel("Flight–habitat effect, TV (pp)"); ax.set_ylabel("Habitat–vivarium effect, TV (pp)")
N["e"] = dict(rho_shared_TV=float(rho_e), n=len(tv_gc), n_outside_30pp=int(((tv_gc > 30) | (tv_cc > 30)).sum()))

for ext in ["pdf", "png"]: fig.savefig(HERE / f"Figure_10_five_views.{ext}", dpi=600)
texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible()]
N["min_font_pt"] = min(t.get_fontsize() for t in texts)
assert N["min_font_pt"] >= 7, N["min_font_pt"]
(HERE / "five_views_numbers.json").write_text(json.dumps(N, indent=1)); print(json.dumps(N, indent=1))
