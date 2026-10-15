#!/usr/bin/env python3
"""Exploratory amended stress benchmark; no target-derived tuning.

Fixed before execution: true gain 8/12, train N24/96, target N24/48,
unchanged/flipped diffuse coupling, 100 independent replicates each.
"""
import itertools
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from statistics_core import Contrast, bh, bounded_gamma, learn_anchor, transform


def main():
    seed, genes, repetitions = 86761014, 1200, 100
    setup = np.random.default_rng(seed)
    sensor = np.arange(31)
    outside = np.arange(31, genes)
    diffuse = np.arange(31,1151)
    coupling_sign = setup.choice([-1.,1.], len(diffuse))
    b = np.zeros(genes)
    b[200:260] = setup.choice([-1.,1.],60)*.8
    positive = abs(b)>1e-12
    sd = setup.uniform(.65,1.15,genes)
    factors = setup.normal(0,.12,(genes,5))
    methods = ["raw","projection","score","bounded_score"]
    configs = list(itertools.product([8.,12.],[24,96],[24,48],["unchanged","flipped"]))
    config = dict(label="Exploratory amended high-gain Gaussian stress benchmark",seed=seed,
                  repetitions_per_configuration=repetitions,genes=genes,sensor_genes=31,
                  diffuse_coupling_genes=len(diffuse),true_score_spectral_gain=[8.,12.],
                  training_samples=[24,96],target_samples=[24,48],coupling=["unchanged","flipped"],
                  ground_nuisance_group_difference=1.5,target_nuisance_group_difference=1.,
                  latent_nuisance_sd=.4,biological_genes=60,absolute_biological_effect=.8,
                  biological_sensor_score=0.,gene_noise_sd_range=[.65,1.15],latent_noise_factors=5,
                  latent_noise_gene_loading_sd=.12,methods=methods,bounded_max_gain=2.,
                  shared_training_target_samples=0,target_tuning=False,
                  flipped_definition="Every outside-sensor nuisance loading reverses sign; sensor loadings stay fixed.",
                  hypothesis_interpretation="FDP compares conditional gene-wise OLS/BH calls to original biological effects; it is a recovery diagnostic rather than formal biological FDR control.",
                  amendment_reason="Root Gaussian benchmark typical learned gain ≈2.15 is far below real score gain ≈8–15; this amendment tests higher diffuse coupling.",
                  limitation="Deliberate Gaussian stress models, not external biological validation or a complete RNA-seq simulation.")
    (HERE/"high_gain_benchmark_config.json").write_text(json.dumps(config,indent=2)+"\n")
    rows=[]
    w = np.zeros(genes)
    w[sensor]=1/np.sqrt(len(sensor))
    for cindex,(true_gain,ntrain,ntarget,coupling) in enumerate(configs):
        u = np.zeros(genes)
        u[sensor]=1.
        u[diffuse]=np.sqrt(len(sensor))*np.sqrt(true_gain**2-1)/np.sqrt(len(diffuse))*coupling_sign
        true_gamma=u/(w@u)
        assert abs(np.linalg.norm(true_gamma)-true_gain)<1e-10
        ut=u.copy()
        if coupling=="flipped":
            ut[outside]*=-1
        gtrain=np.r_[np.zeros(ntrain//2),np.ones(ntrain//2)]
        gtarget=np.r_[np.zeros(ntarget//2),np.ones(ntarget//2)]
        for repeat in range(repetitions):
            rng=np.random.default_rng(np.random.SeedSequence([seed,cindex,repeat]))
            train=10+u[:,None]*(1.5*gtrain+rng.normal(0,.4,ntrain))[None,:]
            train += sd[:,None]*rng.normal(size=(genes,ntrain))+factors@rng.normal(size=(5,ntrain))
            target=10+ut[:,None]*(gtarget+rng.normal(0,.4,ntarget))[None,:]+b[:,None]*gtarget[None,:]
            target += sd[:,None]*rng.normal(size=(genes,ntarget))+factors@rng.normal(size=(5,ntarget))
            anchor=learn_anchor(train,gtrain,None,sensor)
            original,_=Contrast(target).fit(gtarget)
            for method in methods:
                estimate_model=Contrast(transform(target,anchor,method))
                estimate,t=estimate_model.fit(gtarget)
                calls=bh(2*stats.t.sf(abs(t),estimate_model.df))<.05
                ncall=int(calls.sum())
                bio=transform(b[:,None],anchor,method)[:,0]
                gamma=bounded_gamma(anchor[1],anchor[2],2.) if method=="bounded_score" else anchor[2]
                gain=1. if method in ("raw","projection") else float(np.linalg.norm(gamma))
                rows.append(dict(true_gain=true_gain,n_train=ntrain,n_target=ntarget,coupling=coupling,
                                 repetition=repeat,method=method,
                                 learned_unbounded_gain=float(np.linalg.norm(anchor[2])),operator_gain=gain,
                                 biological_effect_RMSE=float(np.sqrt(np.mean((estimate-b)**2))),
                                 outside_sensor_effect_RMSE=float(np.sqrt(np.mean((estimate[outside]-b[outside])**2))),
                                 biological_FDP=float((calls&~positive).sum()/max(ncall,1)),
                                 power=float((calls&positive).sum()/positive.sum()),discoveries=ncall,
                                 biology_operator_relative_distortion=float(np.linalg.norm(bio-b)/np.linalg.norm(b)),
                                 biology_energy_retained=float(np.dot(bio,bio)/np.dot(b,b)),
                                 sensor_effect_energy_removed_pct=float(100*(1-np.dot(estimate[sensor],estimate[sensor])/np.dot(original[sensor],original[sensor]))),
                                 outside_sensor_effect_energy_removed_pct=float(100*(1-np.dot(estimate[outside],estimate[outside])/np.dot(original[outside],original[outside])))))
        print(f"High gain {true_gain:g}; training {ntrain}; target {ntarget}; {coupling}: {repetitions} completed",flush=True)
    frame=pd.DataFrame(rows)
    frame.to_csv(HERE/"high_gain_benchmark_repetitions.csv",index=False)
    groups=["true_gain","n_train","n_target","coupling","method"]
    metrics=[c for c in frame.columns if c not in groups+["repetition"]]
    summaries=[]
    for keys,part in frame.groupby(groups,sort=False):
        row=dict(zip(groups,keys));row["repetitions"]=len(part)
        for metric in metrics:
            values=part[metric].to_numpy()
            row[metric+"_mean"]=float(values.mean())
            row[metric+"_mc_se"]=float(values.std(ddof=1)/np.sqrt(len(values)))
            row[metric+"_lower95"]=float(np.quantile(values,.025))
            row[metric+"_upper95"]=float(np.quantile(values,.975))
        summaries.append(row)
    summary=pd.DataFrame(summaries)
    summary.to_csv(HERE/"high_gain_benchmark_summary.csv",index=False)
    # Four panels retain both coupling directions and all information regimes.
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"svg.fonttype":"none"})
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.3),layout="constrained")
    colors={"raw":"#60666D","projection":"#326C9A","score":"#C75B39","bounded_score":"#277C61"}
    labels={"raw":"Raw","projection":"Projection","score":"Score","bounded_score":"Bounded score (L=2)"}
    for ax,(coupling,metric,title,letter) in zip(axes.flat,[
        ("unchanged","biological_effect_RMSE","Unchanged coupling: effect error","A"),
        ("flipped","biological_effect_RMSE","Flipped coupling: effect error","B"),
        ("unchanged","biological_FDP","Unchanged coupling: biological FDP","C"),
        ("flipped","biological_FDP","Flipped coupling: biological FDP","D")]):
        configs2=list(itertools.product([8.,12.],[24,96],[24,48]))
        x=np.arange(len(configs2))
        for im,method in enumerate(methods):
            part=summary[(summary.coupling==coupling)&(summary.method==method)].set_index(["true_gain","n_train","n_target"]).loc[configs2]
            ax.errorbar(x+(im-1.5)*.13,part[metric+"_mean"],yerr=1.96*part[metric+"_mc_se"],fmt="o",ms=3.6,capsize=2,color=colors[method],label=labels[method],lw=.9)
        ax.set_xticks(x,[f"{g:g}\n{nt}/{ne}" for g,nt,ne in configs2]);ax.set_xlabel("True gain / training N / target N")
        ax.set_ylabel("Mean effect RMSE (log₂ units)" if "RMSE" in metric else "Mean biological FDP")
        ax.set_title(f"{letter}   {title}",loc="left",fontweight="bold")
        ax.grid(axis="y",alpha=.18);ax.spines[["top","right"]].set_visible(False)
        if metric=="biological_FDP":ax.set_ylim(-.015,1.03)
    handles,legend=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,legend,loc="outside lower center",ncol=4,frameon=False,bbox_to_anchor=(.5,-.025))
    fig.suptitle("High-gain stress tests with known sensor-orthogonal biology",fontweight="bold",fontsize=12)
    fig.savefig(HERE/"high_gain_benchmark_panel.png",dpi=600,bbox_inches="tight")
    fig.savefig(HERE/"high_gain_benchmark_panel.svg",bbox_inches="tight")
    plt.close(fig)
    checks=dict(n_configurations=len(configs),n_repetitions=len(frame)//len(methods),n_rows=len(frame),
                all_64_configuration_method_cells_have_100_repetitions=bool((summary.repetitions==100).all() and len(summary)==64),
                maximum_bounded_gain=float(frame.loc[frame.method=="bounded_score","operator_gain"].max()),
                maximum_score_biology_relative_distortion=float(frame.loc[frame.method.isin(["score","bounded_score"]),"biology_operator_relative_distortion"].max()),
                minimum_learned_gain=float(frame.learned_unbounded_gain.min()),maximum_learned_gain=float(frame.learned_unbounded_gain.max()),
                independent_training_target=True,target_tuning=False)
    (HERE/"high_gain_benchmark_checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    print(json.dumps(checks,indent=2))


if __name__=="__main__":main()
