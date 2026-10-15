#!/usr/bin/env python3
"""Prespecified ground-only animal holdouts with training-only normalization."""
import itertools
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from run_revision import load_dataset, indices
from statistics_core import Contrast, bounded_gamma, learn_anchor, transform

KGRID=np.array([0.,.25,.5,.75,1.])


def train_reference(counts,train,test):
    """Filter and normalization reference depend exclusively on training counts."""
    canon=counts.index.astype(str).str.match(r"^ENSMUSG\d+$")
    chosen=canon&(counts.iloc[:,train].mean(axis=1).to_numpy()>=10)
    c=counts.loc[chosen].to_numpy(float)
    positive=(c[:,train]>0).all(axis=1)
    logref=np.log(c[positive][:,train]).mean(axis=1)
    def expression(cols):
        values=c[:,cols]
        factors=[]
        for col in cols:
            nz=c[positive,col]>0
            factors.append(np.exp(np.median(np.log(c[positive,col][nz])-logref[nz])))
        return np.log2(values/np.array(factors)+1),np.array(factors)
    yt,sft=expression(train);ye,sfe=expression(test)
    return counts.index[chosen],yt,ye,sft,sfe,int(positive.sum())


def age_design(meta):
    return np.column_stack([np.ones(len(meta)),(meta.age=="29").to_numpy(float)]) if meta.age.nunique()>1 else np.ones((len(meta),1))


def model(counts,meta,train,test,annotation):
    genes,yt,ye,sft,sfe,nref=train_reference(counts,train,test)
    mt=meta.iloc[train];me=meta.iloc[test]
    xt=age_design(mt)
    # Test gets exactly the training design columns, even if one cell is absent.
    xe=np.column_stack([np.ones(len(me)),(me.age=="29").to_numpy(float)]) if xt.shape[1]==2 else np.ones((len(me),1))
    nuisance=xt[:,1:] if xt.shape[1]>1 else None
    panel=indices(genes,annotation);sensor=panel["Sensor31"]
    outside=np.setdiff1d(np.arange(len(genes)),sensor)
    gt=mt.label.to_numpy(int);ge=me.label.to_numpy(int)
    anchor=learn_anchor(yt,gt,nuisance,sensor)
    B=np.linalg.lstsq(xt,yt.T,rcond=None)[0]
    residual=ye-(xe@B).T
    score=anchor[1]@residual
    return dict(genes=genes,yt=yt,ye=ye,xt=xt,xe=xe,anchor=anchor,sensor=sensor,outside=outside,
                score=score,residual=residual,gt=gt,ge=ge,sft=sft,sfe=sfe,nref=nref)


def tune_shrinkage(counts,meta,train,annotation):
    """Nested ground-training leave-one-animal-out prediction; no outer test use."""
    errors=np.zeros(len(KGRID));baseline=0.;features=0
    for left in train:
        inner=np.array([i for i in train if i!=left])
        d=model(counts,meta,inner,np.array([left]),annotation)
        outside=d["outside"];r=d["residual"][outside]
        predicted=np.outer(d["anchor"][2][outside],d["score"])
        for j,k in enumerate(KGRID):errors[j]+=float(np.sum((r-k*predicted)**2))
        baseline+=float(np.sum(r*r));features+=r.size
    selected=int(np.argmin(errors))
    return float(KGRID[selected]),float(errors[selected]/max(baseline,1e-20)),errors/max(baseline,1e-20)


def stratified_auc(score,meta):
    wins=[]
    for _,m in meta.assign(score=score).groupby("age",sort=True):
        p=m.loc[m.label==1,"score"].to_numpy();n=m.loc[m.label==0,"score"].to_numpy()
        wins.extend((a>b)+.5*(a==b) for a in p for b in n)
    return float(np.mean(wins))


def cosine(a,b):
    return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))


def splits(meta,n):
    if n in ("561","562"):
        cells=[np.asarray(v) for _,v in meta.reset_index(drop=True).groupby(["label","age"],sort=True).indices.items()]
        for held in itertools.product(*cells):
            test=np.array(held);train=np.array([i for i in range(len(meta)) if i not in held]);yield train,test
    else:
        positive=np.flatnonzero(meta.label.to_numpy()==1);negative=np.flatnonzero(meta.label.to_numpy()==0)
        hold=list(itertools.product(itertools.combinations(positive,3),itertools.combinations(negative,3)))
        rng=np.random.default_rng(52520261015)
        for index in np.sort(rng.choice(len(hold),200,replace=False)):
            test=np.array(hold[index][0]+hold[index][1]);train=np.setdiff1d(np.arange(len(meta)),test);yield train,test


def main():
    if '--figures-only' not in sys.argv:
        config=dict(label="Ground-only animal holdout; exploratory stability and prediction analysis",seed_OSD525=52520261015,
                    OSD561_OSD562="All 81 splits: hold out one of three animals from each age×collection cell; train8/test4",
                    OSD525="200 distinct splits sampled without replacement from all 4704 holdouts; train6HGC+5VGC/test3HGC+3VGC",
                    gene_filter="Canonical ENSMUSG and training mean count >=10",normalization="Training-positive genes geometric means; target size factor relative to that training reference using nonzero target entries",
                    no_flight_data=True,outer_test_animals_disjoint=True,training_target_reference_independent=True,
                    shrinkage_grid=KGRID.tolist(),shrinkage_tuning="Nested leave-one-ground-training-animal-out, redoing filtering, normalization, nuisance and sensor fitting each fold; minimize outside-sensor squared prediction error",
                    prediction_target="Held-out expression after subtracting training-fitted intercept and age; sensor score is an available held-out predictor",
                    prediction_baseline="Training-fitted intercept and age only",score_auc="Within-age positive-negative pair discrimination; one-sided orientation fixed by positive training contrast",
                    uncertainty="Split distribution quantiles are robustness ranges, not independent-sample confidence intervals; repeated splits overlap and OSD561/562 share many animals",
                    interpretation="Held-out association and loading generalization within ground collections; no biological validation or proof of true nuisance removal")
        (HERE/"protocol.json").write_text(json.dumps(config,indent=2)+"\n")
        annotation=pd.read_csv(ROOT/"inputs/gene_annotation.csv.gz")
        result=[];predictions=[];roles=[];cvrows=[];normalrows=[]
        for n in ("561","562","525"):
            d=load_dataset(n);meta=d["meta"]
            mask=meta.flight.eq("GC") if n!="525" else meta.flight.isin(["HGC","VGC"])
            meta=meta.loc[mask].copy();counts=d["counts"].loc[:,meta.index]
            meta["label"]=(meta.cohort.eq("ISS_T") if n!="525" else meta.flight.eq("HGC")).astype(int)
            allsplits=list(splits(meta,n))
            for sid,(train,test) in enumerate(allsplits):
                fit=model(counts,meta,train,test,annotation)
                shrink,inner_ratio,inner_errors=tune_shrinkage(counts,meta,train,annotation)
                v,w,gamma=fit["anchor"];outside=fit["outside"];sensor=fit["sensor"]
                ze=fit["xe"][:,1:] if fit["xe"].shape[1]>1 else None
                rawbeta,_=Contrast(fit["ye"],ze).fit(fit["ge"])
                trainbeta,_=Contrast(fit["yt"],fit["xt"][:,1:] if fit["xt"].shape[1]>1 else None).fit(fit["gt"])
                mse0=float(np.mean(fit["residual"][outside]**2))
                score=fit["score"];whole_score=v@fit["residual"]
                common=dict(accession=f"OSD-{n}",split_id=sid,n_training=len(train),n_test=len(test),n_genes=len(fit["genes"]),n_sensor=len(sensor),n_reference_genes=fit["nref"],
                            n_shared_animals=0,sensor_age_stratified_AUC=stratified_auc(score,meta.iloc[test]),whole_genome_age_stratified_AUC=stratified_auc(whole_score,meta.iloc[test]),
                            train_test_contrast_cosine=cosine(trainbeta,rawbeta),train_test_sensor_contrast_cosine=cosine(trainbeta[sensor],rawbeta[sensor]),
                            learned_score_gain=float(np.linalg.norm(gamma)),cv_selected_shrinkage=shrink,inner_cv_selected_error_ratio=inner_ratio,
                            outside_nuisance_only_prediction_MSE=mse0)
                for mode in ("raw","projection","score","bounded_score","training_cv_shrinkage"):
                    if mode=="training_cv_shrinkage":
                        corrected=fit["ye"]-shrink*np.outer(gamma,w@fit["ye"]);predgamma=shrink*gamma
                        gain=float(np.linalg.svd(np.array([[1-shrink,0],[-shrink*np.linalg.norm(gamma-w),1]]),compute_uv=False)[0])
                    elif mode=="projection":
                        corrected=transform(fit["ye"],fit["anchor"],mode)
                        predgamma=v;gain=1.
                    elif mode=="raw":corrected=fit["ye"];predgamma=np.zeros_like(gamma);gain=1.
                    else:
                        corrected=transform(fit["ye"],fit["anchor"],mode)
                        predgamma=bounded_gamma(w,gamma,2.) if mode=="bounded_score" else gamma
                        gain=float(np.linalg.norm(predgamma))
                    beta,_=Contrast(corrected,ze).fit(fit["ge"])
                    if mode=="projection":predicted=np.outer(v,whole_score)
                    else:predicted=np.outer(predgamma,score)
                    mse=float(np.mean((fit["residual"][outside]-predicted[outside])**2))
                    # Whole-genome projection uses the response genes in its score;
                    # its error is reconstruction, not disjoint-predictor prediction.
                    result.append(dict(**common,method=mode,operator_gain=gain,outside_rankone_reconstruction_MSE=mse,
                                       outside_prediction_MSE=np.nan if mode=="projection" else mse,
                                       outside_prediction_MSE_ratio_to_nuisance_only=np.nan if mode=="projection" else mse/mse0,
                                       heldout_sensor_effect_energy_removed_pct=float(100*(1-np.sum(beta[sensor]**2)/np.sum(rawbeta[sensor]**2))),
                                       heldout_outside_effect_energy_removed_pct=float(100*(1-np.sum(beta[outside]**2)/np.sum(rawbeta[outside]**2))),
                                       heldout_total_effect_energy_removed_pct=float(100*(1-np.dot(beta,beta)/np.dot(rawbeta,rawbeta)))))
                for j,i in enumerate(test):predictions.append(dict(accession=f"OSD-{n}",split_id=sid,animal=meta.iloc[i].animal,age=meta.iloc[i].age,label=int(meta.iloc[i].label),sensor_score=float(score[j]),whole_genome_score=float(whole_score[j])))
                for role,draw in [("train",train),("test",test)]:
                    for i in draw:roles.append(dict(accession=f"OSD-{n}",split_id=sid,animal=meta.iloc[i].animal,role=role,age=meta.iloc[i].age,label=int(meta.iloc[i].label)))
                for k,value in zip(KGRID,inner_errors):cvrows.append(dict(accession=f"OSD-{n}",split_id=sid,shrinkage=k,inner_error_ratio=float(value),selected=bool(k==shrink)))
                for role,draw,sf in [("train",train,fit["sft"]),("test",test,fit["sfe"])]:
                    for i,factor in zip(draw,sf):normalrows.append(dict(accession=f"OSD-{n}",split_id=sid,animal=meta.iloc[i].animal,role=role,size_factor=float(factor)))
                if sid%20==0:print(f"OSD-{n}: split {sid+1}/{len(allsplits)} complete",flush=True)
            print(f"OSD-{n}: all {len(allsplits)} holdouts complete",flush=True)
        frame=pd.DataFrame(result);frame.to_csv(HERE/"holdout_split_results.csv",index=False)
        pd.DataFrame(predictions).to_csv(HERE/"heldout_scores.csv",index=False)
        pd.DataFrame(roles).to_csv(HERE/"split_animal_roles.csv",index=False)
        pd.DataFrame(cvrows).to_csv(HERE/"nested_training_cv.csv",index=False)
        pd.DataFrame(normalrows).to_csv(HERE/"normalization_factors.csv",index=False)
        metrics=["sensor_age_stratified_AUC","whole_genome_age_stratified_AUC","train_test_contrast_cosine","train_test_sensor_contrast_cosine","learned_score_gain","cv_selected_shrinkage","outside_prediction_MSE_ratio_to_nuisance_only","heldout_sensor_effect_energy_removed_pct","heldout_outside_effect_energy_removed_pct","heldout_total_effect_energy_removed_pct"]
        summaries=[]
        for (accession,method),part in frame.groupby(["accession","method"],sort=False):
            row=dict(accession=accession,method=method,n_splits=len(part),uncertainty="overlapping-split robustness range")
            for metric in metrics:
                vals=part[metric].to_numpy();row[metric+"_mean"]=float(vals.mean());row[metric+"_median"]=float(np.median(vals));row[metric+"_q025"]=float(np.quantile(vals,.025));row[metric+"_q975"]=float(np.quantile(vals,.975))
            summaries.append(row)
        summary=pd.DataFrame(summaries);summary.to_csv(HERE/"holdout_summary.csv",index=False)
    else:
        config=json.loads((HERE/'protocol.json').read_text())
        frame=pd.read_csv(HERE/'holdout_split_results.csv')
        summary=pd.read_csv(HERE/'holdout_summary.csv')
    # Journal-safe font: explicitly load installed Arial; editable SVG text.
    font_path=Path("/System/Library/Fonts/Supplemental/Arial.ttf")
    family="Arial" if font_path.exists() else "DejaVu Sans"
    if font_path.exists():font_manager.fontManager.addfont(str(font_path))
    plt.rcParams.update({"font.family":family,"font.size":9.5,"svg.fonttype":"none"})
    fig,axes=plt.subplots(2,2,figsize=(7.5,6.5),layout="constrained")
    colors={"OSD-561":"#326C9A","OSD-562":"#277C61","OSD-525":"#A65F29"}
    raw=frame[frame.method=="raw"]
    for x,n in enumerate(["OSD-561","OSD-562","OSD-525"]):
        part=raw[raw.accession==n]
        for off,metric,label,marker in [(-.13,"sensor_age_stratified_AUC","Sensor","o"),(.13,"whole_genome_age_stratified_AUC","Whole genome","s")]:
            vals=part[metric].to_numpy();axes[0,0].plot(x+off,np.median(vals),marker=marker,color=colors[n],ms=5,label=label if x==0 else None)
            axes[0,0].vlines(x+off,np.quantile(vals,.025),np.quantile(vals,.975),color=colors[n],lw=1.3)
    axes[0,0].axhline(.5,ls="--",lw=.8,color="#777777");axes[0,0].set_ylim(-.03,1.03);axes[0,0].set_ylabel("Held-out age-stratified AUC");axes[0,0].set_xticks(range(3),["Cerebellum\nRRRM-2","Hippocampus\nRRRM-2","Cerebellum\nRR-23"]);axes[0,0].legend(frameon=False,fontsize=8.5,loc="lower left")
    methodnames=["raw","projection","score","bounded_score","training_cv_shrinkage"]
    short=["Baseline","Projection","Score","Bounded","CV shrink"]
    for ax,metric,ylabel in [(axes[0,1],"outside_prediction_MSE_ratio_to_nuisance_only","Outside-sensor prediction MSE / baseline"),(axes[1,1],"heldout_outside_effect_energy_removed_pct","Held-out outside-sensor energy removed (%)")]:
        shown=methodnames if "MSE" not in metric else ["raw","score","bounded_score","training_cv_shrinkage"]
        labels=short if "MSE" not in metric else ["Baseline","Score","Bounded","CV shrink"]
        for j,n in enumerate(["OSD-561","OSD-562","OSD-525"]):
            for i,method in enumerate(shown):
                vals=frame[(frame.accession==n)&(frame.method==method)][metric].to_numpy()
                ax.plot(i+(j-1)*.17,np.median(vals),"o",ms=3.5,color=colors[n],label=n if i==0 else None)
                ax.vlines(i+(j-1)*.17,np.quantile(vals,.025),np.quantile(vals,.975),color=colors[n],lw=.9)
        ax.axhline(1 if "MSE" in metric else 0,ls="--",lw=.8,color="#777777");ax.set_xticks(range(len(shown)),labels,rotation=25,ha="right");ax.set_ylabel(ylabel)
    axes[0,1].legend(frameon=False,fontsize=8.5)
    gains=[raw.loc[raw.accession==n,"learned_score_gain"].to_numpy() for n in ["OSD-561","OSD-562","OSD-525"]]
    boxes=axes[1,0].boxplot(gains,tick_labels=["OSD-561","OSD-562","OSD-525"],patch_artist=True,showfliers=False,widths=.45)
    for box,col in zip(boxes["boxes"],colors.values()):box.set_facecolor(col);box.set_alpha(.45)
    axes[1,0].axhline(2,ls="--",lw=.8,color="#777777");axes[1,0].set_ylabel("Training score operator spectral gain")
    for ax,title in zip(axes.flat,["a  Collection-score generalization","b  Loading prediction on held-out animals","c  Small-sample loading amplification","d  Held-out collection contrast retention"]):
        ax.set_title(title,loc="left",fontweight="bold",fontsize=10.5);ax.spines[["top","right"]].set_visible(False);ax.grid(axis="y",alpha=.15)
    fig.savefig(HERE/"ground_holdout_panel.png",dpi=600,bbox_inches="tight");fig.savefig(HERE/"ground_holdout_panel.svg",bbox_inches="tight");plt.close(fig)
    config["figure_font"]=family;(HERE/"protocol.json").write_text(json.dumps(config,indent=2)+"\n")
    print(summary[["accession","method","sensor_age_stratified_AUC_mean","outside_prediction_MSE_ratio_to_nuisance_only_mean","heldout_outside_effect_energy_removed_pct_mean"]].to_string(index=False))


if __name__=="__main__":main()
