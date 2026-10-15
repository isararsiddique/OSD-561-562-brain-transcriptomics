#!/usr/bin/env python3
"""Independent source-count, design, normalization and dense-regression audit."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
roles=pd.read_csv(HERE/"split_animal_roles.csv",dtype={"age":str})
results=pd.read_csv(HERE/"holdout_split_results.csv")
factors=pd.read_csv(HERE/"normalization_factors.csv")
cv=pd.read_csv(HERE/"nested_training_cv.csv")
annotation=pd.read_csv(ROOT/"inputs/gene_annotation.csv.gz").set_index("ensembl_id").gene_symbol
sensor_symbols=set("Fos Fosb Junb Egr1 Egr2 Egr4 Nr4a1 Nr4a3 Dusp1 Dusp5 Dusp6 Arc Npas4 Btg2 Ier2 Gadd45b Gadd45g Klf4 Klf2 Zfp36 Hspa1a Hspa1b Hspa5 Hspb1 Hsph1 Dnajb1 Dnaja1 Hsp90aa1 Bag3 Hspa8 Hspe1".split())
checks=[];errors=[];independent_cases=[]
def check(label,condition,value=None):
    checks.append(dict(check=label,passed=bool(condition),value=value))
    if not condition:errors.append(label)

for acc,nsplit in [("OSD-561",81),("OSD-562",81),("OSD-525",200)]:
    subset=roles[roles.accession==acc]
    check(f"{acc} split count",subset.split_id.nunique()==nsplit,int(subset.split_id.nunique()))
    signatures=[]
    for sid,part in subset.groupby("split_id"):
        train=part[part.role=="train"];test=part[part.role=="test"]
        signatures.append(tuple(sorted(test.animal)))
        check(f"{acc}/{sid} disjoint animals",not bool(set(train.animal)&set(test.animal)))
        check(f"{acc}/{sid} complete source animals",part.animal.nunique()==len(part))
        if acc!="OSD-525":
            check(f"{acc}/{sid} age-cohort allocation",len(train)==8 and len(test)==4 and bool((train.groupby(["age","label"]).size()==2).all()) and bool((test.groupby(["age","label"]).size()==1).all()))
        else:check(f"{acc}/{sid} cohort allocation",dict(train.label.value_counts())=={1:6,0:5} and dict(test.label.value_counts())=={1:3,0:3})
    check(f"{acc} unique holdouts",len(set(signatures))==nsplit)
    n=acc.split("-")[1]
    raw=pd.read_csv(ROOT/f"inputs/OSD-{n}_counts.csv.gz",index_col=0)
    md=pd.read_csv(ROOT/f"inputs/OSD-{n}_metadata.csv.gz",dtype=str).set_index("sample").loc[raw.columns]
    counts=raw.T.groupby(md["Source Name"].to_numpy(),sort=True).sum().T
    for sid in [0,nsplit-1]:
        part=subset[subset.split_id==sid]
        tr=part[part.role=="train"].sort_values("animal");te=part[part.role=="test"].sort_values("animal")
        chosen=counts.index.astype(str).str.match(r"^ENSMUSG\d+$")&(counts[tr.animal].mean(axis=1)>=10)
        c=counts.loc[chosen]
        pos=(c[tr.animal]>0).all(axis=1)
        ref=np.exp(np.log(c.loc[pos,tr.animal]).mean(axis=1)).to_numpy()
        def calc(names):
            out=[]
            for name in names:
                values=c.loc[pos,name].to_numpy(float);nz=values>0
                out.append(np.median(values[nz]/ref[nz]))
            # Median ratio and exponentiated median log ratio differ for an
            # even number of reference genes: production defines log median.
            exact=[]
            for name in names:
                values=c.loc[pos,name].to_numpy(float);nz=values>0
                ratios=values[nz]/ref[nz]
                ratios.sort();j=len(ratios)//2
                exact.append(ratios[j] if len(ratios)%2 else np.sqrt(ratios[j-1]*ratios[j]))
            return np.log2(c[names].to_numpy()/np.array(exact)+1),np.array(exact)
        yt,sft=calc(tr.animal.tolist());ye,sfe=calc(te.animal.tolist())
        altered=c.copy();altered.loc[:,te.animal]=altered.loc[:,te.animal]*3
        altered_ref=np.exp(np.log(altered.loc[pos,tr.animal]).mean(axis=1)).to_numpy()
        check(f"{acc}/{sid} target perturbation leaves training reference unchanged",np.array_equal(altered_ref,ref))
        actual=factors[(factors.accession==acc)&(factors.split_id==sid)].set_index("animal").size_factor
        discrepancy=max(np.max(abs(sft-actual.loc[tr.animal].to_numpy())),np.max(abs(sfe-actual.loc[te.animal].to_numpy())))
        check(f"{acc}/{sid} independent train-reference factors",discrepancy<1e-10,float(discrepancy))
        xt=np.column_stack([np.ones(len(tr)),(tr.age=="29").to_numpy(float)]) if acc!="OSD-525" else np.ones((len(tr),1))
        xe=np.column_stack([np.ones(len(te)),(te.age=="29").to_numpy(float)]) if acc!="OSD-525" else np.ones((len(te),1))
        beta=np.linalg.lstsq(np.column_stack([xt,tr.label.to_numpy()]),yt.T,rcond=None)[0][-1]
        v=beta/np.linalg.norm(beta)
        sensor=np.flatnonzero(annotation.reindex(c.index).isin(sensor_symbols).to_numpy())
        outside=np.setdiff1d(np.arange(len(c)),sensor)
        w=np.zeros(len(c));w[sensor]=beta[sensor];w/=np.linalg.norm(w)
        B=np.linalg.lstsq(xt,yt.T,rcond=None)[0]
        resid_train=yt-(xt@B).T;s=w@resid_train
        gamma=resid_train@s/(s@s)
        resid_test=ye-(xe@B).T
        rawbeta=np.linalg.lstsq(np.column_stack([xe,te.label.to_numpy()]),ye.T,rcond=None)[0][-1]
        rows=results[(results.accession==acc)&(results.split_id==sid)].set_index("method")
        check(f"{acc}/{sid} dense-regression learned gain",abs(np.linalg.norm(gamma)-rows.iloc[0].learned_score_gain)<1e-8)
        check(f"{acc}/{sid} exact sensor identity",abs(w@gamma-1)<1e-10)
        score=w@resid_test
        for mode,row in rows.iterrows():
            if mode=="raw":b=rawbeta;pred=np.zeros_like(gamma)
            elif mode=="projection":b=rawbeta-v*(v@rawbeta);pred=None
            else:
                if mode=="bounded_score":
                    h=gamma-w;pred=w+h*min(1.,np.sqrt(3)/np.linalg.norm(h))
                elif mode=="training_cv_shrinkage":pred=gamma*row.cv_selected_shrinkage
                else:pred=gamma
                b=rawbeta-pred*(w@rawbeta)
            e=100*(1-np.dot(b[outside],b[outside])/np.dot(rawbeta[outside],rawbeta[outside]))
            check(f"{acc}/{sid}/{mode} coefficient-space energy",abs(e-row.heldout_outside_effect_energy_removed_pct)<1e-7,float(abs(e-row.heldout_outside_effect_energy_removed_pct)))
            if pred is not None:
                mse=np.mean((resid_test[outside]-np.outer(pred[outside],score))**2)
                check(f"{acc}/{sid}/{mode} independent prediction MSE",abs(mse-row.outside_prediction_MSE)<1e-8)
        independent_cases.append(dict(accession=acc,split_id=sid,n_genes=len(c),n_training=len(tr),n_test=len(te),normalization_factor_max_discrepancy=float(discrepancy)))
for (acc,sid),part in cv.groupby(["accession","split_id"]):
    part=part.sort_values("shrinkage");best=part.iloc[np.argmin(part.inner_error_ratio.to_numpy())]
    check(f"{acc}/{sid} training CV argmin",part.selected.sum()==1 and bool(best.selected))
check("projection excluded from disjoint-predictor MSE",results.loc[results.method=="projection","outside_prediction_MSE"].isna().all())
report=dict(passed=not errors,n_checks=len(checks),n_failed=len(errors),failures=errors,independently_recomputed_cases=independent_cases,
            unique_source_animals_all_regions=int(roles.animal.nunique()),unique_source_animals_per_accession=roles.groupby("accession").animal.nunique().astype(int).to_dict(),
            scope="Source-count collapse; split disjointness/allocation; normalization factors; dense coefficient estimates and prediction errors; inner CV argmin. Numerical verification does not prove biological validity.")
(HERE/"verification_report.json").write_text(json.dumps(report,indent=2)+"\n")
pd.DataFrame(checks).to_csv(HERE/"verification_checks.csv",index=False)
print(json.dumps(report,indent=2))
if errors:raise SystemExit(1)
