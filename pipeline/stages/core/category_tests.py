from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Animal-level category/null tests, stability and paired-region associations."""
from pathlib import Path
import itertools
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, t

OUT=stage_root('core');TAB=OUT/"tables";PRE=OUT/"prepared"
M=pd.read_csv(TAB/"verified_panel_mapping.csv").fillna("")
SEED=20261016
records=[];points=[];stability=[];Y_by_study={};meta_by_study={}

def assignments(n,k):
    return np.array([[float(j in selected) for j in range(n)] for selected in itertools.combinations(range(n),k)])

for osd in [561,562]:
    meta=pd.read_csv(TAB/f"OSD-{osd}_sample_metadata.csv",index_col=0)
    norm=pd.read_csv(PRE/f"OSD-{osd}_normalized_counts.csv.gz",index_col=0)
    anno=pd.read_csv(PRE/f"OSD-{osd}_annotation.csv").set_index("ENSEMBL")
    sym=anno.SYMBOL.reindex(norm.index)
    # Retain uniquely annotated mapped symbols, with no label-dependent filter.
    eligible=set(M.mouse_symbol_primary)-{""}
    mask=sym.isin(eligible)&~sym.duplicated(keep=False)
    Y=np.log2(norm.loc[mask]+1)
    Y.index=sym.loc[mask]
    Y=Y.loc[:,meta.index]
    sd=Y.std(axis=1,ddof=1)
    Z=Y.loc[sd>0].sub(Y.mean(axis=1),axis=0).div(sd[sd>0],axis=0)
    Y_by_study[osd]=Y;meta_by_study[osd]=meta
    category={r.mouse_symbol_primary:r.category for r in M.itertuples() if r.mouse_symbol_primary}
    categories={cat:[g for g in Z.index if category[g]==cat] for cat in sorted(set(category.values()))}
    categories["Whole corrected panel"]=Z.index.tolist()
    for cat,genes in categories.items():
        if genes:
            score=Z.loc[genes].mean(axis=0)
            for sample,v in score.items():
                points.append({"study":f"OSD-{osd}","sample":sample,"animal":meta.loc[sample,"animal"],
                               "group":meta.loc[sample,"group"],"category":cat,"standardized_score":float(v),"n_genes":len(genes)})
    for cohort in ["ISS_T","LAR"]:
        for age in ["12","29"]:
            subset=meta[(meta.cohort==cohort)&(meta.age.astype(str)==age)]
            samples=subset.index
            isflt=subset.flight.eq("FLT").to_numpy()
            n,k=len(samples),int(isflt.sum())
            assert min(k,n-k)>=2
            A=assignments(n,k)
            W=A/k-(1-A)/(n-k)
            w=isflt/k-(1-isflt)/(n-k)
            for cat,genes in categories.items():
                if not genes:
                    records.append({"study":f"OSD-{osd}","cohort":cohort,"age":age,"category":cat,
                                    "n_genes":0,"n_flight":k,"n_ground":n-k,"effect_standardized":np.nan,"exact_p":1.,"n_assignments":len(A)})
                    continue
                score=Z.loc[genes,samples].mean(axis=0).to_numpy()
                effect=float(w@score);null=W@score
                p=float(np.mean(np.abs(null)>=abs(effect)-1e-12))
                records.append({"study":f"OSD-{osd}","cohort":cohort,"age":age,"category":cat,
                                "n_genes":len(genes),"n_flight":k,"n_ground":n-k,"effect_standardized":effect,
                                "exact_p":p,"n_assignments":len(A)})
            genevals=Z.loc[:,samples].to_numpy()
            effect=genevals@w
            null=genevals@W.T
            observed=float(np.mean(effect**2));nullstat=np.mean(null**2,axis=0)
            p=float(np.mean(nullstat>=observed-1e-12))
            records.append({"study":f"OSD-{osd}","cohort":cohort,"age":age,"category":"Panel direction-free omnibus",
                            "n_genes":len(Z),"n_flight":k,"n_ground":n-k,"effect_standardized":observed,
                            "exact_p":p,"n_assignments":len(A)})
            pd.DataFrame({"null_statistic":nullstat}).to_csv(TAB/f"OSD-{osd}_omnibus_null_{cohort}_age{age}.csv",index=False)
            # Delete-one-animal effect range and direction stability for every eligible gene.
            observed_raw=Y.loc[:,samples].to_numpy()@w
            deleted=[]
            for j in range(n):
                keep=np.arange(n)!=j;fl=isflt[keep];ww=fl/fl.sum()-(1-fl)/(~fl).sum()
                deleted.append(Y.loc[:,samples].to_numpy()[:,keep]@ww)
            deleted=np.array(deleted).T
            for idx,gene in enumerate(Y.index):
                stability.append({"study":f"OSD-{osd}","contrast":f"flight_{cohort}_age{age}","mouse_symbol":gene,
                                  "log_normalized_effect":observed_raw[idx],"loo_effect_min":float(deleted[idx].min()),
                                  "loo_effect_max":float(deleted[idx].max()),
                                  "loo_same_direction_fraction":float(np.mean(np.sign(deleted[idx])==np.sign(observed_raw[idx])))})

res=pd.DataFrame(records)
assert len(res)==2*4*12
res["q_category_full_family"]=false_discovery_control(res.exact_p)
res.to_csv(TAB/"bulk_category_exact_tests.csv",index=False)
pd.DataFrame(points).to_csv(TAB/"bulk_category_animal_scores.csv",index=False)
pd.DataFrame(stability).to_csv(TAB/"bulk_flight_leave_one_out.csv",index=False)

# Shared animals: separate regional normalization followed by within-animal differences.
common=sorted(set(meta_by_study[561].animal)&set(meta_by_study[562].animal))
genes=sorted(set(Y_by_study[561].index)&set(Y_by_study[562].index))
region_values={}
for osd in [561,562]:
    meta=meta_by_study[osd]
    region_values[osd]=Y_by_study[osd].loc[genes].rename(columns=meta.animal.to_dict()).loc[:,common]
per=meta_by_study[561].set_index("animal").loc[common]
diff=region_values[561].to_numpy()-region_values[562].to_numpy()
stratum=per.cohort+"_"+per.age.astype(str)
N=pd.get_dummies(stratum,drop_first=False,dtype=float).to_numpy()
P=np.eye(len(per))-N@np.linalg.pinv(N)
flight=per.flight.eq("FLT").to_numpy(dtype=float)
x=P@flight;Yr=diff@P
beta=Yr@x/(x@x)
df=len(per)-np.linalg.matrix_rank(N)-1
resid=Yr-beta[:,None]*x
se=np.sqrt(np.sum(resid**2,axis=1)/df/(x@x))
pv=2*t.sf(abs(beta/se),df)
paired=pd.DataFrame({"mouse_symbol":genes,"cerebellum_minus_hippocampus_flight_effect":beta,
                     "se":se,"pvalue":pv,"q_secondary_genes":false_discovery_control(pv),"n_paired_animals":len(per)})
paired.to_csv(TAB/"paired_region_flight_associations.csv",index=False)
# Correlated-gene global label null, blocking on age/cohort and moving whole animals.
rng=np.random.default_rng(SEED);B=19999
Xperm=np.empty((len(per),B))
for b in range(B):
    v=flight.copy()
    for level in stratum.unique():
        ix=np.flatnonzero(stratum.to_numpy()==level);v[ix]=rng.permutation(v[ix])
    Xperm[:,b]=P@v
den=np.sum(Xperm**2,axis=0)
null_beta=Yr@Xperm/den[None,:]
observed=float(np.mean(beta**2));null=np.mean(null_beta**2,axis=0)
pd.DataFrame({'mean_squared_adjusted_regional_flight_effect':null}).to_csv(TAB/'paired_region_omnibus_null.csv',index=False)
summary={"n_paired_animals":len(per),"n_genes":len(genes),"omnibus_mean_squared_regional_flight_effect":observed,
         "permutations":B,"seed":SEED,"pvalue":float((1+np.sum(null>=observed-1e-12))/(B+1)),
         "interpretation":"Secondary conditional association; shared animals, no independent replication."}
(TAB/"paired_region_omnibus.json").write_text(json.dumps(summary,indent=2)+"\n")
print("CATEGORY FULL FAMILY q<.05",int(res.q_category_full_family.lt(.05).sum()),flush=True)
print(json.dumps(summary,indent=2),flush=True)
