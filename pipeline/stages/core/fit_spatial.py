from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Four-region GeoMx 0/1 factorial fits with explicit full-scale interaction."""
from pathlib import Path
import gzip
import json
import numpy as np
import pandas as pd
from scipy.stats import t, false_discovery_control

OUT=stage_root('core');ROOT=inputs_root();TAB=OUT/"tables"
mapping=pd.read_csv(TAB/"verified_panel_mapping.csv").fillna("")
vectors={"flight_in_saline":[0,1,0,0],"flight_in_BuOE":[0,1,0,1],
         "BuOE_in_ground":[0,0,1,0],"BuOE_in_flight":[0,0,1,1],
         "flight_balanced_main":[0,1,0,.5],"BuOE_balanced_main":[0,0,1,.5],
         "flight_by_BuOE_interaction":[0,0,0,1]}
allrows=[];genomes=[];metas=[]
for osd,code,region in [(682,"CA","CA1"),(685,"DG","Dentate gyrus"),(698,"FCT","Frontal cortex"),(699,"CT","Cerebral cortex")]:
    path=ROOT/"OSD-682-685-698-699-brain-spatial/data/geo"/f"GSE239336_{code}_GeneExpression_Q3norm.txt.gz"
    data=pd.read_csv(path,sep="\t",index_col=0).drop(index="NegProbe-WTX",errors="ignore")
    assert data.index.is_unique and data.shape[1]==12 and (data.to_numpy()>0).all()
    flight=np.array([float(c.startswith("FLT")) for c in data.columns])
    drug=np.array([float("_TX_" in c) for c in data.columns])
    slides=[c.split(" |")[0] for c in data.columns]
    X=np.column_stack([np.ones(len(flight)),flight,drug,flight*drug])
    assert np.linalg.matrix_rank(X)==4
    for f in [0,1]:
        for d in [0,1]:assert ((flight==f)&(drug==d)).sum()==3
    Y=np.log2(data.to_numpy())
    inv=np.linalg.inv(X.T@X);coef=Y@X@inv;resid=Y-coef@X.T;df=len(flight)-4
    s2=np.sum(resid**2,axis=1)/df
    leverage=np.sum((X@inv)*X,axis=1)
    meta=pd.DataFrame({"roi":data.columns,"study":f"OSD-{osd}","region":region,"flight":flight,
                       "BuOE":drug,"slide":slides,"inference_unit":"section/ROI; animal provenance unresolved"})
    metas.append(meta)
    for tag,vec in vectors.items():
        c=np.array(vec);effect=coef@c;se=np.sqrt(s2*(c@inv@c))
        stat=effect/se;p=2*t.sf(np.abs(stat),df)
        score=X@inv@c
        hcse=np.sqrt(np.sum((resid/(1-leverage))**2*score[None,:]**2,axis=1))
        hcp=2*t.sf(np.abs(effect/hcse),df)
        # Independent numerical identity for the interaction scale.
        if tag=="flight_by_BuOE_interaction":
            means={(f,d):Y[:,(flight==f)&(drug==d)].mean(axis=1) for f in [0,1] for d in [0,1]}
            direct=means[(1,1)]-means[(0,1)]-means[(1,0)]+means[(0,0)]
            assert np.allclose(effect,direct,atol=1e-10)
        genome=pd.DataFrame({"mouse_symbol":data.index,"study":f"OSD-{osd}","region":region,"contrast":tag,
                             "log2FoldChange":effect,"se":se,"pvalue":p,"hc3_se":hcse,"hc3_pvalue":hcp,
                             "q_genome_contrast":false_discovery_control(p)})
        genome.to_csv(TAB/f"OSD-{osd}_{tag}_genome.csv.gz",index=False,compression="gzip")
        genomes.append({"study":f"OSD-{osd}","contrast":tag,"n_genes":len(data),"n_q_genome05":int(genome.q_genome_contrast.lt(.05).sum())})
        indexed=genome.set_index("mouse_symbol")
        for r in mapping.itertuples():
            row={"study":f"OSD-{osd}","region":region,"contrast":tag,"analyte":r.analyte,
                 "human_symbol":r.human_symbol,"mouse_symbol":r.mouse_symbol_primary,"category":r.category,
                 "mapping_status":r.mapping_status,"log2FoldChange":np.nan,"se":np.nan,"pvalue":np.nan,
                 "hc3_se":np.nan,"hc3_pvalue":np.nan,"q_genome_contrast":np.nan,
                 "test_status":"mapping_unavailable" if not r.mouse_symbol_primary else "not_on_assay"}
            if r.mouse_symbol_primary in indexed.index:
                hit=indexed.loc[r.mouse_symbol_primary]
                row.update({k:float(hit[k]) for k in ["log2FoldChange","se","pvalue","hc3_se","hc3_pvalue","q_genome_contrast"]})
                row["test_status"]="section_level_tested_slide_confounded"
            allrows.append(row)
    print("SPATIAL",osd,"genes",len(data),"ROIs",len(flight),flush=True)
full=pd.DataFrame(allrows)
assert len(full)==151*7*4
full["q_panel_primary_full_family"]=false_discovery_control(full.pvalue.fillna(1))
full["q_hc3_panel_primary_full_family"]=false_discovery_control(full.hc3_pvalue.fillna(1))
full["q_panel_within_contrast_exploratory"]=full.groupby(["study","contrast"]).pvalue.transform(lambda v:false_discovery_control(v.fillna(1)))
full["ci95_low"]=full.log2FoldChange-t.ppf(.975,8)*full.se
full["ci95_high"]=full.log2FoldChange+t.ppf(.975,8)*full.se
full.to_csv(TAB/"spatial_panel_complete.csv",index=False)
pd.DataFrame(genomes).to_csv(TAB/"spatial_genome_summary.csv",index=False)
pd.concat(metas).to_csv(TAB/"spatial_sample_registry.csv",index=False)
print("SPATIAL PRIMARY FAMILY",len(full),"q<.05",int(full.q_panel_primary_full_family.lt(.05).sum()),"HC3",int(full.q_hc3_panel_primary_full_family.lt(.05).sum()),flush=True)
