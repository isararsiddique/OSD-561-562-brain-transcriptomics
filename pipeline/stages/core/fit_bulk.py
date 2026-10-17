from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Full endogenous count fits and fixed-family corrected-panel inference."""
from pathlib import Path
import argparse
import json
import pickle
import importlib.metadata
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

OUT=stage_root('core')
ROOT=inputs_root()
TAB=OUT/"tables"
PRE=OUT/"prepared"

def bh(v):
    a=np.asarray(v,dtype=float)
    return false_discovery_control(np.where(np.isfinite(a),a,1.0))

def contrasts(dds):
    # Differences of observed design rows avoid formula reconstruction on reload.
    design=np.asarray(dds.obsm["design_matrix"])
    groups=dds.obs["group"].astype(str).to_numpy()
    def vector(numerator,denominator):
        assert (groups==numerator).any() and (groups==denominator).any()
        return design[groups==numerator].mean(axis=0)-design[groups==denominator].mean(axis=0)
    items=[]
    for cohort in ["ISS_T","LAR"]:
        for age in ["12","29"]:
            items.append((f"flight_{cohort}_age{age}","flight",f"FLT_{cohort}_{age}",f"GC_{cohort}_{age}"))
    for flight in ["GC","FLT"]:
        for cohort in ["ISS_T","LAR"]:
            items.append((f"age_{flight}_{cohort}","age",f"{flight}_{cohort}_29",f"{flight}_{cohort}_12"))
        for age in ["12","29"]:
            items.append((f"cohort_{flight}_age{age}","cohort",f"{flight}_ISS_T_{age}",f"{flight}_LAR_{age}"))
    out=[]
    for tag,kind,num,den in items:
        v=vector(num,den)
        out.append((tag,kind,np.asarray(v),num,den))
    for age in ["12","29"]:
        iss=vector(f"FLT_ISS_T_{age}",f"GC_ISS_T_{age}")
        lar=vector(f"FLT_LAR_{age}",f"GC_LAR_{age}")
        out.append((f"flight_by_cohort_age{age}","interaction",np.asarray(iss-lar),"(FLT-GC)ISS_T","(FLT-GC)LAR"))
    assert len(out)==14
    return out

def run(osd,mode):
    glds={561:556,562:557}[osd]
    meta=pd.read_csv(TAB/f"OSD-{osd}_sample_metadata.csv",index_col=0)
    if mode=="STAR":
        counts=pd.read_csv(PRE/f"OSD-{osd}_STAR_endogenous.csv.gz",index_col=0).T.astype("int64")
    else:
        counts=pd.read_csv(ROOT/f"data/OSD-{osd}/GLDS-{glds}_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",index_col=0).T
        counts=counts.loc[:,counts.columns.str.match(r"^ENSMUSG\d+(?:\.\d+)?$")].round().astype("int64")
    counts=counts.loc[meta.index]
    keep=counts.ge(10).sum(axis=0).ge(3)
    counts=counts.loc[:,keep]
    print("FIT",osd,mode,counts.shape,meta.group.value_counts().to_dict(),flush=True)
    model_path=PRE/f"OSD-{osd}_{mode}_dds.pkl"
    if model_path.exists():
        with model_path.open("rb") as f:dds=pickle.load(f)
        assert pd.Index(dds.obs_names).equals(counts.index)
        assert pd.Index(dds.var_names).equals(counts.columns)
        assert np.array_equal(dds.X,counts.to_numpy()),'Cached model input counts changed'
        assert np.array_equal(dds.obs.group.astype(str),meta.loc[dds.obs_names,'group'].astype(str)),'Cached model group metadata changed'
    else:
        dds=DeseqDataSet(counts=counts,metadata=meta[["group"]].copy(),design="~ group",refit_cooks=True,n_cpus=threads(),quiet=True)
        dds.deseq2()
        with model_path.open("wb") as f:pickle.dump(dds,f)
    if mode=="STAR":
        pd.DataFrame(dds.layers["normed_counts"],index=dds.obs_names,columns=dds.var_names).T.to_csv(PRE/f"OSD-{osd}_normalized_counts.csv.gz",compression="gzip")
        pd.DataFrame({"sample":dds.obs_names,"size_factor":dds.obs["size_factors"].to_numpy()}).to_csv(TAB/f"OSD-{osd}_size_factors.csv",index=False)
    annotation=pd.read_csv(PRE/f"OSD-{osd}_annotation.csv").set_index("ENSEMBL")
    mapping=pd.read_csv(TAB/"verified_panel_mapping.csv").fillna("")
    rows=[];summaries=[]
    all_tests=contrasts(dds)
    if mode!="STAR":all_tests=[c for c in all_tests if c[1]=="flight"]
    for tag,kind,vector,num,den in all_tests:
        stats=DeseqStats(dds,contrast=vector,alpha=.05,quiet=True,n_cpus=threads(),independent_filter=False)
        stats.summary()
        res=stats.results_df.join(annotation)
        res.index.name="ENSEMBL"
        res["genome_q_with_unavailable_p1"]=bh(res.pvalue)
        res.to_csv(TAB/f"OSD-{osd}_{mode}_{tag}_genome.csv.gz",compression="gzip")
        for r in mapping.itertuples():
            symbol=r.mouse_symbol_primary
            hits=res.loc[res.SYMBOL.eq(symbol)] if symbol else res.iloc[:0]
            record={"study":f"OSD-{osd}","count_mode":mode,"contrast":tag,"contrast_kind":kind,
                    "numerator":num,"denominator":den,"analyte":r.analyte,"human_symbol":r.human_symbol,
                    "mouse_symbol":symbol,"category":r.category,"mapping_status":r.mapping_status,
                    "gene_id":"","baseMean":np.nan,"log2FoldChange":np.nan,"lfcSE":np.nan,
                    "stat":np.nan,"pvalue":np.nan,"q_genome_contrast":np.nan,
                    "test_status":"mapping_unavailable" if not symbol else "count_filter_or_annotation_unavailable"}
            if len(hits)==1:
                hit=hits.iloc[0]
                record.update({k:float(hit[k]) for k in ["baseMean","log2FoldChange","lfcSE","stat","pvalue"]})
                record["q_genome_contrast"]=float(hit.genome_q_with_unavailable_p1)
                record["gene_id"]=hits.index[0]
                record["test_status"]="tested" if np.isfinite(hit.pvalue) else "Cooks_or_model_pvalue_unavailable"
            elif len(hits)>1:
                record["test_status"]="multiple_gene_ids_for_symbol"
            rows.append(record)
        summaries.append({"study":f"OSD-{osd}","count_mode":mode,"contrast":tag,"kind":kind,
                          "n_animals":len(meta),"n_filtered_endogenous_genes":len(res),
                          "n_finite_p":int(res.pvalue.notna().sum()),"n_genome_q05":int(res.genome_q_with_unavailable_p1.lt(.05).sum())})
        print("CONTRAST",osd,mode,tag,summaries[-1]["n_genome_q05"],flush=True)
    panel=pd.DataFrame(rows)
    panel["q_panel_within_contrast_exploratory"]=panel.groupby("contrast").pvalue.transform(lambda v:bh(v))
    panel.to_csv(TAB/f"OSD-{osd}_{mode}_panel.csv",index=False)
    pd.DataFrame(summaries).to_csv(TAB/f"OSD-{osd}_{mode}_model_summary.csv",index=False)

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--mode",default="STAR",choices=["STAR","rounded_RSEM"])
    parser.add_argument("--study",type=int,choices=[561,562]);args=parser.parse_args()
    for study in ([args.study] if args.study else [561,562]):run(study,args.mode)
    if all((TAB/f"OSD-{n}_STAR_panel.csv").exists() for n in [561,562]):
        full=pd.concat([pd.read_csv(TAB/f"OSD-{n}_STAR_panel.csv") for n in [561,562]],ignore_index=True)
        assert len(full)==151*14*2
        full["q_panel_primary_full_family"]=bh(full.pvalue)
        full["q_dependency_BY_sensitivity"]=false_discovery_control(full.pvalue.fillna(1),method='by')
        full["ci95_low"]=full.log2FoldChange-1.96*full.lfcSE
        full["ci95_high"]=full.log2FoldChange+1.96*full.lfcSE
        full.to_csv(TAB/"bulk_panel_complete.csv",index=False)
        full[full.q_panel_primary_full_family.lt(.05)].to_csv(TAB/"bulk_panel_primary_q05.csv",index=False)
        print("PRIMARY FULL FAMILY",len(full),"q<.05",int(full.q_panel_primary_full_family.lt(.05).sum()),flush=True)
    (OUT/"logs/model_versions.json").write_text(json.dumps({m:importlib.metadata.version(m) for m in ["pydeseq2","numpy","pandas","scipy","anndata"]},indent=2))
