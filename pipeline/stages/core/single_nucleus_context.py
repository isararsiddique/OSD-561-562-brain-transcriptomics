from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""OSD-352 donor-level panel context from the deposited raw multiome H5.

Reads sparse blocks to keep memory bounded. Author-provided cell labels are
used for context; nuclei are never treated as independent exposure replicates.
"""
from pathlib import Path
from collections import Counter
import itertools
import json
import hashlib
import numpy as np
import pandas as pd
import h5py
from scipy.stats import false_discovery_control

OUT=stage_root('core');ROOT=inputs_root();TAB=OUT/"tables"
BASE=ROOT/"research/brain_integrative_2026/data/OSD-352"
path=BASE/"GLDS-352_snRNA-Seq_filtered_feature_bc_matrix.h5"
meta_path=BASE/"multiomics_brain_metadata_all_cells.csv"
meta=pd.read_csv(meta_path,index_col=0)
mapping=pd.read_csv(TAB/"verified_panel_mapping.csv").fillna("")
assert meta.index.is_unique
origins=meta.groupby("bio_origin").cell_flight_status.unique()
assert all(len(x)==1 for x in origins)
assert (origins.map(lambda x:x[0])=="flight").sum()==3
assert (origins.map(lambda x:x[0])=="ground").sum()==2
groups=sorted(set(zip(meta.bio_origin,meta["predicted.celltype.l2"])))
gindex={g:i for i,g in enumerate(groups)}
n_nuclei=Counter(zip(meta.bio_origin,meta["predicted.celltype.l2"]))
group_for_barcode={b:gindex[(r.bio_origin,r["predicted.celltype.l2"])] for b,r in meta.iterrows()}

with h5py.File(path,"r") as h:
    matrix=h["matrix"]
    barcodes=matrix["barcodes"][:].astype(str)
    names=matrix["features/name"][:].astype(str)
    types=matrix["features/feature_type"][:].astype(str)
    rna=types=="Gene Expression"
    selected=np.flatnonzero(rna & np.isin(names,mapping.mouse_symbol_primary))
    gene_names=names[selected]
    assert len(set(gene_names))==len(gene_names)
    target_index=np.full(len(names),-1,dtype=int);target_index[selected]=np.arange(len(selected))
    indptr=matrix["indptr"][:]
    pb=np.zeros((len(selected),len(groups)),dtype=np.int64)
    library=np.zeros(len(groups),dtype=np.int64)
    positive=np.zeros_like(pb)
    matched=set(meta.index)&set(barcodes)
    assert len(matched)==len(meta),"Some metadata barcodes are absent from H5"
    # A block of 512 barcodes is bounded well below a full ATAC+RNA matrix.
    for start in range(0,len(barcodes),512):
        end=min(start+512,len(barcodes));lo,hi=indptr[start],indptr[end]
        indices=matrix["indices"][lo:hi]
        values=matrix["data"][lo:hi]
        columns=np.repeat(np.arange(start,end),np.diff(indptr[start:end+1]))
        groupcodes=np.array([group_for_barcode.get(b,-1) for b in barcodes[start:end]])
        codes=groupcodes[columns-start]
        valid=(codes>=0)&rna[indices]
        library+=np.bincount(codes[valid],weights=values[valid],minlength=len(groups)).astype(np.int64)
        gi=target_index[indices]
        valid=(codes>=0)&(gi>=0)
        flat=gi[valid]*len(groups)+codes[valid]
        pb+=np.bincount(flat,weights=values[valid],minlength=pb.size).reshape(pb.shape).astype(np.int64)
        positive+=np.bincount(flat,weights=(values[valid]>0).astype(int),minlength=pb.size).reshape(pb.shape).astype(np.int64)
    print("Pseudobulk",pb.shape,"nuclei",len(meta),flush=True)

pbcols=[f"{donor}|{ctype}" for donor,ctype in groups]
pd.DataFrame(pb,index=gene_names,columns=pbcols).to_csv(OUT/"prepared/OSD-352_panel_pseudobulk_counts.csv")
cols=pd.DataFrame([{"bio_origin":d,"cell_type":ct,"n_nuclei":n_nuclei[(d,ct)],
                    "flight_status":origins[d][0],"total_RNA_UMIs":library[j],"column":pbcols[j]}
                   for j,(d,ct) in enumerate(groups)])
cols.to_csv(TAB/"OSD-352_pseudobulk_registry.csv",index=False)
rows=[]
for i,gene in enumerate(gene_names):
    for j,(donor,ct) in enumerate(groups):
        rows.append({"mouse_symbol":gene,"bio_origin":donor,"cell_type":ct,"flight_status":origins[donor][0],
                     "n_nuclei":n_nuclei[(donor,ct)],"raw_count":int(pb[i,j]),
                     "log2_CPM_plus1":float(np.log2(pb[i,j]/library[j]*1e6+1)),
                     "fraction_nuclei_positive":positive[i,j]/n_nuclei[(donor,ct)]})
long=pd.DataFrame(rows)
long.to_csv(TAB/"OSD-352_panel_celltype_context.csv",index=False)
effects=[]
for r in mapping.itertuples():
    for ct in sorted(meta["predicted.celltype.l2"].unique()):
        sub=long[(long.mouse_symbol==r.mouse_symbol_primary)&(long.cell_type==ct)]
        row={"analyte":r.analyte,"mouse_symbol":r.mouse_symbol_primary,"category":r.category,"cell_type":ct,
             "n_flight_animals":0,"n_ground_animals":0,"log2_CPM_effect":np.nan,"exact_p":1.,
             "test_status":"mapping_or_gene_unavailable"}
        if len(sub):
            usable=sub[sub.n_nuclei>=20]
            row["n_flight_animals"]=int(usable.flight_status.eq("flight").sum())
            row["n_ground_animals"]=int(usable.flight_status.eq("ground").sum())
            if len(usable)==5 and row["n_flight_animals"]==3 and row["n_ground_animals"]==2:
                score=usable.log2_CPM_plus1.to_numpy();f=usable.flight_status.eq("flight").to_numpy()
                effect=float(score[f].mean()-score[~f].mean())
                null=[]
                for select in itertools.combinations(range(5),3):
                    mask=np.isin(np.arange(5),select);null.append(score[mask].mean()-score[~mask].mean())
                row.update({"log2_CPM_effect":effect,"exact_p":float(np.mean(np.abs(null)>=abs(effect)-1e-12)),
                            "test_status":"five_animal_exact_context_test"})
            else:row["test_status"]="insufficient_nuclei_per_animal_for_fixed_gate"
        effects.append(row)
effects=pd.DataFrame(effects)
effects["q_external_full_family"]=false_discovery_control(effects.exact_p)
effects.to_csv(TAB/"OSD-352_panel_celltype_exact_tests.csv",index=False)
summary={"n_nuclei_with_author_annotations":len(meta),"n_bio_origin_groups":len(origins),
         "n_flight_animals":3,"n_ground_animals":2,"minimum_nuclei_per_donor_celltype":20,
         "n_primary_mouse_symbols_in_H5":len(gene_names),"n_exact_tested_entries":int(effects.test_status.eq("five_animal_exact_context_test").sum()),
         "n_q05":int(effects.q_external_full_family.lt(.05).sum()),
         "H5_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"metadata_sha256":hashlib.sha256(meta_path.read_bytes()).hexdigest(),
         "source_note":"H5 in NASA OSD-352 file inventory; author metadata extracted from deposited multiomics metadata ZIP. Do not equate bio_origin codes with publication F1/F2/F3 labels without an explicit map.",
         "interpretation":"Separate RR-3 cell-type context; no BuOE arm. Nuclei do not increase the five-animal exposure n."}
(TAB/"OSD-352_context_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2),flush=True)
