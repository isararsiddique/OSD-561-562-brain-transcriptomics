from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Primary-source identity, metadata, gene-mapping and checksum preparation."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import csv
import io
import json
import hashlib
import re
import zipfile
import requests
import pandas as pd
import numpy as np

ROOT = inputs_root()
OUT = stage_root('core')
SRC = OUT / "sources"
TAB = OUT / "tables"
for p in [SRC, TAB, OUT / "prepared", OUT / "figures", OUT / "logs"]:
    p.mkdir(exist_ok=True)


def download(url, name):
    path = SRC / name
    if path.exists() and path.stat().st_size:
        return path
    response = request_get(url, headers={"Accept": "application/json", "User-Agent": "CNS-panel-source-audit/1.0"}, timeout=35)
    response.raise_for_status()
    path.write_bytes(response.content)
    return path


def md5_record(path, sidecar):
    records = dict(line.split("\t")[:2] for line in sidecar.read_text().splitlines() if "\t" in line)
    suffix = path.name.split("_rna_seq_", 1)[-1]
    expected = records.get(suffix)
    actual = hashlib.md5(path.read_bytes()).hexdigest()
    return {"path": str(path.relative_to(ROOT)), "expected_md5": expected,
            "actual_md5": actual, "matches_deposited_md5": expected == actual,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


retrieval = []
urls = [(f"https://osdr.nasa.gov/osdr/data/osd/meta/{n}", f"OSD-{n}_meta.json") for n in [561,562,682,685,698,699,352]]
urls += [(f"https://osdr.nasa.gov/osdr/data/osd/files/{n}", f"OSD-{n}_files.json") for n in [561,562,352]]
with ThreadPoolExecutor(max_workers=4) as pool:
    futs = {pool.submit(download, url, name): (url, name) for url, name in urls}
    for future in as_completed(futs):
        url, name = futs[future]
        try:
            path = future.result()
            retrieval.append({"url": url, "file": str(path.relative_to(OUT)), "status": "retrieved", "date": "2026-10-16"})
        except Exception as exc:
            retrieval.append({"url": url, "file": name, "status": repr(exc), "date": "2026-10-16"})
pd.DataFrame(retrieval).to_csv(TAB / "primary_retrieval_log.csv", index=False)
print("NASA source retrieval", {r["file"]:r["status"] for r in retrieval}, flush=True)

# MGI's primary orthology report gives human-symbol orthology classes and mouse members.
mgi_url = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
mgi = pd.read_csv(download(mgi_url, "MGI_HOM_MouseHumanSequence.rpt"), sep="\t", dtype=str).fillna("")
classes = {key: g for key,g in mgi.groupby("DB Class Key")}
human_mouse = defaultdict(set)
mouse_mgi_ids = defaultdict(set)
for key, g in classes.items():
    humans = g.loc[g["NCBI Taxon ID"] == "9606", "Symbol"]
    mice = g.loc[g["NCBI Taxon ID"] == "10090"]
    for h in humans:
        human_mouse[h].update(mice.Symbol.tolist())
    for _,row in mice.iterrows():
        mouse_mgi_ids[row.Symbol].add(row["Mouse MGI ID"])

old = pd.read_csv(ROOT / "config/cns_targets_mapped.csv").fillna("")
# Current human ENSG resolution: batch lookup saved as a primary record.
lookup_path = SRC / "ensembl_human_lookup_current.json"
if not lookup_path.exists():
    try:
        resp = request_post("https://rest.ensembl.org/lookup/id", json={"ids": old.ensid.tolist()},
                             headers={"Accept":"application/json","Content-Type":"application/json"},timeout=45)
        resp.raise_for_status()
        lookup_path.write_text(json.dumps(resp.json(),indent=1))
    except Exception as exc:
        (SRC / "ensembl_lookup_error.txt").write_text(repr(exc))
lookup = json.loads(lookup_path.read_text()) if lookup_path.exists() else {}
mapped = []
for _, row in old.iterrows():
    current = lookup.get(row.ensid) or {}
    human = current.get("display_name", row.human_symbol_ensembl)
    candidates = sorted(human_mouse.get(human, set()))
    flags = []
    if not current:
        flags.append("Human ENSG resolution uses saved Ensembl cache; current batch lookup unavailable")
    if row.analyte == "HBA1" and human == "HBA2":
        flags.append("Source label HBA1 points to HBA2; retain discrepancy")
    if len(candidates) == 1:
        mouse = candidates[0]
        status = "unique_mouse_member_in_MGI_orthology_class"
    elif not candidates:
        mouse = ""
        status = "no_mouse_member_in_MGI_report"
    else:
        mouse = ""
        status = "multiple_mouse_members_in_MGI_orthology_class"
    mapped.append({"analyte":row.analyte,"human_ensembl":row.ensid,"human_symbol":human,
                   "category":row.category,"mouse_symbol_primary":mouse,"mapping_status":status,
                   "mouse_candidates":";".join(candidates),"previous_mouse_symbol":row.mouse_ortholog_symbol,
                   "mouse_mgi_ids":";".join(sorted(mouse_mgi_ids[mouse])) if mouse else "",
                   "mapping_notes":"; ".join(flags),"orthology_source":mgi_url})
mapping = pd.DataFrame(mapped)
mapping.to_csv(TAB / "verified_panel_mapping.csv",index=False)
print("MGI primary mapping",mapping.mapping_status.value_counts().to_dict(),flush=True)

checks = []
allmeta = []
for osd, glds, region in [(561,556,"Cerebellum"),(562,557,"Hippocampus")]:
    base=ROOT / "data" / f"OSD-{osd}"
    star=base/f"GLDS-{glds}_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv"
    side=base/f"GLDS-{glds}_rna_seq_processed_md5sum_GLbulkRNAseq.tsv"
    check=md5_record(star,side);checks.append(check)
    assert check["matches_deposited_md5"], f"Primary STAR checksum mismatch: {star}"
    counts=pd.read_csv(star,index_col=0)
    assert counts.index.is_unique and counts.columns.is_unique
    assert (counts.to_numpy()>=0).all() and np.equal(counts,np.floor(counts)).all().all()
    mask=counts.index.str.match(r"^ENSMUSG\d+(?:\.\d+)?$")
    counts.loc[~mask].sum(axis=1).rename("total_count").to_csv(TAB/f"OSD-{osd}_excluded_controls.csv")
    counts=counts.loc[mask]
    zip_path=base/f"OSD-{osd}_metadata_OSD-{osd}-ISA.zip"
    with zipfile.ZipFile(zip_path) as z:
        name=next(n for n in z.namelist() if n.split("/")[-1].startswith("s_"))
        isa=pd.read_csv(io.BytesIO(z.read(name)),sep="\t",dtype=str).set_index("Sample Name")
        (SRC/f"OSD-{osd}_sample_isa.tsv").write_bytes(z.read(name))
    missing=set(counts.columns)-set(isa.index)
    assert not missing,f"Missing sample metadata {missing}"
    meta=isa.loc[counts.columns].copy()
    meta.index.name="sample"
    meta["animal"]=meta["Source Name"]
    meta["flight"]=meta["Factor Value[Spaceflight]"].map({"Space Flight":"FLT","Ground Control":"GC"})
    meta["age"]=meta["Factor Value[Age]"].astype(int).astype(str)
    meta["cohort"]=np.where(meta.animal.str.contains("_LAR_"),"LAR","ISS_T")
    meta["group"]=meta.flight+"_"+meta.cohort+"_"+meta.age
    for s,r in meta.iterrows():
        assert ("_LAR_" in s)==(r.cohort=="LAR")
        assert ("_ISS-T_" in s)==(r.cohort=="ISS_T")
    expected=np.where((meta.flight=="FLT")&(meta.cohort=="ISS_T"),"On ISS","On Earth")
    meta["source_expected_location"]=expected
    meta["location_conflict"]=meta["Factor Value[Euthanasia Location]"]!=expected
    meta["region"]=region
    meta["study"]=f"OSD-{osd}"
    meta["endogenous_library_size"]=counts.sum(axis=0)
    meta["detected_genes"]=counts.gt(0).sum(axis=0)
    # Strictly collapse repeated full sources only within this region, if any.
    if not meta.animal.is_unique:
        counts=counts.T.groupby(meta.animal,sort=False).sum().T
        meta=meta.reset_index().drop_duplicates("animal").set_index("animal")
    meta.to_csv(TAB/f"OSD-{osd}_sample_metadata.csv")
    allmeta.append(meta.reset_index())
    counts.to_csv(OUT/"prepared"/f"OSD-{osd}_STAR_endogenous.csv.gz",compression="gzip")
    annotation=pd.read_csv(base/f"GLDS-{glds}_rna_seq_differential_expression_GLbulkRNAseq.csv",
                          usecols=["ENSEMBL","SYMBOL","GENENAME"]).drop_duplicates("ENSEMBL")
    annotation.to_csv(OUT/"prepared"/f"OSD-{osd}_annotation.csv",index=False)
    print(f"Prepared OSD-{osd}",counts.shape,meta.group.value_counts().to_dict(),flush=True)
pd.DataFrame(checks).to_csv(TAB/"primary_checksum_verification.csv",index=False)
meta=pd.concat(allmeta,ignore_index=True)
meta.to_csv(TAB/"bulk_sample_registry.csv",index=False)
meta.loc[meta.location_conflict].to_csv(TAB/"bulk_metadata_conflicts.csv",index=False)
links=meta.groupby("animal").agg(n_regions=("study","nunique"),studies=("study",lambda v:";".join(sorted(v.unique())))).reset_index()
links.to_csv(TAB/"bulk_animal_links.csv",index=False)
summary={"libraries":len(meta),"unique_animal_sources":int(meta.animal.nunique()),
         "paired_animal_sources":int(links.n_regions.eq(2).sum()),"location_conflicts":int(meta.location_conflict.sum()),
         "primary_mapping_counts":mapping.mapping_status.value_counts().to_dict()}
(OUT/"prepared"/"identity_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2),flush=True)
