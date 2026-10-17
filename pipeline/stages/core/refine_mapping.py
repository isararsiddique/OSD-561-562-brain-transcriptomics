from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Resolve coarse MGI multi-member classes with current Ensembl Compara evidence."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import requests
import pandas as pd

OUT=stage_root('core')
SRC=OUT/"sources"
table=OUT/"tables/verified_panel_mapping.csv"
m=pd.read_csv(table).fillna("")

def get_row(row):
    ens=row.human_ensembl
    path=SRC/f"ensembl_orthology_{ens}.json"
    url=f"https://rest.ensembl.org/homology/id/human/{ens}?target_species=mus_musculus&type=orthologues&sequence=none"
    if not path.exists() and path.with_suffix(".error.txt").exists() and not network_allowed():
        return ens,[],url
    if not path.exists():
        r=request_get(url,headers={"Accept":"application/json"},timeout=30)
        if r.status_code==200:
            path.write_text(json.dumps(r.json(),indent=1))
        else:
            path.with_suffix(".error.txt").write_text(f"{r.status_code}\n{r.text[:500]}")
            return ens,[],url
    obj=json.loads(path.read_text())
    return ens,[h for entry in obj.get("data",[]) for h in entry.get("homologies",[])],url

rows=list(m[m.mapping_status!="unique_mouse_member_in_MGI_orthology_class"].itertuples())
results=list(ThreadPoolExecutor(max_workers=4).map(get_row,rows))
ids=sorted({h["target"]["id"] for _,hits,_ in results for h in hits})
lookup={}
lookup_path=SRC/"ensembl_mouse_ambiguity_lookup.json"
if lookup_path.exists() and not network_allowed():
    lookup=json.loads(lookup_path.read_text())
else:
    if ids:
        r=request_post("https://rest.ensembl.org/lookup/id",json={"ids":ids},headers={"Accept":"application/json"},timeout=30)
        if r.status_code==200:
            lookup=r.json()
        else:
            def resolve(mouse_id):
                url=f"https://rest.ensembl.org/lookup/id/{mouse_id}"
                response=request_get(url,headers={"Accept":"application/json"},timeout=30)
                return mouse_id,response.json() if response.status_code==200 else {}
            lookup=dict(ThreadPoolExecutor(max_workers=3).map(resolve,ids))
        (SRC/"ensembl_mouse_ambiguity_lookup.json").write_text(json.dumps(lookup,indent=1))
ledger=[]
for ens,hits,url in results:
    candidates=[{"mouse_id":h["target"]["id"],"mouse_symbol":(lookup.get(h["target"]["id"]) or {}).get("display_name",""),
                 "type":h["type"],"high_confidence":h.get("is_high_confidence") } for h in hits]
    one=[h for h in candidates if h["type"]=="ortholog_one2one" and h["mouse_symbol"]]
    if len(one)==1 and m.loc[m.human_ensembl==ens,"analyte"].iloc[0]!="HBA1":
        m.loc[m.human_ensembl==ens,"mouse_symbol_primary"]=one[0]["mouse_symbol"]
        m.loc[m.human_ensembl==ens,"mapping_status"]="unique_one_to_one_current_Ensembl_ortholog"
        m.loc[m.human_ensembl==ens,"orthology_source"]=url
    for item in candidates:
        ledger.append({"human_ensembl":ens,**item,"source_url":url})
m.to_csv(table,index=False)
pd.DataFrame(ledger).to_csv(OUT/"tables/ensembl_orthology_ambiguity_audit.csv",index=False)
print(m.mapping_status.value_counts().to_dict())
print(m[m.mouse_symbol_primary==""][["analyte","human_symbol","mouse_candidates"]].to_string(index=False))
