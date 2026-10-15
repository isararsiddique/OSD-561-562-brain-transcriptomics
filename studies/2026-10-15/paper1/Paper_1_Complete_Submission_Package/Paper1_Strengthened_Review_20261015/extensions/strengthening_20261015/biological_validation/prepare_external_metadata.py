#!/usr/bin/env python3
"""Standardize audited acquisition metadata without altering external source files."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parent
BASE=OUT.parents[2]
ACQ=OUT.parent/'external_data'
records=[];hashes={};facts=[]
for n in [202,352,374,611]:
    folder=ACQ/f'OSD-{n}';path=folder/'metadata.csv'
    if not path.exists():continue
    snapshot=OUT/f'OSD-{n}_source_metadata_snapshot.csv';shutil.copy2(path,snapshot)
    raw=pd.read_csv(snapshot,dtype=str).fillna('')
    if not raw.sample_id.is_unique:raise ValueError('Duplicate acquired count sample')
    md=raw.copy();md['sample']=raw.sample_id;md['animal']=raw.independent_unit_id
    md['tissue']=raw.tissue.str.lower();md['age_source_value']=raw.age
    md['age']='single age';md['cohort']='matched collection';md['location_conflict']=False
    if n==352:
        md['flight']=raw.spaceflight.map({'Space Flight':'FLT','Ground Control':'GC'})
        if md.flight.isna().any():raise ValueError('Unknown RR-3 exposure value')
        kind='spaceflight';warning='Only 3+3 animals; flight and RIN ranges do not overlap; whole Brain anatomy'
    elif n==202:
        md['flight']='EARTH';kind='Earth radiation unloading factorial'
        warning='Time, ageing and collection cannot be separated; sequencing lanes partly align with design cells'
    else:
        md['flight']='EARTH';kind='Earth radiation'
        warning='Independent units are two-animal pools; 611 pool member identities unavailable' if n==611 else 'Independent units are disjoint two-animal pools'
    if 'primary_analysis_included'not in md:md['primary_analysis_included']='True'
    dest=OUT/f'OSD-{n}_standardized_metadata.csv';md.to_csv(dest,index=False)
    record=dict(accession=f'OSD-{n}',counts=str((folder/'counts.csv.gz').relative_to(BASE)),metadata=str(dest.relative_to(BASE)),
        source_metadata=str(snapshot.relative_to(BASE)),upstream_source_metadata=str(path.relative_to(BASE)),exposure_kind=kind,design_warning=warning)
    records.append(record)
    for p in [path,snapshot,folder/'counts.csv.gz',dest,*list((folder/'raw_sources').glob('*-ISA.zip'))]:
        hashes[str(p.relative_to(BASE))]=hashlib.sha256(p.read_bytes()).hexdigest()
    facts.append(dict(accession=f'OSD-{n}',n_count_columns=len(raw),n_independent_units=raw.independent_unit_id.nunique(),
        primary_included=raw.get('primary_analysis_included',pd.Series('True',index=raw.index)).astype(str).str.lower().eq('true').sum(),
        replicate_units='|'.join(sorted(raw.replicate_unit.unique())),condition_counts=json.dumps(raw.groupby('condition').size().to_dict()),
        tissue='|'.join(sorted(raw.tissue.unique())),warning=warning))
(OUT/'independent_datasets.json').write_text(json.dumps(records,indent=2)+'\n')
(OUT/'external_source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
pd.DataFrame(facts).to_csv(OUT/'external_design_facts.csv',index=False)
print(pd.DataFrame(facts).to_string(index=False))
