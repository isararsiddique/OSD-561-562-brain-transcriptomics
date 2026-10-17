from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Cerebral hemispheres: explicit shared-versus-new animal audit."""
from pathlib import Path
import json,hashlib,zipfile,sys
import requests,numpy as np,pandas as pd
B=stage_root('deep')
old=pd.read_csv(stage_root('extended') / 'tables/validation_sample_registry.csv')
core=pd.read_csv(stage_root('core') / 'tables/bulk_sample_registry.csv')
out=[];provenance=[];overlaps=[];missing=[]
saved=B/'sources/terminal_region_download_provenance.json'
expected={Path(x['path']).name:x['SHA256'] for x in json.loads(saved.read_text())} if saved.exists() else {}
for sid,glds in [(612,588),(613,589)]:
    fs=json.loads((B/'sources'/f'OSD-{sid}_files.json').read_text())['studies'][f'OSD-{sid}']['study_files']
    for f in fs:
        if any(t in f['file_name'] for t in ['ISA.zip','STAR_Unnormalized_Counts','runsheet']):
            url='https://osdr.nasa.gov'+f['remote_url'];p=B/'sources'/f['file_name']
            if not p.exists() or '--refresh' in sys.argv:
                r=request_get(url,timeout=(20,100));r.raise_for_status();assert len(r.content)==f['file_size'];p.write_bytes(r.content)
            content=p.read_bytes();digest=hashlib.sha256(content).hexdigest();assert len(content)==f['file_size']
            if p.name in expected and '--refresh' not in sys.argv:assert digest==expected[p.name],p.name
            provenance.append({'path':record_path(p),'URL':url,'bytes':p.stat().st_size,'SHA256':digest})
    with zipfile.ZipFile(B/'sources'/f'OSD-{sid}_metadata_OSD-{sid}-ISA.zip') as z:
        assert z.testzip() is None;z.extract(f's_OSD-{sid}.txt',B/'sources')
    isa=pd.read_csv(B/'sources'/f's_OSD-{sid}.txt',sep='\t')
    x=pd.read_csv(B/'sources'/f'GLDS-{glds}_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv',index_col=0)
    assert isa['Sample Name'].is_unique and x.columns.isin(isa['Sample Name']).all() and x.index.is_unique
    for _,absent in isa[~isa['Sample Name'].isin(x.columns)].iterrows():
        missing.append({'OSD':f'OSD-{sid}','sample':absent['Sample Name'],'donor':absent['Source Name'],
                        'status':'Primary ISA sample has no deposited bulk STAR count column; not imputed'})
    isa=isa[isa['Sample Name'].isin(x.columns)]
    assert np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all() and np.equal(x.to_numpy(),np.floor(x.to_numpy())).all()
    x=x.loc[x.index.str.startswith('ENSMUSG'),isa['Sample Name']]
    compare=set(old.loc[old.study.isin(['OSD-563','OSD-564']),'animal']) if sid==612 else set(core.animal)
    current=set(isa['Source Name']);overlap=current&compare
    overlaps.append({'OSD':f'OSD-{sid}','n_animals':len(current),'overlap_with_original_region_animals':len(overlap),
                     'overlapping_source_names':sorted(overlap),'independence':'same RR10 animals, new cerebral-region measurement' if sid==612 else 'different RRRM2 animals, same mission/cohort and different region',
                     'new_mice_added':len(current-compare),'source':'Exact primary ISA Source Name identity, no ordinal-ID guessing'})
    for _,s in isa.iterrows():
        agecol='Factor Value[Age]' if sid==613 else 'Characteristics[Age at Launch]'
        age=int(s[agecol]) if sid==613 else np.nan
        flight='FLT' if s['Factor Value[Spaceflight]']=='Space Flight' else 'GC'
        out.append({'OSD':f'OSD-{sid}','key':s['Sample Name'],'donor':s['Source Name'],'flight':flight,'age':age,'launch_age_source':str(s[agecol]),
                    'preservation':s.get('Parameter Value[Carcass Preservation Method]',s.get('Parameter Value[Sample Preservation Method]',''))})
    x.to_csv(B/'prepared'/f'OSD-{sid}_cerebral_mouse_counts.csv.gz')
pd.DataFrame(out).to_csv(B/'tables/terminal_region_mouse_registry.csv',index=False)
(B/'sources/terminal_region_download_provenance.json').write_text(json.dumps(provenance,indent=2))
(B/'tables/terminal_region_animal_independence.json').write_text(json.dumps(overlaps,indent=2))
pd.DataFrame(missing,columns=['OSD','sample','donor','status']).to_csv(B/'tables/terminal_region_missing_bulk_samples.csv',index=False)
print(pd.DataFrame(out).groupby(['OSD','flight','age']).size().to_string(),flush=True)
print(json.dumps(overlaps,indent=2),flush=True)
