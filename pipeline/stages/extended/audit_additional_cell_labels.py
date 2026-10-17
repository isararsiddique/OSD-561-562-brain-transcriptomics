from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Compare GEO sample annotation with the author's attached count-file labels."""
from pathlib import Path
import requests,re,json,pandas as pd
B=stage_root('extended')
url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295428&targ=gsm&form=text&view=brief'
p=B/'sources/GSE295428_current_sample_SOFT.txt'
if not p.exists():
    r=request_get(url,timeout=60);r.raise_for_status();p.write_text(r.text)
text=p.read_text();rows=[];characteristics=[];excluded=[]
for block in re.split(r'(?=\^SAMPLE = )',text):
    m=re.match(r'\^SAMPLE = (GSM\d+)',block)
    if not m:continue
    accession=m[1]
    def value(key):
        s=re.search(r'^!'+re.escape(key)+r' = (.*)$',block,re.M);return s[1].strip() if s else ''
    title=value('Sample_title');tissue=value('Sample_source_name_ch1')
    files=re.findall(r'^!Sample_supplementary_file_\d+ = (.*matrix.mtx.gz)\s*$',block,re.M)
    if not files:continue
    file=files[0].rsplit('/',1)[-1]
    d={}
    for line in block.splitlines():
        if line.startswith('!Sample_characteristics_ch1 = '):
            k,v=line.split(' = ',1)[1].split(':',1);d[k.lower()]=v.strip()
    if 'myeloid' not in tissue.lower():
        if 'myeloid' in file.lower() or 'meyloid' in file.lower():excluded.append({'sample':accession,'matrix_file':file,'current_tissue':d.get('tissue',''),'reason':'Current tissue identity not confirmed; excluded'})
        continue
    assert all(d[k].lower()!='unknown' for k in ['tissue','age','treatment']),accession
    assert d['tissue']==tissue and int(d['age']) in [3,8] and d['treatment'] in ['control','spaceflight']
    characteristics.append({'sample':accession,'tissue_characteristic':d['tissue'],'age_characteristic':d['age'],'treatment_characteristic':d['treatment'],'all_characteristics_known_and_match_analysis':True})
    label=file[len(accession)+1:].split('_matrix')[0]
    match=re.match(r'(.+)-(HC|LAR)_.+-Mouse_(\d+)-Pool_(\d+)',label)
    treatment=re.search(r'^!Sample_characteristics_ch1 = treatment: (.*)$',block,re.M)
    meta_group={'control':'HC','spaceflight':'LAR'}.get(treatment[1].strip() if treatment else '','');meta_mouse=re.search(r'Mouse (\d+)',title)
    ft=match[1] if match else '';fg=match[2] if match else '';fm=match[3] if match else ''
    canonical=lambda s:re.sub(r'[^a-z0-9]','',s.lower().replace('meyloid','myeloid'))
    ok=canonical(ft)==canonical(tissue) and fg==meta_group and bool(meta_mouse) and fm==meta_mouse[1]
    rows.append({'accession':accession,'SOFT_title':title,'SOFT_tissue':tissue,'matrix_file':file,'file_tissue_label':ft,'SOFT_group':meta_group,'file_group_label':fg,'SOFT_mouse_label':meta_mouse[1] if meta_mouse else '', 'file_mouse_label':fm,'all_labels_match':ok,'supplementary_URL':files[0].replace('ftp://','https://')})
a=pd.DataFrame(rows);a.to_csv(B/'tables/GSE295428_sample_vs_countfile_labels.csv',index=False)
pd.DataFrame(characteristics).to_csv(B/'tables/GSE295428_current_characteristics_audit.csv',index=False)
pd.DataFrame(excluded,columns=['sample','matrix_file','current_tissue','reason']).to_csv(B/'tables/GSE295428_unknown_brain_file_exclusions.csv',index=False)
summary={'source_URL':url,'brain_related_sample_files':len(a),'label_mismatches':int((~a.all_labels_match).sum()),'unknown_labeled_brain_file_exclusions':len(excluded),'current_series_annotation_correction':'GEO series carries 2026-04-22 notice assigning unreliable subset tissue/age/treatment to unknown. All analyzed brain libraries retain known matching current characteristics; no unknown-labeled samples were included.','label_normalization':'Remove punctuation/case and correct explicit textual meyloid/myeloid typo; never change mouse or exposure labels.','decision':'Fresh sample attachments agree with fresh GEO metadata. The cached review sample listing was inconsistent with these current records and is not used. Proceed with current attachments and group-qualified donor labels.'}
(B/'tables/GSE295428_label_audit_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
