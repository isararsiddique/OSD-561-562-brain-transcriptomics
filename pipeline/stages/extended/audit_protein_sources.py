from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Preserve and audit RR1 protein values; do not guess corrections or pool assays."""
from pathlib import Path
import json,hashlib,zipfile,xml.etree.ElementTree as E
import numpy as np,pandas as pd
B=stage_root('extended');S=B/'sources/OSD-738'
csv=S/'LSDS-107_western-blot_Roy_Brain_WesternBlot_TRANSFORMED.csv'
xlsx=S/'LSDS-107_western-blot_Roy_Brain_WesternBlot_SUBMITTED.xlsx'
n={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(xlsx) as z:
    strings=[]
    if 'xl/sharedStrings.xml' in z.namelist():strings=[''.join(si.itertext()) for si in E.fromstring(z.read('xl/sharedStrings.xml')).findall('x:si',n)]
    rows=[]
    for row in E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//x:sheetData/x:row',n):
        values=[]
        for c in row.findall('x:c',n):
            v=c.find('x:v',n);it=c.find('x:is',n);t=c.get('t')
            values.append(strings[int(v.text)] if t=='s' else ''.join(it.itertext()) if t=='inlineStr' else float(v.text) if v is not None else None)
        rows.append(values)
    w=pd.DataFrame(rows[1:],columns=rows[0])
d=pd.read_csv(csv).merge(w,left_on='Sample Name',right_on='Sample ID',validate='one_to_one')
assert len(d)==20 and np.allclose(d.asc_relative_density,d.ASC,rtol=1e-8)
d['ASC_source_matches_submitted_workbook']=np.isclose(d.asc_relative_density,d.ASC,rtol=1e-8)
d['ASC_orders_of_magnitude_QC_flag']=d.asc_relative_density.abs()>1000*d.asc_relative_density.abs().median()
d[['Sample Name','Group','asc_relative_density','unit.1','ASC','ASC_source_matches_submitted_workbook','ASC_orders_of_magnitude_QC_flag']].to_csv(B/'tables/OSD738_ASC_source_audit.csv',index=False)
summary={'study':'OSD-738','assay':'Western blot relative density normalized to actin; arbitrary_units, not concentration','submitted_sheet':'Brain (other)','n_records':len(d),'n_ASC_QC_flags':int(d.ASC_orders_of_magnitude_QC_flag.sum()),'flagged_sample_IDs':d.loc[d.ASC_orders_of_magnitude_QC_flag,'Sample Name'].tolist(),'ASC_range':[float(d.ASC.min()),float(d.ASC.max())],'same_values_in_submitted_workbook':True,'decision':'Preserve values, flag scale ambiguity and do not guess a repair or use as clean quantitative mRNA replication. Publication remains prior art.','source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [csv,xlsx]}}
(B/'tables/OSD738_protein_source_audit.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
