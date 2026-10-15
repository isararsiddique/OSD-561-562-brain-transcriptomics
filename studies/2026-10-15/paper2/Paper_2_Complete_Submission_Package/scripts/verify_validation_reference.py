"""One-time online freeze of RRRM2 processing reference metadata."""
from pathlib import Path
import json,hashlib,urllib.request,zipfile,re
B=Path(__file__).resolve().parents[1];I=B/'extensions/independent_validation/inputs';O=B/'config/RRRM2_processing_reference';O.mkdir(exist_ok=True)
inv=json.loads((I/'OSD-561_publisher_inventory.json').read_text())['studies']['OSD-561']['study_files'];side=(I/'GLDS-556_rna_seq_processed_md5sum_GLbulkRNAseq.tsv').read_text();evidence=[]
for row in inv:
 if row['file_name'] not in ['GLDS-556_rna_seq_software_versions_GLbulkRNAseq.md','GLDS-556_rna_seq_processing_info_GLbulkRNAseq.zip']:continue
 name=row['file_name'];url='https://osdr.nasa.gov'+row['remote_url'];path=O/name
 if not path.exists():path.write_bytes(urllib.request.urlopen(url,timeout=60).read())
 raw=path.read_bytes();md5=hashlib.md5(raw).hexdigest();matches=[l for l in side.splitlines() if name in l or name.replace('GLDS-556_rna_seq_','') in l]
 assert len(matches)==1,(name,matches);assert md5 in matches[0],(name,md5,matches)
 evidence.append({'file':name,'url':url,'md5':md5,'publisher_md5_matches':True,'sha256':hashlib.sha256(raw).hexdigest()})
 if path.suffix=='.zip':
  with zipfile.ZipFile(path) as z:
   for member in z.namelist():
    if member.endswith('/') or '__MACOSX' in member:continue
    try:t=z.read(member).decode()
    except UnicodeDecodeError:continue
    lines=[l for l in t.splitlines() if any(k in l for k in ['GRCm','genomeDir','sjdbGTFfile','genomeFasta','STAR version'])]
    if lines:print(member,'\n'+'\n'.join(lines)[:3000]);(O/(Path(member).name+'.reference_extract.txt')).write_text('\n'.join(lines)+'\n')
 else:print(raw.decode()[:3800])
(O/'provenance.json').write_text(json.dumps(evidence,indent=2))
