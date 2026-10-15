"""Build complete accession decisions and compact study evidence from frozen metadata."""
from pathlib import Path
import json,csv,hashlib,re,shutil
import pandas as pd
B=Path(__file__).resolve().parents[1];S=B/'sources';T=B/'tables'
local=Path('/Users/zain/Downloads/OSD-561-562-brain-transcriptomics-main/research/brain_integrative_2026/tables/search_catalogue.csv')
if not (S/'local_search_catalogue_20260910.csv').exists():shutil.copy2(local,S/'local_search_catalogue_20260910.csv')
rows={};terms={};search_summary=[]
for term in ['cerebellum','cerebellar','brain']:
 d=json.loads((S/f'OSDR_search_{term}.json').read_text());hits=d['hits']['hits'];n=0
 for h in hits:
  s=h['_source'];a=s.get('Accession','')
  if a.startswith('OSD-'):rows[a]=s;terms.setdefault(a,set()).add(term);n+=1
 search_summary.append(dict(source='NASA OSDR',query=term,reported_hits=d['hits']['total'],returned=len(hits),osd_records=n))
local=pd.read_csv(S/'local_search_catalogue_20260910.csv').fillna('');sel=local[local.organism.str.contains('Mus musculus') & (local.material.str.contains('brain|cerebell|cortex|hippocamp|neural',case=False)|local.title.str.contains('brain|cerebell',case=False))]
for r in sel.itertuples():
 terms.setdefault(r.accession,set()).add('local mouse brain metadata screen')
 if r.accession not in rows:rows[r.accession]={'Study Title':r.title,'Material Type':r.material,'organism':r.organism,'Mission':{'Name':r.mission},'Study Assay Technology Type':r.assay,'Study Description':r.description}
manual={
525:('existing_discovery','Already analyzed RR23; no additional independent cohort.'),
561:('existing_independent_launch_validation','Already analyzed RRRM2; retain its negative frozen-score result.'),
563:('include_same_launch_external_transport','Distinct animals with matched cerebellar bulk sequencing and publisher-verified SJ files. Shared SpaceX-21 launch; female B6129SF2/J; marker injections; euthanasia and dissection dates each aliased with group. Not independent-launch validation.'),
927:('independent_launch_context_deferred','RRRM1 sorted myeloid cerebellum 3-prime 10x scRNA-seq, 7 samples; no deposited processed SJ matrix. Different cell selection/assay and 5-prime candidate architecture unsuitable for direct gene-count validation. Group/dissection date alias. Raw-only junction reprocessing is a separate assay study.'),
930:('independent_launch_context_deferred','RRRM1 sorted non-myeloid cerebellum 3-prime 10x scRNA-seq, 8 samples; no deposited processed SJ matrix. Shares animals with OSD-927; not another independent cohort. Group/dissection date alias.'),
222:('independent_launch_orthogonal_context','STS135 cerebellum microarray, not junction sequencing. Cannot measure the frozen STAR local junction profiles.'),
352:('independent_launch_other_region_context','RR3 cerebral hemisphere sections at hippocampal level, 3+3 bulk animals. Not separately sampled cerebellum; cannot claim matched-tissue replication.'),
562:('same_mission_other_region_context','RRRM2 hippocampus, same mission/source-animal pool as OSD561; not independent replication.'),
613:('same_mission_other_region_context','RRRM2 left cerebral hemisphere, same mission/source-animal pool as OSD561; not independent replication.'),
564:('same_mission_other_region_context','RR10 hippocampus, same mission/source-animal pool as OSD563; not independent replication.'),
612:('same_mission_other_region_context','RR10 left cerebral hemisphere, same mission/source-animal pool as OSD563; not independent replication.'),
33:('exclude_primary_platform_region_replication','STS128 whole-brain microarray; one surviving flight animal, one ground and one vivarium control; no junction sequencing.'),
536:('independent_launch_orthogonal_context','SpaceX12 brain NanoString/neuroinflammation panel described in linked publication; inventory has no STAR SJ data; not matched cerebellar junction sequencing.'),
916:('exclude_primary_mirna_platform','RRRM1 brain miRNA sequencing; cannot measure coding-transcript splice junctions.'),
}
for n in [682,685,698,699]:manual[n]=('independent_launch_other_region_context','SpaceX24 GeoMx regional gene-expression assay with antioxidant BuOE treatment; not cerebellar bulk junction sequencing. Treatment arms and source animals require explicit modeling for any secondary analysis.')
for n in [928,929,931,932]:manual[n]=('same_mission_other_region_context','RRRM1 other-brain-region sorted-cell 3-prime single-cell assay, same source cohort as OSD927/930; not additional independent cerebellar cohort.')
full=[];evidence=[]
for a,s in sorted(rows.items(),key=lambda x:int(x[0][4:])):
 n=int(a[4:]);title=s.get('Study Title','');mat=s.get('Material Type','');org=s.get('organism','');assay=s.get('Study Assay Technology Type','');desc=s.get('Study Description','');mission=s.get('Mission',{}).get('Name','')
 if n in manual:decision,reason=manual[n]
 elif 'Mus musculus' not in org:decision,reason='exclude_primary_species','Not Mus musculus.'
 elif not re.search('brain|cerebell|cortex|hippocamp|neural|microglia',mat,re.I):decision,reason='exclude_primary_tissue','Material is not separately sampled mouse cerebellum.'
 elif not mission:decision,reason='exclude_primary_exposure','Catalog describes ground analog, irradiation or nonflight study; not an actual independent spaceflight cerebellar cohort.'
 else:decision,reason='exclude_primary_region_or_platform','Different brain region and/or protein, imaging, physiological platform; no matched cerebellar junction sequencing identified.'
 meta=S/f'{a}_meta.json';inv=S/f'{a}_files.json';nr='';ns='';sj='';fastq='';version=''
 if meta.exists():
  md=json.loads(meta.read_text())['study'][a];st=md['studies'][0];nr=len(st['materials'].get('sources',[]));ns=len(st['materials'].get('samples',[]));version=md.get('version','');comments={r['name']:r.get('value','') for r in st.get('comments',[])};mission=comments.get('Mission Name',mission)
  evidence.append(dict(accession=a,title=st['title'],description=st['description'],n_source_records=nr,n_sample_records=ns,mission=mission,project=comments.get('Project Identifier',''),sample_collection='\n'.join(x['description'] for x in st['protocols'] if 'sample collection' in x['name'].lower()),library_construction='\n'.join(x['description'] for x in st['protocols'] if 'library construction' in x['name'].lower()),publisher_source_sha256=hashlib.sha256(meta.read_bytes()).hexdigest()))
 if inv.exists():
  fs=json.loads(inv.read_text())['studies'][a]['study_files'];sj=sum(x['file_name'].endswith('SJ.out.tab') for x in fs);fastq=sum(x['file_size'] for x in fs if x['file_name'].endswith('fastq.gz') and 'raw' in x['file_name'])
 full.append(dict(accession=a,title=title,organism=org,material=mat,assay=assay,mission=mission,search_terms=';'.join(sorted(terms[a])),decision=decision,reason=reason,n_sources=nr,n_samples=ns,star_sj_files=sj,raw_fastq_bytes=fastq,live_metadata_version=version,url=f'https://osdr.nasa.gov/bio/repo/data/studies/{a}',description=desc))
pd.DataFrame(full).to_csv(T/'complete_accession_decisions.csv',index=False);pd.DataFrame(evidence).to_csv(T/'candidate_protocol_evidence.csv',index=False);pd.DataFrame(search_summary).to_csv(T/'search_return_counts.csv',index=False)
# Complete title/abstract screen for the two corrected, date-bounded literature queries.
lit={};litqueries={}
for i in [1,2]:
 d=json.loads((S/f'EuropePMC_supplement_{i}.json').read_text())
 assert d['hitCount']==len(d['resultList']['result'])
 for r in d['resultList']['result']:lit[r['id']]=r;litqueries.setdefault(r['id'],[]).append(i)
out=[]
known={'38862479':'Existing OSD352/561/562/563/564/612/613 multiomics cohorts; no new independent cerebellar bulk cohort.','PPR668612':'Preprint version of existing spatial multiomics publication; not another cohort.','35948598':'OSD536; SpaceX12 targeted NanoString brain assay, no splice-junction counts.','37686374':'OSD682/685/698/699; SpaceX24 GeoMx cortex/hippocampus with BuOE, not cerebellum.','22808101':'OSD33; whole-brain microarray, only one flight survivor.','31443374':'OSD202; ground HLU/radiation analog.','18704384':'OSD32; ground hindlimb unloading microarray.','41644562':'RRRM1 miRNA platform, not local mRNA splice junctions.','40136666':'Protein assays, not junction sequencing.'}
for key,r in lit.items():
 title=r.get('title','');abstract=r.get('abstractText','');low=(title+' '+abstract).lower()
 if key in known:reason=known[key]
 elif any(w in title.lower() for w in ['simulat','unloading','irradiat','radiation','hypergravity','centrifug','clinostat']):reason='Ground analog/radiation/altered-gravity study; no new actual-flight cerebellar bulk cohort identified.'
 elif any(w in title.lower() for w in ['gene expression','neurogenes','bdnf','neurotrophic']):reason='Targeted gene expression, review or non-junction assay; abstract does not identify an additional cerebellar bulk RNA-seq cohort.'
 else:reason='Title/abstract provides no additional matched cerebellar mouse bulk RNA-seq cohort; contextual, physiological, protein or review study. Full text not systematically reviewed.'
 out.append(dict(id=key,source=r.get('source'),pmcid=r.get('pmcid',''),doi=r.get('doi',''),title=title,first_publication_date=r.get('firstPublicationDate',''),queries=';'.join(map(str,litqueries[key])),decision='no_new_eligible_bulk_cohort_identified_at_screen',reason=reason,abstract=abstract))
pd.DataFrame(out).to_csv(T/'complete_literature_title_abstract_screen.csv',index=False)
summary={'date':'2026-10-15','live_OSDR_unique_accessions':len(set(a for a in terms if any(t in terms[a] for t in ['brain','cerebellum','cerebellar']))),'total_accessions_including_prior_local_mouse_brain_screen':len(full),'live_detailed_study_metadata_and_inventory_pairs':len(evidence),'new_distinct_animal_same_launch_bulk_cohorts':1,'new_independently_launched_bulk_cerebellar_cohorts_identified':0,'corrected_date_bounded_literature_records_screened':len(out),'initial_literature_query_2_issue':'Unquoted hyphenated mission terms produced11228 hits; onlyfirst1000 returned. Retained as failed overbroad search, not treated as complete. Corrected quoted mission+spaceflight+date query returned2/2.','scope_limit':'Bounded NASA OSDR accession screen and EuropePMC title/abstract discovery, not systematic full-text review or exhaustive search of all repositories.'}
(T/'search_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
