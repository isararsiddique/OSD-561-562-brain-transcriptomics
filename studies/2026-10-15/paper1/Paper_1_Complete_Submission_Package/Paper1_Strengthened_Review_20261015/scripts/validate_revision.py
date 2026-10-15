#!/usr/bin/env python3
"""Check bundled numerical families, provenance and rendered-document consistency.

This validates a review package; it never certifies causality or acceptance.
--submission-ready additionally requires the recorded author/archive release steps.
"""
import argparse
import hashlib
import json
import math
import re
import sys
import unicodedata
from pathlib import Path
import numpy as np
import pandas as pd
from pypdf import PdfReader
from PIL import Image
from zipfile import ZipFile
import xml.etree.ElementTree as ET
from statistics_core import bh
from run_revision import DATASETS, check_input_integrity
HERE=Path(__file__).resolve().parents[1];R=HERE/'results'
def norm(text):
    return re.sub('[^a-z0-9]','',unicodedata.normalize('NFKD',str(text)).encode('ascii','ignore').decode().lower())
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--submission-ready',action='store_true');args=ap.parse_args()
    checks=[];errors=[]
    def ck(name,ok):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok: errors.append(name)
    def csv(name):return pd.read_csv(HERE/name)
    def js(name):return json.loads((HERE/name).read_text())
    def close(a,b):return np.allclose(a,b,rtol=1e-9,atol=1e-11,equal_nan=True)
    try:check_input_integrity();ck('Input SHA-256 provenance',True)
    except Exception as e:ck('Input SHA-256 provenance: '+str(e),False)
    lib=csv('results/sample_manifest.csv');animals=csv('results/animal_region_manifest.csv')
    ck('129 unique deposited libraries',len(lib)==129 and not lib.duplicated(['accession','sample']).any())
    ck('126 animal-region columns; 83 animals; three cohorts',len(animals)==126 and lib.animal.nunique()==83 and lib.mission.nunique()==3 and animals.n_libraries.sum()==129)
    core=lib[lib.accession.isin(['OSD-561','OSD-562'])]
    ck('53 core libraries from 29 source animals',len(core)==53 and core.animal.nunique()==29)
    ck('Seven location conflicts from six animals',lib.location_conflict.sum()==7 and lib[lib.location_conflict].animal.nunique()==6)
    ident=[]
    for n,(mission,_) in DATASETS.items():
        md=csv(f'inputs/OSD-{n}_metadata.csv.gz')
        ident.extend([{'mission':mission,'source':s,'rfid':r} for s,r in md[['Source Name','Comment[RFID]']].itertuples(index=False,name=None)])
        counts=csv(f'inputs/OSD-{n}_counts.csv.gz').iloc[:,1:].to_numpy()
        ck(f'{n}: raw collapsed-input libraries are nonnegative integer counts',np.isfinite(counts).all() and (counts>=0).all() and np.array_equal(counts,np.round(counts)))
    ident=pd.DataFrame(ident)
    ck('SourceName and RFID agree across regions',(ident.groupby(['mission','rfid']).source.nunique()<=1).all() and (ident.groupby('source').rfid.nunique()<=1).all())
    for rec in js('inputs/upstream_count_provenance.json'):
        n=rec['accession'].split('-')[-1];a=pd.read_csv(HERE/'inputs'/rec['bundled_upstream_path'],index_col=0);b=pd.read_csv(HERE/f'inputs/OSD-{n}_counts.csv.gz',index_col=0)
        ck(f'{n}: every count matches independent NASA STAR download',set(a.index)==set(b.index) and set(a.columns)==set(b.columns) and np.array_equal(a.loc[b.index,b.columns],b))
    for f,expected in [('panel_results.csv',232),('gradient_results.csv',29)]:
        d=csv('results/'+f)
        if 'n_assignments_evaluated' not in d:
            resolution=csv('results/panel_results.csv')[['contrast','sensitivity','n_assignments_evaluated']].drop_duplicates()
            d=d.merge(resolution,on=['contrast','sensitivity'],validate='one_to_one')
        ck(f'{f}: complete family size {expected}',len(d)==expected)
        ck(f'{f}: independent BH recomputation',close(d.q_BH_all_tests_in_family,bh(d.p_animal_permutation)))
        ck(f'{f}: positive exact-enumeration probabilities',(d.p_animal_permutation>0).all() and (d.p_animal_permutation<=1).all() and close(d.p_animal_permutation*d.n_assignments_evaluated,np.round(d.p_animal_permutation*d.n_assignments_evaluated)))
    p=csv('results/panel_results.csv');ck('All 29 primary/sensitivity contrasts present',p[['contrast','sensitivity']].drop_duplicates().shape[0]==29)
    parts=csv('results/OSD525_disjoint_anchor_partitions.csv')
    ck('126 RR23 partitions; 378 operator rows',len(parts)==378 and parts.split.nunique()==126)
    ck('RR23 partition ground roles are disjoint',all(not(set(a.split('|'))&set(b.split('|'))) for a,b in zip(parts.training_HGC_animals,parts.evaluation_HGC_animals)))
    a=csv('results/ground_anchor_sensitivity.csv');apanel=csv('results/ground_anchor_panel_sensitivity.csv');geom=csv('results/operator_geometry.csv')
    ck('15 anchor applications, four modes and eight panels each',len(a)==60 and len(apanel)==480 and len(geom)==15)
    ck('Real score gain identity and bound',close(geom.score_spectral_gain,geom.gamma_l2_norm) and (geom.bounded_score_spectral_gain<=2+1e-10).all() and close(geom.w_dot_gamma,1))
    fits=js('results/count_model_fit_diagnostics.json');nb=csv('results/count_model_panel_comparison.csv')
    ck('Seven integer-count fits and 15 contrasts',len(fits)==7 and sum(x['n_contrasts']for x in fits)==15 and all(x['pydeseq2_version']=='0.5.4'for x in fits))
    ck('All 120 estimator comparisons and 111 direction agreements',len(nb)==120 and int(nb.same_effect_direction.sum())==111 and close(nb.same_effect_direction,np.sign(nb.mean_log2_effect)==np.sign(nb.mean_NB_log2_effect)))
    rep=csv('results/regional_gene_replication.csv.gz')
    ck('Complete 54989 gene IUT family and BH/BY recomputation',len(rep)==54989 and close(rep.q_BH_global_replication_family,bh(rep.replication_max_p)) and close(rep.q_BY_global_replication_family,np.minimum(1,bh(rep.replication_max_p)*np.sum(1/np.arange(1,len(rep)+1)))))
    ck('No cross-cohort gene survives either correction',(rep.q_BH_global_replication_family>.05).all() and (rep.q_BY_global_replication_family>.05).all())
    paired=csv('extensions/paired_biology/paired_panel_interactions.csv')
    ck('Complete 104 paired tests and independent BH',len(paired)==104 and close(paired.q_BH_all_panel_tests,bh(paired.p_inference)) and (paired.q_BH_all_panel_tests>.05).all())
    ip=csv('extensions/flight_collection_interaction/panel_flight_collection_interactions.csv');ig=csv('extensions/flight_collection_interaction/gene_flight_collection_interactions_exploratory.csv.gz')
    ck('Complete 32 interaction panels and independent BH',len(ip)==32 and close(ip.q_BH_complete_panel_family,bh(ip.p_wild_approximate)))
    ck('499999 wild draws and their probability resolution',(ip.n_wild_draws==499999).all() and close(ip.p_wild_approximate*500000,np.round(ip.p_wild_approximate*500000)))
    finite=ig.pvalue.notna();q=np.full(len(ig),np.nan);q[finite]=bh(ig.loc[finite,'pvalue'])
    ck('72023 evaluable interaction gene tests; complete BH',finite.sum()==72023 and close(ig.q_BH_complete_gene_interaction_family,q))
    c=ip[(ip.accession=='OSD-562')&(ip.sensitivity=='all_animals')&(ip.panel=='Circadian')].iloc[0]
    ck('Highlighted direct Circadian interaction matches model table',close(c.interaction_log2,c.FLT_minus_GC_ISST_log2-c.FLT_minus_GC_LAR_log2) and close(c.q_BH_complete_panel_family,.026816))
    candidate=csv('extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv')
    ck('Primary gene discoveries and persistence retained',candidate.groupby('accession').size().to_dict()=={'OSD-561':109,'OSD-562':7} and candidate[candidate.passes_original_global_FDR_in_both].groupby('accession').size().to_dict()=={'OSD-561':90,'OSD-562':2})
    iv=js('extensions/flight_collection_interaction/numeric_model_validation.json');ck('Independent four-group/count/sandwich audit',len(iv['independent_real_data_design_checks'])==4 and all(x['independent_four_group_contrast_max_error']<1e-10 and x['HC3_direct_sandwich_checked_all_8_panels']for x in iv['independent_real_data_design_checks']))
    h=csv('extensions/ground_holdout/holdout_split_results.csv');hs=csv('extensions/ground_holdout/holdout_summary.csv');roles=csv('extensions/ground_holdout/split_animal_roles.csv');cv=csv('extensions/ground_holdout/nested_training_cv.csv')
    ck('362 overlapping ground splits, five methods',len(h)==1810 and h[['accession','split_id']].drop_duplicates().shape[0]==362 and roles.animal.nunique()==30 and (h.n_shared_animals==0).all())
    ck('Ground holdout role allocations disjoint',not roles.duplicated(['accession','split_id','animal']).any())
    ck('Core holdout sensor AUC perfect on all 162 splits',(h[h.accession.isin(['OSD-561','OSD-562'])].sensor_age_stratified_AUC==1).all())
    ck('Nested CV chooses only training argmin',all(x.selected.sum()==1 and close(x.loc[x.selected,'inner_error_ratio'].iloc[0],x.inner_error_ratio.min())for _,x in cv.groupby(['accession','split_id'])))
    ck('Outside prediction ratios independently recomputed',close(h.outside_prediction_MSE_ratio_to_nuisance_only,h.outside_prediction_MSE/h.outside_nuisance_only_prediction_MSE))
    for row in hs.itertuples():
        sel=h[(h.accession==row.accession)&(h.method==row.method)]
        for field in ['sensor_age_stratified_AUC','outside_prediction_MSE_ratio_to_nuisance_only','heldout_outside_effect_energy_removed_pct']:
            ck(f'Holdout summary {row.accession}/{row.method}/{field}',close(getattr(row,field+'_mean'),sel[field].mean()))
    hv=js('extensions/ground_holdout/verification_report.json');ck('1533 independent ground checks incl raw-count reconstructions',hv['passed'] and hv['n_checks']==1533 and hv['n_failed']==0 and len(hv['independently_recomputed_cases'])==6)
    exact=csv('extensions/independent_audit_20261015/exact_ground_sensor_association.csv')
    ck('Exact audit retains all 18 panel/comparator by region tests',len(exact)==18 and exact.groupby('accession').size().to_dict()=={'OSD-561':9,'OSD-562':9} and not exact.duplicated(['accession','panel_or_comparator']).any())
    ep=exact.p_exact_complete_label_enumeration.to_numpy()
    ck('Exact audit independent complete-family BH and BY',close(exact.q_BH_complete_audit_18_tests,bh(ep)) and close(exact.q_BY_complete_audit_18_tests,np.minimum(1,bh(ep)*np.sum(1/np.arange(1,19)))))
    ck('Exact label probabilities reconstructed from all assignments',close(ep,exact.n_null_at_least_observed/exact.n_label_assignments) and (exact.n_label_assignments==400).all() and (exact.n_overlapping_holdouts_per_assignment==81).all() and (exact.n_ground_animals==12).all())
    for accession in ['OSD-561','OSD-562']:
        null=csv(f'extensions/independent_audit_20261015/{accession}_exact_label_null.csv')
        ck(accession+': complete 400-assignment null table',len(null)==400)
        for row in exact[exact.accession==accession].itertuples():
            values=null[row.panel_or_comparator].to_numpy()
            ck(accession+'/'+row.panel_or_comparator+': null means and exceedances independently recomputed',close(values.mean(),row.null_AUC_mean) and close(values.mean(),.5) and int(np.sum(values>=row.mean_crossfitted_age_stratified_AUC-1e-12))==row.n_null_at_least_observed)
    sat=csv('extensions/independent_audit_20261015/saturated_age_panel_interactions.csv')
    ck('Saturated-age complete family remains nonsignificant',len(sat)==64 and close(sat.q_BH_all_64_saturated_estimands,bh(sat.p_HC3_normal_approximate)) and (sat.q_BH_all_64_saturated_estimands>.05).all())
    newverify=js('extensions/independent_audit_20261015/exact_ground_verification_summary.json')
    ck('Independent explicit OLS refits verify cached null algebra',newverify['status']=='PASS' and newverify['n_explicit_refitted_splits']==486 and newverify['maximum_error']<1e-12)
    qc=csv('extensions/real_data_structure/animal_count_qc.csv');qv=js('extensions/real_data_structure/structure_validation.json');ck('All 126 observations in count QC; independent PCA check',len(qc)==126 and qc.accession.nunique()==7 and qv['PCA_eigenspectrum_and_reconstruction_independently_checked'])
    for source,keys,expected,nrows in [('results/known_truth_benchmark_repetitions.csv',['scenario','repeat'],1000,7000),('extensions/mathematics_review/high_gain_benchmark_repetitions.csv',['true_gain','n_train','n_target','coupling','repetition'],1600,6400)]:
        d=csv(source);ck(f'{source}: independent dataset and method accounting',len(d)==nrows and len(d[keys].drop_duplicates())==expected)
        gain='spectral_gain'if 'spectral_gain'in d else 'operator_gain';ck(source+': every bounded gain within cap',(d[d.method=='bounded_score'][gain]<=2+1e-10).all())
    jm=csv('results/joint_bootstrap_summary.csv');ck('Every joint summary uses 500 full refits',(jm.repetitions==500).all() and jm.normalization_refitted.all() and jm.ground_anchor_refitted.all() and jm.target_refitted.all())
    mv=js('extensions/mathematics_review/independent_operator_test_results.json');ck('100 dense operator tests pass',mv['all_numerical_identities_pass'] and mv['random_geometry_trials']==100)
    for row in csv('results/manuscript_numeric_lookups.csv').itertuples():
        path=HERE/row.source_file if row.source_file.startswith('extensions/')else R/row.source_file;f=pd.read_csv(path)
        for k,v in json.loads(row.filters).items():f=f[f[k]==v]
        ck('Manuscript numeric provenance '+row.source_file+'/'+row.field+'/'+row.filters,len(f)==1 and math.isclose(float(f.iloc[0][row.field]),row.value,rel_tol=1e-10,abs_tol=1e-12))
    from validate_strengthening import validate as validate_strengthening
    for entry in validate_strengthening(HERE):ck('Strengthening: '+entry['name'],entry['passed'])
    figure_manifest=js('results/publication_figure_manifest.json');ck('Thirteen main and seven supplementary figures',len(figure_manifest)==20)
    artwork=js('results/journal_artwork_validation.json')
    ck('All 20 panel figures plus optional overview checked at intended print size',len(artwork)==21 and all(x['colourspace']=='srgb' and x['density']>=600 and x['minimum_print_font_points']>=8 and x['minimum_print_stroke_points']>=1 for x in artwork))
    mainblocks=js('results/manuscript_content.json');suppblocks=js('results/supplement_content.json')
    title=next(x['text'] for x in mainblocks if x['type']=='title')
    cited=[]
    for q in mainblocks:
        if q['type'] not in ['table','figure']:
            for m in re.finditer(r'(?<!Supplementary )\b(?:Fig\.|Figure) (\d+)',q.get('text','')):
                if int(m[1]) not in cited:cited.append(int(m[1]))
    ck('Main figures all cited in numerical order',cited==list(range(1,14)))
    ck('Main equations numbered in display order',[int(re.search(r'\((\d+)\)\.?$',q['text'])[1]) for q in mainblocks if q['type']=='equation']==list(range(1,13)))
    ck('Journal title limit and punctuation',len(title.split())<=15 and not re.search(r'[^A-Za-z0-9\s]',title))
    abstract_idx=next(i for i,x in enumerate(mainblocks) if x.get('text')=='Abstract')
    ck('Journal unstructured abstract at most 150 words',len(mainblocks[abstract_idx+1]['text'].split())<=150)
    ck('All display mathematics and Methods retained in main file',sum(x['type']=='equation' for x in mainblocks)==12 and not any(x['type']=='equation' or (x['type']=='heading' and 'method' in x.get('text','').lower()) for x in suppblocks))
    ck('Funding incorporated in Acknowledgements',not any(x['type']=='heading' and x.get('text')=='Funding' for x in mainblocks) and any(x.get('text')=='Acknowledgements' for x in mainblocks))
    for f in figure_manifest:
        path=HERE/'publication_figures'/f['filename'];im=Image.open(path.with_suffix('.png'));svg=path.with_suffix('.svg').read_text()
        ck('Figure 600dpi and editable vector text '+f['filename'],min(im.info.get('dpi',[0,0]))>=599 and '<text'in svg and 'Arial'in svg)
        ck('Figure caption within 350 words '+f['filename'],len(f['caption'].split())<=350)
    pages={}
    for stem,content,expected in [('MANUSCRIPT_COMPLETE','manuscript_content.json',13),('SUPPLEMENTARY_INFORMATION','supplement_content.json',7)]:
        pdf=PdfReader(HERE/(stem+'.pdf'));pages[stem]=len(pdf.pages);strings=[]
        for i,page in enumerate(pdf.pages):
            text=(page.extract_text()or'').replace('Spaceflight brain transcriptomics','');text=re.sub(r'(?m)^\s*\d+\s*$','',text);strings.append(text)
            ck(stem+f': page {i+1} contains text or figure',bool(norm(text)) or bool(page.images))
        full=norm('\n'.join(strings));blocks=js('results/'+content)
        for i,block in enumerate(blocks):
            if block['type']in ['paragraph','title','heading'] and block.get('text')!='Figures':ck(stem+f': rendered text block {i}',norm(block['text'])in full)
            if block['type']=='figure':
                ck(stem+': full caption '+block['filename'],norm(block['caption'])in full)
                ck(stem+': legend shares figure page '+block['filename'],any(norm(block['caption'])in norm(t) and bool(p.images)for p,t in zip(pdf.pages,strings)))
        z=ZipFile(HERE/(stem+'.docx'));doc=ET.fromstring(z.read('word/document.xml'));styles=ET.fromstring(z.read('word/styles.xml'))
        ns={'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        ck(stem+': embedded figure and equation count',len(doc.findall(".//wp:inline",ns))==expected+sum(x['type']=='equation'for x in blocks))
        ck(stem+': no decorative title borders',not styles.findall(".//w:style[@w:styleId='Title']//w:pBdr",ns))
    author=js('authors.json');ck('Verified author affiliation/correspondence/funding/interests',author['corresponding_email']=='isarar@convolity.com' and len(author['authors'])==5 and bool(author['funding_statement']) and bool(author['competing_interests_statement']))
    manifest=HERE/'PACKAGE_FILE_MANIFEST.json'
    if manifest.exists():
        for rec in json.loads(manifest.read_text()):
            path=HERE/rec['path'];ck('Sealed file '+rec['path'],path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==rec['sha256'])
    remaining=['All five authors must review and approve the final manuscript, contributions and related-manuscript disclosure.','Create a permanent public data/code revision archive and insert its verified identifier.','Complete journal portal reporting-summary/editorial-checklist fields and disclose the overlapping Paper 2.']
    report={'computational_checks':len(checks),'computational_errors':errors,'computational_validation_passed':not errors,'submission_ready':False,'pdf_pages':pages,'release_items_remaining':remaining,'scope':'Internal numerical, provenance and rendered-document consistency; not causal certification or acceptance prediction.','checks':checks}
    (R/'validation_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items()if k!='checks'},indent=2))
    if errors or args.submission_ready:sys.exit(1)
if __name__=='__main__':main()
