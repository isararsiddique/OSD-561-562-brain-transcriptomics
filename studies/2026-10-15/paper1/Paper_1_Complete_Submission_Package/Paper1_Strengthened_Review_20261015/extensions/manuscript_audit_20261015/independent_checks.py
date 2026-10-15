#!/usr/bin/env python3
"""Read-only numerical and figure-index audit, independent of production checks."""
from pathlib import Path
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
checks = []
sources = {}
summary = {}

def read(name):
    path = HERE / name
    sources[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)

def check(name, passed, detail=None):
    checks.append(dict(name=name, passed=bool(passed), detail=detail))

def bh(values):
    """Stable rank/reverse-cumulative-min implementation; no production import."""
    p = np.asarray(values, dtype=float)
    finite = np.isfinite(p)
    q = np.full(len(p), np.nan)
    order = np.argsort(p[finite], kind='stable')
    ranked = p[finite][order]
    corrected = np.minimum(1, np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked)+1))[::-1])[::-1])
    back = np.empty(len(order))
    back[order] = corrected
    q[finite] = back
    return q

def one(frame, **filters):
    for field, value in filters.items():
        frame = frame.loc[frame[field] == value]
    assert len(frame) == 1, filters
    return frame.iloc[0]

lib = read('results/sample_manifest.csv')
animal = read('results/animal_region_manifest.csv')
check('129 libraries and 126 collapsed observations', len(lib)==129 and len(animal)==126
      and animal.n_libraries.sum()==129)
check('83 source animals', lib.animal.nunique()==83 and animal.animal.nunique()==83)
core = animal[animal.accession.isin(['OSD-561', 'OSD-562'])]
paired = core.groupby('animal').accession.nunique()
check('29 core animals and 24 paired core animals', core.animal.nunique()==29
      and (paired==2).sum()==24)
rr10 = animal[animal.accession.isin(['OSD-563','OSD-564'])]
check('10 RR-10 paired animals', (rr10.groupby('animal').accession.nunique()==2).sum()==10)
check('7 location entries from 6 animals', lib.location_conflict.sum()==7
      and lib.loc[lib.location_conflict, 'animal'].nunique()==6)

for filename,pfield,qfield,n in [
    ('results/panel_results.csv','p_animal_permutation','q_BH_all_tests_in_family',232),
    ('results/gradient_results.csv','p_animal_permutation','q_BH_all_tests_in_family',29),
    ('extensions/paired_biology/paired_panel_interactions.csv','p_inference','q_BH_all_panel_tests',104),
    ('extensions/flight_collection_interaction/panel_flight_collection_interactions.csv',
     'p_wild_approximate','q_BH_complete_panel_family',32),
    ('results/regional_gene_replication.csv.gz','replication_max_p','q_BH_global_replication_family',54989),
    ('extensions/flight_collection_interaction/gene_flight_collection_interactions_exploratory.csv.gz',
     'pvalue','q_BH_complete_gene_interaction_family',72023),
]:
    data = read(filename)
    check(filename+' finite family size', data[pfield].notna().sum()==n)
    check(filename+' independent BH', np.allclose(bh(data[pfield]), data[qfield],
          equal_nan=True, rtol=1e-10, atol=1e-12))

ip = read('extensions/flight_collection_interaction/panel_flight_collection_interactions.csv')
c = one(ip, accession='OSD-562', sensitivity='all_animals', panel='Circadian')
ce = one(ip, accession='OSD-562', sensitivity='exclude_any_region_metadata_conflict', panel='Circadian')
summary['hippocampal_Circadian'] = {k:float(c[k]) for k in [
    'interaction_log2','bootstrap_lower95','bootstrap_upper95','p_wild_approximate',
    'q_BH_complete_panel_family','p_Monte_Carlo_lower95','p_Monte_Carlo_upper95',
    'mean_NB_interaction_log2']}
summary['hippocampal_Circadian_exclusions'] = {k:float(ce[k]) for k in [
    'interaction_log2','q_BH_complete_panel_family']}
check('Circadian direct difference-in-differences', np.isclose(c.interaction_log2,
      c.FLT_minus_GC_ISST_log2-c.FLT_minus_GC_LAR_log2, atol=1e-12))
check('One interaction panel passes the full family', (ip.q_BH_complete_panel_family<=.05).sum()==1)
check('Wild resolution and 499999 draws', (ip.n_wild_draws==499999).all()
      and np.allclose(ip.p_wild_approximate*500000,
                      np.round(ip.p_wild_approximate*500000), atol=1e-8))
for fname,field in [('leave_one_animal_out_panel_interactions.csv','interaction_without_animal'),
                    ('leave_one_gene_out_panel_interactions.csv','interaction_without_gene_log2')]:
    data = read('extensions/flight_collection_interaction/'+fname)
    sel = data[(data.accession=='OSD-562')&(data.sensitivity=='all_animals')&(data.panel=='Circadian')]
    check('Every full-data Circadian '+fname+' effect negative', (sel[field]<0).all())
    summary[fname] = {'minimum':float(sel[field].min()),'maximum':float(sel[field].max()),'n':len(sel)}

candidate = read('extensions/flight_collection_interaction/all_primary_NB_FDR_candidates_with_sensitivity.csv')
check('All primary discoveries retained, 109/7', candidate.groupby('accession').size().to_dict()
      == {'OSD-561':109,'OSD-562':7})
check('90/2 original-global discoveries persist after exclusions',
      candidate[candidate.passes_original_global_FDR_in_both].groupby('accession').size().to_dict()
      == {'OSD-561':90,'OSD-562':2})
summary['cerebellar_globins'] = candidate[(candidate.accession=='OSD-561')
    & candidate.gene_symbol.isin(['Hbb-bt','Hbb-bs'])][
    ['gene_symbol','log2FoldChange','q_BH_complete_gene_interaction_family',
     'q_BH_complete_gene_interaction_family_metadata_exclusion']].to_dict('records')

hs = read('extensions/ground_holdout/holdout_summary.csv')
h = read('extensions/ground_holdout/holdout_split_results.csv')
roles = read('extensions/ground_holdout/split_animal_roles.csv')
check('362 unique overlapping splits, 30 unique source animals',
      len(h[['accession','split_id']].drop_duplicates())==362 and roles.animal.nunique()==30)
check('Core sensor AUC 1 in every split',
      (h[h.accession.isin(['OSD-561','OSD-562'])].sensor_age_stratified_AUC==1).all())
summary['holdout_score'] = []
for accession in ['OSD-561','OSD-562','OSD-525']:
    row = one(hs, accession=accession, method='score')
    actual = h[(h.accession==accession)&(h.method=='score')]
    fields = ['sensor_age_stratified_AUC','whole_genome_age_stratified_AUC',
              'train_test_contrast_cosine','outside_prediction_MSE_ratio_to_nuisance_only',
              'heldout_outside_effect_energy_removed_pct']
    check(accession+' holdout independent means', all(np.isclose(actual[f].mean(),row[f+'_mean'],
          atol=1e-12) for f in fields))
    check(accession+' per-split MSE ratio', np.allclose(actual.outside_prediction_MSE_ratio_to_nuisance_only,
          actual.outside_prediction_MSE/actual.outside_nuisance_only_prediction_MSE))
    summary['holdout_score'].append(dict(accession=accession,
        n_splits=len(actual),**{f:float(row[f+'_mean']) for f in fields},
        prediction_change_pct=100*(float(row.outside_prediction_MSE_ratio_to_nuisance_only_mean)-1)))

geometry = read('results/operator_geometry.csv')
summary['real_score_gain'] = [float(geometry.score_spectral_gain.min()),float(geometry.score_spectral_gain.max())]
check('15 real geometries; gain equals loading norm and bounded gain <= 2',len(geometry)==15
      and np.allclose(geometry.score_spectral_gain,geometry.gamma_l2_norm)
      and (geometry.bounded_score_spectral_gain<=2+1e-10).all())
joint = read('results/joint_bootstrap_summary.csv')
summary['score_outside_joint_intervals'] = joint[(joint['mode']=='score')
    & (joint.metric=='outside_sensor_energy_removed_pct')][
    ['accession','point_estimate','bootstrap_lower95','bootstrap_upper95']].to_dict('records')
check('All original-score joint outside intervals cross zero', all(
    x['bootstrap_lower95']<0<x['bootstrap_upper95'] for x in summary['score_outside_joint_intervals'])
    and len(summary['score_outside_joint_intervals'])==3)

lookups = read('results/manuscript_numeric_lookups.csv')
for i,row in lookups.iterrows():
    source = row.source_file if row.source_file.startswith('extensions/') else 'results/'+row.source_file
    selected = one(read(source), **json.loads(row.filters))
    check('Recorded manuscript lookup '+str(i), np.isclose(float(selected[row.field]),row.value,
          rtol=1e-10,atol=1e-12), source+' / '+row.field)

manifest_path = HERE/'results/publication_figure_manifest.json'
figs = json.loads(manifest_path.read_text())
check('Final figure order is 1-10 then S1-S4', [re.match(r'figure(S?\d+)',x['filename'])[1]
    for x in figs]==[str(x) for x in range(1,11)]+['S'+str(x) for x in range(1,5)])
panel_inventory = []
for f in figs:
    path = HERE/'publication_figures'/f['filename']
    image = Image.open(path.with_suffix('.png'))
    texts = [''.join(e.itertext()).strip() for e in ET.parse(path.with_suffix('.svg')).getroot().iter()
             if e.tag.rsplit('}',1)[-1]=='text']
    labels = [re.match(r'^([a-j])\s{2}',t)[1] for t in texts if re.match(r'^([a-j])\s{2}',t)]
    panel_inventory.append(dict(figure=f['filename'],labels=labels,caption=f['caption']))
    check(f['filename']+' PNG >=599dpi and editable Arial SVG',
          min(image.info.get('dpi',[0,0]))>=599 and any('Arial' in e.attrib.get('style','')
          for e in ET.parse(path.with_suffix('.svg')).getroot().iter()
          if e.tag.rsplit('}',1)[-1]=='text'))
    number = re.match(r'figure(S?\d+)',f['filename'])[1]
    expected = ('Supplementary Figure '+number[1:]+'.') if number.startswith('S') else 'Figure '+number+'.'
    check(f['filename']+' caption number',f['caption'].startswith(expected))
    check(f['filename']+' caption <=350 words',len(f['caption'].split())<=350)
summary['figure_panel_inventory'] = panel_inventory
summary['obsolete_publication_PNGs'] = sorted(f.stem for f in (HERE/'publication_figures').glob('*.png')
    if f.stem not in {x['filename'] for x in figs} and f.stem!='GRAPHICAL_ABSTRACT')

report = dict(passed=all(x['passed'] for x in checks),n_checks=len(checks),
              failed=[x for x in checks if not x['passed']],summary=summary,
              checks=checks,input_snapshot_sha256=sources,
              scope='Independent read-only saved-result and figure-index audit; no raw-read refitting or render QA.')
(OUT/'independent_checks.json').write_text(json.dumps(report,indent=2,default=str)+'\n')
print(json.dumps({k:report[k] for k in ['passed','n_checks','failed']},indent=2))
