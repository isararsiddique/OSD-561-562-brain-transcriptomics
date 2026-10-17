from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Audit revised Figures 2/6 against complete families and recalculated BH q-values."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

B = stage_root('extended')
V1 = stage_root('core')
checks = []
sources = [V1/'tables/bulk_panel_complete.csv', B/'tables/bulk_method_comparison_panel.csv',
           B/'tables/pathway_CAMERA_complete.csv', B/'tables/sorted_pathway_CAMERA_complete.csv']

def check(label, value):
    checks.append({'check': label, 'passes': bool(value)})

def audit_family(label, source, qcol, displayed, planned):
    pcol = 'PValue' if 'PValue' in source else 'pvalue'
    # Unavailable hypotheses retain p=1, as in the documented complete families.
    q = multipletests(source[pcol].fillna(1).to_numpy(), method='fdr_bh')[1]
    check(label+' complete family', len(source) == planned)
    check(label+' finite stored q-values', np.isfinite(source[qcol]).all())
    check(label+' independently recalculated BH', np.allclose(q, source[qcol], rtol=1e-10, atol=1e-12))
    check(label+' documented discovery count', int((q < .05).sum()) == int(displayed))

old = pd.read_csv(sources[0]); new = pd.read_csv(sources[1])
display = pd.read_csv(B/'tables/figure2_method_counts.csv').set_index('method')
for method, source, qcol in [('DESeq2', old, 'q_panel_primary_full_family'),
                            ('edgeR QL', new[new.method.eq('edgeR_QL_TMM')], 'q_full_family'),
                            ('limma/voom', new[new.method.eq('limma_voom_TMM')], 'q_full_family')]:
    audit_family('Figure 2 '+method, source, qcol, display.loc[method, 'corrected_contrasts'], 4228)
for r in pd.read_csv(B/'tables/figure2_displayed_effects.csv').itertuples():
    if r.method == 'DESeq2':
        s = old[(old.study=='OSD-561') & (old.contrast=='flight_ISS_T_age12') & (old.analyte==r.analyte)]
        effect, q = s.iloc[0][['log2FoldChange', 'q_panel_primary_full_family']]
    else:
        method = {'edgeR QL':'edgeR_QL_TMM', 'limma/voom':'limma_voom_TMM'}[r.method]
        s = new[(new.study=='OSD-561') & (new.contrast=='flight_ISS_T_age12') & (new.analyte==r.analyte) & (new.method==method)]
        effect, q = s.iloc[0][['log2FC', 'q_full_family']]
    check('Figure 2 '+r.analyte+' '+r.method+' displayed values', len(s)==1 and np.allclose([effect,q], [r.log2_fold_change,r.q_full_family]))

ca = pd.read_csv(sources[2]); rr8 = pd.read_csv(sources[3])
family = {'Core bulk':'core_50_sets_28_contrasts', 'Independent bulk':'independent_bulk_50_sets_3_studies',
          'RR3 cell types':'celltypes_50_sets_5_types'}
for r in pd.read_csv(B/'tables/figure6_CAMERA_count_grid.csv').itertuples():
    if r.branch == 'RR8 enriched populations':
        s = rr8[rr8.correlation_model.eq(r.correlation_model)]; qcol = 'q_full_1200_family'
    else:
        s = ca[ca.family.eq(family[r.branch]) & ca.correlation_model.eq(r.correlation_model)]; qcol = 'q_full_family'
    audit_family('Figure 6 '+r.branch+' '+r.correlation_model, s, qcol, r.corrected_contrasts, r.planned_family_size)
for r in pd.read_csv(B/'tables/figure6_example_q_values.csv').itertuples():
    s = ca[ca['set'].eq('HALLMARK_INTERFERON_ALPHA_RESPONSE') & ca.study.eq(r.study) & ca.correlation_model.eq(r.correlation_model)]
    check('Figure 6 '+r.study+' '+r.correlation_model+' displayed q', len(s)==1 and np.isclose(s.iloc[0].q_full_family,r.q_full_family) and np.isclose(-np.log10(r.q_full_family),r.negative_log10_q))

ranked_bulk=pd.read_csv(B/'tables/bulk_ranked_evidence.csv')
for method,source,qcol in [('DESeq2',old,'q_panel_primary_full_family'),('edgeR QL',new[new.method.eq('edgeR_QL_TMM')],'q_full_family'),('limma/voom',new[new.method.eq('limma_voom_TMM')],'q_full_family')]:
    s=ranked_bulk[ranked_bulk.method.eq(method)].sort_values('rank');q=np.sort(source[qcol].to_numpy())
    check('Figure 2 '+method+' every ranked evidence point',len(s)==len(q) and np.array_equal(s['rank'],np.arange(1,len(q)+1)) and np.allclose(q,s.q) and np.allclose(-np.log10(q),s.negative_log10_q))
ranked_pathway=pd.read_csv(B/'tables/pathway_ranked_evidence.csv')
for (branch,model),s in ranked_pathway.groupby(['branch','correlation_model']):
    source=rr8[rr8.correlation_model.eq(model)] if branch=='RR8 enriched populations' else ca[ca.family.eq(family[branch]) & ca.correlation_model.eq(model)]
    q=np.sort(source['q_full_1200_family' if branch=='RR8 enriched populations' else 'q_full_family'].to_numpy());s=s.sort_values('rank')
    check('Figure 6 '+branch+' '+model+' every ranked evidence point',len(s)==len(q) and np.array_equal(s['rank'],np.arange(1,len(q)+1)) and np.allclose(q,s.q) and np.allclose(-np.log10(q),s.negative_log10_q))
donors=pd.read_csv(B/'tables/figure6_interferon_donor_scores.csv');tests=pd.read_csv(B/'tables/figure6_interferon_animal_tests.csv')
for r in tests.itertuples():
    s=donors[donors.study.eq(r.study)];effect=s[s.flight.eq('FLT')].positive_set_singscore.mean()-s[s.flight.eq('GC')].positive_set_singscore.mean()
    check('Figure 6 '+r.study+' individual animals and mean effect',len(s)==r.n_flight+r.n_ground and not s['sample'].duplicated().any() and np.isclose(effect,r.singscore_effect))
dof=pd.read_csv(B/'tables/pathway_camera_degrees_of_freedom.csv')
check('CAMERA degrees of freedom use actual sample and design sizes',len(dof)==16 and (dof.estimated_correlation_df==dof.n_samples-dof.design_rank).all() and (dof.fixed_correlation_df==dof.n_genes-2).all())
correlation=pd.read_csv(B/'tables/pathway_unique_correlation_estimates.csv')
check('Correlation plot retains one index contrast per study and set',not correlation.duplicated(['unit','set']).any() and correlation.loc[correlation.displayed_adequate_coverage].NGenes.ge(10).all())
previous=B/'tables/figure_zero_value_audit.json'
if previous.exists():
    for item in json.loads(previous.read_text())['source_files']:
        check('Underlying analysis unchanged: '+resolve_record(item['path']).name,hashlib.sha256(resolve_record(item['path']).read_bytes()).hexdigest()==item['sha256'])
record = {'checks':len(checks), 'passed':sum(c['passes'] for c in checks), 'details':checks,
          'source_files':[{'path':record_path(p), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],
          'interpretation':'Zero counts are verified no-q<0.05 results in complete planned families. Unavailable hypotheses retain p=1; no expression value was imputed for plotting.'}
(B/'tables/figure_zero_value_audit.json').write_text(json.dumps(record, indent=2))
assert record['passed'] == record['checks'], [c for c in checks if not c['passes']]
print(f"Figure value audit: {record['passed']}/{record['checks']} checks passed")
