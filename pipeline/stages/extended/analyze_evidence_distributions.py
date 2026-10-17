from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Explain threshold-zero results using complete evidence distributions and designs.

No changed testing families, new significance thresholds, or replacement p-values.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

B=stage_root('extended');V1=stage_root('core');T=B/'tables'
old=pd.read_csv(V1/'tables/bulk_panel_complete.csv');new=pd.read_csv(T/'bulk_method_comparison_panel.csv')
bulk=pd.concat([
    old.assign(method='DESeq2',effect=old.log2FoldChange,q=old.q_panel_primary_full_family),
    new.assign(method=new.method.map({'edgeR_QL_TMM':'edgeR QL','limma_voom_TMM':'limma/voom'}),effect=new.log2FC,q=new.q_full_family)
],ignore_index=True)
def rank_rows(s,labels):
    assert np.isfinite(s.q).all() and s.q.between(np.nextafter(0,1),1).all()
    q=np.sort(s.q.to_numpy(),kind='stable')
    return pd.DataFrame({**labels,'rank':np.arange(1,len(q)+1),'q':q,'negative_log10_q':-np.log10(q)})
pd.concat([rank_rows(s,{'method':m}) for m,s in bulk.groupby('method')]).to_csv(T/'bulk_ranked_evidence.csv',index=False)
summary=[]
for method,s in bulk.groupby('method'):
    summary.append({'method':method,'planned_tests':len(s),'estimable_effects':int(s.effect.notna().sum()),'unavailable_effects':int(s.effect.isna().sum()),'raw_p_below_0_05':int(s.pvalue.lt(.05).sum()),'corrected_q_below_0_05':int(s.q.lt(.05).sum()),'minimum_raw_p':s.pvalue.min(),'minimum_q':s.q.min(),'median_absolute_log2_effect':s.effect.abs().median()})
pd.DataFrame(summary).to_csv(T/'bulk_evidence_distribution_summary.csv',index=False)
wide=bulk.pivot(index=['study','contrast','analyte'],columns='method',values='effect');agreement=[]
for a,b in [('DESeq2','edgeR QL'),('DESeq2','limma/voom'),('edgeR QL','limma/voom')]:
    s=wide[[a,b]].dropna();agreement.append({'method_1':a,'method_2':b,'shared_estimable_contrasts':len(s),'effect_Spearman':spearmanr(s[a],s[b]).statistic,'same_effect_direction_fraction':float((np.sign(s[a])==np.sign(s[b])).mean())})
pd.DataFrame(agreement).to_csv(T/'bulk_effect_agreement_summary.csv',index=False)

ca=pd.read_csv(T/'pathway_CAMERA_complete.csv');rr8=pd.read_csv(T/'sorted_pathway_CAMERA_complete.csv')
families={'core_50_sets_28_contrasts':'Core bulk','independent_bulk_50_sets_3_studies':'Independent bulk','celltypes_50_sets_5_types':'RR3 cell types'}
ca=ca.assign(branch=ca.family.map(families),unit=ca.study,q=ca.q_full_family)
rr8=rr8.assign(branch='RR8 enriched populations',unit=rr8.tissue,q=rr8.q_full_1200_family)
pathways=pd.concat([ca,rr8],ignore_index=True)
pd.concat([rank_rows(s,{'branch':branch,'correlation_model':model}) for (branch,model),s in pathways.groupby(['branch','correlation_model'])]).to_csv(T/'pathway_ranked_evidence.csv',index=False)
rows=[]
for (branch,model),s in pathways.groupby(['branch','correlation_model']):
    rows.append({'branch':branch,'correlation_model':model,'planned_tests':len(s),'covered_at_least_10_genes':int(s.NGenes.ge(10).sum()),'raw_p_below_0_05':int(s.PValue.lt(.05).sum()),'corrected_q_below_0_05':int(s.q.lt(.05).sum()),'minimum_q':s.q.min(),'median_q':s.q.median()})
pd.DataFrame(rows).to_csv(T/'pathway_evidence_distribution_summary.csv',index=False)
# Show one index exposure contrast per study/population. Reported weighted
# CAMERA correlations can vary across contrasts; do not assume invariance or
# count repeated contrasts as independent correlation observations.
estimated=pathways[pathways.correlation_model.eq('estimated')]
correlation=[]
for (branch,unit,gene_set),s in estimated.groupby(['branch','unit','set']):
    values=s.Correlation.dropna().to_numpy()
    assert s.NGenes.nunique()==1
    contrast='flight_ISS_T_age12' if branch=='Core bulk' else 'flight_balanced_age' if branch=='RR8 enriched populations' else 'flight_vs_ground'
    selected=s[s.contrast.eq(contrast)];assert len(selected)==1
    r=selected.iloc[0];correlation.append({'branch':branch,'unit':unit,'set':gene_set,'index_contrast':contrast,'NGenes':int(r.NGenes),'estimated_correlation':r.Correlation,'minimum_across_contrasts':values.min() if len(values) else np.nan,'maximum_across_contrasts':values.max() if len(values) else np.nan,'displayed_adequate_coverage':bool(r.NGenes>=10 and np.isfinite(r.Correlation)),'variance_inflation_factor':max(1,1+(r.NGenes-1)*r.Correlation) if np.isfinite(r.Correlation) else np.nan})
pd.DataFrame(correlation).to_csv(T/'pathway_unique_correlation_estimates.csv',index=False)
score=pd.read_csv(T/'pathway_singscore_values.csv')
example=score[score['set'].eq('HALLMARK_INTERFERON_ALPHA_RESPONSE') & score.study.isin(['OSD-352','OSD-563'])].copy()
assert len(example)==16 and not example.duplicated(['study','sample']).any()
example.to_csv(T/'figure6_interferon_donor_scores.csv',index=False)
exact=pd.read_csv(T/'pathway_animal_exact_tests.csv');exact=exact[exact['set'].eq('HALLMARK_INTERFERON_ALPHA_RESPONSE') & exact.study.isin(['OSD-352','OSD-563'])]
for r in exact.itertuples():
    s=example[example.study.eq(r.study)];effect=s[s.flight.eq('FLT')].positive_set_singscore.mean()-s[s.flight.eq('GC')].positive_set_singscore.mean()
    assert np.isclose(effect,r.singscore_effect)
exact.to_csv(T/'figure6_interferon_animal_tests.csv',index=False)
df=pd.read_csv(T/'pathway_camera_degrees_of_freedom.csv')
example_assumptions=pathways[pathways['set'].eq('HALLMARK_INTERFERON_ALPHA_RESPONSE') & pathways.unit.isin(['OSD-352','OSD-563'])].merge(df,on='unit',validate='many_to_one')
example_assumptions['variance_inflation_factor']=np.maximum(1,1+(example_assumptions.NGenes-1)*example_assumptions.Correlation)
example_assumptions['used_camera_df']=np.where(example_assumptions.correlation_model.eq('estimated'),example_assumptions.estimated_correlation_df,example_assumptions.fixed_correlation_df)
example_assumptions.to_csv(T/'interferon_CAMERA_assumption_diagnostic.csv',index=False)

report='''# Why some discovery counts are zero

Zero is a thresholded discovery count, not zero RNA expression, zero estimated effect, or a missing measurement. The underlying effects and evidence remain informative. Figures 2 and 6 now emphasize continuous evidence and individual-animal measurements; the original corrected-discovery counts remain in the tables and captions.

## Bulk methods

All three methods estimate 3,584 of the 4,228 planned target contrasts. The remaining 644 unavailable contrasts retain p=1 in the complete family. Limma/voom has 285 unadjusted p<0.05 contrasts, but its smallest full-family q is 0.0799437. Thus it has no q<0.05 discovery. Unadjusted results are exploratory and are not newly promoted to discoveries.

DESeq2 and voom effects agree in direction for 90.4% of the 3,584 shared estimable contrasts (Spearman 0.894). For young terminal cerebellum, PYCARD is about -1.67 log2 fold change in DESeq2 and -1.59 in voom; the adjusted evidence differs. Different count assumptions, mean-variance weighting and uncertainty estimation affect p-values even when the effects agree. This is method sensitivity, not a missing-data plotting error. The strongest DESeq2 results across the entire family include cohort contrasts; they should not all be presented as flight effects.

The revised full-family rank curves plot every adjusted p-value in order, including unavailable hypotheses at q=1. The dashed q=0.05 threshold is unchanged. Curve height shows evidence strength; it does not show expression magnitude. No significance threshold or correction family was loosened.

## Pathway uncertainty

CAMERA accounts for co-regulation among genes in a set. For a set with m genes and a fixed within-set correlation r, its variance inflation factor is 1+(m-1)r. Increasing fixed r from 0.01 to 0.10 increases uncertainty while preserving the same fixed-correlation degrees-of-freedom rule. For the 86-gene RR3 interferon-alpha set, this factor rises from 1.85 to 9.50.

Estimating correlation changes two assumptions in the installed limma 3.66.0 implementation: it estimates the within-set correlation and uses min(residual sample degrees of freedom, number of genes minus two) for the CAMERA test. Fixed-correlation CAMERA instead uses number of genes minus two here. The actual estimated-correlation designs have only 3-19 residual degrees of freedom. The absence of corrected results under estimated correlation therefore cannot be attributed only to a larger correlation coefficient. Raw estimated correlations can be negative; with the configured default, CAMERA floors the effective variance inflation at one.

The diagnostic saves the installed function source and the actual sample/design dimensions. Revised Figure 6 shows one correlation estimate per study-pathway for the index exposure contrast: young terminal flight for core bulk, age-balanced flight for RR8, and flight versus ground elsewhere. Weighted CAMERA estimates can differ across contrasts; their full range is retained in the correlation table. The figure also shows complete ranked evidence under all three specifications and individual-mouse interferon-alpha singscores. Correlation distributions omit sets with inadequate coverage; these sets remain in their original testing families with p=1. Pathways overlap and their distributions are descriptive, not independent replicates.

## Animal-level biological measurements

Interferon-alpha singscore differences are negative in RR3 (-0.01334; 3 flight and 3 ground mice) and RR10 cerebellum (-0.00895; 5 flight and 5 ground mice). Exact animal-label p-values are 0.300 and 0.302; both have full-family q=0.8535. The direction is visible in individual-animal measurements, but this does not confirm a pathway flight effect. Singscores use within-sample gene ranks; absolute scores across missions and their different gene backgrounds are not a harmonized effect scale.

These analyses clarify the uncertainty rather than manufacture positive results. The focused RR8 AIF1/PYCARD/AIM2 transcript candidate remains subject to its previously documented processing-day and cell-mixture limitations.

## Reproducible outputs

- `audit_camera_assumptions.R` and `logs/CAMERA_implementation_audit.txt`: installed method and actual degrees-of-freedom audit.
- `analyze_evidence_distributions.py`: complete ranked evidence, effect agreement, unique correlations and animal-level score checks.
- `tables/bulk_evidence_distribution_summary.csv` and `tables/pathway_evidence_distribution_summary.csv`: all threshold and continuous summaries.
- `tables/bulk_ranked_evidence.csv` and `tables/pathway_ranked_evidence.csv`: every plotted evidence point.
- `tables/pathway_unique_correlation_estimates.csv`, `tables/figure6_interferon_donor_scores.csv` and `tables/interferon_CAMERA_assumption_diagnostic.csv`: explanatory measurements and assumptions.
'''
(B/'WHY_ZERO_RESULTS.md').write_text(report)
print('Analyzed complete evidence distributions, method agreement, CAMERA assumptions and 16 animal scores')
