#!/usr/bin/env python3
"""Scientific panels generated solely from saved, fully retained analysis outputs."""
import json
import shutil
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parents[1]; R=HERE/'results'; F=HERE/'figures'
plt.rcParams.update({'font.family':'Arial','font.size':14,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
C={'raw':'#667788','projection':'#157c99','score':'#ad517a','bounded_score':'#d18a21'}
N={'raw':'Raw','projection':'Projection','score':'Score','bounded_score':'Bounded (L=2)'}
manifest=json.loads((R/'figure_manifest.json').read_text())[:4]
def save(fig,name,caption):
    fig.savefig(F/(name+'.png'),dpi=600,bbox_inches='tight',facecolor='white')
    fig.savefig(F/(name+'.svg'),bbox_inches='tight',facecolor='white'); plt.close(fig)
    manifest.append(dict(filename=name,caption=caption,png_dpi=600,vector_format='SVG'))

a=pd.read_csv(R/'count_model_panel_comparison.csv'); p=pd.read_csv(R/'panel_results.csv')
fig,ax=plt.subplots(1,2,figsize=(11.8,5),gridspec_kw={'width_ratios':[1,1.2]})
for panel,s in a.groupby('panel'):
    ax[0].scatter(s.mean_log2_effect,s.mean_NB_log2_effect,s=25,alpha=.8,label=panel.replace('_',' / '))
lim=np.max(np.abs(a[['mean_log2_effect','mean_NB_log2_effect']].to_numpy()))*1.12
ax[0].plot([-lim,lim],[-lim,lim],color='#aaaaaa',ls='--');ax[0].set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel='Log-expression panel effect',ylabel='Mean negative-binomial log2 effect')
ax[0].set_title('a  Count-model robustness',loc='left',fontweight='bold')
ax[0].text(.04,.95,f'{a.same_effect_direction.sum()}/{len(a)} directions agree\nr = {a.mean_log2_effect.corr(a.mean_NB_log2_effect):.3f}',transform=ax[0].transAxes,va='top')
ax[0].legend(fontsize=11,frameon=False,ncol=2,loc='lower right')
contrast='OSD562_hippocampus_FLT_vs_GC_within_ISS_T'
rows=[]
for panel in ['SRG','Circadian']:
    for sens in ['design_cohort','exclude_all_regions_of_conflicting_animals']:
        s=p[(p.contrast==contrast)&(p.panel==panel)&(p.sensitivity==sens)].iloc[0]; rows.append(s)
for i,s in enumerate(rows):
    color=C['score'] if s.panel=='SRG' else C['projection']
    ax[1].hlines(i,s.bootstrap_lower95,s.bootstrap_upper95,color=color,lw=2)
    ax[1].scatter(s.mean_log2_effect,i,color=color,s=42,zorder=3)
    ax[1].text(.49,i,f'q={s.q_BH_all_tests_in_family:.3f}',va='center',fontsize=12)
    if s.sensitivity=='design_cohort':
        nb=a[(a.contrast==contrast)&(a.panel==s.panel)].iloc[0]
        ax[1].scatter(nb.mean_NB_log2_effect,i+.13,marker='D',facecolors='none',edgecolors=color,s=45)
ax[1].set_yticks(range(4),['SRG: full','SRG: exclusions','Circadian: full','Circadian: exclusions'])
ax[1].invert_yaxis();ax[1].axvline(0,color='#aaa',lw=.8);ax[1].set_xlim(-.8,.78)
ax[1].set_xlabel('Hippocampal FLT − GC (ISS-T), log2')
ax[1].set_title('b  Biological candidate sensitivity',loc='left',fontweight='bold')
fig.tight_layout(w_pad=3)
save(fig,'figure5_count_model_and_biology','Figure 5. Count-model robustness and biological candidate sensitivity. (a) The 120 panel comparisons summarize 15 regional contrasts fitted using PyDESeq2; agreement is analytical robustness on the same observations, not independent replication. (b) Full and metadata-excluded hippocampal ISS-T flight contrasts for SRG and circadian panels. Circles are log-expression effects with 1,000-animal-bootstrap percentile intervals; hollow diamonds show mean unshrunk negative-binomial effects for the full data. q values use the complete 232-test panel family. Exclusions remove the older FLT ISS-T support; the full SRG result does not establish a region-specific mechanism.')

b=pd.read_csv(R/'known_truth_benchmark_summary.csv'); repetitions=pd.read_csv(R/'known_truth_benchmark_repetitions.csv');scenarios=['pure_null','confounded_no_biology','orthogonal_biology','overlapping_biology','changed_coupling']
labels=['Null','Nuisance only','Separate biology','Overlapping biology','Changed coupling']
fig,axs=plt.subplots(2,2,figsize=(11.8,8.8))
metrics=[('biological_effect_RMSE','Effect RMSE'),('biological_false_discovery_proportion','Biological FDP'),('power','Power'),('true_biology_energy_retained_fraction','Known biology energy retained')]
for k,(metric,label) in enumerate(metrics):
    ax=axs.flat[k]
    for j,m in enumerate(C):
        s=b[b.method==m].set_index('scenario').loc[scenarios];x=np.arange(5)+(j-1.5)*.14
        for position,scenario in zip(x,scenarios):
            values=repetitions.loc[(repetitions.method==m)&(repetitions.scenario==scenario),metric].dropna().to_numpy()
            if len(values)>1 and np.ptp(values)>1e-12:
                violin=ax.violinplot(values,positions=[position],widths=.12,showextrema=False)
                for body in violin['bodies']:body.set_facecolor(C[m]);body.set_edgecolor(C[m]);body.set_alpha(.22);body.set_linewidth(1)
        ax.errorbar(x,s[metric+'_mean'],yerr=1.96*s[metric+'_mc_se'],fmt='o-',color=C[m],lw=1.2,capsize=2,markersize=4,label=N[m])
    ax.set_xticks(range(5),labels,rotation=25,ha='right',fontsize=12);ax.set_ylabel(label)
    ax.set_title(chr(97+k)+'  '+label,loc='left',fontweight='bold');ax.grid(axis='y',alpha=.15)
    if k in [1,2]: ax.set_ylim(-.03,1.04)
    if k==3: ax.axhline(1,color='#aaa',ls='--',lw=.8)
axs[0,0].legend(frameon=False,fontsize=12)
fig.tight_layout(h_pad=2,w_pad=2)
save(fig,'figure6_known_truth_validation','Figure 6. Independent training–target Gaussian mechanism benchmark. Five scenarios each contain 200 independently generated datasets (1,000 datasets; 7,000 method rows including additional oracle-ground and partial-correction references in the supplement). Four operational methods are shown without target tuning. Translucent violins show the full repetition distributions where variable; points and bars are means ±1.96 Monte Carlo standard errors, not real-data confidence intervals. Biological FDP compares discoveries with the original biological coefficient vector and therefore differs from conditional statistical FDR under residual nuisance. Power and biological energy retention are undefined when no biological effects exist. Failure under overlapping biology or changed coupling is retained. The simulator is not a complete RNA-seq generative model.')

geom=pd.read_csv(R/'operator_geometry.csv'); anchors=pd.read_csv(R/'ground_anchor_sensitivity.csv');joint=pd.read_csv(R/'joint_bootstrap_summary.csv')
fig,axs=plt.subplots(1,3,figsize=(12,5.6),gridspec_kw={'width_ratios':[1,.95,1.45]})
for j,n in enumerate(['OSD-561','OSD-562','OSD-525']):
    s=geom[geom.anchor_accession==n].score_spectral_gain
    axs[0].scatter(np.full(len(s),j)+np.linspace(-.12,.12,len(s)),s,color=C['score'],s=28)
axs[0].axhline(2,color=C['bounded_score'],lw=2,label='Bounded gain ≤2')
axs[0].set_xticks(range(3),['Cb 561','Hip 562','Cb 525'],rotation=20);axs[0].set_ylabel('Maximum spectral amplification');axs[0].set_ylim(0,17)
axs[0].set_title('a  Real-data geometry',loc='left',fontweight='bold');axs[0].legend(frameon=False,fontsize=12)
comparisons=['anchor_OSD562_to_OSD562_FLT_recovery','OSD525_shared_HGC_anchor']
# Names are explicit in the output manifest; do not substitute another contrast.
if comparisons[1] not in set(anchors.comparison):
    comparisons[1]=next(n for n in anchors.comparison.unique() if n.startswith('anchor_OSD525_to_OSD525'))
for j,m in enumerate(['projection','score','bounded_score']):
    s=anchors[(anchors['mode']==m)&anchors.comparison.isin(comparisons)].set_index('comparison').loc[comparisons]
    axs[1].plot([0,1],s.sensor_effect_energy_removed_pct,'o-',color=C[m],label=N[m],lw=1.2)
axs[1].set_xticks([0,1],['Hip recovery','RR-23 shared'],rotation=20);axs[1].set_ylabel('Sensor effect energy removed (%)');axs[1].set_ylim(-5,105)
axs[1].set_title('b  Sensor attenuation',loc='left',fontweight='bold');axs[1].legend(frameon=False,fontsize=11)
metric=next(s for s in joint.metric.unique() if 'outside' in s.lower())
for j,n in enumerate(['OSD-561','OSD-562','OSD-525']):
    for k,m in enumerate(['projection','score','bounded_score']):
        s=joint[(joint.accession==n)&(joint['mode']==m)&(joint.metric==metric)].iloc[0];y=j+(k-1)*.18
        axs[2].hlines(y,s.bootstrap_lower95,s.bootstrap_upper95,color=C[m],lw=1.6)
        axs[2].scatter(s.point_estimate,y,color=C[m],s=28,zorder=4)
axs[2].set_yticks(range(3),['Cb recovery','Hip recovery','RR-23 shared HGC']);axs[2].invert_yaxis();axs[2].axvline(0,color='#999',lw=.8)
axs[2].set_xlabel('Outside-sensor effect energy removed (%)');axs[2].set_title('c  Joint sampling uncertainty',loc='left',fontweight='bold')
fig.tight_layout(w_pad=2)
save(fig,'figure7_bounded_operator_and_joint_uncertainty','Figure 7. Bounded score geometry and joint uncertainty. (a) Unbounded score gain for all 15 aligned anchor–target geometries; the bounded operator has gain at most two. (b) Observed sensor-panel effect-energy removal for matched hippocampal recovery and the RR-23 shared-HGC comparison. (c) Observed outside-sensor energy removal and 500-bootstrap percentile intervals, refitting normalization, anchor and target on each resample. Shared RR-23 HGC animals are sampled once for both roles. Intervals remain conditional on the fixed initially retained gene universe and observational sampling assumptions; zero crossings and negative reductions are shown. A point estimate may lie outside a biased percentile interval. Energy removal does not measure biological recovery.')

for source,name,caption in [
('extensions/paired_biology/figure_paired_regional_interactions','figure8_paired_regional_interactions','Figure 8. Paired regional interactions. (a) Source Name matching identifies 24 RRRM-2 and 10 RR-10 cerebellum–hippocampus animal pairs. (b) All eight panel effects in primary regional contrasts; no test has q≤0.05 across the complete 104-test family, including sensitivities. (c) Individual SRG hippocampus-minus-cerebellum scores; outlined points identify metadata conflicts. (d) SRG regional flight coefficients with 2,000-animal-bootstrap percentile intervals and complete-family q values. The approximate three-way interaction uses 9,999 restricted-model wild draws. Metadata exclusion removes all older FLT ISS-T animals, reducing the sensitivity to younger-animal exposure support.'),
('extensions/mathematics_review/high_gain_benchmark_panel','figure9_high_gain_domain_shift','Figure 9. Exploratory high-gain stress benchmark. All 16 combinations of true gain (8 or 12), training size (24 or 96), target size (24 or 48) and unchanged or flipped outside-sensor nuisance coupling are shown. Each uses 100 independent repetitions (1,600 datasets; 6,400 method rows). Known biology is supported outside Sensor31 and is preserved exactly by both score operators. Error bars are ±1.96 Monte Carlo standard errors. Bounding reduces worst score amplification under coupling reversal but leaves substantial nuisance, including under unchanged coupling; neither norm control nor exact sensor-null preservation establishes accurate biological recovery.')]:
    for ext in ['png','svg']: shutil.copy2(HERE/(source+'.'+ext),F/(name+'.'+ext))
    manifest.append(dict(filename=name,caption=caption,png_dpi=600,vector_format='SVG'))

rep=pd.read_csv(R/'regional_gene_replication_summary.csv')
fig,axs=plt.subplots(1,2,figsize=(10.8,4.7))
axs[0].bar(np.arange(3)-.18,rep.n_genes,width=.35,color='#c2cbd1',label='Tested')
axs[0].bar(np.arange(3)+.18,rep.direction_concordant,width=.35,color=C['projection'],label='Same direction')
axs[0].set_xticks(range(3),['Cerebellum','Hippocampus','Left hemisphere']);axs[0].set_ylabel('Region–gene combinations');axs[0].legend(frameon=False,fontsize=12)
axs[0].set_title('a  Cross-cohort direction screen',loc='left',fontweight='bold')
axs[1].bar(np.arange(3)-.18,rep.nominal_max_p05,width=.35,color=C['score'],label='Nominal max-p <0.05')
axs[1].bar(np.arange(3)+.18,rep.BH_FDR05,width=.35,color=C['bounded_score'],label='Global BH q <0.05')
axs[1].set_xticks(range(3),['Cerebellum','Hippocampus','Left hemisphere']);axs[1].set_ylabel('Candidates');axs[1].legend(frameon=False,fontsize=12)
axs[1].set_title('b  No corrected replicated gene',loc='left',fontweight='bold')
fig.tight_layout(w_pad=2)
save(fig,'figureS1_cross_cohort_gene_screen','Supplementary Figure 1. Conservative regional gene reproducibility screen. A candidate must have the same effect direction in every constituent contrast and finite component Wald p values; its intersection–union p value is the maximum component p. Discordant or unavailable effects receive p=1. All 54,989 region–gene combinations enter one global BH and BY family. No gene survives either correction. Component comparisons come from the same negative-binomial fits and can share samples; concordant direction alone is not replication evidence. Ranked nominal candidates are explicitly exploratory.')
(R/'figure_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Saved',len(manifest),'figures')
