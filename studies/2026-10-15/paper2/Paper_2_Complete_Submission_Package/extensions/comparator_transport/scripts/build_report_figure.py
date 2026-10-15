from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
B=Path(__file__).resolve().parents[3];E=B/'extensions/comparator_transport';R=E/'results';F=E/'figures'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':8,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':1,'lines.linewidth':1,'axes.spines.top':False,'axes.spines.right':False})
def panel(ax,l,t):ax.set_title(l+') '+t,loc='left',fontweight='bold',pad=10)
def table(d):return '\n'.join(['| '+' | '.join(d.columns)+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(str(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)])
def main():
 assert json.loads((R/'independent_verification.json').read_text())['status'] in ['PASS','passed']
 g=pd.read_csv(R/'all_comparator_global_transport.csv');primary=g[g.panel=='cohort_measured'];animal=pd.read_csv(R/'all_animal_scores.csv');tc=pd.read_csv(R/'Tcaf1_all_comparators_and_cohorts.csv');training=pd.read_csv(R/'RR23_only_frozen_directions.csv.gz');summary=json.loads((R/'summary.json').read_text());trows=[]
 for tag in ['RRRM2','RR10']:
  z=np.load(R/f'{tag}_frozen_measurements.npz',allow_pickle=True);off=z['offsets'];i=list(z['cluster_ids']).index('6:2:42663893:start');a,b=off[i:i+2]
  for direction in ['FLT_minus_GC','FLT_minus_VIV']:
   tr=training[training.direction==direction].iloc[a:b];score=tr.RR23_unit_vector.to_numpy()@(z['proportions'][a:b]-tr.RR23_comparator_mean.to_numpy()[:,None]);w=np.load(R/f'{tag}_{direction}_bootstrap_weights.npz')['validation_weights'];boot=w@score;rec=tc[(tc.cohort==tag)&(tc.direction==direction)].iloc[0]
   trows.append(dict(cohort=tag,direction=direction,effect_pp=100*rec.effect,conditional_low_pp=100*np.quantile(boot,.025),conditional_high_pp=100*np.quantile(boot,.975),exact_p=rec.exact_p,joint_BY=rec.BY_joint3692))
 tci=pd.DataFrame(trows);tci.to_csv(R/'Tcaf1_conditional_intervals.csv',index=False)
 fig,axs=plt.subplots(2,2,figsize=(183/25.4,162/25.4));fig.subplots_adjust(left=.17,right=.97,top=.93,bottom=.23,hspace=.77,wspace=.65)
 ax=axs[0,0];panel(ax,'a','Three cohorts and measurement overlap');v=[923,281,702];labels=['RR23 discovery\n26 animals','RRRM2 external\n19 animals','RR10 external\n10 animals'];cols=['#8299A6','#0072B2','#D55E00'];ax.barh(np.arange(3),v,color=cols);ax.set_yticks(range(3),labels);ax.invert_yaxis();ax.set_xlim(0,1060);ax.set_xlabel('Measurable frozen clusters')
 for j,n in enumerate(v):ax.text(n+20,j,str(n),va='center')
 ax.text(.0,-.30,'270 clusters measurable in both external cohorts',transform=ax.transAxes,fontsize=8)
 ax=axs[0,1];panel(ax,'b','Global directional transport');labels=[]
 for j,r in enumerate(primary.itertuples()):
  col='#0072B2' if r.cohort=='RRRM2' else '#D55E00';label=r.cohort+' | '+('FLT–GC' if r.direction.endswith('GC') else 'FLT–VIV');labels.append(label);ax.plot([r.conditional_CI95_low_pp,r.conditional_CI95_high_pp],[j,j],c=col,lw=1.5);ax.scatter(r.effect_pp,j,c=col,s=23);ax.text(2.2,j,f'p={r.exact_p:.3f}',va='center',fontsize=8)
 ax.axvline(0,c='#666666',ls='--');ax.set_yticks(range(4),labels);ax.set_ylim(3.6,-.6);ax.set_xlim(-1.4,3.5);ax.set_xlabel('Projected FLT–GC effect (pp)');ax.text(.01,-.29,'All four Holm-adjusted p = 1',transform=ax.transAxes,fontsize=8)
 ax=axs[1,0];panel(ax,'c','Every RR10 animal on the frozen score');f=animal[(animal.cohort=='RR10')&(animal.direction=='FLT_minus_GC')&(animal.panel=='cohort_measured')];center=f.score.mean()
 for j,(ex,col) in enumerate([('GC','#E69F00'),('FLT','#0072B2')]):
  vals=100*(f[f.exposure==ex].score.to_numpy()-center);ax.scatter(j+np.linspace(-.1,.1,len(vals)),vals,c=col,s=25);ax.hlines(vals.mean(),j-.17,j+.17,color='black',lw=1.2)
 ax.axhline(0,c='#CCCCCC',lw=1);ax.set_xticks([0,1],['Ground (n=5)','Flight (n=5)']);ax.set_xlim(-.55,1.55);ax.set_ylabel('Centered frozen global score (pp)');ax.text(.0,-.30,'Distinct animals; shares RR23 launch',transform=ax.transAxes,fontsize=8)
 ax=axs[1,1];panel(ax,'d','Tcaf1 external estimates')
 for j,r in enumerate(tci.itertuples()):
  col='#0072B2' if r.cohort=='RRRM2' else '#D55E00';ax.plot([r.conditional_low_pp,r.conditional_high_pp],[j,j],c=col,lw=1.5);ax.scatter(r.effect_pp,j,c=col,s=23)
 ax.axvline(0,c='#666666',ls='--');ax.set_yticks(range(4),labels);ax.set_ylim(3.6,-.6);ax.set_xlabel('Tcaf1 directed effect (pp)');ax.text(.01,-.30,'All local joint BY q = 1',transform=ax.transAxes,fontsize=8)
 fig.text(.12,.045,'RRRM2 uses an independent launch; RR10 shares SpaceX-21 with RR23.\nDirections use RR23 alone. Points/intervals show effects and conditional 95% animal-bootstrap intervals.',fontsize=8,va='bottom')
 out=F/'Figure_Comparator_Transport'
 for ext in ['pdf','svg','png']:fig.savefig(out.with_suffix('.'+ext),dpi=600,facecolor='white')
 plt.close(fig);im=Image.open(out.with_suffix('.png')).convert('RGB');im.save(out.with_suffix('.png'),dpi=(600,600))
 caption='Comparator-specific external transport in 55 animals from three cohorts and two launches. a) Frozen RR23 multialternative panel and measurable overlap in each external cohort; the common external intersection contains 270 clusters. Gates require at least 20 selected reads and 90% dominance in every included animal. RR10 has 5 FLT and 5 GC distinct source animals but shares the RR23 launch; RRRM2 has 10 FLT and 9 GC animals from an independent launch. b) Equal-cluster scores trained using RR23 FLT–GC or FLT–VIV native-composition unit directions. Points are external FLT–GC effects; segments are conditional 95% percentile intervals from 4,000 whole-animal bootstrap draws with frozen directions. Displayed p values are exact positive-transfer tests; all four Holm-adjusted values are 1. RRRM2 allocations preserve three route/age strata (11,200 assignments); RR10 uses all 252 five-versus-five assignments. c) All ten RR10 animal scores using the FLT–GC-trained direction, centered at the overall RR10 mean for display; black segments show group means. Centering does not change the contrast. d) Tcaf1 projected effects under both frozen comparator directions and the same conditional intervals; all joint BY values are 1 across the full 3,692-slot local family. Local Holm rejection is mathematically unattainable at the available resolution and is not evidence of biological absence. The common-270-cluster sensitivity, two-cohort bootstrap intervals and every unsuccessful test are retained in Supplementary Data. No comparator was selected by external performance. Cohorts differ in sex, strain, collection and other context; exposure/date aliasing persists.'
 assert len(caption.split())<=350;(F/'Figure_Comparator_Transport_caption.txt').write_text(caption+'\n')
 report='''# Comparator-specific transport across three cerebellar cohorts

This extension adds real publisher-verified RR10data from10distinct source animals and tests both discovery ground-comparator directions in both external cohorts. It follows a plan frozen before the new outcomes, while acknowledging that the previous RRRM2 FLT–GC-trained result was already known. RR10 was selected on tissue, source-animal and processed-data criteria, not an observed effect. The total55animals comprise26RR23,19eligibleRRRM2 and10RR10 animals across two launches, not three independent missions.

## Complete global outcomes

'''+table(g.round(7))+'''

No primary or intersection-panel global test passes its four-test Holm family. RR10 measures702/923clusters (76.06%), compared with281/923 (30.44%) in RRRM2. Thus unsuccessful transfer is also observed in a cohort with substantially larger measurement overlap. The RR10 FLT–GC-trained effect is+0.4358pp,95%conditional interval[-1.0133,2.0329],exact positive-transferp=.31349. FLT–VIV-trained direction does not rescue transfer in either cohort. These nonsignificant results do not prove a universal effect is absent or identify why cohorts differ. The intervals preserve uncertainty.

## Design and causal limits

RR10/OSD-563 is female wild-typeB6129SF2/J,5FLT/5GC,14–15weeks at launch, with marker injections shared across arms. RR23 is maleC57BL/6J. Both were launched on SpaceX21; their source animals/RFIDs do not overlap. RR10 flight animals were euthanized02January2021 and controls05January; carcass dissections were12May and14May. These context differences and exposure/date aliasing prevent causal interpretation. All10STARdata files and the processing ZIP match publisher MD5; software-version sidecar has a malformed33-character digest that is preserved as a provenance exception. Verified processing commands and frozen reference metadata establish coordinate compatibility. See the complete cohort-search record for inclusions/exclusions.

RRRM2 is an independent-launch cohort with19eligible animals and the previously documented full-stratum exclusion. Its eligibility and stratification were unchanged. RR10 provides distinct-animal external evidence, not new independent-launch replication.

## Exact frozen analysis

All923original clusters and selected alternatives remain frozen. Label-blind measurement gates are evaluated separately in each external cohort. For comparatorC, d_c is the RR23 FLT−C mean native composition contrast, u_c=d_c/||d_c||. Each validation animal receives u_c'*(p_i−mean_RR23_C); the common reference centering cancels from group contrasts. The score averages equally over the cohort's measurable clusters. No validation outcome enters training weights or comparator choice. Four global one-sided tests form one new retrospective Holm family; the270cluster intersection is a separate complete four-test sensitivity.

RR10 uses252exact allocations. RRRM2 uses11200allocations preserving three route-by-age strata and weights each stratum by its fraction of all19animals. Ties use1e−12tolerance; absolute two-sided tests are also retained. The new local family contains4×923=3692slots, assigningp=1 to unmeasurable clusters, with jointBH/BY/Holm. NeitherBHnorBY yields a local discovery. For Holm, even the finest p minimum1/11200 exceeds.05/3692; the first adjusted minimum is.32964. Holm rejection is therefore impossible throughout the joint family, and that zero count is not biological evidence.

Four thousand bootstrap draws resample biological animals within group/stratum; primary intervals hold RR23directions fixed. A separate two-cohort sensitivity re-estimates directions from RR23animal resamples without fitting to validation outcomes. All bootstrap weights and draws are saved. Intervals are approximate/pointwise and nonsignificance is not equivalence.

## Tcaf1 and complete outputs

'''+table(tci.round(7))+'''

Tcaf1 is measurable in both external cohorts but neither frozen direction supports a transferred association. Its complete fourtest row family remains part of all3692tests, not a selected standalone test. Candidate intervals are descriptive and unadjusted for selection.

## Verification

218independent checks reconstruct every measurement, native direction, global and local exact allocation result, all64000bootstrap effect draws, intervals and multiplicity arithmetic from frozen count matrices/metadata. The verifier does not import analysis functions. All sources, unfavourable outcomes, animal scores, gates and preanalysis decisions are retained. Run scripts/analyse.py, scripts/verify_independently.py, scripts/build_report_figure.py. These checks establish numerical consistency, not causal identification or successful replication.
'''
 (E/'SCIENTIFIC_REPORT.md').write_text(report);(R/'figure_metadata.json').write_text(json.dumps(dict(width_mm=183,height_mm=162,min_font_pt=8,dpi=600,caption_words=len(caption.split())),indent=2));print('Figure/report ready, caption words',len(caption.split()))
if __name__=='__main__':main()
