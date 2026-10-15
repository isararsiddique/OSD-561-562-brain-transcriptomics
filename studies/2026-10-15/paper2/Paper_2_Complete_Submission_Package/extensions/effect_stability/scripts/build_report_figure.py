from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
B=Path(__file__).resolve().parents[3];E=B/'extensions/effect_stability';R=E/'results';F=E/'figures'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':8,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':1,'lines.linewidth':1,'axes.spines.top':False,'axes.spines.right':False})
def panel(ax,l,t):ax.set_title(l+') '+t,loc='left',fontweight='bold',pad=10)
def main():
 e=pd.read_csv(R/'BY_resolution_envelopes.csv');p=pd.read_csv(R/'all_complementary_half_concordance.csv.gz');ef=pd.read_csv(R/'all_frozen_effect_vectors_split_metrics.csv.gz');loo=pd.read_csv(R/'all_animal_leave_one_out_native_effects.csv.gz');s=pd.read_csv(R/'effect_stability_summary.csv');ver=json.loads((R/'verification.json').read_text());assert ver['status']=='PASS'
 fig,axs=plt.subplots(2,2,figsize=(183/25.4,161/25.4));fig.subplots_adjust(left=.10,right=.97,top=.93,bottom=.16,hspace=.70,wspace=.45)
 ax=axs[0,0];panel(ax,'a','Discovery decisions across disjoint halves')
 xs=np.arange(1,11);a=[e.set_index('family').loc[f'split_{i:02d}_A','BY_positive'] for i in xs];b=[e.set_index('family').loc[f'split_{i:02d}_B','BY_positive'] for i in xs]
 ax.bar(xs-.18,a,width=.34,color='#D55E00',label='Original half A');ax.bar(xs+.18,b,width=.34,color='#0072B2',label='Complement B');ax.plot(xs-.18,np.zeros(10),'_',color='#D55E00',ms=5);ax.set_xticks([1,3,5,7,9]);ax.set_ylim(-1,30);ax.set_xlabel('Frozen read split');ax.set_ylabel('BY-positive contrast–cluster tests');ax.legend(frameon=False,loc='upper left',handlelength=1)
 ax=axs[0,1];panel(ax,'b','Split-half effect direction agreement')
 pp=p[p.both_halves_measured];allv=pp.cosine.to_numpy();hits=pp[pp.primary_BY_test].cosine.to_numpy()
 for vals,col,lab in [(allv,'#8299A6','All measured tests'),(hits,'#0072B2','Original BY candidates')]:
  v=np.sort(vals);ax.plot(v,np.arange(1,len(v)+1)/len(v),c=col,label=lab)
 ax.axvline(0,c='#666666',ls='--');ax.set_xlim(-1.03,1.03);ax.set_ylim(0,1.03);ax.set_xlabel('Cosine between disjoint-half effects');ax.set_ylabel('Cumulative fraction of instances');ax.legend(frameon=False,loc='upper left',handlelength=1);ax.text(.04,.49,f'{len(hits)} measured candidate pairs\nAll retain positive direction',transform=ax.transAxes,fontsize=8)
 ax=axs[1,0];panel(ax,'c','Candidate effect magnitude and influence')
 for vals,col,lab in [(ef.loc[ef.primary_BY_test&ef.measurement_valid,'relative_projection'],'#0072B2','Read halves'),(loo.loc[loo.primary_BY_test,'relative_projection'],'#D55E00','One animal omitted')]:
  v=np.sort(vals.to_numpy());ax.plot(v*100,np.arange(1,len(v)+1)/len(v),c=col,label=lab)
 ax.axvline(100,color='#666666',ls='--');ax.set_xlim(25,180);ax.set_ylim(0,1.03);ax.set_xlabel('Projection onto original effect (%)');ax.set_ylabel('Cumulative fraction of instances');ax.legend(frameon=False,loc='lower right',handlelength=1)
 ax=axs[1,1];panel(ax,'d','Discrete p values create a BY boundary')
 xx=np.linspace(10,33,300);scale=2769*sum(1/i for i in range(1,2770))/24310;ax.plot(xx,scale/xx,color='#777777',lw=1,label='Floor block: m H pmin / k');ax.axhline(.05,color='black',ls='--',label='BY 0.05');ax.axvline(20,color='#BBBBBB',ls=':')
 for suffix,col,marker in [('A','#D55E00','o'),('B','#0072B2','s')]:
  r=e[e.family.str.endswith('_'+suffix)];ax.scatter(r.n_at_minimum_p,r.minimum_BY_q,s=22,c=col,marker=marker,zorder=3)
 r=e.iloc[0];ax.scatter(r.n_at_minimum_p,r.minimum_BY_q,s=60,c='#009E73',marker='*',zorder=4);ax.annotate('Original',xy=(30,r.minimum_BY_q),xytext=(25,.028),fontsize=8)
 ax.set_xlim(10,33);ax.set_ylim(.022,.105);ax.set_xlabel('Tests tied at the exact p-value floor');ax.set_ylabel('Minimum BY q');ax.legend(frameon=False,fontsize=8,loc='upper right',handlelength=1)
 fig.text(.10,.05,'All 923 clusters remain in each 2,769-test family. Halves share the same 26 animals.\nCandidate diagnostics are post hoc; repeated splits and omissions are not biological replication.',fontsize=8,va='bottom')
 out=F/'Figure_Effect_Stability'
 for ext in ['pdf','svg','png']:fig.savefig(out.with_suffix('.'+ext),dpi=600,facecolor='white')
 plt.close(fig);im=Image.open(out.with_suffix('.png')).convert('RGB');im.save(out.with_suffix('.png'),dpi=(600,600))
 caption='Read-split repeatability and discrete discovery thresholds. a) All ten pre-existing A halves and all ten complementary B halves, obtained by subtracting A counts from the original matrix; each frozen family retains all 2,769 slots, with p=1 for failed measurements. A halves have zero BY calls; five B halves have 20–23 calls and five have zero. All 60 global comparisons retain within-half Holm significance. b) Native contrast-vector cosine between the two disjoint halves, displayed as empirical cumulative distributions for every jointly measured contrast–cluster instance and the original 43 BY-test subset. All 354 jointly measured candidate instances retain positive orientation; medians are 0.9964 for FLT–GC and 0.9933 for FLT–VIV. c) Projection onto each original native effect, divided by its squared norm, for all 751 measurable candidate half-instances and all 771 candidate animal omissions. These are descriptive repeated-instance distributions, not confidence intervals or independent samples. d) All 20 half-families and the original family at the exact minimum p=1/24,310. The curve is mH_m p_min/k for a tied minimum-p block, m=2,769; 20 floor-level tests are needed for that block to cross 0.05. The complete ordered families were checked; higher ranks can also determine a BY rejection in general. Stable direction does not establish a reliable thresholded biomarker. Halves partition the same observed reads, have dependent conditional sampling noise and share biological animals. Candidate selection is post hoc, and no threshold or family was changed.'
 (F/'Figure_Effect_Stability_caption.txt').write_text(caption+'\n')
 assert len(caption.split())<350
 def tab(d):return d.to_markdown(index=False) if False else '\n'.join(['| '+' | '.join(d.columns)+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(str(v) for v in r)+' |' for r in d.itertuples(index=False,name=None)])
 report='''# Effect stability is distinct from thresholded discovery stability

This retrospective diagnostic uses all923 frozen multialternative clusters and all2769 tests. The ten original read-thinning outcomes were known when this follow-up was planned. Every corresponding complementary half was analyzed without outcome selection; no new favorable thinning seed was chosen.

## Concrete result

The original ten A halves still have zero local BY discoveries. Among the complementary B halves, five have20,20,23,21 and20 calls and five have zero:104 contrast–cluster instances, not104 unique independent discoveries. All60 global comparisons across20halves retain their separate three-contrast Holm corrections. Thus the earlier claim about the original ten thinnings remains true, but a universal claim that no half-depth partition can retain local calls is unsupported.

The original candidate effects remain highly repeatable within these data. Among jointly measurable candidate split pairs, all354 retain positive alignment:325FLT–GC and29FLT–VIV instances. Median pairwise cosines are0.99643 and0.99330. All751 measurable original-to-half candidate instances and all771 single-animal-omission candidate instances retain positive projections. The whole923-cluster panel includes weak effects and lower concordance; all outcomes are preserved, not just selected candidates.

## Meaning and limits

This separates three questions: whether a measurement is retained, whether its observed effect direction/magnitude changes, and whether a multiplicity-corrected decision crosses0.05. None of those questions is independent biological replication. Disjoint read halves share the same animals; conditional on the original count their noises are negatively correlated. Repeated omissions also share animals. The original candidate subset was selected for strong effects and its greater technical repeatability is subject to that selection. No p value or effective biological sample size is assigned to these instance-level distributions. The result supports technical repeatability of candidate point estimates, not a stable externally validated biomarker panel.

## Mathematics

For cluster c, let d_c be the original native mean-composition contrast and h_c its perturbed counterpart. Cosine is d_c'h_c/(||d_c||||h_c||). Relative projection is d_c'h_c/||d_c||²;100% retains the original directed magnitude,0% loses it, and a negative value reverses direction. Native total-variation magnitude is50sum_j|h_cj| percentage points. All metrics use frozen alternatives. Direction and magnitude are evaluated separately from the original square-root exact test.

For ordered raw p values, BY rejects iff max_i [i alpha/(mH_m p_(i))]>=1. With m=2769,H_m=8.50363775 and the minimum exact p=1/24310, the floor-block adjustment is0.96859617/k. Twenty tied minimum-p tests giveq=.04842981; nineteen giveq=.05097875. Every positive B-half decision here is supported by20–23floor-level tests, while A halves have13–19. The complete rank envelopes are supplied. This is not permission to discard hypotheses or relax correction. It describes a discrete decision boundary, not a proof that all nonsignificant candidates are true or a quantitative estimate of power against an unknown truth.

## Full summaries

'''+tab(s.round(6))+'''\n
## Verification and reproduction

Run scripts/analyse.py, scripts/verify.py and scripts/build_report_figure.py. The verifier independently checks all count partitions, selected coordinate reconstructions,405metric/correction/exact-permutation checks, complete dimensions and retained candidate identities. Original-half exact statistics were separately reproduced against the pre-existing tables. All results, ordered p curves, gates, Tcaf1 cases and negative values are retained. The analysis does not infer causal microgravity effects or resolve the exposure/date alias.
'''
 (E/'SCIENTIFIC_REPORT.md').write_text(report)
 (R/'figure_metadata.json').write_text(json.dumps(dict(width_mm=183,height_mm=161,min_font_pt=8,dpi=600,caption_words=len(caption.split()),native_vector=True),indent=2))
 print('Figure/report ready; caption',len(caption.split()),'words')
if __name__=='__main__':main()
