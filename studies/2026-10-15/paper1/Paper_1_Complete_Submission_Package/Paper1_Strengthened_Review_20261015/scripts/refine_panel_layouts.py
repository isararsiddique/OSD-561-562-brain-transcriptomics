#!/usr/bin/env python3
"""Re-render paired and high-gain panels from cached tables at readable type sizes."""
import itertools
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parents[1];F=HERE/'figures'
plt.rcParams.update({'font.family':'Arial','font.size':14,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
res=pd.read_csv(HERE/'extensions/paired_biology/paired_panel_interactions.csv'); scores=pd.read_csv(HERE/'extensions/paired_biology/paired_animal_panel_scores.csv')
fig,axs=plt.subplots(2,2,figsize=(12,10.1),layout='constrained');groups=['GC LAR','FLT LAR','GC ISS-T','FLT ISS-T'];colors=['#4c78a8','#f28e2b','#4c78a8','#f28e2b']
sr=scores[(scores.mission=='RRRM-2')&(scores.panel=='SRG')];ax=axs[0,0];sizes=[len(sr[sr.group==x])for x in groups]
ax.bar(range(4),sizes,color=colors,width=.65)
for i,n in enumerate(sizes):ax.text(i,n+.1,str(n),ha='center')
ax.set_xticks(range(4),groups,rotation=20,ha='right');ax.set_ylim(0,11);ax.set_ylabel('Matched animals');ax.set_title('a  Animal pairing',loc='left',fontweight='bold')
ax.text(.02,.97,'RRRM-2: 24 of 29 animals\nRR-10: 10 pairs',transform=ax.transAxes,va='top',fontsize=13)
labels={'RRRM2_region_by_flight_ISST':'Flight, ISS-T','RRRM2_region_by_flight_LAR':'Flight, LAR','RRRM2_region_by_cohort_GC':'Cohort, GC','RRRM2_region_by_cohort_FLT':'Cohort, FLT','RRRM2_region_by_flight_cohort_adjusted':'Flight, pooled','RR10_region_by_flight':'Flight, RR-10','RRRM2_region_flight_cohort_triple':'Three-way'}
primary=res[res.sensitivity=='all_matched_animals'];cols=['rPRG','dPRG','SRG','HSP','Glia_Complement','BBB_Endothelium','Circadian','Sensor31'];mat=primary.pivot(index='contrast',columns='panel',values='paired_interaction_log2').loc[list(labels),cols];ax=axs[0,1];vm=abs(mat.to_numpy()).max()
im=ax.imshow(mat,cmap='RdBu_r',vmin=-vm,vmax=vm,aspect='auto');ax.set_yticks(range(7),list(labels.values()),fontsize=13);ax.set_xticks(range(8),['rPRG','dPRG','SRG','HSP','Glia','BBB','Clock','Sensor'],rotation=40,ha='right',fontsize=12)
for i in range(7):
 for k in range(8):ax.text(k,i,f'{mat.iloc[i,k]:.2f}',ha='center',va='center',fontsize=12,color='white'if abs(mat.iloc[i,k])>.65*vm else 'black')
ax.set_title('b  Regional effects\nAll 104 tests: q >0.05',loc='left',fontweight='bold',fontsize=14);fig.colorbar(im,ax=ax,shrink=.7,label='Log2 interaction')
ax=axs[1,0];rng=np.random.default_rng(20261014)
for i,g in enumerate(groups):
 s=sr[sr.group==g];ax.scatter(i+rng.uniform(-.13,.13,len(s)),s.paired_score_hip_minus_cb,color=colors[i],s=45,edgecolors=np.where(s.location_conflict,'black',colors[i]),lw=1.3);ax.hlines(s.paired_score_hip_minus_cb.mean(),i-.22,i+.22,color='black',lw=1.8)
ax.set_xticks(range(4),groups,rotation=20,ha='right');ax.set_ylabel('SRG Hip − Cb score');ax.set_title('c  Individual animals',loc='left',fontweight='bold');ax.set_ylim(top=3.75)
ax.text(.02,.98,'Outline: metadata conflict\nLine: group mean',transform=ax.transAxes,va='top',fontsize=12)
ax=axs[1,1]
for i,name in enumerate(['RRRM2_region_by_flight_ISST','RRRM2_region_by_flight_LAR','RR10_region_by_flight']):
 for sens,o,c in [('all_matched_animals',0,'#4c78a8'),('exclude_metadata_conflicts',.22,'#b279a2')]:
  ss=res[(res.contrast==name)&(res.panel=='SRG')&(res.sensitivity==sens)]
  if ss.empty:continue
  s=ss.iloc[0];y=i+o;ax.hlines(y,s.bootstrap_lower95,s.bootstrap_upper95,color=c,lw=1.8);ax.scatter(s.paired_interaction_log2,y,color=c,s=38,label=('Full'if o==0 else'Exclusions')if i==0 else None);ax.text(.53,y,f'q={s.q_BH_all_panel_tests:.3f}',va='center',fontsize=12)
ax.set_yticks(range(3),['ISS-T','LAR','RR-10']);ax.invert_yaxis();ax.axvline(0,color='#aaa',ls='--');ax.set_xlim(-.75,1.12);ax.set_xlabel('SRG regional flight effect');ax.set_title('d  Regional flight contrasts',loc='left',fontweight='bold');ax.legend(frameon=False,fontsize=12,loc='upper left',bbox_to_anchor=(0,-.24),ncol=2)
for ext in ['png','svg']:fig.savefig(F/('figure8_paired_regional_interactions.'+ext),dpi=600,bbox_inches='tight')
plt.close(fig)
s=pd.read_csv(HERE/'extensions/mathematics_review/high_gain_benchmark_summary.csv');repetitions=pd.read_csv(HERE/'extensions/mathematics_review/high_gain_benchmark_repetitions.csv');fig,axes=plt.subplots(2,2,figsize=(12,9.2),layout='constrained');C={'raw':'#667788','projection':'#157c99','score':'#ad517a','bounded_score':'#d18a21'};N={'raw':'Raw','projection':'Projection','score':'Score','bounded_score':'Bounded L=2'}
for ax,(coupling,metric,title,letter) in zip(axes.flat,[('unchanged','biological_effect_RMSE','Stable coupling: effect error','a'),('flipped','biological_effect_RMSE','Flipped coupling: effect error','b'),('unchanged','biological_FDP','Stable coupling: biological FDP','c'),('flipped','biological_FDP','Flipped coupling: biological FDP','d')]):
 combos=list(itertools.product([8.,12.],[24,96],[24,48]));x=np.arange(8)
 for k,m in enumerate(C):
  p=s[(s.coupling==coupling)&(s.method==m)].set_index(['true_gain','n_train','n_target']).loc[combos]
  for position,(gain,nt,ne) in zip(x+(k-1.5)*.14,combos):
   values=repetitions.loc[(repetitions.coupling==coupling)&(repetitions.method==m)&(repetitions.true_gain==gain)&(repetitions.n_train==nt)&(repetitions.n_target==ne),metric].dropna().to_numpy()
   if len(values)>1 and np.ptp(values)>1e-12:
    violin=ax.violinplot(values,positions=[position],widths=.12,showextrema=False)
    for body in violin['bodies']:body.set_facecolor(C[m]);body.set_edgecolor(C[m]);body.set_alpha(.22);body.set_linewidth(1)
  ax.errorbar(x+(k-1.5)*.14,p[metric+'_mean'],yerr=1.96*p[metric+'_mc_se'],fmt='o',ms=4,capsize=2,color=C[m],label=N[m],lw=1)
 ax.set_xticks(x,[f'{g:g}\n{nt}/{ne}'for g,nt,ne in combos],fontsize=12);ax.set_xlabel('True gain\nTraining N / target N',fontsize=12);ax.set_ylabel('Effect RMSE' if 'RMSE'in metric else 'Biological FDP');ax.set_title(letter+'  '+title,loc='left',fontweight='bold',fontsize=14);ax.grid(axis='y',alpha=.15)
 if metric=='biological_FDP':ax.set_ylim(-.02,1.04)
handles,legend_labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,legend_labels,frameon=False,ncol=4,fontsize=12,loc='lower center',bbox_to_anchor=(.5,1.01))
for ext in ['png','svg']:fig.savefig(F/('figure9_high_gain_domain_shift.'+ext),dpi=600,bbox_inches='tight')
plt.close(fig);print('Refined figures 8 and 9 from saved tables.')
