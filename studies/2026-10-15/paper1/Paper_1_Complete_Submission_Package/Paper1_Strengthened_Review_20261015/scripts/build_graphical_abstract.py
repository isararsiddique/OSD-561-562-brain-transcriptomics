#!/usr/bin/env python3
"""Optional graphical overview with numerically verified real-data claims."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'publication_figures';out.mkdir(exist_ok=True)
h=pd.read_csv(ROOT/'extensions/ground_holdout/holdout_summary.csv');p=pd.read_csv(ROOT/'results/panel_results.csv');geom=pd.read_csv(ROOT/'results/operator_geometry.csv')
raw=h[h.method=='raw'].set_index('accession').loc[['OSD-561','OSD-562','OSD-525']];score=h[h.method=='score'].set_index('accession').loc[raw.index]
plt.rcParams.update({'font.family':'Arial','font.size':12,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(11.6,6.4),facecolor='white');fig.text(.055,.96,'Reliable collection scores do not ensure transcriptomic correction',fontweight='bold',fontsize=19,va='top',color='#18384c')
fig.text(.055,.885,'Mouse spaceflight brain reanalysis with paired sampling and bounded oblique projection',fontsize=12,color='#4d626f')
cols=['#326c9a','#277c61','#a65f29'];x=np.arange(3)
ax=fig.add_axes([.06,.36,.23,.4]);ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1);ax.text(0,1.02,'a  Public data and animal identity',fontsize=13,fontweight='bold',va='top')
for i,(name,n,region)in enumerate([('RRRM-2',47,'Core regions + other\nhemisphere animals'),('RR-10',10,'Matched cerebellum\nand hippocampus'),('RR-23',26,'Flight, habitat and\nvivarium controls')]):
 y=.86-i*.29;ax.scatter(.06,y-.075,s=230,color=cols[i]);ax.text(.18,y,name+'  '+str(n)+' source animals',fontweight='bold',fontsize=12,va='top');ax.text(.18,y-.10,region,fontsize=10,color='#4d626f',va='top',linespacing=1.2)
ax.text(0,-.12,'7 deposits  •  129 libraries\n83 animals  •  24 + 10 regional pairs',fontsize=11,fontweight='bold',linespacing=1.6,va='top')
ax1=fig.add_axes([.36,.53,.26,.23]);ax1.bar(x,raw.sensor_age_stratified_AUC_mean,color=cols,width=.58);ax1.set_ylim(0,1.18);ax1.set_yticks([0,.5,1]);ax1.set_xticks(x,['Cb','Hip','RR-23']);ax1.set_ylabel('Held-out sensor AUC',fontsize=10);ax1.set_title('b  Collection-cohort prediction',loc='left',fontsize=13,fontweight='bold');
for i,v in enumerate(raw.sensor_age_stratified_AUC_mean):ax1.text(i,v+.035,f'{v:.3f}',ha='center',fontsize=10)
ax2=fig.add_axes([.36,.2,.26,.23]);ax2.bar(x,score.outside_prediction_MSE_ratio_to_nuisance_only_mean,color=cols,width=.58);ax2.axhline(1,color='#666',ls='--',lw=1);ax2.set_ylim(.85,1.14);ax2.set_yticks([.9,1,1.1]);ax2.set_xticks(x,['Cb','Hip','RR-23']);ax2.set_ylabel('Loading MSE / baseline',fontsize=10);ax2.set_title('c  Outside-sensor prediction',loc='left',fontsize=13,fontweight='bold');
for i,v in enumerate(score.outside_prediction_MSE_ratio_to_nuisance_only_mean):ax2.text(i,v+.008,f'{v:.3f}',ha='center',fontsize=10)
ax3=fig.add_axes([.71,.27,.24,.49]);ax3.axis('off');ax3.set_xlim(0,1);ax3.set_ylim(0,1);ax3.text(0,1.02,'d  A transparent stability bound',va='top',fontsize=13,fontweight='bold')
ax3.text(0,.78,f'Real-data unbounded gain\n{geom.score_spectral_gain.min():.2f}–{geom.score_spectral_gain.max():.2f}',fontsize=14,fontweight='bold',linespacing=1.45,color='#ad517a')
ax3.add_patch(FancyArrowPatch((.25,.6),(.25,.48),arrowstyle='-|>',mutation_scale=18,color='#687a86',lw=1.5))
ax3.text(0,.37,'Bounded operator gain ≤ 2',fontsize=14,fontweight='bold',color='#177e99');ax3.text(0,.17,'Sensor score removal retained\nNuisance removal can be incomplete\nBiological recovery requires validation',fontsize=11,linespacing=1.6)
fig.text(.055,.08,'362 overlapping holdouts from 30 ground animals; training-only filtering and normalization.\nCore refitted-label test: exact p = 0.005; BY q = 0.039 across 18 tests.\nWithin-resource association; bounding limits amplification and does not prove recovered biology.',fontsize=10,color='#4d626f',linespacing=1.5)
for ext in ['png','svg']:fig.savefig(out/('GRAPHICAL_ABSTRACT.'+ext),dpi=600,bbox_inches='tight',facecolor='white')
plt.close(fig)
(out/'GRAPHICAL_ABSTRACT_CAPTION.md').write_text('Optional graphical overview. Collection scores generalize within the available ground resource, while prediction of outside-sensor expression by the fitted score loading can worsen. The 362 overlapping holdouts use 30 unique ground source animals; all filtering, normalization references and loading fitting use training animals only. Bars are descriptive means across retained splits, with no independent-sample confidence claim. Cb: OSD-561 cerebellum; Hip: OSD-562 hippocampus; RR-23: OSD-525 cerebellum. Baseline prediction uses a training-fitted intercept and age where available. Real-data spectral gain is the largest singular value of the oblique score operator on the aligned gene universe. A gain cap of two bounds amplification but can leave nuisance behind. This overview is optional artwork, not a required npj Microgravity display item or a claim of a universally validated brain mechanism.\n')
print('Saved optional graphical overview from real results.')
