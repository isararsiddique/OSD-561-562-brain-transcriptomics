#!/usr/bin/env python3
"""Four-panel static scientific figure derived solely from complete saved tables."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from PIL import Image

OUT=Path(__file__).resolve().parent
font=Path('/System/Library/Fonts/Supplemental/Arial.ttf')
if not font.exists():raise FileNotFoundError('Arial required for publication audit figure')
font_manager.fontManager.addfont(str(font))
plt.rcParams.update({'font.family':'Arial','font.size':8.2,'axes.titlesize':9.3,'axes.labelsize':8.2,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':1,'lines.linewidth':1.2,'xtick.major.width':1,'ytick.major.width':1,'patch.linewidth':1,'svg.fonttype':'none'})
colors={'OSD-561':'#326C9A','OSD-562':'#277C61','OSD-525':'#A65F29'}
names={'OSD-561':'Cerebellum, RRRM-2','OSD-562':'Hippocampus, RRRM-2','OSD-525':'Cerebellum, RR-23'}
fig,axs=plt.subplots(2,2,figsize=(7.5,6.6),layout='constrained')
fig.get_layout_engine().set(w_pad=.075,h_pad=.10,wspace=.085,hspace=.095)
results=pd.read_csv(OUT/'exact_ground_sensor_association.csv')
for ax,n,letter in [(axs[0,0],'OSD-561','a'),(axs[0,1],'OSD-562','b')]:
    tab=pd.read_csv(OUT/f'{n}_exact_label_null.csv')
    row=results[(results.accession==n)&(results.panel_or_comparator=='Sensor31')].iloc[0]
    bins=np.linspace(-.025,1.025,22)
    ax.hist(tab.Sensor31,bins=bins,weights=np.full(len(tab),100/len(tab)),color=colors[n],alpha=.7,edgecolor='white',linewidth=1)
    ax.axvline(row.mean_crossfitted_age_stratified_AUC,color='#A32D2D',linewidth=1.5)
    ax.set_xlim(-.035,1.055);ax.set_ylim(0,20);ax.set_xticks([0,.25,.5,.75,1]);ax.set_xlabel('Mean held-out age-stratified AUC')
    ax.set_ylabel('Label assignments (%)')
    ax.set_title(f'{letter}  '+('Cerebellum' if n=='OSD-561' else 'Hippocampus')+': refitted null',loc='left',fontweight='bold')
    text=f'Observed Sensor31 AUC = {row.mean_crossfitted_age_stratified_AUC:.3f}\nExact p = {row.p_exact_complete_label_enumeration:.3f} (2 / 400)\n18-test BH q = {row.q_BH_complete_audit_18_tests:.5f}\nBY q = {row.q_BY_complete_audit_18_tests:.5f}'
    ax.text(.025,.98,text,transform=ax.transAxes,va='top',ha='left',fontsize=8,linespacing=1.55,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none','pad':3})
tertiles=pd.read_csv(OUT/'training_expression_tertile_prediction_errors.csv')
ax=axs[1,0]
for j,n in enumerate(['OSD-561','OSD-562','OSD-525']):
    means=[]
    for k in [1,2,3]:
        vals=tertiles[(tertiles.accession==n)&(tertiles.training_expression_tertile==k)].score_MSE_ratio.to_numpy()
        x=k+(j-1)*.13;mean=vals.mean();means.append(mean)
        ax.vlines(x,*np.quantile(vals,[.025,.975]),color=colors[n],linewidth=1.2)
        ax.plot(x,mean,'o',color=colors[n],markersize=4,label=names[n] if k==1 else None)
    ax.plot(np.arange(1,4)+(j-1)*.13,means,color=colors[n],linewidth=1,alpha=.55)
ax.axhline(1,color='#777777',linestyle='--',linewidth=1)
ax.set_xticks([1,2,3],['Lower','Middle','Higher']);ax.set_xlabel('Training mean-count tertile')
ax.set_ylabel('Outside-sensor prediction MSE / baseline')
ax.set_ylim(.47,1.52);ax.set_xlim(.65,3.35)
ax.set_title('c  Prediction risk across expression levels',loc='left',fontweight='bold')
ax.legend(frameon=False,fontsize=7,loc='lower left',handletextpad=.45)
ax.text(.025,.98,'Means and overlapping-split 2.5–97.5% ranges',transform=ax.transAxes,va='top',fontsize=7.2)
saturated=pd.read_csv(OUT/'saturated_age_panel_interactions.csv')
ax=axs[1,1]
labels=['12 weeks','29 weeks','Equal-age mean'];estimands=['young','old','equal_age_mean']
for j,n in enumerate(['OSD-561','OSD-562']):
    sub=saturated[(saturated.accession==n)&(saturated.panel=='Circadian')].set_index('estimand').loc[estimands]
    xs=np.arange(3)+(j-.5)*.16
    ax.vlines(xs,sub.HC3_normal_lower95,sub.HC3_normal_upper95,color=colors[n],linewidth=1.2)
    ax.plot(xs,sub.effect_log2,'o',color=colors[n],markersize=4,label='Cerebellum' if n=='OSD-561' else 'Hippocampus')
ax.axhline(0,color='#777777',linestyle='--',linewidth=1)
ax.set_xticks(np.arange(3),labels);ax.set_xlim(-.4,2.4);ax.set_ylim(-1.33,.71)
ax.set_ylabel('Circadian interaction (log2)')
ax.set_title('d  Circadian interaction by age',loc='left',fontweight='bold')
ax.text(.025,.98,'Full-data model; HC3 normal 95% intervals\nNo BH q ≤ 0.05 across 64 estimands',transform=ax.transAxes,va='top',fontsize=7.3,linespacing=1.45)
ax.legend(frameon=False,fontsize=7.3,loc='lower left')
for ax in axs.flat:
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',linewidth=1,alpha=.12)
fig.savefig(OUT/'independent_audit_panel.png',dpi=600,facecolor='white')
fig.savefig(OUT/'independent_audit_panel.svg',facecolor='white')
plt.close(fig)
img=Image.open(OUT/'independent_audit_panel.png');svg=(OUT/'independent_audit_panel.svg').read_text()
checks={'figure_inches':[7.5,6.6],'PNG_pixels':list(img.size),'PNG_dpi':list(img.info['dpi']),'Arial_present_in_SVG':'Arial' in svg,'SVG_has_editable_text':'<text ' in svg,'SVG_has_no_embedded_raster':'<image' not in svg,'panel_mapping':{'a':'Cerebellum Sensor31 exact refitted-label null','b':'Hippocampus Sensor31 exact refitted-label null','c':'Training-expression-tertile prediction MSE ratios across three settings','d':'Full-data Circadian age-specific/equal-age saturated effects, HC3 normal intervals'}}
(OUT/'figure_specification.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
