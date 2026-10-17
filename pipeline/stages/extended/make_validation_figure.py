from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Validation views from complete mouse-level analysis families."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=stage_root('extended')
plt.rcParams.update({'font.family':publication_font(),'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'svg.fonttype':'none','pdf.fonttype':42})
loo=pd.read_csv(B/'tables/sorted_leave_one_mouse_effects.csv');full=pd.read_csv(B/'tables/sorted_cell_panel_complete.csv');ratios=pd.read_csv(B/'tables/sorted_marker_relative_tests.csv');depth=pd.read_csv(B/'tables/sorted_cell_depth_tests.csv')
genes=['AIF1','PYCARD','AIM2','NINJ1'];region='Myeloid_Forebrain';blue='#235A84';red='#B95355';teal='#218A88';gray='#444444'
fig,axes=plt.subplots(1,3,figsize=(15,6.5));fig.suptitle('Forebrain validation across mice, reference transcripts and RNA depth',x=.06,y=.98,ha='left',fontsize=16,weight='bold',color='black')
fig.text(.06,.885,'Eight independent mice; age-balanced flight comparisons. All results come from complete 151-target analysis families.',fontsize=10,color=gray)
for gi,g in enumerate(genes):
    for mi,m in enumerate(['limma_voom','edgeR_QL']):
        a=loo[(loo.tissue==region)&(loo.analyte==g)&(loo.method==m)].log2FC.dropna();r=full[(full.tissue==region)&(full.analyte==g)&(full.QC=='published_gate')&(full.method==m)&(full.contrast=='flight_balanced_age')].iloc[0]
        axes[0].errorbar(r.log2FC,gi+(mi-.5)*.16,xerr=[[r.log2FC-a.min()],[a.max()-r.log2FC]],fmt='o',color=[blue,red][mi],capsize=3,markersize=5,label=['limma/voom','edgeR QL'][mi] if gi==0 else None)
    r=ratios[(ratios.tissue==region)&(ratios.analyte==g)&(ratios.contrast=='flight_balanced_age')&(ratios.error_model=='HC3_t')].iloc[0]
    axes[1].errorbar(r.effect,gi,xerr=[[r.effect-r.ci_low],[r.ci_high-r.effect]],fmt='o',mfc=teal if r.q_full_1812_family<.05 else 'white',mec=teal,ecolor=teal,capsize=3,markersize=6);axes[1].text(1.08,gi,f'q={r.q_full_1812_family:.3f}',va='center',fontsize=8,color=gray)
    r=depth[(depth.tissue==region)&(depth.analyte==g)&(depth.endpoint=='expected_detection_at_500_UMIs')&(depth.error_model=='HC3_t')].iloc[0]
    axes[2].errorbar(100*r.effect,gi,xerr=[[100*(r.effect-r.ci_low)],[100*(r.ci_high-r.effect)]],fmt='o',mfc=teal if r.q_full_1359_family<.05 else 'white',mec=teal,ecolor=teal,capsize=3,markersize=6);axes[2].text(6.2,gi,f'q={r.q_full_1359_family:.3f}',va='center',fontsize=8,color=gray)
for ax in axes:
    ax.set_yticks(range(4),genes);ax.set_ylim(3.5,-.6);ax.axvline(0,c=gray,lw=.8,ls='--')
axes[0].set_title('A  Leave-one-mouse-out refits',loc='left',fontsize=11,weight='bold');axes[0].set_xlabel('Flight log2 fold change');axes[0].legend(loc='upper right',fontsize=8,frameon=False);axes[0].set_xlim(-1.65,1.05)
axes[1].set_title('B  Relative to myeloid transcripts',loc='left',fontsize=11,weight='bold');axes[1].set_xlabel('Flight effect in log2 target/reference ratio');axes[1].set_xlim(-1.75,1.8)
axes[2].set_title('C  Expected detection at 500 UMIs',loc='left',fontsize=11,weight='bold');axes[2].set_xlabel('Flight difference in percentage points');axes[2].set_xlim(-20,12.5)
fig.subplots_adjust(left=.10,right=.97,top=.77,bottom=.25,wspace=.38)
fig.text(.06,.15,'A: bars span eight deletion estimates; they are not confidence intervals. B/C: 95% HC3 model intervals; filled points indicate corrected q < 0.05.',fontsize=9,color=gray)
fig.text(.06,.105,'B: reference = geometric mean of Csf1r, Tmem119, P2ry12, Hexb and Sall1 (biological denominators, not assumed stable housekeeping genes).',fontsize=9,color=gray)
fig.text(.06,.06,'C: analytical rarefaction followed by mouse-level inference. These checks do not remove processing-day confounding or establish protein activity.',fontsize=9,color=gray)
for ext in ['png','pdf','svg']:fig.savefig(B/f'figures/figure7_mouse_and_depth_validation.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
(B/'figures/figure7_source_ledger.json').write_text(json.dumps({'tables':['sorted_leave_one_mouse_effects.csv','sorted_cell_panel_complete.csv','sorted_marker_relative_tests.csv','sorted_cell_depth_tests.csv'],'selection':'Four forebrain genes displayed after full-panel tests; no reduced four-gene correction','reference_family':1812,'depth_family':1359,'font':'Arial','scope':'Deletion stability, relative RNA and standardized detection on the same eight donors; no additional biological replication.'},indent=2))
