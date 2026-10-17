from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Post-result biological illustration; every q comes from the full 3,624-test family."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=stage_root('extended')
plt.rcParams.update({'font.family':publication_font(),'font.size':10,'axes.unicode_minus':False,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
a=pd.read_csv(B/'tables/sorted_cell_panel_complete.csv')
s=a[(a.tissue=='Myeloid_Forebrain')&(a.contrast=='flight_balanced_age')]
q=s.pivot(index='analyte',columns=['QC','method'],values='q_full_3624_family')
effects=s.pivot(index='analyte',columns=['QC','method'],values='log2FC')
consensus=pd.DataFrame({'analyte':q.index,'all_six_model_QC_q_below_0_05':q.lt(.05).all(axis=1).to_numpy(),'minimum_log2FC':effects.min(axis=1).to_numpy(),'maximum_log2FC':effects.max(axis=1).to_numpy(),'maximum_q':q.max(axis=1).to_numpy()})
consensus.to_csv(B/'tables/sorted_forebrain_complete_panel_sensitivity_consensus.csv',index=False)
genes=['AIM2','NLRP3','PYCARD','ICE','GSDMD','IL1B','IL-18','NINJ1','AIF1','ALOX5','ALOX5AP']
order=[(qc,m) for qc in ['published_gate','alternative_gate','marker_positive_gate'] for m in ['limma_voom','edgeR_QL']]
values=effects.reindex(genes).loc[:,order].to_numpy();p=q.reindex(genes).loc[:,order].to_numpy()
fig=plt.figure(figsize=(15,10));gs=fig.add_gridspec(3,2,width_ratios=[1.25,1],left=.14,right=.96,top=.84,bottom=.13,hspace=.5,wspace=.4)
fig.text(.06,.96,'Innate-immune transcript remodeling in forebrain',fontsize=18,weight='bold',color='black')
fig.text(.06,.91,'RR8 myeloid-enriched forebrain | eight exposure donors | age-adjusted models and cell-selection sensitivities',fontsize=11,color='#444444')
ax=fig.add_subplot(gs[:,0]);im=ax.imshow(values,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
ax.set_yticks(range(len(genes)),['AIM2','NLRP3','PYCARD / ASC','CASP1 (panel: ICE)','GSDMD','IL1B','IL18 (panel: IL-18)','NINJ1','AIF1','ALOX5','ALOX5AP'])
ax.set_xticks(range(6),['Main\nvoom','Main\nedgeR','Alternative\nvoom','Alternative\nedgeR','Marker\nvoom','Marker\nedgeR'],fontsize=9)
for i in range(len(genes)):
    for j in range(6):
        ax.text(j,i,f'{values[i,j]:+.2f}'+('\n*' if p[i,j]<.05 else ''),ha='center',va='center',fontsize=9,color='white' if abs(values[i,j])>1 else '#111111')
ax.set_title('A  Point estimates; * = full-family q < 0.05',loc='left',fontsize=11,weight='bold',pad=12)
fig.colorbar(im,ax=ax,orientation='horizontal',pad=.10,fraction=.03,label='Age-balanced flight log2 fold change')
d=pd.read_csv(B/'tables/sorted_target_donor_expression.csv')
for i,g in enumerate(['PYCARD','AIM2','NINJ1']):
    ax=fig.add_subplot(gs[i,1]);v=d[(d.QC=='published_gate')&(d.tissue=='Myeloid_Forebrain')&(d.analyte==g)]
    for x,group in enumerate(['GC','FLT']):
        z=v[v.flight==group].sort_values(['age','donor'])
        for k,r in enumerate(z.itertuples()):ax.scatter(x+(k-1.5)*.065,r.voom_logCPM,c='#218A88' if x==0 else '#B95355',marker='o' if r.age==3 else '^',s=60)
        ax.plot([x-.17,x+.17],[z.voom_logCPM.mean()]*2,c='#444444',lw=1.6)
    ax.set_xticks([0,1],['Ground (4 mice)','Flight (4 mice)']);ax.set_ylabel('Donor pseudobulk log2 CPM');ax.set_title(f'{chr(66+i)}  {g}: individual animals',loc='left',weight='bold',fontsize=11)
fig.text(.06,.065,'Circles: 3 months; triangles: 8 months. Highlighted genes were selected after reviewing the full panel; the full 151-target table is retained.',fontsize=10,color='#444444')
fig.text(.06,.035,'Sensor, adaptor and effector transcripts do not establish protein activation. Processing day, enriched-cell mixtures and ambient RNA remain limitations.',fontsize=10,color='#444444')
for ext in ['png','svg']:fig.savefig(B/f'figures/figureS2_innate_immune_program.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
(B/'figures/figureS2_source_ledger.json').write_text(json.dumps({'tables':['sorted_cell_panel_complete.csv','sorted_target_donor_expression.csv'],'full_inference_family':3624,'selection':'Post-result biological illustration; no new p-values or independent replication from repeated methods/QC.','formats':['PNG','SVG']},indent=2))
