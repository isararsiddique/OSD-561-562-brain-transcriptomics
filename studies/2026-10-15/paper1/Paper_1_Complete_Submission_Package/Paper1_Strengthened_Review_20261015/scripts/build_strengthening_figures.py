#!/usr/bin/env python3
"""Publication panels from complete stored benchmarks; no outcome filtering."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
R=Path(__file__).resolve().parents[1];E=R/'extensions/strengthening_20261015';O=R/'publication_figures'
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':1.0,'lines.linewidth':1.2,'xtick.major.width':1,'ytick.major.width':1})
COL={'original_score':'#a14c60','score_gain2':'#d1954c','ridge_CV':'#587aab','gain2_ridge_CV':'#8b70a8','spectral_constrained_ridge_CV':'#187e72'}
LAB={'original_score':'Original score','score_gain2':'Bounded score','ridge_CV':'Ridge CV','gain2_ridge_CV':'Clipped ridge CV','spectral_constrained_ridge_CV':'Constrained CV'}
M=list(COL);acc=['OSD-561','OSD-562','OSD-525'];names=['RRRM-2 cerebellum','RRRM-2 hippocampus','RR-23 cerebellum']
old=pd.read_csv(E/'method_benchmark/outer_split_results.csv').rename(columns={'outside_prediction_MSE_ratio':'MSE_ratio'})
new=pd.read_csv(E/'method_benchmark/strengthened_outer_results.csv')
matched=pd.read_csv(E/'method_benchmark/matched_grid_outer_results.csv');matched=matched[matched.method.eq('matched_grid_gain2_clipped_ridge_CV')].copy();matched['method']='gain2_ridge_CV'
old=old[~old.method.eq('gain2_ridge_CV')];new=new[~new.method.eq('gain2_ridge_CV')];new=pd.concat([new,matched[matched.scope.eq('independent_exposure_cohort')]],ignore_index=True)
D=pd.concat([old[['accession','split_id','method','MSE_ratio','operator_gain']],new[['accession','split_id','method','MSE_ratio','operator_gain']],matched[matched.scope.eq('original_ground')][['accession','split_id','method','MSE_ratio','operator_gain']]],ignore_index=True)
def finish(fig,name):
 for ext in ['svg','png']:fig.savefig(O/(name+'.'+ext),dpi=600,facecolor='white')
 plt.close(fig)
def title(ax,letter,text):ax.set_title(letter+'  '+text,loc='left',fontweight='bold',pad=9)
rng=np.random.default_rng(71)
fig,axs=plt.subplots(2,3,figsize=(6.9,5.45),gridspec_kw={'height_ratios':[1.3,1]},layout='constrained')
for j,(a,n) in enumerate(zip(acc,names)):
 ax=axs[0,j];d=D[D.accession.eq(a)]
 for k,m in enumerate(M):
  v=d[d.method.eq(m)].MSE_ratio.to_numpy();ax.scatter(k+rng.uniform(-.16,.16,len(v)),v,s=4,color=COL[m],alpha=.24,rasterized=False);ax.plot([k-.27,k+.27],[v.mean()]*2,color=COL[m],lw=2.2)
 ax.axhline(1,color='#555',ls='--');ax.set_xticks(range(5),['Score','Bounded','Ridge','Clipped','Constrained'],rotation=55,ha='right');ax.set_ylim(.75,1.55);title(ax,'abc'[j],n)
 if j==0:ax.set_ylabel('Outside-sensor MSE / baseline')
 ax.text(.98,.96,f'{81 if j<2 else 200} overlapping splits',transform=ax.transAxes,ha='right',va='top',fontsize=8)
 ax=axs[1,j]
 for m in M:
  q=d[d.method.eq(m)];ax.scatter(q.operator_gain.mean(),q.MSE_ratio.mean(),s=45,color=COL[m],edgecolor='white',linewidth=1,zorder=3)
 ax.axhline(1,color='#555',ls='--');ax.axvline(2,color='#187e72',ls=':');ax.set_xlim(.5,35);ax.set_xlabel('Mean operator gain');title(ax,'def'[j],'Error versus gain')
 if j==0:ax.set_ylabel('Mean MSE ratio')
handles=[plt.Line2D([],[],marker='o',ls='',color=COL[m],label=LAB[m]) for m in M]
fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=8.5)
finish(fig,'figure11_constrained_prediction')
fig,axs=plt.subplots(2,2,figsize=(6.9,5.4),layout='constrained')
for ax,a,n,l in zip(axs.flat,['OSD-202','OSD-352','OSD-374','OSD-611'],['Earth radiation and unloading · 40 animals','RR-3 spaceflight · 6 animals','Earth protons · 8 pools','Earth Si28 · 10 primary pools'],'abcd'):
 d=new[new.accession.eq(a)]
 for k,m in enumerate(M):
  q=d[d.method.eq(m)];v=q.MSE_ratio.to_numpy();ax.scatter(k+rng.uniform(-.13,.13,len(v)),v,s=12,color=COL[m],alpha=.5);ax.scatter(k,v.mean(),marker='_',s=250,color=COL[m],linewidth=2.5,zorder=4);ax.scatter(k,q.prediction_SSE.sum()/q.baseline_SSE.sum(),marker='D',s=22,facecolor='white',edgecolor=COL[m],linewidth=1.2,zorder=5)
 ax.axhline(1,color='#555',ls='--');ax.set_xticks(range(5),['Score','Bounded','Ridge','Clipped','Constrained'],rotation=40,ha='right');ax.set_ylabel('Held-out MSE / baseline');title(ax,l,n)
 ax.set_ylim(max(0,float(d[d.method.isin(M)].MSE_ratio.min())-.06),max(1.05,float(d[d.method.isin(M)].MSE_ratio.max())+.06))
fig.legend(handles=[plt.Line2D([],[],marker='o',color='#444',ls='',label='One animal / pool'),plt.Line2D([],[],marker='_',markersize=15,color='#444',ls='',label='Equal-unit mean'),plt.Line2D([],[],marker='D',markerfacecolor='white',color='#444',ls='',label='Pooled SSE ratio')],loc='outside lower center',ncol=3,frameon=False,fontsize=8.5)
finish(fig,'figure12_external_prediction')
B=pd.read_csv(E/'correction_benchmark/benchmark_summary.csv');scens=list(json.loads((E/'correction_benchmark/protocol.json').read_text())['scenarios'])
sl=['Balanced','Correlated','Strongly\ncorrelated','Changed\nloading','Sensor\nbiology','Unequal\ndispersion','Biological\nnull']
ms=['raw','ground_score','ground_cap2','ground_cv_score','ground_ridge_CV','ground_constrained_ridge_CV','known_batch_OLS','neuroCombat_EB','neuroCombat_EB_batch_aware','neuroCombat_noEB','neuroCombat_noEB_batch_aware','RUVg_matrix_k1_oracle_controls','RUVg_covariate_k1_oracle_controls','PCA_covariate_k1']
ml=['Unadjusted','Ground score','Bounded score','Scalar CV','Ridge CV','Constrained CV','Known batch OLS','ComBat EB','ComBat EB + batch','ComBat no EB','ComBat no EB + batch','RUVg matrix*','RUVg covariate*','PCA covariate']
fig,axs=plt.subplots(2,2,figsize=(6.9,6.45),layout='constrained',sharey=True)
for row,n in enumerate([8,24]):
 t=B[B.n_train.eq(n)].set_index(['method','scenario']);raw=t.loc['raw'].reindex(scens).outside_effect_RMSE_mean
 rm=np.array([100*(t.loc[m].reindex(scens).outside_effect_RMSE_mean.to_numpy()/raw.to_numpy()-1) for m in ms]);fd=np.array([t.loc[m].reindex(scens).biological_FDP_mean.to_numpy() for m in ms])
 for col,(mat,cmap,norm) in enumerate([(rm,'RdBu_r',TwoSlopeNorm(vmin=-15,vcenter=0,vmax=110)),(fd,'magma_r',plt.Normalize(0,.8))]):
  ax=axs[row,col];im=ax.imshow(mat,aspect='auto',cmap=cmap,norm=norm);ax.set_xticks(range(7),sl,rotation=55,ha='right',fontsize=8);ax.set_yticks(range(14),ml,fontsize=8);ax.spines['left'].set_visible(False);ax.spines['bottom'].set_visible(False);ax.tick_params(length=0);title(ax,'abcd'[row*2+col],('Effect RMSE change' if col==0 else 'Mean FDP')+f' · n = {n}')
  for yy in range(14):
   for xx in range(7):
    value=mat[yy,xx];color='white' if (col==0 and value>60) or (col==1 and value>.4) else '#152e3b';ax.text(xx,yy,str(int(round(value))) if col==0 else f'{value:.2f}',ha='center',va='center',fontsize=7.8,color=color)
  ax.axhline(5.5,color='white',lw=2)
  if row==1:
   cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.06,pad=.03);cb.set_label('% RMSE change from unadjusted' if col==0 else 'Mean FDP at BH 0.05',fontsize=8)
finish(fig,'figure13_count_comparator_benchmark')
# Optional overview: information, constrained model, observed generalization.
fig=plt.figure(figsize=(6.9,4.6),facecolor='white');fig.text(.04,.96,'Stable prediction across spaceflight and Earth brain data',fontsize=13,fontweight='bold',va='top',color='#183b49')
fig.text(.04,.88,'a  Source data and experimental units',fontweight='bold',fontsize=10)
fig.text(.04,.81,'Original: 7 deposits, 83 animals\n362 overlapping ground holdouts\n\nExternal: 4 datasets, 64 held-out units\n46 individual animals + 18 pools',fontsize=9.5,va='top',linespacing=1.5)
fig.text(.55,.88,'b  Gain constrained regression',fontweight='bold',fontsize=10)
fig.text(.55,.81,'Predict outside-sensor expression\nRefit preprocessing within training folds\nSelect regularization by nested validation\n\n16,084 certified convex fits\nCorrection operator gain ≤ 2',fontsize=9.5,va='top',linespacing=1.5,color='#187e72')
ax=fig.add_axes([.13,.20,.80,.27]);q=new[new.method.eq('spectral_constrained_ridge_CV')].groupby('accession').MSE_ratio.mean().reindex(acc+['OSD-202','OSD-352','OSD-374','OSD-611']);ax.bar(np.arange(7),100*(1-q),color=['#587aab']*3+['#187e72']*4,width=.65);ax.set_xticks(range(7),['Cb','Hip','RR-23','Co57 + HLU','RR-3','Protons','Si28']);ax.set_ylabel('Mean MSE\nreduction (%)');title(ax,'c','Held-out prediction improves in all seven settings');ax.set_ylim(0,9.5)
for i,v in enumerate(100*(1-q)):ax.text(i,v+.3,f'{v:.1f}',ha='center',fontsize=9)
fig.text(.04,.06,'Bounded amplification and prediction gains are established here.\nBiological recovery and cross-context replication remain unproven.',fontsize=9.5,color='#354d58',linespacing=1.4)
finish(fig,'GRAPHICAL_ABSTRACT')
(O/'GRAPHICAL_ABSTRACT_CAPTION.md').write_text('Optional graphical overview. Separate original and external data units are preserved. Bars show descriptive reductions in mean held-out outside-sensor MSE relative to the training-fitted nuisance baseline; no independent-split confidence intervals are implied. External evaluation refits the algorithm within each cohort and does not transfer frozen fitted weights. Cb and Hip denote the original core cerebellum and hippocampus. Co57 + HLU denotes Earth cobalt-57 irradiation and hindlimb unloading. Constrained regression bounds the full correction operator gain at two and all 16,084 fits have relative convex optimality gap below 1e-7. Improved prediction is not evidence that recovered biological effects are causal.\n')
print('Created 3 main panel figures and updated graphical overview.')
