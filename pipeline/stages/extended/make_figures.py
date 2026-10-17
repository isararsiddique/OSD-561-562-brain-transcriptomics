from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Integrated illustrated workflow and publication-sized evidence figures."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle,FancyArrowPatch
from scipy.stats import spearmanr
from PIL import Image
B=stage_root('extended');V1=stage_root('core');F=B/'figures'
plt.rcParams.update({'font.family':publication_font(),'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
blue='#235A84';teal='#218A88';red='#B95355';gray='#444444';colors=[blue,teal,red]
ledger=[]
def save(fig,name,sources):
    fig.savefig(F/(name+'.png'),dpi=210,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    ledger.append({'figure':name,'tables':sources,'authoring_script':__file__})
def title(fig,text,subtitle):
    fig.suptitle(text,x=.055,ha='left',y=.99,fontsize=16,fontweight='bold',color='black')
    fig.text(.055,.99-.50/fig.get_figheight(),subtitle,color=gray,fontsize=10)
import importlib.util
workflow_spec=importlib.util.spec_from_file_location('deep_workflow',code_root()/'stages/deep/workflow_figure.py')
workflow_module=importlib.util.module_from_spec(workflow_spec);workflow_spec.loader.exec_module(workflow_module)
workflow_sources=workflow_module.draw_workflow()
ledger.append({'figure':'figure1_illustrated_study_workflow','tables':workflow_sources,'authoring_script':str(code_root()/'stages/deep/workflow_figure.py')})

# Method-dependent strength, shown alongside whole-genome effect agreement.
old=pd.read_csv(V1/'tables/bulk_panel_complete.csv');new=pd.read_csv(B/'tables/bulk_method_comparison_panel.csv')
fig,axs=plt.subplots(1,3,figsize=(15,5.8),gridspec_kw={'width_ratios':[1.10,1.10,1.2]});title(fig,'Bulk expression effects across statistical methods','Same original samples; every method uses the full 4,228-target contrast family.')
methods=['DESeq2','edgeR QL','limma/voom']
method_counts=[];displayed_effects=[]
for name,source,qcol in [('DESeq2',old,'q_panel_primary_full_family'),('edgeR QL',new[new.method.eq('edgeR_QL_TMM')],'q_full_family'),('limma/voom',new[new.method.eq('limma_voom_TMM')],'q_full_family')]:
    assert len(source)==4228 and np.isfinite(source[qcol]).all()
    method_counts.append({'method':name,'planned_family_size':len(source),'missing_adjusted_p_values':int(source[qcol].isna().sum()),'corrected_contrasts':int(source[qcol].lt(.05).sum())})
pd.DataFrame(method_counts).to_csv(B/'tables/figure2_method_counts.csv',index=False)
for j,g in enumerate(['AIF1','PYCARD']):
    r=old[(old.study=='OSD-561')&(old.contrast=='flight_ISS_T_age12')&(old.analyte==g)].iloc[0]
    vals=[(r.log2FoldChange,r.q_panel_primary_full_family)]
    for method in ['edgeR_QL_TMM','limma_voom_TMM']:
        r=new[(new.study=='OSD-561')&(new.contrast=='flight_ISS_T_age12')&(new.analyte==g)&(new.method==method)].iloc[0];vals.append((r.log2FC,r.q_full_family))
    for i,(effect,q) in enumerate(vals):
        assert np.isfinite(effect) and np.isfinite(q)
        y=j*3+i;axs[0].scatter(effect,y,facecolors=colors[i] if q<.05 else 'white',edgecolors=colors[i],linewidths=1.5,s=65);axs[0].text(effect+.07,y,f'q={q:.3f}',va='center',fontsize=9)
        displayed_effects.append({'analyte':g,'method':methods[i],'log2_fold_change':effect,'q_full_family':q})
pd.DataFrame(displayed_effects).to_csv(B/'tables/figure2_displayed_effects.csv',index=False)
axs[0].set_yticks(range(6),[f'{g} | {m}' for g in ['AIF1','PYCARD'] for m in methods]);axs[0].set_xlim(-2,-.6);axs[0].set_ylim(5.5,-.5);axs[0].set_xlabel('Cerebellum log2 fold change');axs[0].set_title('A  Young terminal flight contrast',loc='left',fontsize=10.5,weight='bold')
ranked_bulk=pd.read_csv(B/'tables/bulk_ranked_evidence.csv')
for i,method in enumerate(methods):
    s=ranked_bulk[ranked_bulk.method.eq(method)].sort_values('rank');assert len(s)==4228
    axs[1].plot(s['rank'],s.negative_log10_q,c=colors[i],ls=['-','--',':'][i],lw=1.8,label=method)
axs[1].axhline(-np.log10(.05),color=gray,ls='--',lw=.8);axs[1].text(.97,-np.log10(.05)+.10,'q = 0.05',transform=axs[1].get_yaxis_transform(),ha='right',fontsize=8)
axs[1].set(xscale='log',xlim=(1,4228),ylim=(-.08,5.95),xlabel='Rank in full family (log scale)',ylabel='-log10 adjusted p-value (q)');axs[1].legend(loc='upper right',fontsize=8,frameon=False);axs[1].set_title('B  Evidence across all targets',loc='left',fontsize=10.5,weight='bold')
ds=pd.read_csv(V1/'tables/OSD-561_STAR_flight_ISS_T_age12_genome.csv.gz',index_col=0);vo=pd.read_csv(B/'tables/OSD-561_voom_flight_ISS_T_age12_genome.csv.gz',index_col=0);z=ds[['log2FoldChange']].join(vo[['logFC']]).dropna()
hb=axs[2].hexbin(z.log2FoldChange,z.logFC,gridsize=55,mincnt=1,cmap='Blues',bins='log');axs[2].plot([-4,4],[-4,4],ls='--',c=red,lw=1);axs[2].set(xlim=(-4,4),ylim=(-4,4),xlabel='DESeq2 log2 fold change',ylabel='limma/voom log2 fold change');axs[2].set_title(f'C  Whole-genome effects | Spearman {spearmanr(z.log2FoldChange,z.logFC).statistic:.2f}',loc='left',fontsize=10,weight='bold');fig.colorbar(hb,ax=axs[2],fraction=.05,pad=.025,label='Genes per bin (log scale)')
inside=z.log2FoldChange.between(-4,4)&z.logFC.between(-4,4)
(B/'tables/figure2_genome_display_summary.json').write_text(json.dumps({'paired_finite_gene_effects':len(z),'inside_display_window':int(inside.sum()),'outside_display_window':int((~inside).sum()),'display_axis_limits':[-4,4],'white_bin_meaning':'No observations in the bin; not a missing gene-effect estimate.'},indent=2))
fig.text(.055,.095,'A: filled circles indicate q < 0.05. B: all 4,228 tests per method; higher curves indicate stronger evidence. Baseline is q=1.',fontsize=9,color=gray)
fig.text(.055,.045,f'C: {len(z):,} shared finite gene effects; white bins contain no observations. Axes show log2 effects from -4 to +4.',fontsize=9,color=gray)
fig.subplots_adjust(top=.80,bottom=.24,wspace=.65,left=.13,right=.965);save(fig,'figure2_method_agreement',['bulk_method_comparison_panel.csv','../core/tables/bulk_panel_complete.csv','../core/tables/OSD-561_STAR_flight_ISS_T_age12_genome.csv.gz','OSD-561_voom_flight_ISS_T_age12_genome.csv.gz','figure2_method_counts.csv','figure2_displayed_effects.csv','figure2_genome_display_summary.json','bulk_ranked_evidence.csv'])

# Reference mixture performance and competing explanation.
cv=pd.read_csv(B/'tables/deconvolution_held_donor_validation.csv');E=pd.read_csv(B/'tables/bulk_RNA_fraction_estimates.csv')
fig,axs=plt.subplots(1,3,figsize=(15,5.5));title(fig,'Estimated immune RNA and conditional gene effects','Reference estimates describe RNA fractions; they do not measure absolute cell abundance or establish mediation.')
s=cv[(cv.method=='NNLS')&(cv.cell_type=='Immune cells')];axs[0].scatter(s.true_RNA_fraction,s.estimated_RNA_fraction,s=10,c=teal,alpha=.25);axs[0].plot([0,1],[0,1],ls='--',c=gray);axs[0].set(xlabel='Known synthetic immune RNA fraction',ylabel='Estimated immune RNA fraction',xlim=(0,1),ylim=(0,1));axs[0].set_title('A  Five held-out donor benchmarks',loc='left',weight='bold',fontsize=11)
s=E[(E.study=='OSD-561')&(E.method=='NNLS')&(E.cell_type=='Immune cells')&E.group.isin(['FLT_ISS_T_12','GC_ISS_T_12'])]
for i,f in enumerate(['GC','FLT']):
    v=s[s.flight==f].estimated_RNA_fraction.to_numpy()*100;axs[1].scatter(i+np.linspace(-.07,.07,len(v)),v,c=teal if i==0 else red,s=55);axs[1].plot([i-.15,i+.15],[v.mean()]*2,c=gray,lw=2)
axs[1].set_xticks([0,1],['Ground (n=3)','Flight (n=5)']);axs[1].set_ylabel('Estimated immune RNA fraction (%)');axs[1].set_title('B  Young terminal cerebellum',loc='left',weight='bold',fontsize=11)
ad=pd.read_csv(B/'tables/composition_adjusted_young_terminal_panel.csv');s=ad[(ad.study=='OSD-561')&(ad.method=='NNLS')&ad.analyte.isin(['AIF1','PYCARD'])]
for i,r in enumerate(s.itertuples()):
    axs[2].scatter(r.flight_effect_unadjusted,i-.12,c=blue,s=55,label='CPM baseline' if i==0 else None);axs[2].errorbar(r.flight_effect_adjusted,i+.12,xerr=[[r.flight_effect_adjusted-r.ci_low],[r.ci_high-r.flight_effect_adjusted]],fmt='o',c=red,label='RNA-covariate OLS' if i==0 else None)
axs[2].set_yticks(range(2),s.analyte);axs[2].set_ylim(-.25,1.5);axs[2].axvline(0,c=gray,ls='--');axs[2].set_xlabel('Flight effect in log2(CPM+1)');axs[2].legend(fontsize=8,loc='lower right');axs[2].set_title('C  Conditional associations\nFlight-immune RNA correlation: -0.94',loc='left',weight='bold',fontsize=10)
fig.subplots_adjust(top=.80,bottom=.17,wspace=.42);save(fig,'figure3_composition_and_identifiability',['deconvolution_held_donor_validation.csv','bulk_RNA_fraction_estimates.csv','composition_adjusted_young_terminal_panel.csv'])

# Matched regulatory context, complete effects and individual donor values.
a=pd.read_csv(B/'tables/RNA_ATAC_promoter_context.csv');d=pd.read_csv(B/'tables/RNA_ATAC_donor_values.csv')
fig,axs=plt.subplots(2,3,figsize=(16,10));title(fig,'Myeloid RNA and matched promoter accessibility','Separate missions: RR3 RNA/ATAC context above; RR8 sorted-cell donor effects below. Filled points: full-family q < 0.05.')
for i,ct in enumerate(['CNS neurons','CNS glia','Immune cells']):
    ax=axs[0,i];s=a[a.cell_type.eq(ct)].dropna(subset=['RNA_log2CPM_effect','ATAC_log2CPM_effect']);ax.scatter(s.RNA_log2CPM_effect,s.ATAC_log2CPM_effect,c=gray,alpha=.25,s=15)
    for g in ['Aif1','Pycard']:
        r=s[s.mouse_symbol.eq(g)].iloc[0];ax.scatter(r.RNA_log2CPM_effect,r.ATAC_log2CPM_effect,c=red if g=='Aif1' else teal,s=70);ax.annotate(g,(r.RNA_log2CPM_effect,r.ATAC_log2CPM_effect),xytext=(5,5),textcoords='offset points',fontsize=10)
    ax.axvline(0,c=gray,lw=.7);ax.axhline(0,c=gray,lw=.7);ax.set(xlabel='RNA flight effect: log2(CPM+1)',ylabel='Promoter ATAC flight effect: log2(CPM+1)');ax.set_title('RR3 | '+ct+' | 5 donors',loc='left',color='black',weight='bold',fontsize=11);ax.text(.03,.96,'0 / 453 corrected promoter tests',transform=ax.transAxes,va='top',fontsize=8)
sorted_panel=pd.read_csv(B/'tables/sorted_cell_panel_complete.csv')
for i,tissue in enumerate(['Myeloid_Cerebellum','Myeloid_Forebrain','Myeloid_Hippocampus']):
    ax=axs[1,i];s=sorted_panel[(sorted_panel.tissue==tissue)&(sorted_panel.contrast=='flight_balanced_age')]
    for gi,gene in enumerate(['AIF1','PYCARD']):
        for qi,qc in enumerate(['published_gate','alternative_gate','marker_positive_gate']):
            for mi,method in enumerate(['limma_voom','edgeR_QL']):
                r=s[(s.analyte==gene)&(s.QC==qc)&(s.method==method)].iloc[0];c=[blue,red][mi]
                ax.scatter(r.log2FC,gi*3+qi+(mi-.5)*.12,s=48,marker=['o','s','^'][qi],facecolors=c if r.q_full_3624_family<.05 else 'white',edgecolors=c,lw=1.3,label=['limma/voom','edgeR QL'][mi] if gi==0 and qi==0 else None)
    ax.set_yticks(range(6),[f'{g} | {q}' for g in ['AIF1','PYCARD'] for q in ['main QC','alternative QC','marker gate']],fontsize=8);ax.invert_yaxis();ax.axvline(0,c=gray,ls='--',lw=.8);ax.set_xlim(-1.6,.1);ax.set_xlabel('Age-balanced flight log2 fold change');ax.set_title('RR8 | '+tissue.replace('Myeloid_','')+f' | n={7 if i==0 else 8} donors',loc='left',color='black',weight='bold',fontsize=11)
    if i==1:ax.legend(frameon=False,loc='lower right',fontsize=8)
fig.text(.055,.035,'RR8: myeloid-enriched populations; 3- and 8-month age strata. Marker gate excludes both lead genes. Processing day remains confounded with flight.',fontsize=10,color=gray)
fig.subplots_adjust(top=.85,bottom=.11,wspace=.45,hspace=.53,left=.11);save(fig,'figure4_RNA_ATAC_context',['RNA_ATAC_promoter_context.csv','RNA_ATAC_donor_values.csv','panel_promoter_peak_links.csv','sorted_cell_panel_complete.csv','sorted_marker_gate_definition.json'])

# Raw spatial layouts, with complete donor identity and no invented region names.
sp=pd.read_csv(B/'prepared/Visium_lead_gene_spots.csv.gz');reg=pd.read_csv(B/'tables/Visium_section_registry.csv');dv=pd.read_csv(B/'tables/Visium_donor_mean_scores.csv')
fig=plt.figure(figsize=(15,9.5));gs=fig.add_gridspec(3,4,top=.85,bottom=.16,wspace=.35,hspace=.42);title(fig,'Spatial gene localization and donor flight effects','Verified coordinates: 29,839 spots, 12 sections, six mice. Displayed sections are the first technical section of each donor.')
for i,(_,r) in enumerate(reg.sort_values('author_section_label').drop_duplicates('donor').iterrows()):
    ax=fig.add_subplot(gs[i//3,i%3]);s=sp[sp.section.eq(r.section)];ax.scatter(s.x,s.y,c='#E1E6EB',s=2);v=s[s.Aif1_raw>0];spatial_colors=ax.scatter(v.x,v.y,c=v.Aif1_log2CPM,cmap='YlOrRd',s=5,vmin=0,vmax=9);ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off');ax.set_title(f'{r.author_section_label} | {r.donor}\nAif1 | {r.flight}',fontsize=9,color='black')
for k,g in enumerate(['AIF1','PYCARD']):
    ax=fig.add_subplot(gs[k,3]);s=dv[dv.analyte.eq(g)]
    for i,f in enumerate(['Ground Control','Space Flight']):
        v=s[s.flight.eq(f)].mean_spatial_score.to_numpy();ax.scatter(i+np.linspace(-.07,.07,3),v,c=teal if i==0 else red,s=50)
    ax.set_xticks([0,1],['Ground','Flight']);ax.set_ylabel('Mean log2(CPM+1)');ax.set_title(f'{g}: each point is one animal',fontsize=10,loc='left')
ax=fig.add_subplot(gs[2,:]);ax.axis('off');s=pd.read_csv(B/'tables/Visium_spatial_structure_tests.csv').query('mouse_symbol in ["Aif1","Pycard"]')
for i,g in enumerate(['Aif1','Pycard']):
    z=s[s.mouse_symbol.eq(g)];ax.text(.02,.8-i*.3,f'{g}: {int(z.q_spatial_structure_family.lt(.05).sum())}/12 sections show corrected spatial structure; sparse sections retain p=1.',fontsize=12,color='black')
ax.text(.02,.15,'Neither lead has a corrected six-animal flight effect. Spatial clustering does not establish a spaceflight-induced mechanism.',fontsize=11,color=gray)
cbax=fig.add_axes([.085,.075,.24,.014]);fig.colorbar(spatial_colors,cax=cbax,orientation='horizontal',label='Aif1: log2(CPM+1), shared scale')
save(fig,'figure5_Visium_localization',['Visium_section_registry.csv','Visium_spatial_structure_tests.csv','Visium_donor_mean_scores.csv','Visium_donor_flight_tests.csv','prepared/Visium_lead_gene_spots.csv.gz'])

# Continuous pathway evidence and independent animal-level score displays.
ca=pd.read_csv(B/'tables/pathway_CAMERA_complete.csv');ex=pd.read_csv(B/'tables/pathway_animal_exact_tests.csv');sorted_ca=pd.read_csv(B/'tables/sorted_pathway_CAMERA_complete.csv')
models=['default_0.01','conservative_0.10','estimated'];model_labels=['Fixed correlation = 0.01','Fixed correlation = 0.10','Estimated correlation']
model_colors=[blue,'#A86120',teal];line_styles=['-','--',':']
families=[('Core bulk','core_50_sets_28_contrasts',1400),('Independent bulk','independent_bulk_50_sets_3_studies',150),('RR3 cell types','celltypes_50_sets_5_types',250),('RR8 enriched populations',None,1200)]
count_rows=[]
for name,fam,planned in families:
    for model in models:
        source=ca[(ca.family==fam)&(ca.correlation_model==model)] if fam else sorted_ca[sorted_ca.correlation_model==model]
        qcol='q_full_family' if fam else 'q_full_1200_family';assert len(source)==planned and np.isfinite(source[qcol]).all()
        count_rows.append({'branch':name,'correlation_model':model,'planned_family_size':planned,'missing_adjusted_p_values':int(source[qcol].isna().sum()),'corrected_contrasts':int(source[qcol].lt(.05).sum())})
pd.DataFrame(count_rows).to_csv(B/'tables/figure6_CAMERA_count_grid.csv',index=False)
z=ca[(ca['set']=='HALLMARK_INTERFERON_ALPHA_RESPONSE')&ca.study.isin(['OSD-352','OSD-563'])]
z[['study','correlation_model','q_full_family']].assign(negative_log10_q=-np.log10(z.q_full_family)).to_csv(B/'tables/figure6_example_q_values.csv',index=False)
ranked=pd.read_csv(B/'tables/pathway_ranked_evidence.csv');correlations=pd.read_csv(B/'tables/pathway_unique_correlation_estimates.csv')
donor_scores=pd.read_csv(B/'tables/figure6_interferon_donor_scores.csv');animal_tests=pd.read_csv(B/'tables/figure6_interferon_animal_tests.csv')
fig,axs=plt.subplots(2,3,figsize=(16,10));title(fig,'Pathway coherence, animal-level scores and statistical uncertainty','Complete evidence distributions under three CAMERA assumptions, with correlations and individual-mouse RNA measurements.')
ax=axs[0,0]
for i,(name,_,_) in enumerate(families):
    v=correlations[correlations.branch.eq(name)&correlations.displayed_adequate_coverage].estimated_correlation.to_numpy()
    ax.scatter(v,i+.085*np.sin(np.arange(len(v))*2.4),s=7,color=teal,alpha=.25,zorder=2)
    quartiles=np.quantile(v,[.10,.25,.5,.75,.90]);ax.plot(quartiles[[0,4]],[i,i],color='black',lw=1.1);ax.plot(quartiles[[1,3]],[i,i],color='black',lw=5);ax.scatter(quartiles[2],i,facecolor='white',edgecolor='black',s=22,zorder=4)
ax.axvline(.01,c=model_colors[0],ls='-',lw=1);ax.axvline(.10,c=model_colors[1],ls='--',lw=1)
ax.set_yticks(range(4),['Core bulk','Independent bulk','RR3 cell types','RR8 enriched\npopulations'],fontsize=9);ax.set(xlim=(-.23,.50),ylim=(3.5,-.55),xlabel='Estimated inter-gene correlation');ax.set_title('A  Within-set correlation',loc='left',fontsize=11,weight='bold')
for letter,ax,(branch,_,planned) in zip(['B','C','E','F'],[axs[0,1],axs[0,2],axs[1,1],axs[1,2]],families):
    s=ranked[ranked.branch.eq(branch)]
    for mi,model in enumerate(models):
        v=s[s.correlation_model.eq(model)].sort_values('rank');assert len(v)==planned
        ax.plot(v['rank'],v.negative_log10_q,c=model_colors[mi],ls=line_styles[mi],lw=1.6,label=model_labels[mi])
    threshold=-np.log10(.05);ax.axhline(threshold,c=gray,ls='--',lw=.7);ax.text(.97,threshold+.12,'q = 0.05',transform=ax.get_yaxis_transform(),ha='right',fontsize=8)
    ax.set(xscale='log',xlim=(1,planned),ylim=(-.18,max(2,s.negative_log10_q.max()+.65)),xlabel='Rank in full family (log scale)',ylabel='-log10 adjusted p-value (q)');ax.tick_params(labelsize=9);ax.set_title(f'{letter}  {branch} | {planned:,} tests',loc='left',fontsize=10.5,weight='bold')
ax=axs[1,0]
for si,sid in enumerate(['OSD-352','OSD-563']):
    s=donor_scores[donor_scores.study.eq(sid)]
    for gi,group in enumerate(['GC','FLT']):
        x=si*3+gi;v=s[s.flight.eq(group)].sort_values('sample').positive_set_singscore.to_numpy();ax.scatter(x+np.linspace(-.12,.12,len(v)),v,c=teal if gi==0 else red,s=35,zorder=3);ax.plot([x-.2,x+.2],[v.mean()]*2,c=gray,lw=1.4)
ax.set_xticks([0,1,3,4],['RR3\nGround\nn=3','RR3\nFlight\nn=3','RR10\nGround\nn=5','RR10\nFlight\nn=5'],fontsize=8);ax.set(xlim=(-.5,4.5),ylim=(-.085,.012),ylabel='Interferon-alpha singscore');ax.text(.03,.96,f'Both full-family exact q={animal_tests.q_full_family.iloc[0]:.3f}',transform=ax.transAxes,fontsize=8,va='top');ax.set_title('D  Individual-animal pathway scores',loc='left',fontsize=10.5,weight='bold')
handles,labels=axs[0,1].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.59,.905),ncol=3,frameon=False,fontsize=10)
fig.text(.055,.095,'A: one index exposure contrast per study/pathway; bars show 10th-90th percentiles, thick bars the interquartile range, white points medians.',fontsize=9,color=gray)
fig.text(.055,.060,'B/C/E/F: every planned test, including unavailable hypotheses at q=1. Y-axis scales differ. D: points are mice; horizontal bars are group means.',fontsize=9,color=gray)
fig.text(.055,.025,'Estimated-correlation CAMERA also uses only 3-19 residual degrees of freedom. No conservative/estimated pathway result passes full-family correction.',fontsize=9,color=gray)
fig.subplots_adjust(top=.83,bottom=.19,left=.11,right=.97,wspace=.43,hspace=.47)
save(fig,'figure6_pathway_correlation_sensitivity',['pathway_CAMERA_complete.csv','pathway_animal_exact_tests.csv','sorted_pathway_CAMERA_complete.csv','sorted_pathway_exact_tests.csv','figure6_CAMERA_count_grid.csv','figure6_example_q_values.csv','pathway_ranked_evidence.csv','pathway_unique_correlation_estimates.csv','figure6_interferon_donor_scores.csv','figure6_interferon_animal_tests.csv','pathway_camera_degrees_of_freedom.csv'])

# Full independent panel, all target rows and missingness retained.
val=pd.read_csv(B/'tables/independent_bulk_panel.csv');rr9=pd.read_csv(B/'tables/OSD536_panel_results.csv');map=pd.read_csv(V1/'tables/verified_panel_mapping.csv')
cols=[];labels=[];qcols=[]
for sid in ['OSD-352','OSD-563','OSD-564']:
    for method in ['edgeR_QL_TMM','limma_voom_TMM']:
        s=val[(val.study==sid)&(val.method==method)].set_index('analyte').reindex(map.analyte);cols.append(s.log2FC.to_numpy());qcols.append(s.q_full_family.to_numpy());labels.append(sid+'\n'+('edgeR' if method.startswith('edgeR') else 'voom'))
s=rr9.set_index('analyte').reindex(map.analyte);cols.append(s.log2_effect.to_numpy());qcols.append(s.q_full_151_exact.to_numpy());labels.append('OSD-536\nnCounter')
for tissue in ['Myeloid_Cerebellum','Myeloid_Forebrain','Myeloid_Hippocampus','Non-myeloid_Cerebellum','Non-myeloid_Forebrain','Non-myeloid_Hippocampus']:
    for method in ['edgeR_QL','limma_voom']:
        s=sorted_panel[(sorted_panel.tissue==tissue)&(sorted_panel.method==method)&(sorted_panel.QC=='published_gate')&(sorted_panel.contrast=='flight_balanced_age')].set_index('analyte').reindex(map.analyte);cols.append(s.log2FC.to_numpy());qcols.append(s.q_full_3624_family.to_numpy());labels.append(tissue.replace('Non-myeloid_','Non-M: ').replace('Myeloid_','M: ')+'\n'+('edgeR' if method.startswith('edgeR') else 'voom'))
mat=np.array(cols).T
fig,ax=plt.subplots(figsize=(16,24));im=ax.imshow(np.ma.masked_invalid(mat),aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set_xticks(range(len(cols)),labels,rotation=90,fontsize=8);ax.set_yticks(range(151),map.analyte,fontsize=7);ax.set_title('Independent mission context | all 151 corrected targets\nM / Non-M = enriched myeloid / non-myeloid RR8 populations; main QC, age-balanced contrast\nGrey = not estimable; black dot = q < 0.05 in that method / assay family',fontsize=15,color='black',pad=20);ax.set_facecolor('#DCE1E7');yi,xi=np.where(np.array(qcols).T<.05);ax.scatter(xi,yi,c='black',s=3);ax.axvline(6.5,c=gray,lw=1.8);fig.colorbar(im,ax=ax,shrink=.2,label='Flight log2 effect (assay-specific scale; clipped at +/-2)');fig.subplots_adjust(left=.15,right=.94,bottom=.11,top=.955)
save(fig,'figureS1_independent_complete_panel',['independent_bulk_panel.csv','OSD536_panel_results.csv','sorted_cell_panel_complete.csv','../core/tables/verified_panel_mapping.csv'])
(F/'figure_source_ledger.json').write_text(json.dumps(ledger,indent=2))
print('Rendered',len(ledger),'figures in PNG, PDF and editable SVG',flush=True)
