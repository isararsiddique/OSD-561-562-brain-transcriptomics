from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Source-linked figures; missing estimates remain missing, never zero-filled."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
B=stage_root('deep');F=B/'figures';R=stage_root('extended')
plt.rcParams.update({'font.family':publication_font(),'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.unicode_minus':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
blue='#245C83';teal='#278A82';red='#B95858';gray='#454545';statecolors=['#3B6B86','#579E94','#D5A464','#9B789B']
ledger=[]
def save(fig,name,sources,pdf=True):
    for ext in (['png','pdf','svg'] if pdf else ['png','svg']):fig.savefig(F/(name+'.'+ext),dpi=210,bbox_inches='tight')
    plt.close(fig);ledger.append({'figure':name,'sources':sources,'script':__file__})
def title(fig,text,subtitle):
    fig.suptitle(text,x=.045,y=.985,ha='left',fontweight='bold',fontsize=19,color='black')
    fig.text(.045,.935,subtitle,color=gray,fontsize=10.5)
def panel(ax,text):ax.set_title(text,loc='left',fontweight='bold',fontsize=12,color='black',pad=10)
def pformat(q):return f'{q:.3g}'
qc=pd.read_csv(B/'tables/state_mouse_QC.csv');dq=pd.read_csv(B/'depth_sensitivity/tables/state_mouse_QC.csv')
reg=qc.drop_duplicates('donor').sort_values(['flight','age','donor']);donors=reg.donor.tolist()
state=pd.read_csv(B/'tables/state_panel_complete.csv');depth=pd.read_csv(B/'depth_sensitivity/tables/state_panel_complete.csv')
fig,axs=plt.subplots(2,2,figsize=(16,11));title(fig,'Forebrain expression states and limits of within-state inference','RR8: 23,036 retained cells from eight mice. State labels are numerical; exposure and age were excluded from feature selection and clustering.')
ax=axs[0,0];s=qc[qc.K.eq(4)].pivot(index='donor',columns='state',values='cell_fraction').reindex(donors)
bottom=np.zeros(8)
for j in range(4):
    ax.bar(np.arange(8),s[j+1]*100,bottom=bottom,color=statecolors[j],width=.72,label=f'State {j+1}');bottom+=s[j+1].to_numpy()*100
ax.set_xticks(range(8),[d.replace('GC_Mouse','G').replace('FLT_Mouse','F')+'\n'+str(int(a))+' mo' for d,a in zip(reg.donor,reg.age)],fontsize=9)
ax.set(ylim=(0,105),ylabel='Fraction of retained cells (%)');ax.axvline(3.5,color=gray,lw=.8)
ax.legend(ncol=4,frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18));panel(ax,'A  Captured state proportions | primary K=4')
ax=axs[0,1];mat=[];labels=[]
for name,q in [('Primary',qc),('Depth sensitivity',dq)]:
    m=q[q.K.eq(4)].pivot(index='state',columns='donor',values='n_cells').reindex(columns=donors)
    mat.extend(m.to_numpy());labels.extend([name+' · S'+str(i) for i in m.index])
mat=np.array(mat);im=ax.imshow(np.log10(mat+1),cmap='Blues',aspect='auto',vmin=0,vmax=3.6)
for i in range(8):
    for j in range(8):ax.text(j,i,str(int(mat[i,j])),ha='center',va='center',fontsize=8.5,color='white' if mat[i,j]>650 else 'black')
    if (mat[i]>=50).all():ax.add_patch(Rectangle((-.5,i-.5),8,1,fill=False,edgecolor=teal,lw=2))
ax.set_xticks(range(8),[d.replace('GC_Mouse','G').replace('FLT_Mouse','F') for d in donors],fontsize=9)
ax.set_yticks(range(8),labels,fontsize=9);ax.axhline(3.5,c='white',lw=3)
fig.colorbar(im,ax=ax,fraction=.045,pad=.025,label='log10(cell count + 1)');panel(ax,'B  Cells from each mouse | K=4')
eligible=[]
for label,data in [('Primary',state),('Depth sensitivity',depth)]:
    z=data[data.contrast.eq('flight_balanced_age')&data.analyte.eq('PYCARD')&data.status.eq('tested')]
    for k,st in z[['K','state']].drop_duplicates().itertuples(index=False,name=None):eligible.append((label,k,st,data))
ax=axs[1,0];display=[]
for i,(label,k,st,data) in enumerate(eligible):
    for mi,method in enumerate(['limma_voom','edgeR_QL']):
        r=data[data.K.eq(k)&data.state.eq(st)&data.method.eq(method)&data.analyte.eq('PYCARD')&data.contrast.eq('flight_balanced_age')].iloc[0]
        y=i+(mi-.5)*.20;c=blue if mi==0 else teal
        if mi==0:ax.plot([r.ci_low,r.ci_high],[y,y],c=c,lw=1.5)
        ax.scatter(r.log2FC,y,s=42,marker='o' if mi==0 else 's',facecolor=c if r.q_full_state_family<.05 else 'white',edgecolor=c,zorder=3)
        ax.text(.04,y,'q='+pformat(r.q_full_state_family),fontsize=8.5,va='center')
        display.append({'analysis':label,'K':k,'state':st,'method':method,'analyte':'PYCARD','log2FC':r.log2FC,'q':r.q_full_state_family})
ax.set_yticks(range(len(eligible)),[f'{l} · K={k} · S{st}' for l,k,st,_ in eligible],fontsize=9)
ax.set(ylim=(len(eligible)-.5,-.6),xlim=(-1.65,.46),xlabel='PYCARD flight log2 fold change');ax.axvline(0,c=gray,ls='--',lw=.8)
panel(ax,'C  PYCARD in all estimable mouse-state models')
ax=axs[1,1];comp=pd.read_csv(B/'tables/state_composition_complete_tests.csv')
focus=['AIF1','PYCARD','AIM2','NINJ1']
for i,g in enumerate(focus):
    for mi,measure in enumerate(['observed_composition','uniform_all_states']):
        r=comp[comp.K.eq(4)&comp.analyte.eq(g)&comp.measure.eq(measure)&comp.contrast.eq('flight_balanced_age')].iloc[0]
        assert r.status=='tested';y=i+(mi-.5)*.24;c=blue if mi==0 else red
        ax.plot([r.ci_low,r.ci_high],[y,y],c=c,lw=1.4)
        ax.scatter(r.effect,y,s=45,facecolor=c if r.HC3_q_full_1812_family<.05 else 'white',edgecolor=c,zorder=3,
                   label=['Observed mixture','Equal state weights'][mi] if i==0 else None)
ax.set_yticks(range(4),focus);ax.set(ylim=(3.5,-.6),xlabel='Flight effect: log2(mean cell CP10K + 1)');ax.axvline(0,c=gray,ls='--',lw=.8)
ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.20),ncol=2,fontsize=9);panel(ax,'D  Fixed-mixture sensitivity | eight mice')
fig.text(.045,.105,'B: outlined rows have at least 50 cells from every mouse. Depth sensitivity regresses cell UMI and mitochondrial fraction from PCs; state numbers do not match across fits.',fontsize=9.5,color=gray)
fig.text(.045,.072,'C: circles = limma/voom with 95% moderated intervals; squares = edgeR QL. Filled = full-family q < 0.05. Finer depth-sensitive states weaken PYCARD support.',fontsize=9.5,color=gray)
fig.text(.045,.039,'D: HC3 mouse-level intervals; filled = q < 0.05 among 1,812 planned tests. Equal weights use sparse observed states (minimum 9 cells). Flight and processing day remain confounded.',fontsize=9.2,color=gray)
fig.subplots_adjust(top=.845,bottom=.22,left=.13,right=.96,wspace=.50,hspace=.48)
pd.DataFrame(display).to_csv(B/'tables/figure8_displayed_count_effects.csv',index=False)
save(fig,'figure8_expression_state_validation',['tables/state_mouse_QC.csv','tables/state_panel_complete.csv','tables/state_composition_complete_tests.csv','depth_sensitivity/tables/state_mouse_QC.csv','depth_sensitivity/tables/state_panel_complete.csv'])

new=pd.read_csv(B/'tables/additional_panel_all_regions.csv');new['q_full_1057_family']=new.q_expanded_1812_family;reg=pd.read_csv(B/'tables/additional_mouse_registry.csv')
ann=pd.read_csv(B/'prepared/additional_gene_annotation.csv');deletion=pd.read_csv(B/'tables/additional_mouse_deletion_expanded_focus.csv');deletion['q_full_1057_family']=deletion.q_expanded_1812_family
fig,axs=plt.subplots(2,3,figsize=(17,11));title(fig,'Additional cohorts reveal control-dependent and glial RNA signals','Primary count reanalysis: RR23 cerebellum (26 mice) and MHU-3 cerebrum (24 mice). Expanded 1,812-test family includes two additional cerebral-region analyses.')
points=[]
def donor_plot(ax,sid,symbol):
    E=pd.read_csv(B/'prepared'/f'{sid}_voom_logCPM.csv.gz',index_col=0)
    ids=ann.loc[ann.SYMBOL.eq(symbol),'ENSEMBL'];ids=ids[ids.isin(E.index)];assert len(ids)==1
    s=reg[reg.OSD.eq(sid)].copy()
    groups=[('GC','WT'),('FLT','WT'),('VIV','WT')] if sid=='OSD-525' else [('GC','WT'),('FLT','WT'),('GC','Nrf2KO'),('FLT','Nrf2KO')]
    labels=['Habitat\nn=9','Flight\nn=9','Vivarium\nn=8'] if sid=='OSD-525' else ['Ground WT\nn=6','Flight WT\nn=6','Ground KO\nn=6','Flight KO\nn=6']
    for i,(f,g) in enumerate(groups):
        ss=s[s.flight.eq(f)&s.genotype.eq(g)].sort_values('donor');y=E.loc[ids.iloc[0],ss.key].to_numpy();c=red if f=='FLT' else teal if f=='GC' else gray
        ax.scatter(i+np.linspace(-.14,.14,len(y)),y,s=35,color=c,zorder=3);ax.plot([i-.19,i+.19],[y.mean()]*2,c='black',lw=1.3)
        for d,yy in zip(ss.donor,y):points.append({'OSD':sid,'gene':symbol,'donor':d,'flight':f,'genotype':g,'voom_logCPM':yy})
    ax.set_xticks(range(len(groups)),labels,fontsize=9);ax.set_ylabel('TMM / voom log2 CPM');ax.set_xlim(-.5,len(groups)-.5)
donor_plot(axs[0,0],'OSD-525','Aim2');panel(axs[0,0],'A  RR23 | AIM2 | individual mice')
donor_plot(axs[0,1],'OSD-525','Pycard');panel(axs[0,1],'B  RR23 | PYCARD | individual mice')
ax=axs[0,2];labels=[];display=[]
contrasts=['flight_vs_habitat_ground','flight_vs_vivarium','habitat_ground_vs_vivarium'];short=['F vs habitat','F vs vivarium','Habitat vs vivarium']
for gi,gene in enumerate(['AIM2','PYCARD']):
    for ci,(c,cl) in enumerate(zip(contrasts,short)):
        i=gi*3+ci;labels.append(gene+' | '+cl)
        for mi,method in enumerate(['limma_voom','edgeR_QL']):
            r=new[new.OSD.eq('OSD-525')&new.analyte.eq(gene)&new.contrast.eq(c)&new.method.eq(method)].iloc[0];y=i+(mi-.5)*.18;color=blue if mi==0 else teal
            if mi==0:ax.plot([r.ci_low,r.ci_high],[y,y],c=color,lw=1.4)
            ax.scatter(r.log2FC,y,s=35,marker='o' if mi==0 else 's',facecolor=color if r.q_full_1057_family<.05 else 'white',edgecolor=color,zorder=3)
            display.append({'OSD':'OSD-525','analyte':gene,'contrast':c,'method':method,'log2FC':r.log2FC,'q':r.q_full_1057_family})
ax.set_yticks(range(6),labels,fontsize=8.5);ax.set(ylim=(5.5,-.5),xlabel='Flight / control log2 fold change');ax.axvline(0,c=gray,ls='--',lw=.8);panel(ax,'C  Control choice changes the inference')
donor_plot(axs[1,0],'OSD-457','Slc1a3');panel(axs[1,0],'D  MHU-3 | SLC1A3 | individual mice')
ax=axs[1,1];genes=['SLC1A3','S100B','AIF1','PYCARD','AIM2']
for i,gene in enumerate(genes):
    r=new[new.OSD.eq('OSD-457')&new.analyte.eq(gene)&new.contrast.eq('flight_by_Nrf2KO')&new.method.eq('limma_voom')].iloc[0]
    ax.plot([r.ci_low,r.ci_high],[i,i],c=blue,lw=1.4);ax.scatter(r.log2FC,i,s=45,facecolor=blue if r.q_full_1057_family<.05 else 'white',edgecolor=blue)
ax.set_yticks(range(5),genes);ax.set(ylim=(4.5,-.5),xlabel='Interaction: flight effect in KO minus WT');ax.axvline(0,c=gray,ls='--',lw=.8);panel(ax,'E  Nrf2 interaction remains uncertain')
ax=axs[1,2];checks=[('OSD-525','AIM2','flight_vs_habitat_ground','RR23 AIM2 | F vs habitat'),('OSD-525','PYCARD','flight_vs_vivarium','RR23 PYCARD | F vs vivarium'),('OSD-457','SLC1A3','flight_WT','MHU-3 SLC1A3 | WT'),('OSD-457','S100B','flight_balanced_genotype','MHU-3 S100B | balanced')]
for i,(sid,g,c,label) in enumerate(checks):
    for mi,method in enumerate(['limma_voom','edgeR_QL']):
        s=deletion[deletion.OSD.eq(sid)&deletion.analyte.eq(g)&deletion.contrast.eq(c)&deletion.method.eq(method)].sort_values('deleted_mouse');color=blue if mi==0 else teal
        y=i+(mi-.5)*.20;ax.plot([s.log2FC.min(),s.log2FC.max()],[y,y],color=color,lw=1.4)
        for j,r in enumerate(s.itertuples()):ax.scatter(r.log2FC,y+.03*np.sin(j*2.3),s=13,marker='o' if mi==0 else 's',facecolor=color if r.q_full_1057_family<.05 else 'white',edgecolor=color,lw=.6,zorder=3)
ax.set_yticks(range(4),[x[3] for x in checks],fontsize=8.5);ax.set(ylim=(3.5,-.5),xlim=(-.9,.9),xlabel='Log2 effect after removing each mouse');ax.axvline(0,c=gray,ls='--',lw=.8);panel(ax,'F  Every mouse deletion | both methods')
fig.text(.045,.110,'A/B/D: each point is one mouse; bars show means. C/E: limma/voom 95% moderated intervals. C/F: circles = limma/voom; squares = edgeR QL; filled = full-family q < 0.05.',fontsize=9.4,color=gray)
fig.text(.045,.075,'F: 26 or 24 deletion estimates per method. AIM2 survives every deletion in the expanded family. SLC1A3, S100B and PYCARD retain direction; corrected support varies with deletion.',fontsize=9.4,color=gray)
fig.text(.045,.040,'RR23 collection date differs by flight/control; narrative and ISA disagree for vivarium. Current MHU-3 has six KO-flight samples, versus five in the paper caption; all KO deletions were tested.',fontsize=9.2,color=gray)
fig.subplots_adjust(top=.845,bottom=.22,left=.075,right=.97,wspace=.80,hspace=.48)
pd.DataFrame(points).to_csv(B/'tables/figure9_displayed_mouse_values.csv',index=False)
pd.DataFrame(display).to_csv(B/'tables/figure9_displayed_contrasts.csv',index=False)
save(fig,'figure9_additional_cohort_validation',['tables/additional_panel_all_regions.csv','tables/additional_mouse_registry.csv','prepared/OSD-525_voom_logCPM.csv.gz','prepared/OSD-457_voom_logCPM.csv.gz','tables/additional_mouse_deletion_expanded_focus.csv'])

# Complete additional panel: quantitative effects, unavailable estimates in grey.
panelmap=pd.read_csv(stage_root('core') / 'tables/verified_panel_mapping.csv');vals=[];qs=[];labels=[]
for sid in ['OSD-525','OSD-457','OSD-612','OSD-613']:
    cs=new[new.OSD.eq(sid)].contrast.drop_duplicates()
    for c in cs:
        for method in ['limma_voom','edgeR_QL']:
            s=new[new.OSD.eq(sid)&new.contrast.eq(c)&new.method.eq(method)].set_index('analyte').reindex(panelmap.analyte)
            vals.append(s.log2FC.to_numpy());qs.append(s.q_full_1057_family.to_numpy());labels.append(sid+'\n'+c.replace('_',' ')+'\n'+('voom' if method=='limma_voom' else 'edgeR'))
fig,ax=plt.subplots(figsize=(17,24));mat=np.array(vals).T;im=ax.imshow(np.ma.masked_invalid(mat),cmap='RdBu_r',vmin=-1.5,vmax=1.5,aspect='auto');ax.set_facecolor('#D9DFE3')
ax.set_xticks(range(len(labels)),labels,rotation=90,fontsize=8);ax.set_yticks(range(151),panelmap.analyte,fontsize=7)
yi,xi=np.where(np.array(qs).T<.05);ax.scatter(xi,yi,c='black',s=5)
ax.set_title('Additional cohorts | all 151 corrected targets\nGrey = unestimable; black point = q < 0.05 in the expanded 1,812-test family per method\nAssay-specific effects shown; no pooled cross-region effect',fontweight='bold',color='black',fontsize=16,pad=20)
fig.colorbar(im,ax=ax,shrink=.23,label='Log2 fold change (display clipped at ±1.5)');fig.subplots_adjust(left=.13,right=.96,bottom=.15,top=.955)
save(fig,'figureS3_additional_complete_panel',['tables/additional_panel_all_regions.csv','../core/tables/verified_panel_mapping.csv'],pdf=False)
new[['OSD','contrast','method','analyte','log2FC','q_expanded_1812_family','status']].to_csv(B/'tables/figureS3_displayed_complete_effects.csv',index=False)

# Different animals from the same mission and the stricter article QC gate.
terminal=pd.read_csv(B/'tables/terminal_region_mouse_registry.csv')
article=pd.read_csv(B/'tables/article_5pct_panel_complete.csv')
fig,axs=plt.subplots(2,2,figsize=(16,10));title(fig,'Regional transfer and the article-defined cell QC check','RRRM-2: 18 different cerebral-region mice. RR10: nine bulk cerebral samples from already analyzed animals. RR8: article 5% versus GEO 50% mitochondrial gate.')
display10=[];mouse10=[]
ax=axs[0,0];fields=[('OSD-613','flight_balanced_age','RRRM-2 cerebral | 18 different mice'),('OSD-612','flight_vs_ground','RR10 cerebral | 9 shared mice'),('OSD-525','flight_vs_habitat_ground','RR23 cerebellum | habitat contrast'),('OSD-457','flight_WT','MHU-3 cerebrum | WT contrast')]
for i,(sid,c,label) in enumerate(fields):
    for mi,method in enumerate(['limma_voom','edgeR_QL']):
        r=new[new.OSD.eq(sid)&new.analyte.eq('PYCARD')&new.contrast.eq(c)&new.method.eq(method)].iloc[0]
        y=i+(mi-.5)*.22;col=blue if mi==0 else teal
        if mi==0:ax.plot([r.ci_low,r.ci_high],[y,y],c=col,lw=1.5)
        ax.scatter(r.log2FC,y,marker='o' if mi==0 else 's',s=40,facecolor=col if r.q_expanded_1812_family<.05 else 'white',edgecolor=col)
        display10.append({'panel':'A','source':'additional_panel_all_regions','OSD':sid,'tissue':'','QC':'','analyte':'PYCARD','contrast':c,'method':method,'log2FC':r.log2FC,'q':r.q_expanded_1812_family})
ax.set_yticks(range(4),[x[2] for x in fields],fontsize=9);ax.set(ylim=(3.5,-.5),xlabel='PYCARD flight log2 fold change');ax.axvline(0,c=gray,ls='--',lw=.8)
panel(ax,'A  PYCARD transfer is not uniformly confirmed')
ax=axs[0,1];E=pd.read_csv(B/'prepared/OSD-613_voom_logCPM.csv.gz',index_col=0);ids=ann.loc[ann.SYMBOL.eq('Pycard'),'ENSEMBL'];ids=ids[ids.isin(E.index)];assert len(ids)==1
s=terminal[terminal.OSD.eq('OSD-613')]
for i,(age,f) in enumerate([(12,'GC'),(12,'FLT'),(29,'GC'),(29,'FLT')]):
    ss=s[s.age.eq(age)&s.flight.eq(f)].sort_values('donor');v=E.loc[ids.iloc[0],ss.key].to_numpy();ax.scatter(i+np.linspace(-.11,.11,len(v)),v,c=teal if f=='GC' else red,s=45);ax.plot([i-.16,i+.16],[v.mean()]*2,c=gray,lw=1.2)
    for donor,val in zip(ss.donor,v):mouse10.append({'OSD':'OSD-613','gene':'Pycard','donor':donor,'voom_logCPM':val})
ax.set_xticks(range(4),['12 wk\nGround n=5','12 wk\nFlight n=4','29 wk\nGround n=4','29 wk\nFlight n=5'],fontsize=9)
ax.set_ylabel('PYCARD: TMM / voom log2 CPM');panel(ax,'B  Different RRRM-2 mice | cerebral region')
ax=axs[1,0];genes=['AIF1','PYCARD','AIM2','NINJ1']
for i,g in enumerate(genes):
    for si,(qc_name,data,col) in enumerate([('GEO 50% gate',pd.read_csv(R/'tables/sorted_cell_panel_complete.csv'),blue),('Article 5% gate',article,teal)]):
        z=data[data.tissue.eq('Myeloid_Forebrain')&data.analyte.eq(g)&data.contrast.eq('flight_balanced_age')&data.method.eq('limma_voom')]
        if si==0:z=z[z.QC.eq('published_gate')]
        r=z.iloc[0];y=i+(si-.5)*.25
        ax.scatter(r.log2FC,y,s=50,facecolor=col if r.q_full_3624_family<.05 else 'white',edgecolor=col,label=qc_name if i==0 else None)
        display10.append({'panel':'C','source':'sorted_cell_panel_complete' if si==0 else 'article_5pct_panel_complete','OSD':'RR8','tissue':'Myeloid_Forebrain','QC':'published_gate' if si==0 else 'article_5pct_gate','analyte':g,'contrast':'flight_balanced_age','method':'limma_voom','log2FC':r.log2FC,'q':r.q_full_3624_family})
ax.set_yticks(range(4),genes);ax.set(ylim=(3.5,-.5),xlabel='Forebrain flight log2 fold change');ax.axvline(0,c=gray,ls='--',lw=.8)
ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.17),ncol=2,fontsize=9);panel(ax,'C  Lead program survives stricter article QC')
ax=axs[1,1];names=['Cerebellum','Forebrain','Hippocampus'];labels=[]
for ti,t in enumerate(names):
    for gi,g in enumerate(['AIF1','PYCARD']):
        i=ti*2+gi;labels.append(t+' | '+g)
        for mi,method in enumerate(['limma_voom','edgeR_QL']):
            r=article[article.tissue.eq('Myeloid_'+t)&article.analyte.eq(g)&article.contrast.eq('flight_balanced_age')&article.method.eq(method)].iloc[0]
            y=i+(mi-.5)*.18;col=blue if mi==0 else teal
            ax.scatter(r.log2FC,y,s=38,marker='o' if mi==0 else 's',facecolor=col if r.q_full_3624_family<.05 else 'white',edgecolor=col)
            display10.append({'panel':'D','source':'article_5pct_panel_complete','OSD':'RR8','tissue':'Myeloid_'+t,'QC':'article_5pct_gate','analyte':g,'contrast':'flight_balanced_age','method':method,'log2FC':r.log2FC,'q':r.q_full_3624_family})
ax.set_yticks(range(6),labels,fontsize=9);ax.set(ylim=(5.5,-.5),xlabel='Article-gated flight log2 fold change');ax.axvline(0,c=gray,ls='--',lw=.8)
panel(ax,'D  Stricter QC preserves regional uncertainty')
fig.text(.045,.105,'A: circles = limma/voom with 95% moderated intervals; squares = edgeR QL. Filled = expanded 1,812-test q < 0.05. B: each point is one mouse; no new animal count from RR10.',fontsize=9.3,color=gray)
fig.text(.045,.070,'C/D: 200–2,500 features and ≤5% mitochondrial RNA retain 62,351 cells across 47 libraries; each method keeps the full 3,624-target contrast family. Points are estimates, not cells.',fontsize=9.3,color=gray)
fig.text(.045,.035,'The source threshold discrepancy is retained explicitly. Whole-mouse inference, processing-day and preservation confounds, and unconfirmed hippocampal PYCARD remain.',fontsize=9.3,color=gray)
fig.subplots_adjust(top=.83,bottom=.23,left=.20,right=.96,wspace=.58,hspace=.47)
pd.DataFrame(display10).to_csv(B/'tables/figure10_displayed_effects.csv',index=False)
pd.DataFrame(mouse10).to_csv(B/'tables/figure10_displayed_mouse_values.csv',index=False)
save(fig,'figure10_regional_and_QC_validation',['tables/additional_panel_all_regions.csv','tables/terminal_region_mouse_registry.csv','prepared/OSD-613_voom_logCPM.csv.gz','tables/article_5pct_panel_complete.csv','tables/article_GEO_QC_discrepancy.json','../extended/tables/sorted_cell_panel_complete.csv'])
(B/'tables/deep_figure_source_ledger.json').write_text(json.dumps(ledger,indent=2))
