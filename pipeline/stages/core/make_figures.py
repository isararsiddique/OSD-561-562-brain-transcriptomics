from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Scientific figures from executed tables; no new tests are created by plots."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
import seaborn as sns
OUT=stage_root('core');TAB=OUT/'tables';PRE=OUT/'prepared';FIG=OUT/'figures'
sns.set_theme(style='ticks',context='paper',font='DejaVu Sans')
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False})
blue='#2166AC';red='#B54343';gray='#717780';purple='#72559A'
bulk=pd.read_csv(TAB/'bulk_panel_complete.csv');spatial=pd.read_csv(TAB/'spatial_panel_complete.csv');mapping=pd.read_csv(TAB/'verified_panel_mapping.csv').fillna('')
meta=pd.read_csv(TAB/'bulk_sample_registry.csv');registry=[]
def save(fig,name,caption,sources):
    fig.savefig(FIG/(name+'.png'),dpi=175,bbox_inches='tight',facecolor='white')
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)
    registry.append({'figure':name,'caption':caption,'source_tables':';'.join(sources),'independent_evidence':'Replotting a table adds no independent replication.'})
    print(name,flush=True)
def point_groups(ax,values,labels,colors):
    for i,(label,v,c) in enumerate(zip(labels,values,colors)):
        off=np.linspace(-.08,.08,len(v));ax.scatter(i+off,v,color=c,s=28,zorder=3)
        ax.plot([i-.18,i+.18],[np.mean(v)]*2,color=c,lw=2)
    ax.set_xticks(range(len(labels)),labels)
def gene_points(osd,gene):
    b=bulk[(bulk.study==f'OSD-{osd}') & (bulk.mouse_symbol==gene)].iloc[0]
    n=pd.read_csv(PRE/f'OSD-{osd}_normalized_counts.csv.gz',index_col=0)
    m=meta[meta.study==f'OSD-{osd}'].copy();m['expression']=np.log2(n.loc[b.gene_id,m['sample']].to_numpy()+1)
    return m

# Figure 1: identity, count inclusion, expression QC, and actual group n.
fig,axes=plt.subplots(2,2,figsize=(12.5,9),layout='constrained')
groups=['GC_ISS_T_12','FLT_ISS_T_12','GC_ISS_T_29','FLT_ISS_T_29','GC_LAR_12','FLT_LAR_12','GC_LAR_29','FLT_LAR_29']
gn=meta.groupby(['study','group']).size().unstack().reindex(index=['OSD-561','OSD-562'],columns=groups)
im=axes[0,0].imshow(gn,aspect='auto',cmap='Blues',vmin=0,vmax=5)
for i in range(2):
    for j in range(8):axes[0,0].text(j,i,str(gn.iloc[i,j]),ha='center',va='center',fontsize=12,color='white' if gn.iloc[i,j]>3 else '#172D45')
axes[0,0].set_xticks(range(8),['G T12','F T12','G T29','F T29','G R12','F R12','G R29','F R29'],rotation=35,ha='right');axes[0,0].set_yticks([0,1],['Cerebellum','Hippocampus'])
axes[0,0].set_title('A  Actual sequenced animals per design cell');axes[0,0].set_xlabel('G = ground; F = flight; T = terminal; R = returned')
cover=[128,127]+[int(spatial[(spatial.study==f'OSD-{x}') & (spatial.contrast=='flight_in_saline')].pvalue.notna().sum()) for x in [682,685,698,699]]
names=['561\nCb','562\nHPC','682\nCA1','685\nDG','698\nFCT','699\nCT']
axes[0,1].bar(names,[151]*6,color='#D9DDE3',label='Corrected analytes: 151');axes[0,1].bar(names,cover,color=blue,label='Primary mapped and testable')
for i,n in enumerate(cover):axes[0,1].text(i,n+2,str(n),ha='center',fontsize=10)
axes[0,1].set_ylim(0,175);axes[0,1].set_ylabel('Analytes');axes[0,1].set_title('B  Mapping and assay coverage');axes[0,1].legend(fontsize=8,loc='upper right')
for osd,ax in zip([561,562],axes[1]):
    norm=pd.read_csv(PRE/f'OSD-{osd}_normalized_counts.csv.gz',index_col=0);log=np.log2(norm+1)
    chosen=log.var(axis=1).nlargest(2000).index;x=log.loc[chosen].T.to_numpy(copy=True);x-=x.mean(axis=0)
    u,s,v=np.linalg.svd(x,full_matrices=False);pcs=u[:,:2]*s[:2];expl=s[:2]**2/np.sum(s**2)
    mm=meta[meta.study==f'OSD-{osd}'].set_index('sample').loc[log.columns]
    pc=pd.DataFrame({'sample':mm.index,'PC1':pcs[:,0],'PC2':pcs[:,1],'flight':mm.flight.to_numpy(),'cohort':mm.cohort.to_numpy(),'age':mm.age.to_numpy()});pc.to_csv(TAB/f'OSD-{osd}_QC_PCA.csv',index=False)
    for f in ['GC','FLT']:
        for c in ['ISS_T','LAR']:
            for age in [12,29]:
                mask=(mm.flight==f)&(mm.cohort==c)&(mm.age==age)
                ax.scatter(pcs[mask,0],pcs[mask,1],c=blue if f=='GC' else red,marker='o' if c=='ISS_T' else '^',s=45 if age==12 else 85,alpha=.85,edgecolor='white',lw=.7)
    ax.set_xlabel(f'PC1 ({expl[0]:.1%})');ax.set_ylabel(f'PC2 ({expl[1]:.1%})');ax.set_title(('C  Cerebellum, n = 26' if osd==561 else 'D  Hippocampus, n = 27'))
axes[1,1].legend(handles=[Line2D([],[],marker='o',ls='',color=blue,label='Ground'),Line2D([],[],marker='o',ls='',color=red,label='Flight'),Line2D([],[],marker='o',ls='',color=gray,label='Terminal'),Line2D([],[],marker='^',ls='',color=gray,label='Returned')],fontsize=8)
fig.suptitle('Corrected CNS panel: verified sample structure and QC\n53 libraries • 29 source animals • 24 animals shared across regions',fontsize=15)
save(fig,'fig01_design_and_QC','Counts are sequenced animals within region. PCA uses centered log2 normalized endogenous expression of the 2,000 highest-variance genes; age 29 uses larger markers. Terminal flight preservation is confounded with exposure.',['bulk_sample_registry.csv','verified_panel_mapping.csv','bulk_panel_complete.csv','spatial_panel_complete.csv'])

# Full retained panel maps: unavailable cells are visible, never zero-filled.
order=mapping.sort_values(['category','human_symbol']).analyte.tolist();human=mapping.set_index('analyte').human_symbol.to_dict()
ctags=['flight_ISS_T_age12','flight_ISS_T_age29','flight_LAR_age12','flight_LAR_age29','age_GC_ISS_T','age_GC_LAR','age_FLT_ISS_T','age_FLT_LAR','cohort_GC_age12','cohort_GC_age29','cohort_FLT_age12','cohort_FLT_age29','flight_by_cohort_age12','flight_by_cohort_age29']
clabs=['F T12','F T29','F R12','F R29','Age G T','Age G R','Age F T','Age F R','T−R G12','T−R G29','T−R F12','T−R F29','F×C 12','F×C 29']
def full_map(frame,studies,tags,labels,name,title,qcol):
    cols=[(s,t) for s in studies for t in tags]
    m=frame.pivot(index='analyte',columns=['study','contrast'],values='log2FoldChange').reindex(index=order,columns=pd.MultiIndex.from_tuples(cols))
    q=frame.pivot(index='analyte',columns=['study','contrast'],values=qcol).reindex(index=order,columns=pd.MultiIndex.from_tuples(cols))
    fig,ax=plt.subplots(figsize=(16.5,23));fig.subplots_adjust(left=.14,right=.92,bottom=.09,top=.94)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#D4D8DD');im=ax.imshow(np.ma.masked_invalid(m),aspect='auto',cmap=cmap,vmin=-3,vmax=3,interpolation='nearest')
    ax.set_yticks(range(len(order)),[human[a] for a in order],fontsize=6.3)
    ax.set_xticks(range(len(cols)),[f'{s}\n{lab}' for s in studies for lab in labels],rotation=90,fontsize=8)
    ii,jj=np.where(q.to_numpy()<.05);ax.scatter(jj,ii,s=9,c='black',marker='o')
    cats=mapping.set_index('analyte').loc[order,'category'];last=0
    for i in range(1,len(order)+1):
        if i==len(order) or cats.iloc[i]!=cats.iloc[last]:
            ax.axhline(i-.5,color='white',lw=1.5);ax.text(len(cols)-.15,(last+i-1)/2,cats.iloc[last],rotation=90,ha='left',va='center',fontsize=7,clip_on=False);last=i
    for i in range(1,len(studies)):ax.axvline(i*len(tags)-.5,color='black',lw=1.3)
    cb=fig.colorbar(im,ax=ax,fraction=.022,pad=.04,shrink=.28);cb.set_label('Unshrunk log2 effect (color clipped ±3)')
    ax.set_title(title+'\nBlack dots: full primary-family q < 0.05; gray: mapping/filter/assay unavailable',fontsize=14,pad=16)
    save(fig,name,'All 151 corrected source analytes are retained. Rows use current human symbols; original labels and orthology ambiguity remain in the tables. Color saturation is display clipping, not changed estimates.',['verified_panel_mapping.csv','bulk_panel_complete.csv' if studies==['OSD-561','OSD-562'] else 'spatial_panel_complete.csv'])
full_map(bulk,['OSD-561','OSD-562'],ctags,clabs,'fig02_complete_bulk_panel','The complete corrected panel across 28 bulk contrasts','q_panel_primary_full_family')

# Figure 3: actual lead-gene animal values plus every flight stratum CI.
fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
for j,gene in enumerate(['Aif1','Pycard']):
    d=bulk[(bulk.mouse_symbol==gene)&(bulk.contrast_kind=='flight')].set_index(['study','contrast']).reindex(pd.MultiIndex.from_product([['OSD-561','OSD-562'],ctags[:4]],names=['study','contrast'])).reset_index()
    ax=axes[0,j];yy=np.arange(len(d))
    for i,r in d.iterrows():ax.errorbar(r.log2FoldChange,i,xerr=1.96*r.lfcSE,fmt='o',c=red if r.q_panel_primary_full_family<.05 else gray,capsize=2)
    ax.axvline(0,c='#ADB3BC',lw=1);ax.set_yticks(yy,[('Cb' if s=='OSD-561' else 'HPC')+' '+t.replace('flight_','').replace('_age',' ') for s,t in zip(d.study,d.contrast)],fontsize=8);ax.invert_yaxis();ax.set_xlabel('Flight − ground log2 fold change; Wald 95% CI');ax.set_title(f'{chr(65+j)}  {gene.upper()}: all eight flight comparisons')
    pts=gene_points(561,gene);pts=pts[(pts.age==12)&(pts.cohort=='ISS_T')];vals=[pts[pts.flight==f].expression.to_numpy() for f in ['GC','FLT']]
    ax=axes[1,j];point_groups(ax,vals,['Ground\nn = 3','Flight\nn = 5'],[blue,red]);ax.set_ylabel('log2(normalized STAR count + 1)')
    hit=d[(d.study=='OSD-561')&(d.contrast==ctags[0])].iloc[0];ax.set_title(f'{chr(67+j)}  Young terminal cerebellum\nlog2FC {hit.log2FoldChange:.2f}; primary q {hit.q_panel_primary_full_family:.3g}')
fig.suptitle('Two primary flight associations: stable direction, conditional interpretation\nRed = primary family q < 0.05; terminal preservation remains confounded',fontsize=14)
save(fig,'fig03_AIF1_PYCARD','Individual animal points and full-model Wald intervals. These intervals are nominal, not simultaneous. Results do not establish age or region specificity. Rounded RSEM and leave-one-out are shared-sample sensitivities; exact-gene sensitivity has no full-family q<0.05.',['bulk_panel_complete.csv','bulk_count_source_sensitivity.csv','bulk_flight_leave_one_out.csv','bulk_gene_exact_sensitivity.csv'])

# Figure 4: collection context, and an animal-paired regional check.
fig,axes=plt.subplots(2,2,figsize=(13.5,9.5),layout='constrained')
for gene,ax,panel in zip(['Fos','Junb'],axes[0],['A','B']):
    pts=gene_points(562,gene);vals=[pts[pts.group==g].expression.to_numpy() for g in groups]
    point_groups(ax,vals,['G T12','F T12','G T29','F T29','G R12','F R12','G R29','F R29'],[blue,red]*4)
    ax.tick_params(axis='x',rotation=40);ax.set_ylabel('log2(normalized STAR count + 1)');ax.set_title(f'{panel}  Hippocampus {gene}: collection-cohort pattern')
paired=pd.read_csv(TAB/'paired_region_flight_associations.csv');lead=paired[paired.mouse_symbol.isin(['Aif1','Pycard','Fos','Junb','Gfap','Cd63'])].sort_values('mouse_symbol')
for i,r in enumerate(lead.itertuples()):axes[1,0].errorbar(r.cerebellum_minus_hippocampus_flight_effect,i,xerr=2.086*r.se,fmt='o',c=gray,capsize=2)
axes[1,0].set_yticks(range(len(lead)),lead.mouse_symbol);axes[1,0].axvline(0,color='#AAB1BA',lw=1);axes[1,0].set_xlabel('Flight effect on within-animal Cb − HPC log expression');axes[1,0].set_title('C  Paired-region secondary association (n = 24)')
null=pd.read_csv(TAB/'paired_region_omnibus_null.csv').iloc[:,0].to_numpy();info=json.loads((TAB/'paired_region_omnibus.json').read_text())
axes[1,1].hist(null,bins=45,color='#C8D7E6',edgecolor='white');axes[1,1].axvline(info['omnibus_mean_squared_regional_flight_effect'],color=red,lw=2,label='Observed');axes[1,1].set_xlabel('Mean squared adjusted regional flight effect');axes[1,1].set_ylabel('Whole-animal permutation count');axes[1,1].set_title(f'D  Paired-region omnibus p = {info["pvalue"]:.3f}');axes[1,1].legend()
fig.suptitle('Collection context is prominent; paired regions do not confirm a regional flight pattern\nTerminal vs returned combines duration, recovery and sample handling',fontsize=14)
save(fig,'fig04_collection_context_and_pairing','Fos and Junb are measured in the same hippocampal animals. Paired analysis is within RRRM-2, not independent replication; per-gene paired q-values and the 127-gene omnibus are nonsignificant.',['bulk_panel_complete.csv','paired_region_flight_associations.csv','paired_region_omnibus.json','paired_region_omnibus_null.csv'])

# Figure 5: factorial spatial effects and uncertainty disagreements.
fig,axes=plt.subplots(1,2,figsize=(13,10),layout='constrained')
chosen=['AIF1','PYCARD','FOS','JUNB','GFAP','GSDMD'];s=spatial[(spatial.analyte.isin(chosen)) & (spatial.contrast=='flight_by_BuOE_interaction')].sort_values(['analyte','study'])
for i,r in enumerate(s.itertuples()):
    if np.isfinite(r.log2FoldChange):axes[0].errorbar(r.log2FoldChange,i,xerr=[[r.log2FoldChange-r.ci95_low],[r.ci95_high-r.log2FoldChange]],fmt='o',c=purple,capsize=2)
    else:axes[0].text(.05,i,'Not on assay',color=gray,fontsize=7,va='center')
axes[0].set_yticks(range(len(s)),[f'{r.analyte} · {r.region}' for r in s.itertuples()],fontsize=8);axes[0].invert_yaxis();axes[0].axvline(0,c='#AAB1BA',lw=1);axes[0].set_xlabel('Full flight × BuOE difference-in-differences\nlog2 expression; OLS 95% CI');axes[0].set_title('A  Explicit interaction estimates')
hits=spatial[spatial.q_hc3_panel_primary_full_family<.05].copy();xx=np.arange(len(hits));width=.25
for k,col,c,label in [(0,'pvalue',gray,'Pooled OLS'),(1,'hc3_pvalue',purple,'HC3, fixed df = 8'),(2,'welch_pvalue',blue,'Welch, group-specific df')]:
    axes[1].barh(xx+(k-1)*width,-np.log10(hits[col].to_numpy()),height=width,color=c,label=label)
axes[1].set_yticks(xx,[f'{r.analyte} · {r.region}\n{r.contrast.replace("_"," ")}' for r in hits.itertuples()],fontsize=8);axes[1].set_xlabel('−log10 nominal p (shared data; separate uncertainty models)');axes[1].set_title('B  All five HC3-only sensitivity calls');axes[1].legend(loc='lower right',fontsize=8)
fig.suptitle('GeoMx: no primary or Welch family-wide calls\n3 sections/ROIs per group; one slide per group; animal provenance unresolved',fontsize=14)
save(fig,'fig05_spatial_factorial_uncertainty','The full interaction is (FLT_BuOE−GC_BuOE)−(FLT_saline−GC_saline). Five fixed-df HC3 calls do not survive the Welch family correction. No model separates slide from group.',['spatial_panel_complete.csv','spatial_sample_registry.csv','spatial_welch_sensitivity.csv'])
stags=['flight_in_saline','flight_in_BuOE','BuOE_in_ground','BuOE_in_flight','flight_balanced_main','BuOE_balanced_main','flight_by_BuOE_interaction']
full_map(spatial,['OSD-682','OSD-685','OSD-698','OSD-699'],stags,['F Sal','F BuOE','BuOE G','BuOE F','F main','BuOE main','F×BuOE'],'figS01_complete_spatial_panel','All four-region spatial contrasts: no primary q < 0.05','q_panel_primary_full_family')

# Figure 6: independent mission context, respecting the five-animal n.
context=pd.read_csv(TAB/'OSD-352_panel_celltype_context.csv');types=['CNS neurons','CNS glia','Immune cells'];typelabel=['Neurons','CNS glia','Immune']
fig,axes=plt.subplots(2,3,figsize=(13,8.5),layout='constrained')
for ax,gene in zip(axes.flat,['Aif1','Pycard','Fos','Junb','Gfap','Casp1']):
    s=context[(context.mouse_symbol==gene)&context.cell_type.isin(types)]
    for i,ct in enumerate(types):
        for flight,off,c in [('ground',-.12,blue),('flight',.12,red)]:
            v=s[(s.cell_type==ct)&(s.flight_status==flight)].sort_values('bio_origin');ax.scatter(i+off+np.linspace(-.03,.03,len(v)),v.log2_CPM_plus1,color=c,s=30)
    ax.set_xticks(range(3),typelabel);ax.set_title(gene);ax.set_ylabel('Donor pseudobulk log2(CPM + 1)');ax.set_ylim(bottom=-.25)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=blue,label='Ground: 2 animals'),Line2D([],[],marker='o',ls='',color=red,label='Flight: 3 animals')],loc='outside lower center',ncols=2)
fig.suptitle('OSD-352: separately assessed RR-3 cell-type context\nAuthor labels; ≥20 nuclei in each donor/type; no BuOE; no corrected confirmation',fontsize=14)
save(fig,'fig06_RR3_celltype_context','21,178 author-annotated nuclei aggregate to five animals. Neurons, CNS glia and immune groups pass the fixed cell-count gate; vascular and neural crest-like groups do not. AIF1/PYCARD effects do not replicate the negative cerebellum direction.',['OSD-352_panel_celltype_context.csv','OSD-352_pseudobulk_registry.csv','OSD-352_panel_celltype_exact_tests.csv'])

# The requested five representations use the exact same genes and estimates.
sense=pd.read_csv(TAB/'bulk_count_source_sensitivity.csv');same=sense[(sense.study=='OSD-561')&(sense.contrast=='flight_ISS_T_age12')].dropna(subset=['log2FoldChange_STAR','log2FoldChange_RSEM'])
x=same.log2FoldChange_STAR.to_numpy();y=same.log2FoldChange_RSEM.to_numpy();lim=[min(x.min(),y.min())-.2,max(x.max(),y.max())+.2]
same.to_csv(TAB/'five_view_source_data.csv',index=False)
fig=plt.figure(figsize=(13,8.5),layout='constrained');gs=fig.add_gridspec(2,3);aa=[fig.add_subplot(gs[0,i]) for i in range(3)]+[fig.add_subplot(gs[1,i]) for i in range(2)]
aa[0].scatter(x,y,s=14,c=blue,alpha=.6)
for gene in ['Aif1','Pycard']:
    p=same[same.mouse_symbol==gene].iloc[0];aa[0].annotate(gene,(p.log2FoldChange_STAR,p.log2FoldChange_RSEM),xytext=(15,20) if gene=='Aif1' else (22,-9),textcoords='offset points',fontsize=8,arrowprops={'arrowstyle':'-','color':gray,'lw':.6})
aa[1].hexbin(x,y,gridsize=15,mincnt=1,cmap='Blues')
grid=np.linspace(*lim,100);xx,yy=np.meshgrid(grid,grid);z=gaussian_kde(np.vstack([x,y]))(np.vstack([xx.ravel(),yy.ravel()])).reshape(xx.shape)
aa[2].contour(xx,yy,z,levels=6,cmap='Blues');aa[2].scatter(x,y,s=3,c=gray,alpha=.25)
aa[3].hist(x,bins=18,density=True,histtype='step',lw=2,color=blue,label='STAR');aa[3].hist(y,bins=18,density=True,histtype='step',lw=2,color=red,label='Rounded RSEM');aa[3].legend(fontsize=8)
aa[4].hist2d(x,y,bins=18,cmap='Blues',cmin=1)
for i,ax in enumerate(aa):
    ax.set_title(['A  Scatter','B  Hexbin counts','C  Smoothed density contours','D  Marginal distributions','E  2D binned density'][i]);ax.set_xlabel('STAR log2 fold change' if i!=3 else 'Log2 fold change')
    if i!=3:ax.plot(lim,lim,'--',c=gray,lw=1);ax.set_xlim(lim);ax.set_ylim(lim);ax.set_ylabel('Rounded RSEM log2 fold change')
    else:ax.set_ylabel('Gene-distribution density')
note=fig.add_subplot(gs[1,2]);note.axis('off');note.text(.02,.95,f'The same {len(same)} testable genes\nappear in all five views.\n\nYoung terminal cerebellum.\nShared animals and sequencing reads.\n\nAgreement is a count-source sensitivity,\nnot an independent replication.\n\nGenes, bins and density estimates\ndo not increase animal n.',va='top',fontsize=11,linespacing=1.5)
fig.suptitle('Five complementary views of one count-source comparison',fontsize=15)
save(fig,'figS02_five_views','Scatter, hexbin, smoothed density contours, marginal distributions and 2D binned density use exactly the same gene estimates. Density is descriptive; shared-sample agreement is not independent validation.',['five_view_source_data.csv','bulk_count_source_sensitivity.csv'])
pd.DataFrame(registry).to_csv(TAB/'figure_source_ledger.csv',index=False)
(FIG/'CAPTIONS.md').write_text('# Figure captions and source tables\n\n'+'\n\n'.join('## '+r['figure']+'\n\n'+r['caption']+'\n\nSources: '+r['source_tables'] for r in registry)+'\n')
