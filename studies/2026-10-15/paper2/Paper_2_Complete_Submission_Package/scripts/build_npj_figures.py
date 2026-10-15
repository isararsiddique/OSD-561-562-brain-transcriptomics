"""Final-size journal figures directly from complete, frozen analysis tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
BASE=Path(__file__).resolve().parents[1];R=BASE/'results';OUT=BASE/'figures/npj';OUT.mkdir(exist_ok=True)
W=183/25.4
font=font_manager.findfont('Arial',fallback_to_default=False)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':1,'lines.linewidth':1.2,'patch.linewidth':1,'xtick.major.width':1,'ytick.major.width':1,'grid.linewidth':1,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
C={'FLT':'#0072B2','GC':'#D55E00','VIV':'#7356A8'};ORDER=['FLT_minus_GC','GC_minus_VIV','FLT_minus_VIV'];LABEL=['FLT–GC','GC–VIV','FLT–VIV'];registry=[]
def panel(ax,letter,title):ax.set_title(letter+')  '+title,loc='left',fontweight='bold',pad=10)
def save(fig,name,caption,sources):
    fig.canvas.draw()
    for ext in ['pdf','svg','png']:fig.savefig(OUT/(name+'.'+ext),dpi=600)
    # Check actual physical size and minimum non-empty text, not raster DPI alone.
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible()]
    assert all(t.get_fontsize()>=8 for t in texts),(name,min(t.get_fontsize() for t in texts))
    registry.append(dict(name=name,width_mm=float(fig.get_figwidth()*25.4),height_mm=float(fig.get_figheight()*25.4),min_font_pt=min(t.get_fontsize() for t in texts),font='Arial',dpi=600,sources=sources,caption=caption))
    (OUT/(name+'_caption.txt')).write_text(caption+'\n');plt.close(fig)
def group_legend(ax,**kw):
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=C[g],label=g,ms=4) for g in C],frameon=False,**kw)
meta=pd.read_csv(BASE/'inputs/animal_metadata.csv').set_index('sample');qc=pd.read_csv(R/'animal_qc.csv');frac=pd.read_csv(R/'all26_animal_fractions.csv.gz',index_col=0)
glob=pd.read_csv(R/'global_contrasts.csv').set_index('contrast').loc[ORDER];pca=pd.read_csv(R/'unsupervised_animal_pca.csv',index_col=0);desc=json.loads((R/'descriptive_summary.json').read_text())
fig,axs=plt.subplots(2,2,figsize=(W,5.55));fig.subplots_adjust(left=.09,right=.98,bottom=.09,top=.94,wspace=.37,hspace=.57)
a=axs[0,0];panel(a,'a','Source animals and collection');a.axis('off')
for k,(g,n,l,date) in enumerate([('FLT',9,10,'14 Jan 2021'),('GC',9,10,'17 Jan 2021'),('VIV',8,9,'17 Jan 2021')]):
 y=.83-k*.24;a.scatter(.02,y,s=42,color=C[g],transform=a.transAxes);a.text(.09,y,f'{g}: {n} animals, {l} libraries\nCollected {date}',va='center',transform=a.transAxes,linespacing=1.5)
a.text(.0,.015,'Technical libraries summed within animal.\nFlight and collection date are aliased.',transform=a.transAxes,linespacing=1.5)
a=axs[0,1];panel(a,'b','All measured binary choices')
for g in C:
 d=pca[pca.group==g];a.scatter(d.PC1,d.PC2,s=22,c=C[g],label=g)
a.set_xlabel(f"PC1 ({100*desc['PCA_PC1_variance_fraction']:.1f}%)");a.set_ylabel(f"PC2 ({100*desc['PCA_PC2_variance_fraction']:.1f}%)");group_legend(a,loc='upper right')
a=axs[1,0];panel(a,'c','Retained junction read depth')
for k,g in enumerate(C):
 y=qc.loc[qc.group==g,'retained_junction_depth'].to_numpy()/1e6;a.scatter(k+np.linspace(-.12,.12,len(y)),y,s=20,color=C[g]);a.hlines(y.mean(),k-.23,k+.23,color='black',lw=1.2)
a.set_xticks(range(3),C);a.set_ylabel('Unique junction reads (millions)');a.set_ylim(0,50)
a=axs[1,1];panel(a,'d','Global binary contrast')
for k,row in enumerate(glob.itertuples()):
 a.scatter(row.RMS_contrast_pp,k,color=list(C.values())[k],s=32);a.text(.15,k-.19,f'Holm p = {row.Holm_three_global_tests:.6f}',fontsize=8)
a.set_yticks(range(3),LABEL);a.set_xlim(0,6.5);a.set_ylim(2.5,-.6);a.set_xlabel('RMS fraction contrast (percentage points)')
save(fig,'Figure_1','Figure 1. Animal-level design and aggregate junction differences. a) RR23 cerebellar source animals and collection dates; technical libraries are not biological replicates. b) Principal components of all 5,249 measured binary choices, centered across 26 animals without selection by test results. c) All animal retained unique-junction depths; segments show group means. d) RMS differences in group mean fractions over the same 5,249 choices. Exact whole-animal allocations enumerate 48,620 FLT–GC and 24,310 assignments in each 9-versus-8 comparison, with Holm correction across three global tests. These are effect summaries, not bars hiding biological variation; the animal distributions appear in b and c. FLT, flight; GC, habitat ground control; VIV, vivarium control. FLT n=9, GC n=9, VIV n=8. Collection-date aliasing limits causal interpretation.', ['results/animal_qc.csv','results/unsupervised_animal_pca.csv','results/global_contrasts.csv'])
gene=pd.read_csv(BASE/'extensions/gene_omnibus/results/gene_omnibus_complete.csv.gz');fg=gene[gene.contrast==ORDER[0]];selected=pd.read_csv(BASE/'extensions/gene_omnibus/results/figure_source_events.csv');loo=pd.read_csv(BASE/'extensions/gene_omnibus/results/descriptive_gene_leave_one_out_summary.csv')
fig,axs=plt.subplots(2,2,figsize=(W,5.8));fig.subplots_adjust(left=.09,right=.98,bottom=.09,top=.94,wspace=.38,hspace=.58)
a=axs[0,0];panel(a,'a','All 1,983 measured genes');hit=fg.q_BY_joint_3G<.05
for mask,col,size in [(~hit,'#A8B3B8',7),(hit,C['FLT'],16)]:a.scatter(fg.loc[mask,'rms_delta_percentage_points'],-np.log10(fg.loc[mask,'q_BY_joint_3G']),s=size,c=col,alpha=.8,rasterized=True)
a.axhline(-np.log10(.05),c='black',ls='--',lw=1);a.set_xlabel('Gene RMS contrast (percentage points)');a.set_ylabel('−log10 joint BY q');a.set_ylim(-.05,1.65);a.text(.96,.93,'47 BY calls\n5,949 joint tests',ha='right',va='top',transform=a.transAxes)
def choices(a,symbol,letter):
 panel(a,letter,symbol+' local choices');ev=selected[selected.gene_symbol==symbol]
 for j,e in enumerate(ev.event_id):
  for k,g in enumerate(C):
   y=100*frac.loc[e,meta.index[meta.group==g]].to_numpy();x=j*4+k;a.scatter(x+np.linspace(-.1,.1,len(y)),y,s=12,c=C[g],alpha=.85);a.hlines(y.mean(),x-.2,x+.2,color='black',lw=1)
 a.set_xticks(np.arange(len(ev))*4+1,[str(j+1) for j in range(len(ev))]);a.set_xlabel('Measured choice');a.set_ylabel('Junction A fraction (%)');a.set_ylim(0,100);group_legend(a,ncol=3,loc='lower center',columnspacing=.7,handletextpad=.25)
choices(axs[0,1],'Hnrnpa2b1','b');choices(axs[1,0],'Stxbp1','c')
a=axs[1,1];panel(a,'d','Primary gene measurement sensitivity')
s=pd.read_csv(BASE/'extensions/gene_omnibus/results/primary_BY_gene_dominance95_followup.csv');print('Gene sensitivity schema',s.columns.tolist())
# Display complete primary list ranked by magnitude, with measured status from fixed-family result.
sens=pd.read_csv(BASE/'extensions/gene_omnibus/results/dominance95_gene_omnibus_complete.csv.gz');m=fg.loc[hit].merge(sens[sens.contrast==ORDER[0]][['gene_id','q_BY_joint_3G','rms_delta_percentage_points']],on='gene_id',suffixes=('_primary','_95'))
a.scatter(m.rms_delta_percentage_points_primary,m.rms_delta_percentage_points_95,c=np.where(m.q_BY_joint_3G_95<.05,C['FLT'],C['GC']),s=22);a.plot([0,24],[0,24],c='black',ls='--',lw=1);a.set_xlim(0,24);a.set_ylim(0,24);a.set_xlabel('Primary gene RMS (pp)');a.set_ylabel('95% dominance gene RMS (pp)');a.text(.04,.96,'45/47 retain joint BY q < 0.05',va='top',transform=a.transAxes);a.text(.96,.05,'Ecpas unmeasured',ha='right',transform=a.transAxes)
save(fig,'Figure_2','Figure 2. Conservative gene localization and animal-level measurements. a) Every measured gene in the FLT–GC comparison; blue points pass joint Benjamini–Yekutieli correction across 5,949 gene–contrast tests. All 47 calls have exact p=2/48,620 and q=0.04825699; the dashed line is q=0.05. b,c) Every measured choice and all animals for post hoc illustrative genes Hnrnpa2b1 and Stxbp1. Each dot is one animal and each segment a group mean; FLT n=9, GC n=9, VIV n=8. Choice labels follow the deposited exact-coordinate table and do not denote transcript isoforms. d) The 46 primary genes remaining measurable under 95% dominance; blue denotes retained BY significance, orange denotes lost significance. Ecpas loses its only choice and is unmeasured rather than assigned a zero effect. The original 5,949-test family is retained, with p=1 where measurements are lost. Diagonal is equal magnitude. This sensitivity does not replace the original candidate list or constitute replication.', ['extensions/gene_omnibus/results/gene_omnibus_complete.csv.gz','extensions/gene_omnibus/results/figure_source_events.csv','results/all26_animal_fractions.csv.gz','extensions/gene_omnibus/results/dominance95_gene_omnibus_complete.csv.gz'])
multi=pd.read_csv(BASE/'extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz');mg=pd.read_csv(BASE/'extensions/multiway_splicing/tables/global_exact_results.csv').set_index('contrast');comp=pd.read_csv(BASE/'extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz');comp=comp.merge(meta.reset_index()[['sample','group']],on='sample');tc=comp[comp.cluster_id=='6:2:42663893:start']
fig=plt.figure(figsize=(W,5.45));gs=fig.add_gridspec(2,2,left=.09,right=.98,bottom=.1,top=.94,wspace=.38,hspace=.66);a=fig.add_subplot(gs[0,0]);panel(a,'a','All 923 local clusters')
m=multi[multi.contrast=='FLT_vs_GC'];hit=m.BY_q_full_3E<.05
for mask,col,size in [(~hit,'#A8B3B8',8),(hit,C['FLT'],17)]:a.scatter(100*m.loc[mask,'max_abs_proportion_difference'],-np.log10(m.loc[mask,'BY_q_full_3E']),s=size,c=col,alpha=.8,rasterized=True)
a.axhline(-np.log10(.05),c='black',ls='--',lw=1);a.set_xlabel('Largest alternative difference (pp)');a.set_ylabel('−log10 joint BY q');a.set_ylim(-.05,1.85);a.text(.97,.95,'FLT–GC\n40 BY calls',ha='right',va='top',transform=a.transAxes)
a=fig.add_subplot(gs[0,1]);panel(a,'b','Global multialternative contrasts')
for k,key in enumerate(['FLT_vs_GC','GC_vs_VIV','FLT_vs_VIV']):
 r=mg.loc[key];a.scatter(r.hellinger_centroid_distance,k,c=list(C.values())[k],s=32);a.text(.003,k-.20,f'Holm p = {r.holm_p_3:.6f}')
a.set_yticks(range(3),LABEL);a.set_xlim(0,.10);a.set_ylim(2.5,-.6);a.set_xlabel('Scaled centroid distance')
a=fig.add_subplot(gs[1,:]);panel(a,'c','Tcaf1 read support at the shared acceptor')
for j,junction in enumerate(sorted(tc.junction_id.unique())):
 for k,g in enumerate(C):
  y=100*tc[(tc.junction_id==junction)&(tc.group==g)].proportion.to_numpy();x=j*4+k;a.scatter(x+np.linspace(-.11,.11,len(y)),y,s=18,c=C[g]);a.hlines(y.mean(),x-.22,x+.22,color='black',lw=1.2)
a.set_xticks(np.arange(3)*4+1,['42,669,866','42,686,462','42,686,892']);a.set_xlabel('Donor-side intron end on chromosome 6 (GRCm39, minus strand)');a.set_ylabel('Read fraction (%)');a.set_ylim(0,100);group_legend(a,ncol=3,loc='upper left')
save(fig,'Figure_3','Figure 3. Multialternative compositions and a coordinate-defined candidate. a) Every FLT–GC cluster in the joint 2,769-test family; blue denotes BY q<0.05 and the dashed line q=0.05. b) Global equal-cluster contrasts from square-root animal compositions; scaled distance is sqrt(T/2), not a percentage or the Hellinger distance between pooled group means. Exact whole-animal allocations and three-test Holm correction are as in Fig. 1. c) All 26 animals at Tcaf1 cluster 6:2:42663893:start; dots are individual animals and segments group means (FLT n=9, GC n=9, VIV n=8). Coordinates denote inclusive STAR intron ends, not exon coordinates; this minus-strand boundary shares an acceptor and varies donor-side support. The native-composition changes are descriptive components of the corrected cluster test, not separately tested alternatives. Tcaf1 passes the primary FLT–GC and FLT–VIV BY tests, but its FLT–VIV result does not survive the 95% dominance sensitivity. The selection of Tcaf1 for illustration is post hoc.', ['extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz','extensions/multiway_splicing/tables/global_exact_results.csv','extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz'])
# Supplementary Figure 1: all 47 candidates, no illustrative subset.
g=fg[fg.q_BY_joint_3G<.05].sort_values('gene_symbol');l=pd.read_csv(BASE/'extensions/gene_omnibus/results/descriptive_all_gene_leave_one_out.csv.gz');print('LOO schema',l.columns.tolist())
fig,axs=plt.subplots(1,2,figsize=(W,6.5));fig.subplots_adjust(left=.14,right=.98,bottom=.09,top=.94,wspace=.65)
for a,chunk,letter in zip(axs,[g.iloc[:24],g.iloc[24:]],['a','b']):
 panel(a,letter,'Primary FLT–GC candidates')
 for k,r in enumerate(chunk.itertuples()):a.scatter(r.rms_delta_percentage_points,k,c=C['FLT'],s=20)
 a.set_yticks(range(len(chunk)),chunk.gene_symbol);a.invert_yaxis();a.set_xlim(0,24);a.set_xlabel('Gene RMS contrast (pp)');a.grid(axis='x',alpha=.12)
save(fig,'Supplementary_Figure_1','Supplementary Figure 1. All 47 primary gene candidates. a,b) The complete alphabetic list of FLT–GC genes passing joint BY correction across all 5,949 tests, with RMS effect magnitude across every measured choice assigned to the gene. These are gene-level effect summaries, not individual-animal distributions or confidence intervals. All calls have exact p=2/48,620 and joint BY q=0.04825699. Source-animal measurements, event counts and all other genes and contrasts remain in Supplementary Data 3; the 95% gate retains 45 of these 47 calls.', ['extensions/gene_omnibus/results/gene_omnibus_complete.csv.gz'])
fig,axs=plt.subplots(1,2,figsize=(W,3.25));fig.subplots_adjust(left=.09,right=.98,bottom=.21,top=.88,wspace=.4)
a=axs[0];panel(a,'a','Shared-control dependence');spl=pd.read_csv(R/'disjoint_control_profile_correlations.csv');pr=json.loads((R/'profile_dependence.json').read_text());a.hist(spl.rho,bins=14,color=C['GC'],alpha=.8);a.axvline(pr['shared_control_rho'],c=C['FLT'],lw=1.5,label='Shared GC');a.axvline(pr['disjoint_rho_median'],c='black',ls='--',lw=1,label='Disjoint median');a.set_xlabel('Effect-profile Spearman correlation');a.set_ylabel('Number of GC partitions');a.legend(frameon=False)
a=axs[1];panel(a,'b','Full-matrix depth thinning');thin=pd.read_csv(R/'depth_thinning_global_sensitivity.csv')
for k,c in enumerate(ORDER):
 d=thin[thin.contrast==c];a.plot(d.replicate,d.RMS_contrast_pp,'o-',color=list(C.values())[k],ms=3,label=LABEL[k]);a.axhline(glob.loc[c,'RMS_contrast_pp'],c=list(C.values())[k],ls='--',lw=1)
a.set_xlabel('Thinning replicate');a.set_ylabel('Global binary RMS (pp)');a.legend(frameon=False,loc='center right')
save(fig,'Supplementary_Figure_2','Supplementary Figure 2. Dependence and read-depth diagnostics for binary choices. a) Effect-profile correlations for all 126 distinct partitions of the nine habitat controls into disjoint sets of four and five. The blue line uses the same controls in both contrasts; the dashed line is the median over disjoint partitions. Partitions reuse animals and are not independent replicates. b) All ten binomial depth thinnings of the full count matrix, with frozen binary definitions; dashed lines show unthinned effects. Target retained depth is 16,249,784 reads. All 30 global tests retain Holm significance. These diagnostics assess mathematical dependence and read sampling; they do not provide new biological replication or remove collection-date aliasing.', ['results/disjoint_control_profile_correlations.csv','results/depth_thinning_global_sensitivity.csv'])
ab=pd.read_csv(BASE/'extensions/abundance_integration/results/abundance_junction_complete.csv.gz');af=ab[ab.contrast==ORDER[0]]
fig,axs=plt.subplots(1,2,figsize=(W,3.3));fig.subplots_adjust(left=.09,right=.98,bottom=.22,top=.88,wspace=.40)
a=axs[0];panel(a,'a','All 1,977 jointly measured genes');a.scatter(af.rms_delta_percentage_points,np.abs(af.mean_log2_difference),s=8,c='#829DA7',alpha=.6,rasterized=True);a.set_xlabel('Junction gene RMS (pp)');a.set_ylabel('Absolute log2 abundance contrast');a.text(.96,.96,'Spearman ρ = 0.0572',ha='right',va='top',transform=a.transAxes)
a=axs[1];panel(a,'b','Abundance equivalence conjunction')
for k,fold in enumerate(['1.2','1.5']):
 y=af['q_BY_conjunction_'+fold].to_numpy();a.scatter(np.arange(len(y)),np.sort(y),s=6,color=[C['FLT'],C['GC']][k],label=fold+'-fold margin',alpha=.65,rasterized=True)
a.axhline(.05,color='black',ls='--',lw=1);a.set_ylim(0,1.05);a.set_xlabel('Genes ranked within each margin');a.set_ylabel('Joint BY conjunction q');a.legend(frameon=False,loc='center')
save(fig,'Supplementary_Figure_3','Supplementary Figure 3. Gene abundance does not establish an independent junction phenotype. a) FLT–GC effects for all 1,977 jointly measurable genes; the rank correlation is descriptive and has no independent-gene significance test. b) All FLT–GC junction-association/abundance-equivalence conjunction q values, sorted separately at each margin. The conjunction p is the maximum of the exact junction p and the Welch equivalence p; BY correction covers 5,931 tests separately for each margin. No conjunction passes BY at either margin. These points summarize gene-level effects and corrected tests; they are not biological replicates. Relative abundance used the publisher-verified STAR matrix, technical-library aggregation, median-ratio normalization and log2(normalized count+0.5).', ['extensions/abundance_integration/results/abundance_junction_complete.csv.gz'])
v=BASE/'extensions/independent_validation/tables';effects=pd.read_csv(v/'complete_common_panel_effects.csv.gz');scores=pd.read_csv(v/'RRRM2_heldout_animal_scores.csv');null=pd.read_csv(v/'RRRM2_exact_score_null.csv.gz');summary=json.loads((v/'summary.json').read_text());print('Null schema',null.columns.tolist())
fig,axs=plt.subplots(1,2,figsize=(W,3.6));fig.subplots_adjust(left=.1,right=.98,bottom=.25,top=.88,wspace=.42)
a=axs[0];panel(a,'a','All 1,512 measurable choices');a.scatter(100*effects.RR23_delta_FLT_minus_GC,100*effects.RRRM2_delta_stratified_FLT_minus_GC,s=7,c='#829DA7',alpha=.5,rasterized=True);a.axhline(0,c='black',lw=1);a.axvline(0,c='black',lw=1);a.set_xlabel('RR23 FLT–GC fraction difference (pp)');a.set_ylabel('RRRM2 stratified difference (pp)');a.text(.03,.96,'ρ = 0.0102',va='top',transform=a.transAxes)
a=axs[1];panel(a,'b','Frozen score in held-out animals')
strata=sorted(scores.stratum.unique())
for j,s in enumerate(strata):
 for k,g in enumerate(['FLT','GC']):
  y=100*scores[(scores.stratum==s)&(scores.exposure==g)].learned_direction_score.to_numpy();x=j*3+k;a.scatter(x+np.linspace(-.08,.08,len(y)),y,s=24,c=C[g]);a.hlines(y.mean(),x-.2,x+.2,color='black',lw=1)
a.set_xticks(np.arange(len(strata))*3+.5,[s.replace('ISS_T_','ISS T / ').replace('LAR_','LAR / ') for s in strata]);a.set_xlabel('Route and design age stratum');a.set_ylabel('RR23 direction score (pp)');a.legend(handles=[Line2D([],[],marker='o',ls='',c=C[g],label=g) for g in ['FLT','GC']],frameon=False,loc='upper right')
save(fig,'Supplementary_Figure_4','Supplementary Figure 4. Independent-mission transport of the frozen binary profile. a) Complete measurable overlap between RR23 and RRRM2; stratum-adjusted RRRM2 effects are compared with RR23 FLT–GC effects. Spearman ρ=0.01019, with two-cohort animal-bootstrap 95% interval −0.06894 to 0.07386. b) All 19 RRRM2 animals in the three eligible strata, coloured by exposure; segments are group means. The score uses RR23-only signs and GC centering, averaged over the same 1,512 choices. Exact allocation of labels within strata yields 11,200 assignments and one-sided p=0.2525. The stratified score effect is 0.1281 pp, with frozen-weight animal-bootstrap 95% interval −0.1372 to 0.4218 pp. No common binary choice passes full-family directional Holm correction. The disputed ISS_T age-29 stratum was excluded in full.', ['extensions/independent_validation/tables/complete_common_panel_effects.csv.gz','extensions/independent_validation/tables/RRRM2_heldout_animal_scores.csv'])
(OUT/'figure_registry.json').write_text(json.dumps(registry,indent=2));print('Built',len(registry),'final-size figures.')
