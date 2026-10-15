#!/usr/bin/env python3
"""Fixed-marker descriptive context; adds no inferential FDR family."""
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import run_flight_collection_interaction as r

MARKERS=['Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Alas2','Slc4a1','Gypa']


def main():
    ann=pd.read_csv(r.BASE/'inputs/gene_annotation.csv.gz')
    genes=pd.read_csv(r.OUT/'gene_flight_collection_interactions_exploratory.csv.gz')
    panels=pd.read_csv(r.OUT/'panel_flight_collection_interactions.csv')
    lookup=ann.set_index('ensembl_id').gene_symbol.to_dict()
    scores=[];avail=[];correlations=[];omissions=[];constituents=[]
    for n in ['561','562']:
        d=r.load_dataset(n);m=d['meta'];syms=np.array([lookup.get(g,'') for g in d['y'].index])
        symbol_scores={}
        for marker in MARKERS:
            ids=d['y'].index[syms==marker]
            avail.append({'accession':f'OSD-{n}','marker':marker,'n_retained_gene_ids':len(ids),
                          'retained_gene_ids':'|'.join(ids),'status':'retained' if len(ids) else 'not available after full-deposit mean-count filter'})
            if len(ids):symbol_scores[marker]=d['y'].loc[ids].mean(axis=0).to_numpy()
        if not symbol_scores:raise ValueError('No fixed markers available')
        blood=np.vstack(list(symbol_scores.values())).mean(0)
        pidx=r.indices(d['y'].index,ann)
        circadian=d['y'].iloc[pidx['Circadian']].mean(axis=0).to_numpy()
        x=r.design(m).to_numpy();resblood=blood-x@np.linalg.lstsq(x,blood,rcond=None)[0]
        rescirc=circadian-x@np.linalg.lstsq(x,circadian,rcond=None)[0]
        correlations.append({'accession':f'OSD-{n}','target':'Circadian_panel','indicator':'fixed erythroid-associated marker mean',
                             'n_animals':len(m),'raw_Spearman':stats.spearmanr(blood,circadian).statistic,
                             'full_design_residual_Pearson':stats.pearsonr(resblood,rescirc).statistic,
                             'interpretation':'descriptive correlation only; full-model residual correlation removes all design-cell means and additive age'})
        for marker,value in symbol_scores.items():
            other=np.vstack([v for s,v in symbol_scores.items() if s!=marker]).mean(0)
            rr=value-x@np.linalg.lstsq(x,value,rcond=None)[0]
            ro=other-x@np.linalg.lstsq(x,other,rcond=None)[0]
            correlations.append({'accession':f'OSD-{n}','target':marker,'indicator':'mean of other fixed markers; target omitted',
                                 'n_animals':len(m),'raw_Spearman':stats.spearmanr(other,value).statistic,
                                 'full_design_residual_Pearson':stats.pearsonr(ro,rr).statistic,
                                 'interpretation':'descriptive; target omitted to avoid mathematical self-correlation'})
        for j,meta in enumerate(m.itertuples()):
            row={'accession':f'OSD-{n}','animal':meta.animal,'flight':meta.flight,'cohort':meta.cohort,'age':meta.age,
                 'group':meta.flight+' '+meta.cohort.replace('ISS_T','ISS-T'),'n_available_marker_symbols':len(symbol_scores),
                 'blood_indicator_mean_log2':blood[j],'circadian_panel_mean_log2':circadian[j],
                 'blood_full_design_residual':resblood[j],'circadian_full_design_residual':rescirc[j]}
            row.update({s:v[j] for s,v in symbol_scores.items()});scores.append(row)
        for sensitivity in ['all_animals','exclude_any_region_metadata_conflict']:
            conflicts=set(json.loads((r.OUT/'analysis_summary.json').read_text())['conflicting_animals_either_region'])
            keep=np.ones(len(m),bool) if sensitivity=='all_animals' else ~m.animal.isin(conflicts).to_numpy()
            xx=x[keep];y=d['y'].loc[:,m.index[keep]].to_numpy()
            for panel,idx in pidx.items():
                full=np.linalg.lstsq(xx,y[idx].mean(axis=0),rcond=None)[0][-1]
                for gene in idx:
                    omit=idx[idx!=gene]
                    effect=np.linalg.lstsq(xx,y[omit].mean(axis=0),rcond=None)[0][-1]
                    omissions.append({'accession':f'OSD-{n}','sensitivity':sensitivity,'panel':panel,
                                      'omitted_ensembl_id':d['y'].index[gene],'omitted_gene_symbol':syms[gene],
                                      'full_interaction_log2':full,'interaction_without_gene_log2':effect})
            sub=genes[(genes.accession==f'OSD-{n}')&(genes.sensitivity==sensitivity)]
            for panel,idx in pidx.items():
                part=sub[sub.ensembl_id.isin(d['y'].index[idx])].copy();part['panel']=panel;constituents.append(part)
    sf=pd.DataFrame(scores);sf.to_csv(r.OUT/'fixed_blood_marker_animal_scores.csv',index=False)
    pd.DataFrame(avail).to_csv(r.OUT/'fixed_blood_marker_availability.csv',index=False)
    pd.DataFrame(correlations).to_csv(r.OUT/'fixed_blood_marker_descriptive_correlations.csv',index=False)
    pd.DataFrame(omissions).to_csv(r.OUT/'leave_one_gene_out_panel_interactions.csv',index=False)
    pd.concat(constituents).to_csv(r.OUT/'all_panel_gene_constituents.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    sf.groupby(['accession','flight','cohort','age']).agg(n_animals=('animal','size'),
        mean_indicator=('blood_indicator_mean_log2','mean'),median_indicator=('blood_indicator_mean_log2','median'),
        min_indicator=('blood_indicator_mean_log2','min'),max_indicator=('blood_indicator_mean_log2','max')).reset_index().to_csv(
            r.OUT/'fixed_blood_marker_group_distributions.csv',index=False)
    primary=genes[(genes.sensitivity=='all_animals')&(genes.q_BH_complete_gene_interaction_family<=.05)].copy()
    sensitivity=genes[genes.sensitivity!='all_animals'][['accession','ensembl_id','log2FoldChange','q_BH_complete_gene_interaction_family']]
    primary=primary.merge(sensitivity,on=['accession','ensembl_id'],suffixes=('','_metadata_exclusion'),validate='one_to_one')
    primary['same_NB_direction_after_exclusion']=np.sign(primary.log2FoldChange)==np.sign(primary.log2FoldChange_metadata_exclusion)
    primary['passes_original_global_FDR_in_both']=primary.q_BH_complete_gene_interaction_family_metadata_exclusion<=.05
    primary['known_fixed_erythroid_marker']=primary.gene_symbol.isin(MARKERS)
    primary.sort_values('q_BH_complete_gene_interaction_family').to_csv(r.OUT/'all_primary_NB_FDR_candidates_with_sensitivity.csv',index=False)
    plt.rcParams.update({'font.family':'Arial','font.size':16,'axes.titlesize':16,'axes.labelsize':16,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(14,11),constrained_layout=True)
    groups=['GC LAR','FLT LAR','GC ISS-T','FLT ISS-T'];colors=['#2878a0','#bd6b2d','#2878a0','#bd6b2d'];rng=np.random.default_rng(561562)
    for ax,accession,letter in [(axes[0,0],'OSD-561','a'),(axes[0,1],'OSD-562','b')]:
        for i,group in enumerate(groups):
            sub=sf[(sf.accession==accession)&(sf.group==group)]
            for age,offset,marker in [('12',-.08,'o'),('29',.08,'^')]:
                vals=sub[sub.age==age].blood_indicator_mean_log2
                ax.scatter(i+offset+rng.uniform(-.06,.06,len(vals)),vals,s=45,color=colors[i],marker=marker,
                           label=f'Age {age}' if i==0 else None)
            ax.hlines(sub.blood_indicator_mean_log2.median(),i-.24,i+.24,color='black',lw=2)
        ax.set_xticks(range(4),groups,rotation=15);ax.set_ylabel('Mean marker expression (log2)')
        region='Cerebellum' if accession=='OSD-561' else 'Hippocampus'
        ax.set_title(f'{letter}  {region}: fixed erythroid markers',loc='left',fontweight='bold',fontsize=16)
        ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15)
    ax=axes[1,0]
    for accession,offset,color in [('OSD-561',-.12,'#2878a0'),('OSD-562',.12,'#bd6b2d')]:
        sub=genes[(genes.accession==accession)&(genes.sensitivity=='all_animals')&genes.gene_symbol.isin(MARKERS)].set_index('gene_symbol').reindex(MARKERS)
        yy=np.arange(len(MARKERS))+offset
        ax.hlines(yy,sub.Wald_lower95,sub.Wald_upper95,color=color,lw=1.8);ax.scatter(sub.log2FoldChange,yy,s=38,color=color,label=accession)
        for i,row in enumerate(sub.itertuples()):
            if row.q_BH_complete_gene_interaction_family<=.05:ax.text(row.log2FoldChange,yy[i]-.19,'*',ha='center',fontsize=15)
    ax.set_yticks(range(len(MARKERS)),MARKERS);ax.invert_yaxis();ax.axvline(0,color='#777777',ls='--',lw=1)
    ax.set_xlabel('Interaction (log2), 95% Wald CI')
    ax.set_title('c  Primary count-model interactions\n* Original complete gene family q ≤ 0.05',loc='left',fontweight='bold',fontsize=16)
    ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
    ax=axes[1,1]
    for accession,color in [('OSD-561','#2878a0'),('OSD-562','#bd6b2d')]:
        sub=sf[sf.accession==accession]
        rho=stats.pearsonr(sub.blood_full_design_residual,sub.circadian_full_design_residual).statistic
        ax.scatter(sub.blood_full_design_residual,sub.circadian_full_design_residual,color=color,s=45,label=f'{accession}: residual r={rho:.2f}')
    ax.axhline(0,color='#777777',lw=.8);ax.axvline(0,color='#777777',lw=.8)
    ax.set_xlabel('Marker score: full-model residual');ax.set_ylabel('Circadian panel: full-model residual')
    ax.set_title('d  Residual association',loc='left',fontweight='bold',fontsize=16)
    ax.legend(frameon=False,loc='upper right');ax.spines[['top','right']].set_visible(False)
    for ext in ['png','svg']:fig.savefig(r.OUT/f'figure_fixed_blood_marker_context.{ext}',dpi=600,bbox_inches='tight')
    plt.close(fig)
    counts=primary.groupby('accession').size().to_dict()
    intersection=primary[primary.passes_original_global_FDR_in_both].groupby('accession').size().to_dict()
    p=panels[(panels.accession=='OSD-562')&(panels.panel=='Circadian')]
    c=p[p.sensitivity=='all_animals'].iloc[0];s=p[p.sensitivity!='all_animals'].iloc[0]
    text='# Biological context and qualified interpretation\n\n'
    text+=f'The primary hippocampal Circadian interaction is {c.interaction_log2:+.3f} log2 (95% animal-bootstrap interval {c.bootstrap_lower95:+.3f} to {c.bootstrap_upper95:+.3f}; approximate wild p={c.p_wild_approximate:.6f}; complete 32-test BH q={c.q_BH_complete_panel_family:.6f}). The Monte Carlo 95% interval for the wild exceedance probability is [{c.p_Monte_Carlo_lower95:.6f},{c.p_Monte_Carlo_upper95:.6f}]. The count-model mean for the same retained panel genes is {c.mean_NB_interaction_log2:+.3f} log2. The metadata sensitivity retains the negative direction ({s.interaction_log2:+.3f}) but does not pass correction (q={s.q_BH_complete_panel_family:.3f}); only one older FLT ISS-T animal remains. This is a qualified cohort-dependent expression candidate, not a robust flight-clock mechanism.\n\n'
    text+=f'Primary gene-model discoveries number {counts.get("OSD-561",0)} in cerebellum and {counts.get("OSD-562",0)} in hippocampus under the original global gene family. Respectively {intersection.get("OSD-561",0)} and {intersection.get("OSD-562",0)} retain global q ≤ 0.05 after exclusions. The much larger exclusion-only cerebellar set is not promoted: that sensitivity has no older FLT ISS-T animals and depends on additive-age extrapolation. The complete candidate table retains every primary discovery, including sensitivity failures.\n\n'
    text+='Hbb-bt and Hbb-bs are supported count-model interactions in full primary cerebellum and keep their direction and global significance after exclusion. Both are beta-globin transcripts; together with the fixed markers Hba-a1, Hba-a2, Alas2, Slc4a1 and Gypa, they motivate a blood/erythroid-associated composition interpretation as one possibility. All seven markers are retained in both regions and their availability, animal distributions, and existing global gene-model inference are exported. Their interaction point estimates are all negative in each primary region, although most do not pass the original global FDR. No new marker-set significance test or favorable FDR family was created. Marker expression cannot prove blood contamination, assign a cell source, or quantify blood fraction; hemoglobin-associated expression and collection-linked biology can also contribute. Correlations use target-omitted marker scores and separately remove the entire full design to avoid simple self-correlation or between-cell means.\n\n'
    text+='Dbp, a circadian output transcription factor, has a negative cerebellar gene interaction passing the original global FDR in full and excluded data. Its hippocampal primary gene q is just above 0.05; a significant aggregate Circadian panel does not imply that every clock gene is individually significant. Other supported primary hippocampal gene candidates include Ccbe1 and Dmkn, which retain direction and original global significance after exclusion; these heterogeneous candidates do not establish a shared neuronal mechanism. All gene inference remains exploratory NB Wald inference; the cohort comparison bundles collection history and other factors and cannot identify microgravity specifically.\n'
    (r.OUT/'biological_context_summary.md').write_text(text)
    print('Fixed markers available:',len(avail),'rows; primary gene discoveries:',counts,'retained after exclusions:',intersection)


if __name__=='__main__':main()
