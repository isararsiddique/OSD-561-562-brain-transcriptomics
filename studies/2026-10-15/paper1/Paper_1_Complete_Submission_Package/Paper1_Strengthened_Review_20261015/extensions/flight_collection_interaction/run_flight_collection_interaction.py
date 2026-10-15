#!/usr/bin/env python3
"""Within-region flight x collection-cohort interactions on collapsed real data.

The direct coefficient is (FLT-GC in ISS-T) - (FLT-GC in LAR).
It is observational effect modification, not a randomized causal flight effect.
"""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / 'scripts'))
from run_revision import load_dataset, indices, PANELS, check_input_integrity
from statistics_core import bh

SEED = 5615622026
COLS = ['Intercept', 'age29', 'flight', 'isst', 'flight_isst']
DEFINITION = '(FLT-GC in ISS-T) minus (FLT-GC in LAR)'


def design(meta):
    if not set(meta.age).issubset({'12', '29'}):
        raise ValueError('Unexpected age coding')
    f = (meta.flight == 'FLT').to_numpy(float)
    c = (meta.cohort == 'ISS_T').to_numpy(float)
    return pd.DataFrame(np.column_stack([np.ones(len(meta)), (meta.age == '29').to_numpy(float), f, c, f*c]),
                        index=meta.index, columns=COLS)


def hc3(y, x):
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Non-identifiable interaction design')
    inv = np.linalg.inv(x.T @ x)
    coefficients = (inv @ x.T @ y.T).T
    residual = y - coefficients @ x.T
    h = np.einsum('ij,ji->i', x, inv @ x.T)
    if np.any(h >= 1-1e-9):
        raise ValueError('Unit-leverage animal: HC3 uncertainty unavailable')
    w = (x @ inv[:, -1]) / (1-h)
    se = np.sqrt(residual**2 @ w**2)
    t = np.divide(coefficients[:, -1], se, out=np.zeros(len(y)), where=se>0)
    return coefficients[:, -1], se, t, coefficients, residual, h


def panel_inference(y, meta, nboot, nwild, seed):
    x = design(meta).to_numpy()
    beta, se, t, coefs, _, _ = hc3(y, x)
    x0 = x[:, :-1]
    fit0 = (np.linalg.lstsq(x0, y.T, rcond=None)[0].T) @ x0.T
    h0 = np.einsum('ij,ji->i', x0, np.linalg.pinv(x0))
    r0 = (y-fit0) / np.sqrt(1-h0)
    rng = np.random.default_rng(seed)
    exceed = np.zeros(len(y), int)
    inv=np.linalg.inv(x.T@x);mapping=inv@x.T
    h=np.einsum('ij,ji->i',x,mapping)
    w=(x@inv[:,-1])/(1-h)
    for start in range(0,nwild,2000):
        # One weight per collapsed animal, shared by every panel in this region.
        n=min(2000,nwild-start)
        weights=rng.choice([-1.,1.],size=(n,len(meta)))
        simulated=fit0[None,:,:]+r0[None,:,:]*weights[:,None,:]
        cs=simulated@mapping.T
        rs=simulated-cs@x.T
        ses=np.sqrt((rs**2)@(w**2))
        ts=np.divide(cs[:,:,-1],ses,out=np.zeros_like(ses),where=ses>0)
        exceed+=(abs(ts)>=abs(t)[None,:]-1e-12).sum(0)
    p = (exceed+1)/(nwild+1)
    mc_lower=np.where(exceed==0,0,stats.beta.ppf(.025,exceed,nwild-exceed+1))
    mc_upper=np.where(exceed==nwild,1,stats.beta.ppf(.975,exceed+1,nwild-exceed))
    frame = meta.reset_index(drop=True)
    strata = [s.index.to_numpy() for _, s in frame.groupby(['flight', 'cohort', 'age'], sort=True)]
    boots = []
    for _ in range(nboot):
        chosen = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
        xb = x[chosen]
        if np.linalg.matrix_rank(xb) == x.shape[1]:
            boots.append(np.linalg.lstsq(xb, y[:, chosen].T, rcond=None)[0][-1])
    if len(boots) < .95*nboot:
        raise ValueError('Too many bootstrap rank failures')
    boots = np.asarray(boots)
    return beta, se, p, np.quantile(boots, [.025, .975], axis=0), coefs, len(boots), mc_lower, mc_upper


def counts_fit(dataset, meta, sensitivity, force):
    """Refit actual integer counts with an explicit numeric full design."""
    accession = meta.accession.iloc[0]
    prefix = OUT / f'{accession}_{sensitivity}'
    destination = Path(str(prefix)+'_NB_interaction.csv.gz')
    diag_path = Path(str(prefix)+'_NB_fit.json')
    x = design(meta)
    counts = dataset['counts'].loc[dataset['y'].index, meta.index].T
    if not np.allclose(counts, np.round(counts)) or (counts<0).any().any():
        raise ValueError('Non-integer or negative counts')
    counts = counts.astype(np.int64)
    digest = hashlib.sha256(counts.to_numpy().tobytes()+x.to_numpy().tobytes()+
                            '|'.join(counts.index).encode()+'|'.join(counts.columns).encode()).hexdigest()
    if destination.exists() and diag_path.exists() and not force:
        diag = json.loads(diag_path.read_text())
        if diag['input_design_sha256'] != digest:
            raise ValueError('Cached fit input/design hash mismatch')
        print(f'{accession} {sensitivity}: verified cached NB fit', flush=True)
        return pd.read_csv(destination), diag
    print(f'{accession} {sensitivity}: NB fit, {len(meta)} animals x {counts.shape[1]} genes', flush=True)
    import pydeseq2
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    dds = DeseqDataSet(counts=counts, metadata=meta, design=x, n_cpus=1, quiet=True,
                       low_memory=True, refit_cooks=True)
    dds.deseq2()
    contrast = np.array([0.,0.,0.,0.,1.])
    st = DeseqStats(dds, contrast=contrast, independent_filter=False, cooks_filter=True,
                    n_cpus=1, quiet=True)
    st.summary()
    frame = st.results_df.rename_axis('ensembl_id').reset_index().drop(columns=['padj'])
    lfc = np.asarray(dds.varm['LFC'])/np.log(2)
    if lfc.shape[1] != len(COLS):
        raise ValueError('Unexpected NB coefficient order')
    assert np.allclose(frame.log2FoldChange, lfc[:, -1], equal_nan=True)
    frame['FLT_minus_GC_LAR_log2'] = lfc[:, 2]
    frame['FLT_minus_GC_ISST_log2'] = lfc[:, 2]+lfc[:, 4]
    assert np.allclose(frame.log2FoldChange, frame.FLT_minus_GC_ISST_log2-frame.FLT_minus_GC_LAR_log2, equal_nan=True)
    frame['Wald_lower95'] = frame.log2FoldChange-1.96*frame.lfcSE
    frame['Wald_upper95'] = frame.log2FoldChange+1.96*frame.lfcSE
    frame['accession'], frame['sensitivity'] = accession, sensitivity
    frame.to_csv(destination, index=False, compression={'method':'gzip','mtime':0})
    pd.DataFrame(lfc, index=counts.columns, columns=COLS).rename_axis('ensembl_id').to_csv(
        Path(str(prefix)+'_NB_coefficients.csv.gz'), compression={'method':'gzip','mtime':0})
    x.assign(animal=meta.animal.to_numpy()).to_csv(Path(str(prefix)+'_design.csv'))
    diag = {'accession':accession,'sensitivity':sensitivity,'input_design_sha256':digest,
            'design_columns':COLS,'numeric_contrast':contrast.tolist(),'coefficient_definition':DEFINITION,
            'n_animals':len(meta),'n_genes':counts.shape[1], 'n_finite_p':int(frame.pvalue.notna().sum()),
            'n_finite_dispersions':int(np.isfinite(dds.var.dispersions).sum()),
            'size_factor_range':[float(dds.obs.size_factors.min()),float(dds.obs.size_factors.max())],
            'pydeseq2':pydeseq2.__version__, 'independent_filter':False,'cooks_filter':True,
            'refit_cooks':True,'coefficient_shrinkage':False,'fit_type':str(getattr(dds,'fit_type','parametric'))}
    diag_path.write_text(json.dumps(diag,indent=2)+'\n')
    print(f'{accession} {sensitivity}: NB complete',flush=True)
    return frame, diag


def cached_contrast_check(frame, accession):
    """Existing full cohort contrasts are an independent parameterization."""
    path = BASE / f'results/count_models/{accession}_gene_results.csv.gz'
    if not path.exists():
        return {'accession':accession,'status':'no existing reference contrast file'}
    old = pd.read_csv(path)
    label = 'OSD561_cerebellum' if accession=='OSD-561' else 'OSD562_hippocampus'
    a = old[old.contrast==f'{label}_FLT_vs_GC_within_ISS_T'].set_index('ensembl_id').log2FoldChange
    b = old[old.contrast==f'{label}_FLT_vs_GC_within_LAR'].set_index('ensembl_id').log2FoldChange
    expected = a-b
    current = frame.set_index('ensembl_id').log2FoldChange
    joined = pd.concat([current.rename('direct'),expected.rename('cached_difference')],axis=1).dropna()
    difference = joined.direct-joined.cached_difference
    return {'accession':accession,'n_compared':len(joined),
            'Spearman_direct_vs_cached_difference':float(stats.spearmanr(joined.direct,joined.cached_difference).statistic),
            'median_absolute_difference':float(abs(difference).median()),'max_absolute_difference':float(abs(difference).max()),
            'comparison':'new explicit interaction vs cached age + four-level condition model; coefficient difference only, no subtraction of p-values'}


def figure(panels, genes, support):
    support=support.copy();support["age"]=support.age.astype(str)
    plt.rcParams.update({'font.family':'Arial','font.size':16,'axes.titlesize':16,'axes.labelsize':16,'svg.fonttype':'none'})
    fig, axes = plt.subplots(2,2,figsize=(14,11),constrained_layout=True)
    labels = [p.replace('_',' / ') for p in PANELS]
    colors = {'all_animals':'#2878a0','exclude_any_region_metadata_conflict':'#bd6b2d'}
    for ax,accession,letter in [(axes[0,0],'OSD-561','a'),(axes[0,1],'OSD-562','b')]:
        for sensitivity,offset in [('all_animals',-.13),('exclude_any_region_metadata_conflict',.13)]:
            sub=panels[(panels.accession==accession)&(panels.sensitivity==sensitivity)].set_index('panel').loc[list(PANELS)]
            y=np.arange(len(sub))+offset
            ax.hlines(y,sub.bootstrap_lower95,sub.bootstrap_upper95,color=colors[sensitivity],lw=2)
            ax.scatter(sub.interaction_log2,y,color=colors[sensitivity],s=40,
                       label='All animals' if sensitivity=='all_animals' else 'Metadata exclusions')
            for i,row in enumerate(sub.itertuples()):
                if row.q_BH_complete_panel_family<=.05:
                    ax.text(row.interaction_log2,y[i]-.17,'*',ha='center',fontsize=15)
        ax.set_yticks(np.arange(len(PANELS)),labels);ax.invert_yaxis();ax.axvline(0,color='#777777',ls='--',lw=1)
        ax.set_xlabel('Interaction (log2), 95% animal bootstrap CI')
        region='Cerebellum' if accession=='OSD-561' else 'Hippocampus'
        qmin=panels[panels.accession==accession].q_BH_complete_panel_family.min()
        ax.set_title(f'{letter}  {region}: direct interaction\nMinimum q={qmin:.3f}; complete family of {len(panels)} tests',loc='left',fontweight='bold')
        ax.grid(axis='x',alpha=.15);ax.spines[['top','right']].set_visible(False)
    axes[0,0].legend(frameon=False,loc='lower right',fontsize=13)
    ax=axes[1,0]
    primary=genes[genes.sensitivity=='all_animals']
    for accession,color in [('OSD-561','#2878a0'),('OSD-562','#bd6b2d')]:
        sub=primary[primary.accession==accession].dropna(subset=['OLS_interaction_log2','log2FoldChange'])
        r=stats.spearmanr(sub.OLS_interaction_log2,sub.log2FoldChange).statistic
        ax.scatter(sub.OLS_interaction_log2,sub.log2FoldChange,s=3,alpha=.13,color=color,
                   label=f'{accession}: Spearman r={r:.3f}')
    lim=np.nanquantile(abs(primary[['OLS_interaction_log2','log2FoldChange']].to_numpy()),.998)
    lim=max(1,lim);ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim);ax.plot([-lim,lim],[-lim,lim],color='#666666',ls='--',lw=1)
    ax.set_xlabel('Gene interaction: normalized-log OLS');ax.set_ylabel('Count-model interaction (log2)')
    ax.set_title('c  Gene interaction estimator agreement',loc='left',fontweight='bold')
    ax.text(.02,.98,'99.8th-percentile display; full tables retained',transform=ax.transAxes,va='top',fontsize=14)
    ax.legend(frameon=False,loc='lower right');ax.spines[['top','right']].set_visible(False)
    ax=axes[1,1];ax.axis('off')
    cols=[('GC','LAR','12'),('GC','LAR','29'),('FLT','LAR','12'),('FLT','LAR','29'),
          ('GC','ISS_T','12'),('GC','ISS_T','29'),('FLT','ISS_T','12'),('FLT','ISS_T','29')]
    rows=[];rowlabels=[]
    for accession in ['OSD-561','OSD-562']:
        for sensitivity in ['all_animals','exclude_any_region_metadata_conflict']:
            s=support[(support.accession==accession)&(support.sensitivity==sensitivity)]
            rows.append([int(s[(s.flight==f)&(s.cohort==c)&(s.age==a)].n_animals.sum()) for f,c,a in cols])
            rowlabels.append(accession+('\nFull' if sensitivity=='all_animals' else '\nExclusions'))
    clabels=[f'{f}\n{c.replace("ISS_T","ISS")}\n{a}' for f,c,a in cols]
    grid=ax.inset_axes([.04,.48,.94,.34])
    from matplotlib.colors import ListedColormap
    grid.imshow(np.array(rows)==0,cmap=ListedColormap(['#f4f4f4','#f6c8be']),vmin=0,vmax=1,aspect='auto')
    grid.set_xticks(range(8),clabels,fontsize=13);grid.xaxis.tick_top()
    grid.set_yticks(range(4),['Cb full','Cb excl','Hip full','Hip excl'],fontsize=14)
    grid.tick_params(length=0)
    for i,row in enumerate(rows):
        for j,value in enumerate(row):grid.text(j,i,str(value),ha='center',va='center',fontsize=16)
    ax.set_title('d  Age support in each design cell',loc='left',fontweight='bold')
    ax.text(.02,.36,'Ages in weeks; ISS denotes ISS-T.\nInteraction = (FLT−GC in ISS-T)\n                  − (FLT−GC in LAR).\nAdditive age; all main effects included.\nWhole conflict animals excluded.\nEmpty cells require age extrapolation.\nWild inference is approximate.\nCollection cohorts are observational.',transform=ax.transAxes,va='top',linespacing=1.35,fontsize=14)
    for ext in ['png','svg']:
        fig.savefig(OUT/f'figure_flight_collection_interaction.{ext}',dpi=600,bbox_inches='tight')
    plt.close(fig)


def main(args):
    check_input_integrity()
    OUT.mkdir(parents=True,exist_ok=True)
    annotation=pd.read_csv(BASE/'inputs/gene_annotation.csv.gz')
    datasets={n:load_dataset(n) for n in ['561','562']}
    conflicts=set(pd.concat([d['meta'] for d in datasets.values()]).query('location_conflict').animal)
    panel_rows=[];gene_frames=[];support_rows=[];manifest=[];diags=[];checks=[];omissions=[]
    for number,d in datasets.items():
        meta=d['meta'].copy();meta['any_region_location_conflict']=meta.animal.isin(conflicts)
        manifest.append(meta)
        panel_index=indices(d['y'].index,annotation)
        if any(len(i)<2 for i in panel_index.values()):raise ValueError('Insufficient panel coverage')
        for sensitivity in ['all_animals','exclude_any_region_metadata_conflict']:
            m=meta if sensitivity=='all_animals' else meta[~meta.any_region_location_conflict]
            x=design(m).to_numpy();y=d['y'].loc[:,m.index].to_numpy()
            # Independent direct least squares and a four-level group parameterization.
            coef=np.linalg.lstsq(x,y.T,rcond=None)[0].T
            cell=np.column_stack([(m.flight.eq(f)&m.cohort.eq(c)).to_numpy(float)
                                  for f,c in [('GC','LAR'),('FLT','LAR'),('GC','ISS_T'),('FLT','ISS_T')]])
            alternate=np.column_stack([(m.age=='29').to_numpy(float),cell])
            alt=np.linalg.lstsq(alternate,y.T,rcond=None)[0].T
            independent=alt[:,4]-alt[:,3]-alt[:,2]+alt[:,1]
            assert np.allclose(coef[:,-1],independent,atol=1e-10)
            checks.append({'accession':f'OSD-{number}','sensitivity':sensitivity,'rank':int(np.linalg.matrix_rank(x)),
                           'independent_four_cell_contrast_max_error':float(abs(coef[:,-1]-independent).max())})
            for (f,c,a),s in m.groupby(['flight','cohort','age'],sort=True):
                support_rows.append({'accession':f'OSD-{number}','sensitivity':sensitivity,'flight':f,'cohort':c,'age':a,'n_animals':len(s)})
            yp=np.vstack([y[idx].mean(0) for idx in panel_index.values()])
            beta,se,p,ci,pc,nb,mc_low,mc_high=panel_inference(yp,m,args.bootstraps,args.wild_draws,SEED+int(number))
            for i,(panel,idx) in enumerate(panel_index.items()):
                panel_rows.append({'accession':f'OSD-{number}','region':m.tissue.iloc[0],'sensitivity':sensitivity,'panel':panel,
                                   'n_animals':len(m),'n_genes':len(idx),'interaction_log2':beta[i],'HC3_se':se[i],
                                   'bootstrap_lower95':ci[0,i],'bootstrap_upper95':ci[1,i],'p_wild_approximate':p[i],
                                   'p_Monte_Carlo_lower95':mc_low[i],'p_Monte_Carlo_upper95':mc_high[i],
                                   'FLT_minus_GC_LAR_log2':pc[i,2],'FLT_minus_GC_ISST_log2':pc[i,2]+pc[i,4],
                                   'n_bootstraps':nb,'n_wild_draws':args.wild_draws,'definition':DEFINITION,
                                   'age_support_warning':('No older FLT ISS-T animals; common additive-age extrapolation'
                                                         if not ((m.flight=='FLT')&(m.cohort=='ISS_T')&(m.age=='29')).any() else
                                                         ('Only one older FLT ISS-T animal' if ((m.flight=='FLT')&(m.cohort=='ISS_T')&(m.age=='29')).sum()==1 else ''))})
            for animal in m.animal:
                keep=m.animal.ne(animal).to_numpy();xb=x[keep]
                if np.linalg.matrix_rank(xb)<x.shape[1]:continue
                b=np.linalg.lstsq(xb,yp[:,keep].T,rcond=None)[0][-1]
                omissions.extend([{'accession':f'OSD-{number}','sensitivity':sensitivity,'panel':panel,'omitted_animal':animal,
                                   'interaction_without_animal':b[i],'full_interaction':beta[i]} for i,panel in enumerate(panel_index)])
            g,diag=counts_fit(d,m,sensitivity,args.force);diags.append(diag)
            g['OLS_interaction_log2']=coef[:,-1]
            if sensitivity=='all_animals':checks.append(cached_contrast_check(g,f'OSD-{number}'))
            gene_frames.append(g)
    panels=pd.DataFrame(panel_rows)
    assert len(panels)==32 and panels.groupby(['accession','sensitivity']).size().eq(8).all()
    panels['q_BH_complete_panel_family']=bh(panels.p_wild_approximate);panels['panel_BH_family_size']=len(panels)
    genes=pd.concat(gene_frames,ignore_index=True).merge(annotation[['ensembl_id','gene_symbol']],on='ensembl_id',how='left',validate='many_to_one')
    genes['q_BH_complete_gene_interaction_family']=bh(genes.pvalue)
    genes['gene_BH_family_size']=int(genes.pvalue.notna().sum());genes['n_planned_gene_effects']=len(genes)
    for i,row in panels.iterrows():
        sub=genes[(genes.accession==row.accession)&(genes.sensitivity==row.sensitivity)]
        idx=indices(sub.ensembl_id,annotation)[row.panel]
        panels.loc[i,'mean_NB_interaction_log2']=sub.iloc[idx].log2FoldChange.mean()
    panels.to_csv(OUT/'panel_flight_collection_interactions.csv',index=False)
    genes.to_csv(OUT/'gene_flight_collection_interactions_exploratory.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    pd.DataFrame(omissions).to_csv(OUT/'leave_one_animal_out_panel_interactions.csv',index=False)
    support=pd.DataFrame(support_rows);support.to_csv(OUT/'age_support_by_design_cell.csv',index=False)
    pd.concat(manifest).to_csv(OUT/'animal_manifest.csv',index=False)
    summary={'definition':DEFINITION,'design_columns':COLS,'seed':SEED,'n_panel_tests':len(panels),
             'n_panel_BH05':int((panels.q_BH_complete_panel_family<=.05).sum()),
             'n_planned_gene_effects':len(genes),'n_evaluable_gene_tests':int(genes.pvalue.notna().sum()),
             'n_gene_BH05':int((genes.q_BH_complete_gene_interaction_family<=.05).sum()),
             'conflicting_animals_either_region':sorted(conflicts),'count_fit_diagnostics':diags,
             'independent_model_checks':checks,'python':platform.python_version(),
             'inference_assumptions':'Approximate fixed-design HC3-studentized restricted-model Rademacher wild bootstrap; independent animals, finite conditional variance, correctly specified common additive age mean model. Percentile animal bootstraps stratified by flight x cohort x age. NB Wald inference is exploratory with empirical-Bayes dispersion estimation and Cook filtering.',
             'limitations':'Collection cohort is observational and bundles collection history and other factors. Interaction tests do not identify microgravity specifically or establish randomized causal flight effects. Metadata exclusions reduce age overlap, especially cerebellum. Normalized bulk counts do not identify cell-type-specific regulation.'}
    (OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    figure(panels,genes,support)
    primary=panels[panels.sensitivity=='all_animals']
    text='# Within-region flight × collection-cohort interaction\n\n'
    text+=f'Each region used all collapsed source animals (26 cerebellum; 27 hippocampus) and the full model intercept + age29 + flight + ISS-T + flight×ISS-T. The direct interaction is {DEFINITION}. Panel uncertainties use {args.bootstraps:,} stratified whole-animal bootstraps and {args.wild_draws:,} HC3-studentized restricted-model wild draws. All 32 region × panel × sensitivity tests share one BH family.\n\n'
    text+=f'{summary["n_panel_BH05"]} panel tests passed BH q ≤ 0.05 (minimum q={panels.q_BH_complete_panel_family.min():.4g}). '
    for row in primary.sort_values('p_wild_approximate').head(4).itertuples():
        text+=f'{row.region} {row.panel}: interaction {row.interaction_log2:+.3f} log2, 95% animal-bootstrap interval [{row.bootstrap_lower95:+.3f}, {row.bootstrap_upper95:+.3f}], approximate p={row.p_wild_approximate:.4g}, q={row.q_BH_complete_panel_family:.4g}. '
    text+='\n\nActual integer counts were refitted in PyDESeq2 with the same explicit five-column design and numeric interaction contrast [0,0,0,0,1]. '
    text+=f'The complete exploratory gene family contains {summary["n_evaluable_gene_tests"]:,} evaluable interaction tests across both regions and sensitivities; {summary["n_gene_BH05"]} passed BH q ≤ 0.05. Cook-filtered unavailable p-values remain unavailable. No favorable post hoc gene family was selected.\n\n'
    text+='Metadata sensitivity excludes any animal with a location disagreement in either regional deposit. This removes all older FLT ISS-T animals from cerebellum and leaves one in hippocampus, reducing the support for an age-adjusted interaction. The coefficient assumes a common additive age effect and cannot resolve age-dependent interactions in the empty design cell. Differences between strata are estimated directly; significance in one stratum and nonsignificance in another are not evidence for interaction. Biological interpretation remains observational, and an interaction cannot isolate microgravity from collection history or other cohort-linked factors.\n'
    (OUT/'manuscript_ready_summary.md').write_text(text)
    print(panels[['accession','sensitivity','panel','interaction_log2','p_wild_approximate','q_BH_complete_panel_family']].to_string(index=False),flush=True)
    print(f'Completed {len(panels)} panel tests and {len(genes)} planned gene effects',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bootstraps',type=int,default=4000)
    parser.add_argument('--wild-draws',type=int,default=499999)
    parser.add_argument('--force',action='store_true')
    parser.add_argument('--figures-only',action='store_true')
    args=parser.parse_args()
    if args.figures_only:
        figure(pd.read_csv(OUT/'panel_flight_collection_interactions.csv'),
               pd.read_csv(OUT/'gene_flight_collection_interactions_exploratory.csv.gz'),
               pd.read_csv(OUT/'age_support_by_design_cell.csv'))
    else:
        main(args)
