#!/usr/bin/env python3
"""Build a publication-scale four-panel figure and data-linked report."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import analyse_gene_omnibus as a

COLORS={'FLT':'#D26926','GC':'#167B82','VIV':'#76619D'}
LABELS={'FLT_minus_GC':'FLT − habitat GC','GC_minus_VIV':'Habitat GC − vivarium','FLT_minus_VIV':'FLT − vivarium'}
ORDER=list(LABELS)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=Path,default=a.DEFAULT_BASE)
    args=parser.parse_args()
    root=a.HERE
    selected,x,metadata,genes=a.load(args.base)
    d=pd.read_csv(root/'results/gene_omnibus_complete.csv.gz')
    summary=pd.read_csv(root/'results/discovery_summary.csv').set_index('contrast').loc[ORDER]
    check=json.loads((root/'results/verification.json').read_text())
    sensitivity_check=json.loads((root/'results/sensitivity_verification.json').read_text())
    manifest=json.loads((root/'results/manifest.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.9),gridspec_kw={'height_ratios':[1,1.18]})
    fig.subplots_adjust(left=.075,right=.965,bottom=.12,top=.895,wspace=.34,hspace=.50)
    fig.suptitle('Gene-level local junction-usage associations in RR23 cerebellum',
                 fontsize=15,fontweight='bold',y=.975)
    fig.text(.5,.938,'1,983 genes · 3,239 measured events · 26 animals · 5,949 joint tests',
             ha='center',fontsize=10,color='#444444')
    ax=axes[0,0]
    positions=np.arange(3)
    for offset,col,color,label in [(-.19,'joint_BH_q_lt_0_05','#6E9BA5','Joint BH'),
                                  (.19,'joint_BY_q_lt_0_05','#C06F36','Joint BY')]:
        vals=summary[col].to_numpy()
        bars=ax.bar(positions+offset,vals,width=.36,color=color,label=label)
        for bar,val in zip(bars,vals):
            ax.text(bar.get_x()+bar.get_width()/2,max(val,0)+7,str(val),ha='center',fontsize=9)
    ax.set_xticks(positions,['FLT–GC','GC–VIV','FLT–VIV'])
    ax.set_ylim(0,340);ax.set_ylabel('Genes with joint q < 0.05')
    ax.set_title('A   Correction across every gene and contrast',loc='left',fontweight='bold',pad=12)
    ax.legend(frameon=False,fontsize=8,loc='upper right')
    ax.text(.5,-.23,'BH: dependence assumptions; BY: arbitrary dependence',
            ha='center',transform=ax.transAxes,fontsize=8,color='#555555')
    ax=axes[0,1]
    f=d[d.contrast=='FLT_minus_GC'].copy()
    passing=f.q_BY_joint_3G<.05
    ax.scatter(f.loc[~passing,'rms_delta_percentage_points'],
               -np.log10(f.loc[~passing,'q_BY_joint_3G']),s=9,c='#A4AFB2',alpha=.6,rasterized=True)
    ax.scatter(f.loc[passing,'rms_delta_percentage_points'],
               -np.log10(f.loc[passing,'q_BY_joint_3G']),s=21,c='#C06F36',edgecolors='white',linewidths=.25,zorder=3)
    ax.axhline(-np.log10(.05),color='#373737',ls='--',lw=.8)
    ax.text(39,-np.log10(.05)+.06,'Joint BY q = 0.05',ha='right',fontsize=8)
    offsets={'Hnrnpa2b1':(-85,25),'Stxbp1':(-10,48),'Mag':(20,25),'Huwe1':(12,43)}
    for symbol,(dx,dy) in offsets.items():
        row=f[f.gene_symbol==symbol].iloc[0]
        ax.annotate(symbol,(row.rms_delta_percentage_points,-np.log10(row.q_BY_joint_3G)),
                    xytext=(dx,dy),textcoords='offset points',fontsize=8,
                    arrowprops={'arrowstyle':'-','lw':.5,'color':'#777777'})
    ax.set_ylim(-.04,1.85);ax.set_xlim(0,max(40,f.rms_delta_percentage_points.max()*1.05))
    ax.set_xlabel('Gene RMS contrast (percentage points)');ax.set_ylabel('−log10(joint BY q)')
    ax.set_title('B   FLT–GC: magnitude and corrected evidence',loc='left',fontweight='bold',pad=12)
    ax.text(.5,-.23,'Labels selected after analysis for illustration',ha='center',
            transform=ax.transAxes,fontsize=8,color='#555555')
    ax=axes[1,0]
    hn=f[f.gene_symbol=='Hnrnpa2b1'].iloc[0]
    rows=genes[hn.gene_id]
    raw=x[rows].T
    im=ax.imshow(raw,aspect='auto',interpolation='nearest',cmap='viridis',vmin=0,vmax=1)
    names=metadata['sample'].str.extract(r'_(F\d+|G\d+|V\d+)')[0]
    ax.set_yticks(np.arange(26),names,fontsize=7)
    ax.set_xticks(np.arange(len(rows)),[f'Choice {i+1}' for i in range(len(rows))],fontsize=8)
    for low,high,group in [(0,9,'FLT'),(9,18,'GC'),(18,26,'VIV')]:
        ax.plot([-.68,-.68],[low-.5,high-.5],lw=5,color=COLORS[group],clip_on=False)
    ax.axhline(8.5,color='white',lw=1.3);ax.axhline(17.5,color='white',lw=1.3)
    ax.set_ylabel('Animals; technical libraries summed first')
    ax.set_title('C   Hnrnpa2b1: all measured choices and animals',loc='left',fontweight='bold',pad=12)
    colorbar=fig.colorbar(im,ax=ax,pad=.025,fraction=.05)
    colorbar.set_label('Local junction-A fraction',fontsize=8)
    ax.text(.5,-.15,f'Post hoc illustration: FLT–GC RMS {hn.rms_delta_percentage_points:.2f} pp; '
            f'joint BY q = {hn.q_BY_joint_3G:.4f}',ha='center',transform=ax.transAxes,fontsize=8)
    ax=axes[1,1]
    st=f[f.gene_symbol=='Stxbp1'].iloc[0]
    rows2=genes[st.gene_id]
    source=[]
    for rr,subplot in [(rows,'C'),(rows2,'D')]:
        for j,row_index in enumerate(rr):
            row=selected.iloc[row_index].to_dict()
            row.update({'panel':subplot,'choice_label':f'Choice {j+1}','posthoc_illustrative':True})
            source.append(row)
    for j,row_index in enumerate(rows2):
        for k,group in enumerate(['FLT','GC','VIV']):
            mask=metadata.group.eq(group).to_numpy()
            vals=x[row_index,mask]*100
            center=j*4+k
            jitter=np.linspace(-.16,.16,len(vals))
            ax.scatter(center+jitter,vals,s=24,color=COLORS[group],alpha=.82,edgecolors='white',linewidths=.35)
            ax.plot([center-.25,center+.25],[vals.mean(),vals.mean()],color='#333333',lw=1.5)
    centers=np.arange(len(rows2))*4+1
    ax.set_xticks(centers,[f'Choice {j+1}' for j in range(len(rows2))])
    ax.set_ylim(0,100);ax.set_ylabel('Local junction-A fraction (%)')
    ax.set_title('D   Stxbp1: individual-animal measurements',loc='left',fontweight='bold',pad=12)
    handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=COLORS[g],
                    markeredgecolor='none',markersize=5,label=g) for g in COLORS]
    ax.legend(handles=handles,frameon=False,loc='lower left',fontsize=8,ncol=3)
    ax.text(.5,-.15,f'Post hoc illustration: FLT–GC RMS {st.rms_delta_percentage_points:.2f} pp; '
            f'joint BY q = {st.q_BY_joint_3G:.4f}',ha='center',transform=ax.transAxes,fontsize=8)
    fig.text(.5,.018,'Association only: flight and collection date are aliased. Local choices do not establish full isoforms or protein changes.',
             ha='center',fontsize=9,color='#444444')
    (root/'figures').mkdir(exist_ok=True)
    for extension in ['png','pdf','svg']:
        fig.savefig(root/f'figures/Figure_Gene_Omnibus.{extension}',dpi=320)
    plt.close(fig)
    pd.DataFrame(source).to_csv(root/'results/figure_source_events.csv',index=False)
    results_table=summary[['joint_BH_q_lt_0_05','joint_BY_q_lt_0_05','within_BH_q_lt_0_05','within_BY_q_lt_0_05']]
    lines=['| Contrast | Joint BH | Joint BY | Within-contrast BH | Within-contrast BY |',
           '|---|---:|---:|---:|---:|']
    for contrast,row in results_table.iterrows():
        lines.append('| '+LABELS[contrast]+' | '+' | '.join(str(int(v)) for v in row)+' |')
    significant=f[f.q_BY_joint_3G<.05].sort_values('gene_symbol')
    sensitivity=pd.read_csv(root/'results/dominance95_sensitivity_summary.csv')
    loo=json.loads((root/'results/descriptive_leave_one_out_summary.json').read_text())
    sf=sensitivity[sensitivity.contrast=='FLT_minus_GC'].iloc[0]
    caption=(f'Figure. Gene-level local junction-usage associations in RR23 cerebellum. '
        f'(A) Gene discoveries at q<0.05 after BH or BY correction across all {len(d):,} gene–contrast tests. '
        f'(B) FLT–GC RMS local-fraction contrasts against joint BY q values; names are illustrative post hoc selections. '
        f'(C) Every measured Hnrnpa2b1 local-choice fraction for all 26 animals, ordered by group. '
        f'(D) Every animal’s fraction at both measured Stxbp1 local choices; horizontal segments are group means. '
        f'FLT n=9, habitat GC n=9, VIV n=8. Technical libraries were summed before fraction calculation. '
        f'Choice coordinates and lexical junction orientation are recorded in figure_source_events.csv. '
        f'Panels C and D are descriptive illustrations of the gene-level tests, not additional independent discoveries. '
        f'Gene RMS values average squared contrast magnitudes over all measured local choices assigned uniquely to that gene. '
        f'FLT and control collection dates are perfectly aliased; the associations therefore cannot identify a flight-only mechanism.')
    (root/'figures/Figure_Gene_Omnibus_caption.txt').write_text(caption+'\n')
    report=f'''# Gene-level local junction-usage extension

This retrospective analysis adds a gene-level observable phenotype to the repaired RR23 all-animal panel. It identifies **47 FLT–habitat GC gene associations under joint arbitrary-dependence BY correction**, with all 5,949 gene–contrast tests retained. This supports a local RNA-processing association in this cohort. It does not separate flight exposure from collection date, demonstrate a full isoform change, or establish a protein, neuronal injury, or causal flight phenotype.

## Universe and measurement

The input panel was constructed from all 26 animals without group labels. The pre-existing primary measurement gate requires ≥20 unique pair-junction reads in every animal and ≥90% dominance of the two retained junctions among all retained junctions sharing their boundary. Both junctions must map to one common, unambiguous gene using exact GRCm39.112 intron annotation. This yields **3,239 local choices from 1,983 genes**, with 1–12 measured choices per gene (median 1). No event p value, event direction, group-specific contrast, or ontology is used to select genes. Genes excluded by annotation or measurement are outside the observable universe, not negative biological findings.

## Mathematics and inference

For gene g, let E_g be its measured choices and x_ei the unique-read junction-A fraction for choice e in animal i. The statistic for groups A and B is

`T_g = mean_{{e in E_g}}[(mean_A x_e − mean_B x_e)^2]`.

The reported effect magnitude is `100 × sqrt(T_g)` percentage points. Its magnitude can combine opposing local-choice directions; the complete table therefore also reports the signed mean contrast and an unbiased squared-distance estimate that subtracts within-group sampling variation. The latter may be negative and is not clamped.

Every exact allocation moves an animal's whole event vector. The test therefore preserves correlated choices and shared junctions within each gene, and does not require independent event coordinates. All 48,620 allocations of 18 animals into 9/9 labels, or 24,310 allocations of 17 animals into 9/8 labels, are enumerated. The two-sided, direction-free evidence is the proportion with squared contrast at least as large as observed, including ties. A fixed 1e−14 absolute tolerance in squared-fraction units handles floating-point ties. For the balanced 9/9 contrast the complementary allocation makes the minimum attainable p value 2/48,620.

The computation uses centered Gram kernels and exact allocation outer products; an independent direct group-mean calculation reproduces 12 selected gene p values across the three contrasts, including the most extreme, largest multi-event, middle, and largest-p genes. Manual BH/BY implementations independently reproduce the SciPy-adjusted results. All **{check['checks_passed']} verification checks passed**.

The primary correction family comprises every one of the **5,949 = 3 × 1,983** gene tests. BH has dependence assumptions; BY controls FDR under arbitrary dependence provided the exact p values are valid under the exchangeability null. Within-contrast corrections are secondary sensitivities and are not substituted for the primary joint correction.

## Results

{chr(10).join(lines)}

Every one of the 47 joint-BY FLT–GC genes attains the exact minimum p=2/48,620, with joint BY q={significant.q_BY_joint_3G.iloc[0]:.9f}. Their gene RMS contrasts range from {significant.rms_delta_percentage_points.min():.2f} to {significant.rms_delta_percentage_points.max():.2f} percentage points. The q value is close to the decision boundary and reflects the discrete exact-null resolution; an arbitrary stronger threshold would remove these calls. GC–VIV and FLT–VIV have no joint-BY gene discoveries, which does not establish the absence of control-context or flight-associated differences. In particular, global all-event tests and gene-wise tests ask different questions and have different multiplicity burdens.

The corrected genes are: {', '.join(significant.gene_symbol)}.

## Measurement and single-animal sensitivities

The stricter 95% all-animal boundary-dominance gate retains 2,710 choices from 1,724 genes. The original 1,983-gene, 5,949-test family is held fixed; the 259 genes losing all measured choices receive p=1. This sensitivity also yields 47 joint-BY FLT–GC calls, with **45/47 primary genes retaining significance**. Rps6kb1 retains one of its two choices but has q=0.070877; Ecpas loses its sole measured choice and receives p=1. Septin4 and Atosa are sensitivity-only calls. They do not replace or enlarge the primary discovery set. The vivarium contrasts again have no joint-BY gene calls. Counts are 251/23/30 joint-BH genes for FLT–GC, GC–VIV, and FLT–VIV, respectively. This is a fixed-universe measurement sensitivity, not a search for a more favorable filter.

Descriptive leave-one-animal-out effects were calculated for every one of the 1,983 genes under each of the 18 FLT or GC omissions, retaining all 35,694 rows. For the 47 primary-BY genes, the least favorable omission retains at least **{100*loo['minimum_rms_ratio_any_primary_gene_omission']:.2f}%** of the full RMS contrast; the largest ratio is {100*loo['maximum_rms_ratio_any_primary_gene_omission']:.2f}%. The median per-gene minimum retained RMS is {100*loo['median_minimum_rms_ratio_primary_genes']:.2f}%. All 18 omissions preserve every measured local-choice direction in 40/47 genes, and preserve the signed mean direction in 42/47. Signed means can cancel opposing lexical local-choice directions and are not a biological direction for the whole gene. No leave-one-out p or q values are computed. Effect stability does not establish independent replication, remove collection-date aliasing, or show that the q values would remain below 0.05 after a smaller-sample test.

Hnrnpa2b1 and Hnrnpc provide RNA-binding follow-up candidates; Thoc2, Fubp1 and Csde1 provide additional RNA-handling candidates. Stxbp1 and Dnm1 provide synaptic trafficking follow-up candidates, and Mag provides a myelin-associated candidate. These are post hoc functional interpretations of gene identities, not tested pathway enrichments, cell-type assignments, or evidence that the associated local choices alter those functions. The biological interpretation requires independent transcript-level confirmation and orthogonal experiments.

## Design limits and next experiments

FLT samples were collected on 14 January 2021, and both habitat GC and vivarium controls on 17 January 2021. Thus exposure and collection date cannot be separately estimated from these animals. Exact label exchangeability is a null assumption, not a correction for this alias. Gene-level significance measures group association under that assumption. Differences between groups could reflect exposure, return handling, sample processing, date, or their combination.

The local choice fractions are read-derived proportions at shared splice boundaries; they do not identify a complete transcript, protein expression, a cell-type-specific phenotype, or neuronal injury. No independent mission replication has been run for this new gene test family. The historical event-level replication uses a different frozen panel and should not be presented as validation of these gene tests.

Useful follow-up would measure the corrected local choices with targeted RT-PCR or long-read RNA in a new cohort whose flight and habitat samples are processed together or randomized across processing batches. The RNA-binding and synaptic trafficking candidates can prioritize assay design without converting annotation into a mechanism.

## Figure

![Gene omnibus figure](figures/Figure_Gene_Omnibus.png)

{caption}

## Reproduction

Run `scripts/analyse_gene_omnibus.py`, `scripts/verify_gene_omnibus.py`, `scripts/sensitivity_dominance95.py`, `scripts/descriptive_leave_one_out.py`, `scripts/verify_sensitivities.py`, and `scripts/build_figure_report.py` with the established project Python environment, in that order. Every script accepts `--base`; the analysis script also accepts `--output`. The default reads only the frozen `paper2_beta_20261014` inputs. `results/manifest.json` records every input SHA-256 hash and the exact allocation family sizes. Complete unsuccessful tests remain in `gene_omnibus_complete.csv.gz`. This is an explicitly retrospective extension constructed after the original all-event results were available, not a preregistered confirmatory analysis.
'''
    (root/'REPORT.md').write_text(report)
    artifacts=['results/gene_omnibus_complete.csv.gz','results/arbitrary_dependence_corrected_genes.csv',
               'results/discovery_summary.csv','results/verification.json','results/figure_source_events.csv',
               'figures/Figure_Gene_Omnibus.png','figures/Figure_Gene_Omnibus.pdf',
               'figures/Figure_Gene_Omnibus.svg','figures/Figure_Gene_Omnibus_caption.txt','REPORT.md',
               'results/dominance95_gene_omnibus_complete.csv.gz','results/dominance95_sensitivity_summary.csv',
               'results/primary_BY_gene_dominance95_followup.csv',
               'results/descriptive_all_gene_leave_one_out.csv.gz','results/descriptive_gene_leave_one_out_summary.csv',
                'results/descriptive_leave_one_out_summary.json','results/sensitivity_verification.json','HISTORY.md']
    manifest['artifact_sha256']={p:a.sha(root/p) for p in artifacts}
    manifest['extension_script_sha256']={p.name:a.sha(p) for p in sorted((root/'scripts').glob('*.py'))}
    manifest['checks_passed_primary_and_sensitivity']=check['checks_passed']+sensitivity_check['checks_passed']
    (root/'results/manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Figure, caption, report, and artifact hashes written.')


if __name__=='__main__':
    main()
