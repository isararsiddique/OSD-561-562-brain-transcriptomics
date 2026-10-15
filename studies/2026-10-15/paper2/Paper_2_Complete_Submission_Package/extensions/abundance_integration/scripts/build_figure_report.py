from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parents[1]
COLORS={'FLT':'#C55841','GC':'#277C8E','VIV':'#8B7EA7'}
def main():
    summary=json.loads((HERE/'results/summary.json').read_text());d=pd.read_csv(HERE/'results/abundance_junction_complete.csv.gz');qc=pd.read_csv(HERE/'results/animal_abundance_qc.csv');x=pd.read_csv(HERE/'results/log2_normalized_gene_counts.csv.gz',index_col=0)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(2,2,figsize=(10.6,8),layout='constrained')
    a=ax[0,0]
    for g in COLORS:
        q=qc[qc.group.eq(g)];a.scatter(q.raw_gene_reads/1e6,q.size_factor,c=COLORS[g],label=g,s=26)
    a.set(xlabel='Summed gene reads per animal (millions)',ylabel='Median-ratio size factor',title='A  Label-invariant abundance normalization');a.legend(frameon=False)
    a=ax[0,1];r=d[d.contrast.eq('FLT_minus_GC')];sig=r.q_BY_junction<.05
    a.scatter(r.mean_log2_difference,r.rms_delta_percentage_points,s=8,c='#a5b3bb',alpha=.5,rasterized=True)
    a.scatter(r.loc[sig,'mean_log2_difference'],r.loc[sig,'rms_delta_percentage_points'],s=22,c=COLORS['FLT'],label='Junction BY q < 0.05')
    a.axvline(0,c='#687984',lw=.7);a.set(xlabel='FLT − GC mean log2 relative abundance',ylabel='Junction-use RMS contrast (pp)',title='B  Different measurements of the same genes');a.legend(frameon=False,fontsize=8)
    a=ax[1,0];contrasts=['FLT_minus_GC','GC_minus_VIV','FLT_minus_VIV'];xx=np.arange(3)
    rows=[next(s for s in summary['summary'] if s['contrast']==c) for c in contrasts]
    for offset,key,label,color in [(-.2,'abundance_BH','Abundance BH','#86a5b0'),(0,'abundance_BY','Abundance BY','#274c62'),(.2,'junction_BY','Junction BY',COLORS['FLT'])]:
        vals=[r[key] for r in rows];a.bar(xx+offset,vals,width=.18,label=label,color=color)
        for k,v in enumerate(vals):a.text(xx[k]+offset,v+1,str(v),ha='center',fontsize=8)
    a.set(xticks=xx,xticklabels=['FLT–GC','GC–VIV','FLT–VIV'],ylabel='Genes at joint q < 0.05',title='C  Full-family abundance and junction tests');a.legend(frameon=False,fontsize=8);a.set_ylim(0,110)
    a=ax[1,1];names=['Hnrnpa2b1','Stxbp1','Opa1'];source=[]
    for j,name in enumerate(names):
        row=r[r.gene_symbol.eq(name)].iloc[0];gid=row.gene_id;vals=x.loc[gid].to_numpy();vals=vals-vals.mean()
        for k,g in enumerate(COLORS):
            v=vals[qc.group.eq(g)];positions=j+(k-1)*.22+np.linspace(-.04,.04,len(v))
            a.scatter(positions,v,s=14,color=COLORS[g],label=g if j==0 else None);a.plot([j+(k-1)*.22-.07,j+(k-1)*.22+.07],[v.mean()]*2,c=COLORS[g],lw=2)
        source.append(row.to_dict())
    a.axhline(0,lw=.7,c='#86959b');a.set(xticks=range(3),xticklabels=names,ylabel='Centered log2 relative abundance',title='D  Abundance of illustrative junction genes');a.legend(frameon=False,fontsize=8)
    for ext in ['png','pdf','svg']:fig.savefig(HERE/f'figures/Figure_Abundance_Integration.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    pd.DataFrame(source).to_csv(HERE/'results/figure_illustrative_genes.csv',index=False)
    caption='Figure. Relative gene abundance and junction use in the same RR23 animals. (A) Summed STAR gene reads and label-invariant median-ratio size factors after technical-library collapse. (B) All 1,977 jointly measurable genes: mean log2 abundance difference versus gene-level junction RMS contrast; 46 of the 47 primary junction BY genes meet the additional gene-abundance coverage gate. (C) Complete-family discovery counts; abundance corrections include all 5,931 gene–contrast tests, and junction BY calls retain their original 5,949-test correction. (D) All animals for three post hoc illustrative genes, centered separately by gene; bars are group means. FLT n=9, GC n=9, VIV n=8. Nonsignificant abundance does not establish abundance equivalence or mechanistic independence. No junction/equivalence conjunction passes joint BY at either declared fold-change margin.'
    (HERE/'figures/Figure_Abundance_Integration_caption.txt').write_text(caption+'\n')
    report=f'''# Gene abundance and local junction usage

The deposited STAR gene-count matrix matched publisher MD5. Twenty-nine libraries were summed into the same 26 source animals used for junction analysis. Median-ratio normalization used {summary['normalization_genes']:,} genes positive in every animal. Of the fixed 1,983 junction genes, 1,977 had at least 20 raw gene reads in every animal; six exclusions are retained in the universe table. This is a restricted observable universe, not genome-wide gene-expression discovery.

Exact two-sided abundance tests were corrected across all 5,931 tests. Joint BH detects 72 FLT–GC, 8 GC–VIV and 24 FLT–VIV genes, whereas joint BY detects none. Global mean-squared log-abundance contrasts also do not pass Holm: adjusted p=0.133032, 0.423200 and 0.423200. These outcomes neither establish unchanged abundance nor prove a difference between the strengths of two statistical tests.

The descriptive FLT–GC correlation between junction RMS and absolute log-abundance contrast is 0.0572 across 1,977 genes. No gene-independence p value is attached. Forty-six primary junction-BY genes meet the abundance gate; Rian does not. The complete six-gene exclusion table should be used to confirm this identity.

Welch TOST sensitivity analyses use operational equivalence intervals ±log2(1.2) and ±log2(1.5), corresponding to geometric-mean ratios 1/1.2–1.2 or 1/1.5–1.5. These margins were declared before this fit but after the junction outcomes were known. None of the junction-association AND abundance-equivalence conjunctions survives joint BY at either margin. Therefore this study does not establish an abundance-independent junction phenotype. Separate abundance equivalence tests do pass for some genes, but selecting their intersection with separately significant junction tests is not substituted for the corrected conjunction family.

The analysis is retrospective and approximate Welch equivalence depends on its distributional assumptions. Exact permutation tests preserve whole-animal vectors and require exchangeability under the null. Flight/date aliasing, normalization assumptions and bulk-tissue composition remain. Relative gene counts and local junction fractions are different measurements from the same libraries and are not independent modalities.

Scripts accept the beta input location via --base. Full measurements, all unsuccessful tests, QC, source provenance and independent arithmetic checks accompany this report. The four-panel figure is descriptive where specified. See RETROSPECTIVE_ANALYSIS_PLAN.md for the declared family and margins.
'''
    excluded=pd.read_csv(HERE/'results/abundance_gene_universe.csv');report=report.replace('Rian does not. The complete six-gene exclusion table should be used to confirm this identity.', 'The excluded genes are '+', '.join(excluded.loc[~excluded.included,'gene_symbol'].fillna(excluded.loc[~excluded.included,'gene_id']).astype(str))+'.')
    (HERE/'REPORT.md').write_text(report)
if __name__=='__main__':main()
