"""Independent arithmetic and inference checks for abundance integration."""
from pathlib import Path
from itertools import combinations
import argparse, json
import numpy as np
import pandas as pd
from scipy.stats import t
HERE=Path(__file__).resolve().parents[1]
def adjust(p,by=False):
    p=np.asarray(p);ix=np.argsort(p);factor=sum(1/np.arange(1,len(p)+1)) if by else 1
    q=np.minimum.accumulate((p[ix]*len(p)*factor/np.arange(1,len(p)+1))[::-1])[::-1];ans=np.empty(len(p));ans[ix]=np.minimum(q,1);return ans
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=HERE.parents[1]/'paper2_beta_20261014');args=ap.parse_args()
    checks=[]
    def check(name,yes):
        assert yes,name
        checks.append(name)
    meta=pd.read_csv(args.base/'inputs/animal_metadata.csv');lib=pd.read_csv(args.base/'inputs/junction_library_manifest.csv')
    raw=pd.read_csv(HERE/'inputs/STAR_gene_counts.csv.gz',index_col=0)
    collapsed=np.column_stack([raw.loc[:,lib.loc[lib['sample'].eq(s),'library']].sum(axis=1) for s in meta['sample']])
    positive=collapsed[(collapsed>0).all(axis=1)]
    sf=np.median(positive/np.exp(np.log(positive).mean(axis=1))[:,None],axis=0);sf/=np.exp(np.log(sf).mean())
    qc=pd.read_csv(HERE/'results/animal_abundance_qc.csv')
    check('Independent technical-library count collapse',np.array_equal(collapsed.sum(axis=0),qc.raw_gene_reads))
    check('Independent median-ratio normalization',np.allclose(sf,qc.size_factor,rtol=1e-12))
    x=pd.read_csv(HERE/'results/log2_normalized_gene_counts.csv.gz',index_col=0)
    idx=raw.index.get_indexer(x.index);expected=np.log2(collapsed[idx]/sf[None,:]+.5)
    check('All normalized animal-gene measurements',np.allclose(x,expected,atol=1e-12))
    result=pd.read_csv(HERE/'results/abundance_junction_complete.csv.gz');glob=pd.read_csv(HERE/'results/global_abundance_contrasts.csv')
    for c,d in result.groupby('contrast',sort=False):
        a,b=c.split('_minus_');ai=np.flatnonzero(meta.group.eq(a));bi=np.flatnonzero(meta.group.eq(b));z=x.to_numpy()[:,np.r_[ai,bi]];na=len(ai);n=z.shape[1]
        order=np.argsort(d.p_exact_abundance.to_numpy());selected=np.unique([order[0],order[len(d)//2],order[-1]])
        v=z[selected];obs=abs(v[:,:na].mean(axis=1)-v[:,na:].mean(axis=1));counts=np.zeros(len(selected),int);total=0
        # Direct enumeration of group means, independent of the fitting weights.
        for chosen in combinations(range(n),na):
            mask=np.zeros(n,bool);mask[list(chosen)]=True
            ds=abs(v[:,mask].mean(axis=1)-v[:,~mask].mean(axis=1))
            counts+=ds>=obs-1e-12;total+=1
        check(c+' independently enumerated exact p values',np.allclose(counts/total,d.iloc[selected].p_exact_abundance,atol=1e-14))
        delta=z[:,:na].mean(axis=1)-z[:,na:].mean(axis=1);va=z[:,:na].var(axis=1,ddof=1)/na;vb=z[:,na:].var(axis=1,ddof=1)/(n-na);se=np.sqrt(va+vb);df=(va+vb)**2/(va**2/(na-1)+vb**2/(n-na-1))
        check(c+' observed contrast and unbiased distance',np.isclose(np.mean(delta**2),glob.set_index('contrast').loc[c,'mean_squared_contrast']) and np.isclose(np.mean(delta**2-va-vb),glob.set_index('contrast').loc[c,'unbiased_squared_distance']))
        for fold in [1.2,1.5]:
            margin=np.log2(fold);pe=np.maximum(1-t.cdf((delta+margin)/se,df),1-t.cdf((margin-delta)/se,df))
            check(c+f' Welch TOST {fold}',np.allclose(pe,d[f'p_equivalence_{fold}'],atol=1e-12))
    for by in [False,True]:
        name='BY' if by else 'BH'
        check(name+' abundance correction',np.allclose(adjust(result.p_exact_abundance,by),result[f'q_{name}_abundance']))
        for fold in [1.2,1.5]:
            check(name+f' equivalence {fold}',np.allclose(adjust(result[f'p_equivalence_{fold}'],by),result[f'q_{name}_equivalence_{fold}']))
            check(name+f' conjunction {fold}',np.allclose(adjust(np.maximum(result.p_exact_junction,result[f'p_equivalence_{fold}']),by),result[f'q_{name}_conjunction_{fold}']))
    check('Complete family with no selected positive subset',len(result)==len(x)*3 and result[['gene_id','contrast']].duplicated().sum()==0)
    report=dict(status='PASS',checks_passed=len(checks),checks=checks)
    (HERE/'results/verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
