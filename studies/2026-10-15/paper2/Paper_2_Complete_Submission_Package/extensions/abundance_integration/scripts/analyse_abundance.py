"""Relative gene abundance and local junction use, with complete test families."""
from pathlib import Path
import argparse, hashlib, importlib.util, json, shutil
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, t, spearmanr

HERE=Path(__file__).resolve().parents[1]
DEFAULT_BASE=HERE.parents[1]/'paper2_beta_20261014'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--base',type=Path,default=DEFAULT_BASE)
    p.add_argument('--gene',type=Path,default=HERE.parent/'gene_omnibus')
    args=p.parse_args(); base=args.base
    for d in ['inputs','results','figures']: (HERE/d).mkdir(exist_ok=True)
    rawpath=HERE/'inputs/STAR_gene_counts.csv.gz'
    if not rawpath.exists():
        source=base.parents[1]/'research/brain_integrative_2026/data/OSD-525/GLDS-525_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv'
        expected='9c729459fe3f89755f89afd9c4c7447d'
        assert hashlib.md5(source.read_bytes()).hexdigest()==expected
        raw=pd.read_csv(source,index_col=0)
        raw.to_csv(rawpath,compression='gzip')
        (HERE/'inputs/provenance.json').write_text(json.dumps(dict(source=str(source),publisher_md5=expected,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),frozen_sha256=hashlib.sha256(rawpath.read_bytes()).hexdigest()),indent=2))
    raw=pd.read_csv(rawpath,index_col=0)
    provenance=json.loads((HERE/'inputs/provenance.json').read_text())
    assert hashlib.sha256(rawpath.read_bytes()).hexdigest()==provenance['frozen_sha256']
    libraries=pd.read_csv(base/'inputs/junction_library_manifest.csv')
    meta=pd.read_csv(base/'inputs/animal_metadata.csv')
    assert set(raw.columns)==set(libraries.library)
    mapping=libraries.set_index('library')['sample']
    collapsed=raw.T.groupby(mapping).sum().T.loc[:,meta['sample']]
    assert len(raw.columns)==29 and collapsed.shape[1]==26 and np.all(collapsed.values>=0)
    positive=collapsed.loc[(collapsed>0).all(axis=1)]
    geom=np.exp(np.log(positive).mean(axis=1))
    sf=positive.div(geom,axis=0).median(axis=0)
    sf=sf/np.exp(np.log(sf).mean())
    normalized=collapsed.div(sf,axis=1)
    genes=pd.read_csv(args.gene/'results/gene_event_universe.csv.gz').drop_duplicates('gene_id').set_index('gene_id')
    eligible=(collapsed.reindex(genes.index).min(axis=1)>=20)&collapsed.reindex(genes.index).notna().all(axis=1)
    universe=genes[['gene_symbol']].copy();universe['minimum_raw_gene_reads']=collapsed.reindex(genes.index).min(axis=1);universe['included']=eligible
    universe.to_csv(HERE/'results/abundance_gene_universe.csv')
    norm=normalized.loc[genes.index[eligible]]
    x=np.log2(norm.to_numpy()+.5)
    pd.DataFrame(x,index=norm.index,columns=norm.columns).to_csv(HERE/'results/log2_normalized_gene_counts.csv.gz')
    qc=meta[['sample','Source Name','group']].copy();qc['raw_gene_reads']=collapsed.sum().values;qc['size_factor']=sf.values;qc['normalization_genes']=len(positive)
    qc.to_csv(HERE/'results/animal_abundance_qc.csv',index=False)
    spec=importlib.util.spec_from_file_location('junction_analyse',base/'scripts/analyse.py'); mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    out=[];glob=[]
    for a,b in mod.CONTRASTS:
        ai=np.flatnonzero(meta.group.eq(a));bi=np.flatnonzero(meta.group.eq(b));ix=np.r_[ai,bi];z=x[:,ix];w=mod.binary_weights(len(ix),len(ai))
        delta,ps,se,df=mod.event_test(z,w,len(ai))
        item=pd.DataFrame(dict(gene_id=norm.index,gene_symbol=genes.loc[norm.index].gene_symbol.values,contrast=f'{a}_minus_{b}',mean_log2_difference=delta,p_exact_abundance=ps,welch_se=se,welch_df=df,CI95_low=delta-t.ppf(.975,df)*se,CI95_high=delta+t.ppf(.975,df)*se))
        for fold in [1.2,1.5]:
            margin=np.log2(fold)
            item[f'p_equivalence_{fold}']=np.maximum(t.sf((delta+margin)/se,df),t.cdf((delta-margin)/se,df))
        out.append(item)
        g=mod.global_test(z,w,len(ai));g.pop('RMS_contrast_pp');g['RMS_log2_contrast']=np.sqrt(g['mean_squared_contrast']);g['contrast']=f'{a}_minus_{b}';glob.append(g)
        print(f'{a}-{b} exact abundance done: {len(item)} genes',flush=True)
    result=pd.concat(out,ignore_index=True)
    junction=pd.read_csv(args.gene/'results/gene_omnibus_complete.csv.gz')
    result=result.merge(junction[['gene_id','contrast','p_exact','q_BY_joint_3G','rms_delta_percentage_points']],on=['gene_id','contrast'],validate='one_to_one').rename(columns={'p_exact':'p_exact_junction','q_BY_joint_3G':'q_BY_junction'})
    for method in ['bh','by']:
        result[f'q_{method.upper()}_abundance']=false_discovery_control(result.p_exact_abundance,method=method)
        for fold in [1.2,1.5]:
            result[f'q_{method.upper()}_equivalence_{fold}']=false_discovery_control(result[f'p_equivalence_{fold}'],method=method)
            result[f'p_conjunction_{fold}']=np.maximum(result.p_exact_junction,result[f'p_equivalence_{fold}'])
            result[f'q_{method.upper()}_conjunction_{fold}']=false_discovery_control(result[f'p_conjunction_{fold}'],method=method)
    result.to_csv(HERE/'results/abundance_junction_complete.csv.gz',index=False)
    globaldf=pd.DataFrame(glob);globaldf['p_Holm_3']=mod.holm(globaldf.exact_p);globaldf.to_csv(HERE/'results/global_abundance_contrasts.csv',index=False)
    summary=[]
    for c,d in result.groupby('contrast',sort=False):
        row=dict(contrast=c,genes=len(d),abundance_BH=int((d.q_BH_abundance<.05).sum()),abundance_BY=int((d.q_BY_abundance<.05).sum()),junction_BY=int((d.q_BY_junction<.05).sum()),junction_AND_abundance_BY=int(((d.q_BY_junction<.05)&(d.q_BY_abundance<.05)).sum()),descriptive_rho_junction_rms_absolute_abundance=float(spearmanr(d.rms_delta_percentage_points,abs(d.mean_log2_difference)).statistic))
        for fold in [1.2,1.5]:
            row[f'equivalence_BY_{fold}']=int((d[f'q_BY_equivalence_{fold}']<.05).sum())
            row[f'conjunction_BY_{fold}']=int((d[f'q_BY_conjunction_{fold}']<.05).sum())
        summary.append(row)
    report=dict(animals=26,libraries=29,normalization_genes=len(positive),normalization='Median of ratios over all-positive genes; geometric-mean size factor 1',junction_genes=len(genes),tested_genes=int(eligible.sum()),excluded_genes=int((~eligible).sum()),joint_tests=len(result),summary=summary,global_results=globaldf.to_dict('records'),retrospective=True)
    (HERE/'results/summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
