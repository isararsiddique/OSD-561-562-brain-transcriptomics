from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Assembly-matched promoter RNA/ATAC context; exact donor-level inference."""
from pathlib import Path
import gzip,re,json,itertools
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control,spearmanr
B=stage_root('extended');V1=stage_root('core')
map=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('')
rna=pd.read_csv(B/'prepared/OSD352_RNA_pseudobulk.csv.gz',index_col=0)
ann=pd.read_csv(B/'prepared/OSD352_RNA_annotation.csv').set_index('ENSEMBL')
rna.index=ann.loc[rna.index,'SYMBOL'];rna=rna.groupby(level=0).sum();rcpm=rna/rna.sum()*1e6
atac=pd.read_csv(B/'prepared/OSD352_ATAC_pseudobulk.csv.gz',index_col=0)
meta=pd.read_csv(B/'tables/OSD352_l2_registry.csv').set_index('column')
coords={}
with gzip.open(B/'sources/Mus_musculus.GRCm38.102.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'):continue
        a=line.rstrip().split('\t')
        if a[2]!='gene':continue
        d=dict(re.findall(r'(\w+) "([^"]+)"',a[8]));s=d.get('gene_name')
        if s not in set(map.mouse_symbol_primary):continue
        tss=int(a[3]) if a[6]=='+' else int(a[4])
        # GTF 1-based inclusive -> half-open BED coordinate.
        coords.setdefault(s,[]).append(('chr'+a[0],tss-1-2000,tss+2000,d['gene_id'],a[6]))
peaks=atac.index.to_series().str.extract(r'^(.*)[:\-](\d+)\-(\d+)$');assert peaks.notna().all().all()
chrom=peaks[0].to_numpy();starts=peaks[1].astype(int).to_numpy();ends=peaks[2].astype(int).to_numpy()
prom=[];links=[]
for r in map.itertuples():
    mask=np.zeros(len(atac),bool)
    for ch,lo,hi,gid,strand in coords.get(r.mouse_symbol_primary,[]):
        mask|=(chrom==ch)&(ends>lo)&(starts<hi)
    for p in atac.index[mask]:links.append({'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'peak':p,'link_rule':'Overlap Ensembl102 GRCm38 gene TSS +/- 2kb; no enhancer link inferred'})
    if mask.any():prom.append(pd.Series(atac.loc[mask].sum(),name=r.mouse_symbol_primary))
pb=pd.DataFrame(prom);pb.to_csv(B/'prepared/OSD352_panel_promoter_counts.csv')
pd.DataFrame(links).to_csv(B/'tables/panel_promoter_peak_links.csv',index=False)
acpm=pb/atac.sum()*1e6
rows=[];donorrows=[]
for ct in ['CNS neurons','CNS glia','Immune cells']:
    md=meta[meta.cell_type.eq(ct)].sort_values('bio_origin');assert len(md)==5 and md.n_nuclei.ge(20).all()
    f=md.flight_status.eq('flight').to_numpy()
    for r in map.itertuples():
        row={'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'cell_type':ct,'promoter_peak_count':sum(x['mouse_symbol']==r.mouse_symbol_primary for x in links),'ATAC_log2CPM_effect':np.nan,'RNA_log2CPM_effect':np.nan,'ATAC_exact_p':1.,'RNA_ATAC_spearman':np.nan,'status':'missing_mapping_or_promoter'}
        if r.mouse_symbol_primary in acpm.index and r.mouse_symbol_primary in rcpm.index:
            av=np.log2(acpm.loc[r.mouse_symbol_primary,md.index]+1).to_numpy();rv=np.log2(rcpm.loc[r.mouse_symbol_primary,md.index]+1).to_numpy()
            effect=av[f].mean()-av[~f].mean();null=[]
            for ids in itertools.combinations(range(5),3):
                m=np.isin(np.arange(5),ids);null.append(av[m].mean()-av[~m].mean())
            row.update(ATAC_log2CPM_effect=effect,RNA_log2CPM_effect=rv[f].mean()-rv[~f].mean(),ATAC_exact_p=np.mean(np.abs(null)>=abs(effect)-1e-12),RNA_ATAC_spearman=spearmanr(av,rv).statistic if np.std(av)>0 and np.std(rv)>0 else np.nan,status='five_donor_descriptive_context')
            for i,(c,m) in enumerate(md.iterrows()):donorrows.append({'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'cell_type':ct,'donor':m.bio_origin,'flight':m.flight_status,'RNA_log2CPM_plus1':rv[i],'promoter_ATAC_log2CPM_plus1':av[i]})
        rows.append(row)
out=pd.DataFrame(rows);out['q_panel_three_celltypes']=false_discovery_control(out.ATAC_exact_p)
out.to_csv(B/'tables/RNA_ATAC_promoter_context.csv',index=False)
pd.DataFrame(donorrows).to_csv(B/'tables/RNA_ATAC_donor_values.csv',index=False)
(B/'tables/promoter_context_summary.json').write_text(json.dumps({'assembly':'GRCm38/mm10','annotation':'Ensembl release102','window':'Gene TSS +/- 2000 bp, every overlapping peak summed once per gene','mapped_panel_genes_with_peaks':len(pb),'planned_family':len(out),'exact_labels':10,'q05':int(out.q_panel_three_celltypes.lt(.05).sum()),'scope':'Same five RR3 donor groups for RNA and ATAC. Accessibility is not inflammasome activation. Overlap establishes promoter proximity, not gene regulation or causal links.'},indent=2))
print(out[out.mouse_symbol.isin(['Aif1','Pycard','Fos','Junb'])].to_string(index=False),flush=True)
