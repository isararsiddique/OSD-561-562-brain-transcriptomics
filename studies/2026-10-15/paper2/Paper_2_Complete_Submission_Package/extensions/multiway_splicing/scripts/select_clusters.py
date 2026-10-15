"""Frozen all26 multiway selection, exact transcript-intron mapping; no labels read."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
from pathlib import Path
from collections import defaultdict, Counter
import gzip, re, json, hashlib, argparse
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[2]
INPUT = ROOT / 'research/paper2_beta_20261014/inputs/all_animal_junction_counts.csv.gz'
GTF = ROOT / 'research/osdr_brain_paper_2026/data/Mus_musculus.GRCm39.112.gtf.gz'

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--base',type=Path,default=INPUT.parent,help='directory containing frozen all_animal_junction_counts.csv.gz and animal_metadata.csv');parser.add_argument('--gtf',type=Path,default=GTF,help='only needed to build the bundled annotation lookup initially');args=parser.parse_args()
    counts_path=args.base/'all_animal_junction_counts.csv.gz'
    assert (BASE/'RETROSPECTIVE_PLAN_2026-10-15.md').exists()
    df=pd.read_csv(counts_path,index_col=0)
    assert df.shape==(805562,26) and df.index.is_unique and df.columns.is_unique
    x=df.to_numpy(dtype=np.int64)
    assert (x>=0).all()
    p=df.index.str.split(':',expand=True)
    info=pd.DataFrame(dict(chrom=p.get_level_values(0),strand=p.get_level_values(1),start=p.get_level_values(2),end=p.get_level_values(3)),index=df.index)
    qualifying=(x>=5).sum(axis=1)>=3
    selected=[]; candidate=[]; funnel=[]
    for boundary in ['start','end']:
        keys=(info.chrom+':'+info.strand+':'+info[boundary]).to_numpy()
        qids=np.flatnonzero(qualifying)
        qgroups=pd.Series(qids).groupby(keys[qids],sort=True).apply(list)
        multi=qgroups[qgroups.map(len)>=3]
        needed=set(multi.index)
        mask=np.array([k in needed for k in keys])
        totals=df.iloc[np.flatnonzero(mask)].groupby(keys[mask],sort=True).sum()
        for key, ids in multi.items():
            y=x[ids]; reads=y.sum(axis=0); allreads=totals.loc[key].to_numpy()
            dom=np.divide(reads,allreads,out=np.zeros(26),where=allreads>0)
            frac=np.divide(y,reads[None,:],out=np.zeros_like(y,dtype=float),where=reads[None,:]>0)
            means=frac.mean(axis=1)
            cov_ok=bool((reads>=20).all()); dom_ok=bool((dom>=.90-1e-12).all())
            comp_ok=bool((means>=.01-1e-12).all() and (means>=.05-1e-12).sum()>=3 and means.max()<=.90+1e-12)
            cid=f'{key}:{boundary}'
            row=dict(cluster_id=cid,boundary_type=boundary,chrom=key.split(':')[0],strand=int(key.split(':')[1]),boundary=int(key.split(':')[2]),n_alternatives=len(ids),min_selected_reads=int(reads.min()),min_dominance=float(dom.min()),min_pooled_proportion=float(means.min()),max_pooled_proportion=float(means.max()),n_pooled_alternatives_ge_05=int((means>=.05-1e-12).sum()),coverage_pass=cov_ok,dominance_pass=dom_ok,composition_pass=comp_ok,selected=cov_ok and dom_ok and comp_ok)
            candidate.append(row)
            if row['selected']: selected.append((row,np.array(ids),frac,reads,dom))
        ca=pd.DataFrame(candidate); ca=ca[ca.boundary_type==boundary]
        funnel.append(dict(boundary_type=boundary,qualifying_boundary_groups=len(qgroups),at_least_three=len(multi),coverage_pass=int(ca.coverage_pass.sum()),coverage_and_dominance_pass=int((ca.coverage_pass&ca.dominance_pass).sum()),selected=int(ca.selected.sum())))
        print(boundary, funnel[-1],flush=True)
    selected.sort(key=lambda z:z[0]['cluster_id'])
    assert selected
    definitions=pd.DataFrame([s[0] for s in selected])
    candidate=pd.DataFrame(candidate).sort_values('cluster_id')
    candidate.to_csv(BASE/'tables/all_candidate_clusters.csv.gz',index=False)
    pd.DataFrame(funnel).to_csv(BASE/'tables/selection_funnel.csv',index=False)
    needed={str(df.index[j]) for _,ids,*_ in selected for j in ids}
    intron_genes=defaultdict(set); symbols={}; exons=defaultdict(list); transcript_info={}
    lookup_path=BASE/'annotation/junction_gene_lookup.csv.gz';source_path=BASE/'annotation/source.json'
    if lookup_path.exists():
        lookup=pd.read_csv(lookup_path).fillna('');assert needed<=set(lookup.junction_id),'frozen lookup does not cover changed input'
        for row in lookup.itertuples():
            intron_genes[row.junction_id]=set(filter(None,row.matching_gene_ids.split(';')))
            for gene,symbol in zip(row.matching_gene_ids.split(';'),row.matching_gene_symbols.split(';')):
                if gene:symbols[gene]=symbol
        annotation_source=json.loads(source_path.read_text())
        assert digest(counts_path)==annotation_source['counts_sha256'],'annotation lookup is frozen to the original input'
    else:
        with gzip.open(args.gtf,'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                z=line.rstrip().split('\t')
                if len(z)<9 or z[2]!='exon':continue
                gene=re.search(r'gene_id "([^\"]+)"',z[8])[1]
                transcript=re.search(r'transcript_id "([^\"]+)"',z[8])[1]
                symbol=re.search(r'gene_name "([^\"]+)"',z[8])
                if symbol:symbols[gene]=symbol[1]
                transcript_info[transcript]=(z[0],1 if z[6]=='+' else 2,gene)
                exons[transcript].append((int(z[3]),int(z[4])))
        for t,coordinates in exons.items():
            chrom,strand,gene=transcript_info[t]; ordered=sorted(set(coordinates))
            for a,b in zip(ordered,ordered[1:]):
                key=f'{chrom}:{strand}:{a[1]+1}:{b[0]-1}'
                if key in needed:intron_genes[key].add(gene)
        lookup_path.parent.mkdir(exist_ok=True)
        pd.DataFrame([dict(junction_id=j,matching_gene_ids=';'.join(sorted(intron_genes[j])),matching_gene_symbols=';'.join(symbols.get(g,'') for g in sorted(intron_genes[j]))) for j in sorted(needed)]).to_csv(lookup_path,index=False)
        annotation_source=dict(source='Mus_musculus.GRCm39.112.gtf.gz',url='https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz',gtf_sha256=digest(args.gtf),counts_sha256=digest(counts_path),mapping='exact chromosome/strand transcript introns, upstream exon end+1 through downstream exon start-1',n_lookup_junctions=len(needed),scope='all retained alternatives; annotation for other junctions is not bundled')
        source_path.write_text(json.dumps(annotation_source,indent=2))
    alternatives=[]; coverage=[]; features=[]; offsets=[0]; overlap=defaultdict(list)
    for row,ids,frac,reads,dom in selected:
        jid=[str(df.index[i]) for i in ids]; common=set.intersection(*[intron_genes[j] for j in jid])
        gene=next(iter(common)) if len(common)==1 else ''
        row.update(gene_mapping_count=len(common),gene_id=gene,gene_symbol=symbols.get(gene,''),gene_candidates=';'.join(sorted(common)))
        features.append(np.sqrt(frac).T); offsets.append(offsets[-1]+len(ids))
        for a,(j,v) in enumerate(zip(jid,frac)):
            overlap[j].append(row['cluster_id'])
            for s,count,proportion in zip(df.columns,x[ids[a]],v):
                alternatives.append(dict(cluster_id=row['cluster_id'],alternative_index=a,junction_id=j,sample=s,count=int(count),proportion=float(proportion),matching_gene_ids=';'.join(sorted(intron_genes[j]))))
        entropy=-(np.where(frac>0,frac*np.log(np.maximum(frac,1e-300)),0)).sum(axis=0)
        for s,tot,d,h in zip(df.columns,reads,dom,entropy):coverage.append(dict(cluster_id=row['cluster_id'],sample=s,selected_reads=int(tot),dominance=float(d),entropy_nats=float(h),normalized_entropy=float(h/np.log(len(ids)))))
    definitions=pd.DataFrame([s[0] for s in selected])
    definitions.to_csv(BASE/'tables/selected_cluster_definitions.csv',index=False)
    pd.DataFrame(alternatives).to_csv(BASE/'tables/animal_alternative_compositions.csv.gz',index=False)
    pd.DataFrame(coverage).to_csv(BASE/'tables/animal_cluster_coverage.csv.gz',index=False)
    junction_overlap=pd.DataFrame([dict(junction_id=j,n_clusters=len(c),cluster_ids=';'.join(c)) for j,c in sorted(overlap.items())])
    junction_overlap.to_csv(BASE/'tables/junction_overlap.csv',index=False)
    pair_counter=Counter()
    import itertools
    for c in overlap.values():
        for a,b in itertools.combinations(sorted(c),2):pair_counter[(a,b)]+=1
    pd.DataFrame([dict(cluster_a=a,cluster_b=b,n_shared_junctions=n) for (a,b),n in pair_counter.items()],columns=['cluster_a','cluster_b','n_shared_junctions']).to_csv(BASE/'tables/overlapping_cluster_pairs.csv',index=False)
    np.savez_compressed(BASE/'tables/frozen_features.npz',H=np.concatenate(features,axis=1),offsets=np.array(offsets),samples=np.array(df.columns,dtype=str),cluster_ids=definitions.cluster_id.to_numpy(dtype=str))
    manifest=dict(plan_sha256=digest(BASE/'RETROSPECTIVE_PLAN_2026-10-15.md'),counts_sha256=digest(counts_path),gtf_sha256=annotation_source['gtf_sha256'],n_animals=26,n_canonical_junctions=len(df),n_qualifying_junctions=int(qualifying.sum()),n_clusters=len(selected),n_alternative_occurrences=int(definitions.n_alternatives.sum()),n_unique_retained_junctions=len(overlap),n_reused_junctions=int((junction_overlap.n_clusters>1).sum()),n_overlapping_cluster_pairs=len(pair_counter),n_unique_mapped_clusters=int((definitions.gene_mapping_count==1).sum()),n_unique_mapped_genes=int(definitions.loc[definitions.gene_mapping_count==1,'gene_id'].nunique()),selection_uses_group_labels=False)
    (BASE/'tables/selection_summary.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':main()
