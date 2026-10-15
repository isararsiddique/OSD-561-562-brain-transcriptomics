"""Exact transcript structure audit of the full pooled observable universe."""
import argparse,json
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from context_common import *

def region(position,rows):
    for feature,label in [('CDS','CDS'),('five_prime_utr','5prime_UTR'),('three_prime_utr','3prime_UTR')]:
        if any(r['feature']==feature and r['start']<=position<=r['end'] for r in rows):return label
    return 'no_CDS_annotation' if not any(r['feature']=='CDS' for r in rows) else 'exonic_nonCDS_unclassified'

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=DEFAULT_BASE);args=p.parse_args()
    multi,alternatives,events,needed=universe(args.base)
    prov=json.loads((HERE/'inputs/provenance.json').read_text())
    gtf=HERE/'inputs/required_gene_context.GRCm39.112.gtf.gz';assert sha(gtf)==prov['frozen_gtf_sha256']
    membership=[];transcripts=[];tcfeatures=[];tc_introns=set()
    for gene,group in gene_groups(gtf):
        bytrans=defaultdict(list)
        for row in group:
            tid=row['attrs'].get('transcript_id')
            if tid:bytrans[tid].append(row)
            if gene=='ENSMUSG00000036667':
                tcfeatures.append({k:row[k] for k in ['chrom','feature','start','end','strand','phase']}|row['attrs'])
        for tid,rows in bytrans.items():
            attr=rows[0]['attrs'];cds=sorted([(r['start'],r['end'],r['phase']) for r in rows if r['feature']=='CDS'])
            transcripts.append(dict(gene_id=gene,gene_symbol=attr.get('gene_name'),transcript_id=tid,
                transcript_version=attr.get('transcript_version'),transcript_name=attr.get('transcript_name'),
                transcript_biotype=attr.get('transcript_biotype'),strand=rows[0]['strand'],
                n_exons=sum(r['feature']=='exon' for r in rows),CDS_bases=sum(b-a+1 for a,b,_ in cds),
                CDS_structure=json.dumps(cds),five_prime_UTR_bases=sum(r['end']-r['start']+1 for r in rows if r['feature']=='five_prime_utr'),
                protein_ids=';'.join(sorted({r['attrs']['protein_id'] for r in rows if 'protein_id' in r['attrs']})),
                CCDS_ids=';'.join(sorted({r['attrs']['ccds_id'] for r in rows if 'ccds_id' in r['attrs']}))))
        for jid,tid,rows,donor,acceptor,ranks in transcript_introns(group):
            if gene=='ENSMUSG00000036667':tc_introns.add(jid)
            if jid not in needed:continue
            attr=rows[0]['attrs'];minus=donor['strand']=='-'
            ds=donor['start'] if minus else donor['end'];ac=acceptor['end'] if minus else acceptor['start']
            membership.append(dict(junction_id=jid,gene_id=gene,gene_symbol=attr.get('gene_name'),
                transcript_id=tid,transcript_name=attr.get('transcript_name'),transcript_biotype=attr.get('transcript_biotype'),
                strand=donor['strand'],donor_exonic_boundary=ds,acceptor_exonic_boundary=ac,
                donor_exon_start=donor['start'],donor_exon_end=donor['end'],acceptor_exon_start=acceptor['start'],
                acceptor_exon_end=acceptor['end'],donor_exon_rank=ranks[(donor['start'],donor['end'])],
                acceptor_exon_rank=ranks[(acceptor['start'],acceptor['end'])],
                donor_boundary_region=region(ds,rows),acceptor_boundary_region=region(ac,rows)))
    m=pd.DataFrame(membership).drop_duplicates();t=pd.DataFrame(transcripts).drop_duplicates('transcript_id')
    m.to_csv(HERE/'tables/junction_transcript_memberships.csv.gz',index=False)
    t.to_csv(HERE/'tables/transcript_context_universe.csv.gz',index=False)
    pd.DataFrame(tcfeatures).to_csv(HERE/'tables/Tcaf1_GTF_features.csv',index=False)
    t[t.gene_symbol.eq('Tcaf1')].to_csv(HERE/'tables/Tcaf1_transcript_CDS_UTR_summary.csv',index=False)
    gm=m.groupby('junction_id').gene_id.apply(set).to_dict()
    tm=m.groupby('junction_id').transcript_id.apply(set).to_dict()
    lookup=alternatives[['cluster_id','junction_id']].drop_duplicates().groupby('cluster_id').junction_id.apply(list)
    corrected=pd.read_csv(args.base/'extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz')
    hits=set(corrected[corrected.BY_q_full_3E<.05].cluster_id)
    audit=[]
    for row in multi.itertuples():
        jids=lookup[row.cluster_id];common=set.intersection(*[gm.get(j,set()) for j in jids])
        ranks=m[m.junction_id.isin(jids)]
        audit.append(dict(cluster_id=row.cluster_id,n_alternatives=len(jids),BY_corrected_any_contrast=row.cluster_id in hits,
            junctions_all_exactly_match_transcript_introns=all(j in tm for j in jids),
            shared_gene_intersection=';'.join(sorted(common)),shared_gene_count=len(common),
            matches_original_gene_intersection=len(common)==row.gene_mapping_count and (not common or ';'.join(sorted(common))==row.gene_candidates),
            alternatives_with_multiple_gene_matches=sum(len(gm.get(j,set()))>1 for j in jids),
            compatible_transcripts=';'.join(sorted(set.union(*[tm.get(j,set()) for j in jids]))),
            compatible_transcript_biotypes=';'.join(sorted(set(ranks.transcript_biotype))),
            first_exon_donor_junctions=len(set(ranks.loc[ranks.donor_exon_rank.eq(1),'junction_id'])),
            junctions_with_any_5prime_UTR_boundary=len(set(ranks.loc[ranks.donor_boundary_region.eq('5prime_UTR')|ranks.acceptor_boundary_region.eq('5prime_UTR'),'junction_id']))))
    ca=pd.DataFrame(audit);ca.to_csv(HERE/'tables/all923_cluster_annotation_audit.csv',index=False)
    ca[ca.BY_corrected_any_contrast].to_csv(HERE/'tables/all41_corrected_cluster_annotation_audit.csv',index=False)
    ea=[]
    for row in events.itertuples():
        common=gm.get(row.junction_a,set())&gm.get(row.junction_b,set())
        ea.append(dict(event_id=row.event_id,gene_id=row.gene_id,gene_symbol=row.gene_symbol,
            exact_annotation_gene_intersection=';'.join(sorted(common)),matches_original=len(common)==1 and row.gene_id in common,
            junction_a_transcript_count=len(tm.get(row.junction_a,set())),junction_b_transcript_count=len(tm.get(row.junction_b,set()))))
    e=pd.DataFrame(ea);e.to_csv(HERE/'tables/all3239_binary_gene_event_annotation_audit.csv.gz',index=False)
    gene_hits=pd.read_csv(args.base/'extensions/gene_omnibus/results/arbitrary_dependence_corrected_genes.csv')
    summary=e.groupby(['gene_id','gene_symbol'],dropna=False).agg(measured_events=('event_id','size'),all_event_annotations_match=('matches_original','all')).reset_index()
    summary['primary_BY_gene']=summary.gene_id.isin(gene_hits.gene_id)
    summary.to_csv(HERE/'tables/all1983_gene_annotation_audit.csv',index=False)
    summary[summary.gene_symbol.isna()].to_csv(HERE/'tables/unnamed_gene_identity_audit.csv',index=False)
    summary[summary.primary_BY_gene].to_csv(HERE/'tables/all47_corrected_gene_annotation_audit.csv',index=False)
    tcj=lookup[TC_CLUSTER];m[m.junction_id.isin(tcj)].to_csv(HERE/'tables/Tcaf1_measured_junction_transcripts.csv',index=False)
    tc=alternatives[alternatives.cluster_id.eq(TC_CLUSTER)].copy()
    meta=pd.read_csv(args.base/'inputs/animal_metadata.csv').set_index('sample')
    tc['group']=tc['sample'].map(meta.group)
    split=tc.junction_id.str.split(':',expand=True);tc['intron_length_bp']=split[3].astype(int)-split[2].astype(int)+1
    tc=tc.sort_values(['sample','intron_length_bp']);tc.to_csv(HERE/'tables/Tcaf1_all26_alternative_read_support.csv',index=False)
    animal=tc.groupby(['sample','group']).apply(lambda d:pd.Series({
        'selected_reads':int(d['count'].sum()),'short5974_reads':int(d.loc[d.intron_length_bp.eq(5974),'count'].sum()),
        'middle22570_reads':int(d.loc[d.intron_length_bp.eq(22570),'count'].sum()),
        'long23000_reads':int(d.loc[d.intron_length_bp.eq(23000),'count'].sum()),
        'long23000_fraction':float(d.loc[d.intron_length_bp.eq(23000),'count'].sum()/d['count'].sum()),
        'long_vs_short_conditional_fraction':float(d.loc[d.intron_length_bp.eq(23000),'count'].sum()/d.loc[d.intron_length_bp.isin([5974,23000]),'count'].sum())}),include_groups=False).reset_index()
    animal.to_csv(HERE/'tables/Tcaf1_all26_long_short_intron_support.csv',index=False)
    groups=animal.groupby('group').agg(n=('sample','size'),minimum_selected_reads=('selected_reads','min'),
        median_selected_reads=('selected_reads','median'),mean_long_fraction=('long23000_fraction','mean'),
        mean_long_vs_short_conditional=('long_vs_short_conditional_fraction','mean')).reset_index()
    groups.to_csv(HERE/'tables/Tcaf1_group_read_support_summary.csv',index=False)
    correlations=[]
    for group,d in [('all26',animal),*animal.groupby('group')]:
        correlations.append(dict(group=group,n=len(d),spearman_rho_fraction_vs_selected_reads=float(spearmanr(d.long23000_fraction,d.selected_reads).statistic),
            interpretation='Descriptive coverage relation; no gene/animal significance assigned and no independent replication'))
    pd.DataFrame(correlations).to_csv(HERE/'tables/Tcaf1_descriptive_fraction_coverage_relation.csv',index=False)
    # Every available annotated Tcaf1 intron, including inclusion partners, retained.
    counts=pd.read_csv(args.base/'inputs/all_animal_junction_counts.csv.gz',index_col=0)
    local=counts.loc[counts.index.isin(tc_introns)].copy();local['junction_id']=local.index
    local.melt(id_vars='junction_id',var_name='sample',value_name='unique_reads').to_csv(HERE/'tables/Tcaf1_all_annotated_intron_read_counts.csv.gz',index=False)
    companion='6:2:42686572:42686892'
    cp=local.loc[companion].drop('junction_id').astype(int) if companion in local.index else pd.Series(dtype=int)
    unique_multiway=set(alternatives.junction_id)
    result=dict(junctions_required_by_full_observable_universe=len(needed),exact_transcript_matched_junctions=m.junction_id.nunique(),
        unique_multiway_junctions=len(unique_multiway),multiway_junctions_without_exact_transcript_match=len(unique_multiway-set(tm)),
        clusters=len(ca),corrected_clusters=int(ca.BY_corrected_any_contrast.sum()),
        clusters_with_every_alternative_exactly_matched=int(ca.junctions_all_exactly_match_transcript_introns.sum()),
        corrected_clusters_with_every_alternative_exactly_matched=int((ca.junctions_all_exactly_match_transcript_introns&ca.BY_corrected_any_contrast).sum()),
        all_cluster_annotation_intersections_match=bool(ca.matches_original_gene_intersection.all()),
        mapped_clusters=int(ca.shared_gene_count.eq(1).sum()),
        corrected_clusters_mapped=int((ca.shared_gene_count.eq(1)&ca.BY_corrected_any_contrast).sum()),
        mapped_clusters_with_multigene_alternative=int((ca.shared_gene_count.eq(1)&ca.alternatives_with_multiple_gene_matches.gt(0)).sum()),
        all_binary_annotations_match=bool(e.matches_original.all()),binary_gene_events=len(e),genes=len(summary),corrected_genes=int(summary.primary_BY_gene.sum()),
        genes_with_blank_symbol=int(summary.gene_symbol.isna().sum()),
        Tcaf1_group_support=groups.to_dict('records'),Tcaf1_companion_junction=companion,
        Tcaf1_companion_available_animals=len(cp),Tcaf1_companion_min_reads=int(cp.min()) if len(cp) else None,
        Tcaf1_companion_max_reads=int(cp.max()) if len(cp) else None,
        no_long_read_or_protein_data=True,reference_names_and_release_match_publisher=True)
    (HERE/'tables/summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
