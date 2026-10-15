"""Independent GTF parser and coordinate/animal audit; no analysis imports."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
from collections import defaultdict
import argparse,gzip,re,json,hashlib,zipfile
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1];PACKAGE=BASE.parents[1];TC='ENSMUSG00000036667';CID='6:2:42663893:start'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=PACKAGE);args=p.parse_args();T=BASE/'tables';prov=json.loads((BASE/'inputs/provenance.json').read_text());checks=[]
    def check(name,value,detail=''):
        assert bool(value),name+' '+str(detail)
        checks.append(dict(check=name,passed=True,detail=str(detail)))
    gtf=BASE/'inputs/required_gene_context.GRCm39.112.gtf.gz';check('immutable frozen GTF hash',sha(gtf)==prov['frozen_gtf_sha256'])
    for f in prov['publisher_files']:
        path=BASE/'inputs'/f['publisher_filename'];check('publisher MD5 and SHA '+path.name,hashlib.md5(path.read_bytes()).hexdigest()==f['publisher_md5'] and sha(path)==f['sha256'])
    with zipfile.ZipFile(BASE/'inputs/processing_info_GLbulkRNAseq.zip') as z:text=z.read('nextflow_processing_info_GLbulkRNAseq.txt').decode()
    software=(BASE/'inputs/software_versions_GLbulkRNAseq.md').read_text()
    check('publisher reference and STAR names',all(v in text for v in ['GRCm39.112','--twopassMode Basic','ERCC92']) and '2.7.11b' in software,'reference names verified; FASTA sequences not reconstructed')
    alt=pd.read_csv(args.base/'extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz');events=pd.read_csv(args.base/'extensions/gene_omnibus/results/gene_event_universe.csv.gz');defs=pd.read_csv(args.base/'extensions/multiway_splicing/tables/selected_cluster_definitions.csv');needed=set(alt.junction_id)|set(events.junction_a)|set(events.junction_b)
    transcript=defaultdict(list);genes=set();tcfeatures=[];records=0
    with gzip.open(gtf,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            z=line.rstrip('\n').split('\t');attrs={}
            for bit in z[8].split(';'):
                m=re.match(r'\s*([^ ]+)\s+"(.*)"\s*$',bit)
                if m:attrs[m[1]]=m[2]
            genes.add(attrs.get('gene_id'));records+=1
            row=(z[0],z[2],int(z[3]),int(z[4]),z[6],z[7],attrs)
            if 'transcript_id' in attrs:transcript[attrs['transcript_id']].append(row)
            if attrs.get('gene_id')==TC:tcfeatures.append(row)
    check('whole-gene frozen scope counts',len(genes)==prov['frozen_gene_records'] and records==prov['frozen_feature_records'],f'{len(genes)} genes / {records} feature records')
    expected=[];transcript_rows=[];tcintrons=set()
    for tid,rows in transcript.items():
        a=rows[0][6];exons=sorted([r for r in rows if r[1]=='exon'],key=lambda r:(r[2],r[3]));cds=sorted((r[2],r[3],r[5]) for r in rows if r[1]=='CDS');uts=[r for r in rows if r[1]=='five_prime_utr'];minus=rows[0][4]=='-';order=exons[::-1] if minus else exons;ranks={(r[2],r[3]):i+1 for i,r in enumerate(order)}
        transcript_rows.append(dict(transcript_id=tid,n_exons=len(exons),CDS_bases=sum(b-a+1 for a,b,_ in cds),CDS_structure=json.dumps(cds),five_prime_UTR_bases=sum(r[3]-r[2]+1 for r in uts),CCDS_ids=';'.join(sorted({r[6]['ccds_id'] for r in rows if 'ccds_id' in r[6]}))))
        def region(pos):
            for feature,label in [('CDS','CDS'),('five_prime_utr','5prime_UTR'),('three_prime_utr','3prime_UTR')]:
                if any(r[1]==feature and r[2]<=pos<=r[3] for r in rows):return label
            return 'no_CDS_annotation' if not cds else 'exonic_nonCDS_unclassified'
        for low,high in zip(exons,exons[1:]):
            start,end=low[3]+1,high[2]-1
            if end<start:continue
            j=f'{low[0]}:{2 if minus else 1}:{start}:{end}'
            if a.get('gene_id')==TC:tcintrons.add(j)
            if j not in needed:continue
            donor,acceptor=(high,low) if minus else (low,high);db=donor[2] if minus else donor[3];ab=acceptor[3] if minus else acceptor[2]
            expected.append(dict(junction_id=j,gene_id=a['gene_id'],gene_symbol=a.get('gene_name'),transcript_id=tid,transcript_name=a.get('transcript_name'),transcript_biotype=a.get('transcript_biotype'),strand=donor[4],donor_exonic_boundary=db,acceptor_exonic_boundary=ab,donor_exon_start=donor[2],donor_exon_end=donor[3],acceptor_exon_start=acceptor[2],acceptor_exon_end=acceptor[3],donor_exon_rank=ranks[(donor[2],donor[3])],acceptor_exon_rank=ranks[(acceptor[2],acceptor[3])],donor_boundary_region=region(db),acceptor_boundary_region=region(ab)))
    expected=pd.DataFrame(expected).drop_duplicates();actual=pd.read_csv(T/'junction_transcript_memberships.csv.gz');cols=list(expected.columns);keys=['junction_id','gene_id','transcript_id']
    pd.testing.assert_frame_equal(expected[cols].fillna('').sort_values(keys).reset_index(drop=True),actual[cols].fillna('').sort_values(keys).reset_index(drop=True),check_dtype=False)
    check('every exact transcript membership and endpoint classification',True,f'{len(expected)} memberships for {expected.junction_id.nunique()} needed junctions')
    expected_tx=pd.DataFrame(transcript_rows).set_index('transcript_id');actual_tx=pd.read_csv(T/'transcript_context_universe.csv.gz').set_index('transcript_id').loc[expected_tx.index]
    pd.testing.assert_frame_equal(expected_tx,actual_tx[expected_tx.columns].fillna(''),check_dtype=False);check('every transcript exon/CDS/phase/UTR summary',True,f'{len(expected_tx)} transcripts')
    gm=expected.groupby('junction_id').gene_id.apply(set).to_dict();tm=expected.groupby('junction_id').transcript_id.apply(set).to_dict();lookup=alt[['cluster_id','junction_id']].drop_duplicates().groupby('cluster_id').junction_id.apply(list).to_dict();ca=pd.read_csv(T/'all923_cluster_annotation_audit.csv').set_index('cluster_id');n_all=0
    for row in defs.itertuples():
        js=lookup[row.cluster_id];common=set.intersection(*[gm.get(j,set()) for j in js]);r=ca.loc[row.cluster_id];allmatch=all(j in tm for j in js);n_all+=allmatch
        assert len(common)==r.shared_gene_count==row.gene_mapping_count and allmatch==r.junctions_all_exactly_match_transcript_introns
        assert ';'.join(sorted(common))==(r.shared_gene_intersection if pd.notna(r.shared_gene_intersection) else '')
    check('all923 exact common-gene intersections and matching gates',len(ca)==923,f'{n_all} clusters have every alternative matched; others remain in analysis')
    original_local=pd.read_csv(args.base/'extensions/multiway_splicing/tables/local_exact_results_complete.csv.gz');corrected=set(original_local.loc[original_local.BY_q_full_3E<=.05,'cluster_id']);check('all41 original corrected cluster flags preserved',len(corrected)==41 and set(ca.loc[ca.BY_corrected_any_contrast].index)==corrected)
    for row in events.itertuples():assert gm.get(row.junction_a,set())&gm.get(row.junction_b,set())=={row.gene_id},row.event_id
    check('all3239 binary common-gene annotations',len(events)==3239)
    ga=pd.read_csv(T/'all1983_gene_annotation_audit.csv');check('all1983 genes including blank names retained',len(ga)==1983 and ga.gene_id.nunique()==1983 and ga.gene_symbol.isna().sum()==3 and set(ga.gene_id)==set(events.gene_id))
    original_genes=pd.read_csv(args.base/'extensions/gene_omnibus/results/arbitrary_dependence_corrected_genes.csv');check('all47 original corrected gene flags preserved',ga.primary_BY_gene.sum()==47 and set(ga.loc[ga.primary_BY_gene,'gene_id'])==set(original_genes.gene_id))
    unnamed=ga[ga.gene_symbol.isna()];unnamed.to_csv(T/'unnamed_gene_identity_audit.csv',index=False)
    status=pd.DataFrame([dict(junction_id=j,exact_transcript_match=j in tm,n_matching_transcripts=len(tm.get(j,set())),n_matching_genes=len(gm.get(j,set())),matching_gene_ids=';'.join(sorted(gm.get(j,set())))) for j in sorted(set(alt.junction_id))]);status.to_csv(T/'multiway_junction_exact_match_status.csv',index=False)
    check('all3248 multiway junctions and1199 unmatched retained',len(status)==3248 and (~status.exact_transcript_match).sum()==1199,f'{status.exact_transcript_match.sum()} exact matches; STAR annotated flag is a different criterion')
    tc=actual[(actual.gene_id==TC)&actual.junction_id.isin(lookup[CID])];pc=actual_tx[actual_tx.transcript_name.isin(['Tcaf1-201','Tcaf1-202','Tcaf1-203'])];check('three coding transcripts same coordinates/phase/CCDS',len(pc)==3 and pc.CDS_structure.nunique()==1 and pc.CDS_bases.eq(2766).all() and pc.CCDS_ids.eq('CCDS20076').all())
    check('205 retained-intron biotype without CDS',actual_tx[actual_tx.transcript_name.eq('Tcaf1-205')].CDS_bases.eq(0).all() and actual_tx[actual_tx.transcript_name.eq('Tcaf1-205')].transcript_biotype.eq('retained_intron').all())
    check('coding transcript measured boundaries are5primeUTR',tc[tc.transcript_name.isin(pc.transcript_name)][['donor_boundary_region','acceptor_boundary_region']].eq('5prime_UTR').all().all())
    check('all matching Tc af1 transcripts retained',set(tc.transcript_name)=={'Tcaf1-201','Tcaf1-202','Tcaf1-203','Tcaf1-205','Tcaf1-206','Tcaf1-208'},'206 retained-intron and208 CDS-not-defined are not ignored')
    counts=pd.read_csv(args.base/'inputs/all_animal_junction_counts.csv.gz',index_col=0);meta=pd.read_csv(args.base/'inputs/animal_metadata.csv').set_index('sample');support=pd.read_csv(T/'Tcaf1_all26_alternative_read_support.csv');js=list(lookup[CID]);selected=counts.loc[js].sum(axis=0)
    for row in support.itertuples():
        assert row.count==counts.loc[row.junction_id,row.sample] and abs(row.proportion-row.count/selected[row.sample])<1e-14 and row.group==meta.loc[row.sample,'group']
    check('all78 Tc af1 animal-alternative counts/fractions',len(support)==78 and set(zip(support.junction_id,support['sample']))=={(j,sample) for j in js for sample in counts.columns})
    animal=pd.read_csv(T/'Tcaf1_all26_long_short_intron_support.csv').set_index('sample');short='6:2:42663893:42669866';middle='6:2:42663893:42686462';long='6:2:42663893:42686892'
    for sample in counts.columns:
        r=animal.loc[sample];sv=counts.loc[short,sample];mv=counts.loc[middle,sample];lv=counts.loc[long,sample];assert r.selected_reads==sv+mv+lv and r.short5974_reads==sv and r.middle22570_reads==mv and r.long23000_reads==lv and abs(r.long_vs_short_conditional_fraction-lv/(lv+sv))<1e-14
    check('all26 total and conditional long/short fractions',len(animal)==26)
    grouped=pd.read_csv(T/'Tcaf1_group_read_support_summary.csv').set_index('group')
    for g in ['FLT','GC','VIV']:
        subset=animal[animal.group==g];row=grouped.loc[g];assert row.n==len(subset) and row.minimum_selected_reads==subset.selected_reads.min() and row.median_selected_reads==subset.selected_reads.median() and abs(row.mean_long_fraction-subset.long23000_fraction.mean())<1e-14 and abs(row.mean_long_vs_short_conditional-subset.long_vs_short_conditional_fraction.mean())<1e-14
    check('all three arm empirical count/fraction summaries',True)
    source=counts.loc[counts.index.isin(tcintrons)];stored=pd.read_csv(T/'Tcaf1_all_annotated_intron_read_counts.csv.gz').pivot(index='junction_id',columns='sample',values='unique_reads').loc[source.index,source.columns];check('every available annotated Tc af1 intron count',np.array_equal(source.to_numpy(),stored.to_numpy()),f'{len(source)} introns ×26 animals')
    cp=counts.loc['6:2:42686572:42686892'];check('companion junction in everyanimal20to169reads',len(cp)==26 and cp.min()==20 and cp.max()==169)
    # Preserve the original snapshot provenance, but do not confuse gzip-byte drift with data drift.
    drift=[]
    for path,oldhash in prov['base_input_sha256'].items():drift.append(dict(package_relative_path=path,original_snapshot_sha256=oldhash,current_sha256=sha(args.base/path),byte_hash_equal=sha(args.base/path)==oldhash,semantic_verification='independent memberships, definitions, gene IDs and raw animal counts checked above; gzip byte equality not required'))
    pd.DataFrame(drift).to_csv(T/'derived_input_snapshot_byte_check.csv',index=False)
    result=dict(all_passed=True,n_checks=len(checks),checks=checks,scope='independent reparse of frozen whole-gene GTF; original full GTF/FASTA sequences not reconstituted',original_provenance_preserved=True);(T/'independent_verification.json').write_text(json.dumps(result,indent=2));pd.DataFrame(checks).to_csv(T/'independent_verification.csv',index=False);print('ALL CHECKS PASS',len(checks))

if __name__=='__main__':main()
