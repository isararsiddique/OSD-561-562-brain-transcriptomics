"""Freeze required gene records and publisher reference provenance; initial step."""
import argparse,gzip,json,re,shutil,hashlib,zipfile
from pathlib import Path
import pandas as pd
from context_common import HERE,DEFAULT_BASE,sha,gene_groups,transcript_introns,universe

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=DEFAULT_BASE)
    p.add_argument('--gtf',type=Path,required=True);p.add_argument('--publisher-dir',type=Path,required=True)
    args=p.parse_args();(HERE/'inputs').mkdir(exist_ok=True)
    source=json.loads((args.base/'extensions/multiway_splicing/annotation/source.json').read_text())
    assert sha(args.gtf)==source['gtf_sha256']
    multi,alternatives,genes,needed=universe(args.base)
    required_genes=set(genes.gene_id)
    count_genes=count_records=0;covered=set()
    with gzip.open(HERE/'inputs/required_gene_context.GRCm39.112.gtf.gz','wt') as out:
        out.write('# Frozen whole-gene records from GRCm39.112 for the pooled observable universe.\n')
        for gene,group in gene_groups(args.gtf):
            matched={jid for jid,*_ in transcript_introns(group) if jid in needed}
            if gene in required_genes or matched:
                count_genes+=1;count_records+=len(group);covered.update(matched)
                out.writelines(row['line'] for row in group)
    print('Frozen genes/records/matching junctions:',count_genes,count_records,len(covered),flush=True)
    sidecar=args.publisher_dir/'GLDS-525_rna_seq_processed_md5sum_GLbulkRNAseq.tsv'
    checksum=pd.read_csv(sidecar,sep='\t',header=None,names=['file','md5']).set_index('file')['md5']
    verified=[]
    for suffix in ['processing_info_GLbulkRNAseq.zip','software_versions_GLbulkRNAseq.md']:
        original=args.publisher_dir/f'GLDS-525_rna_seq_{suffix}'
        expected=checksum.loc[suffix]
        actual=hashlib.md5(original.read_bytes()).hexdigest();assert actual==expected
        shutil.copy2(original,HERE/'inputs'/suffix)
        verified.append({'source_filename':original.name,'publisher_filename':suffix,'publisher_md5':expected,
                         'verified':True,'sha256':sha(original)})
    shutil.copy2(sidecar,HERE/'inputs/publisher_processed_md5sum.tsv')
    with zipfile.ZipFile(HERE/'inputs/processing_info_GLbulkRNAseq.zip') as z:
        processing=z.read('nextflow_processing_info_GLbulkRNAseq.txt').decode()
        command=z.read('nextflow_run_command_GLbulkRNAseq.txt').decode()
    (HERE/'inputs/publisher_processing_reference_extract.txt').write_text('\n'.join(
        line for line in processing.splitlines() if 'BUILD_STAR_INDEX' in line or '--reference_gtf' in line
        or '--sjdbGTFfile' in line or ('STAR --twopassMode' in line))+'\n')
    assert 'GRCm39.112' in processing and '--twopassMode Basic' in processing
    fasta=sorted(set(re.findall(r'Mus_musculus\.GRCm39[^\s]+\.fa',processing)))
    manifest={'status':'Reference assembly and annotation-release names agree with publisher processing',
        'publisher_reference':'Mus_musculus.GRCm39 primary assembly + ERCC92; Ensembl112',
        'publisher_STAR_version':'2.7.11b','publisher_twopass_mode':'Basic','publisher_read_length_max':151,
        'gtf_original_sha256':sha(args.gtf),'gtf_source_url':source['url'],'gtf_source_filename':args.gtf.name,
        'frozen_gene_records':count_genes,'frozen_feature_records':count_records,'needed_junctions':len(needed),
        'junctions_with_exact_transcript_intron_match':len(covered),'publisher_files':verified,
        'annotation_scope':'All1983 measured gene candidates and every gene with an exact match to any junction in all923 multialternative clusters or all3239 measured uniquely gene-mapped binary choices; whole gene records retained.',
        'selection_uses_significance':False,'frozen_gtf_sha256':sha(HERE/'inputs/required_gene_context.GRCm39.112.gtf.gz'),
        'base_input_sha256':{str(f):sha(args.base/f) for f in [
            Path('extensions/multiway_splicing/tables/selected_cluster_definitions.csv'),
            Path('extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz'),
            Path('extensions/gene_omnibus/results/gene_event_universe.csv.gz')]},
        'limits':['Reference names and GTF source hash are verified; source FASTA sequence itself is not independently reconstructed.',
                  'Read counts are reused short-read data; annotation context does not identify a complete isoform or protein consequence.']}
    (HERE/'inputs/provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':main()
