from pathlib import Path
import gzip,hashlib,re
import pandas as pd

HERE=Path(__file__).resolve().parents[1]
DEFAULT_BASE=HERE.parents[1]
TC_CLUSTER='6:2:42663893:start'
ATTR=re.compile(r'(\S+) "([^"]*)"')

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def records(path):
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            z=line.rstrip('\n').split('\t')
            if len(z)!=9:continue
            attrs=dict(ATTR.findall(z[8]))
            yield {'chrom':z[0],'feature':z[2],'start':int(z[3]),'end':int(z[4]),
                   'strand':z[6],'phase':z[7],'attrs':attrs,'line':line}

def gene_groups(path):
    current=None;group=[]
    for row in records(path):
        gene=row['attrs'].get('gene_id')
        if gene!=current and group:
            yield current,group
            group=[]
        current=gene;group.append(row)
    if group:yield current,group

def transcript_introns(group):
    transcripts={}
    for row in group:
        tid=row['attrs'].get('transcript_id')
        if tid:transcripts.setdefault(tid,[]).append(row)
    for tid,rows in transcripts.items():
        exons=sorted([r for r in rows if r['feature']=='exon'],key=lambda r:(r['start'],r['end']))
        if not exons:continue
        minus=exons[0]['strand']=='-'
        ordered=list(reversed(exons)) if minus else exons
        ranks={(e['start'],e['end']):i+1 for i,e in enumerate(ordered)}
        for low,high in zip(exons,exons[1:]):
            first,last=low['end']+1,high['start']-1
            if last<first:continue
            donor,acceptor=(high,low) if minus else (low,high)
            jid=f"{low['chrom']}:{2 if minus else 1}:{first}:{last}"
            yield jid,tid,rows,donor,acceptor,ranks

def universe(base):
    multi=pd.read_csv(base/'extensions/multiway_splicing/tables/selected_cluster_definitions.csv')
    alternatives=pd.read_csv(base/'extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz')
    genes=pd.read_csv(base/'extensions/gene_omnibus/results/gene_event_universe.csv.gz')
    needed=set(alternatives.junction_id)|set(genes.junction_a)|set(genes.junction_b)
    return multi,alternatives,genes,needed
