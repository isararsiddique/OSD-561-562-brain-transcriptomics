#!/usr/bin/env python3
"""All-gene leave-one-animal-out effect stability; no p values or selection."""
import argparse
import json
import numpy as np
import pandas as pd
import analyse_gene_omnibus as a


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=a.Path,default=a.DEFAULT_BASE)
    args=parser.parse_args()
    selected,x,metadata,genes=a.load(args.base)
    ia=np.flatnonzero(metadata.group.to_numpy()=='FLT')
    ib=np.flatnonzero(metadata.group.to_numpy()=='GC')
    full_delta=x[:,ia].mean(axis=1)-x[:,ib].mean(axis=1)
    records=[]
    for animal in np.r_[ia,ib]:
        la=ia[ia!=animal];lb=ib[ib!=animal]
        delta=x[:,la].mean(axis=1)-x[:,lb].mean(axis=1)
        for gene,events in genes.items():
            current=delta[events];original=full_delta[events]
            old_rms=np.sqrt(np.mean(original**2))*100
            new_rms=np.sqrt(np.mean(current**2))*100
            active=np.abs(original)>1e-12
            directions=np.sign(original[active])==np.sign(current[active])
            records.append({'gene_id':gene,'omitted_sample':metadata.iloc[animal]['sample'],
                'omitted_group':metadata.iloc[animal]['group'],'event_count':len(events),
                'full_rms_pp':old_rms,'loo_rms_pp':new_rms,
                'rms_ratio_to_full':new_rms/old_rms if old_rms>0 else np.nan,
                'full_signed_mean_delta_pp':original.mean()*100,
                'loo_signed_mean_delta_pp':current.mean()*100,
                'signed_mean_direction_retained':bool(np.sign(original.mean())==np.sign(current.mean())),
                'fraction_event_directions_retained':directions.mean() if len(directions) else 1.})
    table=pd.DataFrame(records)
    table.to_csv(a.HERE/'results/descriptive_all_gene_leave_one_out.csv.gz',index=False)
    primary=pd.read_csv(a.HERE/'results/gene_omnibus_complete.csv.gz')
    baseline=primary[primary.contrast=='FLT_minus_GC'][['gene_id','gene_symbol','event_count','rms_delta_percentage_points','q_BY_joint_3G']]
    stability=table.groupby('gene_id').agg(
        deletions=('omitted_sample','size'),minimum_loo_rms_pp=('loo_rms_pp','min'),
        maximum_loo_rms_pp=('loo_rms_pp','max'),minimum_rms_ratio=('rms_ratio_to_full','min'),
        maximum_rms_ratio=('rms_ratio_to_full','max'),
        all_signed_mean_directions_retained=('signed_mean_direction_retained','all'),
        minimum_fraction_event_directions_retained=('fraction_event_directions_retained','min'))
    stability=baseline.merge(stability,on='gene_id')
    stability['posthoc_primary_BY_subset']=stability.q_BY_joint_3G<.05
    stability.to_csv(a.HERE/'results/descriptive_gene_leave_one_out_summary.csv',index=False)
    hits=stability[stability.posthoc_primary_BY_subset]
    summary={'scope':'All1983 genes / 18 FLT or GC animal omissions; no p values.',
        'primary_BY_genes':len(hits),'all18_signed_mean_directions_retained':int(hits.all_signed_mean_directions_retained.sum()),
        'all18_all_event_directions_retained':int(hits.minimum_fraction_event_directions_retained.eq(1).sum()),
        'minimum_rms_ratio_any_primary_gene_omission':float(hits.minimum_rms_ratio.min()),
        'maximum_rms_ratio_any_primary_gene_omission':float(hits.maximum_rms_ratio.max()),
        'median_minimum_rms_ratio_primary_genes':float(hits.minimum_rms_ratio.median()),
        'worst_rms_attenuation_gene':hits.sort_values('minimum_rms_ratio').iloc[0][['gene_id','gene_symbol','minimum_rms_ratio']].to_dict(),
        'caution':'Descriptive effect stability is not a recalculated p/q, proof of absence of confounding, or independent replication.'}
    (a.HERE/'results/descriptive_leave_one_out_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    assert len(table)==1983*18 and stability.deletions.eq(18).all()
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
