#!/usr/bin/env python3
"""Preserve the primary 1,983-gene universe under a stricter dominance gate.

Genes losing all events receive p=1; no family reduction or outcome selection.
This sensitivity is secondary and cannot replace the 90% primary analysis.
"""
import argparse
import json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import analyse_gene_omnibus as a


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=a.Path,default=a.DEFAULT_BASE)
    args=parser.parse_args()
    selected,x,metadata,genes=a.load(args.base)
    gates=pd.read_csv(args.base/'results/measurement_gates.csv.gz').set_index('event_id')
    keep=selected.event_id.map(gates.measured_dominance95).to_numpy(bool)
    selected95=selected[keep].reset_index(drop=True)
    x95=x[keep]
    genes95=selected95.groupby('gene_id',sort=True).indices
    primary=pd.read_csv(a.HERE/'results/gene_omnibus_complete.csv.gz')
    all_rows=[]
    for ga,gb in a.CONTRASTS:
        result=a.exact_contrast(x95,metadata,genes95,ga,gb)
        # Empty measured genes are conservative p=1 members of the unchanged family.
        template=primary[primary.contrast==f'{ga}_minus_{gb}'][
            ['gene_id','gene_symbol','contrast','n_group_a','n_group_b','null_allocations']].copy()
        result=template.merge(result,on=['gene_id','contrast','n_group_a','n_group_b','null_allocations'],how='left')
        absent=result.p_exact.isna()
        result.loc[absent,'p_exact']=1
        result.loc[absent,'null_exceedances']=result.loc[absent,'null_allocations']
        result.loc[absent,'event_count']=0
        result['measured_gene_dominance95']=~absent
        all_rows.append(result)
    result=pd.concat(all_rows,ignore_index=True)
    for method in ['bh','by']:
        result[f'q_{method.upper()}_joint_3G']=false_discovery_control(result.p_exact,method=method)
    result.to_csv(a.HERE/'results/dominance95_gene_omnibus_complete.csv.gz',index=False)
    summary=[]
    for contrast,d in result.groupby('contrast',sort=False):
        old=primary[primary.contrast==contrast]
        old_hits=set(old[old.q_BY_joint_3G<.05].gene_id)
        hits=set(d[d.q_BY_joint_3G<.05].gene_id)
        summary.append({'contrast':contrast,'family_genes':len(d),'measured_genes':int(d.measured_gene_dominance95.sum()),
            'measured_events':int(d.event_count.sum()),'joint_BH':int((d.q_BH_joint_3G<.05).sum()),
            'joint_BY':len(hits),'primary_BY_genes_retaining_BY':len(hits&old_hits),
            'primary_BY_genes_total':len(old_hits),'minimum_joint_BY_q':float(d.q_BY_joint_3G.min())})
    pd.DataFrame(summary).to_csv(a.HERE/'results/dominance95_sensitivity_summary.csv',index=False)
    primary_hits=primary[primary.q_BY_joint_3G<.05][['gene_id','gene_symbol','contrast','event_count','q_BY_joint_3G']]
    follow=primary_hits.merge(result[['gene_id','contrast','event_count','p_exact','q_BY_joint_3G','measured_gene_dominance95']],
        on=['gene_id','contrast'],suffixes=('_primary90','_sensitivity95'))
    follow.to_csv(a.HERE/'results/primary_BY_gene_dominance95_followup.csv',index=False)
    assert len(result)==len(primary) and result.p_exact.notna().all()
    manifest=json.loads((a.HERE/'results/manifest.json').read_text())
    manifest['dominance95_sensitivity']={'selection':'Stricter all26 label-invariant dominance; original 1983-gene family held fixed',
        'zero_event_genes':'p=1 retained in joint 5949-test correction','summary':summary,
        'script_sha256':a.sha(a.HERE/'scripts/sensitivity_dominance95.py')}
    (a.HERE/'results/manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
