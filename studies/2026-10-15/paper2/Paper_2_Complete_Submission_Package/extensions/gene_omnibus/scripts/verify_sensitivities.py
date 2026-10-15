#!/usr/bin/env python3
"""Independent all-rows correction and sampled direct sensitivity checks."""
import argparse
import json
import numpy as np
import pandas as pd
import analyse_gene_omnibus as a
from verify_gene_omnibus import adjustment,direct_p


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=a.Path,default=a.DEFAULT_BASE)
    args=parser.parse_args()
    selected,x,metadata,genes=a.load(args.base)
    gates=pd.read_csv(args.base/'results/measurement_gates.csv.gz').set_index('event_id')
    keep=selected.event_id.map(gates.measured_dominance95).to_numpy(bool)
    s95=selected[keep].reset_index(drop=True)
    x95=x[keep]
    g95=s95.groupby('gene_id',sort=True).indices
    result=pd.read_csv(a.HERE/'results/dominance95_gene_omnibus_complete.csv.gz')
    table=pd.read_csv(a.HERE/'results/descriptive_all_gene_leave_one_out.csv.gz')
    checks=[]
    def check(name,condition,detail=None):
        assert condition,name
        checks.append({'name':name,'passed':True,'detail':detail})
    check('95% universe is label-invariant measurement restriction',len(s95)==2710 and len(g95)==1724)
    check('Original 5949-test family retained',len(result)==5949 and result.groupby('contrast').gene_id.nunique().eq(1983).all())
    missing=result[~result.measured_gene_dominance95]
    check('All unmeasured genes retained at p=1',len(missing)==3*(1983-1724) and missing.p_exact.eq(1).all())
    check('Independent manual BH/BY reproduce sensitivity',np.allclose(adjustment(result.p_exact),result.q_BH_joint_3G)
          and np.allclose(adjustment(result.p_exact,True),result.q_BY_joint_3G))
    for contrast,d in result.groupby('contrast'):
        ga,gb=contrast.split('_minus_')
        ia=np.flatnonzero(metadata.group.to_numpy()==ga)
        ib=np.flatnonzero(metadata.group.to_numpy()==gb)
        for symbol in ['Hnrnpa2b1','Rps6kb1']:
            row=d[d.gene_symbol==symbol].iloc[0]
            z=x95[g95[row.gene_id]][:,np.r_[ia,ib]]
            p,count=direct_p(z,len(ia))
            check(f'{contrast} {symbol}: direct sensitivity allocation',
                  count==int(row.null_exceedances) and abs(p-row.p_exact)<1e-12)
    check('All 1983 genes under all 18 animal omissions retained',len(table)==1983*18 and table.groupby('gene_id').size().eq(18).all())
    check('No leave-one-out p or q values computed',not any(c.startswith(('p_','q_')) for c in table.columns))
    # A pandas group-mean path, independent of production's indexed array means.
    frame=pd.DataFrame(x,index=selected.event_id,columns=metadata['sample'])
    chosen=table.sort_values('rms_ratio_to_full').iloc[[0,10,100,1000,10000,-1]]
    for _,row in chosen.iterrows():
        events=selected.iloc[genes[row.gene_id]].event_id
        samples_a=metadata.loc[metadata.group.eq('FLT')&metadata['sample'].ne(row.omitted_sample),'sample']
        samples_b=metadata.loc[metadata.group.eq('GC')&metadata['sample'].ne(row.omitted_sample),'sample']
        delta=frame.loc[events,samples_a].mean(axis=1)-frame.loc[events,samples_b].mean(axis=1)
        rms=float(np.sqrt((delta**2).mean())*100)
        check(f'LOO {row.gene_id}/{row.omitted_sample}: independent pandas means',abs(rms-row.loo_rms_pp)<1e-10)
    output={'checks_passed':len(checks),'checks':checks,
        'sensitivity_result_sha256':a.sha(a.HERE/'results/dominance95_gene_omnibus_complete.csv.gz'),
        'loo_result_sha256':a.sha(a.HERE/'results/descriptive_all_gene_leave_one_out.csv.gz')}
    (a.HERE/'results/sensitivity_verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    main()
