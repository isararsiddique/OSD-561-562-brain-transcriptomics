#!/usr/bin/env python3
"""Independent direct-allocation p values and manual multiple-test adjustment."""
import argparse
import json
from pathlib import Path
from itertools import combinations
import math
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import analyse_gene_omnibus as a


def adjustment(p, by=False):
    p = np.asarray(p)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p)+1)
    if by:
        ranked *= np.sum(1 / np.arange(1, len(p)+1))
    result = np.empty(len(p))
    result[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return result


def direct_p(z, na):
    """No Gram kernel: explicitly allocate and subtract each event's means."""
    n, nb = z.shape[1], z.shape[1] - na
    obs = np.mean((z[:, :na].mean(axis=1) - z[:, na:].mean(axis=1))**2)
    allocations = np.array(list(combinations(range(n), na)), dtype=np.int16)
    total = z.sum(axis=1)
    count = 0
    for offset in range(0, len(allocations), 1024):
        selected = allocations[offset:offset+1024]
        sums_a = z[:, selected].sum(axis=2).T
        delta = sums_a / na - (total[None, :] - sums_a) / nb
        statistic = (delta**2).mean(axis=1)
        count += int((statistic >= obs - a.TOL).sum())
    return count / len(allocations), count


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=Path,default=a.DEFAULT_BASE)
    args=parser.parse_args()
    selected, x, metadata, genes = a.load(args.base)
    result = pd.read_csv(a.HERE/'results/gene_omnibus_complete.csv.gz')
    manifest = json.loads((a.HERE/'results/manifest.json').read_text())
    checks = []
    def check(name, condition, detail=None):
        assert condition, name
        checks.append({'name': name, 'passed': True, 'detail': detail})
    check('Input hashes match manifest', all(a.sha(args.base/p)==h for p,h in manifest['input_sha256'].items()))
    check('Unique gene universe selected without outcome significance', len(genes)==1983 and len(selected)==3239)
    check('Complete three-contrast family and unsuccessful tests retained', len(result)==3*len(genes) and not result.p_exact.isna().any())
    check('Every gene occurs once per contrast', result.groupby('contrast').gene_id.nunique().eq(len(genes)).all())
    check('Event count matches whole pooled observable gene universe',
          all(int(d.event_count.sum())==len(selected) for _,d in result.groupby('contrast')))
    check('Joint BH matches independent manual adjustment', np.allclose(adjustment(result.p_exact),result.q_BH_joint_3G, atol=1e-12))
    check('Joint BY matches independent manual adjustment', np.allclose(adjustment(result.p_exact,True),result.q_BY_joint_3G, atol=1e-12))
    check('SciPy BH and BY separately reproduce saved results',
          all(np.allclose(false_discovery_control(result.p_exact,method=m),result[f'q_{m.upper()}_joint_3G']) for m in ['bh','by']))
    # A one-event gene's squared statistic must agree with the earlier absolute
    # event-contrast statistic, without using event significance for selection.
    event_result=pd.read_csv(args.base/'results/three_arm_event_results.csv.gz')
    single=selected.groupby('gene_id').filter(lambda z:len(z)==1)[['gene_id','event_id']]
    matched=result[result.event_count.eq(1)].merge(single,on='gene_id').merge(
        event_result[['event_id','contrast','exact_p']],on=['event_id','contrast'])
    check('All 3567 single-event gene tests reproduce independently fitted event p values',
          len(matched)==3567 and np.array_equal(matched.p_exact,matched.exact_p),
          {'verification_only_event_result_sha256':a.sha(args.base/'results/three_arm_event_results.csv.gz')})
    for contrast, d in result.groupby('contrast'):
        ga, gb = contrast.split('_minus_')
        na = int(d.n_group_a.iloc[0]); nb = int(d.n_group_b.iloc[0])
        check(f'{contrast}: all exact allocations enumerated', d.null_allocations.eq(math.comb(na+nb,na)).all())
        check(f'{contrast}: within-family adjustments independent',
              np.allclose(adjustment(d.p_exact),d.q_BH_within_G) and np.allclose(adjustment(d.p_exact,True),d.q_BY_within_G))
        ia = np.flatnonzero(metadata.group.to_numpy()==ga)
        ib = np.flatnonzero(metadata.group.to_numpy()==gb)
        # Check smallest p, largest multi-event gene, plus ordinary middle and high p.
        choices = list(dict.fromkeys([d.sort_values('p_exact').gene_id.iloc[0],
                       d.sort_values('event_count',ascending=False).gene_id.iloc[0],
                       d.sort_values('p_exact').gene_id.iloc[len(d)//2],
                       d.sort_values('p_exact').gene_id.iloc[-1]]))
        for gene in choices:
            z = x[genes[gene]][:, np.r_[ia,ib]]
            p, count = direct_p(z,na)
            row = d.set_index('gene_id').loc[gene]
            check(f'{contrast} {gene}: independent direct vector allocation',
                  count==int(row.null_exceedances) and abs(p-row.p_exact)<1e-12,
                  {'events':len(z),'p':p,'exceedances':count})
    # Mathematical invariants on small, exhaustive toy data.
    z = np.array([[.1,.2,.15,.8,.7,.9],[.8,.7,.75,.2,.3,.1]])
    p,_ = direct_p(z,3)
    check('Duplicating all correlated event coordinates preserves omnibus p', direct_p(np.repeat(z,2,axis=0),3)[0]==p)
    check('Uniform scaling preserves exact p', direct_p(z*2,3)[0]==p)
    check('Constant animal vectors yield p=1', direct_p(np.repeat(z[:,:1],6,axis=1),3)[0]==1)
    dates=metadata.groupby('group')['Comment[Euthanasia Date]'].apply(set)
    check('Collection-date alias disclosed and not adjusted away', dates['GC']==dates['VIV'] and dates['FLT'].isdisjoint(dates['GC']))
    output = {'checks_passed':len(checks),'checks':checks,
              'analysis_result_sha256':a.sha(a.HERE/'results/gene_omnibus_complete.csv.gz')}
    (a.HERE/'results/verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    main()
