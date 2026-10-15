#!/usr/bin/env python3
"""Exact, dependence-preserving gene-level local junction usage contrasts.

Selection uses pooled all-animal measurement and unique intron annotation only.
No event significance, direction, group-specific effect, or gene ontology enters
the gene universe. Each permutation moves an animal's whole event vector.
"""
from pathlib import Path
from itertools import combinations
import argparse
import hashlib
import json
import math
import time
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

HERE = Path(__file__).resolve().parents[1]
def resolve_default_base():
    """Support both research/gene extension and beta/extensions/gene_omnibus."""
    relative=Path('inputs/primary_all26_event_definitions.csv.gz')
    for parent in [HERE,*HERE.parents]:
        if (parent/relative).is_file():
            return parent
        sibling=parent/'paper2_beta_20261014'
        if (sibling/relative).is_file():
            return sibling
    raise FileNotFoundError('Frozen beta inputs not found; pass --base to their root.')


DEFAULT_BASE = resolve_default_base()
CONTRASTS = [('FLT', 'GC'), ('GC', 'VIV'), ('FLT', 'VIV')]
TOL = 1e-14  # Absolute numerical tie tolerance in squared fraction units.


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load(base):
    definitions = pd.read_csv(base / 'inputs/primary_all26_event_definitions.csv.gz')
    gates = pd.read_csv(base / 'results/measurement_gates.csv.gz').set_index('event_id')
    fractions = pd.read_csv(base / 'results/all26_animal_fractions.csv.gz').set_index('event_id')
    metadata = pd.read_csv(base / 'inputs/animal_metadata.csv')
    selected = definitions.loc[
        definitions.event_id.map(gates.measured_primary_all26)
        & definitions.gene_mapping_count.eq(1)
        & definitions.gene_id.notna()
    ].copy().sort_values(['gene_id', 'event_id']).reset_index(drop=True)
    x = fractions.loc[selected.event_id, metadata['sample']].to_numpy(float)
    assert np.isfinite(x).all() and ((x >= 0) & (x <= 1)).all()
    genes = selected.groupby('gene_id', sort=True).indices
    return selected, x, metadata, genes


def kernels(x, indices, genes):
    """One centered mean event Gram matrix per gene; preserves event dependence."""
    n = len(indices)
    upper = np.triu_indices(n)
    multiplicity = np.where(upper[0] == upper[1], 1., 2.)
    flat = np.empty((len(upper[0]), len(genes)), dtype=float)
    observed = np.empty(len(genes))
    raw = x[:, indices]
    for j, rows in enumerate(genes.values()):
        z = raw[rows]
        z = z - z.mean(axis=1, keepdims=True)
        k = z.T @ z / len(rows)
        flat[:, j] = k[upper] * multiplicity
    return raw, flat, upper


def exact_contrast(x, metadata, genes, group_a, group_b):
    indices_a = np.flatnonzero(metadata.group.to_numpy() == group_a)
    indices_b = np.flatnonzero(metadata.group.to_numpy() == group_b)
    indices = np.r_[indices_a, indices_b]
    na, nb, n = len(indices_a), len(indices_b), len(indices)
    raw, flat, upper = kernels(x, indices, genes)
    observed = np.array([np.mean((raw[rows, :na].mean(axis=1)
                                  - raw[rows, na:].mean(axis=1)) ** 2)
                         for rows in genes.values()])
    deltas = np.array([np.mean(raw[rows, :na].mean(axis=1)
                                - raw[rows, na:].mean(axis=1))
                       for rows in genes.values()])
    unbiased = []
    for rows in genes.values():
        z = raw[rows]
        difference = z[:, :na].mean(axis=1) - z[:, na:].mean(axis=1)
        unbiased.append(np.mean(difference ** 2 - z[:, :na].var(axis=1, ddof=1) / na
                                - z[:, na:].var(axis=1, ddof=1) / nb))
    allocations = np.array(list(combinations(range(n), na)), dtype=np.int16)
    exceed = np.zeros(len(genes), dtype=np.int64)
    start = time.monotonic()
    for offset in range(0, len(allocations), 512):
        subset = allocations[offset:offset+512]
        weights = np.full((len(subset), n), -1 / nb, dtype=float)
        weights[np.arange(len(subset))[:, None], subset] = 1 / na
        outer = weights[:, upper[0]] * weights[:, upper[1]]
        permuted = outer @ flat
        exceed += (permuted >= observed[None, :] - TOL).sum(axis=0)
    assert (exceed >= 1).all() and (exceed <= len(allocations)).all()
    print(f'{group_a}-{group_b}: {len(genes)} genes, {len(allocations)} allocations, '
          f'{time.monotonic()-start:.1f} seconds', flush=True)
    return pd.DataFrame({
        'gene_id': list(genes), 'contrast': f'{group_a}_minus_{group_b}',
        'n_group_a': na, 'n_group_b': nb,
        'event_count': [len(v) for v in genes.values()],
        'rms_delta_percentage_points': np.sqrt(observed) * 100,
        'signed_mean_delta_percentage_points': deltas * 100,
        'mean_squared_delta_fraction_units': observed,
        'unbiased_mean_squared_distance_fraction_units': unbiased,
        'null_allocations': len(allocations), 'null_exceedances': exceed,
        'p_exact': exceed / len(allocations)
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base', type=Path, default=DEFAULT_BASE)
    parser.add_argument('--output', type=Path, default=HERE)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'results').mkdir(exist_ok=True)
    selected, x, metadata, genes = load(args.base)
    print(f'Gene universe: {len(genes)} genes, {len(selected)} measured events, '
          f'{len(metadata)} animals', flush=True)
    selected.to_csv(args.output/'results/gene_event_universe.csv.gz', index=False)
    results = pd.concat([exact_contrast(x, metadata, genes, a, b)
                         for a, b in CONTRASTS], ignore_index=True)
    symbol = selected.groupby('gene_id').gene_symbol.first()
    results.insert(1, 'gene_symbol', results.gene_id.map(symbol))
    results['q_BH_joint_3G'] = false_discovery_control(results.p_exact, method='bh')
    results['q_BY_joint_3G'] = false_discovery_control(results.p_exact, method='by')
    for method in ['bh', 'by']:
        results[f'q_{method.upper()}_within_G'] = results.groupby('contrast').p_exact.transform(
            lambda p: false_discovery_control(p, method=method))
    results.to_csv(args.output/'results/gene_omnibus_complete.csv.gz', index=False)
    summaries = []
    for contrast, d in results.groupby('contrast', sort=False):
        summaries.append({'contrast': contrast, 'tested_genes': len(d),
                          'measured_events': int(d.event_count.sum()),
                          'joint_BH_q_lt_0_05': int((d.q_BH_joint_3G < .05).sum()),
                          'joint_BY_q_lt_0_05': int((d.q_BY_joint_3G < .05).sum()),
                          'within_BH_q_lt_0_05': int((d.q_BH_within_G < .05).sum()),
                          'within_BY_q_lt_0_05': int((d.q_BY_within_G < .05).sum()),
                          'minimum_joint_BY_q': float(d.q_BY_joint_3G.min()),
                          'median_gene_rms_pp': float(d.rms_delta_percentage_points.median())})
    pd.DataFrame(summaries).to_csv(args.output/'results/discovery_summary.csv', index=False)
    significant = results.loc[results.q_BY_joint_3G < .05].sort_values(
        ['q_BY_joint_3G', 'p_exact', 'gene_id'])
    significant.to_csv(args.output/'results/arbitrary_dependence_corrected_genes.csv', index=False)
    inputs = ['inputs/primary_all26_event_definitions.csv.gz',
              'results/measurement_gates.csv.gz', 'results/all26_animal_fractions.csv.gz',
              'inputs/animal_metadata.csv']
    manifest = {'method': 'Exact animal-vector gene omnibus pairwise contrast',
                'selection': 'All26 pooled measurement gate and unique intron annotation; no event significance',
                'genes': len(genes), 'events': len(selected), 'animals': len(metadata),
                'joint_test_family': len(results), 'seed': None,
                'numerical_tie_tolerance_squared_fraction': TOL,
                'input_sha256': {p: sha(args.base/p) for p in inputs},
                'analysis_script_sha256': sha(Path(__file__)),
                'limitations': ['Flight and collection date perfectly aliased.',
                    'Pairwise labels assumed exchangeable under the measured-association null.',
                    'Measured local junction fractions are not complete isoforms or gene expression.',
                    'Multi-event gene tests average correlated local-choice effects without independence assumptions.',
                    'Gene/event annotation and measurement filters define a restricted observable gene universe.',
                    'New analyses are retrospective, not preregistered.'],
                'summary': summaries}
    (args.output/'results/manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(summaries, indent=2), flush=True)


if __name__ == '__main__':
    main()
