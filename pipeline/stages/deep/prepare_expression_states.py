from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Exploratory, exposure-blind forebrain states from non-panel expression.

Raw integer counts are aggregated by mouse/state. These numerical states are
not validated cell subtypes. No outcome label enters feature selection/PCA.
"""
from pathlib import Path
import gzip, json, hashlib, sys
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score

B = stage_root('deep')
R = stage_root('extended')
depth_adjusted = '--depth-adjusted' in sys.argv
root = B
if depth_adjusted:
    B = root/'depth_sensitivity'
    B.mkdir(exist_ok=True)
for d in ['prepared', 'tables', 'sources', 'figures', 'logs']:
    (B/d).mkdir(exist_ok=True)
md = pd.read_csv(R/'tables/GSE295428_current_sample_registry.csv')
md = md[md.tissue.eq('Myeloid_Forebrain')].sort_values('sample')
panel = pd.read_csv(stage_root('core') / 'tables/verified_panel_mapping.csv')
oldqc = pd.read_csv(R/'tables/GSE295428_sorted_cell_QC.csv').set_index('sample')
xs, cells, hashes = [], [], []
ann = None
for s in md.itertuples():
    folder = R/'sources/GSE295428_sorted_cells'
    paths = {suffix: next(folder.glob(s.sample+'*_'+suffix)) for suffix in
             ['matrix.mtx.gz', 'features.tsv.gz', 'barcodes.tsv.gz']}
    a = pd.read_csv(paths['features.tsv.gz'], sep='\t', header=None,
                    names=['ENSEMBL', 'SYMBOL', 'type'])
    bars = pd.read_csv(paths['barcodes.tsv.gz'], sep='\t', header=None)[0]
    with gzip.open(paths['matrix.mtx.gz'], 'rb') as f:
        x = mmread(f).tocsc()
    assert x.shape == (len(a), len(bars)) and (x.data >= 0).all()
    assert np.equal(x.data, np.floor(x.data)).all()
    if ann is None:
        ann = a
    else:
        assert a.equals(ann)
    lib = np.asarray(x.sum(axis=0)).ravel()
    nfeat = np.diff(x.indptr)
    mt = np.asarray(x[a.SYMBOL.str.lower().str.startswith('mt-').to_numpy()].sum(axis=0)).ravel()/np.maximum(lib, 1)
    keep = (nfeat >= 200) & (nfeat <= 2500) & (mt <= .5)
    assert keep.sum() == oldqc.loc[s.sample, 'n_retained_cells']
    xs.append(x[:, keep].T.tocsr())
    cells.append(pd.DataFrame({'sample': s.sample, 'donor': s.donor,
                              'barcode': bars[keep].to_numpy(), 'library_UMI': lib[keep],
                              'n_features': nfeat[keep], 'mt_fraction': mt[keep]}))
    for p in paths.values():
        hashes.append({'path': record_path(p), 'bytes': p.stat().st_size,
                       'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
    print(s.sample, int(keep.sum()), 'retained', flush=True)
x = sparse.vstack(xs).tocsr()
obs = pd.concat(cells, ignore_index=True)
norm = x.astype(np.float32).multiply((1e4/obs.library_UMI.to_numpy())[:, None]).tocsr()
norm.data = np.log1p(norm.data)
rng = np.random.default_rng(1701)
train = np.concatenate([rng.choice(np.flatnonzero(obs.donor.eq(d)),
                                   min(1500, int(obs.donor.eq(d).sum())), replace=False)
                        for d in sorted(obs.donor.unique())])
train.sort()
z = norm[train]
mean = np.asarray(z.mean(axis=0)).ravel()
var = np.asarray(z.power(2).mean(axis=0)).ravel()-mean**2
detected = np.asarray((z > 0).mean(axis=0)).ravel()
excluded = ann.SYMBOL.isin(panel.mouse_symbol_primary).to_numpy() | ann.SYMBOL.str.match(r'(?i)^(mt-|rpl|rps)').to_numpy()
valid = (~excluded) & (detected >= .01) & (detected <= .95) & (var > 0)
v = pd.DataFrame({'index': np.flatnonzero(valid), 'mean': mean[valid],
                  'variance': var[valid], 'detected_fraction': detected[valid]})
# Rank variance within mean-expression bins, then choose the largest ranks.
v['mean_bin'] = pd.qcut(v['mean'], 20, labels=False, duplicates='drop')
v['dispersion'] = v.variance/np.maximum(v['mean'], 1e-8)
v['within_bin_rank'] = v.groupby('mean_bin').dispersion.rank(pct=True)
hv = v.sort_values(['within_bin_rank', 'dispersion'], ascending=False).head(2000)
features = hv['index'].to_numpy()
assert len(features) == 2000 and not excluded[features].any()
hv = hv.merge(ann.reset_index(), left_on='index', right_on='index')
hv.to_csv(B/'tables/state_nonpanel_variable_features.csv', index=False)
pca = PCA(n_components=30, svd_solver='randomized', random_state=1701)
pca.fit(norm[train][:, features].toarray())
pcs = pca.transform(norm[:, features].toarray()).astype(np.float32)
if depth_adjusted:
    # Exploratory technical sensitivity, prompted by primary PCA/depth audit.
    # Fit the correction on the same balanced training set without exposure.
    technical = np.column_stack([np.ones(len(obs)), np.log10(obs.library_UMI), obs.mt_fraction])
    coeff = np.linalg.lstsq(technical[train], pcs[train], rcond=None)[0]
    pcs = (pcs-technical@coeff).astype(np.float32)
    pd.DataFrame(coeff,index=['intercept','log10_UMI','mt_fraction'],columns=['PC'+str(i+1) for i in range(30)]).to_csv(B/'tables/PCA_technical_regression_coefficients.csv')
np.savez_compressed(B/'prepared/state_PCA_model.npz', components=pca.components_,
                    mean=pca.mean_, selected_feature_indices=features,
                    train_cell_indices=train, explained_variance_ratio=pca.explained_variance_ratio_)
for j in range(3):
    obs['PC'+str(j+1)] = pcs[:, j]
qc, markers, donor_rows, stability = [], [], [], []
context = ['P2ry12', 'Tmem119', 'Sall1', 'Hexb', 'Csf1r', 'Mrc1', 'Pf4',
           'Lyve1', 'Ccr2', 'Lyz2', 'Cd74', 'H2-Ab1', 'Apoe', 'Lpl',
           'Trem2', 'Fos', 'Jun', 'Aif1', 'Pycard', 'Aim2', 'Ninj1']
ann.to_csv(B/'prepared/state_annotation.csv', index=False)
for k in [4, 6]:
    model = MiniBatchKMeans(n_clusters=k, n_init=10, batch_size=2048,
                           random_state=1701, max_iter=200)
    model.fit(pcs[train])
    lab = model.predict(pcs)
    # Numerical naming follows sorted centroid PC1; biological labels are withheld.
    order = np.argsort(model.cluster_centers_[:, 0])
    remap = np.empty(k, dtype=int); remap[order] = np.arange(k)
    lab = remap[lab]
    obs['state_K'+str(k)] = lab+1
    other = MiniBatchKMeans(n_clusters=k, n_init=10, batch_size=2048,
                           random_state=2901, max_iter=200).fit(pcs[train])
    lab2 = other.predict(pcs)
    contingency = pd.crosstab(lab, lab2).reindex(index=range(k), columns=range(k), fill_value=0).to_numpy()
    ii, jj = linear_sum_assignment(-contingency)
    stability.append({'K': k, 'seed1': 1701, 'seed2': 2901,
                      'adjusted_rand_index': adjusted_rand_score(lab, lab2),
                      'matched_assignment_fraction': contingency[ii,jj].sum()/len(lab),
                      'interpretation': 'Algorithm-seed sensitivity in one shared PCA; not independent biological validation'})
    pb = {}
    for state in range(1, k+1):
        for s in md.itertuples():
            mask = obs['sample'].eq(s.sample).to_numpy() & (lab == state-1)
            key = s.sample+'_S'+str(state)
            pb[key] = np.asarray(x[mask].sum(axis=0)).ravel()
            row = {'K': k, 'state': state, 'key': key, 'sample': s.sample,
                   'donor': s.donor, 'flight': s.flight, 'age': s.age,
                   'n_cells': int(mask.sum()), 'donor_total_cells': int(obs['sample'].eq(s.sample).sum()),
                   'cell_fraction': float(mask.sum()/obs['sample'].eq(s.sample).sum()),
                   'raw_UMIs': int(obs.loc[mask, 'library_UMI'].sum()),
                   'median_UMI': float(obs.loc[mask, 'library_UMI'].median()),
                   'median_features': float(obs.loc[mask, 'n_features'].median()),
                   'median_mt_fraction': float(obs.loc[mask, 'mt_fraction'].median())}
            qc.append(row)
            avg = np.asarray(x[mask].astype(float).multiply((1e4/obs.loc[mask, 'library_UMI'].to_numpy())[:,None]).mean(axis=0)).ravel() if mask.any() else np.full(x.shape[1],np.nan)
            for gene in context:
                idx = np.flatnonzero(ann.SYMBOL.eq(gene))
                if len(idx) == 1:
                    markers.append({'K': k, 'state': state, 'donor': s.donor,
                                    'gene': gene, 'mean_cell_CP10K': avg[idx[0]]})
            for p in panel.itertuples():
                idx = np.flatnonzero(ann.SYMBOL.eq(p.mouse_symbol_primary))
                donor_rows.append({'K': k, 'state': state, 'sample': s.sample,
                                   'donor': s.donor, 'flight': s.flight, 'age': s.age,
                                   'analyte': p.analyte, 'mouse_symbol': p.mouse_symbol_primary,
                                   'mean_cell_CP10K': avg[idx[0]] if len(idx) == 1 else np.nan})
    pd.DataFrame(pb, index=ann.ENSEMBL).to_csv(B/'prepared'/f'state_K{k}_pseudobulk.csv.gz')
obs.to_csv(B/'prepared/forebrain_state_cell_registry.csv.gz', index=False)
pd.DataFrame(qc).to_csv(B/'tables/state_mouse_QC.csv', index=False)
pd.DataFrame(markers).to_csv(B/'tables/state_context_markers.csv', index=False)
pd.DataFrame(donor_rows).to_csv(B/'tables/state_target_mouse_expression.csv', index=False)
pd.DataFrame(stability).to_csv(B/'tables/state_algorithm_seed_stability.csv', index=False)
summary = {'dataset': 'GSE295428 / RR8 myeloid-enriched forebrain',
           'n_mice': int(obs.donor.nunique()), 'n_cells': len(obs),
           'training_cells': len(train), 'maximum_training_cells_per_mouse': 1500,
           'n_nonpanel_variable_features': len(features), 'PCA_components': 30,
           'PCA_explained_variance_fraction': float(pca.explained_variance_ratio_.sum()),
           'K_values': [4, 6], 'selection': 'Non-panel, non-mitochondrial, non-ribosomal genes; mean-bin dispersion ranking of balanced training cells; no exposure or age label used',
           'QC': 'Same primary 200-2500 feature and <=50% mitochondrial gate as completed analysis',
           'inference': 'Post-result exploratory; eight mice, not cells, are exposure replicates; state learning uses the same eight mice',
           'remaining_confound': 'Flight and ground dissections occurred on different days; group and day are perfectly confounded',
           'technical_sensitivity': 'Regressed log10 cell UMI and mitochondrial fraction from PCs using balanced training cells, no exposure labels. Post-result sensitivity; may also remove biological RNA-content variation.' if depth_adjusted else 'Primary normalized-expression PCA without technical regression',
           'raw_input_hashes': hashes}
(B/'tables/state_preparation_summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='raw_input_hashes'}, indent=2), flush=True)
