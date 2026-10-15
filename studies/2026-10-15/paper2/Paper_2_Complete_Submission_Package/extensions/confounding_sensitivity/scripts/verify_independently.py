"""Independent direct-feature/pairwise checks; imports no analysis module."""
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
B=Path(__file__).resolve().parents[1];P=B.parents[1];R=B/'results';checks=[]
def check(name,ok,detail=''):
    checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name+': '+str(detail))
z=np.load(P/'extensions/multiway_splicing/tables/frozen_features.npz');samples=z['samples'];H=z['H'];meta=pd.read_csv(P/'inputs/animal_metadata.csv').set_index('sample').loc[samples];g=meta.group.values;groups=['FLT','GC','VIV'];ix=[np.flatnonzero(g==a) for a in groups];W=np.load(R/'bootstrap_weights.npz')['weights'];geometry=pd.read_csv(R/'comparator_geometry.csv');tip=pd.read_csv(R/'bias_tipping_points.csv');q=pd.read_csv(R/'QC_adjustment_sensitivity.csv');qm=pd.read_csv(R/'animal_qc_reconstructed.csv').set_index('sample').loc[samples]
check('unique biological animals',meta['Source Name'].nunique()==26)
check('group sizes',list(map(len,ix))==[9,9,8])
check('flight dissection-date alias',np.array_equal(g=='FLT',meta['Comment[BSP Dissection Date]'].values=='14-Jan-2021'))
check('same-date ground comparison',meta.loc[g!='FLT','Comment[BSP Dissection Date]'].eq('17-Jan-2021').all())
A=np.column_stack([np.ones(26),g=='FLT',g=='GC',meta['Comment[Euthanasia Date]'].eq('14-Jan-2021')]);check('rank exactly three',np.linalg.matrix_rank(A)==3);check('nullspace exact zero',np.array_equal(A@np.array([0,1,0,-1]),np.zeros(26)))
fit=pd.read_csv(R/'identical_fitted_values.csv');check('distinct coefficient sets identical fits',np.array_equal(fit.fit_1,fit.fit_2))
check('bootstrap counts',W.shape==(10000,26) and all(np.all(W[:,i].sum(1)==len(i)) for i in ix))
# Reconstruct H directly from raw per-alternative animal counts, which are frozen
# alongside the package-wide matrix. Every selected junction's counts are then
# independently compared with the canonical all-animal count matrix by chunks.
a=pd.read_csv(P/'extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz');lookup=a[['cluster_id','alternative_index','junction_id']].drop_duplicates();order=[(cid,j) for cid,k in zip(z['cluster_ids'],np.diff(z['offsets'])) for j in range(k)];a2=a.pivot(index=['cluster_id','alternative_index'],columns='sample',values='count').loc[order,samples].to_numpy();tot=np.add.reduceat(a2,z['offsets'][:-1],axis=0);h2=np.sqrt(a2/np.repeat(tot,np.diff(z['offsets']),axis=0)).T
check('all multiway features reconstructed from counts',np.max(abs(h2-H))<2e-15)
mapcount=a.pivot_table(index='junction_id',columns='sample',values='count',aggfunc='first').loc[:,samples];found=[]
for c in pd.read_csv(P/'inputs/all_animal_junction_counts.csv.gz',index_col=0,chunksize=50000):
    use=c.index.intersection(mapcount.index)
    if len(use):check('canonical raw count chunk '+str(len(found)),np.array_equal(c.loc[use,samples],mapcount.loc[use,samples]));found.extend(use)
check('all selected distinct junctions verified',len(set(found))==3248)
f=pd.read_csv(P/'results/all26_animal_fractions.csv.gz',index_col=0);gate=pd.read_csv(P/'results/measurement_gates.csv.gz');f=f.loc[gate.loc[gate.measured_primary_all26,'event_id'],samples];frames={'multiway_Hellinger':h2/np.sqrt(1846),'binary_RMS_pp':100*f.to_numpy().T/np.sqrt(5249)}
for frame,Y in frames.items():
    K=pd.read_csv(R/f'{frame}_animal_kernel.csv',index_col=0).to_numpy();check(frame+' kernel from raw features',np.allclose(Y@Y.T,K,rtol=1e-11,atol=1e-11))
    G=np.array([[np.mean(Y[ii],axis=0)@np.mean(Y[jj],axis=0) for jj in ix] for ii in ix]);Gu=G.copy()
    for k,i in enumerate(ix):
        # Direct averaging over distinct ordered pairs, independent of trace correction.
        v=[Y[u]@Y[v] for u in i for v in i if u!=v];Gu[k,k]=np.mean(v)
    check(frame+' Ustat off-diagonal pairs',np.allclose(Gu,pd.read_csv(R/f'{frame}_Ustat_gram.csv',index_col=0),rtol=1e-10,atol=1e-10))
    mm=np.array([Y[i].mean(0) for i in ix]);d=[mm[0]-mm[1],mm[0]-mm[2],mm[1]-mm[2]];names=['FLT_vs_GC','FLT_vs_VIV','GC_vs_VIV']
    for j,c in enumerate(names[:2]):
        tt=tip[(tip.frame==frame)&(tip.contrast==c)].iloc[0];norm=np.linalg.norm(d[j]);base=np.linalg.norm(d[2]);check(frame+c+' tipping direct norm',abs(norm/base-tt.plugin_tipping_lambda)<1e-10)
        for lam in [0,.5,1,2]:
            bound=lam*base;b=d[j]*min(bound/norm,1);check(frame+c+' attainable lower '+str(lam),abs(np.linalg.norm(d[j]-b)-max(0,norm-bound))<1e-11)
            b=-d[j]*bound/norm;check(frame+c+' attainable upper '+str(lam),abs(np.linalg.norm(d[j]-b)-(norm+bound))<1e-11)
    XX=np.column_stack([np.ones(26),g=='FLT',g=='GC',qm.rin_z,qm.log2_read_depth_z]);coeff=np.linalg.lstsq(XX,Y,rcond=None)[0];adjust=[coeff[1]-coeff[2],coeff[1],coeff[2]]
    for j,c in enumerate(names):
        v=q[(q.frame==frame)&(q.model=='group_RIN_log2_depth')&(q.contrast==c)].iloc[0];check(frame+c+' direct multivariate OLS',abs(np.linalg.norm(adjust[j])-v.adjusted_norm)<1e-10)
    boot=pd.read_csv(R/'comparator_geometry_bootstrap.csv.gz');boot=boot[(boot.frame==frame)&(boot.estimator=='Ustat_noise_corrected')].set_index('replicate')
    for r in [0,1,17,991,9999]:
        yy=[Y[np.repeat(i,W[r,i])] for i in ix];bm=np.array([a.mean(0) for a in yy]);bu=bm@bm.T
        for j,v in enumerate(yy):bu[j,j]=sum(v[a]@v[b] for a in range(len(v)) for b in range(len(v)) if a!=b)/(len(v)*(len(v)-1))
        gc=Gu+bu-G
        expected=gc[0,0]+gc[1,1]-2*gc[0,1];check(frame+' centered Ustat bootstrap '+str(r),abs(expected-boot.loc[r,'FG_squared'])<1e-10)
    for metric in ['FG_squared','FV_squared','GV_squared','FV_dot_GV']:
        row=geometry[(geometry.frame==frame)&(geometry.estimator=='Ustat_noise_corrected')&(geometry.metric==metric)].iloc[0];err=boot[metric]-row.estimate;lo,hi=row.estimate-np.quantile(err,[.975,.025]);check(frame+metric+' basic intervals',np.allclose([lo,hi],[row.interval_low,row.interval_high]))
    er=np.load(R/f'{frame}_centered_errors.npz')['centered_error_norm'];radius=np.quantile(er.max(1),.95);check(frame+' joint radius',np.allclose(tip.loc[tip.frame==frame,'joint_error_radius95'],radius))
    # Verify complete noise-corrected ratio invalid handling, not silent clipping.
    invalid=~boot.valid;check(frame+' invalid ratios missing',boot.loc[invalid,['comparator_coordinate','cosine_FV_GV','orthogonal_fraction_squared','orthogonal_distance']].isna().all().all())
for row in pd.read_csv(R/'input_manifest.csv').itertuples():check('input hash '+row.package_relative_path,hashlib.sha256((P/row.package_relative_path).read_bytes()).hexdigest()==row.sha256)
s=json.loads((R/'summary.json').read_text());check('frozen plan hash',hashlib.sha256((B/'RETROSPECTIVE_PLAN_2026-10-15.md').read_bytes()).hexdigest()==s['plan_sha256']);check('transparent amendment hash',hashlib.sha256((B/'AMENDMENT_USTAT_UNCERTAINTY_2026-10-15.md').read_bytes()).hexdigest()==s['amendment_sha256'])
pd.DataFrame(checks).to_csv(R/'independent_verification.csv',index=False);(R/'independent_verification.json').write_text(json.dumps(dict(status='passed',n_checks=len(checks),failed=0,checks=checks),indent=2));print('Passed',len(checks),'independent checks')
