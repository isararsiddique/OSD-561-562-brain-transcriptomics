"""Frozen identification geometry and explicit, noncausal bias-bound analysis."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
from pathlib import Path
import numpy as np
import pandas as pd
import hashlib,json
BASE=Path(__file__).resolve().parents[1]; P=BASE.parents[1]; R=BASE/'results'
SEED=2026101531; B=10000; GROUPS=['FLT','GC','VIV']; CONTRASTS=['FLT_vs_GC','FLT_vs_VIV','GC_vs_VIV']
L=np.array([[1,-1,0],[1,0,-1],[0,1,-1.]])

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gramnorm(a,K,b=None):return np.einsum('...i,ij,...j->...',a,K,a if b is None else b)
def metrics(G):
    # Gram matrix in order F,G,V; d=F-V, c=G-V.
    C=np.einsum('ai,...ij,bj->...ab',L,G,L)
    df2=C[...,1,1];c2=C[...,2,2];dot=C[...,1,2]
    valid=(df2>0)&(c2>0)&(C[...,0,0]>0)
    t=np.divide(dot,c2,out=np.full_like(dot,np.nan),where=valid)
    residual=df2-np.divide(dot*dot,c2,out=np.full_like(dot,np.nan),where=valid)
    valid=valid&(residual>=0)
    return dict(FG_squared=C[...,0,0],FV_squared=df2,GV_squared=c2,FV_dot_GV=dot,comparator_coordinate=np.where(valid,t,np.nan),cosine_FV_GV=np.where(valid,dot/np.sqrt(np.maximum(df2*c2,1e-300)),np.nan),orthogonal_fraction_squared=np.where(valid,residual/df2,np.nan),orthogonal_distance=np.where(valid,np.sqrt(np.maximum(residual,0)),np.nan),valid=valid)

def main():
    assert (BASE/'RETROSPECTIVE_PLAN_2026-10-15.md').exists()
    z=np.load(P/'extensions/multiway_splicing/tables/frozen_features.npz'); H=z['H'];samples=z['samples'];E=len(z['cluster_ids'])
    m=pd.read_csv(P/'inputs/animal_metadata.csv').set_index('sample').loc[samples]; groups=m.group.to_numpy();ix=[np.flatnonzero(groups==g) for g in GROUPS];n=np.array([len(i) for i in ix]);M=np.array([(groups==g)/len(i) for g,i in zip(GROUPS,ix)])
    design=pd.DataFrame({'sample':samples,'source_animal':m['Source Name'].values,'group':groups,'euthanasia_date':m['Comment[Euthanasia Date]'].values,'dissection_date':m['Comment[BSP Dissection Date]'].values})
    design['intercept']=1.;design['flight_indicator']=(groups=='FLT').astype(float);design['ground_indicator']=(groups=='GC').astype(float);design['early_date_indicator']=(design.euthanasia_date=='14-Jan-2021').astype(float)
    design.to_csv(R/'animal_design.csv',index=False)
    A=design[['intercept','flight_indicator','ground_indicator','early_date_indicator']].to_numpy();null=np.array([0,1,0,-1]);rank=int(np.linalg.matrix_rank(A));assert rank==3 and np.max(abs(A@null))==0
    beta=np.array([1,2,3,4.]);beta2=beta+7*null
    pd.DataFrame({'column':['intercept','flight_indicator','ground_indicator','early_date_indicator'],'null_vector':null,'coefficient_set_1':beta,'coefficient_set_2':beta2}).to_csv(R/'identification_nullspace.csv',index=False)
    pd.DataFrame({'sample':samples,'fit_1':A@beta,'fit_2':A@beta2}).to_csv(R/'identical_fitted_values.csv',index=False)
    idsum=dict(n_animals=len(m),n_columns=A.shape[1],rank=rank,nullity=A.shape[1]-rank,singular_values=np.linalg.svd(A,compute_uv=False).tolist(),flight_date_indicator_identical=bool(np.array_equal(A[:,1],A[:,3])),null_fit_max_difference=float(np.max(abs(A@beta-A@beta2))))
    qc=pd.read_csv(P/'inputs/library_qc.csv');qc['sample']=qc['sample'].str.replace(r'_techrep[0-9]+$','_techrep1',regex=True)
    qc=qc.groupby('sample').agg(read_depth=('read_depth','sum'),rin=('rin','mean')).loc[samples];old=pd.read_csv(P/'results/animal_qc.csv').set_index('sample').loc[samples]
    assert np.allclose(qc[['read_depth','rin']],old[['read_depth','rin']]) and qc.notna().all().all()
    qc['log2_read_depth']=np.log2(qc.read_depth);qc['group']=groups
    for c in ['rin','log2_read_depth']:qc[c+'_z']=(qc[c]-qc[c].mean())/qc[c].std(ddof=0)
    qc.to_csv(R/'animal_qc_reconstructed.csv');qc.groupby('group')[['rin','read_depth','log2_read_depth']].agg(['mean','std','min','max']).to_csv(R/'qc_by_group.csv')
    gates=pd.read_csv(P/'results/measurement_gates.csv.gz');ff=pd.read_csv(P/'results/all26_animal_fractions.csv.gz',index_col=0).loc[gates.loc[gates.measured_primary_all26,'event_id'],samples]
    assert ff.shape==(5249,26)
    frames={'multiway_Hellinger':H/np.sqrt(2*E),'binary_RMS_pp':ff.to_numpy().T*100/np.sqrt(len(ff))}
    rng=np.random.Generator(np.random.PCG64(SEED));W=np.zeros((B,26),dtype=np.int16)
    for i in ix:W[:,i]=rng.multinomial(len(i),np.full(len(i),1/len(i)),size=B)
    np.savez_compressed(R/'bootstrap_weights.npz',weights=W,samples=samples,groups=groups,seed=SEED)
    WM=np.stack([W*(groups==g)[None,:]/len(i) for g,i in zip(GROUPS,ix)],axis=1)
    contrasts=L@M;Db=np.einsum('ag,bgi->bai',L,WM)
    allgeo=[];allboot=[];allsens=[];alltip=[];allqc=[];allqboot=[];allcoef=[];summaries={}
    Xgroup=np.column_stack([np.ones(26),groups=='FLT',groups=='GC']);Z=qc[['rin_z','log2_read_depth_z']].to_numpy();modelcols={'group_only':[],'group_RIN':[0],'group_log2_depth':[1],'group_RIN_log2_depth':[0,1]}
    modeldesign={name:np.column_stack([Xgroup,Z[:,cols]]) for name,cols in modelcols.items()};models={}
    for name,XX in modeldesign.items():
        pinv=np.linalg.pinv(XX);LC=np.zeros((3,XX.shape[1]));LC[:,:3]=[[0,1,-1],[0,1,0],[0,0,1]];weights=LC@pinv;models[name]=(weights,LC)
        for j,row in enumerate(weights):
            for sample,v in zip(samples,row):allcoef.append(dict(model=name,contrast=CONTRASTS[j],sample=sample,linear_coefficient=v))
    pd.DataFrame(allcoef).to_csv(R/'QC_contrast_linear_weights.csv',index=False)
    for frame,Y in frames.items():
        K=Y@Y.T;pd.DataFrame(K,index=samples,columns=samples).to_csv(R/f'{frame}_animal_kernel.csv')
        G=M@K@M.T;noise=np.array([np.trace(K[np.ix_(i,i)])/len(i)-G[j,j] for j,i in enumerate(ix)])/(n-1);Gu=G-np.diag(noise)
        Gb=np.einsum('bgi,ij,bhj->bgh',WM,K,WM);noise_b=np.stack([(W[:,i]@np.diag(K)[i]/len(i)-Gb[:,j,j])/(len(i)-1) for j,i in enumerate(ix)],axis=1);Gub=Gb.copy();Gub[:,np.arange(3),np.arange(3)]-=noise_b
        pd.DataFrame(G,index=GROUPS,columns=GROUPS).to_csv(R/f'{frame}_centroid_gram.csv');pd.DataFrame(Gu,index=GROUPS,columns=GROUPS).to_csv(R/f'{frame}_Ustat_gram.csv')
        geomsummary={}
        for kind,gg,bg in [('plugin',G,Gb),('Ustat_noise_corrected',Gu,Gu[None,:,:]+Gub-G[None,:,:])]:
            mg=metrics(gg);mb=metrics(bg);validfrac=float(np.mean(mb['valid']));geomsummary[kind]={'bootstrap_valid_fraction':validfrac}
            for metric,value in mg.items():
                if metric=='valid':continue
                vals=mb[metric]
                if kind=='plugin':
                    lo,hi=np.nanquantile(vals,[.025,.975]);method='descriptive_resampling_quantiles_not_confidence_interval'
                elif metric in ['FG_squared','FV_squared','GV_squared','FV_dot_GV']:
                    qlo,qhi=np.quantile(vals-float(value),[.025,.975]);lo,hi=float(value)-qhi,float(value)-qlo;method='approximate_basic_analytically_centered_Ustat'
                else:
                    lo=hi=np.nan;method='not_reported_unstable_nonlinear_Ustat_ratio'
                allgeo.append(dict(frame=frame,estimator=kind,metric=metric,estimate=float(value),interval_low=float(lo),interval_high=float(hi),interval_method=method,bootstrap_valid_fraction=validfrac))
                geomsummary[kind][metric]=float(value)
            df=pd.DataFrame({k:v for k,v in mb.items()});df['frame']=frame;df['estimator']=kind;df['replicate']=np.arange(B);allboot.append(df)
        d=np.sqrt(np.maximum(np.diag(L@G@L.T),0));err=np.sqrt(np.maximum(gramnorm(Db-contrasts[None,:,:],K),0));radius=float(np.quantile(err.max(axis=1),.95));benchmark=d[2]
        for ci in range(2):
            alltip.append(dict(frame=frame,contrast=CONTRASTS[ci],observed_norm=d[ci],control_benchmark_norm=benchmark,plugin_tipping_lambda=d[ci]/benchmark,joint_error_radius95=radius,conservative_tipping_lambda=max(0,d[ci]-radius)/(benchmark+radius),causal_effect_identified=False))
            for lam in np.linspace(0,2,81):allsens.append(dict(frame=frame,contrast=CONTRASTS[ci],lambda_benchmark=lam,bias_bound_B=lam*benchmark,plugin_norm_lower=max(0,d[ci]-lam*benchmark),plugin_norm_upper=d[ci]+lam*benchmark,approx_joint_sampling_norm_lower=max(0,d[ci]-radius-lam*(benchmark+radius)),approx_joint_sampling_norm_upper=d[ci]+radius+lam*(benchmark+radius),benchmark_is_date_estimate=False))
        qsummary={}
        for name,XX in modeldesign.items():
            weights,LC=models[name];dd=np.sqrt(np.maximum(gramnorm(weights,K),0));dot=gramnorm(weights,K,contrasts);cos=dot/(dd*d);rankq=int(np.linalg.matrix_rank(XX));lev=np.diag(XX@np.linalg.pinv(XX));bweights=[];valid=[]
            for bb,w in enumerate(W):
                xwx=XX.T@(w[:,None]*XX)
                if np.linalg.matrix_rank(xwx)<XX.shape[1]:continue
                wb=LC@np.linalg.solve(xwx,XX.T*w[None,:]);bweights.append(wb);valid.append(bb)
            bw=np.stack(bweights);bnorm=np.sqrt(np.maximum(gramnorm(bw,K),0));berr=np.sqrt(np.maximum(gramnorm(bw-weights[None,:,:],K),0));rqc=float(np.quantile(berr.max(axis=1),.95))
            originalb=np.sqrt(np.maximum(gramnorm(Db[np.array(valid)],K),0));ratio=bnorm/originalb
            for ci,c in enumerate(CONTRASTS):
                lo,hi=np.quantile(ratio[:,ci],[.025,.975]);allqc.append(dict(frame=frame,model=name,contrast=c,adjusted_norm=dd[ci],unadjusted_norm=d[ci],norm_ratio=dd[ci]/d[ci],cosine_with_unadjusted=cos[ci],ratio_bootstrap_low=lo,ratio_bootstrap_high=hi,joint_error_radius95=rqc,approx_norm_low=max(0,dd[ci]-rqc),approx_norm_high=dd[ci]+rqc,design_rank=rankq,n_parameters=XX.shape[1],residual_df=26-rankq,max_leverage=float(lev.max()),bootstrap_valid=len(valid),bootstrap_invalid=B-len(valid)))
                allqboot.append(pd.DataFrame(dict(frame=frame,model=name,contrast=c,replicate=valid,norm=bnorm[:,ci],norm_ratio=ratio[:,ci],centered_error_norm=berr[:,ci])))
            qsummary[name]=dict(valid_bootstrap=len(valid),rank=rankq,parameters=XX.shape[1],max_leverage=float(lev.max()))
        summaries[frame]=dict(n_features=Y.shape[1],n_units=E if frame=='multiway_Hellinger' else len(ff),joint_error_radius95=radius,geometry=geomsummary,QC_models=qsummary)
        np.savez_compressed(R/f'{frame}_centered_errors.npz',centered_error_norm=err,contrasts=CONTRASTS)
    pd.DataFrame(allgeo).to_csv(R/'comparator_geometry.csv',index=False);pd.concat(allboot).to_csv(R/'comparator_geometry_bootstrap.csv.gz',index=False);pd.DataFrame(allsens).to_csv(R/'bias_bound_sensitivity_curves.csv',index=False);pd.DataFrame(alltip).to_csv(R/'bias_tipping_points.csv',index=False);pd.DataFrame(allqc).to_csv(R/'QC_adjustment_sensitivity.csv',index=False);pd.concat(allqboot).to_csv(R/'QC_bootstrap_complete.csv.gz',index=False)
    files=['inputs/all_animal_junction_counts.csv.gz','extensions/multiway_splicing/tables/animal_alternative_compositions.csv.gz','inputs/animal_metadata.csv','inputs/library_qc.csv','results/animal_qc.csv','results/all26_animal_fractions.csv.gz','results/measurement_gates.csv.gz','extensions/multiway_splicing/tables/frozen_features.npz']
    pd.DataFrame([dict(package_relative_path=f,sha256=sha(P/f),bytes=(P/f).stat().st_size) for f in files]).to_csv(R/'input_manifest.csv',index=False)
    summary=dict(seed=SEED,bootstrap_replicates=B,identification=idsum,frames=summaries,plan_sha256=sha(BASE/'RETROSPECTIVE_PLAN_2026-10-15.md'),amendment_sha256=sha(BASE/'AMENDMENT_USTAT_UNCERTAINTY_2026-10-15.md'),interpretation='All results are observational. GC-VIV is a hypothetical scale benchmark, not an estimate of date bias. No causal rescue.')
    (R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
