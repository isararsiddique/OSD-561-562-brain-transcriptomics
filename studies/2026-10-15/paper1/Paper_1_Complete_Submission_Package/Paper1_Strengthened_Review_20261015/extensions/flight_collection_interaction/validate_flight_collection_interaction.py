#!/usr/bin/env python3
"""Independent numeric/model validation of the finished interaction extension."""
import json
import numpy as np
import pandas as pd
from scipy import stats
from PIL import Image
import run_flight_collection_interaction as r


def main():
    r.check_input_integrity()
    panels=pd.read_csv(r.OUT/'panel_flight_collection_interactions.csv')
    genes=pd.read_csv(r.OUT/'gene_flight_collection_interactions_exploratory.csv.gz')
    assert len(panels)==32 and panels.groupby(['accession','sensitivity']).size().eq(8).all()
    assert np.isfinite(panels.p_wild_approximate).all()
    assert np.allclose(r.bh(panels.p_wild_approximate),panels.q_BH_complete_panel_family)
    assert np.allclose(r.bh(genes.pvalue),genes.q_BH_complete_gene_interaction_family,equal_nan=True)
    finite=genes.pvalue.notna()
    assert np.allclose(2*stats.norm.sf(abs(genes.loc[finite].log2FoldChange/genes.loc[finite].lfcSE)),genes.loc[finite].pvalue)
    assert np.allclose(genes.log2FoldChange,genes.FLT_minus_GC_ISST_log2-genes.FLT_minus_GC_LAR_log2,equal_nan=True)
    summary=json.loads((r.OUT/'analysis_summary.json').read_text())
    conflicts=set(summary['conflicting_animals_either_region'])
    ann=pd.read_csv(r.BASE/'inputs/gene_annotation.csv.gz')
    numeric=[];nb_audits=[]
    for number in ['561','562']:
        d=r.load_dataset(number)
        assert d['meta'].animal.is_unique and d['counts'].columns.equals(d['meta'].index)
        assert np.allclose(d['counts'],np.round(d['counts']))
        for sensitivity in ['all_animals','exclude_any_region_metadata_conflict']:
            m=d['meta'] if sensitivity=='all_animals' else d['meta'][~d['meta'].animal.isin(conflicts)]
            y=d['y'].loc[:,m.index].to_numpy();x=r.design(m).to_numpy()
            groups=np.column_stack([(m.flight.eq(f)&m.cohort.eq(c)).to_numpy(float)
                                    for f,c in [('GC','LAR'),('FLT','LAR'),('GC','ISS_T'),('FLT','ISS_T')]])
            alternate=np.column_stack([groups,(m.age=='29').to_numpy(float)])
            alt=np.linalg.lstsq(alternate,y.T,rcond=None)[0]
            direct=np.linalg.lstsq(x,y.T,rcond=None)[0][-1]
            check=alt[3]-alt[2]-alt[1]+alt[0]
            assert np.allclose(direct,check,atol=1e-10)
            pidx=r.indices(d['y'].index,ann)
            yp=np.vstack([y[i].mean(0) for i in pidx.values()])
            b,se,t,coef,res,h=r.hc3(yp,x)
            inv=np.linalg.inv(x.T@x)
            for i in range(len(yp)):
                covariance=inv@(x.T@np.diag((res[i]/(1-h))**2)@x)@inv
                assert np.isclose(se[i],np.sqrt(covariance[-1,-1]))
            numeric.append({'accession':f'OSD-{number}','sensitivity':sensitivity,'n_animals':len(m),
                            'rank':int(np.linalg.matrix_rank(x)),
                            'independent_four_group_contrast_max_error':float(abs(direct-check).max()),
                            'HC3_direct_sandwich_checked_all_8_panels':True})
        accession=f'OSD-{number}';label='OSD561_cerebellum' if number=='561' else 'OSD562_hippocampus'
        old=pd.read_csv(r.BASE/f'results/count_models/{accession}_gene_results.csv.gz')
        iss=old[old.contrast==label+'_FLT_vs_GC_within_ISS_T'].set_index('ensembl_id')
        lar=old[old.contrast==label+'_FLT_vs_GC_within_LAR'].set_index('ensembl_id')
        expected=(iss.log2FoldChange-lar.log2FoldChange).rename('cached_four_group_interaction')
        current=genes[(genes.accession==accession)&(genes.sensitivity=='all_animals')].set_index('ensembl_id').join(expected)
        current['difference_from_cached_parameterization']=current.log2FoldChange-current.cached_four_group_interaction
        finite=current.pvalue.notna()&current.cached_four_group_interaction.notna()
        discovered=finite & current.q_BH_complete_gene_interaction_family.le(.05)
        max_discovery_error=float(current.loc[discovered].difference_from_cached_parameterization.abs().max())
        assert max_discovery_error<.01
        current.loc[current.difference_from_cached_parameterization.abs()>.01].reset_index().to_csv(
            r.OUT/f'{accession}_NB_parameterization_discrepancies.csv',index=False)
        nb_audits.append({'accession':accession,'n_finite_p_compared':int(finite.sum()),
                          'median_absolute_difference_finite_p':float(current.loc[finite].difference_from_cached_parameterization.abs().median()),
                          'max_absolute_difference_finite_p':float(current.loc[finite].difference_from_cached_parameterization.abs().max()),
                          'n_finite_p_differences_above_0_01':int((current.loc[finite].difference_from_cached_parameterization.abs()>.01).sum()),
                          'primary_global_FDR_discovery_max_difference':max_discovery_error,
                          'note':'Large coefficient differences at unstable/non-significant GLM genes are disclosed; no primary FDR discovery has a discrepancy >0.01 log2.'})
    # Independently reproduce the highlighted panel's 499,999 sign-weight draws.
    d=r.load_dataset('562');m=d['meta'];x=r.design(m).to_numpy()
    idx=r.indices(d['y'].index,ann)['Circadian'];y=d['y'].iloc[idx].mean(axis=0).to_numpy()[None,:]
    observed=r.hc3(y,x)[2][0]
    x0=x[:,:-1];fit=np.linalg.lstsq(x0,y.T,rcond=None)[0].T@x0.T
    h0=np.einsum('ij,ji->i',x0,np.linalg.pinv(x0));res0=(y-fit)/np.sqrt(1-h0)
    inv=np.linalg.inv(x.T@x);mapping=inv@x.T;h=np.einsum('ij,ji->i',x,mapping)
    weight=(x@inv[:,-1])/(1-h)
    row=panels[(panels.accession=='OSD-562')&(panels.sensitivity=='all_animals')&(panels.panel=='Circadian')].iloc[0]
    draws=int(row.n_wild_draws);rng=np.random.default_rng(r.SEED+562);exceeded=0
    for start in range(0,draws,5000):
        n=min(5000,draws-start);signs=rng.choice([-1.,1.],size=(n,len(m)))
        simulated=fit[None,:,:]+res0[None,:,:]*signs[:,None,:]
        coef=simulated@mapping.T;res=simulated-coef@x.T
        se=np.sqrt(res**2@(weight**2));ts=coef[:,:,-1]/se
        if start==0:
            for i in range(100):assert np.allclose(ts[i],r.hc3(simulated[i],x)[2],atol=1e-10)
        exceeded+=int((abs(ts)>=abs(observed)-1e-12).sum())
    recomputed_p=(exceeded+1)/(draws+1)
    assert np.isclose(recomputed_p,row.p_wild_approximate)
    loo=pd.read_csv(r.OUT/'leave_one_animal_out_panel_interactions.csv')
    log=pd.read_csv(r.OUT/'leave_one_gene_out_panel_interactions.csv')
    stability={}
    for sensitivity in ['all_animals','exclude_any_region_metadata_conflict']:
        key=(loo.accession=='OSD-562')&(loo.sensitivity==sensitivity)&(loo.panel=='Circadian')
        a=loo.loc[key].interaction_without_animal
        key=(log.accession=='OSD-562')&(log.sensitivity==sensitivity)&(log.panel=='Circadian')
        g=log.loc[key].interaction_without_gene_log2
        stability[sensitivity]={'animal_omission_effect_range':[float(a.min()),float(a.max())],
                                'gene_omission_effect_range':[float(g.min()),float(g.max())],
                                'interpretation':'Effect-only omissions, no additional significance claims or FDR family'}
    figure_records=[]
    for name in ['figure_flight_collection_interaction','figure_fixed_blood_marker_context']:
        img=Image.open(r.OUT/(name+'.png'))
        assert min(img.info['dpi'])>599
        svg=(r.OUT/(name+'.svg')).read_text()
        assert '<text ' in svg and '<image' not in svg
        figure_records.append({'name':name,'pixels':list(img.size),'dpi':img.info['dpi'],'editable_SVG_text_and_vector_points':True})
    validation={'input_checksums':'passed','animal_collapse_integer_counts_and_unique_identifiers':'passed',
                'complete_32_panel_BH_and_72023_evaluable_gene_BH':'recomputed and passed',
                'NB_Wald_p_values_and_numeric_difference_in_differences':'independently recomputed and passed',
                'independent_real_data_design_checks':numeric,'NB_parameterization_audits':nb_audits,
                'highlighted_Circadian_wild_p_recomputed':recomputed_p,'wild_draws_recomputed':draws,
                'wild_HC3_batch_vs_individual_fit_100_draws':'passed',
                'hippocampal_Circadian_effect_omission_stability':stability,'figure_files':figure_records}
    (r.OUT/'numeric_model_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    summary['NB_parameterization_audit']=nb_audits
    summary['normalization_and_bootstrap_conditioning']='Panel expression uses full-deposit median-ratio normalization and log2(normalized count +1), held fixed under exclusions and animal resampling. NB size factors and dispersions are refitted on each included-animal count matrix. Both use a fixed full-deposit mean-count >=10 gene universe.'
    (r.OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('All numeric/model checks passed; independent highlighted wild p=',recomputed_p)


if __name__=='__main__':main()
