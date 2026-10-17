from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Explicit post-result contrast and uncertainty checks; primary families fixed."""
from pathlib import Path
import pickle,json
import numpy as np,pandas as pd
from scipy.stats import t,false_discovery_control
from pydeseq2.ds import DeseqStats
from fit_bulk import bh
OUT=stage_root('core');TAB=OUT/'tables';PRE=OUT/'prepared';ROOT=inputs_root()
mapping=pd.read_csv(TAB/'verified_panel_mapping.csv').fillna('')
rows=[]
for osd in [561,562]:
    with (PRE/f'OSD-{osd}_STAR_dds.pkl').open('rb') as f: dds=pickle.load(f)
    design=np.asarray(dds.obsm['design_matrix']);groups=dds.obs.group.astype(str).to_numpy()
    ann=pd.read_csv(PRE/f'OSD-{osd}_annotation.csv').set_index('ENSEMBL')
    def v(group):return design[groups==group].mean(axis=0)
    for cohort in ['ISS_T','LAR']:
        tag=f'flight_by_age_{cohort}'
        c=v(f'FLT_{cohort}_29')-v(f'GC_{cohort}_29')-v(f'FLT_{cohort}_12')+v(f'GC_{cohort}_12')
        ds=DeseqStats(dds,contrast=c,alpha=.05,quiet=True,n_cpus=threads(),independent_filter=False);ds.summary()
        res=ds.results_df.join(ann);res.to_csv(TAB/f'OSD-{osd}_{tag}_secondary_genome.csv.gz')
        for r in mapping.itertuples():
            hit=res[res.SYMBOL.eq(r.mouse_symbol_primary)] if r.mouse_symbol_primary else res.iloc[:0]
            row={'study':f'OSD-{osd}','contrast':tag,'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,
                 'log2FoldChange':np.nan,'lfcSE':np.nan,'pvalue':np.nan}
            if len(hit)==1:row.update({k:float(hit.iloc[0][k]) for k in ['log2FoldChange','lfcSE','pvalue']})
            rows.append(row)
age=pd.DataFrame(rows);assert len(age)==604;age['q_secondary_age_family']=bh(age.pvalue)
age.to_csv(TAB/'bulk_secondary_age_interactions.csv',index=False)
print('Secondary age interactions',int(age.q_secondary_age_family.lt(.05).sum()),flush=True)

weights={'flight_in_saline':[-1,0,1,0],'flight_in_BuOE':[0,-1,0,1],
 'BuOE_in_ground':[-1,1,0,0],'BuOE_in_flight':[0,0,-1,1],
 'flight_balanced_main':[-.5,-.5,.5,.5],'BuOE_balanced_main':[-.5,.5,-.5,.5],
 'flight_by_BuOE_interaction':[1,-1,-1,1]}
sp=pd.read_csv(TAB/'spatial_panel_complete.csv');wr=[]
for osd,code in [(682,'CA'),(685,'DG'),(698,'FCT'),(699,'CT')]:
    a=pd.read_csv(ROOT/'OSD-682-685-698-699-brain-spatial/data/geo'/f'GSE239336_{code}_GeneExpression_Q3norm.txt.gz',sep='\t',index_col=0).drop(index='NegProbe-WTX')
    y=np.log2(a.to_numpy());f=np.array([c.startswith('FLT') for c in a.columns]);d=np.array(['_TX_' in c for c in a.columns])
    means=np.column_stack([y[:,(f==ff)&(d==dd)].mean(axis=1) for ff,dd in [(0,0),(0,1),(1,0),(1,1)]])
    vars=np.column_stack([y[:,(f==ff)&(d==dd)].var(axis=1,ddof=1)/3 for ff,dd in [(0,0),(0,1),(1,0),(1,1)]])
    for tag,w in weights.items():
        w=np.array(w);variance=(vars*w**2).sum(axis=1)
        df=variance**2/((vars*w**2)**2/2).sum(axis=1);se=np.sqrt(variance);effect=means@w
        p=2*t.sf(abs(effect/se),df)
        g=pd.DataFrame({'mouse_symbol':a.index,'welch_effect':effect,'welch_se':se,'welch_df':df,'welch_pvalue':p})
        sub=sp[(sp.study==f'OSD-{osd}')&(sp.contrast==tag)][['study','contrast','analyte','mouse_symbol']].merge(g,on='mouse_symbol',how='left')
        wr.append(sub)
welch=pd.concat(wr);welch['q_welch_panel_full_family']=bh(welch.welch_pvalue)
welch.to_csv(TAB/'spatial_welch_sensitivity.csv',index=False)
sp=sp.merge(welch,on=['study','contrast','analyte','mouse_symbol'],how='left',validate='one_to_one')
assert np.allclose(sp.log2FoldChange,sp.welch_effect,equal_nan=True)
sp.to_csv(TAB/'spatial_panel_complete.csv',index=False)
print('Welch sensitivity',int(welch.q_welch_panel_full_family.lt(.05).sum()),flush=True)

star=pd.read_csv(TAB/'bulk_panel_complete.csv');rsem=pd.concat([pd.read_csv(TAB/f'OSD-{osd}_rounded_RSEM_panel.csv') for osd in [561,562]])
rsem['q_RSEM_flight_sensitivity_family']=bh(rsem.pvalue)
sens=star[star.contrast_kind=='flight'].merge(rsem[['study','contrast','analyte','log2FoldChange','lfcSE','pvalue','q_RSEM_flight_sensitivity_family']],on=['study','contrast','analyte'],suffixes=('_STAR','_RSEM'),validate='one_to_one')
sens['same_direction']=np.sign(sens.log2FoldChange_STAR)==np.sign(sens.log2FoldChange_RSEM)
sens.to_csv(TAB/'bulk_count_source_sensitivity.csv',index=False)
meta=pd.read_csv(TAB/'bulk_sample_registry.csv');audit=[]
for osd,g in meta.groupby('study'):
    alt=g.group.copy();mask=g.flight.eq('FLT');alt.loc[mask]=g.loc[mask,'flight']+'_'+np.where(g.loc[mask,'Factor Value[Euthanasia Location]'].eq('On ISS'),'ISS_T','LAR')+'_'+g.loc[mask,'age'].astype(str)
    audit.append({'study':osd,'source_group_n':g.groupby('group').size().to_dict(),'location_flight_group_n':alt.value_counts().to_dict(),
        'n_location_conflicts':int(g.location_conflict.sum()),'source_design_rank':g.group.nunique(),'location_design_rank':alt.nunique(),
        'young_terminal_lead_conflicts':int((g.location_conflict & g.age.eq(12) & g.cohort.eq('ISS_T')).sum()),
        'strain_values':g['Characteristics[Strain]'].unique().tolist()})
(TAB/'metadata_specification_audit.json').write_text(json.dumps(audit,indent=2))
