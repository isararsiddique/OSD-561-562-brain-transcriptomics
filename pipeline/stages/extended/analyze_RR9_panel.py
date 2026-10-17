from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Independent RR9 nCounter counts, control normalization and exact animal tests."""
from pathlib import Path
import gzip,io,re,itertools,json
import numpy as np,pandas as pd
from scipy.stats import gmean,false_discovery_control,ttest_ind
B=stage_root('extended');V1=stage_root('core')
isa=pd.read_csv(B/'sources/OSD-536_sample_ISA.csv').set_index('Sample Name')
df=[];qc=[]
for p in sorted((B/'sources/OSD-536').glob('*.RCC.gz')):
    gsm=re.search(r'GSM\d+',p.name).group();text=gzip.open(p,'rt').read()
    x=pd.read_csv(io.StringIO(text.split('<Code_Summary>')[1].split('</Code_Summary>')[0].strip())).set_index('Name')
    x['sample']=gsm;df.append(x.reset_index())
    n=x[x.CodeClass.eq('Negative')].Count;pos=x[x.CodeClass.eq('Positive')].Count
    qc.append({'sample':gsm,'animal':isa.loc[gsm,'Source Name'],'flight':isa.loc[gsm,'Factor Value[Spaceflight]'],'negative_LOD_mean_plus_2SD':n.mean()+2*n.std(ddof=1),'positive_geomean':gmean(pos),'positive_zero_probes':int(pos.eq(0).sum())})
long=pd.concat(df,ignore_index=True);raw=long.pivot(index='Name',columns='sample',values='Count');classes=long.drop_duplicates('Name').set_index('Name').CodeClass
q=pd.DataFrame(qc).set_index('sample').loc[raw.columns];posfactor=q.positive_geomean.pipe(lambda z:gmean(z)/z)
hk=raw.loc[classes.eq('Housekeeping')];valid=hk.index[hk.gt(0).all(axis=1)];assert len(valid)>=6
positive=raw*posfactor;hkgeomean=gmean(positive.loc[valid].to_numpy(),axis=0);contentfactor=gmean(hkgeomean)/hkgeomean
normalized=positive*contentfactor;normalized.to_csv(B/'prepared/OSD536_nCounter_normalized.csv')
raw.to_csv(B/'prepared/OSD536_nCounter_raw.csv');q['positive_scale']=posfactor;q['HK_scale']=contentfactor;q.to_csv(B/'tables/OSD536_nCounter_QC.csv')
f=q.flight.eq('Space Flight').to_numpy();assert f.sum()==6 and (~f).sum()==6
permutations=np.array([np.isin(np.arange(12),ids) for ids in itertools.combinations(range(12),6)])
mapping=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('');rows=[]
for r in mapping.itertuples():
    row={'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'n_flight':6,'n_ground':6,'detected_above_LOD_samples':0,'log2_effect':np.nan,'welch_p':1.,'exact_p':1.,'status':'not_in_assay_or_unmapped'}
    if r.mouse_symbol_primary in raw.index:
        gene=r.mouse_symbol_primary;detected=(raw.loc[gene]>q.negative_LOD_mean_plus_2SD).sum();row['detected_above_LOD_samples']=int(detected)
        if detected>=6:
            v=np.log2(normalized.loc[gene].to_numpy()+1);effect=v[f].mean()-v[~f].mean();null=permutations@v/6-(~permutations)@v/6
            row.update(log2_effect=effect,welch_p=ttest_ind(v[f],v[~f],equal_var=False).pvalue,exact_p=np.mean(abs(null)>=abs(effect)-1e-12),status='tested_at_least_6_samples_above_negative_probe_LOD')
        else:row['status']='below_fixed_detection_gate'
    rows.append(row)
a=pd.DataFrame(rows);a['q_full_151_exact']=false_discovery_control(a.exact_p);a['q_full_151_welch']=false_discovery_control(a.welch_p.fillna(1));a.to_csv(B/'tables/OSD536_panel_results.csv',index=False)
(B/'tables/OSD536_analysis_summary.json').write_text(json.dumps({'assay':'nCounter Mouse Neuroinflammation panel','endogenous_probes':int(classes.eq('Endogenous').sum()),'HK_probes':13,'all_samples_positive_HK_used':valid.tolist(),'normalization':'Positive-control geometric mean, then all nonzero-housekeeping geometric mean. This is transparent reanalysis, not an exact nSolver reproduction; publication selected HK subset unspecified.','fixed_gate':'Raw count > sample negative-probe mean+2SD in >=6/12 samples','exact_labels':len(permutations),'tested_CNS_targets':int(a.status.str.startswith('tested').sum()),'exact_q05':int(a.q_full_151_exact.lt(.05).sum()),'note':'RR9 male right caudal hemibrain (mid- and hindbrain, publication Methods) post-return data; different sex, region and timing from core study.'},indent=2));print(a[a.mouse_symbol.isin(['Aif1','Pycard','Fos','Junb'])].to_string(index=False),flush=True)
