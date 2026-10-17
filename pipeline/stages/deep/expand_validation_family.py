from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Preserve first-batch correction; additionally audit all new regional tests."""
from pathlib import Path
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
B=stage_root('deep')
first=pd.read_csv(B/'tables/additional_panel_complete.csv')
regions=pd.read_csv(B/'tables/terminal_region_panel_complete.csv')
a=pd.concat([first,regions],ignore_index=True)
for method in ['limma_voom','edgeR_QL']:
    idx=a.method.eq(method);assert idx.sum()==1812
    a.loc[idx,'q_expanded_1812_family']=multipletests(a.loc[idx,'pvalue'],method='fdr_bh')[1]
a.to_csv(B/'tables/additional_panel_all_regions.csv',index=False)
old=pd.read_csv(B/'tables/additional_mouse_deletion_complete.csv.gz')
rows=[];focus=[]
for (drop_sid,drop_mouse,method),s in old.groupby(['deleted_OSD','deleted_mouse','method']):
    extra=regions[regions.method.eq(method)].copy();extra['deleted_OSD']=drop_sid;extra['deleted_mouse']=drop_mouse
    combined=pd.concat([s,extra],ignore_index=True);assert len(combined)==1812
    combined['q_expanded_1812_family']=multipletests(combined.pvalue,method='fdr_bh')[1]
    rows.append(combined)
    focus.append(combined[combined.OSD.eq(drop_sid)&combined.analyte.isin(['AIM2','AIF1','PYCARD','SLC1A3','S100B'])])
pd.concat(rows).to_csv(B/'tables/additional_mouse_deletion_expanded_complete.csv.gz',index=False)
pd.concat(focus).to_csv(B/'tables/additional_mouse_deletion_expanded_focus.csv',index=False)
print(a[a.contrast.isin(['flight_vs_habitat_ground','flight_vs_vivarium','flight_WT','flight_balanced_genotype'])&a.analyte.isin(['AIM2','PYCARD','SLC1A3','S100B'])][['OSD','method','contrast','analyte','log2FC','q_full_1057_family','q_expanded_1812_family']].to_string(index=False))
