#!/usr/bin/env python3
"""All-deposit expression structure and count-derived QC; descriptive, not causal tests."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import svd
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from run_revision import DATASETS,load_dataset,check_input_integrity
check_input_integrity();annotation=pd.read_csv(ROOT/'inputs/gene_annotation.csv.gz').set_index('ensembl_id');mt=set(annotation.index[annotation.gene_symbol.str.match(r'(?i)^mt-',na=False)])
QC=[];pcs=[];summaries=[];loadings=[];model=[]
COL={'FLT_ISS_T':'#b14f80','GC_ISS_T':'#177e99','FLT_LAR':'#d49025','GC_LAR':'#687a86','FLT':'#b14f80','GC':'#177e99','HGC':'#177e99','VGC':'#687a86'}
plt.rcParams.update({'font.family':'Arial','font.size':16,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axes=plt.subplots(3,3,figsize=(13,12),layout='constrained')
for k,n in enumerate(DATASETS):
 d=load_dataset(n);meta=d['meta'].copy();y=d['y'];counts=d['counts'];group=(meta.flight+'_'+meta.cohort)if n in ['561','562']else meta.flight
 total=counts.loc[counts.index.astype(str).str.match(r'^ENSMUSG\d+$')].sum()
 mgenes=counts.index.intersection(list(mt));mitoc=counts.loc[mgenes].sum();detected=(counts>0).sum()
 for ix in meta.index:
  QC.append(dict(accession='OSD-'+n,animal=meta.loc[ix,'animal'],mission=meta.loc[ix,'mission'],tissue=meta.loc[ix,'tissue'],group=group.loc[ix],age=meta.loc[ix,'age'],gene_assigned_count_sum=int(total.loc[ix]),n_detected_genes=int(detected.loc[ix]),mitochondrial_gene_count_fraction=float(mitoc.loc[ix]/total.loc[ix]),n_mitochondrial_symbols_matched=len(mgenes),median_ratio_factor=float(d['sf'].loc[ix]),location_conflict=bool(meta.loc[ix,'location_conflict'])))
 v=y.var(axis=1).sort_values(ascending=False,kind='stable');selected=v.index[:min(2000,len(v))];a=y.loc[selected].to_numpy();center=a.mean(axis=1,keepdims=True);m=a-center
 u,s,vh=svd(m,full_matrices=False);score=(s[:,None]*vh).T;variance=s*s/(s*s).sum()
 # Independent covariance eigenspectrum and exact reconstruction checks.
 eigen=np.linalg.eigvalsh(m.T@m)[::-1];assert np.allclose(eigen,s*s,atol=1e-7,rtol=1e-10)
 assert np.allclose(u@np.diag(s)@vh,m,atol=1e-10)
 summaries.append(dict(accession='OSD-'+n,n_animals=len(meta),n_variable_genes=len(selected),PC1_variance_fraction=variance[0],PC2_variance_fraction=variance[1],PC3_variance_fraction=variance[2],mean_gene_assigned_counts=float(total.mean()),median_gene_assigned_counts=float(total.median()),median_mitochondrial_fraction=float((mitoc/total).median())))
 for i,ix in enumerate(meta.index):
  pcs.append(dict(accession='OSD-'+n,animal=meta.loc[ix,'animal'],group=group.loc[ix],age=meta.loc[ix,'age'],location_conflict=bool(meta.loc[ix,'location_conflict']),PC1=score[i,0],PC2=score[i,1],PC3=score[i,2]))
 for i,gene in enumerate(selected):loadings.append(dict(accession='OSD-'+n,ensembl_id=gene,gene_symbol=annotation.reindex([gene]).gene_symbol.iloc[0],PC1_loading=u[i,0],PC2_loading=u[i,1]))
 xage=pd.get_dummies(meta.age.astype(str),drop_first=True,dtype=float).to_numpy();cond=pd.get_dummies(group,drop_first=True,dtype=float).to_numpy();x0=np.column_stack([np.ones(len(meta)),xage]);x=np.column_stack([x0,cond]);coef=np.linalg.lstsq(x,m.T,rcond=None)[0];e=m.T-x@coef;coef0=np.linalg.lstsq(x0,m.T,rcond=None)[0];e0=m.T-x0@coef0;partial=1-(e*e).sum()/(e0*e0).sum()
 model.append(dict(accession='OSD-'+n,age_rank=int(np.linalg.matrix_rank(x0)),full_rank=int(np.linalg.matrix_rank(x)),condition_df=int(np.linalg.matrix_rank(x)-np.linalg.matrix_rank(x0)),condition_partial_R2_top_variable_genes=float(partial),n_genes=len(selected),inference='Descriptive multivariate R2 conditional on selected log-expression universe; no hypothesis or causal test'))
 ax=axes.flat[k]
 for grp in group.unique():
  for age in meta.age.unique():
   mask=(group==grp)&(meta.age==age);idx=np.where(mask)[0];color=COL[grp];ax.scatter(score[idx,0],score[idx,1],color=color,marker='s'if str(age)=='29'else 'o',s=46,edgecolors=np.where(meta.loc[mask,'location_conflict'],'black',color),lw=1.3,label=grp.replace('_',' ')if str(age)==str(meta.age.unique()[0])else None)
 ax.set_xlabel(f'PC1 ({variance[0]*100:.1f}%)');ax.set_ylabel(f'PC2 ({variance[1]*100:.1f}%)');ax.set_title(chr(97+k)+'  OSD-'+n+' '+d['meta'].tissue.iloc[0].replace('left cerebral hemisphere','Left hemisphere').replace('cerebellum','Cb').replace('hippocampus','Hip'),loc='left',fontweight='bold',fontsize=16)
 ax.legend(frameon=False,fontsize=15,loc='best',ncol=2)
q=pd.DataFrame(QC);s=pd.DataFrame(summaries);pm=pd.DataFrame(model);ax=axes.flat[7]
ax.bar(range(7),pm.condition_partial_R2_top_variable_genes,color='#177e99');ax.set_xticks(range(7),[n.replace('OSD-','')for n in pm.accession]);ax.set_ylim(0,1);ax.set_ylabel('Condition partial R²');ax.set_xlabel('OSD deposit');ax.set_title('h  Descriptive condition structure',loc='left',fontweight='bold',fontsize=16)
ax=axes.flat[8]
for i,acc in enumerate(s.accession):
 part=q[q.accession==acc];ax.scatter(np.full(len(part),i)+np.linspace(-.14,.14,len(part)),part.mitochondrial_gene_count_fraction*100,color='#687a86',s=20)
ax.set_xticks(range(7),[n.replace('OSD-','')for n in s.accession]);ax.set_xlabel('OSD deposit');ax.set_ylabel('Mitochondrial gene counts (%)');ax.set_title('i  Count-composition diagnostic',loc='left',fontweight='bold',fontsize=16)
fig.savefig(OUT/'figureS2_all_deposit_structure.png',dpi=600,bbox_inches='tight');fig.savefig(OUT/'figureS2_all_deposit_structure.svg',bbox_inches='tight');plt.close(fig)
q.to_csv(OUT/'animal_count_qc.csv',index=False);pd.DataFrame(pcs).to_csv(OUT/'animal_pca_scores.csv',index=False);s.to_csv(OUT/'deposit_structure_summary.csv',index=False);pd.DataFrame(loadings).to_csv(OUT/'pca_gene_loadings.csv',index=False);pm.to_csv(OUT/'condition_structure_R2.csv',index=False)
# Count-composition checks from raw counts, not an RNA-integrity assay.
assert q.gene_assigned_count_sum.gt(0).all() and q.mitochondrial_gene_count_fraction.between(0,1).all();assert len(q)==126
(OUT/'structure_validation.json').write_text(json.dumps(dict(all_deposits_included=True,animal_region_rows=len(q),PCA_eigenspectrum_and_reconstruction_independently_checked=True,input_checksums_checked=True,inference='descriptive only',RNA_integrity_numbers='Not present in the bundled metadata; count composition is not a substitute for RIN',font='Arial',dpi=600),indent=2)+'\n')
(OUT/'STRUCTURE_SUMMARY.md').write_text('All seven deposits were examined using gene-assigned count totals, detected genes, mitochondrial count composition and top-2,000-variable-gene PCA. These are descriptive diagnostics, not evidence of mitochondrial injury, RNA integrity or a causal collection mechanism. PC fractions and conditional multivariate condition R² are exported for every deposit. Age is represented in the nuisance design where available, and repeated regions remain distinct observations rather than independent mission replicates. The bundled metadata contain no RNA-integrity-number measurement.\n\n'+s.to_string(index=False)+'\n\n'+pm.to_string(index=False)+'\n')
print(s.to_string(index=False));print(pm.to_string(index=False))
