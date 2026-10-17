from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Raw Visium coordinates and counts; section structure and six-donor context."""
from pathlib import Path
import re,json,itertools
import numpy as np,pandas as pd,h5py
from scipy.sparse import csc_matrix,csr_matrix
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
B=stage_root('extended');V1=stage_root('core');rng=np.random.default_rng(20261016)
isa=pd.read_csv(B/'sources/OSD-352_sample_ISA.csv');panel=pd.read_csv(V1/'tables/verified_panel_mapping.csv').fillna('')
sets={a[0]:a[2:] for line in (B/'sources/mh.all.v2025.1.Mm.symbols.gmt').read_text().splitlines() if (a:=line.split('\t'))}
registry=[];results=[];donorvalues=[];spots=[]
for p in sorted((B/'sources/visium_data').glob('*/filtered_feature_bc_matrix.h5')):
    sample=p.parent.name;token=sample.replace('Sample_','');s=isa[isa['Sample Name'].str.contains('_'+token+'_',regex=False)]
    assert len(s)==1;md=s.iloc[0];donor=md['Source Name'];flight=md['Factor Value[Spaceflight]'];label=md['Comment[Sample Naming in Manuscript]']
    pos=pd.read_csv(p.parent/'tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode')
    with h5py.File(p,'r') as h:
        m=h['matrix'];names=m['features/name'][:].astype(str);bar=m['barcodes'][:].astype(str)
        X=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
    assert set(bar).issubset(pos.index);pos=pos.loc[bar];keep=pos.in_tissue.eq(1).to_numpy();pos=pos.iloc[np.flatnonzero(keep)];X=X[:,keep];bar=bar[keep]
    library=np.asarray(X.sum(axis=0)).ravel();assert (library>0).all();n=len(bar)
    # Six nearest spot neighbours, symmetrized and row-normalized; coordinate evidence only.
    xy=pos[['pixel_col','pixel_row']].to_numpy();dist,ix=cKDTree(xy).query(xy,k=7)
    rr=np.repeat(np.arange(n),6);cc=ix[:,1:].ravel();W=csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(n,n));W=(W+W.T).sign();W=W.multiply(1/np.asarray(W.sum(axis=1)).ravel()[:,None]).tocsr()
    geneidx={g:np.flatnonzero(names==g) for g in set(panel.mouse_symbol_primary) if g}
    scores=[];info=[];raws=[]
    for r in panel.itertuples():
        ids=geneidx.get(r.mouse_symbol_primary,[]);v=np.asarray(X[ids].sum(axis=0)).ravel() if len(ids)==1 else np.zeros(n)
        y=np.log2(v/library*1e6+1);test=len(ids)==1 and v.sum()>=20 and (v>0).sum()>=5
        scores.append(y);raws.append(v);info.append({'analyte':r.analyte,'mouse_symbol':r.mouse_symbol_primary,'kind':'panel_gene','raw_total':int(v.sum()),'positive_spots':int((v>0).sum()),'tested':test})
    for name,genes in sets.items():
        ids=np.flatnonzero(np.isin(names,genes));v=np.asarray(X[ids].sum(axis=0)).ravel();y=np.log2(v/library*1e6+1)
        scores.append(y);info.append({'analyte':name,'mouse_symbol':'','kind':'Hallmark_transcript_sum','raw_total':int(v.sum()),'positive_spots':int((v>0).sum()),'tested':len(ids)>=10 and v.sum()>=20})
    Y=np.asarray(scores).T;Z=Y-Y.mean(axis=0);ss=np.sum(Z*Z,axis=0);tested=np.array([d['tested'] for d in info])&(ss>0)
    I=np.divide(np.sum(Z*(W@Z),axis=0),ss,out=np.zeros_like(ss),where=ss>0)
    # Within-section permutations quantify spatial structure, not independent flight replication.
    center=-1/(n-1);exceed=np.zeros(len(info),int)
    for _ in range(499):
        z=Z[rng.permutation(n)];null=np.divide(np.sum(z*(W@z),axis=0),ss,out=np.zeros_like(ss),where=ss>0)
        exceed+=np.abs(null-center)>=np.abs(I-center)-1e-12
    pv=(exceed+1)/500;pv[~tested]=1
    for j,d in enumerate(info):
        results.append(dict(d,section=sample,donor=donor,flight=flight,author_section_label=label,n_spots=n,Moran_I=I[j] if tested[j] else np.nan,spatial_permutation_p=pv[j],score_definition='log2(CPM+1)' if d['kind']=='panel_gene' else 'log2(sum set transcripts / all RNA library * 1e6 + 1)'))
        donorvalues.append({'section':sample,'donor':donor,'flight':flight,'analyte':d['analyte'],'kind':d['kind'],'mean_spatial_score':Y[:,j].mean(),'Moran_I':I[j] if tested[j] else np.nan})
    leadindices={g:next(j for j,r in enumerate(panel.itertuples()) if r.mouse_symbol_primary==g) for g in ['Aif1','Pycard','Fos','Junb']}
    for i,b in enumerate(bar):
        row={'section':sample,'donor':donor,'flight':flight,'barcode':b,'x':xy[i,0],'y':xy[i,1],'library':int(library[i])}
        for g in ['Aif1','Pycard','Fos','Junb']:
            j=leadindices[g];row[g+'_log2CPM']=Y[i,j];row[g+'_raw']=int(raws[j][i])
        spots.append(row)
    registry.append({'section':sample,'donor':donor,'flight':flight,'author_section_label':label,'ISA_sample':md['Sample Name'],'n_spots':n,'RNA_UMIs':int(library.sum()),'n_genes':len(names)})
    print(sample,donor,n,flush=True)
a=pd.DataFrame(results);a['q_spatial_structure_family']=a.groupby('kind').spatial_permutation_p.transform(false_discovery_control);a.to_csv(B/'tables/Visium_spatial_structure_tests.csv',index=False)
reg=pd.DataFrame(registry);assert len(reg)==12 and reg.donor.nunique()==6 and reg.groupby('donor').size().eq(2).all();reg.to_csv(B/'tables/Visium_section_registry.csv',index=False)
pd.DataFrame(spots).to_csv(B/'prepared/Visium_lead_gene_spots.csv.gz',index=False)
# Average technical section scores first; three flight and three ground animals, 20 exact labels.
dv=pd.DataFrame(donorvalues).groupby(['donor','flight','analyte','kind'],as_index=False)[['mean_spatial_score','Moran_I']].mean();dv.to_csv(B/'tables/Visium_donor_mean_scores.csv',index=False)
tests=[]
for (gene,kind),s in dv.groupby(['analyte','kind']):
    s=s.sort_values('donor');f=s.flight.eq('Space Flight').to_numpy();assert len(s)==6 and f.sum()==3
    v=s.mean_spatial_score.to_numpy();effect=v[f].mean()-v[~f].mean();null=[]
    for ids in itertools.combinations(range(6),3):
        mask=np.isin(np.arange(6),ids);null.append(v[mask].mean()-v[~mask].mean())
    tests.append({'analyte':gene,'kind':kind,'mean_score_flight_minus_ground':effect,'exact_animal_p':np.mean(abs(np.array(null))>=abs(effect)-1e-12),'n_exact_labels':20})
t=pd.DataFrame(tests);t['q_animal_family']=t.groupby('kind').exact_animal_p.transform(false_discovery_control);t.to_csv(B/'tables/Visium_donor_flight_tests.csv',index=False)
(B/'tables/Visium_summary.json').write_text(json.dumps({'sections':len(reg),'animals':reg.donor.nunique(),'flight_animals':3,'ground_animals':3,'technical_sections_per_animal':2,'spots':int(reg.n_spots.sum()),'spatial_permutations_per_section':499,'spatial_gene_q05':int(a.loc[a.kind.eq('panel_gene')].q_spatial_structure_family.lt(.05).sum()),'flight_gene_q05':int(t.loc[t.kind.eq('panel_gene')].q_animal_family.lt(.05).sum()),'note':'Spatial structure and donor flight effects are separate hypotheses. No named anatomical compartment assigned without histological/author annotation. Spots and technical sections never inflate exposure n.'},indent=2))
