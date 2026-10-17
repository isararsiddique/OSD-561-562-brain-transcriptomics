from pipeline.paths import stage_root, results_root, inputs_root, code_root, publication_font, threads, network_allowed, request_get, request_post, record_path, resolve_record
"""Publication workflow: questions, verified units, analyses and evidence."""
from pathlib import Path
import json, hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image
B=stage_root('deep')
R=stage_root('extended')
def draw_workflow():
    plt.rcParams.update({'font.family':publication_font(),'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False})
    fig=plt.figure(figsize=(16,14));ax=fig.add_axes([.025,.025,.95,.95])
    ax.set(xlim=(0,16),ylim=(0,14));ax.axis('off')
    accent='#226B73';ink='#000000';sub='#4D4D4D';used=[]
    def box(x,y,w,h,fill='#F5F7F8',edge='#C7CED2'):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.055,rounding_size=.10',fc=fill,ec=edge,lw=.85))
    def arrow(x1,y1,x2,y2):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=12,color=accent,lw=1.25))
    def txt(x,y,text,size=10.5,weight='normal',color=ink,**kw):
        ax.text(x,y,text,fontsize=size,fontweight=weight,color=color,va='top',**kw)
    txt(.18,13.77,'Spaceflight brain: study workflow',23,'bold')
    txt(.18,13.23,'Biological question: do myeloid RNA changes reflect cell mixture, transcriptional remodeling, or sample handling?',12)
    txt(.18,12.60,'A  Verified data and biological replicates',14,'bold')
    xs=[.20,5.55,10.90]
    for x in xs:box(x,10.28,4.88,1.97)
    txt(.38,12.10,'Original project',12,'bold')
    txt(.38,11.62,'OSD-561 / 562: cerebellum + hippocampus\n53 libraries · 29 mice · 24 mice with paired regions\n\nOSD-682 / 685 / 698 / 699: GeoMx\n4 regions · 48 ROIs · animal links unresolved',10.6)
    txt(5.73,12.10,'Mission comparisons',12,'bold')
    txt(5.73,11.62,'RR8: 47 libraries · 8 mice · 122,084 cells (GEO gate)\nRR10: 10 paired mice · RR9: 12 mice\n\nAdded: RR23 cerebellum · 26 mice\nAdded: MHU-3 cerebrum · 24 mice (WT / Nrf2KO)\nAdded: RRRM-2 cerebral region · 18 different mice',10.6)
    txt(11.08,12.10,'RR3 multimodal context',12,'bold')
    txt(11.08,11.62,'OSD-352: 21,178 nuclei · 5 mice\nRNA + ATAC measured from the same nuclei\n\nVisium: 29,839 spots · 6 mice\n12 technical sections · verified coordinates',10.6)
    box(.20,9.45,15.58,.52,fill='#EAF2F2',edge='#93B7BA')
    txt(7.99,9.81,'151 corrected human targets → 143 primary mouse symbols | ISA identities · count QC · technical-library collapse',11.3,'bold',ha='center')
    for x in [2.64,7.99,13.34]:arrow(x,10.25,x,10.00)
    txt(.18,8.92,'B  Five biological and analytical evidence branches',14,'bold')
    txt(.35,8.48,'QUESTION',10,'bold',color=sub)
    txt(4.05,8.48,'DATA AND ANALYSIS',10,'bold',color=sub)
    txt(10.93,8.48,'CURRENT EVIDENCE',10,'bold',color=sub)
    lanes=[
      ('1  Biological specificity','Which sensor, adaptor and\nexecution transcripts change?',
       'RR8: six enriched populations · whole-genome count models\nedgeR QL + limma/voom · GEO, article and alternative QC\nAll 151 targets × populations × age contrasts retained',
       'Forebrain candidate RNA program:\nAIM2 higher; PYCARD, IL18 and NINJ1 lower.'),
      ('2  State or cell mixture?','Does the signal persist in\ncomparable expression states?',
       'Non-panel reference mixtures: held-donor NNLS / ridge\nNew: balanced, exposure-blind PCA states (K=4 / K=6)\nMouse-state pseudobulk · cell coverage · fixed-composition checks',
       'PYCARD persists in primary shared states.\nTechnical / finer-state sensitivity limits the claim.'),
      ('3  Transfer across missions','Does the association transfer\nacross region, strain and control?',
       'RR3 / RR10 / RR9 context + RR23 / MHU-3 / cerebral counts\nHabitat and vivarium contrasts kept separate\nNrf2 genotype interaction · 50 complete mouse-deletion fits',
       'RR23 AIM2: flight > habitat, deletion-stable.\nPYCARD corrected support varies by control / family.'),
      ('4  Spatial and regulatory support','Is there regional localization\nand promoter agreement?',
       'Original GeoMx contrasts + RR3 donor-linked Visium\nSpatial neighbours separated from exposure tests\nDonor RNA and promoter ATAC · same-mission dependence',
       'Local spatial structure is present.\nCorrected spatial / promoter flight effects unconfirmed.'),
      ('5  Inferential robustness','Which results survive\nmethod and assumption changes?',
       'Three bulk engines · complete planned testing families\nMouse deletion · HC3 intervals · exact animal permutations\n50 Hallmark sets · correlation and residual-df sensitivity',
       'Report effects, uncertainty and failed checks.\nRNA associations remain observational.'),
    ]
    yvals=[7.04,5.77,4.50,3.23,1.96]
    for y,(heading,question,method,result) in zip(yvals,lanes):
        box(.20,y,3.35,1.11,fill='white')
        box(3.89,y,6.55,1.11)
        box(10.78,y,5.00,1.11,fill='#F0F6F4',edge='#B4C7BF')
        txt(.36,y+.97,heading,11.7,'bold')
        txt(.36,y+.58,question,10.6)
        txt(4.05,y+.91,method,10.3)
        txt(10.94,y+.84,result,10.4)
        arrow(3.57,y+.55,3.85,y+.55);arrow(10.46,y+.55,10.74,y+.55)
    box(.20,.60,15.58,.98,fill='white',edge='#7C8E83')
    txt(.38,1.45,'C  Next mechanistic test · proposed laboratory validation',12.5,'bold')
    txt(.38,1.08,'Matched handling + independently identified myeloid cells → multiplex RNA / protein → ASC, caspase and cytokine readouts → perturbation / rescue',10.7)
    txt(.38,.78,'Aim: test sensor–adaptor divergence. Computational evidence does not establish inflammasome activity or a causal protective response.',10.4,color=sub)
    txt(.20,.32,'Limits carried through every branch: terminal bulk preservation; RR8 QC threshold / dissection day; RR23 collection dates; small mouse cohorts and unresolved GeoMx animal identities.',9.6,color=sub)
    out=R/'figures/figure1_illustrated_study_workflow'
    for ext in ['png','pdf','svg']:fig.savefig(out.with_suffix('.'+ext),dpi=210,facecolor='white')
    plt.close(fig)
    (R/'figures/illustration_provenance.json').write_text(json.dumps(used,indent=2))
    sources=['../core/tables/verified_panel_mapping.csv','../core/tables/bulk_sample_registry.csv',
             'tables/sorted_cell_preparation_summary.json','tables/multiome_aggregation_summary.json','tables/Visium_section_registry.csv',
             '../deep/tables/additional_mouse_registry.csv','../deep/tables/state_panel_complete.csv',
             '../deep/tables/additional_panel_all_regions.csv','../deep/tables/additional_mouse_deletion_expanded_focus.csv','../deep/tables/terminal_region_mouse_registry.csv','../deep/tables/article_GEO_QC_discrepancy.json']
    (B/'tables/workflow_source_ledger.json').write_text(json.dumps({'figure':str(out),'sources':sources,'script':__file__,'independent_units':'Mice; cells, sections and technical libraries do not add exposure replicates'},indent=2))
    return sources
if __name__=='__main__':draw_workflow()
