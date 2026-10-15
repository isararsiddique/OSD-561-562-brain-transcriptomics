"""Query definitions exactly as registered in H_PROTOCOL.md section 4-5."""

SPACE = '(spaceflight OR "space flight" OR microgravity OR "International Space Station" OR "Rodent Research" OR RRRM OR "rodent habitat" OR "space station")'
COLLECT = '(euthanasia OR euthanized OR euthanised OR "sample collection" OR "collection technique" OR "collection protocol" OR "dissection" OR "live animal return" OR "on-orbit" OR "on orbit" OR "ISS terminal" OR ketamine OR xylazine OR isoflurane OR anesthesia OR anaesthesia OR "tissue preservation" OR "carcass")'
BRAIN = '(brain OR hippocampus OR hippocampal OR cerebellum OR cerebellar OR cortex OR "nervous system")'
IEG = '("immediate early gene" OR "immediate-early gene" OR "immediate early genes" OR "immediate-early genes" OR "primary response gene" OR "primary response genes" OR "c-fos" OR Fos OR Egr1 OR Npas4 OR Junb OR Fosb OR Nr4a1 OR Arc)'
OMICS = '(transcriptome OR transcriptomic OR transcriptomics OR "RNA-seq" OR "RNA sequencing" OR "gene expression" OR mRNA OR microarray)'
RODENT = '(mouse OR mice OR rat OR rats OR rodent OR rodents OR murine)'
CONF = '(confound OR confounded OR confounding OR artifact OR artefact OR artifactual OR "sample collection" OR "collection technique" OR "collection method" OR dissection OR euthanasia OR "live animal return" OR "on-orbit" OR "tissue preservation")'
EUTH = '(euthanasia OR "euthanasia method" OR "method of euthanasia" OR "killing method" OR "sacrifice" OR "terminal anesthesia" OR "terminal anaesthesia" OR "tissue collection" OR "sample collection" OR "post-mortem interval" OR "postmortem interval" OR "harvest")'
ANAES = '(ketamine OR xylazine OR isoflurane OR anesthesia OR anaesthesia OR anesthetic OR anaesthetic OR pentobarbital OR "carbon dioxide" OR CO2 OR decapitation)'
CORR = '("batch effect" OR "batch effects" OR "batch correction" OR "batch-effect correction" OR "unwanted variation" OR RUV OR "surrogate variable" OR "surrogate variables" OR ComBat OR "confounder adjustment" OR "latent confounder" OR "latent confounders" OR "hidden confounders" OR "data harmonization" OR "removing unwanted")'
GAINN = '("operator norm" OR "spectral norm" OR "spectral gain" OR "amplification factor" OR "noise amplification" OR "error amplification" OR "oblique projection" OR "Lipschitz constant" OR "largest singular value")'
GAINW = GAINN[:-1] + ' OR "variance inflation" OR "condition number" OR "inflates variance" OR "noise inflation")'
OMICS2 = OMICS[:-1] + ' OR genomic OR genomics OR omics OR methylation OR proteomic OR proteomics OR "single-cell")'
FW = '("Frank-Wolfe" OR "Frank Wolfe" OR "conditional gradient" OR "duality gap")'
HOLD = '(holdout OR "hold-out" OR "held-out" OR "held out" OR "negative control samples" OR "control samples")'


def TA(x):
    return '(TITLE:%s OR ABSTRACT:%s)' % (x, x)


C6 = '(euthanasia OR "sample collection" OR "tissue collection" OR sampling OR artifact OR artefact OR artifactual OR confound OR confounding OR "post-mortem" OR postmortem)'

EPMC = {
    'C1': [
        ('C1-1', f'{TA(SPACE)} AND {TA(COLLECT)} AND ({TA(BRAIN)} OR {TA(IEG)})'),
        ('C1-2', f'{TA(SPACE)} AND {TA(CONF)} AND {TA(OMICS)}'),
        ('C1-3', f'{TA(SPACE)} AND {TA(BRAIN)} AND {TA(OMICS)} AND {RODENT}'),
        ('C1-4', f'{TA(SPACE)} AND {IEG} AND {COLLECT} AND {BRAIN} AND {RODENT}'),
        ('C1-5', f'{TA(EUTH)} AND {TA(BRAIN)} AND ({TA(IEG)} OR {TA(OMICS)}) AND {RODENT}'),
        ('C1-6', f'{TA(ANAES)} AND {TA(IEG)} AND {TA(C6)} AND {RODENT}'),
        ('C1-7a', f'SRC:PPR AND {TA(SPACE)} AND {TA(COLLECT)} AND ({TA(BRAIN)} OR {TA(IEG)})'),
        ('C1-7b', f'SRC:PPR AND {TA(SPACE)} AND {TA(CONF)} AND {TA(OMICS)}'),
        ('C1-8', f'METHODS:{SPACE} AND METHODS:(euthanasia OR "live animal return") AND RESULTS:{IEG}'),
    ],
    'C2': [
        ('C2-1', f'{TA(CORR)} AND {GAINN} AND {OMICS2}'),
        ('C2-2', f'{TA(CORR)} AND {TA(GAINW)}'),
        ('C2-3', f'{TA(CORR)} AND {TA("(loading OR loadings)")} AND {TA(chr(40) + "precision OR imprecise OR " + chr(34) + "estimation error" + chr(34) + " OR uncertainty OR " + chr(34) + "sample size" + chr(34) + ")")}'),
        ('C2-4', f'{TA(CORR)} AND {TA(chr(40) + "\"false positive\" OR \"false positives\" OR \"exaggerated confidence\" OR overcorrection OR \"over-correction\" OR inflate OR inflated OR pitfall OR pitfalls)")} AND {TA(OMICS2)}'),
        ('C2-5', f'{TA(chr(40) + "\"unwanted variation\" OR \"confounder adjustment\" OR \"surrogate variable\" OR \"latent confounder\" OR \"latent confounders\" OR \"negative controls\")")} AND {TA(chr(40) + "\"multiple testing\" OR \"multiple hypothesis\" OR identifiability OR identifiable OR \"negative control\")")}'),
        ('C2-6', f'{TA(CORR)} AND (detectability OR removability OR "cannot be removed" OR "not removable") AND (projection OR "regress out" OR "regressing out")'),
    ],
    'C3': [
        ('C3-1a', f'{FW} AND {TA(CORR)}'),
        ('C3-1b', f'{FW} AND {TA(OMICS2)} AND (confounder OR batch)'),
        ('C3-2', f'("spectral norm" OR "operator norm" OR "Lipschitz") AND (constraint OR constrained OR cap OR bound) AND {TA(CORR)}'),
        ('C3-3', f'{TA(SPACE)} AND {TA(chr(40) + "batch OR normalization OR normalisation OR \"cross-validation\" OR holdout OR \"held-out\" OR \"ground control\" OR \"ground controls\" OR confound OR confounding)")} AND {TA(OMICS2)}'),
        ('C3-4', f'"noise ceiling" AND ({TA(CORR)} OR {TA(OMICS)})'),
        ('C3-5', f'{TA(CORR)} AND {TA("(benchmark OR evaluation OR validation OR assessment)")} AND {HOLD} AND ("training set" OR "nested cross-validation" OR "data leakage" OR leakage)'),
    ],
    'C4': [
        ('C4-2', f'("RRRM-2" OR RRRM2 OR "Rodent Research Reference Mission") AND ({BRAIN[1:-1]} OR "nervous system" OR neural)'),
        ('C4-3', f'("RR-10" OR "Rodent Research-10" OR "Rodent Research 10" OR RR10) AND {SPACE} AND {BRAIN}'),
        ('C4-4', f'("RR-23" OR "Rodent Research-23" OR "Rodent Research 23" OR RR23) AND {SPACE} AND ({BRAIN[1:-1]} OR cerebellum)'),
    ],
}
ACC = ['OSD-561', 'OSD-562', 'OSD-563', 'OSD-564', 'OSD-612', 'OSD-613', 'OSD-525',
       'GLDS-556', 'GLDS-557', 'GLDS-568', 'GLDS-569', 'GLDS-588', 'GLDS-589', 'GLDS-525']
for a in ACC:
    EPMC['C4'].append(('C4-1:' + a, '"%s"' % a))

# OpenAlex: same boolean strings, field prefixes removed, in title_and_abstract.search (TA) or fulltext.search (FT)
OA = {
    'C1': [('OA-C1-1', 'TA', f'{SPACE} AND {COLLECT} AND ({BRAIN} OR {IEG})'),
           ('OA-C1-2', 'TA', f'{SPACE} AND {CONF} AND {OMICS}'),
           ('OA-C1-3', 'TA', f'{SPACE} AND {BRAIN} AND {OMICS} AND {RODENT}'),
           ('OA-C1-5', 'TA', f'{EUTH} AND {BRAIN} AND ({IEG} OR {OMICS}) AND {RODENT}'),
           ('OA-C1-6', 'TA', f'{ANAES} AND {IEG} AND {C6} AND {RODENT}')],
    'C2': [('OA-C2-1', 'FT', f'{CORR} AND {GAINN} AND {OMICS2}'),
           ('OA-C2-2', 'TA', f'{CORR} AND {GAINW}'),
           ('OA-C2-3', 'TA', f'{CORR} AND (loading OR loadings) AND (precision OR imprecise OR "estimation error" OR uncertainty OR "sample size")'),
           ('OA-C2-4', 'TA', f'{CORR} AND ("false positive" OR "false positives" OR "exaggerated confidence" OR overcorrection OR "over-correction" OR inflate OR inflated OR pitfall OR pitfalls) AND {OMICS2}'),
           ('OA-C2-5', 'TA', '("unwanted variation" OR "confounder adjustment" OR "surrogate variable" OR "latent confounder" OR "latent confounders" OR "negative controls") AND ("multiple testing" OR "multiple hypothesis" OR identifiability OR identifiable OR "negative control")')],
    'C3': [('OA-C3-1', 'FT', f'{FW} AND {CORR}'),
           ('OA-C3-2', 'FT', f'("spectral norm" OR "operator norm" OR Lipschitz) AND (constraint OR constrained) AND {CORR}'),
           ('OA-C3-3', 'TA', f'{SPACE} AND (batch OR normalization OR normalisation OR "cross-validation" OR holdout OR "held-out" OR "ground control" OR "ground controls" OR confound OR confounding) AND {OMICS2}'),
           ('OA-C3-4', 'FT', f'"noise ceiling" AND ({CORR} OR {OMICS})')],
    'C4': [('OA-C4-2', 'TA', f'("RRRM-2" OR RRRM2 OR "Rodent Research Reference Mission") AND {BRAIN}'),
           ('OA-C4-3', 'TA', f'("RR-10" OR "Rodent Research-10" OR "Rodent Research 10") AND {BRAIN}'),
           ('OA-C4-4', 'TA', f'("RR-23" OR "Rodent Research-23" OR "Rodent Research 23") AND {BRAIN}')]
           + [('OA-C4-1:' + a, 'FT', '"%s"' % a) for a in ACC],
}

CR = {
    'C1': ['spaceflight mouse brain euthanasia immediate early gene expression',
           'sample collection technique confounds rodent spaceflight RNA-seq',
           'ISS terminal versus live animal return mouse tissue transcriptome',
           'on-orbit euthanasia dissection rodent research spaceflight',
           'spaceflight mouse brain hippocampus cerebellum transcriptome',
           'anesthesia for euthanasia influences brain mRNA expression mice',
           'ketamine xylazine immediate early gene expression brain sampling artifact'],
    'C2': ['batch correction operator norm amplification gene expression',
           'removing unwanted variation negative controls gene expression',
           'confounder adjustment multiple hypothesis testing latent factors',
           'batch effect removal exaggerated confidence false positives',
           'oblique projection noise amplification regression adjustment genomics',
           'batch effect detectable but not removable'],
    'C3': ['Frank-Wolfe spectral norm constrained batch correction',
           'Lipschitz constrained batch correction gene expression',
           'validation of batch effect correction held-out control samples',
           'spaceflight transcriptomics batch correction ground control validation'],
    'C4': ['RRRM-2 mouse brain transcriptomics', 'Rodent Research Reference Mission brain',
           'RR-10 mission mouse brain spaceflight', 'RR-23 mouse cerebellum spaceflight',
           'transcriptional profiling of cerebellums from mice flown on the RRRM-2 mission',
           'spatially resolved transcriptional profiling hippocampi mice flown RR-10 mission',
           'single nuclei transcriptomics epigenomics brain tissue mice flown RRRM-2'],
}

CONTROLS = {
    'LaiPolo2020': ('10.1016/j.isci.2020.101733', 'C1'),
    'Grandke2026': ('10.1038/s41467-026-68737-1', 'C1'),
    'Overbey2022': ('10.1016/j.crmeth.2022.100325', 'C1'),
    'Choi2020': ('10.1038/s41598-020-58898-4', 'C1'),
    'Kremsky2023': ('10.3390/ijms241713569', 'C1'),
    'Masarapu2024': ('10.1038/s41467-024-48916-8', 'C1'),
    'Holley2022': ('10.1038/s41526-022-00217-4', 'C1'),
    'StaibLasarzik2014': ('10.1089/neu.2013.3243', 'C1'),
    'Nygaard2016': ('10.1093/biostatistics/kxv027', 'C2'),
    'Zindler2020': ('10.1186/s12859-020-03559-6', 'C2'),
    'Gerard2021': ('10.5705/ss.202018.0345', 'C2'),
    'Wang2017': ('10.1214/16-aos1511', 'C2'),
}
