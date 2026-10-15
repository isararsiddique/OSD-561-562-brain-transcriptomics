"""Build B_closest_prior_art.json and verify every quoted sentence against the cached source text.

Each quote is checked by whitespace-free, HTML-unescaped, case-insensitive substring match against the
file it was read from (abstract record, full text or publisher page). Ellipses '...' split a quote into
fragments that must each match. Entries whose source is metadata only carry quote = null.
"""
import glob
import html
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
C = os.path.join(ROOT, "cache")
PKG = (""
       "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015/extensions/literature_novelty/cache")


def g(pat):
    hits = sorted(glob.glob(os.path.join(C, pat)))
    return hits[0] if hits else None


def text_of(path, key=None):
    raw = open(path, errors="ignore").read()
    if key and path.endswith(".json"):
        d = json.loads(raw)
        raw = d[key]["abstract"] if key in d else raw
    raw = html.unescape(re.sub(r"<[^>]+>", " ", raw))
    return raw


def norm(s):
    s = html.unescape(s).lower()
    s = s.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-").replace("\u00a0", " ")
    return re.sub(r"[\s\"']+", "", s)


E = [
 dict(id="laipolo2020", doi="10.1016/j.isci.2020.101733", pmid="33376967", year=2020,
      title="RNAseq Analysis of Rodent Spaceflight Experiments Is Confounded by Sample Collection Techniques",
      source_read="full text (Europe PMC fullTextXML PMC7756143)", src=("laipolo2020_plain.txt", None),
      showed="RR-1 and follow-on ground experiments in liver and quadriceps: carcass (slow-frozen) versus immediate dissection changes eclipsed flight-ground differences, worsened by polyA selection; ketamine/xylazine versus CO2 euthanasia produced few or no DEGs in liver or quadriceps. Brain was not studied; authors caution against extrapolating to other tissues.",
      quote="Changes due to preservation method eclipsed those between flight and ground samples, making it difficult to identify spaceflight-specific changes. ... These data suggest that the types of euthanasia methods evaluated here do not impact gene expression in select tissues. ... thus extrapolation to other tissues should be done with caution",
      relation="Overturns any 'first to show collection confounding in rodent spaceflight RNA-seq' claim (E1, M1). Narrows E1 to brain, to flight-naive ground animals under the two RRRM-2 protocols, and to a rapid IEG programme; its null euthanasia result (liver/muscle) means the agent cannot be assumed causal in brain.",
      claims=["E1", "M1", "P1"]),
 dict(id="choi2020", doi="10.1038/s41598-020-58898-4", pmid="32047211", year=2020,
      title="Validation of a New Rodent Experimental System to Investigate Consequences of Long Duration Space Habitation",
      source_read="abstract (Europe PMC)", src=("epmc_abs_choi2020_*.json", "abstract_raw"),
      showed="RR-1 validation: female C57BL/6NTac and C57BL/6J mice 21-37 d on ISS, euthanized on orbit; liver and spleen dissected on orbit, remaining tissues frozen in situ; RNA quality acceptable; no standard chronic-stress indices triggered.",
      quote="Mice were euthanized on-orbit, livers and spleens dissected, and remaining tissues frozen in situ for later analyses.",
      relation="Establishes the on-orbit (ISS-terminal) collection capability whose protocol differs from Earth collection; no brain transcriptome or protocol comparison. Context for E1/P1.",
      claims=["E1", "P1"]),
 dict(id="overbey2022", doi="10.1016/j.crmeth.2022.100325", pmid="36452864", year=2022,
      title="Challenges and considerations for single-cell and spatially resolved transcriptomics sample collection during spaceflight",
      source_read="full text (package cache collection_fulltext.xml, PMC9701605)", src=(os.path.join(PKG, "collection_fulltext.xml"), None, True),
      showed="Perspective describing the RRRM-2 design (two ages; ISS-terminal versus LAR with 24 d Earth recovery; matched ground controls) and the trade-offs of in-flight dissection versus LAR, including re-entry stress and circadian disruption for LAR. No data analysis of brain or IEGs.",
      quote="Using LAR, organisms have their circadian rhythms disrupted and experience the stressful conditions of atmospheric re-entry.",
      relation="Already names LAR-specific stressors and collection constraints for RRRM-2 (narrows E1 and C1). Does not quantify a ground-only protocol effect.",
      claims=["E1", "P1"]),
 dict(id="grandke2026", doi="10.1038/s41467-026-68737-1", pmid="41644562", year=2026,
      title="MiRNAs shape mouse age-independent tissue adaptation to spaceflight via ECM and developmental pathways",
      source_read="full text (Europe PMC fullTextXML PMC12876965)", src=("get_ft_grandke2026_*.xml", None),
      showed="RRRM-1 small-RNA atlas: 686 samples, 13 organs incl. brain, two ages, FL/HGC/VGC, LAR and TERM. Extraction cohort (TERM d21 vs LAR d43, preservation) contributed variance; analysis focused on LAR; LAR-TERM comparison framed as reentry stress versus exposure length. Brain had few FL-vs-HGC miRNAs and larger housing effects.",
      quote="This underscores the impact of the extraction timepoint (TERM: day 21, LAR: day 43) and the preservation method (see Methods) on miRNA expression.",
      relation="Closest multi-organ RRRM work: recognises a collection-cohort effect on miRNA but interprets LAR-TERM as reentry/exposure. Different mission (RRRM-1), assay (small RNA) and contrast; does not isolate a protocol effect in flight-naive animals. Narrows E1 wording ('collection cohort affects omics' is not new).",
      claims=["E1", "P1"]),
 dict(id="holley2022", doi="10.1038/s41526-022-00217-4", pmid="35948598", year=2022,
      title="Characterization of gene expression profiles in the mouse brain after 35 days of spaceflight mission",
      source_read="full text (package cache holley_fulltext.xml, PMC9365836)", src=(os.path.join(PKG, "holley_fulltext.xml"), None, True),
      showed="RR-9 male C57BL/6, NanoString neuroinflammation panel; flight mice CO2-euthanized 38 +/- 4 h after splashdown, ground controls by the same method; 54 DEGs. Authors propose on-orbit euthanasia with immediate preservation in future work.",
      quote="In order to evaluate the immune response to spaceflight environment, in the further study, mice will be euthanized in the orbit and brain tissues will be immediately preserved for analysis on the ground.",
      relation="Recognises that return/landing may confound brain expression and that collection should change (narrows E1 'first to recognise'). Different mission, assay and protocol; no ground-only protocol contrast.",
      claims=["E1"]),
 dict(id="kremsky2023", doi="10.3390/ijms241713569", pmid="37686374", year=2023,
      title="Spaceflight-Induced Gene Expression Profiles in the Mouse Brain Are Attenuated by Treatment with the Antioxidant BuOE",
      source_read="full text (Europe PMC fullTextXML PMC10487739)", src=("get_ft_kremsky2023_*.xml", None),
      showed="RR-18 LAR mice, ketamine/xylazine for flight and ground, flight mice processed within 20 h of splashdown; GeoMx spatial profiling of four regions; Arc is among the top five spaceflight DEGs by fold change in frontal cortex and cortex.",
      quote="within 20 h of splashdown. Mice were exsanguinated by closed-cardiac blood collection under deep Ketamine/Xylazine (150/45 mg/kg) anesthesia ... Arc , Ptgds , Mobp , Gm9936 , and Apod for CT",
      relation="Example of an IEG (Arc) reported as a top spaceflight effect when return timing is aliased with flight; supports the relevance of E1 (an IEG shift can arise from collection/return conditions) without showing it.",
      claims=["E1"]),
 dict(id="masarapu2024", doi="10.1038/s41467-024-48916-8", pmid="38862479", year=2024,
      title="Spatially resolved multiomics on the neuronal effects induced by spaceflight in mice",
      source_read="publisher full text (package cache masarapu_publisher.html)", src=(os.path.join(PKG, "masarapu_publisher.html"), None, True),
      showed="RR-3 female BALB/c brain spatial transcriptomics and snMultiome; flight and ground both ketamine/xylazine/acepromazine, carcasses -80 C, dissected on Earth; no IEG findings in the text.",
      quote="At the end of the mission, spaceflight and ground controls were euthanized (using euthanasia chemical ketamine/xylazine/acepromazine) and whole carcasses were stored at −80 °C until the dissection on earth.",
      relation="Protocol-matched flight-ground design in a different mission (RR-3/OSD-352). Not RRRM-2/RR-10/RR-23; no protocol contrast. Context for P1.",
      claims=["P1"]),
 dict(id="staib2014", doi="10.1089/neu.2013.3243", pmid="24945082", year=2014,
      title="Anesthesia for Euthanasia Influences mRNA Expression in Healthy Mice and after Traumatic Brain Injury",
      source_read="abstract (NCBI PMC front matter; full-text API returned HTTP 500 at Europe PMC)", src=("get_ncbi_staib2014_*.xml", None),
      showed="Hippocampal qPCR of 10 transcripts; 1 min isoflurane, 3 min midazolam/fentanyl/medetomidine, or 3 min chloral hydrate before sampling versus no anaesthesia; FosB lower with the injectable regimens in naive mice. Ketamine/xylazine not tested.",
      quote="The data demonstrate that a few minutes of anesthesia before tissue sampling are sufficient to induce immediate mRNA changes, which seem to predominate in the early-regulated gene cluster.",
      relation="Directly overturns any claim that euthanasia-anaesthesia effects on brain early-response genes are new. Bounds E1: our contribution is genome-scale, spaceflight-design-specific quantification, not the phenomenon.",
      claims=["E1"]),
 dict(id="bunting2015", doi="10.3389/fnbeh.2015.00363", pmid="26793081", year=2015,
      title="Influence of Isoflurane on Immediate-Early Gene Expression",
      source_read="abstract (Europe PMC)", src=("epmc_abs_bunting2015_*.json", "abstract_raw"),
      showed="Rat catFISH (Arc, Zif268/Egr1) in hippocampus and cortex: anaesthetic-level isoflurane prevents initiation of IEG transcription but does not stop transcription already started.",
      quote="General anesthesia with isoflurane prevents initiation of learning-related transcription but does not stop ongoing transcription of two plasticity-related IEGs, Arc and Zif268",
      relation="Mechanistic plausibility that inhaled anaesthesia (LAR arm) suppresses new IEG induction; supports direction of E1 but does not test it.",
      claims=["E1"]),
 dict(id="kim2022", doi="10.3389/fnmol.2022.892345", pmid="35706427", year=2022,
      title="Transcriptional Activation, Deactivation and Rebound Patterns in Cortex, Hippocampus and Amygdala in Response to Ketamine Infusion in Rats",
      source_read="abstract (Europe PMC)", src=("epmc_abs_ketamine_ieg_2022_*.json", "abstract_raw"),
      showed="RNA-seq after 1 h or 10 h of anaesthetic-dose ketamine infusion in rat cortex, hippocampus and amygdala: an IEG network (e.g. c-Fos) is induced within 1 h in all three regions.",
      quote="Induction of a shared transcriptional regulatory network occurred within 1 h in all three brain regions consisting of (a) genes involved in stimulus-transcription factor coupling that are induced during altered synaptic activity (immediate early genes, IEGs, such as c-Fos",
      relation="Plausibility that ketamine (ISS-T arm) raises IEG mRNA; timescale (1 h infusion) longer than a euthanasia overdose, so it does not establish causation in RRRM-2.",
      claims=["E1"]),
 dict(id="hu2024", doi="10.7554/eLife.88420", pmid="38512722", year=2024,
      title="Comparative brain-wide mapping of ketamine- and isoflurane-activated nuclei and functional networks in the mouse brain",
      source_read="abstract (Europe PMC)", src=("epmc_abs_hu2024_*.json", "abstract_raw"),
      showed="Brain-wide c-Fos mapping in mice: ketamine activates widespread cortical and subcortical regions; isoflurane mainly hypothalamic regions.",
      quote="Our findings reveal that KET activates a wide range of brain regions, notably in the cortical and subcortical nuclei involved in sensory, motor, emotional, and reward processing",
      relation="Agent-specific IEG activation patterns; consistent with the ISS-T > LAR sign but measured as protein after anaesthesia, not euthanasia.",
      claims=["E1"]),
 dict(id="tyssowski2018", doi="10.1016/j.neuron.2018.04.001", pmid="29681534", year=2018,
      title="Different Neuronal Activity Patterns Induce Different Gene Expression Programs",
      source_read="abstract (Europe PMC)", src=("epmc_abs_tyssowski2018_try_*.json", "abstract_raw"),
      showed="Brief neuronal activity induces only the first temporal wave of activity-regulated genes (the rapid primary response genes), MAPK/ERK-dependent and translation-independent.",
      quote="Brief activity selectively induces a small subset of the activity-regulated gene program that corresponds to the first of three temporal waves of genes induced by sustained activity.",
      relation="Origin of the rPRG panel; supports reading the signal as a minutes-scale programme, which fits a peri-collection window. Does not concern spaceflight.",
      claims=["E1"]),
 dict(id="marsh2022", doi="10.1038/s41593-022-01022-8", pmid="35260865", year=2022,
      title="Dissection of artifactual and confounding glial signatures by single-cell sequencing of mouse and human brain",
      source_read="abstract (Europe PMC) + full text (package cache marsh_fulltext.xml)", src=("abstract_records.json", "marsh2022"),
      showed="Enzymatic dissociation of brain induces an ex vivo signature (IEGs Fos/Jun, Hspa1a, Dusp1, Ccl3/Ccl4), most prominent in microglia; a validated inhibitor protocol prevents it; a similar signature appears in post-mortem human snRNA-seq.",
      quote="Here, we demonstrate that use of enzymatic dissociation on brain tissue induces an aberrant ex vivo gene expression signature, most prominently in microglia, which is prevalent in published literature and can substantially confound downstream analyses.",
      relation="Processing-induced IEG/stress signatures in brain are established; E1 must not claim the class of artefact is new. Different mechanism (dissociation, single-cell).",
      claims=["E1"]),
 dict(id="vandenbrink2017", doi="10.1038/nmeth.4437", pmid="28960196", year=2017,
      title="Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations",
      source_read="publisher page (Correspondence text and Supplementary Fig. 3 legend; no abstract exists)", src=("get_nature_vandenbrink2017_*.html", None),
      showed="Muscle stem cells: tissue dissociation/FACS induces an immediate-early/stress gene subpopulation; smFISH shows Fos absent in intact muscle and induced during isolation.",
      quote="we found a subpopulation that is strongly affected by the widely used dissociation ... Fos expression is induced during the SC isolation procedure",
      relation="Establishes handling-induced IEG expression (muscle). Bounds E1 as a brain/spaceflight instance of a known artefact class.",
      claims=["E1"]),
 dict(id="dachet2021", doi="10.1038/s41598-021-85801-6", pmid="33758256", year=2021,
      title="Selective time-dependent changes in activity and cell-specific gene expression in human postmortem brain",
      source_read="abstract (Europe PMC) + full text (package cache)", src=("abstract_records.json", "dachet2021"),
      showed="Fresh human neocortex held 0-24 h ex vivo: activity-dependent neuronal transcripts fall within hours while housekeeping genes are preserved; glial transcripts rise.",
      quote="Within a few hours, a selective reduction in the number of neuronal activity-dependent transcripts occurred with relative preservation of housekeeping genes commonly used as a reference for RNA normalization.",
      relation="Collection-interval sensitivity of activity genes in brain is known; supports E1 plausibility (and warns that RIN/housekeeping QC will not detect it).",
      claims=["E1"]),
 dict(id="chakkingal2023", doi="10.1186/s12864-023-09794-4", pmid="38093185", year=2023,
      title="The forgotten variable? Does the euthanasia method and sample storage condition influence an organisms transcriptome - a gene expression analysis on multiple tissues in pigs",
      source_read="abstract (Europe PMC)", src=("epmc_abs_chakkingal2023_*.json", "abstract_raw"),
      showed="Piglets: nitrogen anoxia versus T-61 euthanasia gave no protein-coding DEGs in pituitary, hypothalamus, liver or lung (3' mRNA-seq); storage (RNAlater vs LN2) did change expression.",
      quote="On comparison of the 3'mRNA-Seq data generated from pituitary, hypothalamus, liver and lung tissues, no significant differential expression in the protein coding genes were detected between the euthanasia methods.",
      relation="Counter-evidence: euthanasia agent need not alter brain-region transcriptomes. Forces E1 to stay agnostic about the agent and to call the variable 'collection protocol'.",
      claims=["E1"]),
 dict(id="sakilam2026", doi="10.20944/preprints202609.2442.v1", pmid=None, year=2026,
      title="Spaceflight Is Associated with an Age-Conserved Kidney Transcriptomic Signature Linked to Inflammation, Tubular Channelopathy, and Senescence (preprint, not peer reviewed)",
      source_read="abstract (Europe PMC preprint record)", src=("epmc_abs_sakilam2026_*.json", "abstract_raw"),
      showed="RRRM-2 kidney bulk RNA-seq (OSD-771): flight (55-58 d, killed on ISS) versus basal controls, two ages; 7-gene age-conserved signature incl. Cry1 down.",
      quote="We profiled kidney transcriptomes from NASA GeneLab dataset OSD-771/GLDS-674 (Rodent Research Reference Mission-2)",
      relation="Recent RRRM-2 analysis of a different tissue whose flight-versus-basal contrast differs in collection site/protocol; shows the positioning need for E1. Not brain; does not bear on P1 for our deposits.",
      claims=["E1", "P1"]),
 dict(id="szyld2006", doi="10.1007/s11075-006-9046-2", pmid=None, year=2006,
      title="The many proofs of an identity on the norm of oblique projections", source_read="abstract (Springer publisher page)",
      src=("get_springer_szyld2006_*.html", None),
      showed="Reviews the classical identity ||P|| = ||I - P|| for any idempotent P other than 0 and I, noting it has been proven repeatedly.",
      quote="Given an oblique projector P on a Hilbert space, i.e., an operator satisfying P 2=P, which is neither null nor the identity, it holds that ||P|| = ||I –P||.",
      relation="Kills any claim that the oblique-projector form or its spectral gain is a new theorem (T1). The gain of score removal I - g w^T equals ||g||*||w|| by this identity.",
      claims=["T1", "C1"]),
 dict(id="behrens1994", doi="10.1109/78.286957", pmid=None, year=1994,
      title="Signal processing applications of oblique projection operators", source_read="Crossref metadata only (publisher abstract not retrieved: IEEE returned HTTP 202)",
      src=None, showed="NOT READ beyond bibliographic metadata.", quote=None,
      relation="Cite only as the standard reference for oblique projection in structured-interference removal; do not attribute specific content until read.",
      claims=["T1"]),
 dict(id="gagnonbartsch2012", doi="10.1093/biostatistics/kxr034", pmid="22101192", year=2012,
      title="Using control genes to correct for unwanted variation in microarray data", source_read="abstract (Europe PMC)",
      src=("abstract_records.json", "gb2012"),
      showed="RUV-2: restrict factor analysis to negative-control genes to estimate unwanted variation, because separating unwanted from biological variation is the central difficulty.",
      quote="A central problem with this approach is the difficulty in discerning the unwanted variation from the biological variation that is of interest to the researcher.",
      relation="Using a gene subset as a proxy for unwanted variation is established (C1 'sensor' = RUV control-gene idea, with positive rather than negative controls).",
      claims=["T1", "C1"]),
 dict(id="wang2017", doi="10.1214/16-AOS1511", pmid="31439967", year=2017,
      title="Confounder adjustment in multiple hypothesis testing", source_read="abstract (Europe PMC)",
      src=("abstract_records.json", "wang2017"),
      showed="Unifies RUV-4/LEAPP-type confounder adjustment; identification requires either negative controls or sparsity of true effects; strong confounders allow near-oracle power.",
      quote="the first requires a set of \"negative controls\" that are known a priori to follow the null distribution; the second requires the true nonnulls to be sparse.",
      relation="Theory already states the identification conditions a genome-wide correction needs; T1 must be positioned as an empirical, small-n precision result, not new identification theory.",
      claims=["T1", "M1"]),
 dict(id="gerard2021", doi="10.5705/ss.202018.0345", pmid="38148787", year=2021,
      title="Unifying and generalizing methods for removing unwanted variation based on negative controls", source_read="abstract (Europe PMC)",
      src=("abstract_records.json", "gerard2021"),
      showed="RUV* framework unifying RUV variants; RUVB competitive, but reliable calibration across data sets remains difficult.",
      quote="However, providing a consistently reliable calibration among the data sets remains challenging.",
      relation="Calibration failure of control-based correction is known; supports framing gain reporting as a guardrail, not a fix.",
      claims=["T1", "C1"]),
 dict(id="miao2018", doi="10.1093/biomet/asy038", pmid="33343006", year=2018,
      title="Identifying causal effects with proxy variables of an unmeasured confounder", source_read="abstract (Europe PMC)",
      src=("abstract_records_methods.json", "miao2018"),
      showed="With two independent proxies meeting a rank condition the causal effect is identified; with one proxy only a test of the null is available.",
      quote="When only one proxy for the confounder is available, or the required rank condition is not met, we develop a strategy to test the null hypothesis of no causal effect.",
      relation="Formal version of T1's intuition: a single reliable nuisance proxy does not identify the adjustment. Narrows T1 to an empirical illustration.",
      claims=["T1"]),
 dict(id="nygaard2016", doi="10.1093/biostatistics/kxv027", pmid="26272994", year=2016,
      title="Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses", source_read="abstract (Europe PMC)",
      src=("abstract_records.json", "nygaard2016"),
      showed="Two-way ANOVA batch adjustment with unbalanced designs induces incorrect group differences and overconfident downstream inference.",
      quote="Unfortunately, this approach may systematically induce incorrect group differences in downstream analyses when groups are distributed between the batches in an unbalanced manner.",
      relation="Post hoc correction under confounded design is known to mislead; T1/C1 add a gain quantity, not the warning.",
      claims=["T1", "C1", "M1"]),
 dict(id="zindler2020", doi="10.1186/s12859-020-03559-6", pmid="32605541", year=2020,
      title="Simulating ComBat: how batch correction can lead to the systematic introduction of false positive results in DNA methylation microarray studies", source_read="abstract (Europe PMC)",
      src=("abstract_records.json", "zindler2020"),
      showed="ComBat on random simulated DNAm data produced many FDR- and Bonferroni-significant false positives, in unbalanced and balanced designs.",
      quote="Using ComBat to correct for batch effects in randomly generated samples produced alarming numbers of false discovery rate (FDR) and Bonferroni-corrected (BF) false positive results",
      relation="Correction-induced false positives are established; bounds C1.",
      claims=["C1"]),
 dict(id="ye2023", doi="10.1038/s43588-023-00500-8", pmid="38177326", year=2023,
      title="Batch-effect correction with sample remeasurement in highly confounded case-control studies", source_read="abstract (Europe PMC)",
      src=("abstract_records.json", "ye2023"),
      showed="Design-based correction: remeasuring a subset of samples in each batch, with theory and a power calculator; few remeasured samples suffice when between-batch correlation is high.",
      quote="One approach to address the batch effects is repeatedly measuring a subset of samples in each batch.",
      relation="Design-level remedies for confounded batches exist; the paper's 'control by design' recommendation is consistent with, not new relative to, this literature. The spaceflight analogue (both protocols in flight-naive animals) is the specific contribution.",
      claims=["M1", "C1"]),
 dict(id="gerstner2016", doi="10.1186/s12864-016-3065-8", pmid="27801296", year=2016,
      title="Removal of unwanted variation reveals novel patterns of gene expression linked to sleep homeostasis in murine cortex", source_read="abstract (Europe PMC)",
      src=("abstract_records.json", "gerstner2016"),
      showed="RUV normalisation of mouse cortex sleep-deprivation and recovery-sleep expression improves differential-expression detection.",
      quote="RUV improves detection of differential gene expression following sleep deprivation.",
      relation="Control-based correction has been applied to mouse brain activity-state transcriptomes; bounds M1/C1 in brain.",
      claims=["M1", "C1"]),
 dict(id="varma2006", doi="10.1186/1471-2105-7-91", pmid="16504092", year=2006,
      title="Bias in error estimation when using cross-validation for model selection", source_read="abstract (Europe PMC)",
      src=("abstract_records_methods.json", "varma2006"),
      showed="Tuning a classifier by CV and reporting that CV error is biased; nested CV gives an almost unbiased estimate.",
      quote="A nested CV procedure provides an almost unbiased estimate of the true error.",
      relation="Nested selection is standard (M1 is not methodologically new in this respect).",
      claims=["M1"]),
 dict(id="kapoor2023", doi="10.1016/j.patter.2023.100804", pmid="37720327", year=2023,
      title="Leakage and the reproducibility crisis in machine-learning-based science", source_read="abstract (Europe PMC)",
      src=("abstract_records_methods.json", "kapoor2023"),
      showed="Survey finds leakage in 17 fields (294 papers) and gives a taxonomy of eight leakage types with model info sheets.",
      quote="we find 17 fields where leakage has been found, collectively affecting 294 papers",
      relation="Training-only preprocessing to prevent leakage is established practice (M1).",
      claims=["M1"]),
 dict(id="gouk2021", doi="10.1007/s10994-020-05929-w", pmid=None, year=2021,
      title="Regularisation of neural networks by enforcing Lipschitz continuity", source_read="abstract (Springer publisher page)",
      src=("get_springer_gouk2021_*.html", None),
      showed="Upper-bounds a network's Lipschitz constant and trains under that bound by projected stochastic gradient methods.",
      quote="Our technique is then used to formulate training a neural network with a bounded Lipschitz constant as a constrained optimisation problem that can be solved using projected stochastic gradient methods.",
      relation="Gain/Lipschitz-constrained fitting by projection is established in ML (C1/M1 import, not invention).",
      claims=["C1", "M1"]),
 dict(id="zames1966", doi="10.1109/TAC.1966.1098316", pmid=None, year=1966,
      title="On the input-output stability of time-varying nonlinear feedback systems Part one: Conditions derived using concepts of loop gain, conicity, and positivity",
      source_read="Crossref metadata only", src=None, showed="NOT READ beyond bibliographic metadata (classical small-gain reference).", quote=None,
      relation="Cite only as the origin of loop-gain stability conditions for the C1 framing.",
      claims=["C1"]),
 dict(id="rodriguez2026", doi="10.1186/s12974-026-03913-0", pmid=None, year=2026,
      title="Spaceflight suppresses inflammasome signaling and alters dopaminergic regulation in the mouse brain from Rodent Research-1 (RR-1) mission",
      source_read="Crossref metadata only (published 4 Aug 2026; not yet indexed in Europe PMC on 15 Oct 2026)", src=None,
      showed="Title only: RR-1 brain, inflammasome and dopaminergic regulation.", quote=None,
      relation="Recent brain spaceflight paper on a different mission; shows Europe PMC indexing lag, a limit on the P1 absence claim.",
      claims=["P1"]),
]

out, bad = [], []
for e in E:
    rec = {k: e[k] for k in ("id", "doi", "pmid", "title", "year", "source_read")}
    rec["what_it_showed"] = e["showed"]
    rec["quoted_supporting_sentence"] = e["quote"]
    rec["relation_to_our_claim"] = e["relation"]
    rec["claims"] = e["claims"]
    rec["words_in_what_it_showed"] = len(e["showed"].split())
    if e["quote"] and e["src"]:
        spec = e["src"]
        path = spec[0] if (len(spec) > 2 and spec[2]) else (g(spec[0]) or os.path.join(C, spec[0]))
        key = spec[1]
        if key == "abstract_raw":
            d = json.load(open(path))
            rs = d.get("resultList", {}).get("result", [])
            txt = html.unescape(re.sub(r"<[^>]+>", " ", rs[0].get("abstractText", "") if rs else ""))
        elif key:
            txt = html.unescape(re.sub(r"<[^>]+>", " ", json.load(open(path))[key]["abstract"]))
        else:
            txt = text_of(path)  # already tag-stripped and unescaped; do not strip again
        src_norm = norm(txt)
        frags = [f for f in e["quote"].split("...") if f.strip()]
        ok = all(norm(f) in src_norm for f in frags)
        rec["quote_verified_against"] = os.path.relpath(path, ROOT) if path.startswith(ROOT) else path
        rec["quote_verified"] = ok
        if not ok:
            bad.append(e["id"])
    else:
        rec["quote_verified"] = None
    out.append(rec)

with open(os.path.join(ROOT, "B_closest_prior_art.json"), "w") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print(len(out), "entries; quote failures:", bad)
print("over-60-word summaries:", [r["id"] for r in out if r["words_in_what_it_showed"] > 60])
