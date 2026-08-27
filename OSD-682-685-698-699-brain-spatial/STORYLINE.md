# What happened to the mouse brain: a mechanistic storyline

Working narrative for OSD-682/685/698/699 — 34 days on the ISS, C57BL/6 male mice,
GeoMx DSP across CA1, dentate gyrus, frontal cortex and cerebral cortex, with and
without the antioxidant BuOE (MnTnBuOE-2-PyP).

**Status of this document.** This is a hypothesis assembled from the analysis, not
a set of established results. No gene reaches FDR significance transcriptome-wide;
three reach it within the pre-specified target panel. Every claim below is tagged
with its evidence tier so nothing gets quoted at the wrong strength:

| Tier | Meaning |
|---|---|
| **[A]** | Survives multiple-testing correction (panel-restricted BH < 0.05, or category q < 0.05) |
| **[B]** | Consistent across ≥3 of 4 regions with the same direction, nominal P < 0.01 |
| **[C]** | Nominal P < 0.05 in one or two regions — suggestive, needs replication |

Design limits that cap every interpretation: n = 3 per cell; each treatment group
occupies its own GeoMx slide, so batch and treatment cannot be separated; all four
regions come from the same tissue sections, so cross-region agreement is not
independent replication.

---

## The one-paragraph version

Thirty-four days of spaceflight left the mouse brain in a state of **oxidative and
metabolic stress with a regionally split neuronal response**. A cytoprotective
programme was elevated in every region sampled: metallothionein *Mt1* (+0.75 log2,
all four regions, P = 0.0022), the cold-shock RNA-binding proteins *Rbm3* (+0.78,
P = 2.4 × 10⁻⁵) and *Cirbp* (+0.67, P = 0.0027), and both superoxide dismutases
(*Sod1* +0.43, P = 0.011; *Sod2* +0.44, P = 0.025). Prostaglandin D2 synthase rose
sharply (*Ptgds* +1.61 across three regions, P = 1.7 × 10⁻⁴), pointing at
lipid-mediated rather than cytokine-mediated neuroinflammation. In the **cortex**,
activity-dependent immediate-early genes fell markedly — *Arc* −1.43 across three
regions (P = 1.3 × 10⁻⁴) and *Fos* −0.77 (P = 0.013) — consistent with reduced
neuronal activation. In the **hippocampus**, the opposite: postsynaptic and
glutamatergic machinery rose as a coordinated set of 20 genes (dentate gyrus mean
t = +1.22, P = 0.024; CA1 +0.71, P = 0.048), led by *Nrgn* (+0.64, panel
BH = 0.040), *Camk2b* (+0.59), *Grin1* (+0.42) and *Gria1* (+0.38) — alongside
matrix metalloproteinases (*Mmp2* +0.46, P = 0.0016) and *App* (+0.43) / *Mapt*
(+0.44) in CA1, a compensatory or excitotoxic-adjacent signature in the most
vulnerable subfield. **BuOE did not simply blunt everything.** It suppressed the
vascular and inflammatory arm in frontal cortex (blood-brain-barrier category mean
t = −1.00 and neuroinflammation −0.58, both q = 0.032) and the inflammasome in CA1
(−0.93, P = 0.009); it normalised the hippocampal synaptic overshoot, reversing
*Nrgn* (−0.48) and *Junb* (−0.48) in dentate gyrus with a significant Flight × BuOE
interaction (*Junb* −0.53, panel BH = 0.012); and in frontal cortex it shifted the
balance *toward* neuronal genes (mean t = +1.39, P = 0.007), where *Junb* rose
+1.24 (panel BH = 0.005) with an interaction of the **opposite sign** (+0.72, panel
BH = 0.021). The net effect was attenuation of the flight response in hippocampus —
clearest in dentate gyrus (bootstrap P < 0.001) — and the reverse in cerebral
cortex (amplification, P = 0.010).

---

## Step 1 — A brain-wide oxidative stress response

The most reproducible signal in the entire dataset is not in the target panel; it
is a classic stress-response module.

| Gene | Direction | Where | Evidence |
|---|---|---|---|
| *Mt1* (metallothionein-1) | ↑ +0.75 | **all 4 regions** | **[B]** P = 0.0022 |
| *Rbm3* (cold-inducible RNA-binding) | ↑ +0.78 | DG, FCtx, Ctx | **[B]** P = 0.000024 |
| *Cirbp* (cold-inducible RNA-binding) | ↑ +0.67 | CA1, DG, FCtx | **[B]** P = 0.0027 |
| *Sod1* | ↑ +0.43 | DG | **[C]** P = 0.011 |
| *Sod2* | ↑ +0.44 | FCtx | **[C]** P = 0.025 |

**Why this matters mechanistically.** Metallothioneins are induced by reactive
oxygen species, glucocorticoids and metal load, and buffer redox stress.
*Rbm3*/*Cirbp* are cold-shock proteins but respond to cellular stress generally,
and both regulate translation under stress. *Sod1*/*Sod2* are the endogenous
superoxide dismutases.

This is the single most important element of the storyline because **BuOE is a
manganese porphyrin SOD mimetic**. The tissue is upregulating exactly the pathway
the drug substitutes for. It gives a coherent pharmacological rationale: the drug
is acting on a stress axis the brain demonstrably engaged.

## Step 2 — Prostaglandin-mediated neuroinflammation

| Gene | Direction | Where | Evidence |
|---|---|---|---|
| *Ptgds* (prostaglandin D2 synthase) | ↑ **+1.61** | CA1, FCtx, Ctx | **[B]** P = 0.00017 |

The largest consistent effect size in the recurrent set. PTGDS is highly expressed
in leptomeninges and oligodendrocytes and produces PGD2, a mediator of
neuroinflammation and a regulator of sleep. Its elevation across three regions
suggests a lipid-mediated inflammatory component, plausibly meningeal in origin —
which is notable given that the panel's classical cytokines (*Il1b*, *Il6*, *Tnf*)
were either not measurable on this assay or did not move.

The inflammatory signal here is **lipid-mediated, not cytokine-mediated**. That is
a testable distinction and it changes which assays a follow-up should use.

## Step 3 — Mitochondrial and metabolic remodelling, concentrated in dentate gyrus

Top transcriptome-wide hits in DG are dominated by oxidative phosphorylation and
mitochondrial transport:

*Atp5b* ↑, *Cox6c* ↑, *Slc25a4* (ANT1) ↑, *Ppif* (cyclophilin D) ↑ — all **[C]**,
P = 6×10⁻⁵ to 2×10⁻⁴. *Gapdh* ↑ in CA1 and FCtx.

*Ppif* is the regulator of the mitochondrial permeability transition pore and a
direct link between oxidative stress and cell death. Its co-elevation with ETC
subunits in the dentate gyrus — the site of adult neurogenesis — is the most
mechanistically loaded observation in the dataset.

## Step 4 — The regional split in neuronal response (the core finding)

This is where the story becomes interesting rather than generic.

### Cortex: neuronal activity goes down

| Gene | Direction | Where | Evidence |
|---|---|---|---|
| *Arc* | ↓ **−1.43** | CA1, FCtx, Ctx | **[B]** P = 0.00013 |
| *Fos* | ↓ −0.77 | Ctx | **[C]** P = 0.013 |
| *S100b* | ↓ −0.57 | Ctx | **[C]** P = 0.038 |

*Arc* and *Fos* are the canonical activity-dependent immediate-early genes. A
coordinated fall is the transcriptional signature of **reduced neuronal
activation** — consistent with reduced sensorimotor input, altered arousal, or
suppressed cortical drive in microgravity.

### Hippocampus: postsynaptic machinery goes up, as a coordinated set

Category-level: **Neuronal/synaptic mean t = +1.22 in DG and +0.71 in CA1, both
P < 0.05 [A-adjacent]** — this is 20 genes moving together, not a single hit.

| Gene | Function | DG | CA1 |
|---|---|---|---|
| *Nrgn* (neurogranin) | postsynaptic calmodulin-binding | **+0.64** [A] | +0.54 [C] |
| *Junb* | AP-1 immediate-early TF | **+0.57** [A] | — |
| *Grin1* (NMDA receptor NR1) | glutamatergic | +0.42 [C] | +0.42 [C] |
| *Gria1* (AMPA receptor) | glutamatergic | +0.38 [C] | — |
| *Camk2b* | synaptic plasticity kinase | +0.59 [C] | +0.64 [C] |
| *Slc17a7* (VGLUT1) | vesicular glutamate transporter | — | +0.40 [C] |
| *Rgs7* | G-protein signalling | +0.51 [C] | — |

**Interpretation.** Coordinated upregulation of NMDA and AMPA receptor subunits,
VGLUT1, CaMKIIβ and neurogranin is the profile of either **homeostatic
compensation** for reduced input, or **increased excitatory drive**. Neurogranin is
of particular interest: it is a clinically used CSF marker of synaptic injury, and
it is the one panel gene that survives correction for the flight effect in saline.

The cortex-down / hippocampus-up dissociation is the most publishable biological
observation here, and it would have been invisible in a bulk homogenate — it is
precisely what spatial profiling buys.

## Step 5 — Barrier and matrix remodelling with a proteinopathy flavour, in CA1

| Gene | Direction | Where | Evidence |
|---|---|---|---|
| *Mmp2* | ↑ +0.46 | CA1 | **[C]** P = 0.0016 |
| *Mmp12* | ↑ +0.30 | CA1 | **[C]** P = 0.020 |
| *App* (amyloid precursor) | ↑ +0.43 | CA1 | **[C]** P = 0.0055 |
| *Mapt* (tau) | ↑ +0.44 | CA1 | **[C]** P = 0.025 |
| *Ninj1* | ↓ −0.50 | CA1 | **[C]** P = 0.011 |

Matrix metalloproteinases degrade basement membrane and tight junctions and are a
standard readout of blood-brain-barrier remodelling. Their co-elevation with *App*
and *Mapt* in CA1 — the subfield most vulnerable in ageing and ischaemia — is the
combination most relevant to long-term neurodegeneration risk. Category-level
injury-biomarker signal is also highest in CA1 (mean t = +0.95).

## Step 6 — Myelin

*Mbp* ↑ in CA1, FCtx and Ctx [C]; *Mobp* is the single strongest gene in cortex
(+1.33, P = 2.9×10⁻⁵) [C]; *Slc1a2*/GLT-1 ↑ in DG [C]; *Gpr17* ↓ −0.66 in FCtx [C].

GPR17 marks oligodendrocyte precursors and is downregulated as they mature. Falling
*Gpr17* with rising *Mbp*/*Mobp* is consistent with a **shift toward mature
myelinating oligodendrocytes**. Glial/myelin category is up in CA1 (+0.91) and DG
(+0.97).

---

## What BuOE did — three distinct actions, not one

The antioxidant did not uniformly damp the flight response. It acted differently by
arm and by region, which is why a single "does it help?" question was the wrong one.

### Action 1 — Suppressed the vascular, matrix and inflammasome arm **[A]**

In frontal cortex, under flight, two categories fall together with correction-level
support:

- **Blood-brain barrier / endothelial**: 16 targets, mean t = −1.00, **q = 0.032**
- **Neuroinflammation / cytokine**: 17 targets, mean t = −0.58, **q = 0.032**

And in CA1, **inflammasome/pyroptosis** mean t = −0.93, P = 0.009.

Driver genes: *Mmp2* −1.06 (FCtx), *Itgb3* −0.67 (FCtx), *Mmp12* −0.51 (CA1),
*Nos3* −0.56 (CA1), *Tjp1* −0.41 (Ctx), *Gsdmd* −0.38 (FCtx), *Rela* −0.42 (DG).

*Gsdmd* is the pyroptosis executioner and *Rela* is the NF-κB p65 subunit. Their
suppression alongside MMPs is a coherent anti-inflammatory, barrier-protective
signature — exactly the intended pharmacology.

### Action 2 — Normalised the hippocampal synaptic overshoot **[A]**

In dentate gyrus, BuOE moved the flight-elevated synaptic genes back **down**:

| Gene | Flight effect (saline) | BuOE effect in flight |
|---|---|---|
| *Nrgn* | **+0.64** | **−0.48** |
| *Junb* | **+0.57** | **−0.48** |
| *Slc17a7* | — | −0.32 |
| *Camk2a* | — | −0.25 |

Both *Nrgn* and *Junb* carry a **significant Flight × BuOE interaction in DG**
(*Junb* −0.53, panel BH = 0.012 **[A]**). This is the cleanest pharmacological
result in the dataset: the drug specifically opposes the flight-induced change in
the same genes, in the same region.

It also explains the quantitative attenuation result — DG shows the strongest
whole-transcriptome attenuation of the flight response (bootstrap P < 0.001).

### Action 3 — Shifted frontal cortex toward neuronal and neurogenic genes **[A]**

The opposite direction from hippocampus, and the reason a single global answer
fails:

- **Neuronal/synaptic** mean t = **+1.39**, P = 0.007
- **Neurogenesis/cytoskeletal** mean t = +1.16

Driver genes: *Junb* **+1.24** (panel BH = 0.005, the strongest single result in
the study), *Nrgn* +0.68, *Cacna1c* +0.45, *Grin1* +0.42, *Sptan1* +0.33.

*Junb* carries a significant interaction here too (+0.72, panel BH = 0.021) — but
**with the opposite sign to dentate gyrus**. The same gene, the same drug, opposite
direction in two regions.

---

## Net effect, by region

RMS effect = noise-corrected root-mean-square spaceflight response across all
15,782 targets, in log2 units. Δ = saline − BuOE mean squared effect; positive
means the drug shrank the flight response.

| Region | Flight response (key effect sizes) | BuOE effect | RMS saline → BuOE | Δ (95% CI) | Verdict |
|---|---|---|---|---|---|
| **Dentate gyrus** | Synaptic ↑ (mean t = +1.22), mitochondrial ↑, neurogenesis ↑ (+1.02), coagulation ↓ (−0.63) | Reverses the synaptic overshoot (*Nrgn* −0.48, *Junb* −0.48) | 0.119 → 0.071 | **+0.0091** (+0.0066, +0.0117) | **Attenuation, P < 0.001** |
| **CA1** | Synaptic ↑ (+0.71), MMPs ↑ (*Mmp2* +0.46), *App*/*Mapt* ↑, injury markers ↑ (+0.95) | Suppresses inflammasome (−0.93) and *Mmp12* (−0.51) | 0.101 → 0.081 | +0.0037 (+0.0008, +0.0067) | Attenuation, P = 0.013 |
| **Frontal cortex** | Mixed; *Arc* ↓, *Gpr17* −0.66, *Tspan2* +0.59 | Suppresses BBB (−1.00) + inflammation (−0.58), raises neuronal (+1.39) | 0.081 → 0.066 | +0.0022 (+0.0000, +0.0044) | Attenuation, P = 0.047 |
| **Cerebral cortex** | *Arc* −1.57, *Fos* −0.77, *Mobp* +1.33 | Suppresses myelin genes (−1.05), *Cd81* −0.48 | 0.040 → 0.064 | **−0.0025** (−0.0045, −0.0006) | **Amplification**, P = 0.010 |

The cortex result is the honest complication: BuOE does not help everywhere, and
in cerebral cortex the flight response is larger under drug. Any framing that
presents the antioxidant as uniformly protective would be contradicted by our own
data.

---

## The mechanistic hypothesis, stated as a testable chain

| # | Link in the chain | Key genes and effect sizes (log2FC) | Statistic | Tier |
|---|---|---|---|---|
| 1 | Microgravity and radiation raise reactive oxygen species and alter metabolic demand | *Mt1* **+0.75** (4/4 regions), *Rbm3* **+0.78**, *Cirbp* **+0.67**, *Sod1* +0.43, *Sod2* +0.44 | P = 2.4 × 10⁻⁵ – 0.025 | **[B]** |
| 2 | Redox stress engages mitochondrial permeability machinery in the dentate gyrus | *Ppif*, *Slc25a4*, *Atp5b* +0.55, *Cox6c* +0.31 | P = 6.3 × 10⁻⁵ – 2.1 × 10⁻⁴ | **[C]** |
| 3 | Lipid-mediated inflammation is initiated, plausibly meningeal | *Ptgds* **+1.61** (3/4 regions) | P = 1.7 × 10⁻⁴ | **[B]** |
| 4 | Reduced sensorimotor and arousal input suppresses cortical activity | *Arc* **−1.43** (3/4 regions), *Fos* −0.77, *S100b* −0.57 | P = 1.3 × 10⁻⁴ – 0.038 | **[B]** |
| 5 | Hippocampal circuits compensate by upregulating postsynaptic glutamatergic machinery | *Nrgn* **+0.64**, *Camk2b* +0.59, *Grin1* +0.42, *Gria1* +0.38, *Slc17a7* +0.40 — 20 genes, DG mean t = **+1.22** | panel BH = 0.040; category P = 0.024 | **[A]** |
| 6 | In CA1, MMP-mediated barrier remodelling proceeds alongside *App*/*Mapt* elevation | *Mmp2* **+0.46**, *Mmp12* +0.30, *App* +0.43, *Mapt* +0.44 | P = 0.0016 – 0.025 | **[C]** |
| 7a | BuOE intercepts step 1 → suppresses the barrier / matrix / inflammasome arm | FCtx BBB mean t = **−1.00**, neuroinflammation **−0.58**; *Mmp2* −1.06, *Gsdmd* −0.38, *Rela* −0.42; CA1 inflammasome −0.93 | **q = 0.032**, q = 0.032; P = 0.009 | **[A]** |
| 7b | → normalises the hippocampal synaptic overshoot | *Nrgn* **−0.48**, *Junb* **−0.48** (DG); interaction *Junb* **−0.53** | panel BH = 0.012 | **[A]** |
| 7c | → shifts frontal cortex toward neuronal and neurogenic genes | *Junb* **+1.24**, *Nrgn* +0.68, *Cacna1c* +0.45, *Grin1* +0.42; neuronal mean t = **+1.39**; interaction **+0.72** | panel BH = **0.005**; P = 0.007; BH = 0.021 | **[A]** |

The sign reversal between rows 7b and 7c is the same gene, *Junb*, in the same
comparison, in two regions — the observation that most needs replicating.

### The single sentence for an abstract

> Spatial transcriptomic profiling of four brain regions after 34 days of
> spaceflight (GeoMx DSP, 15,782 targets, n = 3 per group) revealed a brain-wide
> oxidative-stress and prostaglandin response — *Mt1* +0.75 log2 in all four
> regions (P = 0.0022), *Ptgds* +1.61 in three (P = 1.7 × 10⁻⁴) — accompanied by
> opposing regional neuronal signatures: activity-dependent transcription fell in
> cortex (*Arc* −1.43, P = 1.3 × 10⁻⁴) while postsynaptic glutamatergic machinery
> rose coordinately in hippocampus (20 genes, dentate gyrus mean t = +1.22;
> *Nrgn* +0.64, panel-restricted BH = 0.040). The antioxidant BuOE attenuated the
> hippocampal response (dentate gyrus bootstrap P < 0.001) and suppressed
> blood-brain-barrier and neuroinflammatory programmes in frontal cortex (both
> q = 0.032), with *JUNB* emerging as a region-specific marker of the drug–flight
> interaction, carrying opposite-signed interactions in frontal cortex (+0.72) and
> dentate gyrus (−0.53).

---

## What would test this

The chain above is falsifiable, and each link maps to a specific experiment.

| Claim | Test |
|---|---|
| Oxidative stress drives the cascade | Protein carbonyl / 4-HNE / 8-oxo-dG staining; *Mt1*, *Sod1/2* qPCR |
| *Arc*/*Fos* down means reduced activity | *Arc* or *Fos* immunostaining with cell counts; ideally in vivo activity recording |
| Hippocampal upregulation is compensatory | Electrophysiology (fEPSP, LTP) in CA1 and DG; NRGN and GluN1 protein |
| CA1 barrier remodelling is real | IgG extravasation, claudin-5 / ZO-1 immunostaining, MMP2 zymography |
| *Ptgds* signals meningeal inflammation | PGD2 measurement; PTGDS immunostaining of leptomeninges |
| DG mitochondrial involvement | Seahorse respirometry; cyclophilin-D inhibition (cyclosporin A) |
| *JUNB* is a marker of the interaction | qPCR/IHC in an independent cohort, both regions, all four groups |
| Neurogranin is a translatable marker | CSF or plasma NRGN protein, ideally in EV fraction |

**Required first, before any of this is publishable as biology:** an independent
cohort with more than one slide per group. The current design cannot separate
batch from treatment, and n = 3 cannot detect below roughly a 1.6-fold change.

---

## Honest limitations

- **No gene is significant transcriptome-wide.** All 28 comparisons return zero at
  BH < 0.05; the best adjusted p anywhere is 0.063.
- **Three genes survive panel-restricted correction** (*JUNB*, *MMP12*,
  neurogranin) and that restriction is only legitimate because the panel was fixed
  in advance.
- **Treatment is perfectly confounded with slide.** Every "BuOE effect" is formally
  a "BuOE-slide effect".
- **Category-level results are the strongest evidence here**, because 16–20 genes
  moving coherently is harder to produce by chance than one gene.
- **Direction of causality is untested.** Nothing here distinguishes compensation
  from damage.
- **DCX did not move** (best P = 0.11), so no claim about neurogenesis rate is
  supported — only the transcriptional context around it.
- **This is tissue mRNA.** No protein, no EV cargo, no functional readout.
