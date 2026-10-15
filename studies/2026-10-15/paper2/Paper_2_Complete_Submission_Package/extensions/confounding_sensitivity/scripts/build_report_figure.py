from pathlib import Path
import json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.colors import ListedColormap
B=Path(__file__).resolve().parents[1];R=B/'results';F=B/'figures'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':1,'lines.linewidth':1.2,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
s=json.loads((R/'summary.json').read_text());t=pd.read_csv(R/'bias_tipping_points.csv');q=pd.read_csv(R/'QC_adjustment_sensitivity.csv');geo=pd.read_csv(R/'comparator_geometry.csv');cur=pd.read_csv(R/'bias_bound_sensitivity_curves.csv');colors=['#0072B2','#D55E00','#009E73'];names=['FLT vs GC','FLT vs VIV','GC vs VIV'];cs=['FLT_vs_GC','FLT_vs_VIV','GC_vs_VIV']
fig,ax=plt.subplots(2,2,figsize=(183/25.4,169/25.4));fig.subplots_adjust(left=.115,right=.98,bottom=.11,top=.92,wspace=.42,hspace=.69)
for i,a in enumerate(ax.flat):a.text(-.22,1.12,'abcd'[i],transform=a.transAxes,fontweight='bold',fontsize=11);a.spines[['right','top']].set_visible(False)
a=ax[0,0];a.set_title('Flight and collection date are aliased',loc='left',pad=12);a.imshow([[1,1,0,1],[1,0,1,0],[1,0,0,0]],cmap=ListedColormap(['#f3f4f4','#0072B2']),vmin=0,vmax=1,aspect='auto')
a.set_xticks(range(4),['Intercept','Flight','Ground','14 Jan']);a.xaxis.tick_top();a.tick_params(length=0);a.set_yticks(range(3),['FLT (9)','GC (9)','VIV (8)']);a.set_ylim(2.5,-.5)
for i,row in enumerate([[1,1,0,1],[1,0,1,0],[1,0,0,0]]):
    for j,v in enumerate(row):a.text(j,i,str(v),ha='center',va='center',color='white' if v else '#333333')
a.text(.5,-.16,'Rank 3 / 4 columns; null vector (0, 1, 0, −1)',ha='center',transform=a.transAxes,fontsize=8);a.text(.5,-.29,'GC and VIV: 17 Jan 2021',ha='center',transform=a.transAxes,fontsize=8)
a=ax[0,1];a.set_title('Unmeasured-bias sensitivity',loc='left',pad=12)
for j,c in enumerate(cs[:2]):
    v=cur[(cur.frame=='multiway_Hellinger')&(cur.contrast==c)];tr=t[(t.frame=='multiway_Hellinger')&(t.contrast==c)].iloc[0]
    a.plot(v.lambda_benchmark,v.plugin_norm_lower/tr.observed_norm,color=colors[j],label=names[j]);a.plot(v.lambda_benchmark,v.approx_joint_sampling_norm_lower/tr.observed_norm,color=colors[j],ls=':')
    a.plot(tr.plugin_tipping_lambda,0,'o',color=colors[j],clip_on=False);a.text(tr.plugin_tipping_lambda,.06,f'{tr.plugin_tipping_lambda:.2f}',rotation=90,ha='center',va='bottom',color=colors[j])
a.set_xlim(0,2);a.set_ylim(-.025,1.05);a.set_xlabel('Assumed bias / GC–VIV norm (λ)');a.set_ylabel('Residual lower bound / observed norm');a.legend(loc='upper right',frameon=False,handlelength=1.5)
a.text(.48,.48,'Solid: plug-in\nDotted: sampling\nallowance',transform=a.transAxes,fontsize=8)
a=ax[1,0];a.set_title('RIN and depth adjustment is small',loc='left',pad=12)
for j,c in enumerate(cs):
    for fi,frame in enumerate(['multiway_Hellinger','binary_RMS_pp']):
        v=q[(q.frame==frame)&(q.model=='group_RIN_log2_depth')&(q.contrast==c)].iloc[0];yy=2-j+(.12 if fi==0 else -.12)
        a.errorbar(v.norm_ratio,yy,xerr=[[v.norm_ratio-v.ratio_bootstrap_low],[v.ratio_bootstrap_high-v.norm_ratio]],fmt='o' if fi==0 else 's',color=colors[j],mfc=colors[j] if fi==0 else 'white',markersize=4,capsize=2)
a.axvline(1,color='#666666',lw=1,ls=':');a.set_yticks(range(3),names[::-1]);a.set_xlabel('Adjusted / unadjusted contrast norm');a.set_xlim(.91,1.22);a.set_ylim(-.5,2.5);a.plot([],[],'o',color='#555555',label='Multiway',ms=4);a.plot([],[],'s',mfc='white',mec='#555555',label='Binary',ms=4);a.legend(loc='lower left',ncol=2,frameon=False,handletextpad=.3,columnspacing=.6,borderpad=0)
a=ax[1,1];a.set_title('Sampling noise affects centroid moments',loc='left',pad=12)
metrics=['FG_squared','FV_squared','GV_squared','FV_dot_GV'];labels=['FLT–GC squared','FLT–VIV squared','GC–VIV squared','Shared-control\ninner product']
for i,m in enumerate(metrics):
    g=geo[(geo.frame=='multiway_Hellinger')&(geo.metric==m)];plug=g[g.estimator=='plugin'].iloc[0];v=g[g.estimator=='Ustat_noise_corrected'].iloc[0]
    a.plot(plug.estimate*1000,3-i+.10,'o',mfc='white',mec='#555555',ms=4)
    a.errorbar(v.estimate*1000,3-i-.10,xerr=[[(v.estimate-v.interval_low)*1000],[(v.interval_high-v.estimate)*1000]],fmt='o',color='#0072B2',ms=4,capsize=2)
a.axvline(0,color='#777777',ls=':',lw=1);a.set_yticks(range(4),labels[::-1]);a.set_ylim(-.65,3.5);a.set_xlabel('Moment × 10³ (Hellinger units squared)');a.set_xlim(-3.9,7.5);a.text(.0,-.28,'○ Plug-in    ● Distinct-animal estimate',transform=a.transAxes,fontsize=8)
# Verify intended text size before saving.
from matplotlib.text import Text
texts=[x for x in fig.findobj(Text) if x.get_text()];assert min(x.get_fontsize() for x in texts)>=8
for ext in ['pdf','svg','png']:fig.savefig(F/f'Figure_Confounding_Sensitivity.{ext}',dpi=600,facecolor='white')
plt.close(fig)
png=F/'Figure_Confounding_Sensitivity.png';im=Image.open(png).convert('RGB');im.save(png,dpi=(600,600))
caption='Figure. Identification limits and measured-QC sensitivity in RR23 cerebellar junction profiles. a, The 26-animal design has identical flight and early-collection indicators and rank three for four columns; no regression can separate their coefficients. b, Primary 923-cluster Hellinger-space sensitivity: the lower bound on a remaining contrast under an assumed additive bias norm no greater than λ times the GC–VIV contrast. Solid curves use observed centroids; dotted curves allow a joint 95% bootstrap error radius across the three contrasts. The same-date GC–VIV magnitude is a hypothetical scale benchmark, not a date-effect estimate. Tipping values label the plug-in bound. c, Contrast norm ratios after simultaneous adjustment for RNA integrity number (RIN) and log2 summed library read depth, relative to the group-only model, for the 923-cluster multiway frame and 5,249 binary choices. Bars are approximate 95% percentile intervals from 10,000 joint within-arm animal bootstrap refits; they quantify measured-QC sensitivity, not causal identification. d, Multiway squared centroid distances and the inner product of FLT–VIV with GC–VIV before and after distinct-animal U-statistic correction for sampling noise. Corrected moment bars are approximate basic bootstrap intervals using analytically centered U-statistic errors; negative lower limits are retained. No inference is made from corrected cosine ratios, because most corrected bootstrap contrast Gram matrices are inadmissible. Animals, rather than reads or overlapping clusters, are the resampling unit. All analyses are retrospective and conditional on the frozen measurement frame.'
(F/'Figure_Confounding_Sensitivity_caption.txt').write_text(caption+'\n');(F/'figure_metadata.json').write_text(json.dumps(dict(width_mm=183,height_mm=169,min_font_pt=min(x.get_fontsize() for x in texts),png_mode='RGB',png_size=Image.open(png).size,dpi=600,caption_words=len(caption.split())),indent=2))
def table(df,cols):
    out=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
    for _,rr in df.iterrows():out.append('| '+' | '.join(f'{rr[c]:.6g}' if isinstance(rr[c],(float,np.floating)) else str(rr[c]) for c in cols)+' |')
    return '\n'.join(out)
report='''# Identification and comparator sensitivity

Author: Isarar Siddique

Affiliation:

Email:

15 October 2026. Retrospective analysis lock plus a transparently logged U-statistic uncertainty amendment. This extension quantifies existing design limitations; it does not claim a causal microgravity effect or biological replication.

## Main findings

The flight indicator and early collection-date indicator are identical across all 26 animals. The design [intercept, flight, ground hardware, early date] has rank 3 for 4 columns, smallest singular value 4.32×10⁻¹⁶ and exact null vector (0,1,0,−1). Adding any constant to the flight coefficient while subtracting it from the date coefficient leaves fitted values exactly unchanged. No increase in algorithmic complexity or precision resolves this lack of identification.

Nevertheless, the observed aggregate profile is insensitive to the two measured quality covariates. Simultaneously adjusting for RIN and log2 summed read depth changes the multiway FLT–GC norm from 0.08313 to 0.08337 (+0.291%) and FLT–VIV from 0.06615 to 0.06662 (+0.707%). The corresponding adjusted/unadjusted vector cosines are 0.99859 and 0.99130. Binary FLT–GC changes by +0.009% and FLT–VIV by −0.053%. Thus measured RIN and library depth do not materially explain the observed centroid pattern under this linear sensitivity model. Adjustment can control measured covariates but cannot separate flight from date; RIN and read depth may also be downstream of treatment or technical processing.

An additive unmeasured-bias vector with norm equal to the observed FLT–GC or FLT–VIV difference can completely erase that difference. Relative to the same-date GC–VIV norm, the plug-in tipping multipliers are 1.269 and 1.010 in the multiway frame; binary values are 1.361 and 0.936. These are quantitative assumption thresholds, not estimates of actual date effects. Under the conservative joint bootstrap sampling allowance, a positive lower norm bound persists only below 0.094× the GC–VIV benchmark for multiway FLT–GC, and it is already zero at λ=0 for FLT–VIV. The latter is uncertainty, not a contradiction of the exact global permutation test, which targets a different null under exchangeability.

Naive shared-control geometry is materially affected by finite-sample noise. In multiway Hellinger coordinates, the FLT–VIV / GC–VIV inner product is 0.000879 before and −0.000894 after distinct-animal moment correction (approximate centered basic interval −0.003323 to 0.000545). The corrected binary inner product is −5.545 percentage-points squared (approximate interval −13.817 to −0.467). These exploratory moment estimates do not establish an antiparallel biological axis: only 45.47% of multiway and 37.39% of binary centered bootstrap Gram matrices admit the required nonnegative two-dimensional geometry. Corrected cosines are ratios of unbiased moment estimates, not unbiased or stable estimators, and are not used as findings.

## Mathematics and estimation

For cluster e and animal i, let p_ie be the frozen vector of junction proportions. Set X_i=concat_e sqrt(p_ie)/sqrt(2E), E=923. Then ||mean_A X−mean_B X|| is the root mean squared Hellinger centroid distance across clusters. For the binary sensitivity frame set X_i=100 f_i/sqrt(5249), giving the RMS difference in percentage points. Samples are 9 FLT, 9 GC and 8 VIV animals; technical replicates are collapsed first. The canonical raw junction matrix was checked independently against every selected multiway junction count.

For group g, the plug-in centroid Gram diagonal is G_gg=||mean_g X||². The distinct-animal estimator replaces this with G_U,gg=[n_g(n_g−1)]⁻¹ sum_(i≠j in g) X_iᵀX_j = G_gg−tr(S_g)/n_g. Cross-group entries remain G_gh=mean_g Xᵀmean_h X. Contrast squared norms and inner products are linear combinations of G_U. These linear moments are unbiased for population-centroid moments under independent animals and a fixed feature frame; the matrix need not be positive semidefinite. The feature frame is conditional on the deposited cohort and its pre-existing measurement gates.

The additive sensitivity model is delta=theta+b, where theta is the residual contrast after an assumed bias vector, not an identified microgravity effect. Under ||b||≤B, the ambient Euclidean bounds are max(0,||delta||−B)≤||theta||≤||delta||+B. They are attained by choosing b parallel or antiparallel to delta. Set B=λ||c|| with c=mean_GC X−mean_VIV X solely as a transparent scenario benchmark. GC and VIV were collected on the same date; their contrast provides no estimate of the three-day flight-versus-ground batch/date bias. Assuming B is small does not remove radiation, habitat, re-entry or other co-exposures from an interpreted flight effect.

Joint within-group resampling uses 10,000 PCG64 count-weight bootstrap draws, seed 2026101531, shared across frames and contrasts. Let r be the 95th percentile of max_k ||delta*_k−delta_k|| across all three contrasts. An approximate conservative lower bound with sampling allowance is max(0,||delta||−r−λ(||c||+r)). It is conditional on the chosen additive bound and the animal-resampling assumptions and is not a causal confidence interval. Existing global permutation results are retained unchanged.

For U-statistic moment uncertainty, ordinary empirical bootstrap percentiles require care: E*(G_U*)=G, the plug-in Gram, exactly. We use centered errors G_U*−G and basic intervals [G_U,contrast−q97.5(error), G_U,contrast−q2.5(error)] for the four linear moments. This was documented in the amendment after the initial diagnostic run. No intervals for corrected nonlinear cosine, position or orthogonal-fraction ratios are interpreted. Raw plug-in bootstrap percentiles are saved only as descriptive resampling ranges; upward norm bias makes them unsuitable for claiming nonzero signal. The correction addresses animal sampling noise and cannot remove batch confounding.

Measured-quality sensitivity fits X_i=alpha+beta_F I_FLT+beta_G I_GC+gamma_R RIN_z+gamma_D log2(depth)_z+epsilon_i by multivariate OLS, alongside group-only and each single-covariate model. RIN is the mean of library RIN values after technical-replicate collapse; depth is their sum. Covariates are standardized in the full sample. The full design has rank 5, 21 residual degrees of freedom and maximum leverage 0.419. All 10,000 within-arm bootstrap refits retain full rank. Output tables include all three contrasts under all four models; no covariate-adjusted significance claim is introduced.

## Complete summary tables

### Bias tipping values

'''
report+=table(t,['frame','contrast','observed_norm','control_benchmark_norm','plugin_tipping_lambda','joint_error_radius95','conservative_tipping_lambda'])
report+='\n\n### All measured-QC sensitivity models\n\n'+table(q,['frame','model','contrast','norm_ratio','cosine_with_unadjusted','ratio_bootstrap_low','ratio_bootstrap_high'])
report+='\n\n### Noise-corrected linear moments\n\n'+table(geo[(geo.estimator=='Ustat_noise_corrected')&geo.metric.isin(metrics)],['frame','metric','estimate','interval_low','interval_high'])
report+='''

## Manuscript-ready compact result

The deposited design was rank deficient when flight and collection date were modeled separately (rank 3 for 4 columns), confirming that date adjustment cannot identify an isolated flight effect. The aggregate junction profile was nevertheless insensitive to measured RIN and sequencing depth: simultaneous adjustment changed the multiway FLT–GC and FLT–VIV centroid norms by +0.29% and +0.71%, with adjusted/unadjusted vector cosines of 0.999 and 0.991. Under an explicit additive bias model, complete erasure of the observed differences required a bias norm equal to 1.269 and 1.010 times the same-date GC–VIV contrast, respectively. GC–VIV was used solely as a hypothetical scale benchmark, not to estimate collection-date effects. Joint animal-bootstrap uncertainty substantially weakened these sensitivity bounds. Distinct-animal centroid moments also showed that shared-control sampling noise affected apparent profile geometry; corrected nonlinear directions were unstable and were not interpreted biologically. These results support technical-covariate robustness of the observed association while quantifying the assumptions and sampling uncertainty that prevent causal attribution.

## Limitations and appropriate use

This is an additional retrospective, cohort-conditional sensitivity family. There is no prospective error guarantee across all analyses, no unmeasured-confounding estimate, no causal identification, no independently replicated junction mechanism and no claim that the bias bound is biologically plausible. Animals may share unrecorded cage-level dependencies; no usable cage randomization factor was available, so animal-bootstrap intervals assume independence at the deposited animal unit. The rank deficiency is exact and assumption-free; bootstrap intervals and noise corrections are conditional on sampling assumptions. These distinctions should be retained when integrating the figure or text into the manuscript.
'''
(B/'SCIENTIFIC_REPORT.md').write_text(report)
(B/'README.md').write_text('''# Confounding/comparator sensitivity extension

Run from the package root using its locked Python environment:

```
python extensions/confounding_sensitivity/scripts/analyse.py
python extensions/confounding_sensitivity/scripts/verify_independently.py
python extensions/confounding_sensitivity/scripts/build_report_figure.py
```

No parent inputs or existing analysis outputs are modified. The dated analysis lock and the separately dated methodological amendment must remain together. Read SCIENTIFIC_REPORT.md for interpretation; results/ contains complete bootstrap, covariate, moment, design and sensitivity tables. figures/ contains the native 183-mm vector PDF/SVG and RGB 600-dpi PNG plus caption and font metadata. Bias benchmarks are scenarios, not date-effect estimates; no causal rescue is claimed.
''')
print('Wrote report and four-panel figure')
