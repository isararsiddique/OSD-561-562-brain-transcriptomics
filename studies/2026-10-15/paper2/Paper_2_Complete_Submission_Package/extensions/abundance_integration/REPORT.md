# Gene abundance and local junction usage

The deposited STAR gene-count matrix matched publisher MD5. Twenty-nine libraries were summed into the same 26 source animals used for junction analysis. Median-ratio normalization used 19,927 genes positive in every animal. Of the fixed 1,983 junction genes, 1,977 had at least 20 raw gene reads in every animal; six exclusions are retained in the universe table. This is a restricted observable universe, not genome-wide gene-expression discovery.

Exact two-sided abundance tests were corrected across all 5,931 tests. Joint BH detects 72 FLT–GC, 8 GC–VIV and 24 FLT–VIV genes, whereas joint BY detects none. Global mean-squared log-abundance contrasts also do not pass Holm: adjusted p=0.133032, 0.423200 and 0.423200. These outcomes neither establish unchanged abundance nor prove a difference between the strengths of two statistical tests.

The descriptive FLT–GC correlation between junction RMS and absolute log-abundance contrast is 0.0572 across 1,977 genes. No gene-independence p value is attached. Forty-six primary junction-BY genes meet the abundance gate; The excluded genes are Mapkapk5, Tomm6, Vps25, Uba52, Gm28557, ENSMUSG00000121470.

Welch TOST sensitivity analyses use operational equivalence intervals ±log2(1.2) and ±log2(1.5), corresponding to geometric-mean ratios 1/1.2–1.2 or 1/1.5–1.5. These margins were declared before this fit but after the junction outcomes were known. None of the junction-association AND abundance-equivalence conjunctions survives joint BY at either margin. Therefore this study does not establish an abundance-independent junction phenotype. Separate abundance equivalence tests do pass for some genes, but selecting their intersection with separately significant junction tests is not substituted for the corrected conjunction family.

The analysis is retrospective and approximate Welch equivalence depends on its distributional assumptions. Exact permutation tests preserve whole-animal vectors and require exchangeability under the null. Flight/date aliasing, normalization assumptions and bulk-tissue composition remain. Relative gene counts and local junction fractions are different measurements from the same libraries and are not independent modalities.

Scripts accept the beta input location via --base. Full measurements, all unsuccessful tests, QC, source provenance and independent arithmetic checks accompany this report. The four-panel figure is descriptive where specified. See RETROSPECTIVE_ANALYSIS_PLAN.md for the declared family and margins.
