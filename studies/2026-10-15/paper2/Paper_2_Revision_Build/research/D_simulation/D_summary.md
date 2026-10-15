> **Correction added by the parent agent, 15 Oct 2026.** `S` is the number of *designated* effects. Separations that
> could not be placed inside the composition bounds are dropped (`split_group_means` returns `dropped`; see
> `n_true_effects` in `mc_replicates_main.csv`), so the realized true-effect count per replicate averages 16.2 (range
> 11-20), 32.6 (22-39) and 65.4 (53-76) for S = 20, 40 and 80. Statements below such as "despite 40 real effects"
> should read "despite a mean of 32.6 real effects". No other number changes; the manuscript uses the corrected wording.

# D_simulation — discrete p-value floor and BY discovery counts

Addresses the credibility-assessment point that the methodological claim in
*Aggregate robustness exceeds local discovery stability* needs simulation
support: that in small exact-permutation designs the number of BY discoveries is
governed by a discrete p-value floor, so read-depth allocation can switch local
discovery counts between 0 and roughly 20 while the underlying effects are
unchanged.

All numbers below come from the outputs in this directory. Nothing is quoted
from the cohort data; the real-data figures are cited only where stated as the
calibration source.

## Files

| File | Contents |
| --- | --- |
| `sim_by_floor.py` | Self-contained seeded simulation (analytic grid, Monte Carlo, figure) |
| `analytic_floor_grid.csv` | Exact k\* over m × design size |
| `sim_by_floor_results_main.json` | All summary statistics, design record, environment |
| `mc_replicates_main.csv` | One row per replicate × S (1,500 rows) |
| `cosines_main.npz` | All 820,000 complementary-half cosine instances |
| `Figure_9_simulation.pdf` / `.png` | Publication figure, 183 mm wide |
| `Figure_9_simulation_caption.txt` | Figure caption |
| `mc_replicates_check{A,B}.csv`, `sim_by_floor_results_check{A,B}.json` | Reproducibility re-runs (20 replicates) |
| `mc_replicates_final25.csv`, `sim_by_floor_results_final25.json` | From-scratch re-run after the final script edits (25 replicates) |

Run: `python3 sim_by_floor.py --reps 500` (243.7 s, single-threaded, peak RSS
208 MB, Python 3.12.13 / numpy 2.5.1, seed 20261015).

## Design

One contrast of 923 three-alternative clusters is simulated for 18 animals split
9 versus 9, and the local family is corrected at the manuscript's m = 2,769
(923 clusters × 3 contrasts). Each replicate is evaluated three times — at full
depth and on the two complementary binomial-thinned halves (A ~ Binom(N, ½),
B = N − A) of the *same* counts, which is how the manuscript's read splits were
built. The generating effects are held identical across the three conditions, so
any difference in discovery count is attributable to read allocation alone.

Calibration sources, stated rather than tuned to an outcome: the real FLT–GC BY
candidates have median total variation 17.4 pp (q75 19.6, q90 22.0), and real
per-animal selected-alternative depths reconstruct from the deposited fraction
denominators to a median near 80–160 reads. The simulation uses 18.0 pp total
variation (13.4 pp RMS fraction difference) and a median depth of 140.

## Part A — analytic, exact arithmetic

For a balanced n versus n comparison the two-sided exact allocation p value
cannot fall below `p_min = 2 / C(2n, n)`. A block of k tests tied at that floor
attains BY adjusted value `q(k) = p_min · m · H_m / k`, which is *decreasing* in
k, so no rejection is possible at level α unless the family contains at least
`k* = ceil(p_min · m · H_m / α)` floor-level tests. k\* depends only on the
design and the family size — never on effect size.

| m | 5v5 | 8v8 | 9v9 | 10v10 |
| --- | --- | --- | --- | --- |
| 500 | 540 † | 11 | 3 | 1 |
| 1,000 | 1,189 † | 24 | 7 | 2 |
| 2,769 | 3,738 † | 74 | **20** | 6 |
| 5,000 | 7,218 † | 142 | 38 | 10 |

† k\* > m: BY can reject **nothing** at 0.05 in a family of this size at that
design, whatever the effects are.

The m = 2,769 / 9v9 cell reproduces the manuscript exactly and is asserted in
the script: H_m = 8.50363774780702, p_min = 2/48,620, k\* = 20,
q(20) = 0.04842980856371377, q(19) = 0.05097874585654081.

## Part B — Monte Carlo, 500 replicates per S (4,500 evaluated families)

**The discovery count is strictly bimodal.** Across all 4,500 families, *zero*
returned between 1 and 19 discoveries. The smallest nonzero count was exactly
20 = k\*; the largest was 49.

| S | condition | mean | median | range | P(0 discoveries) | mean when > 0 | floor ties (mean, range) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | full | 0.00 | 0 | 0–0 | 1.000 | — | 7.2 (2–16) |
| 20 | half A | 0.00 | 0 | 0–0 | 1.000 | — | 4.7 (0–13) |
| 20 | half B | 0.00 | 0 | 0–0 | 1.000 | — | 4.8 (0–10) |
| 40 | full | 0.98 | 0 | 0–23 | 0.952 | 20.5 | 14.0 (6–23) |
| 40 | half A | 0.00 | 0 | 0–0 | 1.000 | — | 9.5 (2–18) |
| 40 | half B | 0.00 | 0 | 0–0 | 1.000 | — | 9.2 (3–18) |
| 80 | full | 30.60 | 29 | 0–49 | 0.018 | 31.2 | 28.3 (17–42) |
| 80 | half A | 9.60 | 0 | 0–29 | 0.572 | 22.4 | 18.8 (6–29) |
| 80 | half B | 9.87 | 0 | 0–28 | 0.556 | 22.2 | 18.8 (7–28) |

**Complementary halves disagree categorically (S = 80).** Of 500 complementary
pairs, 37.2% had exactly one positive half and one with zero; 37.8% had both
zero and 25.0% both positive. The mean within-pair difference in discovery count
was 8.7 and the maximum 28 — between two halves of the same reads from the same
animals carrying the same true effects.

**Aggregate significance is unaffected.** The aggregate permutation test was
Holm-positive within its three-test family in **all 4,500** families, including
every half that returned zero local discoveries. The largest aggregate p value
observed anywhere was 2.057 × 10⁻⁴.

**Effect direction is unaffected.** Across 57,107 true-effect instances the
complementary-half cosine had median 0.9901, 5th percentile 0.9042, and positive
orientation in 99.998% (one negative instance, minimum −0.3804). Null clusters
gave median 0.6403.

**Correctness checks.** Null p values are calibrated: the fraction of null
clusters with p ≤ 0.05 was 0.0499–0.0508 across the nine conditions. False
discoveries were 0.5–0.7% of all calls, so the instability is not a
false-positive artefact. In 151 of 4,500 families the discovery count exceeded
the number of floor ties, because a second discrete level (4/48,620) also
cleared BY once the floor block was large enough — consistent with the
manuscript's statement that higher ranks can also determine a BY rejection.

## Draft Results paragraph

> A simulation confirms that this boundary, rather than effect instability,
> governs local discovery counts at this design size. For a balanced n versus n
> comparison the two-sided exact allocation p value cannot fall below
> 2/C(2n,n), so a block of k tests tied at that floor attains BY adjusted value
> p_min·m·H_m/k, and no rejection is possible at 0.05 unless the family holds at
> least k\* = ceil(p_min·m·H_m/0.05) floor-level tests. For m = 2,769 and nine
> versus nine animals k\* is 20, reproducing the 0.04843 and 0.05098 values
> above; at five versus five animals k\* exceeds m for every family size from 500
> to 5,000, so BY can reject nothing at that design whatever the effects are.
> Dirichlet–multinomial compositional data were then simulated for 18 animals in
> a nine versus nine split, with S = 20, 40 or 80 true cluster effects of 18.0 pp
> total variation (13.4 pp root-mean-square fraction difference) at a median of
> 140 selected-alternative reads per animal–cluster. The generating effects were
> held identical while full depth and both complementary binomial-thinned halves
> were evaluated, enumerating all 48,620 allocations, for 500 replicates per S
> and 4,500 evaluated families. No family returned between one and 19
> discoveries: the smallest nonzero count was exactly 20. With S = 40, full-depth
> families returned zero discoveries in 95.2% of replicates despite 40 real
> effects, and 20.5 calls on average when they rejected at all. With S = 80, full
> depth returned a median of 29 discoveries (mean 30.60, range 0–49) while its own
> complementary halves returned zero in 57.2% and 55.6% of replicates and 22.4
> and 22.2 calls when positive; 37.2% of complementary pairs had exactly one
> positive half, with a mean within-pair difference of 8.7 calls and a maximum of
> 28. Across the same 4,500 families the aggregate Holm-corrected test was
> positive in every one, with largest aggregate p = 2.057 × 10⁻⁴, and the
> complementary-half cosine of the native effect vectors had median 0.9901 and
> 5th percentile 0.9042, retaining positive orientation in 99.998% of 57,107
> true-effect instances. False calls were 0.5–0.7% of all discoveries. Identical
> underlying effects therefore move the local BY count between zero and 20 or
> more through read allocation alone, while aggregate significance and effect
> direction do not change. The simulation uses a parametric count model and a
> single contrast; it supports the correction-boundary mechanism and no
> biological claim.

## Draft Methods paragraph

> **Discrete-floor simulation.** The floor arithmetic was evaluated exactly for
> m = 500, 1,000, 2,769 and 5,000 and for balanced designs of five, eight, nine
> and ten animals per group, using p_min = 2/C(2n,n), H_m = Σ_{i≤m} 1/i and
> k\* = ceil(p_min·m·H_m/0.05). Monte Carlo data comprised 923 three-alternative
> clusters for 18 animals split nine versus nine. Pooled cluster compositions
> were drawn from a Dirichlet(5,3,2) clipped to the selection bounds, and S of
> them received a sum-zero group separation of 18.0 pp total variation; animal
> compositions were Dirichlet with concentration (1−θ)/θ at θ = 0.017 and reads
> were multinomial at depths drawn lognormally about a median of 140 with a floor
> of 20, matching the coverage rule. Complementary halves were A ~ Binomial(N, ½)
> and B = N − A on the same counts, so the two halves share the biological draw
> and differ only in read allocation. The per-cluster statistic was
> ‖Σ_i s_i z_i‖² on within-animal square-root compositions, whose ranking over
> allocations is identical to the squared distance between transformed group
> means; it was evaluated over all 48,620 allocations by mirror-reduced
> enumeration of the 24,310 balanced sign vectors, giving exact two-sided p values
> with a 10⁻¹² tie tolerance, and the aggregate statistic summed it over
> clusters. S + 500 clusters were permuted explicitly; the remaining tests of the
> m = 2,769 family were drawn from their exact discrete null, which is uniform on
> {2/48,620, 4/48,620, …} under label exchangeability. BY was applied to the
> complete 2,769-value family at 0.05, and the aggregate test received Holm
> correction within a three-test family. Null p values were calibrated, with
> 4.99–5.08% below 0.05 across conditions. The script is seeded and reproduces
> byte-identically; 500 replicates per S ran in 244 s single-threaded.

## Limitations

1. **Parametric, not resampled from the cohort.** Counts come from a
   Dirichlet–multinomial with Dirichlet(5,3,2) base compositions and lognormal
   depths. Real junction data have cluster-specific dispersion, correlated
   clusters sharing genomic boundaries, and overlapping clusters within genes;
   none of that dependence is reproduced. The manuscript's own kernel
   construction explicitly preserves such dependence, and the simulation does
   not.
2. **One contrast, three alternatives.** The m = 2,769 family is corrected at
   its real size, but only one contrast of 923 clusters is generated, all true
   effects sit in that contrast, and every cluster has exactly three
   alternatives. The real family spans three contrasts with shared control
   animals and a variable number of alternatives per cluster.
3. **Partly analytic family.** Only S + 500 clusters are permuted from counts;
   the remaining 2,189–2,249 tests are drawn from the exact discrete
   permutation null. This is exact for null tests under exchangeability but
   makes them independent of each other and of depth condition, which real
   tests on shared animals are not. Their contribution is small — the expected
   number of null tests at the floor is 2,769 × 4.11 × 10⁻⁵ ≈ 0.11 — but it is
   an assumption, not a measurement.
4. **The aggregate test is computed on the permuted subfamily** (S true + 500
   null clusters), not on all 923. Diluting it with the remaining null clusters
   would reduce aggregate power, so the reported 100% Holm positivity is an
   upper bound for that statistic under this design.
5. **Effect sizes are homogeneous within a replicate.** All S true effects share
   one total-variation magnitude, whereas real candidate effects are spread
   (real median 17.4 pp, q90 22.0). Heterogeneous effects would broaden the
   floor-tie distribution and therefore the discovery-count distribution.
6. **S is not estimated from the data.** The number of true effects in the real
   cohort is unknown; S = 20, 40 and 80 is a span chosen to bracket k\* = 20, not
   an estimate. Which S the real data resembles is not determined here.
7. **No biological claim.** The simulation speaks only to the behaviour of BY
   under a discrete exact-permutation floor. It does not support any flight,
   isoform or protein conclusion, and does not address the collection-date
   aliasing that limits causal attribution in the cohort.

## Reproducibility

- `--reps 20 --tag checkA` and `--tag checkB` produce **byte-identical** CSVs
  (sha256 `f9853f2251614b2d12eba0c605612f07daab3e48375da81043a915c451e87564`)
  and identical JSON apart from the timestamp and runtime fields.
- Because per-S replicate seeds are spawned by index, the reduced run reproduces
  the 500-replicate run **row for row** on its first 20 replicates (60/60 rows
  matched, 0 differing).
- A from-scratch re-run after the final script edits (`--reps 25 --tag final25`)
  reproduced the 500-replicate run on all 75 shared rows, 0 differing.
- BLAS threading is pinned to one thread inside the script before numpy is
  imported, so results do not depend on the host thread count.
- The PDF measures 183.0 × 126.0 mm and embeds only Arial subsets
  (`ArialMT`, `Arial-BoldMT`, `Arial-ItalicMT`) as Type 42, including the
  mathtext.
