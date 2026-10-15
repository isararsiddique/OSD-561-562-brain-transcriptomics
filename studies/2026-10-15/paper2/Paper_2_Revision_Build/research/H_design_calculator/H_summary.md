# H_design_calculator — how many animals per group before a local discovery is possible

A mission-planning contribution for Paper 2. Part 1 is exact arithmetic: the smallest design at which a
*single* per-cluster (local) finding can ever clear the family correction. Part 2 is a seeded Monte Carlo
power calculation at realistic junction-count parameters. Part 3 is Figure 11 and a recommendation table.

Every number below comes from the outputs in this directory. Nothing is quoted from the cohort data; the
real-data figures appear only where stated as the calibration source.

## Files

| File | Contents |
| --- | --- |
| `design_calculator.py` | Self-contained seeded calculator (attainability, power, figure) |
| `attainability_grid.csv` | Exact k\*, p floor, BY and Holm flags over n = 4–15 × m ∈ {500, 923, 2769, 5000} |
| `power_cells_main.csv` | One row per design cell (48 cells) with every power metric |
| `power_replicates_main.csv` | One row per replicate × cell (19,200 rows) |
| `design_recommendation_main.csv` | The recommendation table |
| `design_calculator_results_main.json` | All of the above plus design record, seeds, environment |
| `Figure_11_design.pdf` / `.png` | Publication figure, 183.0 × 132.0 mm |
| `Figure_11_design_caption.txt` | Figure caption |
| `power_*_check.csv`, `design_calculator_results_check.json` | Reproducibility re-run (40 replicates) |
| `run_main.log`, `run_check.log` | Console logs of both runs |

Run: `python3 design_calculator.py --reps 400` — 801.2 s single-threaded, peak RSS 87.3 MB,
Python 3.12.13 / numpy 2.5.1, macOS-26.4-arm64, seed 20261015. BLAS threading is pinned to one thread
inside the script before numpy is imported.

## Part 1 — attainability (exact, no simulation)

A balanced n vs n two-sided exact allocation test cannot return a p value below
p_min = 2/C(2n, n). Equivalently, mirror-pair reduction leaves R = C(2n−1, n−1) = C(2n, n)/2 distinct
allocation statistics and p = count/R, so p_min = 1/R. A block of k tests tied at that floor attains BY
adjusted value q(k) = p_min · m · H_m / k, which **decreases** in k, so no BY rejection is possible at
level α unless the family holds at least

> **k\*(m, n) = ceil( m · H_m / (R · α) )**

tests at the floor. k\* depends only on the design and the family size — never on the effect size, the
sequencing depth or the biology.

Exact k\* at the manuscript's family m = 2,769 (H_m = 8.50363774780702, α = 0.05):

| n per group | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13–15 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C(2n, n) | 70 | 252 | 924 | 3,432 | 12,870 | 48,620 | 184,756 | 705,432 | 2,704,156 | ≥ 10,400,600 |
| k\* | 13,456 † | 3,738 † | 1,020 | 275 | 74 | **20** | 6 | 2 | **1** | 1 |

† k\* > m: BY can reject **nothing** at 0.05 in a family of this size, whatever the effects are. The same
holds at n = 4 and n = 5 for all four family sizes examined (eight cells in `attainability_grid.csv`).

The m = 2,769 / 9v9 cell reproduces Figure 9 and `../D_simulation/` exactly (k\* = 20,
p_min = 2/48,620 = 4.114 × 10⁻⁵), and the script asserts it.

**Two design thresholds, per family size.**

| m | H_m | k\* at 9v9 | smallest n for one isolated BY discovery | smallest n for one isolated Holm discovery | permutation budget R_min |
| --- | --- | --- | --- | --- | --- |
| 500 | 6.7928 | 3 | 10 (20 animals) | 9 (18 animals) | 67,929 |
| 923 | 7.4054 | 6 | 11 (22 animals) | 9 (18 animals) | 136,704 |
| **2,769** | 8.5036 | **20** | **12 (24 animals)** | **10 (20 animals)** | **470,932** |
| 5,000 | 9.0945 | 38 | 12 (24 animals) | 11 (22 animals) | 909,451 |

Two consequences a mission can act on:

1. **At nine versus nine, m = 2,769, a single isolated cluster can never be reported under BY** — twenty
   tests must sit at the floor simultaneously. Twelve per group (24 animals) is the smallest design at
   which one cluster on its own can clear BY; ten per group suffices under Holm, so BY's H_m factor costs
   two extra animals per group.
2. **The same arithmetic constrains the implementation, not only the design.** A random-permutation test
   with R drawn allocations also has floor 1/R, so the permutation budget must satisfy
   R ≥ R_min(m) = ceil(m · H_m / α) = **470,932** allocations at m = 2,769 before a single test can clear
   BY, however many animals were flown. A study that flies 12 per group and then runs 10,000 permutations
   has thrown the extra animals away.

## Part 2 — power (Monte Carlo, 400 replicates per cell)

One contrast of 923 three-alternative clusters, S = 80 of them carrying a designated group separation,
corrected inside the real family m = 2,769. Allocation p values are **exact by complete mirror-reduced
enumeration for n ≤ 10** (462 / 6,435 / 24,310 / 92,378 levels at n = 6 / 8 / 9 / 10) and estimated from
**480,000 random balanced allocations including the observed one** at n = 12 and n = 15, where enumeration
is infeasible. 480,000 exceeds R_min = 470,932, so k\* = 1 at those designs exactly as in the exact design
and the power curve stays design-limited rather than resampling-limited.

**Grid reduction (stated as required).** The full factorial — 6 design sizes × 3 effect sizes × 3
dispersions × 3 depths — is 162 cells. It was reduced to **48 cells**:

- `band` (18 cells): n ∈ {6, 8, 9, 10, 12, 15} × three θ/depth scenarios at TV = 17.4 pp —
  *favourable* (θ 0.010, depth 280), *calibrated* (θ 0.017, depth 140), *adverse* (θ 0.028, depth 70).
  These are the two extreme corners and the centre of the θ × depth grid.
- `tv` (12 cells): n × TV ∈ {10, 25} pp at the calibrated θ and depth.
- `grid` (18 cells): the **full** θ × depth 3 × 3 at n = 9 and n = 10, where enumeration is cheap. This is
  the check that the three scenarios really do bracket all nine combinations.

The check passes. Per-effect power over the full 3 × 3 is monotone increasing in depth and decreasing in θ
at both n, and spans 0.001–0.859 at n = 9 and 0.213–0.955 at n = 10, against band endpoints 0.001–0.867
and 0.207–0.956. Those endpoints are the *same* parameter combinations run under independent seed
streams, so the 0.001–0.008 gaps give an empirical Monte Carlo standard error on per-effect power of about
0.005.

**Calibrated scenario** — TV 17.4 pp (cohort median candidate total variation), θ 0.017 (cohort median
Dirichlet–multinomial overdispersion), 140 selected-alternative reads per animal–cluster, S = 80
designated (mean 65.6–66.0 realized), m = 2,769:

| n per group | k\* | expected BY discoveries | median (max) | P(≥ 1 discovery) | per-effect power | FDR |
| --- | --- | --- | --- | --- | --- | --- |
| 6 | 1,020 | 0.00 | 0 (0) | 0.000 | 0.000 | — |
| 8 | 74 | 0.00 | 0 (0) | 0.000 | 0.000 | — |
| 9 | 20 | 26.69 (sd 8.81) | 26 (47) | 0.943 | 0.402 | 0.0050 |
| 10 | 6 | 45.05 (sd 5.69) | 46 (60) | 1.000 | 0.683 | 0.0060 |
| 12 | 1 | 58.06 (sd 4.41) | 58 (71) | 1.000 | 0.875 | 0.0054 |
| 15 | 1 | 64.00 (sd 3.87) | 64 (75) | 1.000 | 0.968 | 0.0057 |

**Sensitivity.** Per-effect power at the three θ/depth scenarios and the two off-centre effect sizes:

| n | favourable | calibrated | adverse | TV 10 pp | TV 25 pp |
| --- | --- | --- | --- | --- | --- |
| 9 | 0.867 | 0.402 | 0.001 | 0.000 | 0.888 |
| 10 | 0.956 | 0.683 | 0.207 | 0.012 | 0.972 |
| 12 | 0.993 | 0.875 | 0.504 | 0.169 | 0.996 |
| 15 | 1.000 | 0.968 | 0.747 | 0.382 | 1.000 |

Dispersion and depth move the answer by more than one design step: at the cohort's own 9v9 the same
17.4 pp effects give per-effect power 0.867 at θ 0.010 / depth 280 and 0.001 at θ 0.028 / depth 70.

**False discovery rate.** Across the 31 cells in which a discovery occurred in at least half the
replicates, the mean false discovery proportion was 0.0040–0.0060 and false calls were 0.45–0.64% of all
calls — an order of magnitude below the nominal 0.05, which is the expected conservatism of BY applied to
independent tests. The six marginal cells (0 < P(≥ 1 discovery) < 0.5) show lower values still, because
they rarely reject at all.

## Part 3 — design recommendation table

Family m = 2,769 (923 clusters × 3 contrasts), α = 0.05, S = 80 designated effects. Grid values only:
every interpolated crossing in the n = 9–12 range straddles a change in k\* (20 → 6 → 2 → 1), which is a
genuine discontinuity rather than a slope, so interpolation is not reportable here and the script flags it
(`interp_*_crosses_kstar_step`).

| Requirement | favourable θ 0.010 / depth 280 | **calibrated θ 0.017 / depth 140** | adverse θ 0.028 / depth 70 | TV 10 pp | TV 25 pp |
| --- | --- | --- | --- | --- | --- |
| One isolated cluster attainable under BY (effect-free) | **12** | **12** | **12** | **12** | **12** |
| One isolated cluster attainable under Holm (effect-free) | **10** | **10** | **10** | **10** | **10** |
| ≥ 95% chance of any BY discovery | 9 | **10** | 12 | 12 | 9 |
| ≥ 50% per-effect power | 9 | **10** | 12 | not reached at n ≤ 15 (0.382) | 9 |
| ≥ 80% per-effect power | 9 | **12** | not reached at n ≤ 15 (0.747) | not reached at n ≤ 15 | 9 |
| Permutation allocations needed | 470,932 | 470,932 | 470,932 | 470,932 | 470,932 |

Read off the calibrated column: **ten animals per group to be reasonably sure of returning something, twelve
to have 80% power on an individual cluster, and twelve in any case before a lone cluster can be reported at
all.** The cohort's nine versus nine sits below all three.

## Draft Results paragraph

> The same correction boundary can be inverted into a design calculator, which is what a future mission
> needs. Because k\* = ceil(m·H_m/(R·α)) depends only on the design and the family size, the smallest
> design at which a single isolated cluster can ever clear BY is fixed before any data are collected: for
> m = 2,769 it is twelve animals per group (24 animals), falling to ten per group under Holm, while at nine
> versus nine twenty tests must reach the floor simultaneously and at four or five per group BY can reject
> nothing at 0.05 whatever the effects are. The same arithmetic bounds the implementation as well as the
> design, since a random-permutation test with R allocations also has floor 1/R: at m = 2,769 a single test
> cannot clear BY unless the permutation budget reaches R ≥ 470,932 allocations, so extra animals are
> wasted if the resampling budget does not grow with them. Power was then simulated at the cohort's own
> calibration — Dirichlet–multinomial counts with overdispersion θ = 0.017, a median of 140
> selected-alternative reads per animal–cluster, and 80 designated cluster effects of 17.4 pp total
> variation among 923 clusters (mean 65.9 realized) — with exact enumeration of all allocations for
> n ≤ 10 and 480,000 random allocations at n = 12 and 15, 400 replicates in each of 48 design cells.
> At nine versus nine the expected number of BY discoveries was 26.69 (median 26, maximum 47), the chance
> of at least one discovery 0.943, and per-effect power 0.402; ten per group gave 45.05 expected
> discoveries and 0.683 power, twelve gave 58.06 and 0.875, and fifteen gave 64.00 and 0.968, with the
> mean false discovery proportion between 0.0040 and 0.0060 throughout. Dispersion and depth move the
> answer by more than one design step: the same 17.4 pp effects at nine versus nine give per-effect power
> 0.867 at θ = 0.010 with 280 reads and 0.001 at θ = 0.028 with 70 reads, and a 10 pp effect does not
> reach 50% power anywhere up to fifteen per group. Ten animals per group are therefore needed before a
> local discovery is reasonably certain to appear at all and twelve before an individual cluster carries
> 80% power at the calibrated effect size, against the nine flown here. The calculation uses a parametric
> count model with independent clusters and a single contrast; the attainability thresholds are exact and
> assumption-free, but the power figures are not, and neither supports any biological claim.

## Draft Methods paragraph

> **Design calculator.** The attainability arithmetic was evaluated exactly for m = 500, 923, 2,769 and
> 5,000 and for every balanced design from four to fifteen animals per group, using p_min = 2/C(2n,n),
> H_m = Σ_{i≤m} 1/i and k\* = ceil(p_min·m·H_m/α) at α = 0.05, together with the Holm condition
> p_min·m ≤ α and the permutation-budget condition R ≥ ceil(m·H_m/α). Power was simulated for one
> contrast of 923 three-alternative clusters, of which S = 80 received a sum-zero group separation of
> prescribed total variation; separations that could not be placed inside the composition selection bounds
> were dropped, so the realized effect count is reported alongside every cell. Pooled cluster compositions
> were drawn from a Dirichlet(5,3,2) clipped to the selection bounds, animal compositions were Dirichlet
> with concentration (1−θ)/θ, and reads were multinomial at per-animal–cluster depths drawn lognormally
> (σ_log = 0.70) about the stated median with a floor of 20 reads, matching the coverage rule. The
> per-cluster statistic was ‖Σ_i s_i z_i‖² on within-animal square-root compositions, whose ranking over
> allocations is identical to the squared distance between transformed group means and which is invariant
> to a global sign flip. It was evaluated over all C(2n−1,n−1) mirror-reduced balanced sign vectors where
> C(2n,n) ≤ 200,000 — 462, 6,435, 24,310 and 92,378 levels at six, eight, nine and ten animals per group —
> giving exact two-sided p values with a 10⁻¹² tie tolerance, and over 480,000 uniformly drawn balanced
> allocations including the observed labelling at twelve and fifteen per group, a budget chosen to exceed
> R_min = 470,932 so that the resampling floor is not the binding constraint. The remaining
> m − 80 = 2,689 members of the family were drawn from their exact discrete permutation null, which is
> uniform on {1/R, …, 1} under label exchangeability. Benjamini–Yekutieli was then applied to the complete
> 2,769-value family at 0.05, and expected discoveries, the probability of at least one discovery,
> per-effect power (true rejections divided by realized effects) and the mean false discovery proportion
> were accumulated over 400 replicates in each of 48 cells (19,200 evaluated families). The grid was
> reduced from the full 162-cell factorial to 48 cells: all six design sizes at three θ/depth scenarios
> spanning the corners and centre of the θ × depth grid, all six at effect sizes of 10 and 25 pp, and the
> complete θ × depth 3 × 3 at nine and ten per group, which confirmed that per-effect power is monotone in
> both parameters and that the three scenarios bracket all nine combinations. The script is seeded and
> reproduces exactly; 400 replicates per cell ran in 801 s single-threaded at a peak resident size of
> 87 MB.

## Limitations

1. **Parametric, not resampled from the cohort.** Counts come from a Dirichlet–multinomial with
   Dirichlet(5,3,2) base compositions and lognormal depths. Real junction data have cluster-specific
   dispersion, correlated clusters sharing genomic boundaries, and overlapping clusters within genes. None
   of that is reproduced.
2. **Independence between clusters.** The 80 designated clusters are independent draws and the 2,689
   remaining family members are drawn independently from their exact discrete null. This matters
   asymmetrically: Part 1 is exact arithmetic and holds under arbitrary dependence, because BY itself
   does, but the Part 2 power, discovery-count and false-discovery numbers all assume independence and
   would change under realistic dependence.
3. **One contrast.** The family is corrected at its real size m = 2,769, but only one contrast of 923
   clusters is generated, all designated effects sit in it, and every cluster has exactly three
   alternatives. The real family spans three contrasts with shared control animals and a variable number
   of alternatives per cluster.
4. **Homogeneous effect magnitude within a replicate.** All designated effects in a replicate share one
   total-variation magnitude, whereas real candidate effects are spread (cohort median 17.4 pp). Mixed
   magnitudes would broaden the discovery-count distribution.
5. **The three effect sizes do not hold the realized effect count fixed.** Because a separation that
   cannot be placed inside the composition bounds is dropped, the mean realized count was 74.6 at 10 pp,
   65.9 at 17.4 pp and 55.9 at 25 pp. Per-effect power is normalised by the realized count, but expected
   discovery counts are not comparable across the three effect sizes for that reason.
6. **n = 12 and n = 15 use 480,000 random allocations, not exhaustive enumeration.** That exceeds
   R_min = 470,932, so k\* = 1 exactly as in the exact design, but the p values are estimates with floor
   1/480,000 rather than 2/C(2n,n), making those two columns marginally conservative.
7. **S = 80 is a designated count, not an estimate.** The number of true effects in the real cohort is
   unknown. S = 80 was chosen to match `../D_simulation/`, not inferred from the data, and expected
   discovery counts scale with it. Per-effect power and the attainability thresholds do not.
8. **Grid values only.** Every power threshold crossing in the n = 9–12 range straddles a change in k\*,
   so the interpolated n values in the JSON are flagged and should not be quoted.
9. **The calibration inputs are reconstructions.** The median selected-alternative depth near 140 reads
   was reconstructed from the deposited fraction denominators rather than read from a sequencing report
   (see `../D_simulation/D_summary.md`), and the 17.4 pp median total variation is taken over the real BY
   candidates, which is a selected set.
10. **No biological claim.** The calculator speaks only to the behaviour of a discrete exact-permutation
    test under BY at a given design size. It supports no flight, isoform or protein conclusion, and does
    not address the collection-date aliasing that limits causal attribution in the cohort.

## Reproducibility

- Per-cell replicate seeds are spawned by index from a fixed cell list, so a reduced run reproduces the
  full run row for row. `--reps 40 --tag check` matched the 400-replicate run on **all 1,920 shared rows,
  0 differing**, with identical content hash (`9993ee3a5c1fb80c` over both row sets).
- The script asserts its three load-bearing exact results before any simulation runs: k\*(2769, 9v9) = 20,
  smallest n for a single BY rejection = 12, smallest n for a single Holm rejection = 10, and that the
  Monte Carlo allocation budget exceeds R_min = 470,932.
- The random allocation sets are built from a dedicated seed sequence keyed on (n, R), so they are
  identical across runs and across cells.
- BLAS threading is pinned to one thread before numpy is imported, so results do not depend on the host
  thread count.
- The figure PDF measures 183.0 × 132.0 mm and embeds only Arial subsets (`ArialMT`, `Arial-BoldMT`) as
  Type 42.
