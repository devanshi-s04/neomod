# NEOMOD3 sample count x kNN neighbour count — NEO velocity-density maps

One sky map, one magnitude bin, two NEO source counts, seven k values: fourteen NEO density maps.

**Question.** How do the NEOMOD3 Monte Carlo sample count and the kNN neighbour count change the
smoothness and structure of the NEO velocity-density map?

This is a controlled visual comparison. No preferred k is selected here, and no
calibration, threshold, ROC or F1 is computed. CAL and TEST are not read.

## Fixed

| | |
|---|---|
| sky center | `dlon+000_lat+00` (lon 331.2030, lat 0.0) |
| magnitude bin | `mag24+` (24 <= mag < 25) |
| epoch | 2027-08-25T00:00:00 |
| velocity domain | [-5.0, 5.0] deg/day, step 0.01 |
| patch radius | 30 deg |
| estimator | Bayesian kNN, sealed module `a6de18c2197cbcb9...` |
| Gaussian smoothing | **False** |
| support masking | False |
| quadrature `n_d0_grid` | 8000 (all cases) |

## Source cases

| | rows | n_draws | effective_factor | total physical weight |
|---|---|---|---|---|
| BASE | 27,781 | 100,000,000 | 8.746673 | 3,176.1789 |
| HIGH | 207,059 | 740,000,000 | 64.725384 | 3,199.0386 |

`HIGH = all BASE rows + 179,278 additional independent NEOMOD3 draws`, nested, with
seeds 1,000,000+0..319 (disjoint from the original cache's
42..141). No duplicated,
jittered, resampled or GMM-cloned rows: 0 duplicate orbital records.

**Physical abundance is preserved, not the row count.** Raising the draw budget raises `n_draws`
and therefore `effective_factor = n_draws / total_weight` by the same factor, so the modelled
physical NEO abundance the sample represents is unchanged. The two totals differ by
7.197e-03 (tolerance 2e-02). That residual is
Monte Carlo sampling noise, and it is dominated by BASE: with 27,781 samples BASE carries
~0.60% relative noise on its own, against ~0.22% for HIGH.

## A numerical point that matters for this comparison

The sealed posterior integrand is `p(d0) ~ d0^-k(k+1) exp(-S/d0^2)` with `S = sum_j d_j^2`, so its
peak narrows sharply as k grows. Measured against the exact closed form
`n0 = (k(k+1)/2 - 1/2) / (pi S)`, the production default of 400 quadrature nodes gives max relative
deviation 7.1e-16 at k=10 but **2.6e-04 at k=50 and 1.1e-02 at k=100**. Left at 400, that error
would grow with k and show up as a structural difference between panels. All cases therefore
use 8000 nodes, converged to <= 4e-15 at every k and re-verified at k=150/200/250. Production maps at k=10 are unaffected.

## k-th neighbour distance, central [-1,+1] deg/day window (deg/day; pixel = 0.01)

| source | k | rows | p10 | p50 | p90 | p99 |
|---|---|---|---|---|---|---|
| BASE | 10 | 27,781 | 0.0169 | 0.0310 | 0.0514 | 0.0745 |
| BASE | 25 | 27,781 | 0.0279 | 0.0500 | 0.0796 | 0.1156 |
| BASE | 50 | 27,781 | 0.0399 | 0.0715 | 0.1126 | 0.1609 |
| BASE | 100 | 27,781 | 0.0566 | 0.1011 | 0.1587 | 0.2231 |
| BASE | 150 | 27,781 | 0.0694 | 0.1233 | 0.1933 | 0.2695 |
| BASE | 200 | 27,781 | 0.0799 | 0.1419 | 0.2222 | 0.3068 |
| BASE | 250 | 27,781 | 0.0892 | 0.1584 | 0.2473 | 0.3417 |
| HIGH | 10 | 207,059 | 0.0062 | 0.0113 | 0.0190 | 0.0279 |
| HIGH | 25 | 207,059 | 0.0103 | 0.0183 | 0.0297 | 0.0425 |
| HIGH | 50 | 207,059 | 0.0147 | 0.0260 | 0.0419 | 0.0602 |
| HIGH | 100 | 207,059 | 0.0209 | 0.0369 | 0.0593 | 0.0850 |
| HIGH | 150 | 207,059 | 0.0256 | 0.0452 | 0.0725 | 0.1043 |
| HIGH | 200 | 207,059 | 0.0296 | 0.0522 | 0.0836 | 0.1204 |
| HIGH | 250 | 207,059 | 0.0331 | 0.0584 | 0.0934 | 0.1344 |

## Contents

- `more_neomod_samples_knn_maps.ipynb` — executed notebook (the deliverable to read first)
- `source_base.parquet` — sha256 `03dd6e13a3a4896c1d5d478d3490dc4a7ef515f0a407ccb7b74b46631d04801c`
- `source_high.parquet` — sha256 `40490b3bc4ffaec122919981396168299c1e84a384dd345c46f8a7adb20fc297`
- `maps/density_{BASE,HIGH}_k{010,025,050,100,150,200,250}.npz` — fourteen NEO density maps
- `maps/posterior_{BASE,HIGH}_k{010,025,050,100,150,200,250}.npz` — P(NEO) per case
- `provenance.json`, `acceptance_checks.json`

`P(NEO) = rho_NEO / (rho_NEO + rho_MBA + rho_TNO + rho_Trojan)` uses the fixed production non-NEO
densities, byte-identical across all fourteen cases (hashes in `provenance.json`). Undefined pixels
are NaN, never 0. Max |sum of four class probabilities - 1| over defined pixels:
3.331e-16.

## Acceptance

All checks: **PASS** — see `acceptance_checks.json`.

Reproduce with `neomod/pipeline/knn_sample_experiment.py`
(`export-base` / `draw-high` / `merge-high` / `map` / `posterior` / `finalize` / `bundle`)
plus `neomod/pipeline/slurm/knn_draw_high.sbatch` and `knn_maps.sbatch`.
Code commit `26e88437e7feee3796e04d0d2812bd41bae863e4`.
