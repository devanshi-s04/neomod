# SUPPORT POLICY — DESIGN NOTE (pre-implementation, 2026-08-03)

**Status: design only. No code changed. E0 threshold-10 decision frozen. TEST sealed.**

## §0. What the diagnostics established

| finding | evidence |
|---|---|
| `support_count >= 1` is a **resolution-dependent per-pixel occupancy mask** | halving the pixel creates MBA support holes |
| it causes the catastrophic 0.01-vs-0.005 disagreement | 115 flips, ΔpAUC −0.168, ΔF1 −0.183 |
| disabling it removes the disagreement | **113 of 115 flips gone**; pAUC disagreement −1.68e-01 → −1.20e-03 (**~140×**) |
| the mask also costs accuracy at the production resolution | at 0.01: pAUC 0.8867 → 0.9535, F1 0.8402 → 0.8825 with the mask off |

**Amended conclusion.** The occupancy mask is the **dominant** problem. It is **not** the only one: a
**small residual resolution effect remains** with the mask disabled — **2 flips**, pAUC disagreement
**−1.20e-03**, max |ΔP| **0.177**, p99 3.1e-03. Resolution therefore has a real but second-order
effect, and this note does not claim otherwise. Any final resolution choice must still be justified
separately.

**`support_mask_min=0` is NOT the production answer.** It leaves unmasked kNN density extrapolating
arbitrarily far from any sample, so an out-of-domain observation still receives a confident
probability. The mask was a bad answer to a real question; removing it without replacement
reintroduces the original failure.

---

## §1–§2. Posterior construction

Remove **all** population-specific occupancy masking from posterior construction. For every
population *c*:

```
P(c | X) = rho_c(X) / SUM_c' rho_c'(X)
```

computed from **unmasked** densities with **identical semantics for every population** — same
estimator, same k, same normalisation, same interpolation. No population is exempt, and no
population is independently zeroed. This directly retires `_support_mask_skip`, whose asymmetry
(NEO exempt, everything else maskable) was the visible symptom of the deeper defect.

## §3–§5. One global support statistic, one abstention

Support is a property of **the observation X against the model as a whole**, not of an individual
population. There is exactly one decision:

```
if  S(X) < s*   ->  ABSTAIN  (return NaN for every population)
else            ->  return the full posterior over populations
```

An abstention is **NaN**, never 0 and never 1. Zeroing one population's term while keeping another's
is precisely the bug being retired.

## §6. Which statistic

Three candidates were considered.

**(a) Total mixture/evidence density `SUM_c rho_c(X)`.** Rejected as the primary gate. It conflates
two very different situations: *"the model has no samples here"* (extrapolation — abstain) and
*"the model has samples here and says objects are genuinely rare"* (a real, informative prediction —
do **not** abstain). It also spans many orders of magnitude between the MBA core and the fast-NEO
tail, so any single global threshold would abstain on essentially all |v| > 2 objects — the exact
population the ±5 grid was built to recover.

**(b) kth-neighbour distance to the UNION of model support. ← RECOMMENDED.**
`d_k(X)` = distance in (vlam, vbeta) from X to the k-th nearest sample **pooled across all
populations**. Rationale:

- **grid-independent by construction** — a property of the sample set, not of any raster, which is
  exactly the failure mode being removed;
- **symmetric across populations by construction** — one union, no per-population decisions;
- **physically interpretable** — deg/day, comparable to the velocity scales already in use;
- it is the quantity the kNN estimator is *built on*, so it measures directly where the estimator
  has data rather than inferring it;
- it separates the two cases (a) conflates: sparse-but-sampled regions keep a small `d_k` and remain
  scored, while true extrapolation grows `d_k` without bound.

**(c) Standardised combination of both.** Deferred. Adds a second free parameter and a
standardisation choice for no demonstrated benefit. Revisit only if (b) shows a failure mode.

**Reported alongside, but not gating:** `SUM_c rho_c(X)` and per-population `d_k`, so a future
review can see whether (c) would have helped.

## §7. Choosing `s*` on CAL without laundering metrics through abstention

Selective abstention can improve any metric by discarding hard rows. Controls, **preregistered
before any threshold is fitted**:

1. **Fix coverage first, then read the metric.** Preregister coverage targets
   **{99.9%, 99.5%, 99.0%, 98.0%}**; for each, set `s*` to the CAL quantile achieving it; report
   metrics **at** those coverages. Never search `s*` to maximise a metric.
2. **Report the whole metric-vs-coverage curve**, not a single point. A policy that only looks good
   at one coverage is not a policy.
3. **Random-abstention control.** At each coverage, also abstain on a *random* subset of the same
   size and report the same metrics. **If the support policy does not beat random abstention at
   matched coverage, it is adding nothing** and must not be adopted.
4. **Abstention composition** reported at every coverage, by **truth population**, **magnitude
   bin**, **sky field**, and **velocity band**. Abstention concentrated on one truth class — in
   particular on NEOs, or on |v| > 2 — is a failure regardless of the aggregate metric.
5. `s*` is chosen on **CAL only** and applied unchanged to TEST, consistent with the existing seal
   discipline.

## §8. Can this be tested without rebuilding maps?

**Yes — no rebuild is required.** Every archive already stores, per population and magnitude bin:

| array | use |
|---|---|
| `density_raw` / `density_unsmoothed` | unmasked posterior construction (§1–§2) |
| `support_count` | diagnostic only under the new policy |
| `nearest_dist` | 1-NN distance to that population's samples |

- **Union 1-NN distance** = `min` over populations of `nearest_dist` — directly available.
- **kth-neighbour distance** is recoverable analytically from the stored density without a rebuild:
  the 2-D kNN estimator satisfies `rho ~ k / (pi * d_k^2)`, so `d_k ~ sqrt(k / (pi * rho))` with
  `k = DEFAULT_K_MAP = 10`. This must be **validated against a direct cKDTree query** on one
  representative center+bin before being trusted, since the stored density is a Bayesian posterior
  mean rather than the plug-in estimate.
- Both existing map sets (0.01 and 0.005) are available, so the new policy's **resolution
  independence can be verified immediately** — the decisive check the current mask fails.

**Recommended first experiment (no rebuild):** apply policy (b) to the frozen 15,422 CAL rows at
both resolutions and confirm that the 0.01-vs-0.005 disagreement stays at the mask-off level
(≈2 flips) rather than the mask-on level (115), while abstention composition stays balanced across
truth populations.

---

## Appendix — residual mask-off behaviour

Filled from `SUPPORT_MASK_MECHANISM.json` / `SUPPORT_MASK_FLIP_CORNERS.csv`:

- residual flips with the mask disabled: **2** (of 15,422 rows)
- pAUC disagreement, mask off: **−1.20e-03**; ROC **+1.22e-03**; F1 **−2.52e-03**;
  max |ΔP| **0.177**; p99 **3.1e-03**

### Flip anatomy (the 115 mask-on flips)

| quantity | value |
|---|---|
| truth MBA among the flips | **113 / 115** |
| MBA support > 0 at **all** coarse corners | 12 |
| MBA support == 0 at **>= 1** fine corner | **115 / 115** |
| coarse-supported -> fine hole (both) | 12 (10.4%) |

**Every one of the 115 flips has at least one fine-grid corner with zero MBA support**, which is the
mechanism stated directly: halving the pixel opens MBA occupancy holes, the bilinear stencil picks
one up, the MBA denominator term is zeroed there, and a true MBA is promoted to P(NEO)=1.
Only 12 had MBA support at *every* coarse corner, so most were already marginal at 0.01 — the
coarse grid was masking the fragility, not avoiding it.

### The same rows with the mask disabled

| truth | n | median \|dP\| ON | max ON | median \|dP\| OFF | max OFF | still flip OFF |
|---|---:|---:|---:|---:|---:|---:|
| MBA | 113 | 0.3448 | 0.9995 | **0.000257** | 0.0475 | **0** |
| NEO | 2 | 0.0570 | 0.0684 | 0.0385 | 0.0767 | **0** |

**Median |ΔP| falls ~1,300x and none of the 115 still flips.** The two largest residual mask-off
changes are a NEO (0.8106 -> 0.7339, |ΔP| 0.0767) and an MBA (0.1371 -> 0.1846, |ΔP| 0.0475) —
ordinary interpolation-scale differences, not regime changes.

### Ten largest MBA changes, with and without the mask

Representative rows (full table in `SUPPORT_MASK_FLIP_CORNERS.csv`):

| map | magbin | P(0.01) on | P(0.005) on | Δ on | P(0.01) off | P(0.005) off | Δ off | MBA sup min 0.01 | 0.005 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dlon+050_lat-01 | 18_20 | 0.00052 | **1** | 0.9995 | 0.00052 | 0.00054 | **2.2e-05** | 1 | **0** |
| dlon+050_lat+01 | mag20 | 0.00068 | **1** | 0.9993 | 0.00067 | 0.00083 | **1.6e-04** | 2 | **0** |
| dlon+000_lat-50 | mag23 | 0.00916 | **1** | 0.9908 | 0.00916 | 0.00800 | **1.2e-03** | 1 | **0** |

The `MBA sup min` columns are the tell: `1 or 2 -> 0` across the refinement in every case. With the
mask off the same rows move by ~1e-4, i.e. the densities barely changed — **only the masking
decision did**.

These residuals are the genuine, second-order resolution effect that remains once the occupancy
mask is removed. They are small but non-zero, and are the reason this note does **not** conclude
that resolution is irrelevant.
