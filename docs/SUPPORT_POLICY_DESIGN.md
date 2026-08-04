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

**(b) RAW kth-neighbour distance to the pooled union. REJECTED (amended 2026-08-03).**
Sample spacing differs strongly by population, by map and by magnitude bin — MBA cores are orders
of magnitude denser than the fast-NEO tail. A single distance threshold in deg/day would therefore
reject **sparse but genuinely supported** fast-NEO regions while passing everything near the MBA
core. That is the same population-dependent failure the occupancy mask had, in a different guise.

**(b') STANDARDISED global support score. ← RECOMMENDED.**

```
z_c(X) = d_k,c(X) / q_c(map, mag_bin)
Z(X)   = min_c z_c(X)
```

- `d_k,c(X)` — the **direct** kth-neighbour distance from X to population *c*'s GEN samples.
- `q_c(map, mag_bin)` — a reference spacing for population *c* **in that same map and magnitude
  bin**, derived from that population's own GEN support (leave-one-out kth-neighbour distances).

`z_c` is therefore dimensionless: "how far is X, in units of how far apart this population's own
samples are, here". A sparse fast-NEO region has a large `q_NEO`, so a genuinely supported object
there still gets a small `z_NEO`.

One global decision, no per-population zeroing:

```
Z(X) <= s*  ->  return the FULL unmasked posterior using every population
Z(X) >  s*  ->  ABSTAIN for the observation as a whole (NaN)
```

The `min` over populations means an observation is scored **if at least one modelled population
genuinely supports it** — which is the correct condition, since the posterior is a ratio and needs
only one credible source of density to be meaningful.

**(c) Combining Z with total mixture density.** Deferred. Revisit only if (b') shows a failure mode.

**Reported alongside, but not gating:** `SUM_c rho_c(X)`, per-population `z_c`, and the `argmin_c`
(which population is supplying the support), so a future review can see whether (c) would help.

## §7. Choosing `s*` on CAL without laundering metrics through abstention

Selective abstention can improve any metric by discarding hard rows. Controls, **preregistered
before any threshold is fitted**:

1. **Fix coverage first, then read the metric.** Preregister coverage targets
   **{99.9%, 99.5%, 99.0%, 98.0%}**; for each, set `s*` to the CAL quantile achieving it; report
   metrics **at** those coverages. Never search `s*` to maximise a metric.
2. **Report the whole metric-vs-coverage curve**, not a single point. A policy that only looks good
   at one coverage is not a policy.
3. **Random-abstention control — retained, but NOT the sole adoption criterion** (amended
   2026-08-03). At each coverage, also abstain on a random subset of the same size and report the
   same metrics. Beating random on ROC/F1 is *informative*, not decisive: a correct support policy
   might abstain on genuinely ambiguous rows without improving ranking at all. Adoption is judged
   on the whole of:
   - **calibration / error rate versus Z** — error should rise monotonically with Z;
   - **score stability between 0.01 and 0.005** on retained rows (the check the occupancy mask
     fails: 115 flips vs 2);
   - **abstention composition** by truth population, magnitude bin, sky field, velocity band;
   - **retention of genuine |v| > 2 NEOs** — the population the +-5 grid exists to recover;
   - **whether known extrapolation cases receive large Z**, i.e. the statistic fires where it
     should.
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
- **`d_k` MUST NOT be inferred from the stored density** (amended 2026-08-03). The tempting
  inversion `d_k ~ sqrt(k / (pi * rho))` assumes the simple plug-in kNN estimator, but the stored
  density is a **Bayesian posterior mean** over `d_0` (`estimate_density_full_posterior_2d`), which
  is a different quantity. For the pilot, compute `d_k` **directly from the exact GEN cKDTrees**.
  If the policy is adopted, future map archives must store an explicit **`kth_dist__POP__BIN`**
  array so production scoring needs no tree at query time.
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

---

# §9. PILOT RESULT AND FINAL DECISION (2026-08-03)

## §9.1 Three corrections to the first reading of the pilot

**(1) Error is NOT monotonic in Z.** Across the ten Z deciles the mean absolute error is
0.00246, 0.00282, 0.00300, 0.00230, 0.00230, 0.00199, 0.00246, 0.00386, 0.00635, **0.02740** — it
**fluctuates across the lower deciles** and then **rises strongly in the upper tail**. The correct
statement: *high Z is associated with greater error and sparser support, especially in the final
decile.* The earlier "monotonically increasing" claim is withdrawn.

**(2) High Z is NOT evidence of model extrapolation.** CAL NEOs are **independent NEOMOD3 draws**
and are therefore **in-model by construction**. Z measures **sparse finite-sample representation**,
not out-of-model-ness, and it **cannot distinguish a rare-but-valid observation from a genuinely
out-of-model one**. The earlier framing of Z as detecting "where the model is extrapolating" is
withdrawn.

**(3) The magnitude interpretation was wrong.** At 99% coverage the faintest bin (24–25) has a
**within-bin abstention rate of 1.266%** but supplies **86 of 155 = 55.5% of ALL abstentions**.
Both numbers must be reported: the rate is low, the *share* is a majority. Velocity remains the
strongest concentration (100% of |v|>2 abstained), but **depth may also contribute** and the earlier
dismissal of a depth effect is withdrawn.

## §9.2 Pilot outcome

| truth | median Z | abstained @99% coverage |
|---|---:|---:|
| MBA | 0.262 | 0.27% |
| TNO | 0.379 | 1.36% |
| Trojans | 0.809 | 17.6% |
| **NEO** | **1.029** | **54.7%** |

| velocity band | abstained @99% |
|---|---:|
| ≤0.25 | 0.16% |
| 0.25–0.5 | 1.18% |
| 0.5–1.0 | 16.6% |
| 1.0–2.0 | 85.7% |
| **>2.0** | **100.0%** |

Fast-NEO retention: **58.1%** at 99.9% coverage, **0%** at 99.5% and below. The random-abstention
control was *better* than the policy on ROC/pAUC/F1 at 99.9%.

**Why the standardisation did not rescue it:** per-population normalisation fixes units, but fast
NEOs sit where every population's samples are sparse — including NEO's own relative to typical NEO
spacing — so `min_c z_c` is large there. For this data **low local support is the normal condition
for NEOs, not an error state**, so any gate calibrated to abstain on weak support preferentially
discards them.

## §9.3 FINAL DECISION

1. **REJECT** `support_count` occupancy masking (resolution-dependent, asymmetric, and costly at the
   production resolution: pAUC 0.8867 → 0.9535 with it off).
2. **REJECT** standardized-Z as an automatic abstention gate.
3. **DO NOT** search for another CAL-derived support gate.
4. **Build posterior probabilities from ALL unmasked population densities symmetrically** — same
   estimator, same k, same normalisation, same interpolation, no population exempt, none zeroed.
   `_support_mask_skip` is retired.
5. **Abstain only on explicit TECHNICAL domain failures**: outside map bounds; no sky/magnitude map
   available; invalid observation; non-finite or zero total density. Abstention is NaN, never 0 or 1.
6. **Report Z, nearest/kth-neighbour distance and total density as uncertainty/support METADATA
   only** — surfaced with every score, never gating it.
