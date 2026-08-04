# CALIBRATION STAGE — MANIFEST (frozen 2026-08-04, before anything is fitted)

**CAL-only. TEST sealed. Scores come from frozen interpolation variant C** (protocol v1.3 §11).
No map rebuilding, no interpolation changes, no support-policy changes.

## §1. What is being calibrated, and how the prior stays separable

### §1.1 The raw odds already contain a reference prior — corrected 2026-08-04

With **absolute-count** densities,

```
odds_raw(X) = rho_NEO(X) / SUM_{c != NEO} rho_c(X)
```

is **not** a prior-free likelihood ratio: it already carries the class prior implied by the density
normalisation. The prior-free quantity is

```
LR_raw(X) = odds_raw(X) / odds(pi_ref)          odds(pi) = pi / (1 - pi)
```

**`pi_ref` is defined and computed as the NEO prior implied by the same density normalisation:**

```
pi_ref = INT rho_NEO dA dm  /  INT SUM_c rho_c dA dm
```

integrated over the frozen map domain (all 667 cells x 8 magnitude bins x the +-5 deg/day velocity
grid), using the **same stored density arrays** the scorer reads. The uniform kNN `k/(k-1)` bias
cancels in the ratio. Caveat recorded with the value: the 30-deg patches overlap at 10-deg spacing,
so the integral double-counts sky — acceptable for a ratio only to the extent the overlap is
population-independent, which is reported alongside as the per-map spread of `pi_ref`.

> **`pi_CAL` (the CAL class fraction) MUST NOT be substituted for `pi_ref`** unless the two are
> demonstrably identical. Both are computed and compared explicitly; if they differ, that difference
> is itself reported, not reconciled.

### §1.2 Transferring a calibrated probability to a new prior

Calibration is fitted on `s = logit(P_raw)`. From a calibrated posterior under the CAL prior:

```
logLR_cal(X) = logit(P_cal(X)) - logit(pi_CAL)
P_target(X)  = sigmoid( logLR_cal(X) + logit(pi_target) )
```

**Required assumption — LABEL SHIFT (prior probability shift).** This transfer is valid only if the
class-conditional distributions `p(X | y)` are **unchanged** between CAL and the target, and only the
class prior `p(y)` differs. It is **not** valid under covariate shift, under a different detection or
linking pipeline that changes `p(X | y)`, or if the target population model differs from NEOMOD3.
Any LSST prior update must state that this assumption is being made.

## §2. Candidates — exactly four, no additions

| id | transform | free params |
|---|---|---|
| **R** | **raw C posterior** (no transform) — the reference | 0 |
| **P** | Platt: `sigma(a*s + b)`, **a = exp(alpha) > 0 enforced** | 2 |
| **T** | temperature: `sigma(s / T)`, **T = exp(tau) > 0 enforced** | 1 |
| **I** | isotonic regression on `P` (non-decreasing, non-parametric) | free |

R/P/T are **strictly monotonic**, so ranking is preserved exactly. **Isotonic is non-decreasing and
CREATES TIES**, so small ROC / partial-AUC / F1 changes from tied scores are expected; these are
reported **only as invariance diagnostics**, never as improvement.

## §3. Validation — out-of-fold, never fit-and-score on the same rows

- **5-fold cross-validation within CAL.** Every reported calibrated probability is an
  **out-of-fold** prediction: the transform scoring a row was fitted on the other four folds.
- **Folds are grouped by sky map cell** and stratified to balance NEO counts. Rows in one cell share
  a map, so a random row split would leak map-specific structure across folds.
- **ASSERTED before fitting: no repeated parent object crosses folds.** CAL rows carry one synthetic
  object per row; the assertion checks `ObjID` uniqueness and, if any ObjID repeats, that all its
  rows fall in the same fold. The grouped-by-cell design is kept **only if this assertion passes**.
- The raw C probabilities are **preserved unchanged** and reported beside every calibrated variant.

## §4. Metrics — and what may NOT be claimed

**PRIMARY transform-selection metrics: Brier score and log loss.** Reliability/ECE is reported
alongside as a diagnostic.

**Reliability / ECE calculation — frozen now:**
- **15 equal-width bins on [0, 1]** (fixed rule, fixed count; not quantile bins, whose edges would
  move between candidates and make them non-comparable).
- `ECE = SUM_b (n_b / N) * |obs_b - pred_b|`.
- **Uncertainty:** per-bin Wilson 95% intervals on `obs_b`, and a 500-replicate truth-stratified
  bootstrap 95% CI on ECE itself.
- Bins with `n_b = 0` are excluded from the sum and reported as empty.

**Explicitly NOT claimable:** ROC AUC, standardized partial AUC, or best-F1 improvement. All four
candidates are monotonic, so ranking is unchanged up to ties. These will be **computed and reported
as an invariance check** — if they move by more than tie-breaking noise, that is a bug to
investigate, not a result to report.

## §5. Recorded regardless of outcome

- CAL class prior `pi` (NEO fraction on the evaluated rows)
- fitted parameters per fold (a, b, T; isotonic knot count)
- raw vs calibrated side by side, at every reported quantity
- reliability tables for R and each candidate
- per-fold variation, so a transform that is unstable across folds is visible

## §6. Out of scope here

Operating-threshold selection. That is a **separate** step after a transform is chosen, and it is
also CAL-only. `MODEL_SEAL.json` is written only after both are complete.

**Pause after the comparison and before selecting the final transform.**
