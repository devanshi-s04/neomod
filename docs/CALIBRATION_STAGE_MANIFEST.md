# CALIBRATION STAGE — MANIFEST (frozen 2026-08-04, before anything is fitted)

**CAL-only. TEST sealed. Scores come from frozen interpolation variant C** (protocol v1.3 §11).
No map rebuilding, no interpolation changes, no support-policy changes.

## §1. What is being calibrated, and why the prior stays separable

The frozen posterior is

```
P(NEO | X) = rho_NEO(X) / SUM_c rho_c(X)
```

Because the densities carry **absolute object counts**, this already embeds the model's own class
prior. For future LSST prior updates the prior must remain **swappable**, so calibration is applied
to the **likelihood-ratio side**, not to the posterior directly:

```
LR(X) = rho_NEO(X) / SUM_{c != NEO} rho_c(X)          (prior-free)
P(X)  = pi*LR / (pi*LR + (1 - pi))                    (prior pi applied explicitly)
```

The fitted transform acts on `logit(P)` (equivalently `log LR` at the CAL prior). **The CAL class
prior is recorded** so any later evaluation at a different prior is an explicit re-application of
the formula above, not a refit.

## §2. Candidates — exactly four, no additions

| id | transform | free params |
|---|---|---|
| **R** | **raw C posterior** (no transform) — the reference | 0 |
| **P** | Platt / logistic on `logit(P)`: `sigma(a*logit + b)` | 2 |
| **T** | temperature only: `sigma(logit / T)` | 1 |
| **I** | isotonic regression on `P` (monotonic, non-parametric) | free |

All four are **monotonic**, so ranking is preserved by construction.

## §3. Validation — out-of-fold, never fit-and-score on the same rows

- **5-fold cross-validation within CAL.** Every reported calibrated probability is an
  **out-of-fold** prediction: the transform scoring a row was fitted on the other four folds.
- **Folds are grouped by sky map cell** and stratified to balance NEO counts. Rows in one cell share
  a map, so a random row split would leak map-specific structure across folds.
- The raw C probabilities are **preserved unchanged** and reported beside every calibrated variant.

## §4. Metrics — and what may NOT be claimed

**Judged on:** Brier score · log loss · reliability (predicted vs observed, plus ECE).

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
