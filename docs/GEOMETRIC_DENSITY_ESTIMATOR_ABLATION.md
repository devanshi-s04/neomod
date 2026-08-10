# Geometric Density-Estimator Ablation

**Status:** Proposed experiment. Freeze and hash this protocol before inspecting any new results.

**First evaluation domain:** the controlled geometric benchmark only. Sorcha is explicitly deferred.

## Purpose

The raw VDP probabilities on the geometric benchmark are systematically underconfident at decision-relevant scores. Platt calibration can repair the numerical interpretation of those scores, but it cannot identify why the raw density ratio is underconfident.

This experiment goes back to the density construction itself. It tests whether the Bayesian density estimator, and specifically its treatment of the NEO population, is responsible for the probability distortion.

The experiment compares three map-building cases while holding the population realizations, physical weights, map geometry, scoring rows, interpolation, and all other choices fixed.

## Scientific Question

Which density-estimation choice causes the raw VDP posterior to differ from the empirical NEO frequency?

The experiment must separate:

1. A problem introduced by Bayesian density estimation for all populations.
2. A problem introduced specifically by Bayesian density estimation of the sparse NEO population.
3. A more fundamental mismatch in population normalization, priors, map discretization, or the simulated test distribution.

## Scope and Non-Goals

This is a mechanism test, not a search for the best-looking ROC curve.

- Use the geometric benchmark, where positions and velocities are generated without Sorcha measurement noise.
- Develop and select using an independent geometric `CAL` realization.
- Keep the already inspected geometric `TEST` out of model development.
- After the mechanism and implementation are frozen, evaluate once on a fresh geometric `TEST2` realization.
- Do not use Sorcha in this first experiment.
- Do not apply Platt scaling or any other post-hoc calibration.
- Do not apply Gaussian smoothing or any additional image/map smoothing.
- Do not change the velocity domain, grid spacing, sky cells, magnitude bins, population model, physical weights, interpolation, or support policy between cases.
- Do not assign unsupported or failed rows a score of 0 or 1. Record them as `NaN` with an explicit reason.

## Frozen Population Inputs

All three cases use exactly the same parent realizations and rows:

- **NEO density population:** independent NEOMOD3 `GEN` realization with its frozen seed, orbit cache, physical weights, and provenance hashes.
- **MBA, TNO, and Trojan density populations:** uncloned S3M `GEN` partitions with their frozen split corrections and provenance hashes.
- **Calibration rows:** independent geometric `CAL` parents, with no parent-orbit identity shared with `GEN`.
- **Final evaluation rows:** a fresh geometric `TEST2`, created only after the estimator and all interpretation rules are frozen.

The legacy S3M NEO population must never enter the NEO density tree in any case.

## Frozen Map and Scoring Choices

The following must be identical in all three cases:

- Sky-map centers and nearest-sky-map assignment.
- Magnitude-bin edges.
- Velocity domain and grid spacing.
- Population definitions and physical weights.
- Observer geometry and geometric-tracklet construction.
- Raw-density units.
- Density-first bilinear interpolation in velocity.
- Density-first linear interpolation between magnitude bins.
- Posterior normalization.
- Common CAL rows and common-scorable-row definition.
- Numerical precision and zero-density behavior.

Interpolation remains enabled because it is part of evaluating the same continuous query against a discretized map. It is not a density estimator and is not the smoothing mechanism under test. Population densities must be interpolated first and normalized into a posterior afterward.

No occupancy mask, support mask, or automatic support-based abstention may be introduced in this ablation. Technical failures and zero-total-density queries remain explicit and are reported separately.

## The Three Cases

| Case ID | NEO density | MBA/TNO/Trojan density | Extra Gaussian smoothing |
|---|---|---|---|
| `HIST_ALL` | Direct weighted histogram | Direct weighted histogram | Off |
| `BAYES_ALL` | Current Bayesian density estimator | Current Bayesian density estimator | Off |
| `BAYES_NONNEO` | Direct weighted histogram | Current Bayesian density estimator | Off |

### Case 1: `HIST_ALL`

No Bayesian density estimation is applied to any population. For population `c` and velocity pixel `(i, j)`, compute the direct physical density

```text
rho_c(i, j) = sum(weight_n for samples n in pixel i,j) / pixel_area
```

where

```text
pixel_area = delta_vlambda * delta_vbeta.
```

Use the same grid edges for every population. Do not fill empty pixels, smooth counts, borrow neighboring samples, or add pseudocounts.

This case shows what the finite simulated population directly says before Bayesian density estimation.

### Case 2: `BAYES_ALL`

Apply the latest NEOMOD3 VDP Bayesian density estimator, unchanged, to NEOs, MBAs, TNOs, and Trojans. Disable the separate Gaussian map smoothing.

This isolates the behavior of the current Bayesian estimator without the later Gaussian convolution.

### Case 3: `BAYES_NONNEO`

Use the direct weighted histogram density for NEOs and the unchanged Bayesian density estimator for MBAs, TNOs, and Trojans. Disable Gaussian smoothing.

This is a mechanism diagnostic. It tests whether smoothing the sparse NEO population is what suppresses the raw NEO odds. Because it mixes estimator types, it is interpretable only if the density-integral and unit checks below pass.

## Posterior Calculation

For every query tracklet, interpolate the four raw physical densities using the same frozen interpolation path and calculate

```text
rho_nonNEO = rho_MBA + rho_TNO + rho_Trojans
P_NEO = rho_NEO / (rho_NEO + rho_nonNEO)
```

If the denominator is zero or any required density is invalid, return `NaN` and a reason code. Never silently substitute zero.

No Platt transformation is applied. The experiment evaluates the raw density-derived posterior.

## Implementation Architecture

Do **not** create three independent full copies of the VDP implementation. Three copies could drift and create differences unrelated to the estimator.

Instead:

1. Preserve `src/velocity_density_pipeline_neomod_clone_only.py`, all sealed maps, and both existing seals unchanged.
2. Create one dedicated ablation implementation from the latest NEOMOD3 version, for example:

   `src/velocity_density_pipeline_neomod_density_ablation.py`

3. Add one explicit, validated mode parameter:

   `density_mode = hist_all | bayes_all | bayes_nonneo`

4. Route each population to either the direct histogram or the existing Bayesian estimator according to the table above.
5. Write three separate immutable map products and manifests. These are the three experimental variants; the source logic remains shared.

Thin mode-specific launch configurations are encouraged. Full duplicated source modules are not.

## Pre-Run Validation Gates

The following checks must pass before interpreting any classifier result.

### 1. Population and split provenance

- NEO density rows come exclusively from NEOMOD3 `GEN`.
- Non-NEO density rows come exclusively from S3M `GEN`.
- No `CAL`, old `TEST`, or future `TEST2` parent identity enters a density map.
- Report input row counts, unique parent counts, seed/cache hashes, and source-file SHA256 values.

### 2. No extra smoothing

- Assert Gaussian sigma is zero/disabled for every population and case.
- Assert no alternate KDE, pseudocount, hole filling, or post-map filter is invoked.
- Store this assertion in every map manifest.

### 3. Direct-histogram correctness

For a synthetic unit-test sample with known coordinates and weights:

- Verify exact pixel assignment, including edge behavior.
- Verify the density integral equals the sum of physical sample weights.
- Verify empty pixels remain exactly zero.

### 4. Bayesian-path fidelity

- `BAYES_ALL` and the non-NEO branches of `BAYES_NONNEO` must call the exact current Bayesian density routine with identical arguments.
- Do not modify its normalization as part of this ablation.
- Record the routine name, source commit, parameters, and code hash.

### 5. Density units and integrals

For every population, sky map, magnitude bin, and case, report:

```text
integral_c = sum(rho_c * pixel_area)
expected_c = sum(physical input weights represented by that map/bin)
ratio_c = integral_c / expected_c
```

The estimator cases must use compatible physical-density units. A large or systematic integral mismatch is an implementation/normalization finding and blocks posterior interpretation.

### 6. Scoring identity

- Score exactly the same CAL tracklet UIDs in all cases.
- Assert identical positions, velocities, magnitudes, map choices, and interpolation inputs.
- Compare scientific metrics on the common finite-score subset.
- Report coverage and failure composition separately so a case cannot win by dropping difficult rows.

## Execution Stages

### Stage A: Read-Only Inventory

Before editing or running:

- Identify the exact latest NEOMOD3 VDP source commit and map-generation entry point.
- Identify the Bayesian density function and the separate Gaussian-smoothing call.
- Identify frozen GEN and geometric CAL artifacts and their hashes.
- Identify the existing 16 representative E0 sky centers and all magnitude bins.
- Produce a short implementation plan and file list.

### Stage B: Unit Tests and Representative Pilot

Run the three cases on the same 16 E0 sky centers and all magnitude bins.

The pilot is for implementation validation only. It must establish:

- The three estimator modes dispatch correctly.
- Gaussian smoothing is absent.
- Density integrals and units pass.
- The same CAL rows are scored.
- Empty-density and invalid-query behavior is explicit.
- Representative density and posterior maps differ only where the estimator choice predicts.

Do not use this limited pilot to declare a scientific winner.

### Stage C: Full Geometric CAL Comparison

After reviewing the pilot gates, build all frozen sky maps for all three cases and score the full geometric CAL set.

Use CAL to diagnose the mechanism and, if necessary, choose one estimator policy. Do not inspect the old TEST again and do not create TEST2 yet.

### Stage D: Freeze and Final Geometric Evaluation

Once the estimator policy and interpretation are fixed:

- Write a new ablation/model seal with source hashes, map hashes, configuration, and decision rationale.
- Generate a fresh independent geometric `TEST2` parent realization.
- Score `TEST2` once with the frozen cases or selected policy and deterministic digest2.
- Report the final result without further tuning.

Sorcha transfer testing begins only after this geometric mechanism test is closed.

## Required Metrics

Compute all metrics on CAL for diagnosis and on TEST2 only for final evaluation:

- ROC AUC.
- Standardized partial AUC for the frozen low-FPR range.
- Completeness versus contamination.
- Brier score and log loss for raw probabilities.
- Reliability table with fixed, preregistered score bins and Wilson intervals.
- Calibration slope and intercept as diagnostics only; do not apply them to scores.
- Coverage, zero-density rate, non-finite rate, and reason codes.
- Metrics by population, velocity band, magnitude bin, and sky direction.
- Per-object score differences among the three cases.

The comparison must include both all valid rows per case and the common finite-score subset.

## Mechanism Diagnostics

For each score bin and selected phase-space region, compare empirical odds with density odds:

```text
observed_odds = observed_NEO_fraction / (1 - observed_NEO_fraction)
density_odds = rho_NEO / rho_nonNEO
delta_log_odds = log(observed_odds) - log(density_odds)
```

Also produce:

- Per-population density-integral tables.
- Raw NEO and non-NEO density maps for representative fields and magnitude bins.
- Raw posterior maps for all three cases using identical color scales.
- One-dimensional slices through representative dense and sparse velocity regions.
- Histograms of per-object score changes: `HIST_ALL -> BAYES_ALL` and `BAYES_ALL -> BAYES_NONNEO`.
- Reliability curves for all three raw cases on one plot.
- Lists of objects with the largest score changes, including velocity, magnitude, sky map, and all four population densities.

## Interpretation Rules

- If `HIST_ALL` is calibrated, `BAYES_ALL` is underconfident, and `BAYES_NONNEO` returns toward calibration, the NEO Bayesian estimator is the leading mechanism.
- If `BAYES_ALL` and `BAYES_NONNEO` behave similarly, smoothing the NEO density is not the principal cause.
- If `HIST_ALL` is also underconfident, investigate physical normalization, priors, discretization, and population mismatch before calibration.
- If `BAYES_NONNEO` changes the posterior but its density integrals or units differ from the other cases, treat the result as an estimator-unit artifact, not evidence about NEO smoothing.
- If a Bayesian case improves ranking but worsens probability calibration, report those as separate effects.
- If the effect changes strongly by velocity, magnitude, or sky direction, do not summarize it as one global correction.

## Deliverables

Write all outputs under:

`outputs/geometric_density_estimator_ablation/`

Required products:

- `protocol_snapshot.md` and SHA256.
- `inventory.json`.
- `maps_hist_all/`, `maps_bayes_all/`, and `maps_bayes_nonneo/`.
- One manifest per map variant with source/config/input hashes.
- `density_integral_audit.parquet`.
- `cal_scores_hist_all.parquet`.
- `cal_scores_bayes_all.parquet`.
- `cal_scores_bayes_nonneo.parquet`.
- `cal_common_rows.parquet`.
- `cal_metrics.csv`.
- `cal_reliability_bins.csv`.
- `cal_coverage_and_failures.csv`.
- `cal_per_object_deltas.parquet`.
- Diagnostic figures described above.
- `PILOT_REPORT.md` followed by `CAL_REPORT.md`.

Large map products should remain on Hyak. Commit code, launch files, manifests, compact tables, figures, and reports to Git; do not commit large map archives or scored parquet files unless already permitted by repository policy.

## Stop Rule

After the three raw cases are run and explained, pause for advisor review.

Do not add Platt scaling, reintroduce Gaussian smoothing, tune another support rule, alter grid limits, start Sorcha optimization, or inspect a final TEST2 until the density-estimator mechanism has been interpreted and the next step is explicitly approved.
