# GMM Normalisation Bug — Diagnostic Note for Hyak Claude
# Written by Arnor Claude, 2026-06-03

## Summary

The advisor identified a bug: GMM P(NEO) is much weaker than S3M kNN P(NEO),
especially at vlam ≈ −0.5 where only NEOs exist and MBAs don't.
We ran a diagnostic on Arnor and pinpointed the exact problem.

---

## Diagnostic run (on Arnor, notebook: sorcha_gmm_s3m_singleepoch_comparison.ipynb)

Maps used:
- S3M: `prob_maps/prob_maps_2026-05-09T22_neocp.npz`  (center lon=229°, grid ±0.8)
- GMM: `prob_maps_gmm/prob_maps_2026-05-01_antisun.npz` (center lon=220°, grid ±2.0)

Mag bin tested: `mag22`

### Table 1 — Raw density at vbeta=0

| vlam  | S3M_MBA    | GMM_MBA    | S3M_NEO    | GMM_NEO    |
|-------|-----------|-----------|-----------|-----------|
| -0.50 | 4.830e+00 | 4.309e+01 | 1.054e+02 | 3.923e+01 |
| -0.30 | 8.797e+02 | 5.254e+03 | 1.865e+03 | 5.224e+02 |
| -0.20 | 1.321e+07 | 1.688e+07 | 1.255e+03 | 3.989e+02 |
|  0.00 | 5.697e+01 | 2.343e+02 | 9.807e+00 | 6.285e+00 |
|  0.20 | 3.343e+00 | 2.412e+01 | 2.995e+00 | 6.284e-01 |

### Table 2 — Grid integrals (total weight = sum × cell area)

| pop     | S3M integral | GMM integral | GMM/S3M ratio |
|---------|-------------|-------------|---------------|
| MBA     | 9.8783e+04  | 9.8277e+04  | **0.995**     |
| NEO     | 7.9556e+01  | 2.5642e+01  | **0.322**     |
| TNO     | 1.4886e+04  | 1.0030e+04  | 0.674         |
| Trojans | 5.3493e+02  | 2.4245e+03  | 4.532         |

### Table 3 — P(NEO) at vbeta=0

| vlam  | P(NEO) S3M | P(NEO) GMM | ratio S3M/GMM |
|-------|-----------|-----------|---------------|
| -0.50 | 0.9543    | 0.4568    | 2.09          |
| -0.30 | 0.6793    | 0.0902    | 7.53          |
| -0.20 | 0.0001    | 0.0000    | 4.02          |
|  0.00 | 0.0448    | 0.0005    | 95.86         |
|  0.20 | 0.4324    | 0.0191    | 22.61         |

---

## Two bugs identified

### Bug 1 — GMM NEO density is 3× underscaled (primary bug)

The integrals table is the smoking gun:
- MBA integral ratio: **0.995** — GMM MBA is correctly scaled relative to S3M
- NEO integral ratio: **0.322** — GMM NEO density is only 1/3 what S3M produces

The MBA training produces the correct total weight; the NEO training does not.
GMM NEO density needs to be ~3.1× larger to match S3M's relative scaling.

Effect at vlam=−0.5:
- S3M NEO/MBA density ratio = 105 / 4.8 = **22** → P(NEO) ≈ 0.95 (correct)
- GMM NEO/MBA density ratio = 39 / 43 = **0.91** → P(NEO) ≈ 0.46 (wrong)

### Bug 2 — GMM MBA Gaussian tails bleed into vlam=−0.5 (secondary bug)

The GMM is a smooth Gaussian mixture. At vlam=−0.5 where no MBA clones exist in S3M:
- S3M MBA density: 4.8 (near zero — no clones there)
- GMM MBA density: **43** (9× higher — tail of the Gaussian)

This is a fundamental property of the GMM (smoothness). The tails raise MBA density
in the NEO-only region, further suppressing P(NEO). Bug 1 is fixable; Bug 2 is
partially intrinsic to GMM but made worse by Bug 1.

---

## Where to look in the code

File: `src/velocity_density_pipeline_gmm.py`

Look for how `density_maps_downweighted_raw` is computed per population. Specifically:

1. **Check clone_factor per population** — is the same clone_factor applied to NEO
   training objects as to MBA training objects? If NEO uses a lower clone_factor, 
   the total NEO weight will be lower.

2. **Check per-population normalization** — is any division by N_training_objects or
   N_components applied AFTER the density is computed? If so, check whether this
   is applied consistently across populations.

3. **Quick diagnostic to run on Hyak** — add these print statements to the map
   generation code and re-run on one map to see the raw integrals:

```python
for pop in ["MBA", "NEO", "TNO", "Trojans"]:
    d = results["mag22"]["density_maps_downweighted_raw"][pop]
    print(f"{pop}: sum={d.sum():.4e}  max={d.max():.4e}  nonzero={np.count_nonzero(d)}")
```

The expected result (if normalisation is correct) is that the ratio
`NEO.sum() / MBA.sum()` should match the S3M survey population ratio
(approx 268K NEOs / 13.9M MBAs ≈ 0.019). If GMM gives a much lower ratio,
the NEO density is being produced at the wrong scale.

---

## The fix

After computing each population's density map, ensure the total weight
(integral over the velocity grid) is consistent with the survey-based
population priors — i.e. the S3M N_pop values, not the training set sizes.

The simplest normalisation fix:

```python
# After computing density_maps_downweighted_raw for each population:
# Normalise each map to integrate to 1, then rescale by the desired prior weight.

N_survey = {"MBA": 13_900_000, "NEO": 268_000, "TNO": 49_000, "Trojans": 180_000}

for pop in pops:
    d = density_maps_downweighted_raw[pop]
    total = d.sum()
    if total > 0:
        density_maps_downweighted_raw[pop] = d / total * N_survey[pop]
```

This ensures the MBA:NEO density ratio reflects the actual survey population
ratio, regardless of how many training objects or GMM components were used
per population.

---

## Expected outcome after fix

At vlam=−0.5 with corrected NEO density (×3.1) and same MBA density:
- GMM NEO/MBA ratio → 39×3.1 / 43 ≈ 2.8  (was 0.91)
- GMM P(NEO) at vlam=−0.5 → ~0.73  (was 0.46, S3M is 0.95)

Won't fully close the gap with S3M (Bug 2, MBA tails, remains) but should
substantially improve P(NEO) in the NEO-only velocity region and raise F1
significantly above the current 0.837.

---

## Current science status (as of 2026-06-03)

- S3M kNN: F1=0.847 (single epoch), F1=0.740 (Sorcha 2yr)
- GMM (with bug): F1=0.842 (single epoch, same objects), F1=0.837 (Sorcha 2yr)
- digest2: F1=0.655 (single epoch), F1=0.837 (Sorcha 2yr)
- Target after fix: F1 > 0.857 on Sorcha 2yr

Good luck! — Arnor Claude
