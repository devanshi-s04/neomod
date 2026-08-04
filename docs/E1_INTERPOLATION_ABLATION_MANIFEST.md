# E1 INTERPOLATION ABLATION — MANIFEST (frozen before any variant is scored)

**CAL-only. TEST sealed. No map rebuilding — every candidate uses the sealed 667-map grid
`prob_maps_grid_neomod3_GEN_final` (MAP_BUILD_SEAL `c30cc29d24ead70d…`).**

`MAP_BUILD_SEAL` freezes **map generation**, not the scoring/interpolation implementation, and
`MODEL_SEAL` is not yet written — so density-first interpolation is permitted at this stage. Every
candidate is scored **end-to-end**; no term is substituted into another candidate's path.

## Variants (exactly four — no additions)

| id | sky cell | magnitude bin | velocity | posterior formed |
|---|---|---|---|---|
| **A** | nearest | nearest | existing sealed within-map interpolation | **pixelwise, then interpolated** (sealed scorer) |
| **B** | nearest | nearest | density-first bilinear | after interpolation |
| **C** | nearest | **density-first linear between bins** | density-first bilinear | after interpolation |
| **D** | **density-first bilinear across neighbouring cells** | density-first linear between bins | density-first bilinear | after interpolation |

A is `ProbMapSet.score_visible` unchanged. B/C/D interpolate each population's density and then
normalise, so `P = rho_NEO / SUM_c rho_c` is formed **after** interpolation.

## Rows and mask

- **Full sealed CAL v2** (`CAL_DATASET_SEAL` `847bb72265479154…`), all 667 cells — not the 16-cell
  E0 subset, which existed only because only 16 maps did. This raises the NEO count from 181 to
  ~5,000 and removes the sample-size limitation that made the E0 per-center rule undefined.
- **Identical rows for every candidate.** The evaluation mask is the **intersection of all four
  candidates' technical-valid masks**, so no candidate can benefit by abstaining on hard rows.
- Technical-valid per `MAP_BUILD_SEAL.FROZEN_POLICY`: in-bounds, finite, non-zero total density.
- **Class-conditional coverage reported** (overall and NEO separately), since every E0 abstention
  was a NEO.

## Metrics (all raw posterior — NO calibration fitted at this stage)

ROC AUC · standardized partial ROC AUC (FPR ≤ 0.01) · best F1 · Brier · log loss · reliability
(predicted vs observed by bin) · fast-NEO recovery (|v| > 2) · runtime.

Selection uses the rule already frozen in `EVALUATION_PROTOCOL.md`: compare at **matched
completeness or matched FPR**, never at self-selected thresholds; prevalence-dependent quantities
reported with the prior stated.

## Explicitly out of scope here

Calibration fitting or selection — a **separate** bounded CAL-only stage afterwards, preserving the
likelihood-ratio/prior separation and always reporting raw and calibrated side by side.
