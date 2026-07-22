# mag245 tests and LSSTCam field analysis

**Updated:** 2026-07-21
**Scope:** n-body benchmark/Sorcha comparisons restricted to the mag < 24.5 products, with the night-61642 identity test and the exact LSSTCam field experiment.

This note records the work completed in the mag245 notebooks so the analysis can be resumed without reconstructing the notebook history from cached figures.

## Inputs and conventions

The primary local files are:

| File | Use |
|---|---|
| `outputs/phase2_benchmark_s3m_nbody_mag245/benchmark_comparison_s3m_nbody_mag245.parquet` | n-body benchmark tracklets with benchmark VDP and digest2 scores |
| `sorcha_comparison_case1_nbody_Vband.parquet` | Sorcha case 1 tracklets with V-band VDP and digest2 scores |
| `../baseline_v5.0.0_2yrs.db` | Rubin pointing database used to recover visit centers, rotations, and visit IDs |
| `/opt/anaconda3/lib/python3.13/site-packages/sorcha/modules/data/LSST_detector_corners_100123.csv` | Bundled LSSTCam detector geometry, loaded through Sorcha's `Footprint()` class |

The pointing database is intentionally one directory above `neomod/`:

```text
/Users/devanshisingh/Downloads/research/NEO_probability/baseline_v5.0.0_2yrs.db
```

The benchmark and Sorcha comparison use n-body kinematics. VDP scores are compared in Johnson V: benchmark `P_NEO_vdp` is compared with Sorcha `P_NEO_vdp_Vband`. Sorcha's original `mean_mag` and `P_NEO_vdp` columns remain available for the before/after magnitude-band comparison.

Rates use `vlam` and `vbeta` in deg/day, with no extra cos(beta) factor. The identity key is benchmark `s3m_objid ==` Sorcha `ObjID`.

## Main notebooks

| Notebook | Purpose |
|---|---|
| `mag245_nbody_benchmark_vs_sorcha_roc.ipynb` | Full-sky benchmark-versus-Sorcha VDP/digest2 operating curves, conventional ROC/AUC, elongation and magnitude breakdowns, and same-object score comparisons |
| `mag245_nbody_scoring_consistency_night61642.ipynb` | Same-night scoring consistency and per-object benchmark/Sorcha score diagnostics |
| `nbody_tracklet_velocity_benchmark_night61642_identity.ipynb` | N-body velocity and offset-vector comparisons for case 1, case 2, and other available Sorcha cases |
| `LSSTcam_tracklets_scoring.ipynb` | Exact LSSTCam visit-pair selection, Field A visualization, dense-field selection, and field-level score curves |
| `lsstcam_tracklets_scoring.py` | Reusable loader, visit matcher, LSSTCam geometry, field-selection, plotting, and score-curve code used by the LSSTCam notebook |

## Magnitude-band correction

The VDP maps are V-band maps, while the original Sorcha comparison used the observed LSST-band tracklet magnitude for the VDP magnitude-bin lookup. That caused a systematic approximately 0.6 mag offset and changed many magnitude-bin assignments.

The corrected Sorcha output is add-only. It preserves the original columns and adds:

- `mean_mag_V` — converted Johnson-V tracklet magnitude;
- `mag_bin_label_Vband` — V-band magnitude bin;
- `P_NEO_vdp_Vband` — VDP score evaluated with the V-band magnitude.

The conversion was validated against the benchmark: the median Sorcha-minus-benchmark magnitude residual changed from approximately `-0.606` mag to `+0.007` mag. The residual scatter is photometric scatter, not a remaining band offset. See `SORCHA_MAG_BAND_FIX_APPLIED.md` for the derivation and production-side details.

## Night-61642 identity test

Night 61642 is the busiest Sorcha night in the case products. The scoring consistency notebook reports:

- 1,762 identity pairs;
- 794 benchmark-truth NEOs;
- benchmark/Sorcha VDP score difference median approximately `+0.0000`;
- absolute VDP score difference median `0.0002`;
- 6.0% of pairs have `|Delta P| > 0.05` in the reported comparison;
- the mag245 ROC notebook reports 794 same-night, benchmark-truth NEOs in the matched score table;
- VDP absolute score agreement within 0.1 is 90.1%; digest2 agreement within 0.1 is 88.4% in that same table.

The fair interpretation is that the n-body benchmark and Sorcha produce nearly identical scores for the same object when the geometry, epoch, and magnitude convention are aligned. Population-level curves answer a different question because the benchmark and Sorcha detected sets have different base rates and selection functions.

## Plot and notebook fixes

The night-61642 identity notebook contains the following fixes:

1. Quiver panels now set explicit symmetric limits. Matplotlib's quiver autoscale sees the common `(0, 0)` anchor points but not the arrow tips, so relying on autoscale clipped the vectors into a tiny origin window.
2. Section 4 offset-arrow plots use limits from the actual `dvl`/`dvb` data.
3. RMS labels are computed from the points inside the displayed plot limits. The notebook reports in-view component RMS and vector RMS, and separately reports the explicitly defined core RMS where that zoom is intentional.
4. The retained notebook no longer includes the removed test1 magnitude sections or their figure.
5. The historical `ROC` label is retained where necessary for continuity, but the completeness-versus-contamination plots are identified as precision-recall-derived operating curves. Conventional TPR-versus-FPR ROC and AUC are plotted separately in the mag245 comparison notebook.

## Exact LSSTCam experiment

The original circular search idea was replaced by the actual camera geometry. Sorcha's built-in footprint has 189 detector polygons, including chip gaps. The executed notebook reports:

- active focal-plane area: 9.49 deg^2 in the tangent-plane calculation;
- camera width: 3.51 degrees;
- all 4,316 night-61642 case-1 Sorcha tracklets matched uniquely to visits;
- both detections of all retained tracklets pass their respective rotated detector footprints;
- visit matching uses `observationStartMJD + visitTime/2`, which agrees with the Sorcha detection midpoint to better than 0.01 seconds and preserves the observing filter.

### Field A

Fields are actual visit pairs `(observationId0, observationId1)`, not arbitrary RA/Dec cells. Candidate pairs are ranked by total tracklet count only. Field A is the densest pair:

```text
observationId0 = 373765
observationId1 = 373817
center         = (RA 315.682518 deg, Dec -13.277114 deg)
rotSkyPos      = 144.966521 deg / 120.543191 deg
Sorcha tracklets = 25
```

Field A contains 5 NEOs, 19 MBAs, and 1 TNO in the Sorcha population labels. Its benchmark identity matches contain 5 NEOs and 4 MBAs, so Field A is a descriptive crowded-field case study rather than a sufficiently powered standalone ROC sample.

### Dense-field aggregate

Starting with Field A, fields are added in density order when their centers are at least 3.5 degrees from all previously selected centers. This avoids overlapping camera fields while keeping selection independent of scores and labels. The prespecified target is at least 75 matched NEOs and 75 matched non-NEOs.

The exact-camera run reaches the target with:

```text
27 fields
547 Sorcha tracklets
199 identity-matched objects
75 benchmark-truth NEOs
124 benchmark-truth non-NEOs
```

The earlier 20-field estimate came from a circular proxy and is superseded by these exact visit-pair counts. The 1.75-degree spherical circle remains useful only as a sensitivity check, not as the primary camera model.

The exact-camera figures are saved under `Figures/`:

- `LSSTcam_01_geometry.png`
- `LSSTcam_02_field_A_tracklets.png`
- `LSSTcam_03_selected_fields.png`
- `LSSTcam_04_dense_fields_precision_recall.png`
- `LSSTcam_05_dense_fields_roc.png`
- `LSSTcam_06_zero_diagnostics.png`

The selected-field tables are:

- `outputs/LSSTcam_top30_field_candidates.csv`
- `outputs/LSSTcam_selected_fields.csv`
- `outputs/LSSTcam_selected_tracklets.csv`

### Dense-field score results

For the 199 matched objects in the selected fields:

| Dataset | Score | Best F1 | Conventional ROC AUC |
|---|---|---:|---:|
| Benchmark | VDP | 0.847 | 0.879 |
| Benchmark | digest2 | 0.837 | 0.919 |
| Sorcha | VDP | 0.847 | 0.888 |
| Sorcha | digest2 | 0.837 | 0.914 |

The completeness-versus-contamination curves use benchmark population truth in both panels and field-cluster bootstrap bands. Field A itself is not used for threshold tuning.

### VDP-zero cohort

Within the selected dense-field sample, three objects receive Sorcha VDP exactly zero. All three are benchmark-truth NEOs and have Sorcha digest2 scores near one. This is the first focused cohort for the next investigation: compare their velocity location, magnitude bin, probability-map assignment, focal-plane position, detector edge distance, and benchmark/Sorcha score differences.

## Reproducibility checks

The executed LSSTCam notebook passed these checks:

- no notebook cell execution errors;
- all code cells executed;
- exact database path resolved to the project-level copy;
- visit midpoint matching remained below 0.01 seconds;
- both detections passed the exact camera footprint;
- selected field labels and membership tables reconcile with the plotted sample counts;
- no duplicate tracklet or identity object entered the selected aggregate.

## Next tests

The immediate next analysis is the VDP-zero/digest2-positive cohort, beginning with the three dense-field NEOs and then expanding to the full night-61642 matched set. After that, repeat the exact-camera field selection for Sorcha case 2 and case 3, keeping the same density-only selection rule and benchmark-truth scoring convention. Any comparison across cases should report both total tracklets and identity-matched objects because the cases have different detection and population mixes.

Related context: `NBODY_REFERENCE.md`, `SORCHA_MAG_BAND_FIX_APPLIED.md`, `pipeline_paper.md`, and `SORCHA_V5_PIPELINE.md`.
