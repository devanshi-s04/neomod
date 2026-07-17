# Hyak Agent Handoff: Pure-S3M Blackbox Products

## Mission

Produce the full-grid evidence bundle for the pure-S3M advisor deck.

Hyak has the full data products. The task is not to write the final presentation.
The task is to generate and verify the full set of comparable products needed by
the Mac-side presentation notebook.

The scientific question:

> When the same VDP blackbox is given pure-S3M benchmark tracklets versus
> pure-S3M Sorcha tracklets, why does performance change, and what do we need
> to fix?

## Experiment Matrix

Produce products for both cells (GMM branch only; kNN branch dropped):

| Input type | Map branch | Comparator |
| --- | --- | --- |
| Benchmark pure-S3M | GMM 667 maps | Digest2 |
| Sorcha pure-S3M | GMM 667 maps | Digest2 |

Use the same direction bins everywhere:

- full sky
- 0-20 deg from antisun
- 20-40 deg from antisun
- 40-70 deg from antisun
- 70-110 deg from antisun
- 110-141 deg from antisun

Also compute magnitude-bin versions wherever the columns/maps support it.

## Required Inputs To Verify

Canonical Hyak path:

```text
/mmfs1/gscratch/dirac/ds2004/sorcha
```

Verify these exist or identify their canonical replacements:

```text
neomod/
S3Mdata/ or raw S*.s3m files
prob_maps_grid_s3m/              # pure-S3M GMM 667 maps (complete)
outputs/phase2_s3m/sorcha_comparison_s3m.parquet
outputs/benchmark_tracklets_s3m/tracklets_benchmark.parquet
outputs/phase2_benchmark_s3m/   # benchmark VDP + digest2 shards
digest2/digest2                  # digest2 binary
```

kNN/K|M maps (prob_maps_grid_knn_s3m/) are not used for this deck.

## Grid Geometry

GMM grid (only branch used for this deck):

```text
reference epoch: 2026-01-01T00:00:00
longitude step: 10 deg
sun exclusion: 40 deg
latitude base: 0,1,2,3,4,5,8,12,18,25,35,50
total maps: 29 longitudes x 23 latitudes = 667
velocity grid: (-2.0, 2.0) deg/day, 401x401
```

Manifest (grid_manifest_gmm.csv):

```text
index
delta_lon_from_antisun_deg
lat_deg
filename
map_branch
```

Validation:

- exactly 667 map files in prob_maps_grid_s3m/
- no zero-byte files
- filenames and metadata agree

## Benchmark Products

Benchmark pure-S3M: 475,027 synthetic tracklets (S3M objects propagated to
2026-01-01, 30-min 2-det X05 tracklets) scored against the GMM 667-map grid.

Pipeline:
  gen_benchmark_tracklets_s3m.py (COMPLETE — outputs/benchmark_tracklets_s3m/)
  sorcha_phase2.py score-vdp     (GMM maps, support-mask-min 1)
  sorcha_phase2.py run-digest2   (96-task array)
  sorcha_phase2.py combine       -> outputs/phase2_benchmark_s3m/benchmark_comparison_s3m.parquet

Produce:

```text
benchmark_metrics_by_direction_bin.csv
benchmark_metrics_by_mag_bin.csv
benchmark_roc_gmm_<bin>.png
```

Required metrics:

```text
N_total
N_NEO
threshold_best_f1
completeness_best_f1
contamination_best_f1
F1_best
optional AUC
```

Digest2 must be evaluated on the same benchmark tracklets used by VDP.

## Sorcha Products

Sorcha pure-S3M: outputs/phase2_s3m/sorcha_comparison_s3m.parquet (648,908 rows,
GMM branch, COMPLETE). Contains P_NEO_vdp and P_NEO_d2.

Produce:

```text
sorcha_metrics_by_direction_bin.csv
sorcha_metrics_by_mag_bin.csv
sorcha_roc_gmm_<bin>.png
```

## Map Diagnostic Figures

For GMM maps, generate panels for the standard direction bins and representative
mag bins:

- log density maps in ecliptic rate space
- P(NEO) maps
- mask OFF panels
- mask ON panels (support-count mask, support_mask_min=1)

Suggested output names:

```text
fig_density_gmm_<direction_bin>_<mag_bin>.png
fig_pneo_mask_off_gmm_<direction_bin>_<mag_bin>.png
fig_pneo_mask_on_gmm_<direction_bin>_<mag_bin>.png
```

## Sanity Checks

Produce one compact validation note:

```text
pure_s3m_blackbox_hyak_validation.txt
```

It should include:

- grid counts for both branches
- map directories used
- Sorcha parquet row count
- benchmark row counts
- whether Digest2 scores are reused or rerun
- known failures or caveats
- exact commands/scripts used

Also produce object-level rate sanity material if available:

```text
jpl_horizons_rate_check.csv
jpl_horizons_rate_check.png
```

This should compare selected object ecliptic rates from the pipeline to JPL
Horizons or another trusted reference.

## Bundle For Mac

Create a transferable bundle directory:

```text
outputs/presentation_pure_s3m_blackbox/
```

Minimum contents:

```text
README.txt
pure_s3m_blackbox_hyak_validation.txt
grid_manifest_gmm.csv
benchmark_metrics_by_direction_bin.csv
benchmark_metrics_by_mag_bin.csv
sorcha_metrics_by_direction_bin.csv
sorcha_metrics_by_mag_bin.csv
figures/
```

If feasible, include compact scored parquets. If the full scored parquets are
too large, include metrics and figures first.

## Do Not Do Yet

- Do not add hybrid catalog results.
- Do not change the scientific conclusion to make VDP look better.
- Do not mix exact old offsets with direction-bin claims unless labelled as
  historical context.
- Do not hide the Sorcha underperformance relative to Digest2.

The whole point is to make the remaining gap clear enough for Mario/Zeljko to
help fix.

