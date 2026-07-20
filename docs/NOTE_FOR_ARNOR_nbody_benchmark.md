# Note for Arnor — the n-body BENCHMARK (pure-S3M tracklets), full + night-61642

**From:** Hyak Claude · **Date:** 2026-07-19 · **Re:** fixing_integrator.md §13.6–§13.8

These are the **benchmark** parquets you asked for — our pipeline's own pure-S3M tracklets, propagated
with the new ASSIST n-body integrator, for object-by-object / distribution comparison against Sorcha.
(Distinct from the earlier `arnor_case1_nbody_*` files, which were Sorcha tracklets + our predictions
joined — the *Sorcha* side.)

## Two files

| file | rows | contents |
|---|---:|---|
| `benchmark_night61642.parquet` | 11,612 | night-61642 object set (Sorcha's detected union across case1/2/3), **velocity/position layer** (no P_NEO), carries `s3m_objid` for the identity join to Sorcha |
| `benchmark_comparison_s3m_nbody.parquet` | 670,500 | full v3-scale benchmark (caps NEO 12.9k / MBA 650k / TNO 1.6k / Trojans 6k), **VDP + digest2 scored** (`P_NEO_vdp`, `P_NEO_d2`) |

Same tracklet schema as the old two-body `benchmark_comparison_s3m_v3.parquet`, so it's a drop-in A/B.
The two-body versions are preserved at `benchmark_night61642_twobody/` and `phase2_benchmark_s3m_v3/`.

## How they were built (and why not by re-running the generator)

The benchmark generator (`gen_benchmark_*`) is **doubly broken for n-body**: (1) it loads the hybrid
catalog (`t_0=60065`, 2023) whose elements propagate ~0.47° (max 77°) off the validated pure-S3M cache;
(2) it puts `adam_core_stub` on `sys.path`, shadowing the real `adam_core` that sorcha/assist need, so
`score_orbital_df` silently ran two-body. So both benchmarks are built **directly from the Stage 0
n-body cache** (the validated state that matched Sorcha to 0.03° in acceptance #3), reformatted into the
benchmark tracklet schema. Details in §13.6.

## Key results (already checked here)

- **Night-61642 vs Sorcha** (shared `s3m_objid`): |Δv| median **0.0078** deg/day, raw sky-sep median
  **0.034°** (collapsed from the two-body benchmark's 1.50°), **0% of objects >5°** (was 15%). The
  n-body benchmark now agrees with Sorcha, as expected.
- **Full benchmark VDP (scorable subset): AUC 0.878, F1 0.662** — beats two-body v3 (0.873 / 0.589,
  contamination 22% vs 31%) and **matches Sorcha's n-body VDP** (0.881 full-2yr / 0.893 night-61642).

## IMPORTANT for the VDP comparison — read the metric on the scorable subset

We deliberately left the benchmark in its **full, unfiltered state** (no detection cut). Consequence:
**68% of objects are fainter than mag 25**, outside the maps' 14–25 range, so they are unscorable
(`P_NEO_vdp = 0`, `mag_bin_label` is null). Sorcha never *detects* those, so the raw full-set VDP AUC
(0.074) is meaningless — it's dominated by faint unscorable objects. For any VDP comparison, **filter to
`mag_bin_label.notna()`** (equivalently `14 <= mean_mag < 25`); that's the 0.878/0.662 above. The same
is true of the old v3 (67% faint), so it's apples-to-apples. digest2 (`P_NEO_d2`) is orbit-based and
scores everything (~0.9 AUC). If you'd prefer the on-disk set restricted to detectable objects
(mag < 24.5) to match Sorcha's detection directly, say so — it's a one-line rebuild.

## Columns

`s3m_objid` (true S3M identity — join to Sorcha `ObjID`) · `population` · `mean_ra/mean_dec`,
`ra0/dec0/mjd0_utc`, `ra1/dec1/mjd1_utc` (2 synthetic detections, `ra1 = ra0 + dra·30min`) ·
`mean_dra/mean_ddec`, `vlam/vbeta` (ecliptic rates, deg/day, `vlam = dλ/dt`, no cosβ) · `mean_mag`
(Johnson V, from S3M H) · `lam_deg/beta_deg`, `dlon_from_antisun_deg`, `prob_map/prob_map_file` ·
`H`, `e` · (full only) `P_NEO_vdp`, `P_NEO_d2`, `mag_bin_label`, `digest2_id`.

Epoch: MJD 61642 (2027-08-25), the busiest Sorcha night; grid geometry antisun-relative at that epoch.
Full write-up: fixing_integrator.md §13.6–§13.8.
