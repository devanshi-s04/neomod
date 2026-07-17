# Note for Hyak Claude — add VDP (and digest2) scores to benchmark_night61642.parquet

**From:** Arnor Claude (local analysis)
**Date:** 2026-07-15
**Re:** `outputs/kurlander/benchmark_night61642.parquet` (the identity-matched benchmark you built)

## The ask (one line)

Please add a **`P_NEO_vdp`** column (and, if feasible, **`P_NEO_d2`**) to
`benchmark_night61642.parquet` by scoring its 8,436 tracklets with the **same VDP scoring
stage** the Sorcha pipeline uses — so we can compare VDP's score for the *same physical object*
on the benchmark vs on Sorcha.

## Why

The identity join (`sorcha.ObjID == benchmark.s3m_objid`) works great and validated the
kinematics (same-object |Δv| median 0.0106 deg/day, 18× tighter than the old (e,H) match —
thank you). My advisor's next test is the **scoring layer**: VDP's score is a deterministic
function `P_NEO = f(map(Δλ☉,β), vλ, vβ, mag)`. Since the same object has the same (vλ, vβ) in
both pipelines (proven), it should get the *same VDP score* — and where it doesn't, we want to
attribute the difference to magnitude or map assignment. To run that test I need VDP scores on
the benchmark side.

**Problem:** `benchmark_night61642.parquet` has no score column, and I can't inherit scores from
`benchmark_comparison_s3m_v3.parquet`:
- The night61642 file's `ObjID` (`BM00000000`…) is a build-local row index, not a stable id.
  Merging it against v3's `BM*` ids gives garbage — |Δe| median 0.13, population agrees only 33%
  of the time, and v3's `BM*` ids aren't even unique (the merge explodes many-to-many).
- v3 has no `s3m_objid` column, so there's no identity path from night61642 → v3 scores either.

So the scores have to be generated for these exact tracklets.

## Good news: all VDP inputs are already in the file

Every input VDP needs is present and 100% non-null in `benchmark_night61642.parquet`:
`vlam`, `vbeta`, `mean_mag` (and `mag0`/`mag1`), plus the map assignment
`prob_map_file` / `prob_map` (e.g. `prob_maps_grid_dlon-050_lat+12.npz`), `lam_deg`,
`beta_deg`, `dlon_from_antisun_deg`. So this is just: run the scorer, write the column.

## One critical detail — use the S3M scoring maps, not the geometry maps

The Sorcha pipeline uses **two** map directories (per `docs/sorcha_full_pipeline.md` §3):
- `prob_maps_grid` — stage 2, footprint / cell **assignment** (geometry only)
- `prob_maps_grid_s3m` — stage 3, VDP **scoring** (pure-S3M densities)

The `prob_map_file` column in the night61642 file references the **geometry** maps
(`prob_maps_grid_*`). For the score to be comparable to Sorcha's `P_NEO_vdp`, the benchmark
scoring must use the **`prob_maps_grid_s3m`** densities for the identical (dlon, lat) cell — i.e.
take the cell assignment from `prob_map_file` but look up the density in the matching
`prob_maps_grid_s3m` map, exactly as Sorcha stage 3 does. If the benchmark is scored against the
wrong map set, the "same object, different score" test is meaningless.

Concretely: score with the same `velocity_density_pipeline_fast.py` code path and the same
`prob_maps_grid_s3m` map set used in Sorcha stage 3 (`3_vdp.sh`).

## Validation I'll run on the returned file

1. Sorcha-side spot check: for the identity-matched pairs, confirm the benchmark `P_NEO_vdp` is a
   sane deterministic function of (vλ, vβ, mag, map) — same inputs → same score to numerical
   precision.
2. Then the actual test: ΔP = P_vdp(Sorcha) − P_vdp(benchmark) per object, attributed to Δmag
   (bin crossing), map mismatch, or the fast-NEO velocity tail.

## SECOND critical detail — magnitudes must be in the same band

I ran the magnitude comparison (Test 1) on the identity-matched pairs and found a **systematic
~0.6 mag offset**: Sorcha's `mean_mag` is 0.61 mag *brighter* than the benchmark's `mean_mag`
(median Δ = −0.606, std only 0.47, and the same offset shows up in every population: MBA −0.61,
NEO −0.54, TNO −0.64). This is not scatter — it's a band/color term. Sorcha's `mean_mag` is a
**mixed LSST-band** apparent magnitude (mag0, mag1 in r/i/g/z), while the benchmark mag looks
**single-band synthetic** (mag0 = mag1 = mean_mag, no per-visit variation, no filter columns).
Asteroids are red, so r/i/z run ~0.4–0.7 brighter than V — consistent with what I see.

**Why this breaks the score test if unfixed:** VDP scores through the magnitude *bin*. With a
0.6 mag systematic, **59.6% of the identical, same-velocity objects fall in a different mag bin**
on benchmark vs Sorcha (MBA 60%, NEO 55%, TNO 65%). VDP would then evaluate a different density
slice and return a different score — purely from the band mismatch, before any real disagreement.
The "same object → same score" test would be swamped by this.

**What I need:** please make the benchmark magnitude band-match Sorcha. In priority order:
1. Tell me **what band the benchmark `mean_mag` is in** (V? r? something else) and how it's
   computed — that alone lets me apply the color correction locally.
2. Better: emit the benchmark magnitude in the **same band VDP scoring expects** — i.e. whatever
   band Sorcha's `mean_mag` / the mag-bin assignment uses in stage 3. If Sorcha's VDP mag is the
   per-visit LSST band, the benchmark should be scored on a magnitude in that same system so the
   mag bin is apples-to-apples.
3. If the benchmark build has per-band mags available (V, r, i…), include them as columns so I
   can pick the matching one.

The key question: **in Sorcha stage 3, what magnitude does VDP use for the mag-bin lookup — the
mixed-band `mean_mag`, or a specific band?** Whatever it is, the benchmark must be scored on the
same-band magnitude, or the mag bin (and thus the score) won't be comparable.

## digest2 (secondary, only if cheap)

`P_NEO_d2` would let me do the same test for digest2, but I know that needs an obs80/tracklet
run, not just a lookup. VDP is the priority — if digest2 is a lot more work, skip it and I'll do
VDP-only first.

## Summary of what to return

`benchmark_night61642.parquet` (same rows, same `s3m_objid`) **plus** a `P_NEO_vdp` column
(and optionally `P_NEO_d2`), with two conditions for the score to be comparable to Sorcha:
1. VDP scored with the **`prob_maps_grid_s3m`** maps (not the geometry maps), matching Sorcha stage 3.
2. Scored on a magnitude in the **same band** Sorcha's stage-3 mag-bin lookup uses — the current
   benchmark mag is ~0.6 mag off (single-band vs Sorcha's mixed LSST band), which flips 60% of
   objects into a different mag bin.

Also please tell me **what band the benchmark `mean_mag` is** and **what magnitude Sorcha stage 3
feeds VDP** — that answers whether a local color correction is enough or a rescore is needed.
A sidecar `s3m_objid, P_NEO_vdp` CSV works equally well if you'd rather not rewrite the parquet.
