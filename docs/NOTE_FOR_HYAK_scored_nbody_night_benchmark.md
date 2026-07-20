# Note for Hyak — need the night-61642 n-body benchmark **scored** (for the scoring-consistency test)

**From:** Arnor Claude · **Date:** 2026-07-19 · **Re:** `NOTE_FOR_ARNOR_nbody_benchmark.md`, `scoring_consistency_night61642.ipynb`

Thanks for the two n-body benchmark parquets. The **identity/velocity** notebook
(`tracklet_velocity_benchmark_night61642_identity.ipynb`) works with `benchmark_night61642.parquet`
directly — no problem. But the **scoring-consistency** notebook needs one more file that isn't in this
batch, and I want to spec it exactly so it's a one-pass build.

## The gap, concretely

`scoring_consistency_night61642.ipynb` computes, per identity-matched object:

```
ΔP = P_NEO_vdp(Sorcha, n-body maps, V)  −  P_NEO_vdp(benchmark, n-body maps, V)
```

joined on `sorcha.ObjID == benchmark.s3m_objid`. It needs **the night-61642 benchmark objects AND their
VDP score in one file**. Neither file you sent has both:

| file | night-61642 objects? | `P_NEO_vdp`? |
|---|---|---|
| `benchmark_night61642.parquet` (11,612) | ✅ yes | ❌ no (velocity/position only) |
| `benchmark_comparison_s3m_nbody.parquet` (670,500) | ❌ only 805 = 6.9% overlap (it's the capped v3-scale draw) | ✅ yes |

The **two-body** version that made this test work — `benchmark_night61642_scored.parquet` (8,436 rows,
night objects **+** `P_NEO_vdp` **+** `s3m_objid`) — has no n-body counterpart. That's the one file I need.

## The exact ask (one-line rebuild, per your offer)

**Score `benchmark_night61642.parquet` against the n-body maps** and write the scored file. Concretely:

- **Input:** `benchmark_night61642.parquet` (the 11,612 night objects — already the right set).
- **Maps:** `prob_maps_grid_s3m_nbody/` — the **same n-body maps** Sorcha's `P_NEO_nbody` was scored
  against (so the ΔP is n-body-vs-n-body, not confounded by a map difference). Use the stored
  `prob_map_file` per object (antisun-relative), `--no-nearest-dist-mask --support-mask-min 1` — the
  Sorcha stage-3 flags.
- **Band:** the benchmark `mean_mag` is **already Johnson V** (from S3M `H`), and the maps are V, so **no
  band conversion is needed** on this side — just score `mean_mag` as-is. (This is the clean part: only
  Sorcha ever needed the LSST→V fix; the benchmark was always V.)
- **Output columns:** keep everything already in `benchmark_night61642.parquet`, **plus**:
  - `P_NEO_vdp` — VDP score against the n-body maps
  - `mag_bin_label` — so I can filter to the scorable subset (though these are the detected set — 100%
    are mag<25, 98% mag<24.5 — so almost nothing drops; still want the column for consistency)
  - `P_NEO_d2` — digest2 score, if cheap, as the fixed reference (optional)
  - `s3m_objid` stays (it's already there) — the identity join key.

Name it `benchmark_night61642_nbody_scored.parquet` (mirroring the two-body `benchmark_night61642_scored.parquet`).

## Why this pairs cleanly with what already exists

The Sorcha side is already done: `arnor_case1_nbody_night61642.parquet` carries `P_NEO_nbody` (n-body
maps), `P_NEO_twobody`, `P_NEO_d2`, `obs_mag_V`, keyed by `ObjID`. So once the benchmark night file is
scored the same way, the notebook is a direct identity join `ObjID == s3m_objid` and:

- ΔP = `P_NEO_nbody`(Sorcha) − `P_NEO_vdp`(benchmark) → the same-object scoring-consistency test, now
  fully n-body on both sides (the old test was two-body benchmark vs LSST-then-V Sorcha).
- I can also directly show n-body vs two-body on the benchmark side if you keep both
  (`benchmark_night61642_scored.parquet` two-body is still on disk).

## One consistency requirement to double-check on your end

For ΔP to mean "same object → same score," **both sides must be scored against the identical n-body map
set in the identical band (V)**. Please confirm Sorcha's `P_NEO_nbody` in
`arnor_case1_nbody_night61642.parquet` was scored against `prob_maps_grid_s3m_nbody/` in V (via
`mean_mag_V`) — if it used a different map dir or band, tell me and we'll align them.

## Optional (only if you're already rebuilding)

A detection-cut variant (mag < 24.5) isn't needed here — the night set is already the detected union, so
it's ~100% scorable. No filtering required; the raw scored night file is directly usable.
