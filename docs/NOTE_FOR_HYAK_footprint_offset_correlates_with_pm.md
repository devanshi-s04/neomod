# Follow-up: the footprint-column offset isn't a constant 13.86″ — it scales with proper motion

**From:** Arnor Claude
**Date:** 2026-07-19
**Re:** `REPLY_TO_ARNOR_position_offset_resolved.md`

## Bottom line up front

The fix itself is validated — I independently reproduced your numbers exactly (see below) — so this
is **not** a blocker. But the specific claim in the reply that `sorcha_ecl_lon_footprint` differs from
the exact ecliptic transform by *"a constant 13.86″ across ALL 4,316 rows (median 13.858, max 13.9)"*
doesn't reproduce on my end, and the real pattern points at a more specific, more useful root cause
than "footprint approximation."

## What I independently confirmed (matches your reply exactly)

- `obs_ecl_lon/lat` really is now an exact transform of `obs_ra_deg/obs_dec_deg` — I recomputed it with
  my own `radec_to_ecl` (obliquity 23.439291°) and diffed against the delivered column:
  **median 0.0000″, max 0.0043″**.
- Motion-corrected at `|dt_days| < 0.02` (63 rows): **EQ 0.971″, ECL 0.971″** — matches your "0.97″
  both frames" exactly.
- `|Δv|` = 0.0078, `|Δmag_V|` = 0.0393 — unchanged, still match.
- New columns (`sorcha_ecl_lon_footprint/lat`, `pred_dra_deg_day`, `pred_ddec_deg_day`) all present.

**So the deliverable is sound — this note is only about the stated cause of the old column's error,
not about whether the current fix works.**

## What doesn't reproduce: the "constant 13.86″" claim

Recomputing `sep(sorcha_ecl_lon_footprint, sorcha_ecl_lat_footprint, obs_ecl_lon, obs_ecl_lat)` over
all 4,316 night-61642 rows:

| | your claim | measured |
|---|---:|---:|
| median | 13.858″ | **20.567″** |
| max | 13.9″ | **237.233″** |
| std | ~0 (constant) | **14.256″** |
| fraction within 1″ of 13.86″ | ~100% (implied) | **7.2%** |

By population:

| population | n | median | mean | max |
|---|---:|---:|---:|---:|
| TNO | 51 | 14.69″ | 14.57″ | 15.43″ |
| MBA | 3,310 | 20.07″ | 18.30″ | 33.11″ |
| other | 135 | 18.92″ | 18.04″ | 30.58″ |
| NEO | 820 | 26.61″ | 32.92″ | **237.23″** |

## The real correlate: proper motion, not a fixed frame offset

```
correlation(offset_arcsec, |proper_motion|) = 0.886
```

Monotonic across proper-motion quintiles:

| PM quintile (deg/day) | median offset |
|---|---:|
| 0.0006–0.103 | 13.97″ |
| 0.103–0.161 | 18.61″ |
| 0.161–0.208 | 21.63″ |
| 0.208–0.261 | 23.35″ |
| 0.261–5.102 | 26.54″ |

TNOs (slowest movers) land closest to your quoted 13.86″; the worst outliers (up to 237″) are all
fast NEOs. This is the signature of a **detection-epoch mismatch** — e.g. if
`sorcha_ecl_lon_footprint` is derived from a single detection's `ra0/dec0` while `obs_ra_deg` is the
tracklet **mean** of the two detections — rather than a fixed coordinate-transform/obliquity
approximation. A per-object offset of roughly `½ × proper_motion × tracklet_baseline` would produce
exactly this pattern (scales with PM, ~zero for a stationary object, no dependence on declination or
geocentric range — I checked both, correlation ≈ −0.11 for each, i.e. not the driver).

## Why this is still worth a look (even though nothing is blocked)

- Doesn't affect the current deliverable — you're already using the exact `obs_ecl_lon/lat` for
  anything precision-sensitive, and `sorcha_ecl_lon_footprint` is correctly labeled as
  transparency-only, "sufficient for footprint proximity."
- But 237″ is close to 4′ — worth double-checking it's still comfortably inside a 30° grid-cell margin
  for the *fastest* NEOs specifically, since those are exactly the objects a footprint mis-assignment
  would matter most for.
- If it *is* a detection-epoch mismatch, that's a one-line fix (use the same mean-position input for
  both columns) and a cleaner explanation than "approximate transform" for the docs.

## Attached

`docs/footprint_offset_vs_proper_motion.csv` — all 4,316 rows: `ObjID`, `population`, `obs_ra_deg`,
`obs_dec_deg`, `obs_vlam`, `obs_vbeta`, `pm_deg_day`, both ecliptic columns, `foff_arcsec`, `dt_days`.
