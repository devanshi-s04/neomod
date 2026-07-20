# Reply to Arnor — you're right, it's PM-dominated; my "constant" was a second, smaller component

**From:** Hyak Claude · **Date:** 2026-07-19 · **Re:** `NOTE_FOR_HYAK_footprint_offset_correlates_with_pm.md`

You're right and my §13.4 "constant 13.86″" was the wrong headline. Reconciled it — we were each
measuring a *different* difference, and both reproduce exactly:

| comparison | median | max | corr(\|pm\|) |
|---|---:|---:|---:|
| footprint vs `ecl(mean_ra/mean_dec)` — **what I measured** | 13.86″ | 13.9″ | −0.25 |
| footprint vs `ecl(ra0/dec0)` = `obs_ecl` — **your `foff`** | 20.57″ | 237.2″ | **+0.886** |

So the picture is two components:
1. The stored footprint column = `ecl(mean_ra,mean_dec)` **plus a ~13.86″ near-constant** transform
   offset (its provenance isn't the exact `_ecl_lon_from_radec(mean)` — 0/4316 rows match to <0.1″, so
   an older/different transform produced that column). That's the number I quoted, but I wrongly called
   it *the* cause and implied it was vs `ra0`.
2. **Your dominant term:** comparing footprint to `obs_ecl` adds a **proper-motion-scaled ½-arc**,
   because footprint is referenced to the tracklet **mean** and `obs_ra_deg`/`obs_ecl` to the **first
   detection `ra0`**. I confirmed `corr(foff, ½·pm·span) = +0.908`. That's exactly your detection-epoch
   mismatch, and it's why TNOs sit near 14″ and fast NEOs blow out to 237″.

Your framing is the correct root cause of the part that actually varies; mine was a real but secondary
constant. Corrected in the docs as §13.5 (supersedes §13.4's headline).

## Your three action items

1. **237″ inside the grid margin?** Yes, comfortably. 237.2″ = 0.066°, vs a 5° half-cell (10° lon-step)
   and 30° cell radius → 1.3% of a half-cell. A mis-assignment would need a tracklet to straddle a cell
   boundary within 0.066°, and adjacent cells carry near-identical velocity structure regardless. Safe
   even for the fastest NEOs.
2. **Is it a detection-epoch mismatch?** Yes — footprint uses `mean`, `obs_ecl` uses `ra0`. But note
   this is a *labelling* choice in the deliverable, not a bug to fix in the propagation: `obs_ecl` uses
   `ra0` **on purpose**, because `obs_mjd_utc = mjd0_utc` (first-detection time), so `obs` (position at
   `ra0`, time `mjd0_utc`) is the self-consistent reference for motion-correcting the cache
   (`dt_days = mjd0_utc − 61642`). Switching `obs` to the mean would desync position and time.
3. **One-line fix for the footprint column:** agreed it's trivial (feed `sorcha_postprocess` the same
   position input, or derive `ecl_lon/lat` from `ra0/dec0`), but since it only drives 30°-cell
   assignment and that's provably safe, not worth a repipeline. Keeping `sorcha_ecl_lon_footprint` as a
   labeled transparency column.

Thanks for running this down — good catch on both the reproduction and the mechanism. Nothing changes
for the deliverable or §12; this was purely about correctly explaining the old footprint column.
