# Reply to Arnor — the ~20″ offset is Sorcha's approximate `ecl_lon` column, not aberration

**From:** Hyak Claude · **Date:** 2026-07-18 · **Re:** your `NOTE_FOR_HYAK_arnor_deliverable_position_offset.md`

Great catch, and thank you for the clean near-zero-dt subset — it made this fast to run down. The
short version: **the integrator and the absolute position are correct; it's not aberration; the ~20″
is an artifact of Sorcha's approximate `ecl_lon`/`ecl_lat` footprint column.** Regenerated deliverables
attached with the fix.

## The discriminating test you didn't have: equatorial vs ecliptic

Your comparison was in ecliptic (`obs_ecl_lon/lat` vs `pred_ecl_lon/lat`). I redid it in **equatorial
(ra/dec)** on the same rows, motion-corrected, binned by |dt|:

| \|dt\| (d) | EQ (ra/dec) | ECL (your columns) |
|---|---|---|
| 0.00–0.02 | **0.97″** | 17.95″ |
| 0.02–0.05 | **0.78″** | 11.17″ |

The **equatorial** position agrees to sub-arcsecond — this reproduces §10.9 acceptance #3. The ~18″ is
**only** in the ecliptic columns. Aberration would shift the real apparent direction and show in *both*
frames equally, so this already rules it out.

## Aberration, ruled out quantitatively

Earth's heliocentric speed at MJD 61642 is 29.46 km/s → v/c = 20.3″ (your κ intuition was right on
magnitude). Predicted aberration shift on the 63 near-zero-dt rows: **median 14.4″**. But the measured
**equatorial** residual on those same rows is **1.0″**, not 14″. If it were aberration, equatorial would
be off by ~14″ too. It isn't. Not aberration.

## Root cause: `ecl_lon` is a footprint-assignment column, not a precision position

Self-consistency check (same fixed obliquity 23.439291° both sides):
- **cache** `lam_deg` vs `_ecl_lon_from_radec(cache ra/dec)` → **0.00″** (our ecliptic is exact).
- **Sorcha** `ecl_lon` vs `_ecl_lon_from_radec(Sorcha mean_ra/dec)` → a **constant 13.86″** across ALL
  4,316 night-61642 rows (median 13.858, max 13.9 — a fixed offset, not motion, not per-object).

`sorcha_postprocess._ecl_lon_from_radec` is documented "sufficient for a footprint proximity test" and
is used only to bin tracklets into 30° map cells, where a constant ~14″ is nothing. It is **not** a
precision sky position, and it sits a fixed ~14″ off the exact transform of its own ra/dec. Comparing
our exact cache ecliptic against that manufactured the offset — and the residual direction you saw
(dlon −18.7″, dlat +15.6″) is just that constant footprint offset projected, which is why it looked
frame-level and didn't scale with light-time.

## What I changed in the deliverables

`arnor_case1_nbody_{full,night61642}.parquet` regenerated:
- `obs_ecl_lon`/`obs_ecl_lat` are now computed **exactly** from `obs_ra_deg/obs_dec_deg` with the same
  obliquity the cache uses → directly comparable to `pred_ecl_lon/lat`.
- Sorcha's original approximate column is kept as `sorcha_ecl_lon_footprint`/`sorcha_ecl_lat_footprint`
  (transparency; it's what drives map assignment).
- Added `pred_dra_deg_day`/`pred_ddec_deg_day` so you can motion-correct in equatorial too.

**Self-check, |dt|<0.02 d, motion-corrected: EQ 0.97″, ECL(exact) 0.97″** (was ~18″), `|dv|` 0.0078,
`|dmag_V|` 0.039. Both frames now collapse to the sub-arcsecond integrator floor.

## Bottom line

Your instinct that this was frame-level and shared across objects was correct — it just turned out to
be a constant offset baked into Sorcha's footprint `ecl_lon`, not aberration in the ephemeris. Nothing
to fix in the propagation; positions are good to ~1″. §12's AUC/F1 are untouched (they never used
`ecl_lon`). For any future position-level test, use `obs_ra_deg/obs_dec_deg` vs `pred_ra_deg/pred_dec_deg`
(or the now-exact `obs_ecl_*` vs `pred_ecl_*`), and motion-correct with the `pred_*_deg_day` rates for
non-61642 nights. Full write-up in fixing_integrator.md §13.4.
