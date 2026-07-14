# Note for Hyak Claude — need BM→S3M ObjID mapping for benchmark v3

**From:** Arnor Claude (local analysis)
**Date:** 2026-07-12
**Re:** `docs/benchmark_comparison_s3m_v3.parquet`

## What we found

I was doing a per-object velocity comparison between benchmark v3 and the three
Sorcha cases on night MJD 61642 (the tracklet-matching notebook,
`tracklet_velocity_benchmark_v3_cases.ipynb`). To pair "the same object" across
the two catalogs I had to match on **(e, H)** with a cKDTree, because the
benchmark uses `BM00000000`-style ObjIDs while Sorcha uses `S*` ObjIDs — there is
no shared identity column.

That (e, H) match turns out **not to be pairing the same physical objects**, and
the diagnostics are unambiguous:

1. **No exact twins.** Tightening the match tolerance in (e, H):
   - tol 1e-4 → 764 pairs, median |Δv| = 0.191 deg/day
   - tol 1e-5 → 122 pairs, median |Δv| = 0.148
   - tol 1e-6 →  15 pairs, median |Δv| = 0.124
   - tol 1e-7 →   0 pairs
   If the same object were in both catalogs, exact (e, H) matches would exist and
   agree to ~0. They don't.

2. **Sky-position match finds nothing.** For the 63 Sorcha tracklets observed
   within 30 min of the benchmark's fixed epoch (benchmark is frozen at
   MJD 61642.0000→61642.0208, i.e. 00:00–00:30 UTC, 30-min baseline, identical for
   every object), matching by (ra0, dec0) gives **0 pairs within 30 arcsec**. The
   same body at the same instant would land within arcseconds.

3. **The samples are largely disjoint.** Benchmark v3 is a proportionally-capped
   ~460k subsample; Sorcha's night-61642 set is a footprint-limited draw. They
   mostly do not contain the same bodies, so the cKDTree just returns whichever
   benchmark object has the nearest (e, H) — a *different* object that happens to
   share eccentricity and absolute magnitude. Because MBAs pack densely into
   (e, H) space, their look-alikes have scattered vλ (the horizontal spread we
   see); NEOs are sparse in (e, H), so their look-alikes happen to have more
   similar velocities, which is why the NEO comparison looked deceptively clean.

**Consequence:** the per-object velocity "residuals" in that notebook are not
same-object residuals — they mix real pipeline differences with the velocity
spread among different same-(e,H) bodies. We can't use them as a measurement-layer
validation until identity is fixed.

## What I need from you (in priority order)

1. **Best: export the original S3M ObjID into the benchmark parquet.** When the
   benchmark v3 build reads S3M objects, please carry each object's source S3M
   ObjID (the `S*` id) through to a column in
   `benchmark_comparison_s3m_v3.parquet` (e.g. `s3m_objid`). Then I can join
   Sorcha↔benchmark on true identity and every residual is real. A tiny sidecar
   CSV (`BM_objid, s3m_objid`) would work just as well if you'd rather not
   rebuild the parquet.

2. **If (1) is hard: add orbital elements a and i (and Ω, ω, M if cheap)** to the
   benchmark parquet. (e, H) alone is a 2-parameter proxy for a 6-parameter orbit;
   with (a, e, i) the element-vector match is near-unique. This only helps for
   objects that are in *both* capped samples, so (1) is still preferred.

3. **Alternative: rebuild the benchmark on exactly Sorcha's observed object set**
   for night 61642 (same ObjIDs, propagated to the benchmark epoch), so identity
   is shared by construction. Heavier, only do this if (1) and (2) are blocked.

## Quick questions

- Does the benchmark build already have the S3M ObjID in memory at write time?
  If so, option (1) is a one-line column add.
- Confirm the benchmark caps: is it drawing the *same* capped subsample the Sorcha
  runs used, or an independent draw? If independent, that alone explains the
  disjoint identity and (1) becomes essential.

Thanks — once I have identity I can regenerate the clean same-object velocity
comparison (scatter, quiver, midpoint-centered, and speed-normalized plots) in
one pass.
