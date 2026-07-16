# Finding: Sorcha VDP scoring uses LSST-band magnitudes against V-band maps

**Date:** 2026-07-16
**Found via:** the benchmark_night61642 identity-matched scoring test (Arnor's magnitude
comparison, `NOTE_FOR_HYAK_benchmark_night61642_scores.md`), traced to source on Hyak.
**Severity:** real, systematic (~0.6 mag), affects every mag-dependent VDP result
(case1/2/3, v5, and the benchmark comparisons). Not a bug in one script — a
band-consistency gap between how the maps are built and how Sorcha tracklets are scored.

---

## TL;DR

The VDP probability maps' magnitude bins are in **Johnson V**. Sorcha's VDP scoring feeds
each tracklet's **observed LSST-filter magnitude** (r/i/z/g `PSFMag`) into those V-band
bins. For red asteroids the LSST bands run ~0.5–0.7 mag brighter than V, so **Sorcha
tracklets are mag-binned ~0.6 mag too bright**. Measured on 3,205 identity-matched
objects: median Δ(Sorcha − benchmark) = **−0.606 mag**, and **~60% of the same physical
objects land in a different VDP mag bin** on Sorcha vs on the (correctly V-binned)
benchmark. The benchmark is the self-consistent one; Sorcha is the mis-binned one.

---

## What `mean_mag` actually is on each side

Neither side is a raw Sorcha output column — both are *derived*, from **different
magnitude systems**:

| Side | Where `mean_mag` comes from | Band |
|---|---|---|
| **VDP maps** (`prob_maps_grid_s3m`) | built by `compute_apparent_magnitude_for_population` = `H + 5log₁₀(r·Δ) − 2.5log₁₀(Φ)`, using the S3M `.s3m` `H` | **Johnson V** |
| **Benchmark** (`gen_benchmark_*.py`) | same `compute_apparent_magnitude_for_population` (`pop_df["mag_app"]`) | **Johnson V** |
| **Sorcha** (`sorcha_postprocess.py:393`) | `(first.PSFMag + last.PSFMag) / 2` — mean of Sorcha's per-visit `PSFMag` | **observed LSST filter** (r/i/z/g, mixed) |

The S3M `.s3m` `H` is the standard Johnson-V absolute magnitude, so the maps and the
benchmark are both in V. Sorcha, by contrast, reports each detection's apparent magnitude
in whatever LSST filter that visit used, and our postprocess averages the two.

### The input photometry system (why they differ, exactly)
From `prep_s3m_sorcha_inputs.py`, the Sorcha physical-parameter file was built as:
- `H_r = H_V − (Johnson V − LSST r)`  → Sorcha's absolute magnitude is **LSST r**,
- per-object LSST colors `u-r, g-r, i-r, z-r, y-r` sampled from `CDS_colors.parquet` (seed 42).

So Sorcha computes `PSFMag` in the observed filter from `H_r` + colors + phase, while the
maps/benchmark compute apparent mag from `H_V`. The two are related, per object, by the
sampled V→LSST color — a *known* quantity, but not applied anywhere in the scoring path.

---

## Evidence (3,205 identity-matched night-61642 objects)

- Δmag (Sorcha − benchmark): **median −0.606, std 0.466**; per population MBA −0.62,
  NEO −0.54, TNO −0.64 — a systematic offset, not scatter.
- Sorcha's night-61642 filters: **i 3,503 · r 2,747 · z 1,555 · g 827** — almost all red
  bands, exactly where V→filter colors are ~0.4–0.9 mag. The offset direction and size
  match asteroid colors.
- Consequence for VDP: **~60% of the identical, same-velocity objects fall in a different
  1-mag bin** on Sorcha vs benchmark (MBA 60% / NEO 55% / TNO 65%). VDP then reads a
  different density slice and returns a different `P_NEO` — purely from the band mismatch,
  before any real disagreement.

The maps are V-band by construction, so the **benchmark's V-mag scoring is correct** and
**Sorcha's LSST-mag scoring is the inconsistent one**.

---

## Why this matters beyond the night-61642 test

Every Sorcha VDP result to date scored tracklets with the observed-LSST-band `mean_mag`
against V-band map bins: the case1/2/3 linking runs, the v5 ROC study, and the
benchmark-vs-Sorcha comparisons. Anything that depends on the magnitude bin is affected by
a ~0.6 mag systematic:
- objects are assigned to a bin ~0.6 mag too bright;
- the fainter mag bins (where NEO/MBA velocity structure matters most for the Rubin
  stream) are the most affected;
- the **absolute** F1/AUC impact is probably modest (a systematic shift, not scatter), but
  any per-mag-bin claim, and any benchmark↔Sorcha score comparison, is contaminated until
  fixed.

Not affected: the kinematics (v_λ, v_β) — those are band-independent and were validated to
median |Δv| = 0.0106 deg/day on the same identity match. This is *only* the magnitude axis.

---

## First-order check of the fix (and why it isn't one line)

The reduced magnitude `PSFMag − H_filter` is filter-independent (pure geometry), so in
principle `V ≈ PSFMag − H_filter + H_V`. Applying that with the benchmark's V-band `H`
removes most of the offset (scatter 0.52 → 0.34 mag) but leaves a **+0.22 mag median
residual**. The residual is second-order: phase-function/`G` treatment differences between
the synthetic benchmark mag and Sorcha's model, trailing-loss in `PSFMag` (moving-source
photometry) vs the benchmark's static synthetic mag, and per-visit photometric noise
(`PSFMagSigma`). So the ~0.6 mag *band* term is cleanly correctable; the last ~0.2 mag is a
genuine synthetic-vs-observed photometry difference that a color correction alone won't
erase. The fix plan (`SORCHA_MAG_BAND_FIX_PLAN.md`) addresses both.
