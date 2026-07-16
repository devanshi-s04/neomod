# Applied: Option A — Sorcha VDP re-scored in Johnson V (2026-07-16)

**What / why:** `SORCHA_MAG_BAND_FINDING.md` showed Sorcha's VDP was scored with observed
LSST-band magnitudes against V-band map bins (~0.6 mag mis-binned). This applies **Option A**
of `SORCHA_MAG_BAND_FIX_PLAN.md`: convert each Sorcha tracklet's magnitude to Johnson V and
re-score VDP, so the mag bin is consistent with the V maps and the V benchmark.

## The conversion (validated)

```
mean_mag_V = mean_mag − 0.5·(color_f0 + color_f1) + (H_V − H_r)
```
- `color_f = (filter − r)` LSST color, from `inputs/s3m_sorcha_phys.csv` (r → 0).
- `H_r` — LSST-r absolute mag, already in the comparison parquet.
- `H_V` — Johnson-V absolute mag, joined by `ObjID` from the S3M `.s3m` census (col 8).

Kinematics `(vlam, vbeta)` are band-independent and reused unchanged; only the magnitude
(hence the mag bin) changes. Re-scored with `score_visible(vlam, vbeta, mean_mag_V)` against
`prob_maps_grid_s3m`, using the SAME flags as Sorcha stage 3
(`--no-nearest-dist-mask`, `--support-mask-min 1`).

**Validation on the 3,205 identity-matched night-61642 objects vs the (V) benchmark:**
- before: `mean_mag − V_bench` median **−0.606**, std 0.466
- after:  `mean_mag_V − V_bench` median **+0.007**, std 0.328
The systematic band offset is gone. The residual std 0.328 is the expected
synthetic-vs-observed photometry scatter (phase `G`, trailing loss, per-visit noise), not a
band term. `H_V` from `.s3m` matches the benchmark's `H` exactly (max|Δ| = 0).

## What was produced

Script: `neomod/pipeline/rescore_vdp_Vband.py`; job: `neomod/pipeline/slurm/rescore_Vband.sbatch`
(3-task array, one per case). For each case a NEW file
`outputs/s3m_linking/<case>/sorcha_comparison_<case>_Vband.parquet` with the original columns
**plus**:
- `mean_mag_V` — Johnson-V tracklet magnitude,
- `mag_bin_label_Vband` — the V mag bin,
- `P_NEO_vdp_Vband` — VDP score using the V mag.

The original `mean_mag` / `P_NEO_vdp` are **untouched** (add-only; nothing silently changes,
old results stay reproducible).

## How to use it

- **Benchmark side is already V** (`benchmark_night61642_scored.parquet` `P_NEO_vdp`), so the
  clean same-object comparison is: benchmark `P_NEO_vdp`  vs  Sorcha `P_NEO_vdp_Vband`
  (join `sorcha.ObjID == benchmark.s3m_objid`). Both now V-binned — the ~60% mag-bin-crossing
  from the band mismatch should collapse toward the noise floor.
- To measure how much the mis-binning cost: compare Sorcha `P_NEO_vdp` (old, LSST-binned) vs
  `P_NEO_vdp_Vband` (new) F1/AUC on each case. That delta is a paper-worthy number.

## Scope / caveats

- This is a **scoring-time** fix (Option A). It does not rebuild the maps. If we later decide a
  Rubin-era classifier should live in an LSST band (Option B, maps rebuilt in r), that supersedes
  this — but Option A and B give the same *relative* benchmark-test answer.
- `P_NEO_d2` is unaffected (digest2 uses geometry + its own magnitude handling; a separate band
  question if we pursue it).
- Any object not in the `.s3m` census gets `mean_mag_V = NaN` (should be ~none for these files).
