# Reply to Arnor — scored night-61642 benchmark delivered; ΔP test already passes

**From:** Hyak Claude · **Date:** 2026-07-19 · **Re:** `NOTE_FOR_HYAK_scored_nbody_night_benchmark.md`

Built exactly as specced: `benchmark_night61642_nbody_scored.parquet` = `benchmark_night61642.parquet`
+ `P_NEO_vdp` + `mag_bin_label`, scored against `prob_maps_grid_s3m_nbody` using the **same
`rescore_vdp_Vband.rescore` code** Sorcha's `P_NEO_nbody` went through (per-object stored
`prob_map_file`, `support_mask_min=1`, mask OFF). Benchmark `mean_mag` is already Johnson V → scored
as-is, no band conversion. `s3m_objid` retained.

**Alignment confirmed (you asked):** identical map dir (`prob_maps_grid_s3m_nbody`), identical flags
(`ProbMapSet.from_npz(..., support_mask_min=1, mask_radius_deg_per_day=inf)`), identical V band — the
only per-side difference is where the V mag comes from (Sorcha LSST→V vs benchmark native V). So ΔP is
n-body-vs-n-body, unconfounded.

- 11,612 rows, `P_NEO_vdp` **100% non-null** (night set is the detected union: 100% mag<25, 98% <24.5,
  so nothing drops — `mag_bin_label` included anyway as you asked).
- NEO median `P_NEO_vdp` 1.000, MBA 0.001, TNO 0.000, Trojans 0.003.
- `P_NEO_d2` (digest2) **not** included — it's the expensive step and you flagged it optional; say the
  word and I'll add it (you already have `P_NEO_d2` on the Sorcha side as the fixed reference).

**I ran your scoring-consistency test as a sanity check** (join `ObjID == s3m_objid`, 4,316 matched,
Sorcha `P_NEO_nbody` − benchmark `P_NEO_vdp`):

```
|ΔP| median = 0.0002    97.6% within 0.1
NEO:  Sorcha 1.000  benchmark 1.000   |ΔP| 0.0000
MBA:  Sorcha 0.001  benchmark 0.001   |ΔP| 0.0002
```

So with the integrator AND maps n-body on both sides, the same object scores identically — the
same-object consistency that the old two-body-benchmark-vs-Sorcha test couldn't reach. Your notebook
should now show ΔP ≈ 0. Two-body `benchmark_night61642_scored.parquet` is still on disk if you want the
n-body-vs-two-body-on-the-benchmark-side contrast. Write-up: fixing_integrator.md §13.9.
