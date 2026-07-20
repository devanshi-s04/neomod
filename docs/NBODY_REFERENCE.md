# N-body benchmark/Sorcha reference — parquets, conventions, and the exchange

**Purpose of this doc:** one place that covers every n-body file we now have, what it's for, and
the three issues Hyak and I found/resolved while validating them. Consolidates:
`NOTE_FOR_ARNOR_nbody_benchmark.md`, `NOTE_FOR_HYAK_arnor_deliverable_position_offset.md`,
`NOTE_FOR_HYAK_footprint_offset_correlates_with_pm.md`, `NOTE_FOR_HYAK_scored_nbody_night_benchmark.md`,
`NOTE_FOR_HYAK_sorcha_nbody_Vband_for_scoring.md`, `NOTE_FOR_HYAK_v5_nbody_map_files.md`,
`REPLY_TO_ARNOR_footprint_pm_reconciled.md`, `REPLY_TO_ARNOR_position_offset_resolved.md`,
`REPLY_TO_ARNOR_scored_night_benchmark.md`.

**Getting this doc onto your Mac** — run from your Mac's terminal (not from Arnor):

```bash
scp ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/docs/NBODY_REFERENCE.md ~/Desktop/
```

To also pull the underlying parquets/notebooks referenced below, same pattern — e.g.:

```bash
scp ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/benchmark_night61642_nbody_scored.parquet ~/Desktop/
scp ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/nbody_benchmark_vs_sorcha_roc.ipynb ~/Desktop/
```

(Swap the remote filename and local destination as needed. This is Arnor → Mac; the original
Hyak → Arnor transfers that built these files are documented per-file in §1.)

---

## TL;DR

The integrator went from two-body Kepler to ASSIST n-body (see `fixing_integrator.md` for the
full derivation). Once that landed, Hyak rebuilt every benchmark/Sorcha data product with it, and
across several rounds we (a) got every file we needed on Arnor, (b) chased down two apparent
~20″-scale position bugs that both turned out to be **measurement artifacts in an approximate
column, not integrator errors** — real n-body positions agree to sub-arcsecond. Nothing in §12's
AUC/F1 results changed; the whole exchange was about trusting position-level comparisons.

---

## 1. The parquet files — what we have, on Arnor

All at `/astro/users/ds2004/vdp/` unless noted. Epoch for all n-body products: **MJD 61642
(2027-08-25 UTC)** — the busiest night across Sorcha case1/2/3. Grid is antisun-relative at that
epoch. Magnitude convention: **Johnson V** everywhere scores are computed (Sorcha's raw
`mean_mag` is mixed LSST-band; converted via `mean_mag_V`). Rate convention: `vlam = dλ/dt`,
`vbeta = dβ/dt`, deg/day, **no cosβ factor**.

### Benchmark side (pure S3M, our pipeline)

| file | rows | scored? | purpose |
|---|---:|---|---|
| `benchmark_night61642.parquet` | 11,612 | no (velocity/position only) | night-61642 object set (union of what Sorcha detected across case1/2/3); carries `s3m_objid` for identity join |
| `benchmark_night61642_nbody_scored.parquet` | 11,612 | **yes** — `P_NEO_vdp` + `P_NEO_d2` | same set, scored against `prob_maps_grid_s3m_nbody/`; ~100% scorable (it's the detected set, 98% mag<24.5) |
| `benchmark_comparison_s3m_nbody.parquet` | 670,500 | **yes** — `P_NEO_vdp` + `P_NEO_d2` | full v3-scale benchmark, proportional caps (NEO 12.9k / MBA 650k / TNO 1.6k / Trojan 6k); only **32% (215,645) scorable** — filter `mag_bin_label.notna()`, 68% are fainter than mag 25 |

**Columns (both):** `s3m_objid` (identity key) · `population` · `mean_ra/mean_dec`, `ra0/dec0/mjd0_utc`,
`ra1/dec1/mjd1_utc` (synthetic 2-detection tracklet, `ra1 = ra0 + dra·30min`) · `mean_dra/mean_ddec`,
`vlam/vbeta` · `mean_mag` (already Johnson V — **no band conversion needed on this side**) ·
`lam_deg/beta_deg`, `dlon_from_antisun_deg`, `prob_map/prob_map_file` · `H`, `e` · (scored files)
`P_NEO_vdp`, `P_NEO_d2`, `mag_bin_label`, `digest2_id`.

**Why not just re-run the generator:** `gen_benchmark_*` is doubly broken for n-body — (1) it loads
the hybrid catalog (`t_0=60065`, 2023) whose elements propagate ~0.47° (max 77°) off the validated
pure-S3M cache; (2) it puts `adam_core_stub` on `sys.path`, shadowing the real `adam_core`, so
`score_orbital_df` silently ran two-body. Both benchmarks were instead built **directly from the
Stage 0 n-body cache** (the state validated to 0.03° against Sorcha in `fixing_integrator.md`
acceptance #3), reformatted into the benchmark tracklet schema.

### Sorcha side (case1)

| file | rows | schema | scores |
|---|---:|---|---|
| `sorcha_comparison_case1_nbody_Vband.parquet` | 648,769 | original `mean_*` schema (add-only rescore) | `P_NEO_vdp_Vband` (n-body maps, V) + `P_NEO_d2`; `mean_mag` still raw LSST intact (V−LSST ≈ +0.6, confirmed), `mean_mag_V` the conversion |
| `arnor_case1_nbody_full.parquet` | 648,769 | `obs_*` (Sorcha) / `pred_*` (our n-body prediction) side by side | `P_NEO_nbody`, `P_NEO_twobody`, `P_NEO_d2` |
| `arnor_case1_nbody_night61642.parquet` | 4,316 | same `obs_*`/`pred_*` | same three, plus `sorcha_ecl_lon_footprint`/`lat` (transparency column, see §3) and `pred_dra_deg_day`/`pred_ddec_deg_day` |

**Use `sorcha_comparison_case1_nbody_Vband.parquet`** for anything needing the original schema
(mag before/after band-fix panels, sky maps, ROC by RA/Dec cell) — it's the drop-in replacement for
the old two-body `sorcha_comparison_case1_Vband.parquet`. **Use the `arnor_case1_nbody_*` files**
for position-level work against the Stage-0 cache prediction (`pred_*` columns).

### Maps

`prob_maps_grid_s3m_nbody/` — only 6 of the 667 centers are on Arnor (the ones the v5-normalisation
notebook actually loads): `dlon ∈ {+000, −030, −040, −050, −090, −120}`, all `lat+00`, ~24–33 MB
each. Full grid lives only on Hyak. All scoring already happened server-side — these 6 are for
direct map-image visualization only.

### Evidence CSVs (from the debugging exchange, §3)

- `docs/near_zero_dt_rows_for_hyak.csv` — 63 rows with `|dt_days| < 0.02` (essentially zero elapsed
  time), used to isolate the first position-offset puzzle.
- `docs/footprint_offset_vs_proper_motion.csv` — all 4,316 night-61642 rows, used to show the
  footprint-offset correlates with proper motion (r=0.886).

### Two-body counterparts (kept for A/B)

`outputs/kurlander/benchmark_night61642_scored.parquet`, `docs/benchmark_comparison_s3m_v3.parquet`,
`outputs/kurlander/sorcha_comparison_case1_Vband.parquet` — untouched, still on disk.

---

## 2. Validated results (numbers I independently reproduced, not just claims)

**Kinematic agreement, night 61642 (n-body benchmark vs Sorcha, identity join):**
`|Δv|` median **0.0078** deg/day; sky-separation collapsed **1.50° → 0.034°**; 0% of objects >5°
apart (was 15% under two-body).

**Scoring agreement, night 61642 (same map set, both V):**
ΔP = Sorcha − benchmark, `|ΔP|` median **0.0002**, **97.6%** within 0.1.

**Full benchmark VDP (scorable subset, `mag_bin_label.notna()`):** AUC **0.878**, F1 **0.662** —
beats two-body v3 (0.873 / 0.589, contamination 22% vs 31%) and matches Sorcha's own n-body VDP
(0.881 full-2yr / 0.893 night-61642).

**Same-object, same-epoch NEO score comparison** (the 61 benchmark NEOs with a Sorcha tracklet on
night 61642 — the fair test, since population-level ROC compares different NEO sets: 1,051 capped
benchmark NEOs vs 148,634 Sorcha NEOs):
`|ΔP|` median **0.0000** for both VDP and digest2; 91.8% (VDP) / 86.9% (digest2) agree within 0.1.
→ At matched geometry the two pipelines score the same NEO identically. The population-level F1
gap is a base-rate artifact (0.5% vs 23% NEO), not a measurement inconsistency.

**digest2 still leads VDP** at every elongation bin except 110–141° (near-tied) — same ranking as
two-body. The integrator fix corrected positions/scores; it didn't flip which classifier wins.

---

## 3. The exchange: two ~20″ position puzzles, both resolved as measurement artifacts

### 3.1 First puzzle — "is this aberration?" → No, it's an approximate column

**Found:** comparing `arnor_case1_nbody_night61642.parquet`'s `obs_ecl_lon/lat` (Sorcha) to
`pred_ecl_lon/lat` (our cache) gave ~18–25″ residuals that did **not** shrink to zero even at
`|dt_days| → 0` (same object, same instant) — ruling out motion/curvature as the cause. Direction
was consistent across many unrelated objects (dlon median −18.7″, dlat +15.6″), which pointed at
something frame-level rather than per-object integrator error. Magnitude (~20″) was suspiciously
close to the annual aberration constant (κ≈20.5″) — flagged as the leading hypothesis, not confirmed.

**Hyak's diagnosis:** aberration was ruled out quantitatively — predicted aberration shift on the
same 63 rows was 14.4″, but the measured **equatorial** (RA/Dec) residual was only 1.0″; aberration
would show equally in both frames, and it didn't. The actual cause: **`ecl_lon/lat` in Sorcha's
output is an approximate footprint-assignment column** (`sorcha_postprocess._ecl_lon_from_radec`,
documented as "sufficient for a footprint proximity test"), not a precision transform of its own
RA/Dec — it sat ~14″ off the exact transform of its own coordinates, by construction, because it's
only used for coarse 30° map-cell binning.

**Fix delivered:** `obs_ecl_lon/lat` regenerated as the **exact** transform of `obs_ra_deg/obs_dec_deg`
(same obliquity as the cache). Old approximate column kept as `sorcha_ecl_lon_footprint/lat` for
transparency. Result: **EQ 0.97″, ECL(exact) 0.97″** — both frames collapse to the sub-arcsecond
integrator floor. I independently reproduced this exactly (median diff of my own recomputed
transform vs the delivered column: 0.0000″, max 0.0043″).

### 3.2 Second puzzle — the "constant 13.86″" claim didn't reproduce → it's proper-motion-scaled

**Found:** Hyak's diagnosis said the footprint column differs from the exact transform by "a
constant 13.86″ across all 4,316 rows." Recomputing that exact comparison gave median **20.57″**,
max **237.2″**, only 7.2% of rows within 1″ of 13.86″ — not remotely constant. The real driver:
`correlation(offset, |proper_motion|) = 0.886`, monotonic across PM quintiles (14″ slowest → 26.5″
fastest), with the worst outliers (up to 237″) all fast NEOs.

**Reconciliation:** both measurements were real, just answering different questions. Hyak's 13.86″
was `footprint vs ecl(mean_ra/mean_dec)` — a genuine small constant transform-provenance offset
(corr with PM only −0.25). My 20.57″/237″ was `footprint vs obs_ecl` (= `ecl(ra0/dec0)`) — and
`obs_ecl` uses the **first detection** (`ra0`) while `footprint` uses the tracklet **mean** — a
proper-motion-scaled term (confirmed `corr(offset, ½·pm·span) = +0.908`). This is intentional, not
a bug: `obs_ecl` uses `ra0` on purpose because `obs_mjd_utc = mjd0_utc`, so position and time stay
synced for motion-correcting the cache; using the mean would desync them.

**Resolution:** 237″ = 0.066° is safely inside a 30° grid-cell / 5° half-cell margin (1.3% of a
half-cell) — footprint assignment is robust even for the fastest NEOs. No fix needed; documented as
§13.5 (supersedes §13.4's headline) in `fixing_integrator.md`.

### 3.3 Net effect on everything else

Neither puzzle touched §12's AUC/F1 results (never used `ecl_lon`) or the integrator itself.
**For any future position-level work: use `obs_ra_deg/obs_dec_deg` vs `pred_ra_deg/pred_dec_deg`,
or the now-exact `obs_ecl_*` vs `pred_ecl_*`** — never the `*_footprint` columns for precision work.

---

## 4. Notebooks built on these files

| notebook | uses |
|---|---|
| `nbody_tracklet_velocity_benchmark_night61642_identity.ipynb` | `benchmark_night61642.parquet` + Sorcha case1/2/3 |
| `nbody_scoring_consistency_night61642.ipynb` | `benchmark_night61642_nbody_scored.parquet` + `sorcha_comparison_case1_nbody_Vband.parquet` |
| `nbody_benchmark_v5_normalisation_s3m.ipynb` | `benchmark_comparison_s3m_nbody.parquet` + the 6 map files |
| `nbody_benchmark_vs_sorcha_roc.ipynb` | `benchmark_comparison_s3m_nbody.parquet` + `sorcha_comparison_case1_nbody_Vband.parquet` |

All four run clean end-to-end (verified by execution, not just inspection).
