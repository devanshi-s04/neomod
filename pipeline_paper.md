# NEO Probability Pipeline — Complete Methods Reference (for paper writing)

**Purpose.** A single, exhaustive reference for the methods section: every input, parameter,
script, threshold, coordinate convention, and design decision in the Sorcha → VDP/digest2
comparison pipeline. Where a number comes from the earlier 2-year run (baseline v3.3 window)
vs the current **v5.0** run, it is labelled. Results for v5.0 are pending (the 667-map grid is
generating as of 2026-06-16); prior-run results are given for context and are expected to be
reproduced/updated by the v5.0 run.

**Working dir (Hyak/Klone):** `/mmfs1/gscratch/dirac/ds2004/sorcha/`
(the old `gscratch/astro/...` path is a symlink → `dirac`).
**Repo:** `neomod/` (GitHub `devanshi-s04/neomod`). Paper `.tex` is in Overleaf, *not* this repo.
**Python:** always the project env — `conda_prep/bin/python`, `conda_prep/bin/sorcha-run`
(the shell shows `(base)`, which is NOT the active env — ignore it).

---

## 0. Goal and one-paragraph summary

Reproduce the Tom Wagg et al. (arXiv:2408.12517) LSST/Rubin NEOCP simulation workflow with
Sorcha, then score the resulting nightly tracklets with two classifiers — the Velocity
Density Pipeline (VDP, our method) and digest2 (the MPC's NEOCP tool) — for a head-to-head
ROC comparison on realistic Rubin-cadence simulated data. Sorcha simulates Rubin detections
of a 14.4-million-object hybrid S3M+MPCORB Solar System population over 2 years of the Rubin
baseline cadence; post-processing builds Wagg-style nightly tracklets and scores them.

**Our contribution vs Wagg:** Wagg predicts NEOCP *traffic* (uses difi for self-linkage, no
classifier comparison). We instead compare VDP vs digest2 as NEO classifiers via full ROC.
We do **not** implement difi. Observatory code X05 (Rubin) in both.

---

## 1. Input population: hybrid S3M + MPCORB catalog

Source: `hybrid.h5`, key `/df`, shape **14,444,912 × 10**: `id, x, y, z, vx, vy, vz, t_0, H, g`.
- Heliocentric **ecliptic J2000** Cartesian state vectors (AU, AU/day) at epoch `t_0` (MJD TDB,
  = 60065 for the catalog).
- `H` = absolute magnitude (V-band); `g` = phase slope (≈0.15).
- Composition: **S3M** synthetic (TNOs, Centaurs, MBAs, Trojans) + **MPCORB** real known
  NEOs/MBAs. S-prefix = S3M synthetic (13,112,612; 90.8%); non-S = real MPCORB (1,332,300;
  9.2%). The real MPCORB NEOs (~33k) carry wild e/i that S3M does not synthesize — these fill
  VDP calibration gaps.

### Sorcha input CSVs (from `hybrid.h5`)
- `inputs/hybrid_sorcha_orbits.csv` (~2.0 GB, 14,444,912 rows):
  `ObjID, FORMAT=CART, x, y, z, xdot, ydot, zdot, epochMJD_TDB`
- `inputs/hybrid_sorcha_phys.csv` (~1.8 GB): `ObjID, H_r, GS, u-r, g-r, i-r, z-r, y-r`
  - Colors sampled from `CDS_colors.parquet` (LSST color distribution, 534 rows).
  - Magnitude conversion: `H_r = H_V − (Johnson V − LSST r)`.
- Note: Sorcha ObjIDs are `S`-prefixed base-N strings (e.g. `S00000dra`) — `S` is Sorcha's
  encoding for ALL objects, NOT a population label. Population is classified post-hoc from the
  state vectors (Section 6).

### Keplerian element catalog for VDP training
`hybrid_catalog_prep.py` → `hybrid_elements.parquet` (0.93 GB, 14,444,912 rows):
`OID, q, e, i, node, argperi, t_p, H, t_0, a, population`.
- Cartesian → Keplerian (heliocentric two-body, μ = 4π²/365.25² AU³/day²): specific energy
  `eps=v²/2−μ/r`, `a=−μ/2eps`, `h=r×v`, `i=arccos(h_z/|h|)`, `e=|e_vec|` (Laplace-Runge-Lenz),
  `q=a(1−e)`; node, argperi, t_p from standard relations.
- Built with chunked reads (500k rows/chunk, 29 chunks) to avoid login-node OOM
  (`pd.read_hdf(chunksize=)` fails on fixed-format HDF5; flat `block0_values[:]` also OOMs).
- Population labels in this catalog: NEO `q<1.3`; MBA `1.7≤a<4.1, q≥1.3`;
  Trojan `4.7<a<5.9, e<0.3`; TNO `a>30`; other = rest. Counts: MBA 13,883,703; NEO 295,040
  (≈261k S3M + ≈33k real MPCORB); Trojan 185,417; TNO 63,381; other 17,371.
- Frame check (`sanity_check_frame.py`): propagated S3M Keplerian elements to t=60065 vs
  hybrid.h5 vector → residual |dr|=1.1e-3 AU, |dv|=9 m/s (physical N-body drift over 16 yr,
  not a frame error). Confirms heliocentric ecliptic J2000.

---

## 2. Rubin pointing database

`baseline_v5.0.0_2yrs.db` — SQLite3, **148 MB, 414,488 observations**, MJD **60980.0–61709.4**
(**2025-Nov-04 to 2027-Nov-08**, a 730-day window).
- Trimmed from the full 10-year `baseline_v5.0.0_10yrs.db` (728 MB):
  ```sql
  ATTACH DATABASE 'baseline_v5.0.0_2yrs.db' AS new;
  CREATE TABLE new.observations AS
    SELECT * FROM observations WHERE observationStartMJD < 61710.0;
  CREATE TABLE new.info   AS SELECT * FROM info;
  CREATE TABLE new.events AS SELECT * FROM events;
  ```
- **v3.3 vs v5.0 windows:** v3.3 = 2025-May–2027-Apr (`outputs/production_2yr/`); v5.0 =
  2025-Nov–2027-Nov (`outputs/production_2yr_v5/`). The 6-month offset changes orbital phases,
  which is why v5.0 surfaced 2 *new* pathological objects (Section 5).

---

## 3. Sorcha configuration

File: `neomod/pipeline/config/Rubin_full_footprint_detections.ini`. Exact settings:
```
[INPUT]
  ephemerides_type = ar          # ASSIST+REBOUND N-body ephemeris
  eph_format = csv ; aux_format = csv ; size_serial_chunk = 20000
  pointing_sql_query = SELECT observationId, observationStartMJD as observationStartMJD_TAI,
    visitTime, visitExposureTime, filter, seeingFwhmGeom as seeingFwhmGeom_arcsec,
    seeingFwhmEff as seeingFwhmEff_arcsec, fiveSigmaDepth as fieldFiveSigmaDepth_mag,
    fieldRA as fieldRA_deg, fieldDec as fieldDec_deg, rotSkyPos as fieldRotSkyPos_deg
    FROM observations order by observationId
[SIMULATION]
  ar_ang_fov = 2.06 ; ar_fov_buffer = 0.2 ; ar_picket = 1
  ar_obs_code = X05              # Vera C. Rubin Observatory
  ar_healpix_order = 6           # FOV pre-filtering discretisation
[FILTERS] observing_filters = r,g,i,z,u,y
[SATURATION] bright_limit = 16.0 # saturated sources removed
[PHASECURVES] phase_function = HG
[FOV] camera_model = footprint   # full LSSTCam CCD footprint incl. chip gaps
      footprint_edge_threshold = 2.   # arcsec; exclude detections within 2" of detector edge
[FADINGFUNCTION] fading_function_width = 0.1 ; fading_function_peak_efficiency = 1.
[OUTPUT] output_format = hdf5 ; output_columns = all   # 56 columns
[LIGHTCURVE] lc_model = none
[ACTIVITY] comet_activity = none
```
- **No `[LINKINGFILTER]` section** — all individual detections kept (no linking/tracklet
  filtering at the Sorcha stage; tracklets are built in Phase 1).
- SPICE/ephemeris kernels pre-cached in `sorcha_cache_2025-07-06/` (passed via `--ar`).

### Sorcha output schema (56 columns; key ones)
`ObjID, fieldMJD_TAI, RA_deg, Dec_deg, RATrue_deg, DecTrue_deg, RARateCosDec_deg_day,
DecRate_deg_day, PSFMag, PSFMagTrue, H_r, H_filter, optFilter, detectorID, Range_LTC_km,
Obj_Sun_LTC_km, seeingFwhmGeom_arcsec, fieldFiveSigmaDepth_mag, SNR, x,y,z,xdot,ydot,zdot,
epochMJD_TDB` (the state-vector columns are the *input* vectors preserved per detection).
- HDF5 key: `/sorcha_results` (read with `pd.read_hdf(f, key="sorcha_results")`).
- **Column choice rules** (critical):
  - Use **observed** quantities for scoring: `RA_deg`, `Dec_deg`, `PSFMag` (astrometric/photometric
    noise applied). NOT `RATrue_deg`/`DecTrue_deg`/`PSFMagTrue`.
  - Rates `RARateCosDec_deg_day`, `DecRate_deg_day` are noiseless instantaneous ephemeris rates.
  - **Rate conversion for VDP:** `dra_deg_day = RARateCosDec_deg_day / cos(Dec_deg)`;
    `ddec_deg_day = DecRate_deg_day`.

---

## 4. HPC infrastructure (Hyak/Klone)

- Partition `ckpt` (checkpoint/preemptible), `--account=astro`. `--partition=ckpt` needs no
  explicit account for basic use.
- **NEVER use `--requeue` on ckpt:** preempted requeued jobs fail to re-fetch
  `--export=ALL,...` env → stuck in `user env retrieval failed requeued held`. Instead use
  **skip-existing** logic + manual resubmit for stragglers.
- **QOS limit:** `QOSMaxSubmitJobPerUserLimit = 2000` tasks (pending+running). Group ckpt cap
  ~377 jobs (`AssocGrpJobsLimit`).
- VDP import needs `adam_core`; satisfied by a minimal stub at `neomod/adam_core_stub/`
  (`score_observation` and `ProbMapSet.from_npz` don't call adam_core at runtime — only orbit
  propagation does). Map-generation scripts `os.chdir(neomod/)` so `s3m_loader` finds `S3Mdata/`.
- Production benchmarks (16 cores): 16,000 objects → 21–24 min, MaxRSS 53–56 GB.

---

## 5. Step 1 — Sorcha production run

**Worker:** `production_run/multi_sorcha_production.py`.
**Slurm:** `neomod/pipeline/slurm/multi_sorcha_production_v5.sh`
(`--cpus-per-task=16 --mem=128G --time=04:00:00 --array=0-903%16`,
`--chunksize 16000 --norbits 1000 --cores 16`).

### Chunking
14,444,912 rows → **904 instances** (array 0–903). Each instance reads a 16,000-row slice and
splits it into **16 parts of ≤1,000 objects**, run simultaneously via
`multiprocessing.Pool(16)`; each part calls `sorcha-run` → one `inst{NNNNN}_part{NNN}.h5`.
Per-part orbit/physical CSV slices land in `work/production_2yr_v5/instance_{NNNNN}/`.

### Output
`outputs/production_2yr_v5/inst{NNNNN}_part{NNN}.h5` — **14,445 files, ~232 GB**.
(Instances 0–901 → 16 parts each; 902 → 13 parts; 903 → 0 rows.)
**Skip-existing:** `h5_is_good()` checks the file exists and has `/sorcha_results` before each
`sorcha-run`; safe to resubmit any array index without `--requeue`.

### Step 1b — Pathological objects (integrator hangs)
A handful of objects cause ASSIST+REBOUND to take astronomically small adaptive steps (close
planetary encounters) and hang for tens of minutes to indefinitely. Only `inst00820_part003`
was affected. **5 objects excluded total:**
- `A804 RA`, `A854 OA`, `A868 TA` — known from the v3.3 run.
- `A854 RA`, `A868 WA` — **new in v5.0** (different orbital phase from the 6-month window shift).
**Found by 3-round bisection:** 16×~63-obj batches (15-min timeout) → 16×8-obj → individual
objects; timed-out tasks contain the bad object.
**Recovery:** `neomod/pipeline/slurm/split820v5_part003_final.sh` runs `sorcha-run` directly on
`work/production_2yr_v5/instance_00820/orbits_003_skip5.csv` (995 objects) →
`inst00820_part003.h5`.

---

## 6. Step 2 — Phase 1: tracklet construction

**Worker:** `neomod/pipeline/sorcha_postprocess.py`.
**Slurm (v5.0, grid maps):** `neomod/pipeline/slurm/sorcha_postprocess_v5_grid.sh`
(`--cpus-per-task=1 --mem=8G --time=01:00:00 --array=0-999%64`),
`--indir outputs/production_2yr_v5 --outdir outputs/tracklets_v5_grid
--prob_maps_dir prob_maps_grid --overwrite`.
**Key insight:** all detections of an ObjID are within one h5 file → tracklets are
self-contained per file → process each file independently in a Slurm array.

### Batching under the QOS cap
14,445 files > 2,000-task cap → 15 batches via `INDEX_OFFSET` (file_index = array_id +
INDEX_OFFSET), 2 batches at a time:
```
sbatch --export=ALL,INDEX_OFFSET=0     --array=0-999%64 sorcha_postprocess_v5_grid.sh
... INDEX_OFFSET=1000, 2000, ... 13000 (--array=0-999)
    INDEX_OFFSET=14000 --array=0-444     (final partial batch)
```

### Tracklet logic
For each h5: read `/sorcha_results`; drop non-finite RA-rate/Dec-rate/mag; group by
`(ObjID, night)` where `night = floor(fieldMJD_TAI)`; take first & last detection of the night;
keep groups with **n_det ≥ 2** AND **3 min ≤ dt ≤ 90 min** (Wagg-style pairing window).
Per-tracklet quantities: mean RA/Dec/rates/mag, `mjd0/mjd1` (TAI and UTC), ecliptic lon/lat,
orbital elements (a,e,q), population, and prob-map assignment. Atomic write (tmp→rename);
`--skip-existing` for safe resubmit, `--overwrite` for the re-run.

### Output columns (44)
```
tracklet_id, source_file, source_file_index, ObjID, night, n_det, night_span_min, dt_min,
mean_ra, mean_dec, mean_dra, mean_ddec, mean_mag,          <- VDP inputs (rates from Sorcha)
ra0, dec0, mjd0_tai, mjd0_utc, mag0, filter0, snr0,        <- digest2 obs 1
ra1, dec1, mjd1_tai, mjd1_utc, mag1, filter1, snr1,        <- digest2 obs 2
H_r, H_filter, range_km, obj_sun_km, ecl_lon, ecl_lat,
prob_map, prob_map_file, prob_map_obstime_str,
prob_map_center_lon_deg, prob_map_center_lat_deg, prob_map_dist_deg,
population, a_au, e, q_au, n_det_per_night
```
- **`n_det_per_night` (col 44, v5.0 addition)** = total raw Sorcha detections for that
  (ObjID, night) before first/last pairing. Lets ≥3/≥4-detection re-analysis (Wagg uses ≥3)
  happen post-hoc without re-running Phase 1.
- v5.0 tracklet count (monthly-map run, in `outputs/tracklets_v5/`, before the grid re-run):
  **29,844,550** tracklets, 14,445 parquets, 7.85 GB.

### Population classification (`classify_population`, from state vectors)
`MU = 0.01720209895²`. NEO `q<1.3`; MBA `2.0<a<3.3, e<0.3`; Trojan `4.8<a<5.6`; TNO `a>30`;
other = rest. (Slightly different bins than the *training* catalog in §1 — this one labels
tracklets; treated as non-NEO unless NEO for the ROC.) "other" = real MPCORB designations that
don't fit a clean bin.

### Tracklet → map assignment (`assign_probability_maps`)
Two modes, auto-selected:
- **Monthly antisun maps:** footprint = within `max_sep_deg` (30°) of the **current antisun**
  ecliptic longitude at observation time (NOT the map's historical sky position — this fix
  prevented off-axis MBAs near their stationary point from being scored as NEOs). Among
  in-footprint maps, pick nearest epoch.
- **Grid maps (v5.0)** — auto-detected by the `delta_lon_from_antisun_deg` key: place each
  tracklet in antisun-relative coords (Δlon = ecl_lon − antisun_lon wrapped to ±180°; ecl_lat)
  and assign to the **nearest grid center** in 2-D (Δlon, lat). Tracklets within 40° of the Sun
  are unassigned.
- Helpers: `_antisun_ecl_lon(mjd)` (low-precision Sun formula, ~1°), `_ecl_lon_from_radec`,
  `_ecl_lat_from_radec`, `_angular_sep_ecl_lon` (all fixed-obliquity, vectorised).

---

## 7. VDP — Velocity Density Pipeline (probability maps)

Core idea: at a given sky patch and epoch, NEOs and MBAs separate cleanly in apparent-velocity
space (vlam, vbeta) — especially near the antisun, where NEOs show fast retrograde motion
(vlam ≈ −0.4 to −0.6 deg/day) while MBAs are near-stationary. A VDP "map" is, per magnitude
bin, the probability P(NEO | vlam, vbeta) on a velocity grid.

**Code:** `neomod/src/velocity_density_pipeline_gmm.py` (the current pipeline; older variants
`_fast.py` (S3M kNN) and `_hybrid.py` exist and are untouched). Loader: `hybrid_loader.py`
(drop-in for `s3m_loader.py`; keeps function names `define_s3m`, `s3m_array`).

### Velocity grid & per-map params
- Grid: `DEFAULT_GRID_LIM = (−2.0, +2.0)` deg/day, `DEFAULT_GRID_STEP = 0.01` → **401×401**
  in (vlam, vbeta).
- `DEFAULT_MAX_SEP_DEG = 30` (training sky-cut radius around the map center).
- Magnitude bins (8): `14_16, 16_18, 18_20, mag20, mag21, mag22, mag23, mag24+`.

### Density estimation
- **NEO: GMM cloner.** Feature vector per object: `log(a), q, sin/cos(i), sin/cos(node),
  sin/cos(argperi), sin/cos(M_obs)` (mean anomaly at obstime). `StandardScaler` → Gaussian
  Mixture, **n_components = 80**, `covariance_type='full'`, `reg_covar=1e-6`,
  `random_state=42`. Trained on **visible** source NEOs per mag bin (same sky cut as K|M).
  Sample `clone_factor × N_vis` clones from the full NEO orbit distribution, re-apply the sky
  cut; **acceptance ≈ 0.556**. `H` sampled empirically from the source NEO H-distribution (NOT
  a GMM feature — preserves the faint-end rise). `e` derived from (a, q), never sampled
  independently. K|M fallback if GMM fails or <2 visible clones (happens in very sparse bins).
- **MBA / TNO / Trojan: conditional K|M cloner** (no GMM).
- **Clone factors** (`DEFAULT_POPULATION_SETTINGS`): NEO **80**, MBA **5** *(see §7.1)*,
  TNO 10, Trojans 5. (Code default for MBA is 1.)
- Density evaluator parallelised with joblib (`--n-jobs`).

### Smoothing & masking
- Smoothing applied to **NEO only** (`DEFAULT_SMOOTH_POPULATION_NAMES = ('NEO',)`),
  support-aware Gaussian: passes at support thresholds `[3, 10]`, sigma_pixels `[1.5, 4.0]`,
  truncate `[5, 5]`; `smooth_support_scale_by_clone_factor = True`.
- **nearest_dist mask:** cells whose nearest cloned point is > 0.2 deg/day away → P=0.
  Stored always; toggled at scoring via `--no-nearest-dist-mask`. For the GMM pipeline the mask
  is **OFF** in Phase 2 (`sorcha_phase2_vdp_gmm.sh` passes `--no-nearest-dist-mask`).
- `.npz` stores raw ingredients (`density_raw__POP__BIN`, `support_count__POP__BIN`,
  `nearest_dist__POP__BIN`, `magcut_count__POP__BIN`, `x_grid`, `y_grid`, center/obstime/
  max_sep/labels, smoothing params; overlays only if `--save-overlays`). The final P(NEO) is
  computed at load time by `ProbMapSet.from_npz` as density_NEO / Σ density_pop.

### 7.1 GMM normalisation fix (committed 2026-06-16)
Diagnosed 2026-06-07. The downweighting used `density_clone_map / f` (intended clone_factor),
but GMM samples from the full orbit space then sky-cuts, so only `acceptance ≈ 0.32–0.56`
survive — underscaling NEO density by ~3×. **Fix:** for the NEO GMM path,
`effective_factor = n_visible_clones_gmm / n_source` (and pass it to the support map &
smoothing), else `effective_factor = f`. Lives in `build_cloned_maps_for_center_magbin`.
This fix was uncommitted until 2026-06-16 (commit `b4cf62b`).

### 7.2 MBA clone_factor = 5 (the F1≈0.837 config)
The maps behind the documented F1≈0.837 result used **MBA clone_factor = 5** (confirmed:
`support_count__MBA` is exactly 5× the cf=1 maps; MBA is the only population that differs;
NEO/TNO/Trojan match). The code default is 1. The v5.0 grid uses cf=5 to match, applied via
`--mba-clone-factor 5` (default) in `sorcha_gen_maps_grid.py`, which builds a custom
`population_settings` so the global `DEFAULT_POPULATION_SETTINGS` (and the monthly/hybrid
pipelines) stay untouched. This was a *deferred* small win, adopted on instruction — not a
previously forbidden change.

---

## 8. Step 3 — the 667-map antisun-relative sky grid (v5.0)

See `docs/SORCHA_V5_PIPELINE.md` §3 for the full as-built block. Essentials:
- **Coords:** (Δlon-from-antisun, ecliptic-lat) — time-independent, reusable every year.
- **Longitude:** 10° steps, |Δlon| ≤ 140° (≥40° from Sun) → **29 usable** (−140°…+140°).
- **Latitude (non-uniform):** symmetric `0,1,2,3,4,5,8,12,18,25,35,50` → **23 values**
  (0, ±1…±5 at 1° steps near the ecliptic, coarsening to ±50°).
- **Total: 29 × 23 = 667 maps.** Reference epoch 2026-01-01T00:00:00. MBA cf=5.
- **Script:** `neomod/pipeline/sorcha_gen_maps_grid.py` (configurable `--lon-step`,
  `--lat-base`/`--lat-points`, `--sun-exclusion`, `--ref-obstime`, `--mba-clone-factor`,
  `--task-id`, `--list-only`/`--map-grid-file`, `--n-jobs`, `--save-overlays`). Each `.npz`
  augmented with `delta_lon_from_antisun_deg`, `grid_lat_deg`, `ref_obstime_str`.
- **Slurm:** `sorcha_gen_maps_grid_slurm.sh` (`--cpus-per-task=16 --mem=32G --array=0-666%48`,
  no `--overwrite` → skip-existing). **Test:** `sorcha_gen_maps_grid_test.sh` (idx 333 antisun,
  652 lat50). **Manifest:** `slurm/grid_map_manifest.csv` (667 rows).
- **Timing:** ~3.6 min (sparse high-lat) to ~10 min (dense antisun); full grid ~2–3 h at %48.
- **Validation (2026-06-16):** antisun map P(NEO) peaks ≈0.93 at vlam≈−0.54 (NEO locus); MBA
  stationary/retrograde ≈0; with cf=5 it reproduces the old monthly antisun map (support ratio
  5→1, mean P(NEO) matches per bin). Per-cell scatter = known joblib non-determinism (~0.02 F1).
- Output: `prob_maps_grid/`.

---

## 9. Step 4 — Phase 2: scoring

**Script:** `neomod/pipeline/sorcha_phase2.py` (commands: `audit`, `score-vdp`, `sample`,
`run-digest2`, `combine`). v5.0 Slurm wrappers to be created:
`sorcha_phase2_vdp_v5.sh`, `sorcha_digest2_v5.sh` (model on the existing
`sorcha_phase2_vdp.sh`/`sorcha_digest2_slurm.sh`, repointed to `outputs/tracklets_v5_grid/`,
`prob_maps_grid/`, `outputs/phase2_v5/`).

### 9a. VDP scoring (`score-vdp`)
Reads tracklet parquets, loads each tracklet's assigned map `.npz`, scores
P(NEO) by **bilinear interpolation** on the (vlam, vbeta) grid for the tracklet's mag bin,
using `mean_dra`/`mean_ddec` (rates from Sorcha, converted by `/cos(Dec)`). Run as a Slurm
array of shards (prior run: 113 shards × 128 parquets, `--array=0-112%32`). Adds
`P_NEO_vdp, vlam, vbeta, mag_bin_label`. Output: `outputs/phase2_*/vdp_shards/`.
GMM run uses `--no-nearest-dist-mask`.

### 9b. digest2 scoring (`run-digest2`)
Formats each tracklet as a pair of MPC 80-column observations (`ra0/dec0/mjd0_utc` and
`ra1/dec1/mjd1_utc`, observed positions, UTC) and runs the `mpcdev-digest2` binary.
- **Designation bug (fixed):** must use `"     D{i:06d}"` (5 spaces + D + 6 digits). The old
  `D{i:011d}` put D at col 1, so digest2's 5-char truncation collapsed all tracklets in a chunk
  to one id → they merged into one object → digest2 hung. Lookup key = `D{i:06d}`.
- Slurm: chunks of ~5,000 tracklets; slow ckpt nodes can time out → retry with
  `--digest2-chunk-tracklets 1000`. Output: `outputs/phase2_*/digest2_shards/`.
- **Threshold:** Wagg submits at digest2 ≥ 65 (0–100 scale) = 0.65 on our 0–1 scale.

### 9c. combine
```
conda_prep/bin/python neomod/pipeline/sorcha_phase2.py combine \
  --work-dir outputs/phase2_v5 --output outputs/phase2_v5/sorcha_comparison_v5.parquet
```
**CRITICAL gotcha:** digest2 shards embed VDP columns from whatever run produced them. When
combining a *new* VDP run, drop stale `P_NEO_vdp,vlam,vbeta,mag_bin_label` from the d2 shards,
then left-join the fresh `vdp_shards` on `tracklet_id`. Final parquet → SCP to Arnor.

### digest2 binary validation (2026-06-08)
`validate_digest2_neocp.py` vs MPC published NEOCP scores (real 80-col tracklets through
`mpcdev-digest2`): Pearson r=0.627, Spearman 0.586, MAE 0.101, 12/13 high-MPC(≥0.90) scored
≥0.50. Correlation limited by input differences (MPC full arc vs our 1-hr synthetic). Binary
confirmed functioning. Used synthetic 5-char keys `VA000`–`VA0NN` to avoid truncation collisions.

---

## 10. Step 5 — ROC analysis (on Arnor)

Notebook `sorcha_roc_comparison.ipynb`. NaN P(NEO) (outside footprint) **dropped** for a fair
in-footprint comparison (digest2 has no footprint). Sweep both classifiers; report AUC, best
F1, completeness, contamination, per-population fractions. Fill Section 4.8 of
`neomod/paper/NEOrocks.tex` (Overleaf) with v5.0 numbers.

### Prior-run results (context; v5.0 pending)
Full 2-yr (v3.3 window), in-footprint, GMM + mask OFF + grid ±2.0:
| Classifier | F1 | Completeness | Contamination |
|---|---|---|---|
| VDP (GMM) | **0.837** | 78.3% | 10.2% |
| digest2 | 0.836 | 77.1% | 8.7% |
Per-population "above threshold": NEO 78.3% vs 77.3% (VDP completeness slightly *exceeds*
digest2); MBA 1.5% vs 0.2% (VDP's remaining gap is MBA contamination, 7.5× more false MBAs);
TNO 0.0% vs 69.5% (digest2 can't separate slow TNOs from NEOs; VDP can); Trojan 0.1% vs 7.6%.
F1 progression (S3M kNN → fixes): 0.342 → … → 0.740 (S3M) → 0.803 (hybrid kNN) → **0.837**
(GMM, tied/slightly beating digest2). The headline scientific point: VDP, calibrated to the
scored population, matches or beats digest2 while using richer 2-D velocity structure.

---

## 11. Conventions, constants, reproducibility

- **Frame:** heliocentric ecliptic J2000 throughout. Obliquity 23.439291° (fixed) in the fast
  ecliptic helpers.
- **μ_sun:** `0.01720209895²` AU³/day² (Gaussian gravitational constant squared).
- **Velocity components:** vlam = ecliptic-longitude rate, vbeta = ecliptic-latitude rate
  (deg/day), in the antisun-relative frame.
- **Reproducibility:** GMM `random_state=42`, pipeline `seed=42`. BUT joblib parallel density
  evaluation is **not** bit-reproducible (thread-scheduling-dependent float summation) → F1
  drifts ~0.02 between regenerations of the *same* map. **For final paper maps, regenerate with
  `--n-jobs 1`** (≈16× slower) if exact reproducibility is required. The comparison run uses
  `--n-jobs 16`.
- **NEOMOD3 augmentation:** attempted (adds debiased NEO orbits to GMM training), **zero net
  effect** (code left in place, harmless) — unusual velocities require specific orbital phase at
  opposition, rare regardless of training diversity.

---

## 12. File / directory inventory (v5.0)

```
inputs/hybrid_sorcha_orbits.csv, hybrid_sorcha_phys.csv     # Sorcha inputs
CDS_colors.parquet, hybrid.h5, hybrid_elements.parquet      # catalog + colors + Keplerian
baseline_v5.0.0_2yrs.db                                     # pointing DB (148 MB)
sorcha_cache_2025-07-06/                                    # SPICE cache
outputs/production_2yr_v5/   inst*_part*.h5  (14,445; ~232 GB)
outputs/tracklets_v5/        monthly-map Phase-1 parquets (29,844,550 tracklets)
outputs/tracklets_v5_grid/   grid-map Phase-1 re-run (target of Step 5)
prob_maps_gmm/               24 monthly antisun GMM maps (superseded by grid)
prob_maps_grid/              667 antisun-relative grid maps (Step 3)
outputs/phase2_v5/           VDP + digest2 shards + sorcha_comparison_v5.parquet

neomod/
  pipeline/
    sorcha_postprocess.py            # Phase 1 (grid-aware assignment, n_det_per_night)
    sorcha_gen_maps_grid.py          # 667-map grid generator (NEW)
    sorcha_gen_map_gmm.py            # single monthly map (GMM)
    sorcha_phase2.py                 # Phase 2 (score-vdp/run-digest2/combine)
    validate_digest2_neocp.py        # digest2 binary validation
    hybrid_catalog_prep.py           # Cartesian->Keplerian
    config/Rubin_full_footprint_detections.ini
    slurm/
      multi_sorcha_production_v5.sh
      split820v5_part003_final.sh    # + bisect820v5*.sh (pathological-object hunt)
      sorcha_postprocess_v5.sh, sorcha_postprocess_v5_grid.sh
      sorcha_gen_maps_grid_slurm.sh, sorcha_gen_maps_grid_test.sh
      grid_map_manifest.csv
      sorcha_gen_maps_gmm_slurm.sh, sorcha_phase2_vdp*.sh, sorcha_digest2_*.sh
  src/
    velocity_density_pipeline_gmm.py # GMM pipeline (NEO=GMM, others=K|M; + norm fix)
    velocity_density_pipeline_fast.py / _hybrid.py   # earlier variants (untouched)
    hybrid_loader.py, s3m_loader.py, neoscore.py
  adam_core_stub/                    # minimal stub for vdp import
  docs/                              # context docs (this file's companions)
  paper/  -> Overleaf (gitignored)   # NEOrocks.tex lives here, NOT in git
```

---

*Compiled 2026-06-16 for paper methods writing. Keep in sync with `docs/SORCHA_V5_PIPELINE.md`
(operational reference) and `docs/WAGG_SORCHA_HYAK_CONTEXT.md` (full history). v5.0 Phase 2
results pending the 667-map grid run.*
