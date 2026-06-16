# Sorcha v5.0 Simulation Pipeline
**Working directory:** `/mmfs1/gscratch/dirac/ds2004/sorcha/`
**Cadence database:** Rubin Operations Baseline v5.0.0, trimmed to 2 years (MJD 60980–61710, 2025-Nov to 2027-Nov)
**Cluster:** Hyak/Klone HPC, partition `ckpt`, account `astro`

---

## Input Files

| File | Size | Description |
|------|------|-------------|
| `inputs/hybrid_sorcha_orbits.csv` | 2.0 GB | 14,444,912 objects, Cartesian state vectors (AU, AU/day) at `epochMJD_TDB`, columns: `ObjID, FORMAT, x, y, z, xdot, ydot, zdot, epochMJD_TDB` |
| `inputs/hybrid_sorcha_phys.csv` | 1.8 GB | 14,444,912 rows, physical parameters: `ObjID, H_r, GS, u-r, g-r, i-r, z-r, y-r` |
| `baseline_v5.0.0_2yrs.db` | 148 MB | SQLite3 pointing database, 414,488 observations, MJD 60980.0–61709.4. Trimmed from the full 10-year `baseline_v5.0.0_10yrs.db` (728 MB) via `CREATE TABLE AS SELECT WHERE observationStartMJD < 61710.0` |
| `neomod/pipeline/config/Rubin_full_footprint_detections.ini` | — | Sorcha configuration (see details below) |
| `sorcha_cache_2025-07-06/` | — | ASSIST+REBOUND SPICE kernel cache, pre-fetched to avoid repeated downloads |

### Population composition (`inputs/hybrid_sorcha_orbits.csv`)
The hybrid catalog combines:
- **S3M** (Synthetic Solar System Model): TNOs, Centaurs, MBAs, Trojans
- **MPCORB**: real known NEOs and MBAs from the MPC catalog

Objects are classified post-hoc in Phase 1 by semi-major axis and eccentricity into: NEO (q < 1.3 AU), MBA (2.0 < a < 3.3 AU, e < 0.3), Trojan (4.8 < a < 5.6 AU), TNO (a > 30 AU), other.

### Sorcha config key settings (`Rubin_full_footprint_detections.ini`)
- `ephemerides_type = ar` — ASSIST+REBOUND integrator for ephemeris generation
- `ar_obs_code = X05` — Vera C. Rubin Observatory MPC code
- `ar_healpix_order = 6` — spatial discretisation for field-of-view pre-filtering
- `camera_model = footprint` — full LSSTCam CCD footprint with chip gaps (not circular approximation)
- `footprint_edge_threshold = 2.` arcsec — exclude detections within 2 arcsec of detector edges
- `output_format = hdf5` — one HDF5 file per chunk-part (key: `sorcha_results`)
- `output_columns = all` — all 56 Sorcha output columns retained
- `observing_filters = r,g,i,z,u,y`
- `bright_limit = 16.0` — saturated sources removed
- `phase_function = HG`
- No `[LINKINGFILTER]` section — all individual detections kept (no linking/tracklet filtering at this stage)

---

## Step 1 — Sorcha Production Run

### Script
**`production_run/multi_sorcha_production.py`** — Python parallel worker  
**`neomod/pipeline/slurm/multi_sorcha_production_v5.sh`** — Slurm submission script

### Chunking strategy
The 14,444,912-row input is processed in 904 instances (array tasks 0–903):
- Each instance reads a 16,000-row slice: `rows = instance_id * 16000 : (instance_id+1) * 16000`
- Each 16,000-row slice is split into 16 parts of ≤1,000 objects each
- All 16 parts run simultaneously via Python `multiprocessing.Pool(processes=16)`
- Each part calls `sorcha-run` as a subprocess → writes one `inst{NNNNN}_part{NNN}.h5`

### Slurm parameters
```
--partition=ckpt  --account=astro
--cpus-per-task=16  --mem=128G  --time=04:00:00
--array=0-903%16   (16 instances running simultaneously)
```

### Output
```
outputs/production_2yr_v5/inst{NNNNN}_part{NNN}.h5
```
- **Total files:** 14,445 (instances 0–901 produce 16 parts each = 14,432; instance 902 produces 13 parts = 12,912 objects; instance 903 exits immediately with 0 rows)
- **Format:** HDF5, key `/sorcha_results`, 56 columns including RA/Dec, rates, PSFMag, heliocentric state vectors, filter, SNR
- **Total size:** ~232 GB

### Skip-existing logic
`multi_sorcha_production.py` checks `h5_is_good()` before each `sorcha-run` call. If the output HDF5 exists and contains the `sorcha_results` key it is skipped. This allows safe resubmission after preemption without `--requeue`.

### Work directory
```
work/production_2yr_v5/instance_{NNNNN}/
    orbits_{NNN}.csv     # per-part orbit slice
    physical_{NNN}.csv   # per-part physical slice
```
These CSV slices are written by `multi_sorcha_production.py` before spawning the pool.

---

## Step 1b — Recovery: Pathological Objects in inst00820_part003

Instance 820, part 003 (1,000 objects) hung the ASSIST+REBOUND integrator indefinitely. This is caused by objects with close planetary encounters that make the adaptive step-size integrator take astronomically small steps.

### Excluded objects (5 total from original part003)
| ObjID | Reason |
|-------|--------|
| `A804 RA` | Integrator hung — found in v3.3 run |
| `A854 OA` | Integrator hung — found in v3.3 run |
| `A868 TA` | Integrator hung — found in v3.3 run |
| `A854 RA` | Integrator hung — found via bisection in v5.0 run |
| `A868 WA` | Integrator hung — found via bisection in v5.0 run |

The new objects (`A854 RA`, `A868 WA`) were not problematic in v3.3 — the 6-month offset in the v5.0 simulation window (Nov 2025 vs May 2025 start) places them at different orbital phases, likely triggering close encounters that did not occur in the v3.3 window.

### Bisection method
1. Split 997 remaining objects (after removing v3.3 known bad objects) into 16 batches of ~63 each. Submit as Slurm array with 15-min limit. Batches that timeout contain bad objects.
2. Split each failing batch into 8-object sub-batches. Repeat.
3. Run individual objects from remaining candidates. Timed-out single-object runs are the bad objects.

### Recovery script
**`neomod/pipeline/slurm/split820v5_part003_final.sh`**  
Reads `work/production_2yr_v5/instance_00820/orbits_003_skip5.csv` (995 objects) and runs `sorcha-run` directly, outputting `outputs/production_2yr_v5/inst00820_part003.h5`.

**Filtered input files:**
- `work/production_2yr_v5/instance_00820/orbits_003_skip5.csv` (995 objects)
- `work/production_2yr_v5/instance_00820/physical_003_skip5.csv` (995 objects)

---

## Step 2 — Phase 1: Tracklet Building

Converts raw Sorcha HDF5 detections into nightly two-detection tracklets. One parquet file per input HDF5.

### Script
**`neomod/pipeline/sorcha_postprocess.py`** — Python worker  
**`neomod/pipeline/slurm/sorcha_postprocess_v5.sh`** — Slurm submission script

### Tracklet construction logic
For each `inst*_part*.h5` file:
1. Read all detections from `/sorcha_results` key
2. Filter: remove detections with non-finite RA rate, Dec rate, or magnitude
3. Group by `(ObjID, night)` where `night = floor(fieldMJD_TAI)`
4. For each group, extract the first and last detection of the night
5. Keep groups where: `n_det >= 2` AND `3 min ≤ dt ≤ 90 min` (Wagg-style pairing window)
6. Compute per-tracklet quantities: mean RA/Dec, mean rates, `mjd0/mjd1`, ecliptic coordinates, orbital elements, population classification
7. Assign to nearest VDP probability map by epoch proximity (antisun footprint check)

### Output columns (44 total)
`tracklet_id, source_file, source_file_index, ObjID, night, n_det, night_span_min, dt_min, mean_ra, mean_dec, mean_dra, mean_ddec, mean_mag, ra0, dec0, mjd0_tai, mjd0_utc, mag0, filter0, snr0, ra1, dec1, mjd1_tai, mjd1_utc, mag1, filter1, snr1, H_r, H_filter, range_km, obj_sun_km, ecl_lon, ecl_lat, prob_map, prob_map_file, prob_map_obstime_str, prob_map_center_lon_deg, prob_map_center_lat_deg, prob_map_dist_deg, population, a_au, e, q_au, n_det_per_night`

**`n_det_per_night`** (column 44): total raw Sorcha detections for that (ObjID, night) pair before the first/last pairing. Allows re-analysis at ≥3 or ≥4 detection thresholds (Wagg uses ≥3) without rerunning Phase 1.

### Map assignment
Each tracklet is assigned to the VDP probability map whose epoch is nearest to the observation time, subject to the tracklet being within `max_sep_deg` of the **current antisun direction** at observation time (not the historical map center position). Tracklets outside all map footprints are dropped (unless `--keep-unmapped`).

Probability map directory: `prob_maps/` (24 monthly antisun maps, May 2025 – Apr 2027)
→ *Will be replaced by the ~500-map antisun-relative grid (Step 3) before Phase 2 scoring.*

### Slurm parameters
```
--partition=ckpt  --account=astro
--cpus-per-task=1  --mem=8G  --time=01:00:00
--array=0-999%64   (64 tasks running simultaneously per batch)
```

### Submission batching
14,445 files exceed the per-user QOS limit of 2,000 pending tasks. Submitted in 15 batches using `INDEX_OFFSET`:

| Batch | `INDEX_OFFSET` | `--array` | Files covered |
|-------|---------------|-----------|---------------|
| 0 | 0 | 0-999 | 0–999 |
| 1 | 1000 | 0-999 | 1000–1999 |
| 2 | 2000 | 0-999 | 2000–2999 |
| 3 | 3000 | 0-999 | 3000–3999 |
| 4 | 4000 | 0-999 | 4000–4999 |
| 5 | 5000 | 0-999 | 5000–5999 |
| 6 | 6000 | 0-999 | 6000–6999 |
| 7 | 7000 | 0-999 | 7000–7999 |
| 8 | 8000 | 0-999 | 8000–8999 |
| 9 | 9000 | 0-999 | 9000–9999 |
| 10 | 10000 | 0-999 | 10000–10999 |
| 11 | 11000 | 0-999 | 11000–11999 |
| 12 | 12000 | 0-999 | 12000–12999 |
| 13 | 13000 | 0-999 | 13000–13999 |
| 14 | 14000 | 0-444 | 14000–14444 |

Submit 2 batches at a time (= 2,000 tasks) and add more as the queue drains.

### Output
```
outputs/tracklets_v5/tracklets_inst{NNNNN}_part{NNN}.parquet
```
One parquet per input HDF5. Concatenated across all 14,445 files in Phase 2.

---

## Step 3 — Antisun-Relative Sky Grid (667 maps) — AS BUILT (2026-06-16)

### Motivation
The previous 24 monthly maps all sit at ecliptic lat=0, covering only the antisun
direction along the ecliptic. A full sky grid in **antisun-relative (Δlon, lat)
coordinates** lets any LSST tracklet — regardless of sky position — be scored with
the nearest map. The grid is time-independent: NEO/MBA velocity structure is set by
geometry relative to the Sun, so a map built at offset (Δlon, lat) from the antisun
is reusable at any epoch as the Sun moves.

### Grid specification (as built)
- **Coordinate system:** (Δlon from current antisun, ecliptic lat) — antisun-relative ecliptic
- **Longitude:** 10° steps; centers kept where |Δlon| ≤ 140° (i.e. ≥ 40° from the Sun,
  which sits at |Δlon| = 180°) → **29 usable longitudes** (Δlon = −140°…+140°)
- **Latitude (non-uniform, fine near the ecliptic):** symmetric expansion of
  `0, 1, 2, 3, 4, 5, 8, 12, 18, 25, 35, 50` → **23 values**
  (0, ±1, ±2, ±3, ±4, ±5, ±8, ±12, ±18, ±25, ±35, ±50°). 1° steps at the ecliptic
  where most NEOs/MBAs lie, coarsening to 10–15° steps at high latitude.
- **Sun exclusion:** 40° (no centers within 40° of the Sun ecliptic longitude)
- **Total maps:** **29 × 23 = 667**
- **Reference epoch:** 2026-01-01T00:00:00 (only fixes where patches land on the real
  sky for this generation; velocity statistics are epoch-independent in the antisun frame)
- **MBA clone_factor = 5** (matches the config behind the documented F1≈0.837 result;
  applied via `--mba-clone-factor 5`, the script default, NOT by editing the global
  `DEFAULT_POPULATION_SETTINGS`). NEO=80, TNO=10, Trojans=5 are the code defaults.
- **Velocity grid:** (−2.0, +2.0) deg/day, step 0.01 (401×401), per `velocity_density_pipeline_gmm.py`
- **Training sky cut:** max_sep_deg = 30° per map center

### Script
`neomod/pipeline/sorcha_gen_maps_grid.py` — single configurable generator:
`--lon-step`, `--lat-base`/`--lat-points`, `--sun-exclusion`, `--ref-obstime`,
`--mba-clone-factor` (default 5), `--task-id N` (generate the Nth grid map for a
Slurm array), `--list-only`/`--map-grid-file` (write the grid manifest), `--n-jobs`,
`--save-overlays`. Same GMM map-generation path as `sorcha_gen_map_gmm.py`, but at the
specified (Δlon, lat) offset from the antisun; afterwards each `.npz` is augmented with
exact `delta_lon_from_antisun_deg` and `grid_lat_deg` keys for assignment.

**Slurm:** `neomod/pipeline/slurm/sorcha_gen_maps_grid_slurm.sh` —
`--cpus-per-task=16 --mem=32G --time=06:00:00 --array=0-666%48`, calls the script with
`--task-id $SLURM_ARRAY_TASK_ID`. No `--overwrite` → skip-existing, so resubmitting
after preemption/maintenance only fills unfinished maps.
**Test script:** `sorcha_gen_maps_grid_test.sh` (indices 333 = antisun, 652 = lat 50;
`--save-overlays`).
**Manifest:** `neomod/pipeline/slurm/grid_map_manifest.csv` (667 rows: index, Δlon, lat, filename).

### Measured performance
- Per-map runtime: ~3.6 min (sparse high-lat) to ~10 min (dense antisun).
- Full grid at %48 concurrency: ~2–3 h wall (667 tasks, under the 2000 QOS cap and ~377 group ckpt limit).

### Validation (test maps, 2026-06-16)
- Map 333 (antisun) P(NEO) is physically correct: peak ≈0.93 at vlam≈−0.54 (fast
  retrograde NEO locus); MBA stationary (0,0) and retrograde (−0.18,0) score ≈0.
- With MBA cf=5, map 333 reproduces the old `prob_maps_gmm/2026-01-01_antisun.npz`
  (MBA support ratio 5.0→1.0; mean P(NEO) agrees per mag bin). Confirms the grid
  pipeline matches the F1≈0.837 configuration. Residual per-cell scatter is the known
  joblib parallel non-determinism (~0.02 F1 between regenerations).

### Map-to-tracklet assignment (in `sorcha_postprocess.py`)
Grid maps are auto-detected by the `delta_lon_from_antisun_deg` key. For each tracklet:
compute the antisun ecliptic longitude at observation time, take Δlon = ecl_lon −
antisun_lon (wrapped to ±180°) and ecliptic lat, then assign to the **nearest grid
center** in 2-D (Δlon, lat) space. Tracklets within 40° of the Sun are left unassigned.

### Required sanity plots per batch (advisor)
- Heliocentric distance histogram (AU, split by population: NEO/MBA/TNO/Trojan)
- Heliocentric x-y scatter (1–5 AU range, showing inner belt structure)
- Velocity (vlam, vbeta) coverage plot — derivable from `--save-overlays` outputs

Output directory: `prob_maps_grid/`

---

## Step 4 — Phase 2: Scoring

*On hold until Step 3 (500-map grid) is complete.*

### Phase 2a — VDP scoring
**Script:** `neomod/pipeline/sorcha_phase2.py score-vdp`  
**Slurm:** new v5 version of `neomod/pipeline/slurm/sorcha_phase2_vdp.sh` (to be created)

Reads all tracklet parquets from `outputs/tracklets_v5/`, looks up each tracklet's assigned map file, loads the VDP probability grid, and assigns `P(NEO)` via bilinear interpolation on the (vlam, vbeta) grid.

Output: `outputs/phase2_v5/vdp_shards/`

### Phase 2b — digest2 scoring
**Script:** `neomod/pipeline/sorcha_phase2.py run-digest2`  
**Slurm:** new v5 version of `neomod/pipeline/slurm/sorcha_digest2_slurm.sh` (to be created)

Formats each tracklet as a pair of MPC 80-column observations, runs the `mpcdev-digest2` binary, and stores the raw digest2 score.

Output: `outputs/phase2_v5/digest2_shards/`

### Phase 2c — Combine
```bash
conda_prep/bin/python neomod/pipeline/sorcha_phase2.py combine \
  --work-dir outputs/phase2_v5 \
  --output outputs/phase2_v5/sorcha_comparison_v5.parquet
```
Merges VDP and digest2 shards into a single parquet for ROC analysis on Arnor.

---

## Key Technical Notes

### Why no `--requeue`
On Hyak `ckpt` partition, preempted jobs with `--requeue` fail to re-fetch `--export=ALL` environment variables. They get stuck in state `"user env retrieval failed requeued held"`. Instead: use skip-existing logic (`h5_is_good()` in the Python worker) and resubmit manually.

### Cadence database trimming
`baseline_v5.0.0_2yrs.db` was created from the full 10-year database:
```bash
sqlite3 baseline_v5.0.0_10yrs.db "
  ATTACH DATABASE 'baseline_v5.0.0_2yrs.db' AS new;
  CREATE TABLE new.observations AS
    SELECT * FROM observations WHERE observationStartMJD < 61710.0;
  CREATE TABLE new.info AS SELECT * FROM info;
  CREATE TABLE new.events AS SELECT * FROM events;
"
```
MJD 60980 = 2025-Nov-04, MJD 61710 = 2027-Nov-08. This 730-day window is 6 months later than the v3.3 window (May 2025 – Apr 2027) because baseline_v5.0.0 starts its simulated survey in November 2025.

### v3.3 vs v5.0 simulation windows
| Run | Cadence DB | Window | Notes |
|-----|-----------|--------|-------|
| v3.3 | `baseline_v3.3_2yrs.db` | 2025-May – 2027-Apr | Previous run, outputs in `outputs/production_2yr/` |
| v5.0 | `baseline_v5.0.0_2yrs.db` | 2025-Nov – 2027-Nov | Current run, outputs in `outputs/production_2yr_v5/` |

### Conda environment
All Python/Sorcha commands use the project-local environment:
```
conda_prep/bin/python
conda_prep/bin/sorcha-run
```
The shell always shows `(base)` — this is NOT the active environment. Always use the explicit path above.

### QOS limit
`QOSMaxSubmitJobPerUserLimit = 2000` tasks (pending + running combined) on `ckpt`. Each Phase 1 batch = 1,000 tasks. Maximum 2 batches submitted simultaneously.
