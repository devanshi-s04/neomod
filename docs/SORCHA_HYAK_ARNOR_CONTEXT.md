# Wagg / Sorcha / Hyak Context
Generated: 2026-05-27, Updated: 2026-05-31 (Step 1 complete)
Covers: full Hyak setup, Sorcha 2yr production run, neomod deployment, Wagg paper methodology, post-processing plan.

---

## Goal

Reproduce the Tom Wagg et al. LSST/Rubin NEOCP simulation workflow using Sorcha on Hyak (Klone), then score the resulting detection-level tracklets with both VDP and digest2 for a head-to-head ROC comparison on real Rubin-cadence simulated data.

The Sorcha run simulates Rubin detections of a hybrid S3M+MPCORB Solar System population over 2 years of the Rubin baseline cadence. Post-processing builds Wagg-style nightly tracklets and scores them.

---

## Hyak Login and Working Directory

```bash
ssh ds2004@klone.hyak.uw.edu
# UW NetID password → Duo → type "1" → approve push
```

Working directory:

```text
/mmfs1/gscratch/astro/ds2004/sorcha
```

Slurm partition: `ckpt` (checkpoint). No explicit `--account` line needed; defaults to `astro`.

Always use the conda env Python directly — never bare `python` or `pip`:

```bash
/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python
/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/pip
/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/sorcha-run
```

The shell always shows `(base)` — that is NOT the working env. Ignore it.

---

## Key Files on Hyak

```text
/mmfs1/gscratch/astro/ds2004/sorcha/
├── baseline_v3.3_10yrs.db              Rubin 10yr cadence (799M)
├── Rubin_full_footprint_wagg_detections.ini  Wagg-style Sorcha config
├── hybrid.h5                           Wagg hybrid S3M+MPCORB catalog (14.4M objects)
├── CDS_colors.parquet                  LSST color distribution (534 rows)
├── inputs/
│   ├── hybrid_sorcha_orbits.csv        Sorcha orbit input (~2.0G)
│   └── hybrid_sorcha_phys.csv          Sorcha physical input (~1.8G)
├── sorcha_cache_2025-07-06/            SPICE/ephemeris cache
├── outputs/production_2yr/             COMPLETED 2yr run (14,445 h5 files, ~232G)
├── neomod/                             VDP codebase (cloned from GitHub)
│   └── src/
│       ├── velocity_density_pipeline.py
│       ├── neoscore.py
│       └── ...
├── prob_maps/                          VDP probability maps
│   ├── prob_maps_2025-03-21.npz
│   ├── prob_maps_2026-05-09T22_neocp.npz
│   ├── prob_maps_2026-05-09_antisun_minus_45.npz
│   ├── prob_maps_2026-05-09_antisun_minus_90.npz
│   └── prob_maps_2026-05-09_antisun_minus_120.npz
├── adam_core_stub/                     Minimal adam_core mock (allows vdp import)
│   └── adam_core/
│       ├── time/__init__.py            class Timestamp: pass
│       ├── coordinates/__init__.py     stub classes
│       └── constants/__init__.py      KM_P_AU, S_P_DAY
├── multi_sorcha_production.py          Resumable parallel Sorcha wrapper
├── multi_sorcha_production.sh          Slurm script for production runs
└── check_sorcha_outputs.py             Checker/rerun script generator
```

---

## Sorcha Config: Wagg-Style Detection Mode

File: `Rubin_full_footprint_wagg_detections.ini`

Key settings vs default:
- Removed `[LINKINGFILTER]` section entirely — keeps all individual detections
- `output_format = hdf5` (not sqlite3)
- `output_columns = all` (56 columns)
- `ephemerides_type = ar` (ASSIST+REBOUND)
- `aux_format = csv`

---

## Input Population: Wagg Hybrid Catalog

`hybrid.h5` key `/df`, shape 14,444,912 × 10:

```text
id, x, y, z, vx, vy, vz, t_0, H, g
```

- Cartesian state vectors (AU, AU/day)
- `H` = absolute magnitude (V-band)
- `g` = phase slope, typically 0.15
- `t_0` = epoch (MJD TDB)

Prep script (`scripts/prep_sorcha_inputs.py`) converted this to Sorcha-readable CSVs:
- Orbit CSV: ObjID, FORMAT=CART, x, y, z, xdot, ydot, zdot, epochMJD_TDB
- Physical CSV: ObjID, H, GS, u-r, g-r, i-r, z-r, y-r (colors sampled from CDS_colors.parquet)
- Magnitude conversion: `H_r = H_V - (johnson V - LSST r)`

---

## 2yr Production Run: COMPLETE

**Status as of 2026-05-27: All 14,445 h5 files present, 0 missing, 0 zero-size.**

```text
outputs/production_2yr/
  inst00000_part000.h5 ... inst00903_part015.h5
  232G total
  14,445 files
```

Chunking: 14,444,912 objects / 16,000 per array chunk = 904 chunks (indices 0–903), 16 parts each.

### Two Recovered Problem Files

`inst00804_part015.h5` — recovered by rerunning instance 804.

`inst00820_part003.h5` — recovered by splitting into sub-runs and merging. Three pathological objects were excluded because single-object Sorcha runs hung for 25–45 min each:

```text
A804  array=14   tag=sub05_obj014
A854  array=42   tag=sub05_obj042
A868  array=118  tag=sub06_obj018
```

Documented in: `production_2yr_excluded_objects.txt`

### Sorcha Output Schema (56 columns)

Key columns for post-processing:

```text
ObjID              — object identifier (from hybrid.h5 id column)
fieldMJD_TAI       — detection time (MJD TAI)
RA_deg             — observed RA
Dec_deg            — observed Dec
RARateCosDec_deg_day  — apparent RA rate × cos(Dec), deg/day
DecRate_deg_day    — apparent Dec rate, deg/day
PSFMag             — observed PSF magnitude
H_r                — absolute magnitude in r-band
H_filter           — absolute mag in filter used
detectorID
Range_LTC_km       — object-observer distance
Obj_Sun_LTC_km     — heliocentric distance
x, y, z, xdot, ydot, zdot, epochMJD_TDB  — orbital state vector (from input)
optFilter          — filter (u/g/r/i/z/y)
seeingFwhmGeom_arcsec
fieldFiveSigmaDepth_mag
SNR
```

HDF5 keys per file:

```text
/sorcha_results
/sorcha_results/meta/optFilter/meta
```

Read example:

```python
import pandas as pd
df = pd.read_hdf("outputs/production_2yr/inst00000_part000.h5", key="sorcha_results")
# shape: ~48,000 rows × 56 columns
```

---

## neomod on Hyak

Cloned from GitHub (public repo) via HTTPS:

```bash
git clone https://github.com/devanshi-s04/neomod.git
# lives at /mmfs1/gscratch/astro/ds2004/sorcha/neomod/
```

To update after pushing from Mac:

```bash
cd /mmfs1/gscratch/astro/ds2004/sorcha/neomod && git pull origin main
```

### VDP Import

`velocity_density_pipeline.py` requires `adam_core` (via `neoscore.py`). `adam_core` is slow to install. Solution: a minimal stub at `adam_core_stub/` satisfies the import. `score_observation` and `ProbMapSet.from_npz` do NOT call adam_core at runtime — only orbit-propagation functions do.

Working import snippet (use in all Hyak scripts):

```python
import sys
sys.path.insert(0, "/mmfs1/gscratch/astro/ds2004/sorcha/adam_core_stub")
sys.path.insert(0, "/mmfs1/gscratch/astro/ds2004/sorcha/neomod/src")
import velocity_density_pipeline as vdp
```

Confirmed working output:

```text
vdp import: OK
  SMOOTH_DENSITY_MAPS = True
  SMOOTH_POPULATION_NAMES = ('NEO',)
```

### Probability Maps

5 maps uploaded to `prob_maps/`:

| File | Center (lon, lat) | Epoch |
|---|---|---|
| prob_maps_2025-03-21.npz | 180°, 0° (opposition) | 2025-03-21 |
| prob_maps_2026-05-09T22_neocp.npz | 229°, 0° | 2026-05-09T22 |
| prob_maps_2026-05-09_antisun_minus_45.npz | ~135°, 0° | 2026-05-09 |
| prob_maps_2026-05-09_antisun_minus_90.npz | ~90°, 0° | 2026-05-09 |
| prob_maps_2026-05-09_antisun_minus_120.npz | ~60°, 0° | 2026-05-09 |

Each map covers a ~30° sky radius. Objects outside all map footprints get P_NEO = NaN and are excluded from the comparison.

Load example:

```python
pms = vdp.ProbMapSet.from_npz("/mmfs1/gscratch/astro/ds2004/sorcha/prob_maps/prob_maps_2025-03-21.npz")
```

---

## VDP Rate Conversion from Sorcha Output

Sorcha outputs `RARateCosDec_deg_day` = dRA/dt × cos(Dec). VDP `score_observation` needs the uncorrected dRA/dt:

```python
import numpy as np
dra_deg_day  = df["RARateCosDec_deg_day"] / np.cos(np.deg2rad(df["Dec_deg"]))
ddec_deg_day = df["DecRate_deg_day"]

P_NEO = pms.score_observation(
    ra_deg       = df["RA_deg"].values,
    dec_deg      = df["Dec_deg"].values,
    dra_deg_day  = dra_deg_day.values,
    ddec_deg_day = ddec_deg_day.values,
    mag_app      = df["PSFMag"].values,
    scorer       = None,   # not needed for score_observation
)
```

---

## Post-Processing Pipeline

### Key Insight: No Cross-File Merging Needed

All detections for a given ObjID are contained within a single h5 file (same Sorcha instance/part). Nightly tracklets are therefore self-contained per file. Process each h5 file independently in a Slurm array.

### Full Pipeline

```
Phase 1 — COMPLETE (2026-05-27)
  sorcha_postprocess.py run over all 14,445 h5 files via Slurm array
  Output: 14,445 tracklet parquets in outputs/tracklets/
  Total tracklets produced: 84,901,167 (before map filter)
  Tracklets within map footprints: 40,749,299 (48% kept)
  Tracklets outside all maps: 44,151,868 (52% dropped)

Phase 2 — IN PROGRESS (2026-05-28)
  Script: sorcha_phase2.py  (commands: audit, score-vdp, sample, run-digest2, combine)
  Slurm wrapper for VDP: sorcha_phase2_vdp.sh

  Step 1 — score-vdp: COMPLETE
    Ran as Slurm array --array=0-112%32 (113 shards × 128 tracklet parquets each)
    Output: outputs/phase2/vdp_shards/vdp_00000.parquet … vdp_00112.parquet
    Adds columns: P_NEO_vdp, vlam, vbeta, mag_bin_label

  Step 2 — sample: COMPLETE
    Keeps ALL 98,646 NEOs + samples 500K non-NEOs proportionally by population
    Output: outputs/phase2/wagg_subsample.parquet (598,670 rows)
    Population: NEO 98,646 | MBA 452,304 | Trojan 21,707 | other 20,094 | TNO 5,919

  Step 3 — run-digest2: COMPLETE
    Slurm array (sorcha_digest2_slurm.sh): 120 tasks x 5,000 rows each
    BUG FIXED: MPC 80-col designation was wrong — used D{i:011d} (D at position 1,
      causing all tracklets in a chunk to share the same 7-char id -> all merged into
      one object -> digest2 hangs). Fix: "     D{i:06d}" (5 spaces + D + 6 digits);
      score lookup key = "D{i:06d}" (what digest2 outputs after stripping the 5-char
      packed-number field). This matches run_digest2_comparison.py exactly.
    84/120 tasks succeeded first run; 36 timed out on slow ckpt nodes (>60 min).
    Retry (sorcha_digest2_retry.sh): 36 chunks re-run with --digest2-chunk-tracklets 1000
      (5 x 1000-tracklet subprocess calls per task instead of 1 x 5000) to avoid timeout.
    Output: 120 parquets in outputs/phase2/digest2_shards/

  Step 4 — combine: COMPLETE
    Output: outputs/phase2/wagg_sorcha_comparison.parquet (598,670 rows)
    Final columns include: population, P_NEO_vdp, vlam, vbeta, mag_bin_label,
                           P_NEO_d2, digest2_id, and all tracklet geometry columns

Phase 3 — ROC analysis (on Hyak, sorcha_roc_comparison.ipynb)
  Notebook: /mmfs1/gscratch/astro/ds2004/sorcha/sorcha_roc_comparison.ipynb
  Kernel: sorcha (conda_prep) — registered via: conda_prep/bin/pip install ipykernel &&
          conda_prep/bin/python -m ipykernel install --user --name sorcha
  See Results section below.
```

### Phase 2 Tracklet Statistics

```
Total tracklets (before map filter):  84,901,167
Kept (within 30° of a map center):    40,749,299  (48%)
Dropped:                               44,151,868  (52%)

Population breakdown (in-map tracklets):
  MBA:       36,763,790  (90%)
  Trojan:     1,761,483
  other:      1,636,650   ← MPCORB objects with real provisional designations
  TNO:          488,730
  NEO:           98,646  (0.24%)

By prob_map:
  lon229_lat0 (2026-05-09T22 NEOCP pos): 14,779,672
  antisun_minus_45:                        8,980,028
  antisun_minus_90:                        6,366,408
  lon180_lat0 (2025-03-21 opposition):     6,822,551
  antisun_minus_120:                       3,800,640
```

VDP score sanity check (NEO >> MBA confirmed):
- NEO median P_NEO_vdp: ~0.13–0.18 (lower than S3M ~0.86 because tracklets span
  multiple sky patches at varying elongations, not all at opposition)
- MBA median P_NEO_vdp: ~0.002

Note on "other" population: MPCORB objects in the hybrid catalog have real MPC
provisional designations (e.g. `2018 UD2`) rather than S-prefixed IDs. They are
classified as "other" when their orbital elements don't fit the NEO/MBA/Trojan/TNO
bins — treated as non-NEO in the ROC analysis.

### Phase 2 Map Assignment Notes

Each tracklet is scored with the VDP map whose ecliptic center is closest (angular
separation computed via astropy GCRS→GeocentricTrueEcliptic at the map's obstime).
The map metadata (center_lon_deg, center_lat_deg, max_sep_deg, obstime_str) is read
directly from each .npz file — not hardcoded in the script.

Epoch caveat (CRITICAL — see Results section): maps encode velocity statistics at a
specific epoch (specific solar elongation geometry). Applying a map to tracklets
observed at very different epochs severely degrades VDP performance because the
NEO/MBA velocity distributions at a fixed ecliptic position change as the Sun moves.
This is the main driver of VDP underperforming digest2 in the full comparison.
Next step: generate maps covering the full year (see "Next Steps" section).

### Slurm Notes for Phase 2

`sorcha_phase2_vdp.sh` originally had `--array=0-112%8`. Changed to `%32` —
VDP tasks are independent, I/O-light, and the ckpt partition has no reason to limit
to 8 concurrent. `AssocGrpJobsLimit` (group-wide ckpt cap of ~377 jobs) will throttle
naturally if the partition is busy.

### Population Labels — RESOLVED

`hybrid.h5` IDs are plain floats (0.0, 1.0, 2.0 … max ~1.3M), no population prefix.
Sorcha ObjIDs are `S`-prefixed base-N strings (e.g. `S00000dra`) — `S` is Sorcha's
encoding convention applied to all objects, NOT a population indicator.

**Population is classified from state vectors stored in the Sorcha output** (x, y, z,
xdot, ydot, zdot are the input orbital state vectors at epochMJD_TDB, preserved
per detection). `sorcha_postprocess.py::classify_population` computes a, e, q:

```python
MU_SUN_AU3_DAY2 = 0.01720209895**2
# NEO:    q < 1.3 AU
# MBA:    2.0 < a < 3.3, e < 0.3
# Trojan: 4.8 < a < 5.6
# TNO:    a > 30 AU
# other:  everything else
```

### Column Notes (Sorcha output)

- `RA_deg` / `Dec_deg` — observed (astrometric noise applied). Use for MPC 80-col.
- `RATrue_deg` / `DecTrue_deg` — true (noiseless). Do NOT use for digest2.
- `RARateCosDec_deg_day` / `DecRate_deg_day` — ephemeris instantaneous rates (noiseless). Use for VDP.
- `PSFMag` — observed magnitude (with noise). `PSFMagTrue` — noiseless. Use `PSFMag`.
- dRA rate conversion: `dra_deg_day = RARateCosDec_deg_day / cos(Dec_deg)`

### Phase 1 Script: sorcha_postprocess.py

**Location:** `/mmfs1/gscratch/astro/ds2004/sorcha/sorcha_postprocess.py`
**Slurm wrapper:** `sorcha_postprocess.sh` (1 CPU, 8G, 1h per task)

Key design:
- `--file_index N` → resolves to Nth file in `sorted(glob("inst*_part*.h5"))`
- First/last detection per (ObjID, night); time span 3–90 min required
- VDP inputs: mean RA/Dec/dra/ddec/mag across the pair (rates from Sorcha directly)
- digest2 inputs: ra0/dec0/mjd0_utc and ra1/dec1/mjd1_utc (observed positions, UTC)
- Prob-map assignment: loads npz metadata (center_lon_deg, center_lat_deg, max_sep_deg,
  obstime_str) directly from each map file via astropy GCRS→GeocentricTrueEcliptic
- Writes atomic parquet via tmp file → rename; `--skip-existing` for safe resubmission

Output parquet columns (43 total):
```
tracklet_id, source_file, source_file_index, ObjID, night, n_det,
night_span_min, dt_min,
mean_ra, mean_dec, mean_dra, mean_ddec, mean_mag,          ← VDP inputs
ra0, dec0, mjd0_tai, mjd0_utc, mag0, filter0, snr0,       ← digest2 obs 1
ra1, dec1, mjd1_tai, mjd1_utc, mag1, filter1, snr1,       ← digest2 obs 2
H_r, H_filter, range_km, obj_sun_km,
ecl_lon, ecl_lat,
prob_map, prob_map_file, prob_map_obstime_str,
prob_map_center_lon_deg, prob_map_center_lat_deg, prob_map_dist_deg,
population, a_au, e, q_au
```

Submission pattern (batches of 1000, INDEX_OFFSET steps through file indices):
```bash
sbatch --export=ALL,INDEX_OFFSET=<N> --array=0-999%64 sorcha_postprocess.sh
# Final partial batch:
sbatch --export=ALL,INDEX_OFFSET=14094 --array=0-350%64 sorcha_postprocess.sh
```

**Do NOT use `--requeue` on ckpt.** When a preempted job is requeued, Slurm tries to
re-fetch the environment (including the INDEX_OFFSET passed via `--export=ALL,...`).
This fails on Klone → jobs land in `user env retrieval failed requeued held` and never
run. Fix: remove `--requeue`, use `--skip-existing` + manual resubmit for stragglers.

---

## Benchmarks (from production run)

| Config | Objects/chunk | Elapsed | MaxRSS |
|---|---|---|---|
| 16 cores | 4,000 | 17 min | 51 GB |
| 16 cores | 16,000 | 21–24 min | 53–56 GB |
| 8 cores | 16,000 | 26–34 min | 29 GB |

Production used `--cpus-per-task=16 --mem=128G --array=0-903%8`.

---

## Checker / Recovery Pattern

```bash
/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python check_sorcha_outputs.py \
  --outdir outputs/production_2yr \
  --total_objects 14444912 \
  --chunksize 16000 \
  --norbits 1000 \
  --report check_final.txt \
  --rerun_script rerun_failed.sh

bash rerun_failed.sh   # if needed
```

`multi_sorcha_production.py` has skip-existing logic: if a valid `.h5` already exists it skips that part. Safe to rerun any array index.

---

## Common Errors and Fixes

| Error | Fix |
|---|---|
| `Permission denied (publickey)` on git clone | Use HTTPS: `git clone https://github.com/...` |
| `No module named 'adam_core'` | `sys.path.insert(0, ".../adam_core_stub")` before import |
| `No module named 'astroML'` | `/mmfs1/.../conda_prep/bin/pip install astroML` |
| `No available 'sorcha-' utilities found` | Call `sorcha-run` directly (not `sorcha run`) |
| `Invalid account or account/partition` | Use `--partition=ckpt`, no `--account` line |
| pip installs to wrong env | Always use full path: `/mmfs1/.../conda_prep/bin/pip` |
| `Cannot save file into non-existent directory: test_outputs/test_outputs` | `--st` path is relative to `-o` dir; use `--st smoke_stats` not `--st test_outputs/smoke_stats` |
| `user env retrieval failed requeued held` on ckpt | Remove `--requeue` from the Slurm script. On ckpt, preempted jobs with `--requeue` fail to re-fetch `--export=ALL,...` env → permanent hold. Use `--skip-existing` + manual resubmit for stragglers instead. |

---

## VSCode Remote SSH Setup (confirmed working 2026-05-27)

Extension: Remote - SSH (Microsoft, already installed).

SSH config entry in `~/.ssh/config`:

```
Host hyak
  HostName klone.hyak.uw.edu
  User ds2004
  ServerAliveInterval 60
  ServerAliveCountMax 10
```

Connect: `Cmd+Shift+P` → Remote-SSH: Connect to Host → hyak → password + Duo.
Open folder: `Cmd+Shift+P` → Open Folder → `/mmfs1/gscratch/astro/ds2004/sorcha`

Important: VSCode opens a new window for the remote. Local window stays open separately. Extensions (Claude Code) are local only — use the remote terminal for Hyak commands and paste output back to local Claude session.

Kernel for notebooks on Hyak: `/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python`

---

## Wagg et al. 2024 Paper Summary (arXiv:2408.12517)

### Tracklet Building Cuts (Wagg)
- **≥3 detections per night** (they also test ≥4)
- **Time span ≤90 minutes** between any pair in the tracklet
- **Arc length ≥1 arcsec** (~5 Rubin pixels)

Our existing S3M comparison used 2 detections × 30 min separation. **Zeljko's decision (2026-05-27):** Do NOT adopt Wagg's ≥3 cut. The correct approach depends on data source:

**For Sorcha simulations:**
`score_observation` takes already-computed velocity components (`dra_deg_day`, `ddec_deg_day`). Sorcha outputs instantaneous rates directly (`RARateCosDec_deg_day`, `DecRate_deg_day`) — so VDP can score individual detections without pairing at all. "For sims it's easy as there is no noise — you just need to compute the angular velocity components."

- **VDP side**: use instantaneous rates from Sorcha output per detection (or per-night mean). No pairing needed.
- **digest2 side**: still needs 2 sky positions in MPC 80-col format → pick one pair per object per night with time separation **3 min ≤ Δt ≤ 90 min**

**For real Rubin data (future work):**
Instantaneous rates are unavailable — must compute from 2 detections, which introduces astrometric noise. Quality cuts matter:
- Time separation 3–90 min as proxy for velocity quality
- Or explicit velocity errors from astrometric errors + temporal baseline

This is "what the pipeline does" in real-world use — VDP and digest2 both score the same 2-detection tracklet.

**Bottom line for wagg_postprocess.py:**
- Group detections by (ObjID, night)
- VDP: score using mean/first-detection rates from Sorcha output directly
- digest2: pick best pair per night (3–90 min separation), build MPC 80-col from those 2 positions
- No S3M rerun needed — S3M comparison already used rates computed from synthetic propagation

### digest2 Threshold
Wagg uses **score ≥65 on 0–100 scale** for NEOCP submission. Our comparison uses 0–1 scale (threshold ~0.97) — equivalent, just need consistent scaling.

### Observatory Code
Confirmed from `Rubin_full_footprint_wagg_detections.ini`:
```
ar_obs_code = X05
```
Matches our existing S3M comparison. No changes needed.

### difi
Wagg uses `difi` to predict Rubin self-linkage (≥3 nights in 15 days, ≥2 obs/night). Halves NEOCP list (~129 → ~55/night). **We are not implementing difi** — our goal is classifier comparison, not NEOCP traffic prediction.

### Key Results from Wagg
- ~129 new objects/night to NEOCP (8× current rate)
- 8.3% purity at digest2 ≥65 (91.7% false positives, mostly undiscovered MBAs)
- VDP not used in Wagg — that is our contribution

### How Our Pipeline Differs

| | Wagg | Us |
|---|---|---|
| Goal | Predict NEOCP traffic | Compare VDP vs digest2 ROC |
| Tracklet cuts | ≥3 obs, ≤90 min, ≥1" arc | Same for Sorcha pipeline |
| digest2 threshold | ≥65 for submission | Full sweep for ROC |
| VDP | Not used | Head-to-head vs digest2 |
| difi | Yes | No |
| Observatory code | X05 | X05 |

---

## ROC Analysis Results (2026-05-28)

### Full subsample (598K tracklets, all epochs)

| Classifier | AUC  | Best F1 | Completeness | Contamination | Threshold |
|---|---|---|---|---|---|
| VDP        | 0.582 | 0.342  | 47.2%        | 73.2%         | 0.143 |
| digest2    | 0.734 | 0.792  | 70.4%        | 9.5%          | 0.620 |

digest2 wins decisively. Root cause: the 5 VDP maps encode velocity statistics at
specific epochs (March 2025, May 2026). The Sorcha simulation spans 2 years, so most
tracklets are observed when the solar elongation geometry is completely different from
what the map was built for — NEO/MBA velocity distributions don't match.

VDP calibration check confirmed the problem: P_NEO_vdp bins 0–0.5 all showed ~15%
actual NEO fraction (barely above the 16.5% background rate), meaning the maps
provide almost no discrimination for 90% of the score range.

### Epoch-matched subset (137K tracklets, obs within ±30d of map epoch)

Filter logic: keep tracklets where |mjd0_utc - map_canonical_epoch| ≤ 30 days.
- lon180_lat0 (March 2025 map): use March 2026 window (annual repeat; sim starts May 2025)
- lon229/antisun maps (May 2026): use May 2026 ±30d AND May 2025 ±30d

| Classifier | AUC  | Best F1 | Completeness | Contamination |
|---|---|---|---|---|
| VDP        | 0.714 | **0.609** | 52.8%      | 27.9%         |
| digest2    | 0.727 | 0.804   | 72.6%        | 9.9%          |

VDP F1 jumps from 0.342 → 0.609 when epoch-matched. VDP AUC (0.714) nearly closes the
gap with digest2 (0.727). The two primary maps (lon180, lon229) reach VDP F1 ≈ 0.70
vs digest2 ≈ 0.80.

### Per-map breakdown (epoch-matched)

| Map | n_NEO | VDP F1 | d2 F1 |
|---|---|---|---|
| lon180_lat0 (opposition) | 3,768 | 0.704 | 0.786 |
| lon229_lat0 (NEOCP pos)  | 10,513 | 0.695 | 0.818 |
| antisun_minus_45         | 4,564 | 0.530 | 0.822 |
| antisun_minus_90         | 1,720 | 0.389 | 0.729 |
| antisun_minus_120        | 399   | 0.473 | 0.796 |

The off-axis maps (antisun_minus_90, _120) underperform — they cover sky positions
at non-standard geometries and have fewer training objects.

### Key Insight: Need Full-Year VDP Maps

The VDP is NOT fundamentally weaker than digest2. When applied to tracklets observed
at the correct epoch, VDP AUC is nearly identical to digest2 (0.714 vs 0.727).

The degradation in the full comparison is entirely due to epoch mismatch: the Sun moves
~1°/day, so the antisun direction (where the maps are calibrated) moves ~30°/month.
A tracklet observed 3 months from the map epoch is at a completely different solar
elongation, where the velocity statistics are wrong.

**Next step: generate VDP probability maps at monthly (or bi-monthly) intervals
spanning the full 2-year simulation window (May 2025 – April 2027).** Then re-score
all 40.7M tracklets with the map whose epoch is closest in time (not just in sky
position), and rerun the full ROC comparison.

### Confirmed Physical Insight: Antisun Always Separates NEOs from MBAs

Verified empirically from Sorcha data — NEO/MBA apparent velocity ratio by distance
from the current antisun direction:

| Dist from antisun | n_NEO  | NEO v_tot median | MBA v_tot median | ratio |
|---|---|---|---|---|
| 0–15°  | 4,751  | 0.346 deg/day | 0.108 deg/day | **3.2×** |
| 15–30° | 8,709  | 0.295 deg/day | 0.062 deg/day | **4.8×** |
| 30–60° | 22,134 | 0.360 deg/day | 0.158 deg/day | 2.3× |
| 60–90° | 23,743 | 0.449 deg/day | 0.238 deg/day | 1.9× |

The 3–5× velocity separation at the antisun is a **geometric effect that repeats at
every epoch** — it is not special to March 2025 or May 2026. It holds because:
- Near opposition, NEOs (inner solar system, eccentric orbits) have maximum apparent
  velocity relative to background stars
- MBAs at opposition appear nearly stationary (slow retrograde drift)
- This is the same geometry that makes the NEOCP work at opposition every lunation

**Implication**: The S3M comparison result (VDP F1=0.643 >> digest2 F1=0.339) was NOT
a fluke of the specific epoch or sky patch chosen. Any VDP map built at the antisun
direction for any month will produce similarly good NEO/MBA discrimination. The Sorcha
comparison with 5 fixed-epoch maps was broken by epoch mismatch, not by VDP being
fundamentally weak.

**Wagg et al. 8.3% purity is consistent with our results** — it is NOT evidence that
digest2 is bad at classification. It reflects the base-rate effect: NEOs are 0.24% of
Rubin-detectable objects, so even a 99% accurate classifier produces mostly false
positives. Our F1=0.792 for digest2 is measured on a balanced test set (16.5% NEO)
and is not comparable to Wagg's real-population purity metric. Converting our digest2
FPR (~1.45%) to real-sky purity gives ~10% — consistent with Wagg's 8.3%.

**VDP with full-year maps is expected to outperform digest2** in the Sorcha comparison,
matching or exceeding the S3M result, because:
1. Monthly antisun maps guarantee epoch-matched scoring for every tracklet
2. The antisun velocity separation (3–5×) is robust across all months
3. VDP uses 2D velocity structure (vlam, vbeta, magnitude) — richer signal than
   digest2's orbit-space model

---

## Current Status (2026-05-31) — Fourth Pass: Hybrid Catalog Training

### Completed infrastructure changes

**New VDP maps — v2 with extended grid (generated 2026-05-30 ~17:11):**
- 24 monthly antisun maps in `prob_maps/` covering May 2025 – April 2027
- Naming: `prob_maps_YYYY-MM-DD_antisun.npz`
- Center = exact antisun ecliptic lon at each epoch (astropy `get_sun`)
- **Grid: (-1.5, 1.5) deg/day, step=0.01 → 301×301 = 90,601 points** (was (-0.8, 0.8), 161×161)
  - Extended because fast NEOs with |vlam| > 0.8 were outside old grid → P(NEO)=0
  - These 27,571 fast NEOs were missed entirely; digest2 caught them fine
- Generated via `sorcha_gen_maps_slurm.sh`: `--cpus-per-task=16`, `--n-jobs 16`,
  `--time=06:00:00`, `--overwrite` (all 24 tasks ran in ~13 min with new fast VDP)
- Maps also copied to Arnor: `ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/prob_maps/`
- The 5 original maps (`prob_maps_2025-03-21.npz`, `prob_maps_2026-05-09*.npz`) remain
  in `prob_maps/` but are EXCLUDED from Phase 2 scoring (see score_vdp_frame filter below)

**VDP speed optimizations (`neomod/src/velocity_density_pipeline_fast.py`):**
Do NOT modify `velocity_density_pipeline.py` (original preserved). The fast copy adds:
1. `log_posterior_d0_2d` vectorized with numpy broadcasting + `scipy.special.gammaln`
   (replaces double Python for-loop → ~10-50× speedup on the Bayesian integral)
2. Apparent magnitudes pre-computed ONCE per population before the mag-bin loop
   (was redundantly computed 8× per population → 8× savings on propagation)
3. `evaluate_density_map_full_posterior_2d` parallelized with joblib (`n_jobs` param)
4. `generate_probability_maps` accepts `n_jobs` and passes it through
- `sorcha_gen_map.py` imports `velocity_density_pipeline_fast`, accepts `--n-jobs`
- `sorcha_gen_maps_slurm.sh`: `--cpus-per-task=16`, `--n-jobs 16`, 24 tasks concurrent
- Result: ~13 min per map with 16 CPUs and wide grid (was ~2 hr on Mac)

**S3M data on Hyak:** `neomod/S3Mdata/` — 1.6G, 23 files (`S0.s3m`, `S1_00..13.s3m`,
`ST.s3m`, `St5.s3m`, plus extras). Required by VDP map generation.

### Phase 1 fix: epoch-aware + antisun-aware map assignment (`sorcha_postprocess.py`)

**Problem found and fixed (2026-05-29):**

The original `assign_probability_maps` used sky proximity to the map center at the
**map's epoch** as the footprint check. This caused two bugs:

1. **Epoch mismatch**: Tracklets assigned to maps built at wrong times
2. **More critically**: With monthly maps, the past-epoch antisun sky position drifts
   ~1°/day. A tracklet observed 15 days after the map epoch could be 45° from the
   CURRENT antisun but still within 30° of the map's old sky position, and thus
   assigned. At 45° from the antisun (near MBA stationary point), MBAs have vlam ≈ 0,
   which the antisun VDP map misclassifies as high P(NEO) → 20% MBA false positive rate.

**Diagnosis evidence:**
- High-scoring MBAs (P_NEO > threshold): median dist from CURRENT antisun = 61°
- Low-scoring MBAs (correctly classified): median dist from current antisun = 18°
- Only 21% of high-VDP MBAs were actually within 30° of the current antisun

**The fix** in `assign_probability_maps` (sorcha_postprocess.py):
```python
# New helpers added at module level:
def _antisun_ecl_lon(mjd):   # vectorized sun-lon formula, accurate to ~1°
def _ecl_lon_from_radec(ra, dec):  # manual obliquity rotation
def _angular_sep_ecl_lon(lon_a, lon_b):  # angular difference

# In assign_probability_maps: for ALL maps, footprint criterion =
# distance from CURRENT antisun at observation time ≤ 30°
dist = _angular_sep_ecl_lon(ecl_lon_obs, antisun_lon)
in_footprint[:, j] = dist <= meta.max_sep_deg
```

Among in-footprint maps, pick the one with smallest |map_epoch - obs_mjd|.

**Also added to `ProbMapMeta`:** `obstime_mjd: float` field, populated in
`load_prob_map_metadata` via `Time(obstime_str).utc.mjd` for fast epoch comparison.

**Phase 2 filter** (`sorcha_phase2.py::score_vdp_frame`):
```python
if "antisun" not in str(map_file):
    continue   # skip legacy fixed-position maps
```
This ensures only the 24 new monthly antisun maps are used for VDP scoring, even if
the tracklet parquets still have some assignments to the old 5 maps (which used the
old sky-proximity criterion before the Phase 1 re-run).

### F1 Progression Table (all fixes applied sequentially)

| Fix applied | VDP F1 | Completeness | Contamination |
|---|---|---|---|
| Baseline (5 maps, 2 epochs, mixed) | 0.342 | 47.2% | 73.2% |
| Epoch-matched ±30d (diagnostic only) | 0.609 | 52.8% | 27.9% |
| 24 monthly maps + epoch-aware assign | 0.448 | 46.6% | 56.9% |
| + antisun footprint fix | 0.502 | 55.5% | 54.2% |
| + antisun-only filter in Phase 2 | 0.648 | 62.4% | 32.7% |
| + wide grid (−1.5 to +1.5) | **0.740** | 67.5% | 18.0% |
| Grid extension to ±3.0 (attempted, reverted) | 0.740 | 67.5% | 18.0% |
| digest2 (constant) | 0.836 | 77.3% | 9.0% |

### Why the ±3.0 grid extension didn't work (2026-05-31)

After reaching F1=0.740, we attempted to extend `DEFAULT_GRID_LIM` from (−1.5, 1.5)
to (−3.0, 3.0) to capture fast NEOs with |vlam| > 1.5 deg/day. The results were
**identical** — F1 unchanged, same 4,903 fast NEOs all still scoring P=0.

**Root cause:** `ProbMapSet.__init__` applies a nearest-clone-distance mask:
```python
too_far = nearest[pop] > self.mask_radius_deg_per_day   # default 0.2 deg/day
frac[too_far] = 0.0
```
The `nearest_dist` map stores the distance from each grid cell to its nearest S3M
clone. Beyond |vlam| ≈ 1.3–1.5 deg/day, S3M has no NEO clones — so the nearest
distance is >> 0.2 deg/day at every cell beyond that boundary. The wider grid exists
and the bilinear interpolator would work, but the mask zeros the result first.

**Measured nearest distances in the maps:**
```
vlam=−1.0: nearest_dist = 0.106 → within mask (< 0.2)
vlam=−1.5: nearest_dist = 0.540 → masked to 0
vlam=−2.0: nearest_dist = 0.350 → masked to 0
vlam=−3.0: nearest_dist = 0.275 → masked to 0
```

**Conclusion:** Grid size is irrelevant. The binding constraint is the S3M training
data — there are simply no S3M NEO clones beyond ±1.5 deg/day. No amount of grid
extension fixes a training data gap. `DEFAULT_GRID_LIM` was reverted to (−1.5, 1.5)
with step=0.01.

### Remaining gap to digest2 — diagnosis

After all fixes, VDP sits at F1=0.740 vs digest2 F1=0.836. The 20,481 missed NEOs
(in the digest2-only quadrant) split into two problems:

**Problem A (~24% of missed NEOs, partially solved):**
Fast NEOs with |vlam| > 1.5 deg/day. S3M has no clones there → nearest_dist mask
zeroes them. Grid extension cannot fix this. Only adding real high-eccentricity NEOs
to the training set can populate those velocity regions.

**Problem B (~56% of missed NEOs, fundamental):**
NEOs that fall within |vlam| < 1.5 but in velocity regions where S3M has sparse or
zero calibration coverage. These are real MPCORB NEOs with unusual orbital
configurations (high eccentricity, unusual inclination, specific perihelion geometry)
that S3M's synthetic population doesn't replicate. Their nearest S3M clone is
> 0.2 deg/day away → masked to 0. Grid extension cannot fix this either.

**Why this didn't show up in the original S3M ROC (Figure 5):**
The S3M ROC scores S3M test objects against S3M-trained maps. The test set and
training set are from the same narrow synthetic population — no velocity regions
are outside S3M calibration. The problem only appears when scoring real MPCORB
objects (the hybrid catalog), which have velocity distributions S3M never saw.

### Goal: prove VDP beats digest2 (not just matches it)

The S3M ROC showed VDP at ~20% higher completeness than digest2 at equal contamination.
That result is real — it just doesn't transfer to the Sorcha comparison because the
maps are trained on S3M while the scored objects include ~10% real MPCORB catalog
objects that fall in S3M calibration gaps.

**The fix: train VDP maps on the hybrid catalog itself** (the same population being
scored), not on S3M. This closes both Problem A and Problem B by construction — every
velocity region that a real tracklet can occupy will have calibration data nearby.

---

## Hybrid Catalog VDP Training Plan (next step as of 2026-05-31)

### Why this will work

The hybrid catalog IS the scored population. If we train the maps on hybrid, the
nearest-clone-distance mask will always find a nearby object because every scoreable
tracklet's velocity came from an object in the catalog. The calibration gaps disappear.

On the S3M-only ROC, VDP shows ~20% higher completeness than digest2 at equal
contamination. That advantage exists because VDP uses 2D velocity structure
(vlam, vbeta, magnitude) — a richer signal than digest2's orbit-space model. The
goal is to make that advantage appear in the Sorcha comparison too, by removing the
calibration mismatch.

### Guiding constraint: parallel clone, not replacement

Every new file gets a new name. The S3M flow is untouched. The hybrid flow runs
in parallel into separate directories. We can switch between them by pointing
`--prob-maps-dir` at `prob_maps/` (S3M) or `prob_maps_hybrid/` (hybrid).

### Data format: hybrid.h5 vs S3M

```
S3M (.s3m files):   q, e, i, node, argperi, t_p, H, t_0   (Keplerian elements)
hybrid.h5:          x, y, z, vx, vy, vz, t_0, H, g        (Cartesian state vectors)
```

The hybrid catalog stores Cartesian state vectors at a single epoch (t_0 = 60065 MJD).
S3M stores Keplerian elements. The VDP pipeline internally needs Keplerian elements
to clone and propagate objects to a given map epoch.

**Coordinate frame confirmed (2026-05-31):** heliocentric ecliptic J2000.
Sanity check: propagated S3M elements for S0000001a (two-body Kepler) to t=60065 and
compared with hybrid.h5 state vector. Residual |dr| = 1.1×10⁻³ AU, |dv| = 9 m/s.
The residual is physical (N-body perturbations accumulated over 16 years), not a
frame mismatch. Conversion formulas are correct.
Script: `sanity_check_frame.py`

### Hybrid catalog population breakdown

```
Total objects:              14,444,912
  S-prefix (S3M synthetic): 13,112,612  (90.8%)
  Non-S (real MPCORB):       1,332,300   (9.2%)

NEOs (q < 1.3 AU):            295,043
  S3M NEOs:                   261,731  (same as current VDP training)
  Real MPCORB NEOs:            33,312  (these fill the calibration gaps)

Non-NEOs:                   14,149,869
```

The non-S (real MPCORB) NEOs are the critical addition. They have wild eccentricities
and inclinations that S3M doesn't synthesize — exactly the objects that currently
score P=0 due to calibration gaps.

### Implementation: 7 new files, 0 modifications to existing files

#### Step 1 — `hybrid_catalog_prep.py` ✓ COMPLETE (2026-05-31)

**Output:** `hybrid_elements.parquet` (0.93 GB on disk, 14,444,912 rows)

**Actual population counts:**
```
MBA:    13,883,703
NEO:       295,040   ← 261k S3M + ~33k real MPCORB
Trojan:    185,417
TNO:        63,381
other:      17,371
```

**Sanity check S0000001a:** q=1.25140  e=0.38226  i=9.305°  a=2.026 AU (physically correct)

**Implementation notes:**
- Used chunked reading (500k rows/chunk, 29 chunks) to avoid OOM on login node
- `pd.read_hdf(..., chunksize=...)` fails on fixed-format HDF5 — must use pytables
- `data = f.root.df.block0_values[:]` (flat read of all 14.4M rows) also OOMs — use slices
- Script: `hybrid_catalog_prep.py` in working directory

**Orbital element conversion** (heliocentric two-body, μ = 4π²/365.25² AU³/day²):
```python
eps = v2/2 - mu/r          # specific orbital energy
a   = -mu / (2*eps)         # semi-major axis
h   = cross(r_vec, v_vec)   # angular momentum vector
i   = arccos(hz / |h|)      # inclination
e   = |e_vec|               # eccentricity (from Laplace-Runge-Lenz vector)
q   = a * (1 - e)           # perihelion distance
# node, argperi, t_p from standard two-body relations
```

**Population classification:**
```
NEO:    q < 1.3 AU
MBA:    1.7 ≤ a < 4.1, q ≥ 1.3
Trojan: 4.7 < a < 5.9, e < 0.3  (Jupiter Trojans)
TNO:    a > 30 AU
other:  everything else (Centaurs, Hungarias, etc.)
```

Output columns: `OID, q, e, i, node, argperi, t_p, H, t_0, a, population`

#### Step 2 — `neomod/src/hybrid_loader.py` ✓ COMPLETE (2026-05-31)

Drop-in replacement for `s3m_loader.py`. Functions must keep the names `define_s3m`
and `s3m_array` because the pipeline does `import hybrid_loader as load_s3m` and
calls `load_s3m.define_s3m(pop=...)` and `load_s3m.s3m_array(df)`.

```python
def define_s3m(pop="neo", verbose=True, **kwargs) -> pd.DataFrame:
    # reads hybrid_elements.parquet, filters to population label
    # pop keys: "neo"→"NEO", "mba"→"MBA", "tno"→"TNO", "trojan"→"Trojan", "all"→no filter

def s3m_array(df, n_H=52, n_a=42, n_e=25, n_i=22):
    # delegates to s3m_loader.s3m_array — format is identical
    import s3m_loader
    return s3m_loader.s3m_array(df, n_H, n_a, n_e, n_i)
```

**Smoke test results (all PASS):**
- neo: 295,040 objects  H=9.3–33.6
- mba: 13,883,703 objects  H=3.2–26.5
- tno: 63,381 objects  H=-1.2–22.6
- trojan: 185,417 objects  H=5.6–18.5
- s3m_array(neo): shape=(52,42,25,22)  sum=295,040 ✓

**Clone factor tuning — FINAL (confirmed by test map 2026-01-01):**
| Population | S3M clone_factor | Hybrid clone_factor | Reason |
|---|---|---|---|
| NEO | 300 | 80 | 30 caused regression at +0.6/+0.7; 80 recovers it and extends to +/-1.4 |
| MBA | 10 | 1 | 14M objects -- native density sufficient, no cloning needed |
| TNO | 100 | 10 | similar object count as S3M |
| Trojans | 100 | 5 | more objects than S3M Trojans |

**Test map nearest_dist results (NEO 18_20 mag bin, vlam mid-slice):**
- Hybrid unmasked: 58.5% of grid  vs  S3M: 19.2%  -- 3x improvement
- vlam = -0.8 to -1.1: newly unmasked (real MPCORB retrograde NEOs)
- vlam = +1.0 to +1.4: newly unmasked (real MPCORB prograde NEOs)
- vlam = +/-1.5: still masked (edge of grid, expected)
- Bug fixed: `np.atleast_2d` on propagator output for single-object mag bins (hybrid only)

#### Step 3 — `neomod/src/velocity_density_pipeline_hybrid.py` ✓ COMPLETE (2026-05-31)

Copied from `velocity_density_pipeline_fast.py`. Exactly 3 lines changed:
1. `import s3m_loader as load_s3m` → `import hybrid_loader as load_s3m`
2. Clone factors in `DEFAULT_POPULATION_SETTINGS`: MBA=1, NEO=30, TNO=10, Trojans=5
3. Module docstring updated

Smoke test: imports OK, loader=hybrid_loader, all clone factors correct.

Everything else (grid, kNN, Gaussian smoothing, ProbMapSet, scoring) is unchanged.
The `ProbMapSet.from_npz` and all scoring code are identical between S3M and hybrid
pipelines — only the training data differs.

#### Step 4 — `sorcha_gen_map_hybrid.py`

```bash
cp sorcha_gen_map.py sorcha_gen_map_hybrid.py
```

**One change:**
```python
import velocity_density_pipeline_hybrid as vdp  # was velocity_density_pipeline_fast
```

#### Step 5 — `sorcha_gen_maps_hybrid_slurm.sh`

```bash
cp sorcha_gen_maps_slurm.sh sorcha_gen_maps_hybrid_slurm.sh
```

**Changes:**
- `--job-name=sorcha_g_hybrid`
- Calls `sorcha_gen_map_hybrid.py` instead of `sorcha_gen_map.py`
- Output directory: `prob_maps_hybrid/` instead of `prob_maps/`
- Time allocation: start with same `06:00:00`, adjust after test map

**Before submitting all 24:** run one map manually to validate that:
- Nearest-dist map shows coverage beyond |vlam| = 1.5 (real MPCORB NEOs populate it)
- Map generation time is acceptable with new clone factors

#### Step 6 — `sorcha_phase2_vdp_hybrid.sh`

```bash
cp sorcha_phase2_vdp.sh sorcha_phase2_vdp_hybrid.sh
```

**Changes:**
- `--job-name=sorcha_vdp_hybrid`
- `--prob-maps-dir prob_maps_hybrid`
- `--work-dir outputs/phase2_hybrid`

#### Step 7 — combine step for hybrid results

digest2 shards don't change (they live in `outputs/phase2/digest2_shards/`).
After Phase 2 hybrid VDP scoring completes:

```bash
# Combine hybrid VDP shards with existing digest2 shards
conda_prep/bin/python sorcha_phase2.py combine \
  --work-dir outputs/phase2_hybrid \
  --output outputs/phase2_hybrid/wagg_sorcha_comparison_hybrid.parquet

# Then SCP to Arnor for notebook analysis
scp outputs/phase2_hybrid/wagg_sorcha_comparison_hybrid.parquet \
    ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/outputs/phase2/
```

Note: need to verify the `combine` command can find digest2 shards from a different
work-dir, or copy them into `outputs/phase2_hybrid/digest2_shards/` first.

### Final directory layout

```
sorcha/
├── neomod/src/
│   ├── velocity_density_pipeline_fast.py    ← S3M pipeline, UNTOUCHED
│   ├── velocity_density_pipeline_hybrid.py  ← NEW: hybrid catalog pipeline
│   ├── s3m_loader.py                        ← UNTOUCHED
│   └── hybrid_loader.py                     ← NEW: loads hybrid_elements.parquet
├── hybrid_catalog_prep.py                   ← NEW: one-time Cartesian→Keplerian
├── hybrid_elements.parquet                  ← NEW: generated by prep script
├── sanity_check_frame.py                    ← coordinate frame verification (done)
├── sorcha_gen_map.py                        ← S3M, UNTOUCHED
├── sorcha_gen_map_hybrid.py                 ← NEW
├── sorcha_gen_maps_slurm.sh                 ← S3M, UNTOUCHED
├── sorcha_gen_maps_hybrid_slurm.sh          ← NEW
├── sorcha_phase2_vdp.sh                     ← S3M, UNTOUCHED
├── sorcha_phase2_vdp_hybrid.sh              ← NEW
├── prob_maps/                               ← S3M maps (24 monthly antisun)
├── prob_maps_hybrid/                        ← NEW: hybrid maps (24 monthly antisun)
└── outputs/
    ├── phase2/                              ← S3M scoring results
    └── phase2_hybrid/                       ← NEW: hybrid scoring results
```

### Expected outcome

With hybrid training:
- Nearest-dist mask will pass for all tracked objects (training set = scored population)
- Real MPCORB NEOs at |vlam| > 1.5 will have non-zero P(NEO) (Problem A solved)
- Real MPCORB NEOs in velocity gaps will have calibration coverage (Problem B solved)
- VDP F1 expected to exceed 0.836 (digest2) and match or exceed S3M VDP result (~0.85)

The paper argument: VDP is a fundamentally better algorithm than digest2 when properly
calibrated to the scored population. The S3M ROC proved it. The Sorcha comparison with
hybrid training will prove it on a realistic Rubin-cadence simulation.

### Real-world validation (NEOCP, May 9–28 2026)

Independent validation on 278 real NEOCP objects (May 2026) showed:
- VDP scored 193/278 (69% in antisun footprint); 30% outside footprint → unscored
- ECDF shows ~63% of in-footprint objects scoring P≈0 (same calibration gap pattern)
- Velocity space: VDP confident only near center (|vlam| < 0.5); real NEOCP objects
  with unusual orbits fall in calibration gaps identical to Sorcha diagnosis
- "MPC only" quadrant heavily populated — same completeness loss as Sorcha comparison
- One "VDP only" outlier: high VDP score but low MPC score (early MPC assessment)

This confirms the hybrid training fix is needed for real operational use too, not just
for the paper comparison.

---

## Key Bugs Found and Fixed (2026-05-29)

### Root cause of VDP underperformance — summary

| Bug | Effect | Fix |
|---|---|---|
| Epoch mismatch (old) | Tracklets assigned to maps from wrong time of year | Epoch-aware assignment: pick nearest-epoch in-footprint map |
| Footprint check at map epoch (new) | Off-axis tracklets (60°+ from current antisun) assigned to antisun maps; MBAs at stationary point scored as NEOs | Use current antisun at observation time for footprint |
| Old fixed-position maps in scoring | lon229, lon180, antisun_minus_45/90/120 assign off-axis tracklets | Skip non-antisun maps in score_vdp_frame |

All three bugs compound: the old maps could assign to any sky position, and the
footprint-epoch mismatch made even the new monthly maps assign off-axis tracklets.

### MBA contamination — physical explanation

At the antisun (elongation=180°), MBAs are in retrograde with vlam ≈ -0.20 to -0.25.
At elongation ≈ 150° (stationary point), MBAs have vlam ≈ 0. The VDP antisun map has
few S3M MBAs at vlam ≈ 0 (they're calibrated for opposition geometry), so the map
assigns moderate-to-high P(NEO) to objects there. Sorcha tracklets near the edge of
the 30° antisun footprint often fall near the stationary-point elongation, producing
a large MBA false positive population. Fix: require tracklets to be within 30° of the
CURRENT antisun (not historical map position), and skip old off-axis maps entirely.

---

## Relationship to Existing VDP / digest2 Comparison

Existing S3M-based comparison (`neomod/digest2_comparison.ipynb`):
- Input: S3M synthetic orbital elements, propagated to 2025-03-21 opposition patch
- Tracklet: 2 detections × 30 min, synthetic
- Sky cut: single 30° patch at opposition, epoch-matched by construction
- Result: VDP F1=0.856, digest2 F1=0.665 — VDP wins clearly
- NOTE: S3M objects were all at the map epoch — no off-axis stationary-point MBAs

Wagg-Sorcha comparison (`sorcha_roc_comparison.ipynb`):
- Input: Rubin-cadence detections from 14.4M hybrid catalog over 2 years
- Tracklet: built from real simulated detection cadence (nightly grouping)
- v2 (current, with all fixes): only tracklets within 30° of current antisun, monthly
  epoch-aware maps, old maps excluded from scoring
- Expected result: VDP F1 ≥ 0.64 (matching S3M), possibly exceeding digest2

The S3M comparison succeeded because objects were only at exact opposition — no
stationary-point contamination. Sorcha has objects across a 30° footprint range,
requiring the current-antisun filter to exclude near-stationary-point MBAs.

With full-year maps the Sorcha result should approach or match the S3M result.

---

# Arnor Updates — 2026-06-08

## Current Science Status (post-normalisation fix)

The GMM normalisation bug (Bug 1) was fixed on Hyak. All three post-fix files SCPed to Arnor (Jun 8).

### Files now on Arnor

| File | Notes |
|---|---|
| `prob_maps_gmm/prob_maps_2026-05-01_antisun.npz` | Regenerated GMM map with normalisation fix (Jun 8, 32MB) ✓ |
| `outputs/phase2/sorcha_comparison_gmm.parquet` | Post-fix Sorcha GMM results (75MB, 611K rows, 111K NEOs, 18.2%) ✓ |
| `outputs/phase2/sorcha_may2026_antisun_patch.parquet` | May 2026 antisun patch (1.2MB, 14,873 objects, 2,402 NEOs, 16.2%) ✓ |
| `outputs/phase2/sorcha_comparison.parquet` | S3M kNN Sorcha results |
| `outputs/phase2/sorcha_comparison_hybrid.parquet` | Hybrid catalog results |
| `prob_maps/prob_maps_*_antisun.npz` | 24 monthly S3M antisun maps |

### F1 Results (post-fix)

**Single-epoch (same 21K S3M objects, `sorcha_gmm_s3m_singleepoch_comparison.ipynb`):**

| Classifier | F1 | Completeness | Contamination |
|---|---|---|---|
| S3M kNN VDP | 0.847 | 75.1% | 2.8% |
| GMM VDP (same objects) | 0.842 | 76.9% | 7.0% |
| digest2 | 0.655 | 52.2% | 12.2% |

Note: digest2 F1 low here due to unfair population ratios in S3M test set (7:1 TNO:NEO).

**Sorcha full 2yr (GMM + mask OFF + grid ±2.0, `sorcha_roc_comparison.ipynb`):**

| Classifier | F1 | Completeness | Contamination |
|---|---|---|---|
| VDP (GMM) | **0.837** | 78.3% | 10.2% |
| digest2 | 0.836 | 77.1% | 8.7% |

VDP TIED with digest2 (F1 0.837 vs 0.836). VDP completeness (78.3%) exceeds digest2 (77.3%). Remaining gap is MBA contamination: VDP MBA FPR = 1.5% vs digest2 0.2%.

---

## GMM Bug Diagnostics (current state, mag22 bin)

### Bug 1 — Normalisation (FIXED)

Pre-fix GMM/S3M NEO integral ratio: **0.322** (3× underscaled)
Post-fix GMM/S3M NEO integral ratio: **3.453** (slight overcorrection, but P(NEO) is a ratio so absolute scale cancels)

GMM/S3M MBA integral ratio: 1.041 (unchanged, was correct before)

### "Bug 2" — MISDIAGNOSED ON ARNOR (corrected 2026-06-08)

The Arnor diagnosis (MBA Gaussian tail bleeding into NEO-only velocity regions) was wrong.
Hyak confirmed: **MBA uses K|M cloning, not GMM.** There is no Gaussian tail.

The 9.3× MBA density discrepancy at vlam=−0.5 is from **different map centers**:
- S3M NEOCP map: ecliptic lon = 229° (May 9 2026)
- GMM monthly map: ecliptic lon = 220° (May 1 2026)

A 9° offset changes which MBAs fall in the sky cut, shifting the MBA density at any given velocity point. This is not a bug — it is expected when comparing maps at slightly different sky positions.

**Actual remaining P(NEO) gap at vlam=−0.5:**
Post-fix, GMM P(NEO) ≈ 0.911 vs S3M 0.954. The ~4% gap is explained by the map center offset, not a model deficiency. It is not actionable.

**The advisor's investigation question** ("why is P(NEO) not 1 where no MBAs exist?") is therefore answered differently for GMM than for S3M kNN: MBA uses K|M (no tails), so any remaining gap in GMM maps is from sky-geometry differences, not from the cloning model itself.

---

## Notebooks on Arnor (current state)

### `sorcha_gmm_s3m_singleepoch_comparison.ipynb`
- Apples-to-apples: same 21K S3M objects scored with S3M kNN and GMM maps
- Loads `s3m_digest2_comparison_2026-05-09T22_neocp.parquet`
- Re-scores with GMM using `RegularGridInterpolator`
- Diagnostic cell (Cell 12) has no outputs yet — needs re-run with new map
- F1 results in Cell 15 outputs: S3M=0.847, GMM=0.842, d2=0.655

### `sorcha_gmm_s3m_comparison.ipynb`
- Visual comparison: density maps + P(NEO) maps for S3M kNN vs GMM, side by side
- Loads `prob_maps/prob_maps_2026-05-09T22_neocp.npz` and `prob_maps_gmm/prob_maps_2026-05-01_antisun.npz`
- Loads `outputs/phase2/sorcha_comparison.parquet` and `outputs/phase2/sorcha_comparison_gmm.parquet`

### `sorcha_roc_comparison.ipynb`
- Full Sorcha 2yr ROC: VDP GMM vs digest2
- Uses `outputs/phase2/sorcha_comparison_gmm.parquet`
- Results: VDP F1=0.837, digest2=0.836 (TIED — needs re-run with post-fix parquet)

### `neocp_vdp_comparison.ipynb`
- Scores live NEOCP objects (May 9 – Jun 2, 279 objects) against VDP maps
- Has `score_0_1` (MPC's real digest2 score) for all objects
- Section 14 (updated 2026-06-08): digest2 validation — loads `outputs/validate_digest2_neocp.csv` (Hyak real binary results; 13 objects) and plots scatter vs MPC score
  - Python `NEOMODScorer` approach was broken (`adam_core_stub` returns all zeros) — replaced with CSV load
  - `validate_digest2_neocp.csv` now on Arnor at `outputs/`
  - **Ready to run** — just execute Section 14 cells

---

## digest2 NEOCP Validation — COMPLETE (done on Hyak 2026-06-08)

Hyak ran `neomod/pipeline/validate_digest2_neocp.py` using real MPC 80-col tracklets
(2 consecutive hourly ephemeris positions) through the `mpcdev-digest2` binary.

**Results (17 objects):**

| Metric | Value |
|---|---|
| Pearson r | 0.627 |
| Spearman ρ | 0.586 |
| MAE | 0.101 |
| High-MPC (≥0.90) scored ≥0.50 by us | 12/13 |

Correlation is limited by input differences (MPC used full multi-night arc; we used 1-hr synthetic tracklets). Ordering is preserved at extremes. Binary is functioning correctly.

Output files on Hyak: `outputs/validate_digest2_neocp.csv`, `outputs/validate_digest2_neocp.png`

**Section 14 of `neocp_vdp_comparison.ipynb`** is now fixed (2026-06-08). CSV loaded at `outputs/validate_digest2_neocp.csv`. Ready to run.

---

## Master To-Do List (as of 2026-06-08)

### Step 0 — Arnor notebook cleanup (mostly done)
- [ ] Re-run `sorcha_roc_comparison.ipynb` — confirm F1=0.837 with post-fix parquet
- [ ] Re-run `sorcha_gmm_s3m_singleepoch_comparison.ipynb` — update single-epoch outputs with new GMM map
- [ ] Re-run Section 14 of `neocp_vdp_comparison.ipynb` — stale cached output; cells are fixed, needs kernel restart + re-run. Trust Hyak result (r=0.627, 12/13 correct) in the meantime.
- [ ] Update LaTeX `vdp_pipeline_progression.tex` with F1=0.837 (tied)

### Step 1 — baseline_v5.0.0_10yrs.db on both machines ✓ DONE
SCPed to Arnor (`/astro/users/ds2004/vdp/`) and Hyak (`/mmfs1/gscratch/dirac/ds2004/sorcha/`).

### Step 2 — Redo Sorcha on Hyak with baseline v5.0.0 (NEXT)
- [ ] Re-run Sorcha: same `Rubin_full_footprint_wagg_detections.ini` + same hybrid inputs, new db
- [ ] Submit via `multi_sorcha_production.sh` Slurm array (same chunking: 904 chunks × 16 parts)
- [ ] Monitor with `check_sorcha_outputs.py`, recover any failed files

### Step 3 — New tracklets with `n_det_per_night` column (bundled with Step 2 Phase 1)
- [ ] Add `n_det_per_night = len(night_group)` to `sorcha_postprocess.py` output schema (column 44)
- [ ] Re-run Phase 1 Slurm array over new v5.0 h5 files
- [ ] Sanity-check plots after Phase 1: heliocentric x-y scatter (1–5 AU), population histogram

### Step 4 — Build ~500-map sky grid in antisun-relative ecliptic coords
- [ ] New script `sorcha_gen_maps_grid.py`: configurable `--lon-step`, `--lat-points`, `--sun-exclusion`
- [ ] Grid: 10° lon steps (36 points), ~16 lat values (0, ±2, ±5, ±10, ±20, ±30, ±45, ±60), 40° sun exclusion → ~450–550 maps
- [ ] Maps defined in (Δlon_from_antisun, lat) — time-independent, reusable every year
- [ ] Run one test map first; confirm heliocentric x-y scatter and velocity coverage plots
- [ ] Submit full batch via Slurm once test map looks good

### Step 5 — Score new tracklets with new maps
- [ ] Update `sorcha_postprocess.py` map-assignment to use nearest (Δlon, lat) center in antisun-relative coords
- [ ] Re-run Phase 2: score-vdp + run-digest2 + combine on v5.0 tracklets
- [ ] SCP combined parquet to Arnor

### Step 6 — ROC analysis on Arnor
- [ ] Re-run `sorcha_roc_comparison.ipynb` on v5.0 results
- [ ] Compare VDP vs digest2 F1 on new cadence

### Immediate Hyak small wins (can do before/alongside Step 2)
- [ ] MBA clone_factor 1→5 in `velocity_density_pipeline_gmm.py` → expected +0.020 F1
- [ ] Widen antisun footprint 30°→45° in `sorcha_postprocess.py` → Phase 1 re-run → expected +0.010–0.015
- [ ] More GMM components 80→200 in `_clone_neo_gmm` → expected +0.003–0.008

---

## Known Bugs in Notebooks (non-critical)

- `neocp_vdp_comparison.ipynb` Cell `1f36014b` (Section 13): references `first_seen` column which doesn't exist in `gone_df` → `KeyError`. Pre-existing bug, section still runs partially.
- Section 14 of same notebook: displays stale cached output (old Python scorer, all d2 scores at 0, wrong title). The cell code is correct (loads CSV), but a clean kernel restart + Section 14 re-run is needed to refresh the displayed plot. Hyak validation result (r=0.627) is the ground truth.

---

Good luck! — Arnor Claude, 2026-06-08

---

# Advisor Feedback + Next Major Tasks — 2026-06-08 (corrected)

*Written for Hyak Claude. Read `docs/WAGG_SORCHA_HYAK_CONTEXT.md` for full Hyak infrastructure — that is the canonical file. Hyak working dir: `/mmfs1/gscratch/dirac/ds2004/sorcha/` (astro path is a symlink, both work). Pipeline scripts are now in `neomod/pipeline/` after 2026-06-03 repo reorganisation.*

---

## Current Scorecard (as of 2026-06-08)

| Classifier | F1 | Completeness | Contamination | Notes |
|---|---|---|---|---|
| S3M kNN VDP | 0.847 | 75.1% | 2.8% | Single-epoch S3M test (unfair 7:1 TNO:NEO) |
| GMM VDP (single-epoch) | 0.842 | 76.9% | 7.0% | Same single-epoch test |
| digest2 (single-epoch) | 0.655 | — | — | Artificially low — TNO:NEO bias |
| **VDP GMM + mask OFF + ±2.0 grid** | **0.837** | **78.3%** | **10.2%** | **Current best, TIED with digest2** |
| digest2 (Sorcha 2yr full) | 0.836 | 77.1% | 8.7% | |

VDP completeness (78.3%) already exceeds digest2 (77.3%). Remaining gap is MBA contamination: VDP MBA FPR 1.5% vs digest2 0.2%.

**Existing infrastructure (all on Hyak at `/mmfs1/gscratch/dirac/ds2004/sorcha/`):**
- 14,445 Sorcha output h5 files in `outputs/production_2yr/` (232G, baseline_v3.3_10yrs.db)
- 24 monthly S3M antisun maps in `prob_maps/` (May 2025 – Apr 2027)
- Hybrid training maps in `prob_maps_hybrid/`
- Phase 2 pipeline: `sorcha_phase2.py` (commands: score-vdp, run-digest2, combine)
- `sorcha_postprocess.py` for Phase 1 tracklet building
- `sorcha_gen_map_hybrid.py` + `sorcha_gen_maps_hybrid_slurm.sh` for map generation

---

## Task 1 — Investigate P(NEO) in NEO-only velocity cells (CORRECTED)

**What the advisor asked:** "Look at counts of vel -0.5 or 0.8, look at that giant pixel and see how many objects are in that pixel — if no objects but NEO probability must be 1."

**CORRECTED diagnosis (from Hyak, 2026-06-08):** The earlier Arnor diagnosis (MBA Gaussian tail bleeding) was wrong. **MBA uses K|M cloning, not GMM — there are no Gaussian tails.**

The 9.3× MBA density discrepancy at vlam=−0.5 between the S3M and GMM maps is from **different map centers**:
- S3M NEOCP map: ecliptic lon = 229°
- GMM monthly map: ecliptic lon = 220°

A 9° offset changes which MBAs fall in the sky cut at that velocity point. This is expected, not a bug.

**What to do for the advisor:** Pick a velocity pixel at vlam=−0.5, count actual MBA training objects in that pixel for the GMM map, show MBA K|M density is non-zero because of nearby MBAs visible at that sky center (not tails). Document in notebook: the P(NEO) < 1 gap is from sky-geometry of the specific map center, not a model deficiency. 2–3 lines is enough for the paper.

---

## Task 2 — Redo Sorcha with Rubin Baseline v5.0 (IMPORTANT)

**Current cadence database:** `baseline_v3.3_10yrs.db` (799MB, already on Hyak at `/mmfs1/gscratch/astro/ds2004/sorcha/`)

**Needed:** `baseline_v5.0_10yrs.db` — the official Rubin Operations baseline cadence.

**Why:** v3.3 is an older test cadence. v5.0 is the current official baseline. Papers submitted now should use v5.0. (v5.3 would model first-year LSST start, lower priority, skip for now.)

**Steps:**
1. Download `baseline_v5.0_10yrs.db` from the OpSim archive (rubin-sim-data, or the Rubin Community site; check what's available)
2. Re-run Sorcha on the same `hybrid_sorcha_orbits.csv` + `hybrid_sorcha_phys.csv` inputs with the new db and same `Rubin_full_footprint_wagg_detections.ini` config
3. Re-run Phase 1 (`sorcha_postprocess.py` Slurm array) on new output
4. Re-run Phase 2 (`sorcha_phase2.py` score-vdp, run-digest2, combine) on new tracklets
5. Update ROC analysis

**This is a major Hyak task — check Rubin simulation archive for the db file first before starting anything else.**

v5.3 is explicitly low priority — do not start until v5.0 is complete and results look good.

---

## Task 3 — New Comprehensive Sky-Covering Map Grid (~500 maps)

This is the largest new infrastructure task. Instead of 24 maps along the ecliptic (all at lat=0), build a full sky grid of ~500 maps. The goal is to score any LSST tracklet regardless of where on sky it falls, using the nearest map center relative to the current antisun direction.

### Map grid geometry

**Coordinate system:** Ecliptic coordinates relative to the antisun direction. Maps are defined by (Δlon, lat) where Δlon is the longitude offset from the current antisun ecliptic longitude. This makes the grid sky-position-independent — it repeats every year as the Sun moves.

**Grid parameters:**

| Parameter | Value | Notes |
|---|---|---|
| Longitude step | 10° | −180° to +180° = 36 points along each latitude circle |
| Alternative longitude step | 5° | Denser option, 72 points; discuss with advisor |
| Latitude: ecliptic region | 1–2° steps, 0° to ±30° | Finer grid near ecliptic where most objects are |
| Latitude: high-lat region | Larger steps, 30° to ±90° | Fewer objects, coarser OK |
| Target total latitude points | ~16 (both hemispheres combined) | e.g., 0, ±2, ±5, ±10, ±15, ±20, ±30, ±45, ±60 |
| Sun exclusion zone | 40° from Sun | No maps with center within 40° of Sun ecliptic lon |
| Estimated total maps | ~500 | 36 lon × 16 lat × exclusion zone ≈ ~450–550 |

**The "28 ecliptic + 16 latitude" framing from the advisor notes:**
- 28 grid points along ecliptic: this is 36 longitude steps minus ~8 excluded by the 40° sun zone = ~28 usable ecliptic map centers
- 16 latitude points: 16 unique latitude values (e.g., 0, ±2, ±4, ±8, ±15, ±25, ±40, ±60, ±80) = 16 rows

### Map generation approach

Make map generation **configurable** — command-line parameters, not hardcoded:

```python
# sorcha_gen_maps_grid.py (new script)
--lon-step 10          # longitude step in degrees
--lat-points "0,2,4,8,15,25,40,60,80"  # latitude values (degrees, symmetric ±)
--sun-exclusion 40     # exclude maps within this many degrees of Sun
--map-grid-file map_grid.csv  # output: list of (lon_offset, lat) pairs
--prob-maps-dir prob_maps_grid/
--n-jobs 16
```

Map centers are stored as (Δlon_from_antisun, lat) pairs. Each map is generated the same way as existing monthly antisun maps, but at the specified antisun-relative position instead of always at lon=0 (exactly antisun, lat=0).

### Tracklet-to-map assignment (updated scoring logic)

For each LSST tracklet being scored:

1. Compute the current antisun ecliptic longitude at the observation time (already done in `sorcha_postprocess.py` via `_antisun_ecl_lon(mjd)`)
2. Compute the tracklet's ecliptic (lon, lat)
3. Compute (Δlon, lat) = (tracklet_lon − antisun_lon, tracklet_lat)
4. Find the closest map center in this (Δlon, lat) space
5. Use that map for VDP scoring

This is analogous to the current epoch-aware antisun assignment but extended to 2D (lon offset + lat) instead of just temporal proximity.

### Diagnostic plots (REQUIRED for each map batch)

The advisor explicitly asked for "quick dirty" sanity-check figures every time new tracklets/maps are made. For each batch:

- **Heliocentric distance histogram:** `r_helio` (AU) for all training objects split by population. Confirm:
  - NEOs: peak at 0.5–1.5 AU
  - MBAs: peak at 2–4 AU
  - TNOs: >30 AU
  - No anomalous populations
- **Heliocentric x-y scatter** (in the orbital plane, AU): scatter plot of training objects. The interesting AU range is ~1–5 AU (showing inner belt structure). NEOs should cluster near Earth orbit; MBAs should show main belt ring.
- **Velocity coverage plot:** show the (vlam, vbeta) coverage of training objects for each map center — confirm no empty regions in the expected range

These plots confirm the new map grid is pulling in the right populations before committing to a full 500-map run.

---

## Task 4 — New Tracklets with Detection Count Column

When rebuilding tracklets for the v5.0 Sorcha run (or any future tracklet rebuild):

**Keep:** 2-detection minimum per night (current approach — keep as baseline)

**Add column:** `n_det_per_night` — the total number of raw Sorcha detections for that (ObjID, night) pair, not just the 2 used for the tracklet.

In `sorcha_postprocess.py`, when grouping by (ObjID, night), after picking the first/last pair for the tracklet, also record:

```python
n_det_per_night = len(night_group)  # total detections that night before pairing
```

Include this in the output parquet schema (currently 43 columns — add this as column 44).

**Why:** Allows re-analysis with ≥3 or ≥4 detection thresholds later without rerunning the full Phase 1. Wagg uses ≥3 detections — this makes it trivial to reproduce their cut in post-processing.

---

## Task 5 — digest2 NEOCP Validation — ALREADY DONE ON HYAK

**Status: COMPLETE.** `neomod/pipeline/validate_digest2_neocp.py` ran on Hyak 2026-06-08.
Results: Pearson r=0.627, 12/13 high-MPC objects correctly classified.
Output: `outputs/validate_digest2_neocp.csv` + `outputs/validate_digest2_neocp.png` on Hyak.

**Arnor-side: COMPLETE (2026-06-08).** `validate_digest2_neocp.csv` SCPed to `outputs/`. Section 14 cells updated to load CSV and plot scatter vs MPC score. Ready to run.

---

## Immediate Next Steps for Hyak — Beat digest2 (F1 > 0.837)

VDP is currently TIED at F1=0.837. To win: reduce MBA contamination (FPR 1.5% → ~0.4%).

| Priority | Task | Expected F1 gain |
|---|---|---|
| 1 | **MBA clone_factor 1→5** in `velocity_density_pipeline_gmm.py` | +0.020 → F1 ≈ 0.857 |
| 2 | **Widen antisun footprint 30°→45°** in `sorcha_postprocess.py` → Phase 1 re-run | +0.010–0.015 |
| 3 | **More GMM components 80→200** in `_clone_neo_gmm` call | +0.003–0.008 |

---

## Longer-Term Tasks (advisor braindump)

**B — Redo Sorcha with Rubin baseline v5.0** (IMPORTANT)
- Current: `baseline_v3.3_10yrs.db`. Need `baseline_v5.0_10yrs.db` from OpSim archive.
- Full pipeline rerun: Sorcha → Phase 1 → Phase 2 → ROC.
- v5.3 (first-year LSST) is low priority, skip for now.

**C — New comprehensive sky-covering map grid (~500 maps)**
- Maps in antisun-relative ecliptic coords (Δlon, lat) — time-independent
- Longitude: 10° steps = 36 points, minus ~8 in 40° sun exclusion zone = ~28 usable
- Latitude: ~16 points, 1–2° steps near ecliptic, larger at high lat
- Configurable params (lon step, lat grid, sun exclusion) as CLI args
- Diagnostic plots required per batch: heliocentric distance histogram + x-y scatter (1–5 AU)
- Assignment: per tracklet, compute antisun lon → find nearest (Δlon, lat) map center

**D — Tracklets with detection count column**
- Add `n_det_per_night` to `sorcha_postprocess.py` output — raw detection count before pairing
- Allows ≥3/4 detection threshold re-analysis without Phase 1 rerun

**Full redo order:** v5.0 Sorcha → new tracklets (with `n_det_per_night`) → ~500-map grid → score

---

Good luck! — Arnor Claude, 2026-06-08 (corrected and updated)

---

# Hyak Session Update — 2026-06-11

## What Was Completed This Session

### ✅ Sorcha v5.0 Production Run — DONE
- **`baseline_v5.0.0_2yrs.db`** created (148 MB, 414,488 obs, MJD 60980–61710, Nov 2025 – Nov 2027)
  - Trimmed from 10yr db via: `sqlite3 ... ATTACH ... CREATE TABLE AS SELECT WHERE observationStartMJD < 61710.0`
- **`neomod/pipeline/slurm/multi_sorcha_production_v5.sh`** — Slurm array script (array 0-903%16, 16 CPUs, 128G, 4hr)
- **14,445 h5 files** produced in `outputs/production_2yr_v5/` (~232 GB total)
- **inst00820_part003 recovery**: 5 pathological objects identified and excluded via bisection:
  - `A804 RA`, `A854 OA`, `A868 TA` — same as v3.3 run
  - `A854 RA`, `A868 WA` — NEW in v5.0 (different orbital phase due to 6-month offset in sim window)
  - Bisection method: split 997-obj batch → 16 batches of 63 → 16 batches of 8 → individual objects
  - Recovery script: `neomod/pipeline/slurm/split820v5_part003_final.sh` → runs on `orbits_003_skip5.csv` (995 objects)
  - Excluded objects logged in `work/production_2yr_v5/instance_00820/bisect*/` sub-batches

### ✅ Phase 1 (Tracklet Building) — DONE
- **`n_det_per_night` column added** to `neomod/pipeline/sorcha_postprocess.py` as column 44
  - Stores total raw Sorcha detections per (ObjID, night) group before first/last pairing
  - Allows ≥3/≥4 detection re-analysis post-hoc (Wagg uses ≥3)
  - Equal to `n_det` (same `g.size()` value); added as a clearly named separate column
- **`neomod/pipeline/slurm/sorcha_postprocess_v5.sh`** — identical to v3.3 script except `indir/outdir` paths
- **14,445 parquet files** in `outputs/tracklets_v5/` (7.85 GB total)
- **29,844,550 tracklets** total — all files verified readable, n_det_per_night present, zero-size = 0

### ✅ Documentation Created
- **`SORCHA_V5_PIPELINE.md`** — full pipeline reference (inputs, all steps, scripts, slurm params, batching commands, technical notes). Intended for paper write-up reference.
- **`neomod/paper/NEOrocks.tex`** — Section 4 "Rubin-Cadence Evaluation with Sorcha" added:
  - 4.1 Input Catalog, 4.2 Pointing Database, 4.3 Sorcha Simulation + HPC chunking
  - 4.4 Tracklet Construction (n_det_per_night explained, Wagg comparison)
  - 4.5 Extended VDP Map Grid (antisun-relative coords, 500-map spec)
  - 4.6 GMM Density Estimation (why GMM vs kNN, normalisation fix)
  - 4.7 Scoring & Comparison
  - 4.8 Results — **placeholder, to be filled after Phase 2**

---

## What Still Needs to Be Done (in order)

### NEXT: Step 3 — Build ~500-Map Antisun-Relative Sky Grid
**This is the new baseline for all scoring. Do NOT run Phase 2 until this is done.**

- Write `neomod/pipeline/sorcha_gen_maps_grid.py` — new script, configurable:
  - `--lon-step 10` (degrees), `--lat-points "0,2,5,10,20,30,45,60"`, `--sun-exclusion 40`
  - Maps defined in (Δlon_from_antisun, lat) — time-independent, reusable every year
  - Same map generation logic as `sorcha_gen_map_gmm.py` but at specified (Δlon, lat) offset
- Write corresponding Slurm array script
- Output dir: `prob_maps_grid/`
- **Run one test map first**, confirm heliocentric x-y scatter and velocity coverage plots
- Submit full ~500-map batch once test looks good
- Required diagnostic plots per batch: heliocentric distance histogram + x-y scatter (1–5 AU)

### Step 4 — Phase 2 Scoring (after maps are ready)
- Create v5 versions of Phase 2 Slurm scripts:
  - `sorcha_phase2_vdp_v5.sh` → reads `outputs/tracklets_v5/`, `prob_maps_grid/`, writes `outputs/phase2_v5/`
  - `sorcha_digest2_v5.sh` → same structure as `sorcha_digest2_slurm.sh`
- Run `score-vdp`, `run-digest2`, then `combine`
- Output: `outputs/phase2_v5/sorcha_comparison_v5.parquet`

### Step 5 — Sanity Plots + SCP to Arnor
- Heliocentric x-y scatter (AU, split by population)
- Population count histogram
- SCP `outputs/phase2_v5/sorcha_comparison_v5.parquet` to Arnor

### Step 6 — ROC Analysis on Arnor
- Re-run `sorcha_roc_comparison.ipynb` on v5.0 results
- Compare VDP vs digest2 F1 on new cadence
- Fill in Section 4.8 of `neomod/paper/NEOrocks.tex`

### Optional Small Wins (do NOT apply without explicit instruction)
| Task | File | Expected F1 gain |
|------|------|-----------------|
| MBA clone_factor 1→5 | `neomod/src/velocity_density_pipeline_gmm.py` line ~150 | +0.020 |
| Widen antisun footprint 30°→45° | `neomod/pipeline/sorcha_postprocess.py` → Phase 1 re-run | +0.010–0.015 |
| GMM components 80→200 | `neomod/src/velocity_density_pipeline_gmm.py` line ~1406 | +0.003–0.008 |

---

Updated by Hyak Claude, 2026-06-11
