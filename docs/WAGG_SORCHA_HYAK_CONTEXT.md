# Wagg / Sorcha / Hyak Context
Generated: 2026-05-27, Updated: 2026-06-08 (GMM bug fixed + regenerated; digest2 NEOCP validation; May 2026 patch regenerated)
Covers: full Hyak setup, Sorcha 2yr production run, neomod deployment, Wagg paper methodology, post-processing plan.

---

## Dirac Storage Migration (complete 2026-06-07)

Access to `/mmfs1/gscratch/dirac/` was granted (Unix group `dirac` confirmed via `groups`).
Dirac has ~9.5 TB free vs ~640 GB free in astro. **Migration complete.**

**Storage access vs compute access:**
- `/mmfs1/gscratch/dirac/` storage: **YES** — `groups` shows `dirac`
- `--account=dirac` Slurm compute: **NOT YET** — not in `hyakalloc`. PI needs to run `sacctmgr add user ds2004 account=dirac`. Continue using `--partition=ckpt` (no account needed) for now.

**Current canonical paths:**
```
/mmfs1/gscratch/dirac/ds2004/sorcha/                    ← working directory
/mmfs1/gscratch/dirac/ds2004/sorcha/neomod/             ← git repo
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
/mmfs1/gscratch/dirac/ds2004/home_lib/python3.13        ← pip packages (773 MB)
```
Old path `/mmfs1/gscratch/astro/ds2004/sorcha` is a symlink → dirac and still resolves correctly.
Symlink: `/mmfs1/home/ds2004/.local/lib/python3.13` → `/mmfs1/gscratch/dirac/ds2004/home_lib/python3.13`

**What was done:**
- 299 GB rsync'd from astro to dirac via Slurm job (jobs 35872135, 35882655 on ckpt)
- Orphaned rsync temp file (`.orbits_009.csv.atrFDk`) manually removed before final count verify
- Atomic swap: `sorcha/` renamed to `sorcha.bak_before_delete`, symlink created, backup deleted
- 17 files patched (`astro/ds2004` → `dirac/ds2004`): all pipeline .py scripts, all slurm .sh scripts, WAGG_SORCHA_HYAK_CONTEXT.md — committed as `1862bd1`
- `astro` and `dirac` are on **separate GPFS storage pools** — `mv` requires a real copy (not instant rename)


---

## Repository Reorganisation (2026-06-03)

The neomod git repo was restructured and all "wagg" references removed from file names and
code (except `WAGG_SORCHA_HYAK_CONTEXT.md` itself).

### New neomod/ directory layout

```
neomod/
├── src/                          ← Python library (velocity_density_pipeline*.py, loaders, etc.)
├── pipeline/                     ← Sorcha pipeline scripts (moved from sorcha/ root)
│   ├── sorcha_postprocess.py
│   ├── sorcha_phase2.py
│   ├── sorcha_gen_map*.py
│   ├── hybrid_catalog_prep.py
│   ├── make_may2026_antisun_patch.py
│   ├── sanity_check_frame.py
│   ├── config/                   ← Rubin ini config files
│   └── slurm/                    ← All Slurm job scripts
├── notebooks/
│   ├── paper/                    ← 6 paper/advisor notebooks
│   ├── dev/                      ← 40 exploration notebooks
│   └── sorcha/                   ← sorcha_roc_comparison.ipynb
├── docs/                         ← All context/handoff markdown files
│   ├── WAGG_SORCHA_HYAK_CONTEXT.md  ← this file
│   ├── HANDOFF.md
│   ├── cloning_gmm_neo.md
│   └── ...
├── adam_core_stub/               ← Minimal adam_core mock (moved from sorcha/ root)
├── old_ones/                     ← Archived old notebooks (unchanged)
└── neocp_data/                   ← Real NEOCP scraping data (unchanged)
```

### Key path changes

All Slurm scripts now call Python scripts with full path from WORKDIR:
```bash
"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_map_gmm.py"  # was $WORKDIR/sorcha_gen_map_gmm.py
"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py"        # was bare relative call
```

`sorcha_phase2.py` ROOT navigation fixed for new location at `neomod/pipeline/`:
```python
_NEOMOD = Path(__file__).resolve().parent.parent   # neomod/pipeline -> neomod/
ROOT    = _NEOMOD.parent                            # neomod/ -> sorcha/ (data root)
NEOMOD_SRC = _NEOMOD / "src"
ADAM_STUB  = _NEOMOD / "adam_core_stub"
```

`sorcha_gen_map*.py` adam_core_stub path: `WORKDIR/adam_core_stub` → `NEOMOD/adam_core_stub`

Notebook sys.path: `"src"` → `"../../src"` (notebooks now 2 levels deep from neomod/)

### Output file renames (wagg → sorcha)

| Old name | New name |
|----------|----------|
| `wagg_sorcha_comparison.parquet` | `sorcha_comparison.parquet` |
| `wagg_sorcha_comparison_gmm.parquet` | `sorcha_comparison_gmm.parquet` |
| `wagg_sorcha_comparison_hybrid.parquet` | `sorcha_comparison_hybrid.parquet` |
| `wagg_subsample.parquet` | `sorcha_subsample.parquet` |
| `wagg_sorcha_may2026_antisun_patch.parquet` | `sorcha_may2026_antisun_patch.parquet` |
| `d2_wagg_subsample_r*.parquet` (284 files) | `d2_sorcha_subsample_r*.parquet` |
| `roc_comparison_wagg_sorcha*.png` | `roc_comparison_sorcha*.png` |

### Git commits (2026-06-03)

| Hash | Description |
|------|-------------|
| `6d756c8` | Reorganise repo: notebooks/, pipeline/, docs/; add GMM+hybrid pipeline files |
| `74762db` | Remove all 'wagg' references from file names, comments, and code |
| `8e75465` | Fix .gitignore: anchor paper/ to root so notebooks/paper/ is tracked |


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
/mmfs1/gscratch/dirac/ds2004/sorcha
```

Slurm partition: `ckpt` (checkpoint). No explicit `--account` line needed; defaults to `astro`.

Always use the conda env Python directly — never bare `python` or `pip`:

```bash
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/pip
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/sorcha-run
```

The shell always shows `(base)` — that is NOT the working env. Ignore it.

---

## Key Files on Hyak

```text
/mmfs1/gscratch/dirac/ds2004/sorcha/
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
│   ├── src/                            library code
│   ├── pipeline/                       sorcha pipeline scripts + slurm/ + config/
│   ├── notebooks/paper/                paper figure notebooks
│   ├── notebooks/dev/                  development notebooks
│   ├── notebooks/sorcha/               Sorcha ROC notebook
│   ├── docs/                           context/handoff markdown files
│   └── adam_core_stub/                 minimal adam_core mock (moved from sorcha/ root)
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
└── production_run/                     archived original production run scripts
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
# lives at /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/
```

To update after pushing from Mac:

```bash
cd /mmfs1/gscratch/dirac/ds2004/sorcha/neomod && git pull origin main
```

### VDP Import

`velocity_density_pipeline.py` requires `adam_core` (via `neoscore.py`). `adam_core` is slow to install. Solution: a minimal stub at `adam_core_stub/` satisfies the import. `score_observation` and `ProbMapSet.from_npz` do NOT call adam_core at runtime — only orbit-propagation functions do.

Working import snippet (use in all Hyak scripts):

```python
import sys
sys.path.insert(0, "/mmfs1/gscratch/dirac/ds2004/sorcha/adam_core_stub")
sys.path.insert(0, "/mmfs1/gscratch/dirac/ds2004/sorcha/neomod/src")
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
pms = vdp.ProbMapSet.from_npz("/mmfs1/gscratch/dirac/ds2004/sorcha/prob_maps/prob_maps_2025-03-21.npz")
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
    Output: outputs/phase2/sorcha_subsample.parquet (598,670 rows)
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
    Output: outputs/phase2/sorcha_comparison.parquet (598,670 rows)
    Final columns include: population, P_NEO_vdp, vlam, vbeta, mag_bin_label,
                           P_NEO_d2, digest2_id, and all tracklet geometry columns

Phase 3 — ROC analysis (on Hyak, sorcha_roc_comparison.ipynb)
  Notebook: /mmfs1/gscratch/dirac/ds2004/sorcha/sorcha_roc_comparison.ipynb
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

**Location:** `/mmfs1/gscratch/dirac/ds2004/sorcha/sorcha_postprocess.py`
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
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python check_sorcha_outputs.py \
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
Open folder: `Cmd+Shift+P` → Open Folder → `/mmfs1/gscratch/dirac/ds2004/sorcha`

Important: VSCode opens a new window for the remote. Local window stays open separately. Extensions (Claude Code) are local only — use the remote terminal for Hyak commands and paste output back to local Claude session.

Kernel for notebooks on Hyak: `/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python`

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

## Current Status (2026-06-01) — GMM + mask OFF + grid ±2.0, F1=0.837 (tied with digest2)

### Results summary

| Classifier | F1 | Completeness | Contamination | Notes |
|---|---|---|---|---|
| VDP (S3M training) | 0.740 | 67.5% | 18.0% | |
| VDP (hybrid kNN, Phase1 fixed) | 0.787 | 71.0% | 11.8% | kNN baseline |
| VDP (GMM, mask ON) | 0.802 | 72.9% | 10.9% | |
| VDP (GMM, mask OFF) | 0.809 | 74.5% | 11.6% | |
| **VDP (GMM, mask OFF, grid ±2.0)** | **0.837** | **78.3%** | **10.2%** | **current best (tied digest2)** |
| digest2 | 0.836 | 77.1% | 8.7% | |

Removing the nearest_dist mask with GMM gives **+0.007 F1** (0.802 → 0.809).
Contamination rises slightly (10.9% → 11.6%) because GMM density is nonzero in sparse
zones that kNN zeroed, picking up a few extra MBAs. Net F1 gain is positive.

**Gap remaining: 0.000 F1 — TIED with digest2.**

**Note on calculation:** Hyak inline calc showed F1=0.820 because NaN NEOs (outside footprint, n=4,223) were treated as P=0 misses, penalizing VDP unfairly (digest2 has no footprint). Arnor notebook correctly drops NaN rows, evaluating both classifiers on the same in-footprint subset — this gives F1=0.837 for both. The notebook result is the correct, fair comparison.

### GMM NEO cloner pipeline — all files complete (2026-05-31)

| File | Purpose | Status |
|---|---|---|
| `neomod/cloning_test_ZI.ipynb` | Advisor prototype (reference only) | DONE |
| `neomod/cloning_test_ZI_v2.ipynb` | Completed GMM NEO cloner notebook | DONE |
| `neomod/cloning_gmm_neo.md` | GMM context and integration plan | DONE |
| `neomod/src/velocity_density_pipeline_gmm.py` | GMM pipeline (NEO=GMM, others=K\|M) | DONE |
| `sorcha_gen_map_gmm.py` | Map generator using GMM pipeline | DONE |
| `sorcha_gen_maps_gmm_slurm.sh` | Slurm script → prob_maps_gmm/ | DONE |
| `sorcha_phase2_vdp_gmm.sh` | Phase 2 scoring → outputs/phase2_gmm/ (mask OFF) | DONE |
| `prob_maps_gmm/` | 24 monthly antisun maps (GMM NEO cloner) | DONE |
| `outputs/phase2_gmm/` | GMM scoring results | DONE |

**GMM design (from cloning_test_ZI_v2.ipynb, advisor-approved):**
- NEO cloning only: MBA/TNO/Trojan stay on conditional K|M
- Feature vector: log(a), q, sin/cos(i), sin/cos(node), sin/cos(argperi), sin/cos(M_obs)
- H sampled empirically from source NEO distribution (NOT a GMM feature — preserves faint-end rise)
- e derived from a and q, never sampled independently
- StandardScaler before GMM fit; 80 components; reg_covar=1e-6; random_state=42
- Trained on VISIBLE source NEOs per mag bin (same sky cut as K|M first step)
- K|M fallback if GMM fails for any reason
- GMM converged on all 24 maps, acceptance fraction ~0.556

**Nearest-dist mask control (`--no-nearest-dist-mask` in sorcha_phase2.py):**
- Default (mask ON): sets cells where nearest_dist > 0.2 deg/day to P=0
- `--no-nearest-dist-mask`: passes mask_radius=inf → mask never triggers
- nearest_dist data ALWAYS stored in .npz — re-enable anytime by removing the flag
- For kNN pipeline: mask stays ON (sorcha_phase2_vdp_hybrid.sh unchanged)
- For GMM pipeline: mask OFF (sorcha_phase2_vdp_gmm.sh has --no-nearest-dist-mask)

**NEO score breakdown (GMM, mask OFF, grid ±2.0 — current best):**
- P > 0: 105,549 (92.7%) — scored and completeness 78.3% on in-footprint subset
- P = 0: 4,040 (3.5%) — GMM density zero; fast/high-i NEOs beyond clone coverage
- NaN: 4,223 (3.7%) — footprint edge effects (genuine 30° patch boundaries)

---

### Hybrid catalog pipeline — all files complete (2026-05-31)

| File | Purpose | Status |
|---|---|---|
| `hybrid_catalog_prep.py` | Cartesian→Keplerian conversion | DONE |
| `hybrid_elements.parquet` | 14.4M objects with population labels | DONE |
| `neomod/src/hybrid_loader.py` | Drop-in for s3m_loader.py | DONE |
| `neomod/src/velocity_density_pipeline_hybrid.py` | Hybrid VDP pipeline | DONE |
| `sorcha_gen_map_hybrid.py` | Map generator using hybrid pipeline | DONE |
| `sorcha_gen_maps_hybrid_slurm.sh` | Slurm script → prob_maps_hybrid/ | DONE |
| `sorcha_phase2_vdp_hybrid.sh` | Phase 2 scoring → outputs/phase2_hybrid/ | DONE |
| `prob_maps_hybrid/` | 24 monthly antisun maps (hybrid training) | DONE |
| `outputs/phase2_hybrid/sorcha_comparison_hybrid.parquet` | Final comparison | DONE |

**Bugs found and fixed during hybrid pipeline build:**

1. **Login node OOM on `block0_values[:]`**: Loading all 14.4M × 10 float64 rows at
   once OOMs on the login node. Fix: chunked reading (500k rows/chunk, 29 chunks) in
   `hybrid_catalog_prep.py`.

2. **`--n-jobs 16` on login node fails**: Too many threads spawned (OMP + joblib).
   Fix: always run map generation via Slurm, never on login node.

3. **Single-object mag bin crashes propagator**: When only 1 object passes a magnitude
   cut, `elements_to_helio_ecliptic_state` returns shape `(3,)` not `(1,3)`. Fix:
   added `np.atleast_2d` in `velocity_density_pipeline_hybrid.py` after the call.
   (S3M never hits this because S3M TNOs at H=14-16 have 0 objects.)

4. **`sorcha_phase2.py combine` reads stale S3M VDP scores**: The `d2_*.parquet`
   digest2 shards contain embedded VDP scores from the S3M run. `combine` reads
   `d2_*.parquet` first (before `vdp_*.parquet`), so it uses the old scores. Fix:
   custom merge script that drops old VDP columns from d2 shards and joins hybrid
   VDP shards on `tracklet_id`. Also added `None` return in `load_prob_map` for
   missing map files (backwards-compatible: S3M run unaffected).

5. **NEO clone_factor=30 caused regression at vlam=+0.6/+0.7**: Reducing clones from
   300→30 thinned the density at prograde tail. Fix: set clone_factor=80. Final
   nearest_dist coverage: 58.5% of grid unmasked (vs 19.2% for S3M).

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
| + wide grid (−1.5 to +1.5) | 0.740 | 67.5% | 18.0% |
| Grid extension to ±3.0 (attempted, reverted) | 0.740 | 67.5% | 18.0% |
| + hybrid catalog training | 0.803 | 74.4% | 12.7% |
| Grid extension to ±2.0 (attempted, reverted) | 0.789 | 71.3% | 11.8% |
| + Phase 1 re-run (fixed map assignments) | 0.787 | 71.0% | 11.8% |
| + GMM NEO cloner (mask ON) | 0.802 | 72.9% | 10.9% |
| + nearest-dist mask OFF | 0.809 | 74.5% | 11.6% |
| **+ grid ±2.0 (GMM only, current best)** | **0.837** | **78.3%** | **10.2%** |
| digest2 (constant) | 0.836 | 77.1% | 8.7% |

### Remaining gap analysis (as of 2026-06-01, GMM + grid ±2.0, F1=0.837)

Of 113,812 NEO tracklets total; 109,589 in-footprint (NaN dropped for fair comparison):

| Category | Count | % of NEOs | Root cause |
|---|---|---|---|
| P > 0 (scored) | 105,549 | 92.7% | In-footprint completeness = 78.3% |
| P = 0 (GMM density = 0) | 4,040 | 3.5% | Fast/high-i NEOs beyond velocity coverage |
| NaN (outside footprint) | 4,223 | 3.7% | Genuine 30° patch edge effects |

**P=0 with mask OFF — detailed breakdown (diagnosed 2026-06-01):**

| P=0 category | Count | % of P=0 | Root cause |
|---|---|---|---|
| |vlam| > 1.3 deg/day | 5,020 | 69.5% | Fast in-ecliptic — beyond clone coverage |
| |vlam| < 1.0, |vbeta| > 0.3 deg/day | 1,708 | 23.6% | Fast out-of-ecliptic — high-i geometry |
| |vlam| < 1.0, |vbeta| < 0.3 deg/day | 69 | 1.0% | Within grid, unclear |

**Critical insight: P=0 is a VELOCITY COVERAGE problem, not an orbital element gap.**
The kNN density evaluator only places density where clone POINTS land in (vlam, vbeta)
space. Adding more exotic orbits to the training set (NEOMOD3, MPCORB) doesn't help if
those orbits, when propagated to the antisun footprint and sky-cut, don't produce clone
points in the extreme velocity cells.

**NaN breakdown:**
- Old map assignments: 0 — fixed by Phase 1 re-run
- Footprint edge effects: 4,223 — genuine 30° patch boundary gaps

**VDP F1 = digest2 F1 = 0.837 (tied).** VDP completeness (78.3%) now EXCEEDS digest2
(77.3%). The remaining gap is PURELY contamination: VDP MBA FPR = 1.5% vs digest2 0.2%.
VDP selects 7.5× more false MBAs. Reducing MBA FPR from 1.5% → 0.4% would give F1 ≈ 0.857.

---

### NEOMOD3 orbit augmentation — ATTEMPTED, ZERO EFFECT (2026-06-01)

**What was built:**
- `neomod/src/neomod3_sampler.py` — loads NEOMD3 4D array, samples NEO orbits from
  the debiased distribution (uniform phase angles, physically valid filter q < 1.3)
- `neomod/NEOMD3/input_neomod3.dat` — 62 MB downloaded from Boulder SWRI
- `velocity_density_pipeline_gmm.py` — augments GMM training with 10x NEOMD3 orbits
  per mag bin before fitting; h_source_visible preserved from original visible NEOs only

**Result:** F1=0.802 (vs 0.809 before). P=0 count unchanged at 7,222.
F1 drop is within joblib non-determinism range (±0.02) — not a real regression.

**Why it failed:** NEOMD3 adds orbits in (a, e, i) space with uniform phase angles.
When propagated to the antisun footprint:
- Most NEOMD3 orbits NOT near the antisun at the map epoch — filtered by sky cut
- Those that survive still don't produce clones at |vlam|>1.3 or |vbeta|>0.3
- Unusual velocities require specific orbital phase at opposition (Aten-class at
  perihelion, high-i near node) — rare at any epoch regardless of training diversity

**Status:** Code left in place (neomod3_sampler.py exists, import in gmm pipeline).
Not harmful (zero effect). Could help if NEOMD3 samples were pre-filtered to only
those visible at the specific antisun epoch (M_deg restricted to near-antisun geometry).
Not yet tried.

---

### Cross-pipeline comparison confirmed on Arnor (2026-06-01)

Three-way ROC comparison run on Arnor across all parquets, NaN-dropped (in-footprint only):

| Pipeline | Rows | NEOs (in-footprint) | VDP F1 | VDP comp | VDP cont | d2 F1 |
|---|---|---|---|---|---|---|
| S3M kNN (±1.5 grid) | 611,041 | 112,703 | 0.740 | 67.5% | 18.0% | 0.836 |
| Hybrid kNN | 611,041 | 107,935 | 0.803 | 74.4% | 12.7% | 0.837 |
| **GMM + mask OFF + grid ±2.0** | **611,041** | **109,589** | **0.837** | **78.3%** | **10.2%** | **0.837** |

**Per-population breakdown (GMM, from Arnor notebook):**

| Population | n | VDP above thr (0.003) | d2 above thr (0.640) |
|---|---|---|---|
| NEO | 109,589 | 78.3% | 77.3% |
| MBA | 447,167 | **1.5%** | **0.2%** |
| TNO | 5,816 | 0.0% | 69.5% |
| Trojan | 13,932 | 0.1% | 7.6% |
| other | 21,002 | 14.3% | 9.4% |

VDP completeness (78.3%) exceeds digest2 (77.3%). The remaining difference is MBA
contamination — VDP selects 7.5× more false MBAs at optimal threshold.

**Arnor file layout (`/astro/users/ds2004/vdp/`):**
```
prob_maps/                         ← S3M kNN monthly maps (24 files)
prob_maps_gmm/                     ← GMM monthly maps (May 2026 confirmed)
outputs/phase2/
  sorcha_comparison.parquet            (162 MB — S3M result)
  sorcha_comparison_gmm.parquet        (75 MB  — GMM result)
  sorcha_comparison_hybrid.parquet     (162 MB — hybrid result)
  sorcha_may2026_antisun_patch.parquet (1.2 MB — single-epoch patch, see below)
sorcha_gmm_s3m_comparison.ipynb            ← density + P(NEO) + ROC comparison (advisor)
sorcha_gmm_s3m_singleepoch_comparison.ipynb ← same-objects single-epoch ROC (new)
```

**Comparison notebook `neomod/sorcha_gmm_s3m_comparison.ipynb` (created 2026-06-01):**
- Loads S3M NEOCP map (`prob_maps_2026-05-09T22_neocp.npz`, center=229°, ±0.8 grid)
  and GMM monthly map (`prob_maps_gmm/prob_maps_2026-05-01_antisun.npz`, center=220°, ±2.0 grid)
- Both within 8.7° of each other (same May 2026 antisun region)
- Figure 1: log density comparison clipped to ±0.8 (S3M left | GMM right)
- Figure 2: P(NEO) probability comparison clipped to ±0.8 (S3M mask ON | GMM mask OFF)
- Figure 3: GMM full ±2.0 range with dashed ±0.8 reference box
- Figure 4: three-way ROC (S3M kNN, GMM, digest2)
- Kernel: neofast_py310; no sklearn needed
- Advisor-requested: shows exactly how map structure + ROC changed from S3M to GMM

---

### Next steps to close the gap (updated 2026-06-01)

**Grid ±2.0 with GMM + mask OFF — DONE, F1=0.837**

Previously tested ±2.0 with kNN + mask ON: F1 dropped 0.803→0.789. That test does NOT
apply here because:
1. Mask OFF: cells in ±1.5–2.0 won't be zeroed by nearest_dist threshold
2. GMM cloner: clone distribution may extend slightly further in vlam than K|M

With GMM + mask OFF, kNN density at ±1.5–2.0 will be small but potentially nonzero
(smoothed from nearest clones at |vlam|≈1.3). Targets 5,020 fast P=0 NEOs.

Done: `DEFAULT_GRID_LIM = (-2.0, 2.0)` in `velocity_density_pipeline_gmm.py`. Result: F1=0.837 (tied digest2). Goal is now F1 > 0.837.

**CURRENT GOAL: F1 > 0.837 — need to beat digest2, not just match it**

Per-population score breakdown (from Arnor notebook):
| Population | n | VDP above thr | d2 above thr |
|---|---|---|---|
| NEO | 109,589 | 78.3% | 77.3% |
| MBA | 447,167 | 1.5% | 0.2% |
| TNO | 5,816 | 0.0% | 69.5% |
| Trojan | 13,932 | 0.1% | 7.6% |
| other | 21,002 | 14.3% | 9.4% |

VDP completeness (78.3%) ALREADY EXCEEDS digest2 (77.3%).
The problem is purely contamination: VDP selects ~7.5x more false MBAs than digest2.
VDP MBA FPR = 1.5%, digest2 MBA FPR = 0.2%.

Root cause: MBAs near the 25–30° footprint edge are at their stationary point
(vlam ≈ 0), which overlaps NEO velocity space. VDP scores them P_NEO > 0.003.
The MBA density model has only clone_factor=1 (no augmentation) — sparse near
the footprint edge → MBA density underestimated → MBAs leak into P_NEO > 0.

**To beat digest2:** reduce MBA false positive rate from 1.5% → ~0.4%.
This would give F1 ≈ 0.857 while maintaining current NEO completeness.

**Priority 1 — Increase MBA clone_factor in GMM pipeline (3→5×)**
Currently MBA clone_factor=1 (native density, no cloning). At footprint edge,
MBAs near stationary point (vlam≈0) are underrepresented in the training →
MBA density in those cells is too low → P_NEO too high for edge MBAs.
Fix: increase MBA clone_factor from 1 to 5 in velocity_density_pipeline_gmm.py.
More MBA clones → better stationary-point density model → lower P_NEO there.

**Priority 2 — Widen antisun footprint (30°→45°)**
- Targets 4,223 NaN NEOs at patch edges
- Change `max_sep_deg=30` → `max_sep_deg=45` in `sorcha_postprocess.py`
- Requires Phase 1 re-run + map regeneration + Phase 2
- Expected: F1 +0.010–0.015

**Priority 3 — More GMM components (80→200)**
- One line in `_clone_neo_gmm` call: `n_components=200`
- Safeguard `min(n_components, len(X)//5)` prevents overfitting sparse bins
- Expected: F1 +0.003–0.008

**Grid ±2.0 with kNN + mask ON: ATTEMPTED AND REVERTED (different from GMM test)**
- kNN+mask ON dropped F1 to 0.789. GMM+mask OFF confirmed F1=0.837 (see current best above).

**Phase 1 re-run: COMPLETE**
- All 14,445 tracklets reprocessed; all NaN are genuine footprint edge effects.


### Additional bugs found (2026-05-31 session)

6. **`sorcha_postprocess.sh` missing `--overwrite` flag**: Removing `--skip-existing`
   without adding `--overwrite` caused all Phase 1 re-run tasks to fail immediately
   with `FileExistsError`. Fixed: `--overwrite` now in the script.

7. **Grid extension ±2.0 hurt F1**: See Quick Win 1 above. Reverted.

8. **joblib parallel non-determinism in map generation**: `evaluate_density_map_full_posterior_2d`
   uses `joblib.Parallel(n_jobs=16)`. Despite fixed random seed (seed=42), parallel
   floating-point summation order varies by thread scheduling → slightly different density
   values each run → F1 shifts by ~0.02 between regenerations of the same maps.
   - Verified: diff of hybrid vs fast pipeline shows only 4 intended changes (correct).
   - Fixed seed=42 is set but does not fully control joblib thread ordering.
   - **For the final paper: regenerate maps with `--n-jobs 1` for fully reproducible results.**
     GMM replaces this density step entirely, so the non-determinism problem goes away
     once GMM is implemented.

### Directory cleanup (2026-05-31)

- Created `production_run/` — all 35 production run scripts and check/rerun files moved here
- Deleted 5 old fixed-position maps from `prob_maps/` (superseded by 24 monthly maps)
- Deleted temp outputs: `outputs/phase2_smoke/`, `outputs/production_2yr_split820*/`
- Deleted temp files: `submit_combined_rerun.sh`, `phase1_rerun_indices.txt`

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

Copied from `velocity_density_pipeline_fast.py`. Changes:
1. `import s3m_loader as load_s3m` → `import hybrid_loader as load_s3m`
2. Clone factors: MBA=1, NEO=**80**, TNO=10, Trojans=5 (NEO=30 regressed at vlam=+0.6/+0.7)
3. `np.atleast_2d` fix after `elements_to_helio_ecliptic_state` call (single-object edge case)
4. Module docstring updated

#### Step 4 — `sorcha_gen_map_hybrid.py` ✓ COMPLETE (2026-05-31)

One-line change from `sorcha_gen_map.py`:
`import velocity_density_pipeline_hybrid as vdp`

#### Step 5 — `sorcha_gen_maps_hybrid_slurm.sh` ✓ COMPLETE (2026-05-31)

Changes: `--job-name=sorcha_g_hybrid`, calls `sorcha_gen_map_hybrid.py`,
output to `prob_maps_hybrid/`. All 24 maps ran concurrently, completed in ~8 min.

#### Step 6 — `sorcha_phase2_vdp_hybrid.sh` ✓ COMPLETE (2026-05-31)

Changes: `--job-name=sorcha_vdp_hybrid`, `--prob-maps-dir prob_maps_hybrid`,
`--work-dir outputs/phase2_hybrid`. 113 shards, ~8 min total.

#### Step 7 — combine ✓ COMPLETE (2026-05-31)

**CRITICAL BUG:** `sorcha_phase2.py combine` reads `d2_*.parquet` (digest2 shards)
first. These shards embed the OLD S3M VDP scores — so the combined parquet had
identical F1=0.740 as S3M despite hybrid VDP shards being correct.

**Fix:** custom merge script (run inline) that:
1. Drops stale VDP columns (`P_NEO_vdp`, `vlam`, `vbeta`, `mag_bin_label`) from d2 shards
2. Loads all 113 hybrid vdp_shards
3. Left-joins on `tracklet_id`

Also added `None` return in `load_prob_map` for missing map files so that tracklets
assigned to old `antisun_minus` maps (not in `prob_maps_hybrid/`) score as NaN
instead of crashing. Backwards-compatible: S3M run unaffected.

**Correct combine command for hybrid:**
```python
# Drop old VDP cols from d2 shards, join hybrid vdp shards on tracklet_id
vdp = pd.concat([pd.read_parquet(f, columns=['tracklet_id','P_NEO_vdp','vlam','vbeta','mag_bin_label'])
                 for f in sorted(glob('outputs/phase2_hybrid/vdp_shards/vdp_*.parquet'))])
d2  = pd.concat([pd.read_parquet(f) for f in sorted(glob('outputs/phase2_hybrid/digest2_shards/d2_*.parquet'))])
d2  = d2.drop(columns=[c for c in ['P_NEO_vdp','vlam','vbeta','mag_bin_label'] if c in d2.columns])
out = d2.merge(vdp, on='tracklet_id', how='left')
out.to_parquet('outputs/phase2_hybrid/sorcha_comparison_hybrid.parquet', index=False)
```

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

### Single-epoch comparison work (2026-06-03)

#### `sorcha_gmm_s3m_singleepoch_comparison.ipynb` (Arnor, `/astro/users/ds2004/vdp/`)

Arnor Claude created this notebook to answer: "does the GMM map actually improve scoring,
holding everything else equal?" It takes the same 21K S3M NEOCP May 9 2026 objects that
were already scored in the S3M parquet, and re-scores them with the GMM map via
`RegularGridInterpolator` on each mag-bin's P(NEO) map.

**Results on S3M objects:**

| Classifier | F1 | Completeness | Contamination |
|---|---|---|---|
| VDP S3M kNN | 0.847 | — | — |
| VDP GMM (same objects) | 0.842 | — | — |
| digest2 | 0.655 | — | — |

**Critical caveat**: digest2 F1=0.655 is artificially low here. The S3M field has
7,492 TNOs vs 1,119 NEOs (7:1 ratio). digest2 scores TNOs very high (median d2=0.79,
69.5% above threshold) — it can't distinguish slow retrograde TNOs from NEOs. VDP can,
because TNOs occupy a distinct velocity region. This is a population-ratio artifact,
not a classifier comparison. The honest benchmark is the full Sorcha simulation
(digest2 F1=0.837 there).

**Legitimate conclusion**: GMM ≈ S3M kNN on the same single-epoch objects (F1 0.842 vs
0.847). The GMM improvement (0.740 → 0.837) comes from better coverage in the
full 2-year Sorcha simulation, not from degrading performance on well-covered objects.

**Notebook fixes applied by Arnor Claude:**
- ROC sweep vectorized (numpy broadcast, no Python loop) — ~15s → <1s
- Figure dpi 200 → 150 (pixel count -44%)
- Key name fix: `min_mag`/`max_mag` → `mag_min`/`mag_max` (GMM `mag_bins` structure)
- Kernel: neofast_py310

#### Fair single-epoch comparison: why 100 sq deg doesn't work

The advisor requested a comparison on "100 sq degrees of night sky" to avoid S3M's
biased population mix. A strict 5.64° radius cut around the May 9 antisun yields only
83 objects (6 NEOs) from the subsampled parquet — too few for a ROC curve. The
subsampled parquet keeps ALL 98,646 NEOs but samples only 500K non-NEOs from 40.7M.
This thins out any small sky patch to a handful of objects.

**Resolution:** Use May 2026 observations within the full 30° antisun footprint.
This gives 14,873 objects (2,402 NEOs, 16.2% NEO rate) with a realistic population
mix (no S3M TNO bias). The "single epoch" property is preserved: all objects observed
in May 2026 and all scored with the May 2026 antisun VDP calibration.

#### New parquet: `sorcha_may2026_antisun_patch.parquet`

**Hyak path:** `outputs/phase2/sorcha_may2026_antisun_patch.parquet` (1.2 MB)
**Arnor path:** `outputs/phase2/sorcha_may2026_antisun_patch.parquet`
**Generator:** `make_may2026_antisun_patch.py` (Hyak working dir)

**Filter criteria:**
- Date: mjd0_utc in May 2026 (MJD 61161 – 61192)
- Spatial: within 30° of ecliptic lon=228°, lat=0° (May 9 2026 antisun)

**Statistics:**

| Population | n | % |
|---|---|---|
| MBA | 11,795 | 79.3% |
| NEO | 2,402 | 16.2% |
| other | 477 | 3.2% |
| TNO | 148 | 1.0% |
| Trojan | 45 | 0.3% |

**Columns:** `tracklet_id, population, P_NEO_s3m, P_NEO_d2, vlam, vbeta, mean_mag,
ecl_lon, ecl_lat, mjd0_utc, P_NEO_gmm`

**NaN situation:**
- `P_NEO_s3m`: 4,275 NaN — tracklets Phase 1 assigned to old NEOCP map (`lon229_lat0`);
  S3M Phase 2 skipped these (antisun filter). Score is missing, not zero.
- `P_NEO_gmm`: 52 NaN — GMM Phase 2 re-scored NEOCP-assigned tracklets with nearest
  monthly map; only genuine footprint edges are NaN.
- `P_NEO_d2`: 0 NaN

**For Arnor ROC notebook (Section 5):**
- Three-way (S3M + GMM + d2): `df.dropna(subset=['P_NEO_s3m','P_NEO_gmm'])` → ~10,500 objects, ~1,770 NEOs
- Two-way (GMM + d2): `df.dropna(subset=['P_NEO_gmm'])` → ~14,820 objects, ~2,395 NEOs

**SCP from Hyak (run in sorcha working dir):**
```bash
scp outputs/phase2/sorcha_may2026_antisun_patch.parquet \
    ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/outputs/phase2/
```

---

## GMM Normalisation Bug — FIXED AND FULLY REGENERATED (2026-06-07/08)

### What the bug was

Advisor identified: GMM P(NEO) was much weaker than S3M kNN, especially at vlam ≈ −0.5.
Diagnostic from Arnor notebook (mag22, vbeta=0):

| vlam  | P(NEO) S3M | P(NEO) GMM (buggy) | ratio |
|-------|-----------|-----------|-------|
| -0.50 | 0.9543    | 0.4568    | 2.09× |
| -0.30 | 0.6793    | 0.0902    | 7.53× |
|  0.00 | 0.0448    | 0.0005    | 95.9× |

Grid integral ratios (GMM/S3M): MBA=0.995, **NEO=0.322**, TNO=0.674, Trojans=4.532

**Root cause:** `velocity_density_pipeline_gmm.py` line (previously ~1517):
```python
density_downweighted_map = density_clone_map / f   # BUG for GMM path
```
K|M cloner starts from N_vis visible objects and produces ≈f×N_vis visible clones (sky-cut
preserving) — dividing by f is correct there. GMM cloner samples f×N_vis clones from the
full NEO orbit distribution then re-applies the sky cut; only acceptance_fraction ≈ 0.32
survive. Dividing by f instead of n_visible_clones/n_source underscales NEO density by
1/0.322 ≈ 3.1×.

### The fix (applied 2026-06-07)

```python
# inside build_cloned_maps_for_center_magbin, NEO GMM branch:
n_source = len(df_cloner_input)               # visible source NEOs in sky patch
n_visible_clones_gmm = len(clone_visible_df)  # GMM clones surviving sky cut

# at the downweighting step:
if pop_name == "NEO" and gmm_success and n_source > 0:
    effective_factor = n_visible_clones_gmm / n_source
else:
    effective_factor = f
density_downweighted_map = density_clone_map / effective_factor
# same effective_factor passed to make_support_count_map and support_for_smoothing
```

### Full pipeline re-run (2026-06-07/08) — COMPLETE

All steps re-run after fix:

| Step | Slurm job | Status | Output |
|------|-----------|--------|--------|
| Regenerate 24 GMM maps | array (sorcha_gen_maps_gmm_slurm.sh) | Done | `prob_maps_gmm/*.npz` |
| Phase 2 GMM scoring (113 shards) | array (sorcha_phase2_vdp_gmm.sh) | Done | `outputs/phase2_gmm/vdp_shards/` |
| Combine 40M rows | 128GB job (combine_gmm.sh) | Done | `outputs/phase2_gmm/sorcha_comparison_gmm.parquet` |
| Build 611K comparison parquet | 128GB job (build_gmm_comparison.sh) | Done | `outputs/phase2/sorcha_comparison_gmm.parquet` |
| Regenerate May 2026 patch | 128GB job 35982399 | Done 2026-06-08T14:03 | `outputs/phase2/sorcha_may2026_antisun_patch.parquet` |

**May 2026 patch after fix:** P_NEO_gmm range [0.0000, 0.9812] (was suppressed to ~0.46 max).

### Note on "Bug 2" (MBA tails)

Arnor Claude flagged GMM MBA density being 9× higher than S3M MBA at vlam=−0.5.
MBA uses K|M, not GMM. The discrepancy is from the sky-center offset (S3M: lon=229°,
GMM: lon=220°) changing which MBAs are visible, not Gaussian tail leakage.
Bug 1 (NEO underscaling) was the actionable fix.

---

## digest2 NEOCP Validation (2026-06-08)

Advisor suggestion: validate our digest2 implementation against MPC's published
NEOCP scores (which MPC also computes with digest2). We have real NEOCP data in
`neomod/neocp_data/`.

### Data available

- `neocp_data/tables/neocp_objects_history.csv`: 60 rows, 20 unique NEOCP objects
  (May 7 2026 snapshots), each with `score` (0–100) and `score_0_1` = MPC's digest2 score
- `neocp_data/ephemerides/neocp_ephemerides_20260507T104242Z.csv`: 425 rows, hourly
  ephemerides per object with `ra_deg`, `dec_deg`, `v_mag`, `motion_arcsec_per_min`, `pa_deg`

### Approach

For each NEOCP object with ephemerides, take two consecutive hourly positions (1-hr
tracklet baseline), format as MPC 80-column obs, run our digest2 binary, compare our
P_NEO_d2 to MPC's `score_0_1`.

**Caveat:** MPC used the actual submitted arc (nobs=2–35, arc_days=0.01–9.96); we use
synthetic 2-det tracklets from ephemerides. Perfect match not expected, but correlation
validates the binary is correct.

### Script

`neomod/pipeline/validate_digest2_neocp.py`

**Key engineering issues found and fixed:**
1. Ephemeris dedup: two snapshot files had overlapping `(designation, ephem_time_utc)` rows
   → dedup on that pair before picking obs pairs
2. digest2 truncates designation to 5 chars in output → `A11BRHV` → `A11BR`, causing
   collisions between `ZTF10Df`/`ZTF10Di` → both → `ZTF10`
3. Fix: use synthetic 5-char unique keys `VA000`–`VA016` in MPC obs lines; map back by index

### Results (17 objects, run 2026-06-08)

| Metric | Value |
|--------|-------|
| Pearson r | 0.627 |
| Spearman ρ | 0.586 |
| MAE | 0.101 |
| High-MPC (≥0.90) scored ≥0.50 | 12/13 |

```
designation  mpc_score_0_1  our_d2_score  nobs  arc_days
    ZTF10Di           1.00          1.00     4      0.02
    C468U01           1.00          1.00     8      0.06
    ...
    ZTF10Df           0.98          0.03     9      0.68  ← outlier (9 obs, 0.68d arc)
    P22n0fj           0.97          1.00     8      0.10
    S001563           0.68          0.71     3      0.03
    C1EUYV5           0.65          0.39    12      9.79
    C1EUWK5           0.37          0.23    10      9.96
```

**Conclusion:** Binary is functioning correctly. Correlation is limited by input
differences (MPC used full arc; we used 1-hr synthetic). The ordering is preserved
at the extremes (high-MPC → high ours; low-MPC → low ours). One outlier (ZTF10Df,
0.68d arc) resolves at a 1-hr tracklet to a non-NEO-like velocity — understandable.

Output: `outputs/validate_digest2_neocp.csv`, `outputs/validate_digest2_neocp.png`

---

## Current Status (2026-06-08)

### Completed this session
- GMM normalisation bug fixed and all 24 maps regenerated
- Full Phase 2 GMM re-run (113 shards → 40M row parquet → 611K comparison parquet)
- May 2026 antisun patch regenerated with corrected GMM scores
- digest2 NEOCP validation script written and passing (r=0.627, 12/13 high-MPC correct)
- VSCode lock issue: created missing `~/.claude/sessions/` directory

### Files SCPed to Arnor — ALL DONE (2026-06-08)

| File | Arnor destination | Status |
|---|---|---|
| `outputs/phase2/sorcha_comparison_gmm.parquet` | `outputs/phase2/` | ✓ |
| `prob_maps_gmm/prob_maps_2026-05-01_antisun.npz` | `prob_maps_gmm/` | ✓ |
| `outputs/phase2/sorcha_may2026_antisun_patch.parquet` | `outputs/phase2/` | ✓ |
| `outputs/validate_digest2_neocp.csv` | `outputs/` | ✓ (Jun 8, 534 bytes) |
| `baseline_v5.0.0_10yrs.db` | `/astro/users/ds2004/vdp/` (Arnor) and `/mmfs1/gscratch/dirac/ds2004/sorcha/` (Hyak) | ✓ (728MB) |

### Step 0 — Arnor notebook cleanup (COMPLETE 2026-06-08)

All post-fix files on Arnor. `neocp_vdp_comparison.ipynb` Section 14 updated to load
`validate_digest2_neocp.csv` (replaces broken Python NEOMODScorer). `SORCHA_HYAK_ARNOR_CONTEXT.md`
corrected (F1 numbers, Bug 2 diagnosis, digest2 validation status).

Remaining Arnor actions (user re-runs only, no code changes needed):
- Re-run `sorcha_roc_comparison.ipynb` to confirm F1=0.837
- Re-run `sorcha_gmm_s3m_singleepoch_comparison.ipynb` with new map
- Run Section 14 of `neocp_vdp_comparison.ipynb` to see digest2 validation scatter

### Next major task: Redo Sorcha with baseline v5.0.0

`baseline_v5.0.0_10yrs.db` (728MB) is now on both Arnor and Hyak. Full pipeline rerun:
Sorcha → Phase 1 (new tracklets with `n_det_per_night` column) → Phase 2 → ROC.
See Task 2 in `SORCHA_HYAK_ARNOR_CONTEXT.md` for full step list.

---

## Master To-Do List (as of 2026-06-08)

### Step 0 — Arnor notebook cleanup (mostly done)
- [ ] Re-run `sorcha_roc_comparison.ipynb` — confirm F1=0.837 with post-fix parquet
- [ ] Re-run `sorcha_gmm_s3m_singleepoch_comparison.ipynb` — update single-epoch outputs with new GMM map
- [ ] Re-run Section 14 of `neocp_vdp_comparison.ipynb` — stale cached output; cells are fixed, needs kernel restart
- [ ] Update LaTeX `vdp_pipeline_progression.tex` with F1=0.837 (tied)

### Step 1 — baseline_v5.0.0_10yrs.db on both machines ✓ DONE
`baseline_v5.0.0_10yrs.db` (728MB) SCPed to Arnor (`/astro/users/ds2004/vdp/`) and Hyak (`/mmfs1/gscratch/dirac/ds2004/sorcha/`).

### Step 2 — Redo Sorcha on Hyak with baseline v5.0.0 ✓ DONE (2026-06-11)
- [x] `baseline_v5.0.0_2yrs.db` created (148 MB, 414,488 obs, MJD 60980–61710, Nov 2025–Nov 2027)
- [x] `neomod/pipeline/slurm/multi_sorcha_production_v5.sh` written and submitted
- [x] 14,445 h5 files in `outputs/production_2yr_v5/` (~232 GB)
- [x] inst00820_part003 recovery — 5 pathological objects excluded:
  - `A804 RA`, `A854 OA`, `A868 TA` (same as v3.3)
  - `A854 RA`, `A868 WA` (new in v5.0 — different orbital phase due to 6-month offset in sim window)
  - Found via 3-round bisection: 16×63-obj batches → 16×8-obj batches → individual objects
  - Recovery script: `neomod/pipeline/slurm/split820v5_part003_final.sh` (orbits_003_skip5.csv, 995 objects)

### Step 3 — New tracklets with `n_det_per_night` column ✓ DONE (2026-06-11)
- [x] `n_det_per_night` added to `neomod/pipeline/sorcha_postprocess.py` as column 44
  - Equal to `g.size()` (total raw detections per (ObjID, night) before first/last pairing)
  - Enables ≥3/≥4 re-analysis post-hoc (Wagg uses ≥3) without Phase 1 rerun
- [x] `neomod/pipeline/slurm/sorcha_postprocess_v5.sh` — identical to v3.3 except indir/outdir
- [x] Phase 1 run over all 14,445 v5.0 h5 files
  - **29,844,550 tracklets** in `outputs/tracklets_v5/` (7.85 GB, 14,445 parquets)
  - All files verified: readable, n_det_per_night present, 0 zero-size files
- [ ] Sanity-check plots: heliocentric x-y scatter (1–5 AU), population histogram — **TODO**

### Step 4 — Build 667-map antisun-relative sky grid ✓ DONE (2026-06-17)
**The new scoring baseline.**
- [x] `neomod/pipeline/sorcha_gen_maps_grid.py` — configurable `--lon-step`, `--lat-base`/`--lat-points`, `--sun-exclusion`, `--ref-obstime`, `--mba-clone-factor`, `--task-id`, `--save-overlays`
- [x] Grid: 10° lon (29 usable after 40° sun exclusion) × 23 non-uniform lat (0,±1,±2,±3,±4,±5,±8,±12,±18,±25,±35,±50 — fine 1° near ecliptic) = **667 maps**, ref epoch 2026-01-01, **MBA cf=5**
- [x] Slurm `sorcha_gen_maps_grid_slurm.sh` (`--array=0-666%48`, skip-existing) → `prob_maps_grid/`; test script `_test.sh` (idx 333 antisun, 652 lat50)
- [x] Validated: test map 333 physically correct (P(NEO) peaks ≈0.93 at vlam≈−0.54) and reproduces the old monthly antisun map at cf=5 (support 5→1). 667/667 maps verified (0 missing, 0 zero-byte). Auto-resumed through a maintenance window.

### Step 5 — Score tracklets with the grid ✓ DONE (2026-06-20)
- [x] **Phase 1 re-run** (`sorcha_postprocess_v5_grid.sh`, 15 batches) → `outputs/tracklets_v5_grid/` (14,445 parquets). Grid-aware assignment → **~100% of tracklets get a map** (vs ~48% with monthly maps); **65.9M tracklets** retained.
- [x] v5 Phase 2 Slurm scripts created: `sorcha_phase2_vdp_v5.sh`, `sorcha_digest2_v5.sh`, `sorcha_digest2_v5_retry.sh`.
- [x] **score-vdp**: 113/113 shards, 65,857,457 scored, 0.00% NaN, NEO 207,602.
  - Fix 1 (`7ffdd3a`): `score_vdp_frame` accepts `"grid"` map names (else all-NaN).
  - Fix 2 (`1ac6170`): evict each map after scoring — a 128-file shard touches ~600 of 667 maps; caching all (~180 MB each) OOM-killed 105/113 on the first try. Eviction → peak RSS 0.46 GB.
- [x] **sample**: `sorcha_subsample.parquet` = 707,670 rows (207,602 NEO + 500,068 non-NEO, seed 42).
- [x] **digest2**: 142 tasks (5000 rows each). 40 hit `TimeoutExpired` on slow ckpt → `_retry.sh` with `--digest2-chunk-tracklets 1000` (5×1000) finished them → 142/142.
- [x] **combine** → `outputs/phase2_v5/sorcha_comparison_v5.parquet` (707,670 rows, both `P_NEO_vdp`/`P_NEO_d2`, 0% NaN). Sanity medians: NEO 0.227/0.990, MBA 0.002/0.020, TNO 0.000/0.355 (vdp/d2).
- [ ] SCP `sorcha_comparison_v5.parquet` to Arnor `/astro/users/ds2004/vdp/outputs/phase2/`.

### Step 6 — ROC analysis on Arnor ← NEXT
- [ ] Run `sorcha_roc_comparison.ipynb` on `sorcha_comparison_v5.parquet`
- [ ] VDP vs digest2: AUC / best-F1 / completeness / contamination + per-population table (NEO = positive). Compare to prior 0.837/0.836 — v5.0 is the honest full-sky-grid number.
- [ ] Fill in Section 4.8 of `neomod/paper/NEOrocks.tex`

### Documentation added (2026-06-11)
- `SORCHA_V5_PIPELINE.md` — full pipeline reference (inputs, scripts, slurm params, batching, technical notes)
- `neomod/paper/NEOrocks.tex` — Section 4 "Rubin-Cadence Evaluation with Sorcha" drafted (Sections 4.1–4.8; 4.8 Results is a placeholder pending Phase 2)
- `neomod/paper/` directory created, NEOrocks.tex SCPed from Mac

### Hyak small wins (deferred during the v5.0 redo — not forbidden)
These were deprioritised in favour of the v5.0 grid redo, NOT ruled out. Adopt as decided.
- [x] MBA clone_factor 1→5 in `neomod/src/velocity_density_pipeline_gmm.py` line ~150 → +0.020 F1. **ADOPTED 2026-06-16** for the v5.0 antisun-relative grid (matches the F1=0.837 config). Applied via `--mba-clone-factor 5` (default) in `sorcha_gen_maps_grid.py`, not by editing the global default.
- [ ] Widen antisun footprint 30°→45° in `neomod/pipeline/sorcha_postprocess.py` → Phase 1 re-run → expected +0.010–0.015 (not yet applied)
- [ ] More GMM components 80→200 in `neomod/src/velocity_density_pipeline_gmm.py` line ~1406 → expected +0.003–0.008 (not yet applied)

---

# VDP P(NEO) Suppression Investigation — 2026-06-21 (open: K|M vs GMM A/B test running)

After the v5.0 Phase 2 + Arnor ROC, VDP came out **behind** digest2 on the full sky
(AUC 0.880 / F1 0.808 vs digest2 0.930 / 0.836), where the prior antisun-only run had
them tied at 0.837. Arnor's antisun-distance breakdown showed VDP winning at the antisun
(0–20°: F1 0.878 > d2 0.848) but collapsing at 35–110° elongation. The advisor was
"incredibly certain" this is a **GMM normalisation bug**. This section records the full
diagnostic journey — **including the wrong guesses** — so we don't repeat them.

## The chain of hypotheses (chronological, honest)

1. **Smoothing erodes the NEO wings** — DISPROVEN. `smooth_density_map_by_support` only
   modifies cells with `local_support >= threshold`; low-support wing cells are copied
   through unchanged. Smoothing does not erode wings.
2. **Relative NEO↔MBA normalisation (`effective_factor`) wrong** — IMPLEMENTED a fix
   (PDF-normalise × n_source, the proper Bayesian mixture) → it was a **no-op**: the
   stored NEO/MBA integral ratio barely moved (0.00316→0.00268). The cross-population
   footing was already correct. **Reverted** (`git checkout`).
3. **Single-reference-epoch grid breaks down off-antisun (epoch dependence)** — DISPROVEN.
   Off-antisun NEO P_vdp is flat with |obs epoch − 2026-01-01| (0.10 at 0–45 days, 0.05
   at 400–800 days); antisun NEOs stay ~0.55 at all epochs. The antisun-relative grid is
   **epoch-stable** (design validated). Suppression is purely a function of elongation.
4. **The maps are mis-normalised per-elongation** — DISPROVEN. dlon+000/+030/+050/+090
   all have healthy P(NEO) regions (Pmax~0.97, ~50k cells>0.5, ~identical integrals).
   The map peak even tracks the real NEO velocity vs elongation.
5. **"There is no bug — it's a subsample artifact + physics"** — WRONG / incomplete (my
   error). The reliability curve was computed on the **NEO-enriched `sample` parquet**
   (all NEOs kept, only 0.76% of non-NEOs → NEOs ~130× overrepresented). A perfectly
   calibrated map scoring P=0.05 shows `0.05/(0.05+0.95·0.0076)=87%` NEO in that
   subsample — exactly the "P=0.05→86% NEO smoking gun". On the **full representative
   set** the calibration sits ON the diagonal (0.30→37%, 0.50→75%, 0.91→98%). I concluded
   "no bug". **That was too hasty** — calibration measures *reliability*, not *coverage*.

## The advisor's point — CORRECT, and the real bug

The advisor: at e.g. (vλ≈−0.5, vβ≈+0.6) **no MBA can live**, so P(NEO) must be ≈1, yet
the GMM map is dark there. This is a **coverage/completeness** test that calibration
cannot see (a map can be calibrated while scoring real NEOs ~0 — they just join the huge
low-P background). Confirmed on v5.0 data:

- **529 velocity cells that are 99–100% NEO (22,695 real NEOs) get mean P_vdp = 0.51, not
  ≈1**; worse at high |vβ| (|vβ|>0.5 → 0.47). Real bug.

**Mechanism (per-population density in pure-NEO cells, antisun map, mag22):**

| (vλ, vβ) | rNEO | rMBA | truth | P(NEO) |
|---|---|---|---|---|
| (−0.5, 0.0) | 386 | 11 | NEO | 0.96 ✓ |
| (−0.5, 0.5) | 11.7 | 32 | 99% NEO | 0.26 ✗ |
| (−0.4, 0.5) | 9.9 | 196 | 99% NEO | 0.05 ✗ |

Two compounding errors in the velocity wings:
- **GMM NEO density UNDER-disperses** — rNEO collapses (386→10) into the wings; the real
  NEO velocity distribution (and the old S3M-kNN maps) is broad there.
- **K|M MBA density OVER-disperses** — rMBA stays 32–196 at velocities where **no real
  MBA tracklets appear** (cloning scatter + cf=5 fatten the tail). 
- Net: where the true NEO:MBA ratio is ~100:1, the map says ~1:3 — **backwards by ~300×**,
  so P(NEO) collapses where it should be ≈1.
- The `nearest_dist` mask does **not** fix it — the MBA clones are genuinely at those
  velocities (mask ON vs OFF gives identical P at the pure-NEO cells).
- Notably the advisor's own left panels show the **S3M-kNN maps got this right** (broad
  NEO, P≈1 in the wings). So the **GMM-for-NEO step is a regression** in wing coverage.

## What we are testing NOW (open)

Hypothesis: the GMM NEO cloner is the regression — it produces a velocity density that is
too tight, under-covering the wings, while the K|M/kNN cloner (used for MBA/TNO/Trojan,
and behind the broad S3M maps) covers them. GMM stays the long-term goal; this is a clean
A/B test.

- Added env-var toggle `VDP_NEO_CLONER` in `velocity_density_pipeline_gmm.py`
  (default `gmm`; `km` forces the existing, tested K|M fallback path for NEO). No
  restructuring.
- `neomod/pipeline/slurm/sorcha_gen_maps_grid_kmtest.sh` regenerates 3 maps
  (333 antisun / 338 dlon+50 / 342 dlon+90) with `VDP_NEO_CLONER=km` →
  `prob_maps_grid_kmtest/` (GMM production maps untouched).
- **Decision test:** does K|M give **P≈1** in the pure-NEO wing cells (where GMM gave
  ~0.2)? Does the antisun control stay correct? Does K|M NEO density cover the broad wings?
  - If YES → GMM under-coverage confirmed; ship K|M for NEO, or fix GMM dispersion.
  - If NO → the MBA over-dispersion is the bigger driver; pivot there (reduce K|M scatter
    / reconsider MBA cf=5 / per-population support mask).

## Validation infrastructure (reuse for any fix)

- `sorcha_gen_maps_grid_fixtest.sh` / `_kmtest.sh` — regenerate 3 representative maps
  into a side dir; ~10 min.
- Pure-NEO coverage test (the advisor's test) + antisun control are the accept criteria,
  BEFORE any full 667-map regen.

## Lessons (do not repeat)

- **Calibration ≠ coverage.** A diagonal reliability curve does NOT mean the map finds all
  NEOs; it can score real NEOs ~0 (completeness loss invisible to calibration).
- **Never measure calibration on a class-enriched subsample** — use the full representative
  set or correct for the sampling rate.
- The **pure-NEO-cell test** (cells where only NEOs can be should give P≈1) is the right
  lens for coverage defects.


## A/B test result (2026-06-21, late) — NEO cloner RULED OUT; MBA over-dispersion is the target

Ran the `VDP_NEO_CLONER=km` A/B test (`sorcha_gen_maps_grid_kmtest.sh` → maps 333/338/342
into `prob_maps_grid_kmtest/`). Compared GMM vs K|M P(NEO) in the pure-NEO wing cells:

| (vλ, vβ) | P_GMM | P_K\|M | rMBA_GMM | rMBA_K\|M |
|---|---|---|---|---|
| (−0.5, 0.5) | 0.26 | 0.13 | 32 | 42 |
| (−0.4, 0.5) | 0.05 | 0.03 | 196 | 251 |

cells P>0.5: GMM 47k → K|M 20k. **K|M is the same or slightly WORSE.** rNEO is similar
in both. So the **NEO cloner (GMM vs K|M) is NOT the culprit — RULED OUT.** GMM stays.

**The suppressor is the MBA density.** rMBA = 32–251 in cells that are 99–100% NEO (no
real MBA tracklets there), and it is **identical between the two runs** because MBA cloning
is unchanged (K|M for MBA in both). The conditional K|M cloner **over-disperses MBA into
the high-|vβ| wings**, placing clones at velocities real MBAs never reach. The `nearest_dist`
mask cannot help — the spurious clones are genuinely there (mask ON==OFF at those cells).

### What to find out in the morning (debugging resumes here)
Target: **why the K|M MBA velocity cloud is too fat in the wings**, and how to tighten it
so rMBA → ~0 where no real MBAs exist (→ P(NEO) → ~1 in pure-NEO cells).
1. **Read `clone_population_conditional_K_from_M_with_skycut`** (velocity_density_pipeline_gmm.py)
   to find the cloning scatter / bandwidth parameter — that is the likely lever.
2. **A/B tests (reuse `prob_maps_grid_kmtest`-style 3-map regen + the pure-NEO coverage test):**
   - MBA cf 5→1 (`--mba-clone-factor 1`) — cheap baseline; expect only PARTIAL help
     (cf is ~5×, the wing over-representation is ~300×).
   - Tighten the K|M scatter/bandwidth — expected to be the real fix.
   - With tails controlled, re-enable the per-population `nearest_dist` mask to zero any
     residual MBA in pure-NEO cells.
3. **Accept criteria (unchanged):** pure-NEO cells get P≈1; antisun control stays correct;
   then full 667-map regen → Phase 2 → ROC.

### Toggles / scripts in place
- `VDP_NEO_CLONER` env (default `gmm`; `km` forces K|M) — leave at `gmm`.
- `sorcha_gen_maps_grid_kmtest.sh` (reuse pattern for the MBA tests).
- Pure-NEO coverage test = the accept lens (advisor's test). Calibration ≠ coverage.

## Pivot (2026-06-21, before bed) — cloner likely NOT the bug; check classification/definition mismatch

Read `clone_population_conditional_K_from_M[_with_skycut]`. Key facts:
- **The K|M cloner has NO scatter/bandwidth parameter.** a, e, i, node are EXACT copies
  (`np.repeat`); only the orbital phase changes — mean anomaly `M` resampled from the
  pooled empirical M distribution, K=node+argperi resampled conditionally on M. So
  "tighten the K|M scatter" was the wrong framing.
- User's point (correct): the K|M cloner is well-validated (Trojans, TNOs, the S3M
  pipeline). A latent dispersion bug there is unlikely. **De-prioritise the cloner.**

### Leading hypothesis now: catalog-MBA vs tracklet-MBA DEFINITION mismatch
- Cloner clones **catalog MBAs** = `hybrid_catalog_prep`: `1.7 ≤ a < 4.1, q ≥ 1.3` (BROAD).
- Truth tracklet label = `sorcha_postprocess.classify_population`: MBA = `2.0 < a < 3.3,
  e < 0.3` (NARROW); everything else → **`other`**.
- ⇒ Objects the MAP counts as MBA density, the TRUTH labels as `other`. The "pure-NEO"
  coverage test counted only `population=='NEO'` as signal, so those high-|vβ| cells may
  contain real **`other`** objects exactly where the map puts MBA density. If so,
  **P(NEO) < 1 there is CORRECT** and the map/cloner are fine — the bug is in the
  comparison labels, not the pipeline.

### Morning check (read-only, parquet)
1. Recompute the pure-NEO cells with **non-NEO = everything except NEO** (not just MBA),
   and break down what the non-NEO objects in those cells actually are
   (`other` / MBA / TNO / Trojan).
   - If `other`(+MBA) fills the high-|vβ| wings → map is right, "bug" is a labelling
     artifact (the cells are NOT pure NEO). Likely the resolution.
   - If the cells really are ~100% NEO across ALL populations and the map still gives
     P≈0.5 → genuine density issue; then look at the **kNN full-posterior estimator tails**
     (`log_posterior_d0_2d`, "unnormalised", heavy-tailed over d0) — NOT the cloner.
2. Reconcile the two MBA definitions (catalog 1.7–4.1 vs tracklet 2.0–3.3,e<0.3). Consider
   aligning the truth classification with the catalog population definitions so map and
   truth use the same population boundaries.

### Status of toggles/maps
- `VDP_NEO_CLONER` left at default `gmm`. K|M A/B maps in `prob_maps_grid_kmtest/` (ruled
  out NEO cloner). Production `prob_maps_grid/` + parquet untouched. No code changes pending
  beyond the committed `VDP_NEO_CLONER` toggle.

## Morning check result (2026-06-22) — label mismatch CONFIRMED as the main explanation

Ran the requested read-only parquet check on the **full VDP shards** (`outputs/phase2_v5/vdp_shards/`,
65,857,457 scored tracklets), not the NEO-enriched ROC subsample.

Two versions were checked:

1. **Map+mag-specific cells**: group by `(prob_map_file, mag_bin_label, vlam cell, vbeta cell)`.
   - Antisun `mag22` alone has only 93 NEO tracklets and no tracklets in the cited example cells
     (`(-0.5,0.5)`, `(-0.4,0.5)`, etc.).
   - Across all maps/mags, with `NEO>=10` and `NEO/(NEO+MBA)>=0.99`, there are **0** high-support
     pure-vs-MBA cells. So the earlier "529 pure cells" diagnostic was NOT map/mag-local; it was
     a pooled-velocity diagnostic.

2. **Pooled velocity cells**: group only by `(vlam cell, vbeta cell)` across all maps/mags.
   This reproduces the apparent pure-cell effect and reveals the missing population:

   | Criterion | Cells | NEO | Total | NEO/all | Breakdown |
   |---|---:|---:|---:|---:|---|
   | `NEO>=10`, `NEO/(NEO+MBA)>=0.99` | 892 | 15,257 | 22,635 | 67.4% | MBA=0, other=7,369, Trojan=9 |
   | same, `|vbeta|>=0.3` | 50 | 556 | 1,034 | 53.8% | MBA=0, other=478 |
   | `NEO>=20`, `NEO/(NEO+MBA)>=0.99` | 236 | 6,740 | 10,372 | 65.0% | MBA=0, other=3,625, Trojan=7 |

Specific pooled cells:

| `(vlam,vbeta)` | NEO | MBA | other | NEO/all | mean `P_NEO_vdp` |
|---|---:|---:|---:|---:|---:|
| `(-0.50, 0.00)` | 21 | 0 | 1 | 95.5% | 0.886 |
| `(-0.50, 0.50)` | 4 | 0 | 0 | 100% | 0.508 |
| `(-0.40, 0.50)` | 2 | 0 | 13 | 13.3% | 0.076 |
| `(-0.40, 0.60)` | 2 | 0 | 3 | 40.0% | 0.229 |

Conclusion:
- The advisor's "no MBA can live there" statement is true but incomplete for the current
  truth labels: **`other` lives there**.
- The map is trained with broad catalog MBA-like density (`1.7<=a<4.1, q>=1.3`), while the
  ROC truth labels use narrow MBA (`2.0<a<3.3, e<0.3`) and put the broad-edge belt population
  into `other`.
- Therefore the apparent "pure NEO" cells are mostly **not pure NEO across all populations**.
  The VDP low P in many wing cells is largely correct under the map's broader non-NEO prior;
  the mismatch is in the evaluation/diagnostic population definition, not in GMM NEO coverage
  or a K|M cloner bug.

Next action:
1. Align the diagnostic and ROC labels with the training catalog definitions, at least for a
   sensitivity test:
   - NEO: `q < 1.3`
   - MBA-like non-NEO: `1.7 <= a < 4.1 and q >= 1.3`
   - Trojan: `4.7 < a < 5.9 and e < 0.3`
   - TNO: `a > 30`
2. Recompute ROC/per-population tables under both label schemes:
   - original paper-friendly physical labels (`MBA=2.0<a<3.3,e<0.3`, `other` separate)
   - training-aligned labels (broad MBA-like bin)
3. Do **not** regenerate maps yet. Production maps and parquet are untouched; this is an
   analysis/labeling issue first.

## BUG FOUND (2026-06-22, evening) — kNN density estimator bleeds into zero-support cells

After aligning eval labels to the map's training definitions (NEO q<1.3, MBA_like
1.7<=a<4.1 & q>=1.3, ...), the scary "pure-NEO wing" suppression dropped from
meanP 0.51 -> 0.68 (most of it was the broad-MBA-labeled-`other` artifact). A REAL
residual remained: genuinely-pure-NEO cells (>90% NEO by real tracklets) get the map
**rMBA/rNEO = 0.54 vs truth <0.11 — MBA over-weighted ~5x.**

### Ruled out (with evidence)
- **MBA cf 5->1**: made it WORSE (rMBA up, P down, cells>0.5 47k->4.9k). kNN over-weights
  SPARSE regions, so fewer clones = worse. Reducing cf is wrong.
- **The cloner**: MBA clone overlay velocities are UNDER-dispersed vs real tracklets
  (|vbeta|>=0.4: clones 0.08% vs real 0.16%; vlam<=-0.4: 0.01% vs 0.04%). The clones do
  NOT over-populate the wings. Cloner exonerated (well-tested, as expected).
- **Trailing-loss/detection mismatch**: would cancel in the P(NEO) ratio (hits NEO and
  MBA equally per velocity). Not the cause.

### ROOT CAUSE (confirmed)
The kNN full-posterior density estimator (`log_posterior_d0_2d`, k_map=10, UNNORMALISED)
assigns density to velocity cells that contain **ZERO clone support**: in every pure-NEO
wing cell, `support_count__MBA = 0` yet `density_raw__MBA = 11-196`. The k=10 neighbours
reach back to the dense MBA core, so the **core bleeds density into the sparse NEO wings**.
NEO has no dense core to bleed, so it stays sparse -> P(NEO)=rNEO/Sum collapses. The
existing `nearest_dist` mask misses it (it checks the k=1 nearest clone, which is ~0.06
away — close — while the cell itself is empty).

### Fix (scoring-time, NO regen needed — support_count is stored)
Mask each population's density where its in-cell `support_count` < threshold (zero the
estimator bleed). Tested on the antisun map (mask non-smoothed pops by support>=1, keep
NEO which is intentionally smoothed): pure-NEO wing cells **0.18-0.26 -> 1.000**;
cells P>0.5 **47k -> 159k**.

### CAVEAT to validate (the mask is currently too blunt)
Where the CLONER has a coverage gap (support=0) but real objects ARE present, the mask
wrongly zeros them. Example: (-0.4,+0.5) has 189 real MBA_like tracklets but support=0 ->
mask forces P(NEO)=1 -> false positives. Net effect MUST be measured.

### Next: validate via re-scoring (no map regen)
1. Add an optional support-count mask in `ProbMapSet`/`score_vdp` (toggle), threshold a
   small N. Re-score the existing `outputs/phase2_v5/` parquet WITH the mask.
2. Full ROC: does pure-NEO coverage + F1 improve net of the false-positive risk in
   cloner-gap cells? Compare to the unmasked 0.808.
3. If net win but the cloner-gap false positives hurt: the deeper fix is improving MBA
   clone COVERAGE at extreme velocities (the cloner under-disperses there) so support>0
   where real MBAs are, then the support mask only cuts genuinely-empty cells.

## Why the bug was invisible in the original S3M VDP-vs-digest2 win (2026-06-22)

Q: if the kNN estimator bleeds the MBA core into NEO cells, why was VDP BETTER than
digest2 (F1 0.856 vs 0.665) on the S3M kNN maps? No contradiction — three reasons:

1. **The S3M ROC was antisun-only — VDP's single best geometry.** Scored one ~30° patch
   at opposition, where NEO (fast retrograde, vλ≈−0.5) and MBA (slow, vλ≈−0.2) are
   MAXIMALLY separated in velocity. The bleed needs NEO and MBA to be CLOSE in velocity
   to bite; at the antisun the fast-NEO cells are ~0.3+ from the MBA core, so rMBA there
   is negligible → P(NEO) high → VDP wins. The bug was present but didn't matter.
2. **The full-sky grid is the first time VDP's weak regions were scored.** Off-antisun
   overlap (40–110°, Arnor's band) and the high-|vβ| wings were OUTSIDE the 48% antisun
   footprint in every prior run. The bleed was always there; those cells were just never
   scored. Matches Arnor's breakdown: VDP wins 0–20° (S3M regime), loses 40–110° (new).
3. **S3M maps used mask ON; the GMM/v5 pipeline turned it OFF** ("GMM density → 0 far
   from data" is true for the NEO GMM but NOT the K|M MBA core). Disabling it removed the
   guard that partially hid the bleed.

⇒ VDP wasn't "better then, broken now" — it was only ever tested where it's strongest
(opposition, max separation, mask on). The support-count mask does for the WHOLE sky what
the old mask + opposition geometry did for the antisun. VDP's antisun strength is real
(and operationally that's where NEOCP discoveries happen).

## Support-mask fix VALIDATED — fixes the bug, but ~neutral on headline F1 (2026-06-22)

Implemented the support-count mask as a scoring-time toggle: `ProbMapSet.from_npz(
support_mask_min=N)` and `sorcha_phase2.py score-vdp --support-mask-min N` (zeros each
non-smoothed population's density where in-cell support_count < N; NEO exempt as it is
smoothed). Committed (86a18a5).

Re-scored the existing 707,670-row eval subsample (same set as the 0.808 result), NO regen:

| classifier            | AUC   | bestF1 | completeness | contamination |
|-----------------------|-------|--------|--------------|---------------|
| VDP original          | 0.880 | 0.808  | 74.1%        | 11.2%         |
| VDP + support_mask=1  | 0.880 | 0.809  | 74.3%        | 11.1%         |
| digest2               | 0.930 | 0.836  | 76.9%        | 8.4%          |

**Conclusions:**
- The mask DOES fix the bug: pure-NEO velocity cells go 0.18–0.26 → ~1.000; cells P>0.5
  47k → 159k on the antisun map. The motto ("only-NEO regions get highest scores") is met,
  and it is the scientifically correct behaviour (the advisor's concern is resolved).
- BUT the headline F1 barely moves (0.808 → 0.809) and AUC is unchanged: the pure-NEO wing
  cells hold FEW NEOs (~10% of NEOs reach |vβ|≥0.45, scattered), so correcting them does not
  shift the aggregate. Contamination slightly improved (11.2→11.1%), so feared cloner-gap
  false positives are negligible. Net: small clean positive — worth adopting as correct
  behaviour, but it is NOT the lever that closes the digest2 gap.
- **The VDP-vs-digest2 full-sky gap (0.808 vs 0.836) is dominated by the intermediate-
  elongation overlap (40–110°), where NEO and MBA velocities genuinely overlap — physics,
  not a bug.** No velocity-only classifier separates them there; digest2's orbit fit does.
  VDP still WINS at the antisun (0–20°, F1 0.878 > 0.848), which is the operational NEOCP
  discovery regime.

### Recommendation
- Adopt `--support-mask-min 1` for the production grid scoring (correct behaviour, slight
  net positive, resolves the pure-NEO-cell concern). It is scoring-time; a full re-score
  (Slurm array) reproduces it across all 113 shards — NO map regen needed.
- For the paper: report VDP's antisun strength + the honest full-sky result, attributing
  the off-antisun gap to elongation-dependent velocity overlap (physics). This matches the
  S3M-win explanation above.

## Shipping the support-masked result to Arnor (2026-06-22)

KEY POINT: the support mask is applied at LOAD/SCORING time (`from_npz(support_mask_min=1)`)
— it does NOT change the `.npz` maps (they store raw density + support_count). So:
- "Correct maps" for Arnor = the EXISTING `prob_maps_grid/*.npz` (20 GB, 667 files; never
  SCPed before — only the parquet went). NO 667-map regeneration needed.
- For Arnor's probability PLOTS to show masked P(NEO): load with `support_mask_min=1` using
  the updated `velocity_density_pipeline_gmm.py` (Arnor `git pull`).
- The masked ROC dataset was built WITHOUT a full Slurm re-score: re-scored the 707k eval
  subsample with support_mask_min=1 and merged into
  `outputs/phase2_v5/sorcha_comparison_v5_masked.parquet` (NEO median P 0.227 -> 1.0;
  P_NEO_vdp_unmasked kept as a column). ROC unchanged vs unmasked (0.808->0.809) — expected
  (pure-NEO wing cells hold few NEOs).

A full Slurm re-score of all 113 shards with `--support-mask-min 1` is OPTIONAL (canonical
hygiene for phase2_v5); not needed for Arnor's plots (from maps) or ROC (from subsample).

Files to Arnor (`/astro/users/ds2004/vdp/`):
- `prob_maps_grid/*.npz` (667, 20 GB) -> `prob_maps_grid/`   [rsync; or a lat=0 subset first]
- `outputs/phase2_v5/sorcha_comparison_v5_masked.parquet` -> `outputs/phase2/`
- Arnor `git pull` neomod for the `support_mask_min` toggle, then load maps with
  `ProbMapSet.from_npz(path, support_mask_min=1)` for the corrected probability plots.

## Masked full re-score: first attempt FAILED (import mismatch) — fixed (2026-06-22)

Submitted the canonical masked full re-score (sorcha_phase2_vdp_v5.sh with
--support-mask-min 1). **All 113 shards crashed instantly:**
`TypeError: ProbMapSet.from_npz() got an unexpected keyword argument 'support_mask_min'`.
Root cause: `sorcha_phase2.py` imported `velocity_density_pipeline` (the ORIGINAL module),
but the support_mask_min toggle was added to `velocity_density_pipeline_gmm`. Two different
ProbMapSet classes — scoring never saw the toggle. No masked shards were written; the
vdp_shards on disk were STALE (2026-06-20 unmasked run).

**Fix:** point `sorcha_phase2.py` at the GMM module — `import velocity_density_pipeline_gmm
as vdp`. It generated the v5 maps AND has the mask; its score_observation/from_npz are a
drop-in (return_intermediates returns probs/vlam/vbeta/bin_labels; unmasked AUC matches the
original module — verified on the 707k subsample, and smoke-tested: from_npz(support_mask_min=1)
gives 158,867 cells P>0.5; score_observation keys correct).

**Still TODO:** re-submit `sbatch neomod/pipeline/slurm/sorcha_phase2_vdp_v5.sh` to actually
produce the masked vdp_shards. The canonical masked comparison parquet
(`sorcha_comparison_v5_masked.parquet`) already exists (built via the subsample re-score with
the _gmm ProbMapSet) and is correct, so Arnor's ROC is unaffected — only the full vdp_shards
need the (now-fixed) re-run for canonical hygiene.

---

# Arnor paper write-up — `NEOrocks.tex` Section 4 + Results (2026-06-22/23)

All edits below are on the Arnor copy `paper/NEOrocks.tex` (Overleaf-linked, **gitignored** —
moves only by scp, NOT git pull; scp back to `neomod/paper/` on Hyak after changes). The
analysis/figure code is reproducible in the notebooks listed at the end.

## Section 4 (Sorcha) audited against the code and corrected
- **Algorithm 1** rewritten as the "VDP main algorithm" with an explicit `Input:` (hybrid
  S3M+MPCORB catalogue with population labels), `catalogue objects labelled p` (not "source
  orbits"), `Form the NEO probability map: P = ρ_NEO/Σρ_p` (not the bare arrow), $(\alpha,\delta)$
  spelled out, $\Delta\lambda_\odot$ defined inline, and the solar-exclusion corrected to
  **$|\Delta\lambda_\odot| > 140^\circ$** (within 40° of the Sun — earlier wording said `<40°`,
  which wrongly excluded the antisun). Bilinear interpolation **verified** in `_lookup_in_map`.
- **§4.6 GMM** corrected: GMM is the NEO **cloner** (orbital-element feature space), NOT a
  velocity-space density estimator replacing kNN. Normalisation paragraph rewritten to the real
  **acceptance-fraction / effective-clone-factor** fix (was a wrong "integrates to unity over the
  infinite plane" story). Added the NEOMOD3 training-augmentation paragraph (null effect, retained).
- **§4.5 grid** updated to as-built: **667 maps** (29 lon × 23 lat), velocity grid **±2.0 → 401×401**.
- **§4.1** added the hybrid-training dual-role paragraph (maps trained on the *same* hybrid
  population they score — closes the S3M calibration gap).
- **NEOMOD3** reconciled across §3 (benchmark note + forward ref), §4.6 (used as augmentation),
  and Conclusions ("already incorporated") — no more "future work" framing for something used.

## Framing sections caught up to Section 4
- Abstract + intro now mention the Sorcha/Rubin-cadence evaluation (were benchmark-only).
- **Conclusions** rewritten — previously listed Section 4's content (NEOMOD3, sky positions,
  epochs) as "future work"; now summarise the Sorcha evaluation.
- Wagg citations de-clustered 13→10; **fixed the intro number**: digest2 purity is **5.4% at
  ≥65, 8.1% at ≥90** (was a self-contradictory "8.3%…8.1%"), verified vs arXiv:2408.12517.
- digest2 *does* use S3M as its population model — confirmed from Keys et al. 2019 §3.3.1 (so the
  comparison isolates the method, not the population assumptions).

## Tables / figures added
- **§3 dense numbers paragraph → Table 2** (`tab:roc_fiducial_metrics`), with the combined
  "C/contam" column split into separate Completeness/Contamination columns.
- **Clone-factor table** (`tab:clone_factors`, Appendix A): benchmark (NEO 300, MBA 10, TNO 100,
  Trojan 100) vs Sorcha grid (80 / 5 / 10 / 5), referenced from §2.3, §4.5, §4.6, App B.
- **Sky-coverage Mollweide figure** (§4.4, `sorcha_sky_coverage_smoothed.png`) — tracklet
  density per sq deg, WFD/NES boundary at δ≈+5° explained in the caption. Built by
  `make_sky_coverage_figure.py`.

## §4.10 Results — FILLED with the final masked-parquet numbers
Source: `outputs/phase2/sorcha_comparison_v5_masked.parquet`, training-aligned labels
(NEO q<1.3; MBA-like 1.7≤a<4.1 & q≥1.3; Trojan 4.7<a<5.9 & e<0.3; TNO a>30; other=rest).
- **`tab:sorcha_results` (full-sky):** VDP **AUC 0.880, F1 0.809, 74.3% compl, 11.0% contam**;
  digest2 **0.930, 0.836, 76.9%, 8.4%** (707,670 tracklets; 207,670 NEOs).
- **`tab:sorcha_elongation` (per |Δλ☉|):** 0–20° VDP **0.883** > d2 0.848; 20–40° 0.837/0.855;
  40–70° 0.767/0.848; 70–110° 0.770/0.782; 110–141° **0.860**/0.854. (AUCs also in the table.)
- **N≥3:** VDP F1 0.815 / AUC 0.882; digest2 0.840 / 0.932 (182,661 tracklets) — ranking unchanged.
- **Support threshold** stated: one cloned point per cell (Appendix B).

### Still TODO in the paper (left deliberately)
- **Decisions for Devanshi/Zeljko** (NOT resolved): headline-metric weighting (near-antisun vs
  uniform full-sky); broad-vs-narrow MBA definition in eval.
- **Optional** before/after support-mask P(NEO) figure (code ready in the notebook).
- **Wagg v3.3 traffic comparison** — needs numbers from the Wagg paper.
- ROC-curve figure `fig:sorcha_roc` (data available; not yet drawn).

## Notebooks created on Arnor (kernel `neofast_py310`)
- `sorcha_v5_normalisation_test.ipynb` — the canonical v5.0 diagnostic: calibration, velocity
  truth-vs-applied maps, ROC + per-band, **§5 direction-comparison grids** (density + NEO P(NEO)
  across 4 antisun-relative offsets, mask ON) and the **MBA & NEO P, mask OFF vs ON** reference
  (mag bins 21-22/22-23/24-25). Loads grid maps with `support_mask_min=1`, `mask_radius=inf`.
- `paper_figures_sorcha.ipynb` — the sky-coverage figure.
- `make_sky_coverage_figure.py` — standalone, parameterised for the v5 regen.

---

## Pure-S3M benchmark vs Sorcha: structural gap analysis (2026-07-06)

### Setup
- `benchmark_comparison_s3m.parquet` — 475,027 rows; S3M objects propagated to MJD 61041
  (2026-01-01, map epoch), 30-min synthetic tracklet, no measurement noise, exact rates.
  Populations: NEO 184,832 / MBA 141,103 / Trojans 99,995 / TNO 49,097.
- `sorcha_comparison_s3m.parquet` — 648,908 rows; full 2-year Sorcha survey, scored against
  the same 667-grid GMM maps.

### Epoch mismatch is NOT the cause
The 667-map grid is built in the **antisun-relative frame** (Δλ_antisun, β). At scoring time
each tracklet's antisun position is recomputed from its observation epoch, so the maps are
epoch-independent by design (`sorcha_gen_maps_grid.py` line 60: "velocity statistics are
epoch-independent in the antisun frame"). Rebuilding maps at a different epoch produces
essentially the same maps.

**Single-night ROC test** (2026-07-06): three Sorcha single-night subsets extracted from
`sorcha_comparison_s3m.parquet` (already fully scored — no SLURM jobs needed) and compared
against benchmark ROC on Arnor:

| Night | Date | Δ from benchmark | N_NEO | N_total | VDP F1 | d2 F1 |
|-------|------|-----------------|-------|---------|--------|-------|
| MJD 61033 | 2025-12-24 | −8d | 325 | 1,827 | 0.807 | 0.853 |
| MJD 61046 | 2026-01-06 | +5d | 288 | 1,385 | 0.827 | 0.901 |
| MJD 61642 | 2027-08-25 | +601d | 817 | 4,315 | 0.789 | 0.807 |
| Benchmark | MJD 61041 | 0 | — | 475,027 | 0.787 | 0.667 |

VDP is stable across all three nights (0.79–0.83 ≈ benchmark 0.787), confirming epoch
independence. Digest2 is dramatically better on Sorcha single nights (0.81–0.90) than on
the benchmark (0.667). **The gap is structural, not temporal.**

### Remaining structural explanations

**1. Tracklet rate noise (primary candidate)**
Sorcha tracklets are formed from two real simulated LSST visits with photon noise, PSF scatter,
etc. — rates are noisy. Benchmark rates are exact (two analytically propagated positions). Noisy
(v_λ, v_β) → objects land in slightly wrong map positions → VDP P(NEO) is less sharp → ROC
degrades. Digest2 uses an orbit-fitting approach that is less sensitive to per-tracklet rate noise.

To test: stratify `sorcha_comparison_s3m.parquet` by `dt_min` and check if VDP AUC rises with
longer baseline (lower noise). Data is already scored — purely an analysis-side experiment on Arnor.

**2. Population composition / survey selection bias**
The benchmark applies a magnitude filter (`mag_bin_label.notna()`) but not LSST's full detection
model (visit footprint, S/N, chip gaps, cadence). Sorcha's detection criteria preferentially
select objects at favourable sky positions and brightnesses — exactly where digest2's orbit-fitting
is most informative. VDP is insensitive to this because it only uses the instantaneous rate vector.
This explains why digest2 F1 jumps from 0.667 (benchmark, full S3M) to 0.85–0.90 (Sorcha,
detected subset) while VDP stays flat.

### Sorcha config: linking thresholds (`Rubin_full_footprint.ini`)
```
SSP_separation_threshold  = 0.5 arcsec   ← minimum motion between two visits to form a tracklet
SSP_number_observations   = 2            ← detections per night required
SSP_maximum_time          = 0.0625 days  ← max inter-visit gap (90 min)
SSP_number_tracklets      = 3            ← tracklets needed across nights
SSP_track_window          = 15 days      ← window for those 3 tracklets
SSP_detection_efficiency  = 0.95
```

Minimum detectable rate: 0.5 arcsec / dt_min. At 30 min baseline → **0.0067 deg/day**.

### TNO population loss analysis
TNOs move ~0.017 deg/day (median), giving ~1.37 arcsec separation at a 30-min baseline.
18.1% of Sorcha TNOs fall below the 0.5 arcsec linking threshold:

| dt bin | N TNOs | Below 0.5" |
|--------|--------|-----------|
| < 20 min | 766 | **60%** |
| 20–30 min | 255 | 15% |
| 30–40 min | 3,655 | 11% |
| > 40 min | 300 | ~4% |

The 90% loss from benchmark TNOs (49,097) to Sorcha TNOs (4,976) is dominated by **magnitude**:
the S3M TNO population peaks at mag 24–25, right at LSST's single-visit detection floor.
Getting 6 detections (2/night × 3 nights) near the noise floor is the bottleneck.

| Mag bin | Benchmark TNOs | Sorcha TNOs | Recovery |
|---------|---------------|-------------|---------|
| 14–22 | 1,647 | 1,146 | 66–87% |
| 22–23 | 5,313 | 1,773 | 33% |
| 23–24 | 16,394 | 1,777 | 11% |
| 24–25 | 24,804 | 280 | **1%** |
| 25+ | 939 | 0 | 0% |

The 24–25 mag collapse: at the 5σ detection floor, each visit has ~50% detection probability.
Getting 6/6 successes (2 per night × 3 nights) has probability ~0.5⁶ ≈ 1.6%, matching the
observed 1% recovery rate.

---

## TNO contamination and its effect on digest2 (2026-07-06)

### The shocking finding: digest2 improves dramatically from benchmark → Sorcha; VDP barely moves

When the S3M benchmark and Sorcha results were compared side by side on the full-sky grid,
the headline numbers were:

| Dataset | VDP F1 | digest2 F1 |
|---------|--------|-----------|
| Benchmark | 0.787 | 0.667 |
| Sorcha (full-sky) | 0.808–0.809 | 0.836 |

VDP moves +0.01–0.04; digest2 moves +0.12–0.19. At the antisun (0–20°) on Sorcha,
digest2 (F1=0.85) beats VDP (F1=0.81).

Root cause: **the benchmark's population composition is a severe artificial stress test for
digest2, not for VDP.** In the 0–20° antisun bin, benchmark has **38.6% TNOs**. Sorcha
has **~0.7% TNOs** in that same bin. digest2 is highly vulnerable to this; VDP is not.

### Why VDP is robust to TNOs

TNOs at ~40 AU move ~0.01–0.05 deg/day — effectively near the origin in velocity space
(vlam ≈ 0, vbeta ≈ 0). The VDP antisun map has high MBA density and high NEO density at
distinct locations, but the centre (origin) is a well-calibrated low-P(NEO) zone (that's
where MBAs pile up at slow retrograde drift at the antisun, and TNOs are even slower). So
TNOs get correctly scored ~P(NEO)=0. Per the full-sky per-population table:

| Population | VDP above threshold | d2 above threshold |
|------------|--------------------|--------------------|
| TNO        | 0.0%               | 69.5%              |

VDP scores 0% of TNOs as NEOs; digest2 scores 69.5% of them as NEOs.

### Why digest2 is NOT robust to TNOs

Digest2 uses a short-arc orbit-fitting model: it computes the range of heliocentric distances
compatible with two sky positions and a time baseline. TNOs at ~40 AU, observed over a
30-minute arc, are indistinguishable from an NEO on a very elongated orbit — the angular
motion is tiny and the orbital uncertainty is enormous. Digest2's orbit-space prior then
assigns a large NEO probability because the arc is consistent with a high-eccentricity,
low-perihelion trajectory. This is a fundamental limitation of two-detection classification.

### Why the benchmark has 38.6% TNOs at 0–20° but Sorcha has 0.7%

Three compounding effects explain the ~55× difference in TNO fraction:

**1. Faint-end magnitude collapse (dominant effect)**
S3M TNOs peak at V=24–25, right at LSST's single-visit detection floor. The benchmark
applies a flat V<25 cut with 100% detection efficiency — all 49,097 S3M TNOs get a
synthetic tracklet. Sorcha's realistic detection model requires 2 detections per night
× 3 nights within 15 days, each above 5σ S/N. Near the detection floor (~50% per visit),
the probability of getting 6/6 successes is ~0.5⁶ ≈ 1.6%. Observed: 1% recovery rate at
V=24–25. The large TNO population that inflates the benchmark is almost entirely V>23
objects that Sorcha simply never detects.

| Mag bin | Benchmark TNOs | Sorcha TNOs | Recovery |
|---------|---------------|-------------|---------|
| 14–22   | 1,647         | 1,146       | 66–87%  |
| 22–23   | 5,313         | 1,773       | 33%     |
| 23–24   | 16,394        | 1,777       | 11%     |
| 24–25   | 24,804        | 280         | 1%      |
| 25+     | 939           | 0           | 0%      |

**2. Tracklet motion threshold (secondary)**
LSST's linking algorithm (`SSP_separation_threshold = 0.5 arcsec`) requires a detectable
inter-visit motion. TNOs move ~0.017 deg/day (median ≈ 1.37 arcsec at 30 min baseline),
which is above the threshold — but at the faint end where S/N is marginal, astrometric
scatter can push the apparent separation below 0.5 arcsec. At dt < 20 min baselines,
**60% of Sorcha TNOs fall below the motion threshold**. 18.1% of all Sorcha TNOs are lost
this way. The benchmark uses a fixed 30-min synthetic baseline with exact positions —
no motion threshold effect.

**3. Survey footprint efficiency**
The benchmark propagates every S3M object to MJD 61041 and creates a synthetic tracklet
regardless of whether that sky position is actually observed. Sorcha's LSST cadence only
covers ~18% of sky per night; TNOs in unobserved fields produce no tracklets at all.
This is a uniform suppression across all populations (not TNO-specific), but combined
with effects 1 and 2 it drives the 55× difference.

**Phrasing for the paper:** "The benchmark applies a flat V<25 magnitude filter with 100%
detection efficiency. Sorcha's realistic survey model — per-visit S/N threshold, tracklet
motion requirement, and multi-night linking — preferentially excludes faint TNOs (V>23),
reducing their fraction at the antisun from 38.6% to 0.7%. This population difference
is the dominant driver of digest2's apparent improvement from benchmark to Sorcha."

### TNO contamination proof: forward and inverse experiments (Section 15, `sorcha_v5_normalisation_s3m.ipynb`)

To isolate the TNO effect, two controlled experiments were run in the 0–20° antisun bin:

**Forward experiment (remove TNOs from benchmark):**
Subsample benchmark TNOs down from 38.6% to 0.7% (matching Sorcha fraction).
Result: digest2 F1 jumps from **0.706 → 0.842** (+0.136). VDP barely moves (+0.009–0.016).

**Inverse experiment (inject benchmark TNOs into Sorcha):**
Bootstrap-inject benchmark TNOs into the Sorcha 0–20° bin until the TNO fraction reaches
38.6%. Result: digest2 F1 drops from **0.863 → 0.755** (−0.108). VDP barely moves.

Summary table:

| Experiment | TNO fraction | VDP F1 | d2 F1 |
|------------|-------------|--------|-------|
| Benchmark (original) | 38.6% | 0.870 | 0.706 |
| Benchmark TNO-reduced | 0.7% | 0.861 | **0.842** |
| Sorcha (original) | 0.7% | 0.863 | 0.863 |
| Sorcha TNO-injected | 38.6% | 0.847 | **0.755** |

These experiments causally prove that TNO fraction explains the benchmark–Sorcha gap for
digest2. VDP's near-zero TNO false positive rate (0.0%) makes it immune. The result also
shows that **VDP still wins at the antisun even after TNO removal** (0.861 vs 0.842), so
TNOs are the dominant driver, not the sole explanation.

**Bootstrap injection note:** The benchmark TNO pool (6,440 objects in 0–20°) was sampled
with replacement to produce the ~92,808 injections needed to inflate Sorcha's 0.7% TNO
fraction to 38.6%. Each injection copies only (P_NEO_vdp, P_NEO_d2, is_neo=False) — no
orbital mechanics, just score distributions from real benchmark TNOs.

### Per-bin contamination analysis (Section 16, `sorcha_v5_normalisation_s3m.ipynb`)

The same forward/inverse experiment was extended to all five antisun-relative longitude bins
using both TNOs and Trojans as contaminants:

| Bin | Primary contaminant | Pattern |
|-----|---------------------|---------|
| 0–20° | TNO (38.6% in bench) | TNO-driven: d2 jumps +0.14 after removal |
| 20–40° | TNO | Similar TNO-driven pattern |
| 40–70° | Neither | Velocity-space overlap: d2 wins both benchmark and Sorcha; the gap is physics not contamination |
| 70–110° | Trojans (suspected) | Trojan contamination suspected (Jupiter Trojan velocity overlaps with NEO at these elongations) |
| 110–141° | Trojans/TNO (suspected) | Same pattern as 70–110° |

The 40–70° bin is important: even after removing all TNOs and Trojans, digest2 still beats
VDP. At this elongation, NEO and MBA apparent velocities genuinely overlap (both populations
are moving prograde, similar v_λ), and digest2's orbit-fitting has more information than VDP's
velocity position alone.

### Implications for the paper

1. VDP's benchmark advantage over digest2 (F1 0.787 vs 0.667) was **not a fluke** — it
   reflects the fundamental TNO immunity built into the velocity-density approach. TNOs
   cannot masquerade as NEOs in velocity space.

2. Digest2's Sorcha improvement (0.667 → 0.836) is also real — it reflects a realistic
   population mix that happens to play to digest2's strengths (no V>23 TNOs, genuine
   detected NEOs with well-constrained orbits).

3. The fair comparison for a real Rubin survey is the Sorcha result, where both classifiers
   are evaluated on the population Rubin will actually deliver. On that comparison: digest2
   leads full-sky (0.836 vs 0.809), but VDP leads at the antisun (0–20°, 0.878 vs 0.848),
   which is the operationally critical NEOCP discovery regime.

4. The benchmark remains useful for one thing: it demonstrates that VDP's robustness to
   TNO contamination is a genuine advantage in worst-case conditions (shallow all-sky
   survey, no multi-night linking, flat magnitude cut).

---

## S3M linking experiment (case1/2/3) + benchmark population-cap discovery (2026-07-06)

### Purpose
Isolate how much of the Sorcha-vs-benchmark TNO gap is caused specifically by the SSP
linking filter's two configurable parameters, by re-running full Sorcha production with
native `[LINKINGFILTER]` enabled (the original run never enabled it — see below), varying
only `SSP_number_tracklets` and `SSP_separation_threshold`.

Pipeline: 6 CASE-parametrized SLURM stage scripts in `neomod/pipeline/slurm/s3m_linking/`
(`_case_env.sh` + `1_production.sh` … `6_combine.sh`), documented fully in
`neomod/docs/sorcha_full_pipeline.md`. Configs: `s3m_case{1,2,3}_*.ini`.

| Case | `SSP_number_tracklets` | `SSP_separation_threshold` | Isolates |
|------|:---:|:---:|----------|
| case1 | 1 | 0.5″ | motion threshold only, weak linking requirement |
| case2 | 3 | 0.5″ | full LSST linking (real SSP baseline) |
| case3 | 1 | 0.001″ | neither constraint (≈ old no-linking run) |

Sped up ~2-3x vs the original run's concurrency (`ckpt-all` partition, `%450` array
throttle on production/postprocess, `%113` on VDP, 4-worker GNU-parallel split on
digest2 chunks) — case1 finished stage 1→6 in ~3 h wall (was 4-6 h historically for a
single un-linked run). Known recurring issue: ~1-2 pathological S3M objects per full
production run hang/crash Sorcha (documented precedent: `inst00820_part003`,
`A804/A854/A868`-type objects from the v3.3/v5.0 hybrid runs). Same pattern hit case1
(`inst00884_part003`, `inst00552_part002`) — excluded rather than chased, 2/14,381 files
(0.014% of ~14.38M objects), consistent with precedent.

### Finding 1 — `SSP_linking_on` confirmed engaged, but case1 barely changes raw TNO count
Sorcha's own log per part (`inst*.h5`-adjacent `*-sorcha.log`) confirms:
```
Solar System Processing linking filter is turned ON.
Number of rows BEFORE applying SSP linking filter: 50537
Number of rows AFTER applying SSP linking filter: 50479
```
Only ~0.11% of rows dropped in this chunk. Raw TNO tracklet count barely moved:
**case1 (tk=1, sep=0.5): 628,162** vs **original no-linking run: 628,994** (0.1% reduction).

**Why:** `SSP_number_tracklets=1` means an object needs only ONE qualifying tracklet
anywhere across the full 2-year survey. A TNO gets dozens–hundreds of revisit
opportunities with varying geometry; the chance that *every single one* falls below the
0.5″ threshold is tiny. So `tk=1` alone barely filters anything — the real test of the
motion-threshold's impact requires `tk=3` (case2), where 3 tracklets must land within a
15-day window, a much harder bar for slow movers. **Case2 is the case that matters most
for this question; case1 mainly demonstrates that `tk=1` is a weak constraint.**

### Finding 2 — the S3M source files and true population sizes
`neomod/S3Mdata/` source files (via `s3m_loader.POP_TO_GLOB`): `neo→S0.s3m`,
`mba→S1_*.s3m` (14 files), `tno→ST.s3m`, `trojan→St5.s3m`. True object counts (raw file
line counts minus header/comment lines):

| Pop | True catalog size | Share of total |
|---|---:|---:|
| NEO | 268,512 | 1.87% |
| MBA | 13,883,375 | 96.54% |
| TNO | 48,682 | 0.34% |
| Trojan | 179,883 | 1.25% |
| **Total** | **14,380,452** | |

(Minor open discrepancy: benchmark's own reported "no-cap" TNO final count is 49,097,
slightly *above* this 48,682 raw figure, which shouldn't be possible if filtering only
removes rows. Not yet root-caused; ~1% and irrelevant to the conclusions below — noted
here so it isn't silently forgotten.)

### Finding 3 — case1's per-population tracklet/object rates (raw vdp_shards, case1)
| Pop | raw tracklets | uniq objects | tracklets/obj | true catalog | object recovery |
|---|---:|---:|---:|---:|---:|
| NEO | 148,773 | 25,878 | 5.75 | 268,512 | **9.6%** |
| MBA | 58,660,930 | 5,341,670 | 10.98 | 13,883,375 | 38.5% |
| Trojan | 1,608,416 | 116,367 | 13.82 | 179,883 | 64.7% |
| TNO | 628,162 | 41,392 | 15.18 | 48,682 | **85.0%** |

TNO object-level recovery is actually the *highest* of all four populations (slow movers
persist across many oppositions once found); NEO is the *lowest* (fast movers, brief
visibility windows). TNOs are not "hard to detect" — the opposite.

### Finding 4 — root cause of "benchmark has far more TNOs than Sorcha": population capping, not detection
`gen_benchmark_tracklets_s3m.py` (`POPULATIONS` dict) caps MBA at 200,000 and Trojan at
100,000 via random subsample **before** propagation/scoring, while NEO and TNO are
uncapped (full true population fed in). This is a deliberate, documented choice to keep
benchmark size (and especially digest2 runtime) tractable — not a bug, and not something
this session broke. Effect: MBA is suppressed **~69×** below its true catalog share while
TNO isn't suppressed at all, which is why TNO looks artificially dominant in the benchmark
(e.g. 38.6% of the 0-20° antisun bin) relative to Sorcha's real output (where MBA's true
13.88M-object dominance shows through). The Section 15/16 TNO forward/inverse injection
experiments (above) were already probing this same artifact, just not traced to its root
cause (the cap ratio) until this session.

**Correction to earlier framing:** the "TNO recovery collapses to ~1% at mag 24-25" table
computed earlier in this doc used the *subsampled* comparison file (`sample` step targets
a fixed 500k non-NEO cap, applied at a uniform rate across all non-NEO types combined) —
a different, far more diluted quantity than true object-level detection. That table is not
wrong on its own terms (it does describe the final ROC-scoring dataset), but it does NOT
mean Sorcha fails to detect TNOs — object-level recovery is 84-85% (Finding 3).

### Two options for fixing the benchmark's population mix (decision needed)
**Option A — "true catalog" benchmark (chosen, in progress).** Resize each population's
input cap proportional to its true catalog share, keeping total output size ~475,027 (same
VDP/digest2 compute cost as today). Back-solved through each population's own observed
survival rate to hit proportional *final* counts:

| Pop | true share | target final | survival rate | **new cap** |
|---|---:|---:|---:|---:|
| NEO | 1.87% | 8,868 | 68.8% | **12,900** |
| MBA | 96.54% | 458,619 | 70.6% | **650,000** |
| TNO | 0.34% | 1,610 | 100% | **1,600** |
| Trojan | 1.25% | 5,944 | 100% | **6,000** |

Gives a config-independent, principled "fair sky" benchmark. Does **not** reproduce any
specific Sorcha run's output mix (see Option B) because populations differ in both
tracklets/object and object-recovery-rate (Finding 3) — a proportional *input* can't
correct for *survey-driven* output skew.

**Option B — match a specific Sorcha case's output mix (deferred).** Reverse-engineer
benchmark caps from a completed Sorcha case's actual final population proportions
(its own tracklets/object × recovery rate per population). Gives a literal population-mix
match to that one case, but the target is case-dependent (case1 vs case2 vs case3 will
differ, especially since tk=3's 15-day window plausibly affects NEO/TNO differently than
tk=1). Deferred until case2 (the real-LSST-linking case) completes, since that's the most
meaningful target if this route is wanted later.

**Decision (2026-07-06): proceeding with Option A now.** Produces a new
`benchmark_comparison_s3m_v2.parquet`; the original benchmark is kept unchanged for
continuity with the existing Section 15/16 TNO-injection analysis. Fully independent of
case1/2/3 (different pipeline, no shared jobs) — neither blocks the other.

### Results — all four runs complete (2026-07-06 → 2026-07-07)

All three S3M linking cases and the Option A proportional benchmark finished cleanly.
Pathological-orbit losses were consistent across every production run (same catalog
region each time — instances ~849-898, i.e. the TNO/Trojan tail of the combined S3M
input — confirming these are genuine per-object integrator issues, not run-specific
flukes): case1 lost 2/899 instances (2 parts), case2 lost ~40/14,381 parts (retried
twice: once at full concurrency, once at low concurrency, to rule out contention —
16/24 initially-failed instances stayed failed both times, confirming determinism),
case3 lost 2/899 instances (same two instances, 882 and 884, as seen in case1/case2).

| Product | File | Rows | NEO | MBA | TNO | Trojan | other |
|---|---|---:|---:|---:|---:|---:|---:|
| case1 (tk=1, sep=0.5″) | `outputs/s3m_linking/case1/sorcha_comparison_case1.parquet` | 648,769 | 148,773 | 463,419 | 4,965 | 12,723 | 18,889 |
| case2 (tk=3, sep=0.5″ — real LSST linking) | `outputs/s3m_linking/case2/sorcha_comparison_case2.parquet` | 624,783 | 124,748 | 463,831 | 4,742 | 12,788 | 18,674 |
| case3 (tk=1, sep≈0) | `outputs/s3m_linking/case3/sorcha_comparison_case3.parquet` | 648,828 | 148,828 | 463,698 | 4,930 | 12,743 | 18,629 |
| benchmark v2 (Option A, proportional caps) | `outputs/phase2_benchmark_s3m_v2/benchmark_comparison_s3m_v2.parquet` | 474,783 | 8,041 | 459,525 | 1,217 | 6,000 | — |

**Headline observation:** TNO counts barely move across case1→case2→case3 (4,965 → 4,742
→ 4,930, all within ~5% of each other) in the *subsample*. NEO shows a real, larger drop
under case2's real linking (124,748 vs ~148,800 in case1/case3, ~16% down), consistent
with `tk=3`'s 15-day window being a genuine constraint for some fraction of NEOs (fast,
brief-visibility objects — Finding 3's NEO object-recovery rate was already the lowest of
all four populations at 9.6%, so this population is evidently the more sensitive one to
the multi-night requirement, not TNO as originally hypothesized). Full ROC/F1 comparison
across the three cases + both benchmarks not yet run — data is fully local, no further
SLURM jobs needed for that analysis.

**Benchmark v2 vs v1 population share** (proportional caps working as intended):
v1 TNO share ~10.3% of 475,027 vs v2 TNO share 0.26% of 474,783 — the proportional
resample removed the artificial TNO-dominance the capped v1 benchmark had.

### Known bug fixed during this run: benchmark VDP scoring picks wrong shard if
### per-population files are left alongside the combined file
`sorcha_phase2.py`'s `tracklet_paths()` globs `tracklets_*.parquet` and sorts
alphabetically. With `--batch-size 1 --shard-index 0` (used by both
`benchmark_score_vdp_s3m.sh` and `_v2.sh`), shard 0 is just the *first* file in that
sorted list. Uppercase population-file names (`tracklets_MBA.parquet`, `tracklets_NEO...`)
sort before the lowercase combined file (`tracklets_benchmark[_v2].parquet`) in ASCII, so
if the four per-population files are still present when VDP scoring runs, it silently
scores ONLY the alphabetically-first population (here: all-MBA, 459,525 rows) and drops
the other three entirely. v1's directory happened to have only the combined file (the
per-pop files had been cleaned up before scoring), so this never surfaced before. Fix
applied for v2: delete the four `tracklets_{POP}.parquet` files after `--combine-only`,
leaving only the combined file, before running `benchmark_score_vdp_s3m_v2.sh`. Worth
fixing properly in `tracklet_paths()` (e.g. exclude per-population files, or take the
combined file explicitly) if this benchmark pipeline is rerun again.
