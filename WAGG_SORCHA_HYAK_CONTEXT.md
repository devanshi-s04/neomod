# Wagg / Sorcha / Hyak Context
Generated: 2026-05-27, Updated: 2026-05-27
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

## Post-Processing Pipeline (NEXT STEPS)

### Key Insight: No Cross-File Merging Needed

All detections for a given ObjID are contained within a single h5 file (same Sorcha instance/part). Nightly tracklets are therefore self-contained per file. Process each h5 file independently in a Slurm array.

### Full Pipeline

```
Phase 1 — Slurm array job over all 14,445 h5 files
  wagg_postprocess.py --infile inst*.h5 --outfile tracklets_inst*.parquet
  For each h5 file:
    1. Read sorcha_results
    2. Group by (ObjID, night = floor(fieldMJD_TAI))
    3. Keep nights with ≥2 detections, time span ≤90 min
    4. Per tracklet: mean RA/Dec/rates/mag, first+last detection for MPC obs
    5. Assign to nearest prob_map within 30°
  Output: small tracklet parquet per file

Phase 2 — single combine job
  1. Concatenate all tracklet parquets → wagg_tracklets_all.parquet
  2. VDP score_observation on all tracklets (vectorized, fast)
  3. Build MPC 80-col observations → run digest2 in 5000-tracklet chunks
  4. Merge VDP + digest2 scores → wagg_sorcha_comparison.parquet

Phase 3 — local analysis (download ~few-GB parquet to Mac)
  ROC curves, F1 metrics, population breakdowns — same notebook as digest2_comparison.ipynb
```

### Open Question: Population Labels

`hybrid.h5` has no explicit population column (only id, x, y, z, vx, vy, vz, t_0, H, g). Two options:

**Option A** — classify from state vectors in Sorcha output:
Sorcha preserves x, y, z, xdot, ydot, zdot per detection. Convert to orbital elements, then:
- q = a(1−e) < 1.3 AU → NEO
- 2.0 < a < 3.3, e < 0.3 → MBA
- a > 30 AU → TNO
- a ≈ 5.2 AU, librating → Trojan

**Option B** — check hybrid.h5 id column format (may encode population):

```bash
/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python - <<'PY'
import pandas as pd
df = pd.read_hdf("hybrid.h5", key="df")
print(df["id"].head(30).tolist())
PY
```

Run this first — if IDs have population prefixes, Option B is simpler.

### Slurm Script Pattern for Post-Processing

```bash
#!/bin/bash
#SBATCH --job-name=wagg_postprocess
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-14444%32
#SBATCH --output=logs/postprocess_%A_%a.out
#SBATCH --error=logs/postprocess_%A_%a.err

/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python wagg_postprocess.py \
  --file_index ${SLURM_ARRAY_TASK_ID} \
  --indir outputs/production_2yr \
  --outdir outputs/tracklets \
  --prob_maps_dir prob_maps
```

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

Our existing S3M comparison used 2 detections × 30 min separation. Whether to unify tracklet definition across S3M and Sorcha comparisons is an open advisor question.

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

## Current Status (2026-05-27)

Completed:
- 2yr Sorcha production run: 14,445 h5 files, 0 missing, ~232G
- neomod cloned to Hyak at `/mmfs1/gscratch/astro/ds2004/sorcha/neomod/`
- VDP import working with adam_core stub
- Probability maps uploaded to `prob_maps/` (5 maps covering 4 sky directions)
- Observatory code confirmed X05
- VSCode Remote SSH connected and sorcha folder open

Pending before writing `wagg_postprocess.py`:

1. **Population label check** — run on Hyak to see if ObjID encodes population:

```bash
/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python - <<'PY'
import pandas as pd
df = pd.read_hdf("hybrid.h5", key="df")
print(df.shape)
print(df.columns.tolist())
print(df["id"].head(30).tolist())
PY
```

If IDs have population prefixes → use directly. If not → classify from state vectors stored in Sorcha output (q = a(1−e) < 1.3 AU → NEO, 2.0 < a < 3.3 → MBA, a > 30 → TNO, a ≈ 5.2 → Trojan).

2. **Advisor decision** on whether to redo S3M comparison with 3-detection tracklets for consistency with Wagg, or keep 2-detection for S3M and 3-detection for Sorcha.

---

## Relationship to Existing VDP / digest2 Comparison

Existing S3M-based comparison (`neomod/digest2_comparison.ipynb`):
- Input: S3M synthetic orbital elements, propagated to 2025-03-21 opposition patch
- Tracklet: 2 detections × 30 min, synthetic
- Result: VDP F1=0.856, digest2 F1=0.665

Wagg-Sorcha comparison (this pipeline):
- Input: Rubin-cadence detections from 14.4M hybrid catalog over 2 years
- Tracklet: built from real simulated detection cadence (nightly grouping)
- Same VDP maps, same digest2 binary, same scoring logic
- This is the more realistic validation — actual Rubin observing patterns, real detection noise/limits, full population mix

The Wagg comparison is the upgrade: same classifier, harder test.
