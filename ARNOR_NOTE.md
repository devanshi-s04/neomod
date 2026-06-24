# Note for Arnor Claude — 2026-06-03

Hi! This is a handoff note from Hyak Claude. A lot happened on Hyak today. Here's what you need to know before touching anything on Arnor.

---

## What happened on Hyak today

### 1. Git repo fully reorganised (neomod on GitHub)

The neomod repo (`git@github.com:devanshi-s04/neomod.git`) was restructured in 4 commits today (`6d756c8` → `74762db` → `8e75465` → `34fdd50`). New layout:

```
neomod/
├── src/                 ← library code (5 new files added: fast, gmm, hybrid pipelines)
├── pipeline/            ← sorcha pipeline .py scripts (moved from sorcha/ root)
│   ├── slurm/           ← all Slurm .sh scripts
│   └── config/          ← Rubin .ini config files
├── notebooks/
│   ├── paper/           ← 6 paper/advisor notebooks (including sorcha_gmm_s3m_comparison.ipynb)
│   ├── dev/             ← 40 exploration notebooks
│   └── sorcha/          ← sorcha_roc_comparison.ipynb
├── docs/                ← all context/handoff .md files (including this file)
└── adam_core_stub/      ← minimal adam_core mock (was at sorcha/ root on Hyak)
```

### 2. All "wagg" references removed from file names and code

Every output file was renamed. The key ones you'll care about on Arnor:

| Old name | New name |
|----------|----------|
| `wagg_sorcha_comparison.parquet` | `sorcha_comparison.parquet` |
| `wagg_sorcha_comparison_gmm.parquet` | `sorcha_comparison_gmm.parquet` |
| `wagg_sorcha_comparison_hybrid.parquet` | `sorcha_comparison_hybrid.parquet` |
| `wagg_subsample.parquet` | `sorcha_subsample.parquet` |
| `wagg_sorcha_may2026_antisun_patch.parquet` | `sorcha_may2026_antisun_patch.parquet` |

**You need to rename these on Arnor too** (they only got renamed on Hyak):
```bash
cd /astro/users/ds2004/vdp/outputs/phase2
mv wagg_sorcha_comparison.parquet sorcha_comparison.parquet
mv wagg_sorcha_comparison_gmm.parquet sorcha_comparison_gmm.parquet
mv wagg_sorcha_comparison_hybrid.parquet sorcha_comparison_hybrid.parquet
mv wagg_subsample.parquet sorcha_subsample.parquet        # if it exists
mv wagg_sorcha_may2026_antisun_patch.parquet sorcha_may2026_antisun_patch.parquet
```

### 3. New parquet scped to Arnor today

`sorcha_may2026_antisun_patch.parquet` (1.2 MB) should have arrived at:
```
/astro/users/ds2004/vdp/outputs/phase2/sorcha_may2026_antisun_patch.parquet
```
This is the single-epoch May 2026 antisun patch (14,873 objects, 2,402 NEOs, 16.2% NEO rate).
Columns: `tracklet_id, population, P_NEO_s3m, P_NEO_d2, vlam, vbeta, mean_mag, ecl_lon, ecl_lat, mjd0_utc, P_NEO_gmm`

---

## What Arnor has that is NOT yet in git — action needed

### Priority 1: `sorcha_gmm_s3m_singleepoch_comparison.ipynb`

This notebook was **created entirely on Arnor** and does not exist in the git repo at all. It's the single-epoch same-objects comparison (re-scores S3M objects with GMM map using RegularGridInterpolator). It should go in `notebooks/paper/` in the repo.

**Action:** After `git pull`, copy it in and commit:
```bash
cd /astro/users/ds2004/vdp    # wherever it lives on Arnor
git pull origin main
cp sorcha_gmm_s3m_singleepoch_comparison.ipynb neomod/notebooks/paper/
cd neomod
git add notebooks/paper/sorcha_gmm_s3m_singleepoch_comparison.ipynb
git commit -m "Add single-epoch same-objects comparison notebook (created on Arnor)"
git push origin main
```

### Priority 2: `sorcha_gmm_s3m_comparison.ipynb` — your version vs git version

The version in git (`notebooks/paper/sorcha_gmm_s3m_comparison.ipynb`) is the Hyak original.
Your Arnor version has fixes applied by you (path corrections for Arnor layout, dpi 200→150, mag_min/max key fix). The Hyak version in git has the path references pointing at `../outputs/phase2/sorcha_comparison.parquet` etc. (already updated for the wagg rename).

**Action:** Check if your Arnor version has all the wagg renames applied. If yes, your version is probably better (has your fixes). You can overwrite the git version with yours after confirming paths match:
```bash
# On Arnor, check the parquet paths in the notebook
grep -o '"[^"]*\.parquet"' sorcha_gmm_s3m_comparison.ipynb | sort -u
# Should show sorcha_comparison*, not wagg_sorcha_comparison*
```
If paths still say "wagg", update them first, then commit your version over the git version.

---

## What is and isn't in git (summary)

### IN git (pull to get):
- All Python source files (`src/velocity_density_pipeline*.py`, loaders, etc.)
- All pipeline scripts (`pipeline/sorcha_phase2.py`, `sorcha_postprocess.py`, `sorcha_gen_map*.py`)
- All Slurm scripts (`pipeline/slurm/*.sh`)
- All notebooks in `notebooks/paper/`, `notebooks/dev/`, `notebooks/sorcha/`
- All docs in `docs/` (including `WAGG_SORCHA_HYAK_CONTEXT.md`)
- `adam_core_stub/`

### NOT in git (large data, Hyak-only):
- All `outputs/` parquets (`sorcha_comparison.parquet`, shards, etc.)
- All `prob_maps/` and `prob_maps_gmm/` .npz files
- `NEOMD3/input_neomd3.dat` (62 MB, gitignored)
- `S3Mdata/` (large S3M input files)

### On Arnor only (not in git, not on Hyak):
- `sorcha_gmm_s3m_singleepoch_comparison.ipynb` ← needs to go to git
- The output parquets that live in `outputs/phase2/` on Arnor

---

## Current science status

VDP F1 = **0.837**, digest2 F1 = **0.837** — tied. Goal is F1 > 0.837.

**Next priority (not done yet):** Increase MBA `clone_factor` from 1 → 5 in `velocity_density_pipeline_gmm.py` to reduce MBA false positive rate (1.5% → ~0.4%), targeting F1 ≈ 0.857. This needs to happen on Hyak (map regeneration via Slurm).

See `docs/WAGG_SORCHA_HYAK_CONTEXT.md` for full context.

---

Good luck! — Hyak Claude 🙂
