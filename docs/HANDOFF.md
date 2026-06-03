# Sorcha / VDP / digest2 Pipeline — Session Handoff
Generated: 2026-06-01 (v7 — super-detailed for post-compact resume)

Read this document top to bottom before touching any code. Then read
`neomod/WAGG_SORCHA_HYAK_CONTEXT.md` for the full history.

---

## 1. Current state (what we have, what we want)

**Goal:** VDP F1 > digest2 F1 on the Sorcha/Wagg Rubin simulation.

| Classifier | F1 | Completeness | Contamination |
|---|---|---|---|
| VDP (S3M kNN, original) | 0.740 | 67.5% | 18.0% |
| VDP (hybrid kNN) | 0.787 | 71.0% | 11.8% |
| VDP (GMM NEO, mask ON) | 0.802 | 72.9% | 10.9% |
| **VDP (GMM NEO, mask OFF — current best)** | **0.809** | **74.5%** | **11.6%** |
| digest2 | 0.836 | 77.1% | 8.7% |

**Gap: 0.027 F1.** 113,812 NEO tracklets in subsample. Current breakdown:
- P > 0 (scored): 102,367 (89.9%)
- P = 0 (GMM density genuinely zero): 7,222 (6.3%) — unusual NEO orbits
- NaN (unscored, outside footprint): 4,223 (3.7%) — 30° patch edge

---

## 2. The NEOMOD3 augmentation plan (NEXT THING TO BUILD)

### Background — why this exists

During the original Mario-era work (`neomod/old_ones/`), the approach was to use
NEOMOD3 directly as a VDP scoring function:
- Given a tracklet's (ra, dec, dra, ddec), try many distances d and radial velocities ddot
- For each (d, ddot) pair, invert to get (a, e, i) and look up the NEOMOD3 4D table
- Marginalize over d and ddot to get P(NEO)

This FAILED for two reasons:
1. Marginalization over d and ddot was ~9 minutes per grid cell — completely impractical
   (see comment in `old_ones/antisunneo.ipynb`: "DO NOT RUN AGAIN OR ITS AN UNSKIPPABLE CUTSCENE")
2. The discrete bins (coarse 52×42×25×22 grid) produced blocky velocity maps that
   Zeljko rejected — he wanted smooth (vlam, vbeta) space maps, which led to switching
   to S3M + kNN for smooth continuous density estimation.

### Why it works now (completely different use)

We are NOT using NEOMOD3 as a scoring function. We use it to **sample training orbits
for the GMM**. The GMM smooths over all bin-boundary artifacts. The 9-minute
marginalization problem is completely irrelevant because we never look up a tracklet
in the table — we only sample from the table.

### What NEOMOD3 is

NEOMOD3 is Nesvorny et al.'s debiased NEO orbit population model (Boulder SWRI).
It is a 4D probability table: `array4D[iH, ia, ie, ii]` with bins:

```
n_H=52, H_min=15.0, H_max=28.0  →  dH = 0.25 mag
n_a=42, a_min=0.0,  a_max=4.2   →  da = 0.10 AU
n_e=25, e_min=0.0,  e_max=1.0   →  de = 0.04
n_i=22, i_min=0.0,  i_max=88.0  →  di = 4.0°
```

Bin centers:
```python
H_center = 15.0 + (np.arange(52) + 0.5) * 0.25
a_center = 0.0  + (np.arange(42) + 0.5) * 0.1
e_center = 0.0  + (np.arange(25) + 0.5) * 0.04
i_center = 0.0  + (np.arange(22) + 0.5) * 4.0
```

The values in the table are relative probabilities (debiased model weights).
Loading code (from old_ones/neomod-scratch.ipynb):
```python
import re, pandas as pd, numpy as np
from urllib.request import urlopen

url = "https://www2.boulder.swri.edu/~davidn/NEOMD_Simulator/input_neomod3.dat"
text = urlopen(url).read().decode("utf-8", errors="ignore")
lines = text.splitlines()
header_idx = next(i for i, ln in enumerate(lines)
                  if re.match(r"^\s*iH\s+ia\s+ie\s+ii\s+binned\s+mod\.\s+nu6\b", ln))
# parse header — "binned mod." is two tokens that become one column
# then parse data rows into df with columns: iH, ia, ie, ii, binned_mod, nu6, ...
# assemble array4D:
array4D = np.zeros((52, 42, 25, 22))
idx = df[['iH','ia','ie','ii']].astype(int).values - 1  # 0-based
array4D[idx[:,0], idx[:,1], idx[:,2], idx[:,3]] = df['binned_mod'].values
```

**Why this fills the 7,222 P=0 gap:** NEOMD3 is DEBIASED — it represents the TRUE
underlying NEO distribution, including:
- High-inclination NEOs (i > 40°) underrepresented in MPCORB (biased toward low-i)
- Aten-class objects (a < 1 AU) with unusual fast apparent motion
- Very high eccentricity objects (e > 0.8) rarely caught at opposition
- Objects currently far from Earth that never appear visible at training epoch

These are exactly the orbital configurations the current GMM misses because it was
trained on VISIBLE source NEOs per mag bin — a biased sample.

### Exact implementation plan

**New file:** `neomod/src/neomod3_sampler.py`

This module does one thing: given the NEOMD3 4D array, sample N orbital elements.

```python
import numpy as np
import pandas as pd
from pathlib import Path

NEOMOD3_URL = "https://www2.boulder.swri.edu/~davidn/NEOMD_Simulator/input_neomod3.dat"
NEOMOD3_LOCAL = Path("/mmfs1/gscratch/astro/ds2004/sorcha/neomod/NEOMOD3/input_neomod3.dat")

# Bin definitions (must match old_ones notebooks exactly)
N_H, H_MIN, H_MAX = 52, 15.0, 28.0
N_A, A_MIN, A_MAX = 42, 0.0, 4.2
N_E, E_MIN, E_MAX = 25, 0.0, 1.0
N_I, I_MIN, I_MAX = 22, 0.0, 88.0

def load_neomod3_array():
    """Load NEOMD3 4D probability table. Downloads from URL if local copy absent."""
    import re
    if NEOMOD3_LOCAL.exists():
        text = NEOMOD3_LOCAL.read_text(errors="ignore")
    else:
        from urllib.request import urlopen
        text = urlopen(NEOMOD3_URL).read().decode("utf-8", errors="ignore")
        NEOMOD3_LOCAL.parent.mkdir(parents=True, exist_ok=True)
        NEOMOD3_LOCAL.write_text(text)
    lines = text.splitlines()
    header_idx = next(i for i, ln in enumerate(lines)
                      if re.match(r"^\s*iH\s+ia\s+ie\s+ii\s+binned\s+mod\.", ln))
    rows = []
    for ln in lines[header_idx+1:]:
        parts = ln.split()
        if len(parts) >= 5:
            try:
                rows.append([int(parts[0]), int(parts[1]),
                              int(parts[2]), int(parts[3]), float(parts[4])])
            except ValueError:
                continue
    df = pd.DataFrame(rows, columns=['iH','ia','ie','ii','weight'])
    da = (A_MAX - A_MIN) / N_A
    de = (E_MAX - E_MIN) / N_E
    di = (I_MAX - I_MIN) / N_I
    dH = (H_MAX - H_MIN) / N_H
    array4D = np.zeros((N_H, N_A, N_E, N_I))
    idx = df[['iH','ia','ie','ii']].astype(int).values - 1
    array4D[idx[:,0], idx[:,1], idx[:,2], idx[:,3]] = df['weight'].values
    centers = {
        'H': H_MIN + (np.arange(N_H) + 0.5) * dH,
        'a': A_MIN + (np.arange(N_A) + 0.5) * da,
        'e': E_MIN + (np.arange(N_E) + 0.5) * de,
        'i': I_MIN + (np.arange(N_I) + 0.5) * di,
    }
    widths = {'H': dH, 'a': da, 'e': de, 'i': di}
    return array4D, centers, widths


def sample_neomod3_orbits(n_samples, obstime_str, rng=None,
                           array4D=None, centers=None, widths=None):
    """Sample n_samples NEO orbits from the NEOMD3 debiased distribution.

    Returns a DataFrame with columns compatible with the VDP pipeline:
        a, e, i, node, argperi, t_p, H, q
    where node, argperi, M are drawn uniformly (phase angles have no preferred
    value in the debiased population model).
    t_p is reconstructed from M at obstime_str.
    """
    from astropy.time import Time
    from neomod.src.velocity_density_pipeline_gmm import MU_SUN, AU_KM
    # allow lazy loading
    if array4D is None:
        array4D, centers, widths = load_neomod3_array()
    if rng is None:
        rng = np.random.default_rng(42)

    # Flatten 4D → probability vector, sample bin indices
    flat = array4D.flatten()
    flat = np.maximum(flat, 0.0)
    flat /= flat.sum()
    flat_idx = rng.choice(len(flat), size=n_samples, p=flat)
    iH, ia, ie, ii = np.unravel_index(flat_idx, array4D.shape)

    # Perturb within bin (uniform within bin width)
    H = centers['H'][iH] + rng.uniform(-widths['H']/2, widths['H']/2, n_samples)
    a = centers['a'][ia] + rng.uniform(-widths['a']/2, widths['a']/2, n_samples)
    e = centers['e'][ie] + rng.uniform(-widths['e']/2, widths['e']/2, n_samples)
    i = centers['i'][ii] + rng.uniform(-widths['i']/2, widths['i']/2, n_samples)

    # Phase angles: uniform (debiased model has no angular preference)
    node   = rng.uniform(0.0, 360.0, n_samples)
    argperi = rng.uniform(0.0, 360.0, n_samples)
    M_deg  = rng.uniform(0.0, 360.0, n_samples)

    # Reconstruct t_p from M at obstime
    t_obs = Time(obstime_str, scale='tdb').mjd
    a_km = a * AU_KM
    n_rad_day = np.sqrt(MU_SUN / a_km**3) * 86400.0
    M_rad = np.deg2rad(M_deg)
    t_p = t_obs - M_rad / n_rad_day

    # Physical validity filter
    q = a * (1.0 - e)
    valid = (
        (a > 0.05) & (a < 4.2) &
        (e >= 0.0) & (e < 1.0) &
        (q >= 0.0) & (q < 1.3) &
        (i >= 0.0) & (i <= 90.0) &
        np.isfinite(H) & np.isfinite(a) & np.isfinite(e) & np.isfinite(i)
    )

    return pd.DataFrame({
        'a': a[valid], 'e': e[valid], 'i': i[valid],
        'node': node[valid], 'argperi': argperi[valid],
        't_p': t_p[valid], 'H': H[valid], 'q': q[valid],
    })
```

**Integration point in `velocity_density_pipeline_gmm.py`:**

In `_clone_neo_gmm`, after building the training DataFrame from the visible source
NEOs (`df` = df_cloner_input), augment it with NEOMD3 samples BEFORE fitting the GMM:

```python
# In _clone_neo_gmm, after train = df.loc[valid].copy():

from neomod3_sampler import sample_neomod3_orbits, load_neomod3_array

# Lazy-load NEOMD3 once per map generation (cache via module-level variable)
_NEOMOD3_CACHE = None
def _get_neomod3():
    global _NEOMOD3_CACHE
    if _NEOMD3_CACHE is None:
        _NEOMD3_CACHE = load_neomod3_array()
    return _NEOMD3_CACHE

# Sample NEOMD3 orbits — oversample relative to visible training size
n_neomod3 = max(500, len(train) * 10)   # 10x augmentation
array4D, centers, widths = _get_neomod3()
neomod3_df = sample_neomod3_orbits(n_neomod3, obstime_str, rng=rng,
                                    array4D=array4D, centers=centers,
                                    widths=widths)

# Augment training set
train_augmented = pd.concat([train, neomod3_df], ignore_index=True)
print(f"  NEOMD3 augmentation: {len(train)} visible → {len(train_augmented)} total training orbits")

# Then use train_augmented instead of train for GMM feature matrix:
X, feature_names = _gmm_feature_matrix(train_augmented, obstime_str)
# ... rest of _clone_neo_gmm unchanged
```

**Important:** The NEOMD3 orbits don't have `M_obs_deg` or `q` as existing columns
the way `df_cloner_input` does. The `sample_neomod3_orbits` function ADDS q and
reconstructs t_p from M, so the returned DataFrame is compatible with
`_gmm_feature_matrix` which calls `_gmm_as_q(df)` and `_gmm_mean_anomaly_rad(df, obstime_str)`.

**Tunable parameter:** `n_neomod3 = max(500, len(train) * 10)` — controls how much
NEOMD3 data augments each mag bin. Higher values = better coverage of unusual orbits
but slower GMM fitting. Start with 10× and increase if P=0 count doesn't drop enough.

**Expected gain:** The 7,222 P=0 NEOs that are in unusual orbital configurations
should now have training orbits nearby in the GMM. Expected F1 +0.010–0.020.

---

## 3. Other planned improvements

### Widen antisun footprint to 45°

**Current state:** `max_sep_deg=30` in `sorcha_postprocess.py::assign_probability_maps`
and in `build_visible_subset_dataframe` calls.

**Change needed:**
1. In `sorcha_postprocess.py`, change `MAX_SEP_DEG = 30` (or wherever the constant is)
   to `45`. This affects Phase 1 map assignment.
2. Regenerate all 24 GMM maps (with `max_sep_deg=45` in the VDP pipeline).
3. Re-run Phase 1 with `--overwrite` (already in sorcha_postprocess.sh).
4. Re-run Phase 2 GMM.

**Risk:** At 40–45° from antisun, MBA velocities shift (less retrograde, more variable),
potentially increasing contamination. Test with a single map first.

**Expected gain:** ~3,000–4,000 of the 4,223 NaN NEOs recovered, F1 +0.010–0.015.

### Increase GMM components (80 → 200)

In `velocity_density_pipeline_gmm.py`, in the `_clone_neo_gmm` call at line ~1370:
```python
# Change:
gmm_clone_df, _, _, gmm_diag = _clone_neo_gmm(
    df=df_cloner_input, n_clones=n_clones_target, obstime_str=obstime_str,
    n_components=80,   # ← change to 200
    random_state=42,
)
```

The safeguard inside `_clone_neo_gmm`:
```python
n_comp = min(int(n_components), max(1, len(X_scaled) // 5), len(X_scaled))
```
automatically reduces components for sparse bins. Safe to increase to 200.

---

## 4. Active pipeline — exact commands

### Generate GMM maps
```bash
cd /mmfs1/gscratch/astro/ds2004/sorcha
sbatch --array=0-23%24 sorcha_gen_maps_gmm_slurm.sh
# Writes: prob_maps_gmm/prob_maps_YYYY-MM-DD_antisun.npz (24 files)
# Time: ~20–40 min all concurrent
```

### Run Phase 2 GMM scoring
```bash
rm -rf outputs/phase2_gmm/vdp_shards
sbatch sorcha_phase2_vdp_gmm.sh
# Has --no-nearest-dist-mask (mask disabled for GMM)
# Writes: outputs/phase2_gmm/vdp_shards/vdp_00000.parquet … vdp_00112.parquet
# Time: ~15 min
```

### Combine and compute F1 (memory-safe)
```python
# Run on login node — load only needed columns to avoid OOM
import pandas as pd, numpy as np
from glob import glob
from sklearn.metrics import precision_recall_curve

vdp = pd.concat([pd.read_parquet(f, columns=['tracklet_id','P_NEO_vdp'])
                 for f in sorted(glob('outputs/phase2_gmm/vdp_shards/vdp_*.parquet'))])
base = pd.read_parquet('outputs/phase2_hybrid/wagg_sorcha_comparison_hybrid.parquet',
                       columns=['tracklet_id','population','P_NEO_d2'])
out = base.merge(vdp, on='tracklet_id', how='left')

neo = (out['population'] == 'NEO')
y_true = neo.astype(int).values

def best_f1(scores, y_true):
    scores = np.nan_to_num(scores, nan=0.0)
    prec, rec, _ = precision_recall_curve(y_true, scores)
    f1s = 2*prec*rec/(prec+rec+1e-9)
    idx = np.argmax(f1s)
    return f1s[idx], rec[idx], 1-prec[idx]

f1_vdp, comp_vdp, cont_vdp = best_f1(out['P_NEO_vdp'].values, y_true)
f1_d2, comp_d2, cont_d2   = best_f1(out['P_NEO_d2'].values, y_true)
print(f"VDP GMM  F1={f1_vdp:.3f}  comp={comp_vdp*100:.1f}%  cont={cont_vdp*100:.1f}%")
print(f"digest2  F1={f1_d2:.3f}   comp={comp_d2*100:.1f}%  cont={cont_d2*100:.1f}%")
print(f"NEO P>0: {(neo & (out['P_NEO_vdp']>0)).sum():,}")
print(f"NEO P=0: {(neo & (out['P_NEO_vdp']==0)).sum():,}")
print(f"NEO NaN: {(neo & out['P_NEO_vdp'].isna()).sum():,}")
```

---

## 5. Switching between pipelines

| Pipeline | Maps dir | Phase 2 script | Mask | Notes |
|---|---|---|---|---|
| S3M kNN | `prob_maps/` | `sorcha_phase2_vdp.sh` | ON | original |
| Hybrid kNN | `prob_maps_hybrid/` | `sorcha_phase2_vdp_hybrid.sh` | ON | baseline |
| GMM (current) | `prob_maps_gmm/` | `sorcha_phase2_vdp_gmm.sh` | OFF | active |

To re-enable nearest_dist mask for GMM: remove `--no-nearest-dist-mask` from
`sorcha_phase2_vdp_gmm.sh`.

---

## 6. Key file locations

| File | Role | Status |
|---|---|---|
| `hybrid.h5` | 14.4M Cartesian state vectors | READ ONLY |
| `hybrid_elements.parquet` | 14.4M Keplerian + population | READ ONLY |
| `neomod/src/hybrid_loader.py` | Loads hybrid_elements.parquet | READ ONLY |
| `neomod/src/velocity_density_pipeline_hybrid.py` | Hybrid kNN pipeline | READ ONLY |
| `neomod/src/velocity_density_pipeline_gmm.py` | GMM pipeline (active) | MODIFY for NEOMD3 |
| `neomod/src/neomod3_sampler.py` | TO BE CREATED — NEOMD3 sampler | BUILD NEXT |
| `neomod/NEOMOD3/input_neomod3.dat` | TO BE DOWNLOADED | Download on first use |
| `neomod/old_ones/neomod-scratch.ipynb` | NEOMD3 loading reference code | READ ONLY |
| `neomod/old_ones/antisunneo.ipynb` | NEOMD3 scoring prototype (history) | READ ONLY |
| `neomod/cloning_test_ZI_v2.ipynb` | GMM cloner reference notebook | READ ONLY |
| `sorcha_gen_maps_gmm_slurm.sh` | Slurm: generate 24 GMM maps | MODIFY if pipeline changes |
| `sorcha_phase2_vdp_gmm.sh` | Slurm: Phase 2 GMM scoring | Has --no-nearest-dist-mask |
| `sorcha_phase2.py` | Shared Phase 2 pipeline | Has --no-nearest-dist-mask flag |
| `sorcha_postprocess.sh` | Phase 1 Slurm (--overwrite) | Ready for footprint change |
| `prob_maps_gmm/` | 24 monthly GMM maps (active) | Regen after pipeline changes |
| `outputs/phase2_hybrid/wagg_sorcha_comparison_hybrid.parquet` | kNN baseline result | READ ONLY (used as d2 source) |
| `outputs/phase2_gmm/` | GMM scoring results | Active output directory |
| `production_run/` | Archived production scripts | READ ONLY |

---

## 7. Python environment

```bash
/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python
```

Always use this full path. Shell shows `(base)` — that is NOT this environment. Ignore it.

Slurm partition: `ckpt` (default). `astro` account also has `cpu-g2` (778 free CPUs,
can use with `--partition=cpu-g2 --account=astro` to bypass ckpt congestion).

QOS limit on ckpt: max 2 job arrays at once. Check with `hyakalloc`.

---

## 8. What NOT to modify

- `velocity_density_pipeline.py` — original, preserved
- `velocity_density_pipeline_fast.py` — S3M fast pipeline, preserved
- `velocity_density_pipeline_hybrid.py` — hybrid kNN, preserved as clean baseline
- `s3m_loader.py` — S3M loader, preserved
- `prob_maps/` — 24 S3M monthly maps
- `prob_maps_hybrid/` — 24 hybrid kNN maps
- `outputs/phase2/` — S3M scoring results
- `outputs/phase2_hybrid/wagg_sorcha_comparison_hybrid.parquet` — kNN result (used as source for d2 columns in combine)

---

## 9. F1 progression history

```
0.342  baseline (5 maps, wrong epochs)
0.448  24 monthly maps + epoch-aware assignment
0.502  + antisun footprint fix
0.648  + antisun-only filter in Phase 2
0.740  + wide grid (-1.5 to +1.5)
0.787  + hybrid catalog training + Phase 1 re-run
0.802  + GMM NEO cloner (mask ON)
0.809  + nearest-dist mask OFF   ← current best
------
0.836  digest2 (target to beat)
```

---

## 10. Paper argument (for context)

The S3M-only comparison showed VDP F1=0.856 vs digest2 F1=0.665. That result stands.

The Sorcha comparison is harder because ~10% of objects are real MPCORB objects with
unusual orbits not in S3M. The NEOMD3 augmentation closes this gap by training on the
TRUE debiased NEO distribution rather than the biased observed sample.

The narrative: VDP was broken by epoch mismatch (F1=0.342) and we've brought it to
0.809 through calibration. NEOMD3 augmentation should push it past digest2.

---

## 11. Step-by-step: building NEOMD3 augmentation

1. **Download NEOMD3 data** (test on login node):
   ```bash
   mkdir -p /mmfs1/gscratch/astro/ds2004/sorcha/neomod/NEOMOD3
   wget "https://www2.boulder.swri.edu/~davidn/NEOMD_Simulator/input_neomod3.dat" \
        -O /mmfs1/gscratch/astro/ds2004/sorcha/neomod/NEOMOD3/input_neomod3.dat
   ```

2. **Verify loading** using reference code from `old_ones/neomod-scratch.ipynb` cells 1-2.
   The file has a header line matching `iH  ia  ie  ii  binned mod.  nu6`.
   After loading: `array4D.sum()` should be nonzero, `array4D.shape = (52, 42, 25, 22)`.

3. **Create `neomod/src/neomod3_sampler.py`** with the code in Section 2 above.

4. **Modify `velocity_density_pipeline_gmm.py`**:
   - Add `import neomod3_sampler` at top
   - In `_clone_neo_gmm`, augment `train` with NEOMD3 samples before GMM fit
   - See exact code in Section 2 above

5. **Test with a single map** (one month, watch the logs):
   ```bash
   # Test on login node first (will be slow, just to check it runs):
   /mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python sorcha_gen_map_gmm.py \
     --obstime 2026-01-01T00:00:00 \
     --center-lon <antisun_lon> --center-lat 0.0 \
     --center-label antisun \
     --output /tmp/test_neomod3_map.npz --n-jobs 1 --overwrite
   ```
   Check log output for "NEOMD3 augmentation: X visible → Y total training orbits"

6. **Check velocity coverage** of the test map:
   ```python
   import numpy as np
   m = np.load('/tmp/test_neomod3_map.npz', allow_pickle=True)
   nd = m['nearest_dist__NEO__18_20']
   xg = m['x_grid']
   mid = len(xg)//2
   for vlam in [-1.5, -1.3, -1.0, -0.5, 0.0, 0.5, 1.0, 1.3, 1.5]:
       idx = min(np.searchsorted(xg, vlam), len(xg)-1)
       print(f"vlam={vlam:+.1f}  nd={nd[idx,mid]:.4f}")
   # With NEOMD3 aug, expect coverage to extend further (lower nd at extremes)
   ```

7. **Submit all 24 maps** once test looks good:
   ```bash
   rm prob_maps_gmm/*.npz
   sbatch --array=0-23%24 sorcha_gen_maps_gmm_slurm.sh
   ```

8. **Run Phase 2 and combine** (see Section 4 above).

9. **Evaluate** — compare P=0 count before/after NEOMD3 augmentation:
   - Before: 7,222 P=0 NEOs
   - After: should be significantly lower if NEOMD3 fills orbital gaps
