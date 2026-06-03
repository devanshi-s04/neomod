# Project Context: Solar System Population Velocity-Space Classifier

## What this project does

Builds ecliptic velocity-space probability maps — P(population | v_λ, v_β) — for four
S3M synthetic solar system populations (MBA, NEO, TNO, Trojan). These maps are used to
classify observed solar system objects and compute ROC curves for threshold selection.

The core deliverable is `velocity_density_pipeline.py`, a self-contained Python module
that handles everything from S3M loading through map generation, persistence, lookup,
and ROC-curve scoring.

---

## Populations

| Label    | S3M key    | Size (rows) | Notes                        |
|----------|------------|-------------|------------------------------|
| MBA      | `"mba"`    | ~13.9M      | Dominates runtime            |
| NEO      | `"neo"`    | ~269k       |                              |
| TNO      | `"tno"`    | ~49k        |                              |
| Trojans  | `"trojan"` | ~180k       |                              |

S3M DataFrames have 15 columns; the pipeline uses: `a, e, i, node, argperi, t_p, H`.

---

## Key files

```
velocity_density_pipeline.py   ← the entire pipeline (single file, ~1950 lines)
neoscore.py                    ← user's package: propagation, observer state, scoring
s3m_loader.py                  ← user's package: loads S3M population DataFrames
```

`velocity_density_pipeline.py` imports both packages at module level:
```python
import neoscore as nsc
import s3m_loader as load_s3m
```

---

## Public API surface

### Generation (slow, run once per obstime)
```python
import velocity_density_pipeline as vdp

vdp.generate_probability_maps(
    obstime_str="2025-03-21T00:00:00",   # ISO, TDB scale
    output_path="prob_maps_2025-03-21.npz",
    # optional overrides:
    center_lon_deg=180.0,   # ecliptic longitude of sky-patch center
    center_lat_deg=0.0,
    max_sep_deg=30.0,       # sky-patch radius
    chunk=500_000,          # Newton propagation chunk size
    seed=42,
)
```

Internally: loads all 4 populations via `load_all_populations()`, iterates over 8
magnitude bins, runs conditional K|M cloner + Bayesian kNN density estimator per
population per bin, saves to one `.npz` file.

### Loading (fast, use from any notebook)
```python
pms = vdp.ProbMapSet.from_npz("prob_maps_2025-03-21.npz")
```

### Scoring — Case 1: single observation (no S3M needed)
```python
probs = pms.score_observation(
    ra_deg=183.2, dec_deg=-1.5,
    dra_deg_day=-0.35, ddec_deg_day=0.04,
    mag_app=21.7,
)
print(probs["NEO"])   # array([0.243...])
```
Uses pure-math manual equatorial→ecliptic rotation. No scorer/propagation needed.
`apply_mask=True` by default (0.2 deg/day nearest-object mask applied — correct for
single-object scoring).

### Scoring — Case 2: bulk S3M for ROC analysis
```python
df, scorer = vdp.load_s3m_population("neo")
probs, visible_df = pms.score_orbital_df(
    df, scorer=scorer,
    chunk=500_000,
    apply_mask=False,    # ← MUST be False for ROC; True zeros out sparse regions
)
# visible_df columns: a, e, i, node, argperi, t_p, H,
#                     ra_deg, dec_deg, sky_sep_deg, vlam, vbeta,
#                     mag_app, mag_bin_label,
#                     P_NEO, P_MBA, P_TNO, P_Trojans
```

`apply_mask=False` is the default for `score_orbital_df`. It returns raw density
ratios (positive wherever kNN put any density). `apply_mask=True` uses the 0.2
deg/day mask (same as the plot), which zeroes out sparsely-sampled regions and makes
TPR start below 1.0 at X=0 — wrong for ROC.

### Plotting (probability maps with 0.2 deg/day mask)
```python
fig, axes = vdp.plot_probability_maps(pms, mag_bin_label="mag21", mask_radius=0.2)
```

---

## Magnitude bins

8 apparent-magnitude bins (not H bins):

| Label    | mag_min | mag_max |
|----------|---------|---------|
| 14_16    | 14.0    | 16.0    |
| 16_18    | 16.0    | 18.0    |
| 18_20    | 18.0    | 20.0    |
| mag20    | 20.0    | 21.0    |
| mag21    | 21.0    | 22.0    |
| mag22    | 22.0    | 23.0    |
| mag23    | 23.0    | 24.0    |
| mag24+   | 24.0    | 26.0    |

An object's **apparent magnitude** (H + geometry) determines which bin's map is
consulted. Objects outside all bins get P=0 for all populations.

---

## Velocity grid

- v_λ, v_β ∈ [−0.8, 0.8] deg/day
- Step: 0.01 deg/day → 161 × 161 = 25,921 grid points
- Lookup: bilinear interpolation; out-of-grid → 0

---

## Cloning strategy

**Active path for all 4 populations:** conditional K|M cloner
(`clone_population_conditional_K_from_M_with_skycut`).

Key design:
- Sky cut is applied *before* cloning (to learn the correct M–K anti-correlation)
- Jointly samples M (mean anomaly at obstime) and K = Ω + ω from their empirical
  conditional distribution
- Preserves the M–K anti-correlation that keeps mean longitude L ≈ patch longitude
- This is critical for Trojans and TNOs where standard uniform-M cloning produces
  wrong velocity distributions (U-shaped M histogram)

**Dead code path:** uniform mean-anomaly cloner (`use_conditional_cloner=False`),
retained for backward compatibility; requires `nsc.clone_population_uniform_mean_anomaly`.

---

## Density estimation

Bayesian kNN posterior (Appendix B 2D) from the original notebook:
- k=5 nearest neighbors (k_map=5; previously 10, reduced to suppress MBA artifacts)
- d0 grid: 400 log-spaced points
- Posterior mean of 1/(π d0²) used as density estimate
- Per-population `nn_threshold_deg_per_day` masking at 0.2 deg/day for plots

---

## Ecliptic rotation: critical bug fix

**Do NOT use `astropy GeocentricTrueEcliptic` frame transforms for rates.**
Astropy drops differentials in that frame transform, giving wrong v_λ/v_β.

The fix (carried through the entire pipeline): manual rotation via obliquity matrix.
```python
eps = np.deg2rad(23.439291)
c, s = np.cos(eps), np.sin(eps)
# equatorial -> ecliptic:
x_ecl  =  x_eq
y_ecl  =  c * y_eq + s * z_eq
z_ecl  = -s * y_eq + c * z_eq
# same rotation applied to velocity components
```

`radec_rates_to_ecliptic_rates_manual()` is the public function for this.
`elements_to_vlam_vbeta_d_h()` (legacy uniform-cloner path) still calls the broken
`nsc.radec_rates_to_ecliptic_rates_at_obstime` — do not use it on the active path.

---

## NPZ file format

One file per obstime + sky center. Arrays stored:

```
x_grid, y_grid                            # 1D coordinate axes (deg/day)
obstime_str, center_lon_deg, center_lat_deg, max_sep_deg, center_label
population_names, mag_bin_labels, mag_bin_mins, mag_bin_maxs
density_raw__{POP}__{LABEL}               # (ny, nx) float32 downweighted density
nearest_dist__{POP}__{LABEL}              # (ny, nx) float32 nearest-clone distance
magcut_count__{POP}__{LABEL}              # int64 scalar
```

---

## ROC curve workflow

```python
# Step 1: score all 4 populations, tag true label, save
import velocity_density_pipeline as vdp
import pandas as pd

pms = vdp.ProbMapSet.from_npz("prob_maps_2025-03-21.npz")

scored_dfs = []
for pop_label, s3m_key in [("NEO","neo"), ("MBA","mba"), ("TNO","tno"), ("Trojans","trojan")]:
    df, scorer = vdp.load_s3m_population(s3m_key)
    probs, visible_df = pms.score_orbital_df(df, scorer=scorer, chunk=500_000, apply_mask=False)
    visible_df["true_population"] = pop_label
    scored_dfs.append(visible_df)

all_scored = pd.concat(scored_dfs, ignore_index=True)
all_scored.to_parquet("s3m_scored_2025-03-21.parquet", index=False)

# Step 2: ROC (can run from any notebook after parquet exists)
import numpy as np, pandas as pd, matplotlib.pyplot as plt

df = pd.read_parquet("s3m_scored_2025-03-21.parquet")
is_neo = df["true_population"] == "NEO"
p_neo  = df["P_NEO"].to_numpy()

thresholds = np.linspace(0, 1, 1001)
TPR = np.array([(p_neo[is_neo]  > X).sum() / is_neo.sum()  for X in thresholds])
FPR = np.array([(p_neo[~is_neo] > X).sum() / (~is_neo).sum() for X in thresholds])
purity = np.array([is_neo[p_neo > X].sum() / max((p_neo > X).sum(), 1) for X in thresholds])

# optimal threshold by max F1
f1 = 2 * (purity * TPR) / np.maximum(purity + TPR, 1e-9)
X_star = thresholds[np.argmax(f1)]
```

---

## Common gotchas

| Symptom | Cause | Fix |
|---------|-------|-----|
| TPR starts at ~0.15 at X=0 on ROC | `apply_mask=True` zeroes out sparse NEO regions | Use `apply_mask=False` in `score_orbital_df` |
| `NameError: scorer_neo not defined` | Copied `score_observation` example that had stale `scorer=` param | `score_observation` doesn't need a scorer |
| `FileNotFoundError: s3m_scored.parquet` | Ran ROC cell before saving parquet | Run the scoring/saving cell first |
| `AttributeError: jax has no attribute version` | JAX/jaxlib version mismatch triggered by `%autoreload 2` | Restart kernel cleanly; no autoreload |
| Wrong v_λ/v_β | Using astropy `GeocentricTrueEcliptic` for rate transforms | Use `radec_rates_to_ecliptic_rates_manual()` |

---

## Environment notes

- Python 3.11, conda env `lsdb_latest`
- Key packages: numpy, pandas, astropy, scipy, astroML, tqdm, matplotlib, pyarrow
- `neoscore` and `s3m_loader` are local packages in `~/digest3/`
- `neoscore` imports `adam_core` which uses JAX — sensitive to version mismatches
- Observed S3M obstime: `"2025-03-21T00:00:00"` (TDB scale)
- Default sky center: ecliptic lon=180°, lat=0°, radius=30°
