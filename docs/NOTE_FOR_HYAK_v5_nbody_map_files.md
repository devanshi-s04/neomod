# Note for Hyak — need 6 n-body map files for the v5-normalisation notebook

**From:** Arnor Claude · **Date:** 2026-07-19 · **Re:** `nbody_benchmark_v5_normalisation_s3m.ipynb`

Built the n-body copy of `benchmark_v5_normalisation_s3m.ipynb`. The **ROC sections already work** —
they only need `benchmark_comparison_s3m_nbody.parquet`, which I have — and the full-sky number
reproduces your earlier note exactly: **AUC 0.878, best-F1 0.662** (compl 57.5%, contam 22.0%,
digest2 F1 0.725 for comparison). Per-elongation-bin table also ran clean (attached below).

## What's blocked: the raw map-image sections

This notebook also renders the actual probability-density maps (not just ROC), and those cells call
`vdp.ProbMapSet.from_npz(...)` directly on the `.npz` files — not the scored parquet. I only have the
**two-body** maps locally (`prob_maps_grid_s3m/`, 667 files, 20 GB); the n-body maps
(`prob_maps_grid_s3m_nbody/`) only exist on Hyak.

I don't need all 667 — the notebook only ever loads **6 specific centers** (all `lat+00`):

```
prob_maps_grid_s3m_nbody/prob_maps_grid_dlon+000_lat+00.npz
prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-030_lat+00.npz
prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-040_lat+00.npz
prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-050_lat+00.npz
prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-090_lat+00.npz
prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-120_lat+00.npz
```

At ~31 MB/file (matching the two-body size) that's **~185 MB total** — one scp, not the full grid.

## The ask

```bash
scp /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/.../prob_maps_grid_s3m_nbody/prob_maps_grid_dlon+000_lat+00.npz \
    /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/.../prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-030_lat+00.npz \
    /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/.../prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-040_lat+00.npz \
    /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/.../prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-050_lat+00.npz \
    /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/.../prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-090_lat+00.npz \
    /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/.../prob_maps_grid_s3m_nbody/prob_maps_grid_dlon-120_lat+00.npz \
    ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/prob_maps_grid_s3m_nbody/
```
(fix the source path prefix to wherever `prob_maps_grid_s3m_nbody/` actually lives — I don't have that
path on my side to fill in exactly.)

## What each file unlocks

| file | section |
|---|---|
| `dlon+000` (antisun) | §1 support-mask on/off maps; §1b was dropped (map-method comparison, not n-body relevant) |
| `dlon+000, -040, -090, -120` | §5/§10 direction-comparison grids (density + P(NEO), all 4 pops) |
| `dlon+000, -030, -040, -050, -090, -120` | §7 mega-grid across elongation, mag 24–25 |

## What's already done, no blocker

- Section 3 (ROC/calibration) — reproduces your AUC 0.878 / F1 0.662 exactly
- Section 6 (ROC by elongation, 6-panel figure) — rendered, saved to `Figures/nbody_benchmark_roc_by_elongation.png`
- Section 8 (ROC grid, elongation × magnitude) — rendered
- Section 9 (rate-space tracklet distribution) — needs only the parquet, rendered, saved

## Per-elongation-bin ROC table (n-body, already have this)

```
             Bin  N_NEO  N_total  VDP F1  VDP compl  VDP contam  d2 F1  d2 compl  d2 contam
        Full sky   1051   215645   0.662      57.5%        22.0%  0.725     63.5%      15.6%
   |Δλ|=0--20°       125    38557   0.721      74.4%        30.1%  0.804     72.0%       9.1%
  |Δλ|=20--40°       116    33700   0.575      52.6%        36.5%  0.717     61.2%      13.4%
  |Δλ|=40--70°       195    44910   0.565      47.7%        30.6%  0.767     65.1%       6.6%
 |Δλ|=70--110°       281    55054   0.632      52.7%        20.9%  0.685     53.0%       3.2%
|Δλ|=110--141°       334    43424   0.760      65.0%         8.4%  0.762     67.7%      12.7%
```

digest2 leads VDP at every elongation except 110–141° (where they're essentially tied, 0.760 vs 0.762)
— consistent with the two-body pattern, so the integrator fix didn't flip the VDP-vs-digest2 ranking,
just corrected the underlying positions/scores it's computed from.
