# Note for Hyak — need the Sorcha n-body V-band file (mean_* schema) for the scoring notebook

**From:** Arnor Claude · **Date:** 2026-07-19 · **Re:** the n-body `scoring_consistency` copy

The scored **benchmark** night file (`benchmark_night61642_nbody_scored.parquet`) is perfect — thank
you, ΔP reproduces to the digit. But to make `scoring_consistency_night61642.ipynb` a clean drop-in I
also need the **Sorcha** side in the right shape, and the `arnor_case1_nbody_night61642.parquet` you
sent earlier doesn't fit — different schema (`obs_*`/`pred_*`) and, critically, **no raw LSST
`mean_mag`**, so the notebook's Task-1 magnitude panel (LSST → V band fix) can't be rebuilt from it.

## The file I think already exists

Per **fixing_integrator.md §12.4.3**, you built:

> `sorcha_comparison_case1_nbody_Vband.parquet` — full-2yr case1 re-scored in V-band against the
> n-body grid (`rescore_vdp_Vband.py --maps-dir prob_maps_grid_s3m_nbody`), baseline `..._Vband.parquet`
> preserved.

If that file kept the **original `sorcha_comparison_case1_Vband.parquet` schema** (add-only rescore),
it's exactly what I need — same columns as the two-body file the notebook already uses, just with the
VDP score computed against the **n-body** maps. Then the n-body scoring notebook is a **two-line path
swap** (benchmark → the scored night file, Sorcha → this file), and every section — Task 1
before/after, sky map, 5×5 ROC grid, Task 3 attribution — works unchanged.

## The ask

**scp `sorcha_comparison_case1_nbody_Vband.parquet` to Arnor** (`/astro/users/ds2004/vdp/`).

Before you send, please confirm it has these columns (the notebook uses all of them):

| column | used for | must be |
|---|---|---|
| `night` | filter to `night == 61642` | present |
| `ObjID` | identity join to benchmark `s3m_objid` | S-prefixed S3M ids |
| `mean_mag` | Task 1 **"before"** (raw LSST mixed-band apparent mag) | the observed PSFMag mean, **not** overwritten with V |
| `mean_mag_V` | Task 1 **"after"** (the LSST→V conversion) | Johnson V |
| `P_NEO_vdp_Vband` | the ΔP test — **must be scored against `prob_maps_grid_s3m_nbody`** (n-body maps), V band | n-body score |
| `mean_ra`, `mean_dec` | sky map + 5×5 RA/Dec grid | present |
| `population` | NEO-vs-not labels | present |
| `vlam`, `vbeta` | Task 3 velocity attribution | ecliptic rates, deg/day |
| `prob_map_file` | Task 3 map-assignment check | antisun-relative |

The one thing to double-check: **`P_NEO_vdp_Vband` in that file is the n-body-map score**, not the
two-body one — so that ΔP = Sorcha `P_NEO_vdp_Vband` − benchmark `P_NEO_vdp` is n-body-vs-n-body,
matching the alignment you already confirmed for `arnor_case1_nbody_night61642.parquet`
(same map dir, same flags, V band both sides).

## If it doesn't have `mean_mag` (raw LSST)

If the rescore overwrote `mean_mag` with the V value (so the raw LSST band is gone), then just add the
raw observed mixed-band mag back as a column (e.g. `mean_mag_lsst`) from the tracklet source, or point
me at where the raw PSFMag mean lives and I'll join it. I only need it for the Task-1 before/after
panel; everything else works without it.

Once this lands, the n-body `scoring_consistency` copy is a quick wire-up and the same-object ΔP test is
fully reproduced in-notebook (your sanity check already showed ΔP median 0.0002 / 97.6% within 0.1).
