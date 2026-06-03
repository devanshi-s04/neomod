# NEO-only GMM cloning context

Date: 2026-05-31

This file summarizes the work in `cloning_test_ZI_v2.ipynb` and the intended next step for integrating it into the VDP probability-map pipeline.

## Goal

The current production VDP map generation uses conditional K|M cloning. That works, but the NEO velocity support can be too sparse or too narrow in parts of projected velocity space. The goal of this notebook is to test a **NEO-only Gaussian Mixture Model cloner** that can replace the current NEO clone generator while leaving MBA/TNO/Trojan cloning unchanged.

In other words:

- Use GMM cloning only for the NEO population.
- Keep the existing conditional K|M path for all non-NEO populations.
- Project both clone sets through the same VDP visibility/projection path.
- Compare orbital distributions, visible velocity support, and density maps before touching production.

## Files read / context used

- `cloning_test_ZI.ipynb`: advisor-edited notebook with an initial GMM prototype.
- `TODAY_CONTEXT_DIGEST2_VDP.md`: context for digest2 vs VDP, Sorcha scoring, and why NEO velocity support matters.
- `src/velocity_density_pipeline.py`: production VDP map generation and cloning functions.

Relevant production code:

- `clone_population_conditional_K_from_M(...)`
- `clone_population_conditional_K_from_M_with_skycut(...)`
- map generation loop that builds visible clones and density maps
- NEO-only smoothing behavior already exists via `DEFAULT_SMOOTH_POPULATION_NAMES = ("NEO",)`

## What the advisor prototype did

The original advisor prototype fit a `GaussianMixture` to transformed NEO orbital features. It included angular variables using sine/cosine transforms, but it was not yet production-ready because:

- It included `H` directly in the GMM feature vector.
- It sampled some orbital quantities in ways that could become physically inconsistent.
- It returned raw sampled arrays rather than a clean VDP-ready orbital-element DataFrame.
- It was a proof of concept, not wired into the existing VDP projection diagnostics.

## What `cloning_test_ZI_v2.ipynb` adds

The new notebook turns the prototype into a testable NEO-only cloner.

The GMM feature vector is now:

```text
log(a),
q,
sin(i), cos(i),
sin(node), cos(node),
sin(argperi), cos(argperi),
sin(M_obs), cos(M_obs)
```

where:

- `q = a(1-e)` is used as the NEO-defining perihelion coordinate.
- `e` is reconstructed from `e = 1 - q/a`, rather than sampled independently.
- `M_obs` is the mean anomaly at the map observation epoch.
- `t_p` is reconstructed from `M_obs` and `a`.
- angular variables are sampled through sine/cosine pairs to avoid artificial edges at 0/360 degrees.

The clone output is a normal VDP-compatible DataFrame:

```text
a, e, i, node, argperi, t_p, H, q, M_obs_deg
```

## Important latest change: H is not in the GMM

Zeljko suggested not using `H` in the GMM because the GMM smooths the observed absolute-magnitude distribution and does not preserve the count increase toward the faint end.

This is now fixed.

Current behavior:

- `H` is removed from the GMM feature matrix.
- After orbital clones are sampled and validated, each clone gets an `H` value drawn with replacement from the empirical/source NEO `H` distribution.
- Diagnostics report:

```text
H_sampling: empirical_source_distribution
```

Why this matters:

- The GMM now learns orbital/phase-space structure only.
- The observed faint-end `H` pileup is preserved directly.
- The 1D `H` sanity plot should now match the source distribution much more closely, especially near the faint end.

## Current qualitative result

From the latest notebook plots:

- The `H` distribution is better after empirical sampling. The faint-end rise is preserved instead of smoothed away by the GMM.
- The GMM still smooths some details in `q`, especially the source spike near `q ~ 1.1`, but Zeljko suggested leaving that alone for now.
- The GMM visible velocity support is broader than the current K|M visible clone support, especially into the more negative `v_lambda` region.
- The GMM raw and smoothed density maps look broadly similar to K|M, but with somewhat broader NEO support.

Interpretation:

- This is promising for recovering NEOs that VDP currently misses because NEO support is too sparse.
- It is not proven better until it is inserted into map generation and tested on the Sorcha/digest2 comparison parquet.
- The main risk is that broader NEO density support could increase MBA contamination, so the full ROC/F1 comparison is the deciding test.

## What is similar to the advisor notebook

- Still uses `sklearn.mixture.GaussianMixture`.
- Still uses transformed angular features.
- Still tests source NEOs vs GMM clones with 1D/2D plots.
- Still focuses on cloning the NEO population, not replacing the full population model.

## What changed from the advisor notebook

- GMM output is converted into a physically consistent orbital-element table.
- `q` and `a` are sampled, then `e` is derived from them.
- `M_obs` is sampled and `t_p` is reconstructed at the observation epoch.
- `H` is now sampled empirically outside the GMM.
- Invalid sampled clones are rejected.
- The clones are projected through `vdp.build_visible_subset_dataframe(...)`.
- The notebook compares K|M vs GMM in VDP velocity space and density-map space.

## Production integration idea

Add a production option like:

```python
cloner = "conditional_km"  # existing default
cloner = "neo_gmm"         # new NEO-only option
```

Then in map generation:

- If `population_name == "NEO"` and `cloner == "neo_gmm"`, use the GMM cloner.
- Otherwise use the existing conditional K|M cloner.

The production GMM cloner should:

- fit only on the NEO source table for the current map epoch / source subset,
- use the same feature representation as `cloning_test_ZI_v2.ipynb`,
- sample `H` empirically from source NEOs,
- return the same columns expected by the existing VDP projection functions.

## Suggested next test on Hyak

1. Copy `cloning_test_ZI_v2.ipynb` and this context file to Hyak.
2. Implement a NEO-only GMM cloner in `src/velocity_density_pipeline.py` or a small side module first.
3. Generate a small test probability map with:
   - NEO = GMM clones
   - MBA/TNO/Trojan = existing conditional K|M clones
4. Score the existing Sorcha comparison sample.
5. Compare against the current best VDP run:
   - VDP wide-grid result: F1 about 0.740, completeness about 67.5%, contamination about 18.0%.
   - digest2 reference: F1 about 0.836, completeness about 77.3%, contamination about 9.0%.

The useful question is not only "does GMM fill more NEO support?" It is:

```text
Does NEO-GMM increase completeness without giving back too much contamination?
```

## Quick scp targets

From the local machine, likely copy:

```bash
scp neomod/cloning_test_ZI_v2.ipynb neomod/cloning_gmm_neo.md ds2004@klone.hyak.uw.edu:/mmfs1/gscratch/astro/ds2004/sorcha/neomod/
```

Adjust the remote path if the Hyak repo lives somewhere else.
