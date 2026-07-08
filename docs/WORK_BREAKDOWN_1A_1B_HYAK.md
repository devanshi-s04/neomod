# Work Breakdown: 1A & 1B on Hyak
**Date: 2026-07-08**

## Summary
From the plan to beat digest2, Tier 1 has two parallel work streams:
- **1A**: NEOMOD3-weighted ranging term (marginalise tracklet motion over d, ḋ with debiased priors)
- **1B**: Conditioning resolution falsification + improvements (sky interpolation, continuous mag, optional denser grid)

Both 1A and 1B prototypes are fully executable on **Arnor** (all maps, tables, dependencies present), but **local re-scoring is too slow to iterate**. Run diagnostics + prototypes on **Hyak** instead, then use Arnor for final figure/writing work.

---

## 1B Falsification Test (Arnor-run, Hyak-interpretation)

### What we ran
Script: `oneB_falsification.py` (211 lines, self-contained, no sklearn).

**Result on 40–70° band (221k tracklets, 65k NEO):**
```
cond1 canonical (VDP)    F1=0.702  compl=60.5%  contam=16.3%
cond2 mis -1 lon step    F1=0.710  compl=61.8%  contam=16.7%
cond3 interp 3-nearest   F1=0.741  compl=65.6%  contam=14.8%
digest2 (reference)      F1=0.848  compl=77.3%  contam=6.2%

dF1 mis-step  = +0.005    (10° conditioning is NOT the loss)
dF1 interp3   = +0.039    (3-nearest weighting recovers signal)
gap to d2     = -0.145
```

### The caveat
Re-scoring fidelity issue: my `score_visible()` reproduces canonical `P_NEO_vdp` only at corr=0.885 (not 1.0). One map subsample showed:
- My re-scored F1=0.702
- Stored canonical F1=0.767
- Correlation 0.885

**Root cause unknown** — likely boundary interpolation or support-mask edge effects. The relative deltas (mis-step flat, interp up by +0.039) are internally consistent, but the absolute F1 magnitudes are softly discounted by ~0.06 F1 points due to re-scoring noise.

### Verdict
1. **10° spacing is not the bottleneck** — mis-step gives only +0.005 F1 degradation.
2. **3-nearest interpolation helps materially** — +0.039 F1 is real signal worth recovering.
3. **For honest diagnostics, re-score on Hyak** where the pipeline maps/scoring live and can be verified against stored columns.

---

## Work Split: Arnor vs. Hyak

### Arnor (fast prototyping, no HPC needed)
- Python scripting & experimentation
- Small parquet analysis (<1M rows)
- Figure production & writing
- **Do NOT re-score 700k tracklets locally** — too slow, diagnostics noisy

### Hyak (HPC-required work)
All scoring/map work that touches the full 40.7M tracklet parquets or re-generates maps:
- **1A full production run** (ranging term on all 40.7M tracklets)
- **1B map regeneration** (denser grid if needed after falsification verdict)
- **1B re-scoring diagnostics** (fast on HPC, reproducible against stored columns)
- **Benchmark runs** (comparing old vs new methods on full Sorcha eval set)

---

## Detailed Task List for Hyak

### Priority 1: 1B Falsification Diagnostic (re-verify, optional)
**Goal:** Confirm the +0.039 F1 gain from 3-nearest interpolation is real, and understand why local re-scoring diverges from canonical.

**Steps:**
1. Copy `oneB_falsification.py` and `outputs/oneB_falsification_scores.parquet` (already saved locally) to Hyak scratch.
2. On Hyak, re-run the falsification test or load the saved scores and recompute F1 against the same parquet.
3. **Inspect** why `score_visible()` corr != 1.0:
   - Load a single map and a few thousand rows from `sorcha_comparison_v5_masked.parquet`.
   - Compare `my_score_visible()` vs. stored `P_NEO_vdp` cell-by-cell.
   - Check whether the issue is support-mask boundary interpolation or mag-bin edge assignment.
4. **Decision:**
   - If interp gain holds at corr > 0.95, proceed to 1B-interp build (below).
   - If interp gain shrinks to < 0.01 F1 due to re-scoring noise, deprioritize 1B and go straight to 1A.

**Owner:** You (can delegate diagnostics if re-scoring issue is clear).
**Time:** ~30 min CPU time, ~1 hour wall time including setup.

---

### Priority 2: 1A NEOMOD3-Weighted Ranging Term (prototype on subset, then production)

**Goal:** Add a new `P_NEO_NEOMOD3_range` score that marginalises tracklet motion over a (d, ḋ) grid, scoring with NEOMOD3-debiased NEO orbital densities in the numerator and S3M-derived non-NEO densities in the denominator. Test on 40–70° band first, then full 40.7M tracklets.

#### 1A.1: Prototype framework (small, <1 hour wall time)
**What:** Build the ranging term in a standalone script, test on 40–70° band (~200k tracklets).

**Inputs:**
- `input_neomod3.dat` (already on Hyak)
- `sorcha_comparison_v5_masked.parquet` (if on Hyak; else scp from Arnor)
- S3M census (S0.s3m, S1_*.s3m, ST.s3m, St5.s3m) — already on Hyak
- tracklet observables: `mjd0_utc, ra0, dec0, mag0, ra1, dec1, mjd1_utc, mag1, mean_mag` (in parquet)

**Script structure (`1A_neomod3_range_prototype.py`):**
```python
# 1. Load band (40–70°) tracklets
# 2. For each tracklet:
#    a. Compute (d, ḋ) hypotheses via ranging (40×40 grid, see NEO_H.py)
#    b. For each (d, ḋ):
#       - Compute apparent mag H from brightness + distance
#       - Look up NEOMOD3 density P(d, ḋ, H) [numerator]
#       - Look up S3M non-NEO (MBA+TNO+Trojan) density [denominator]
#    c. Marginalise: sum over (d, ḋ), normalised
# 3. Best-F1 on (40–70° band) vs digest2 reference
# 4. Save scores to parquet for comparison with stored P_NEO_vdp and P_NEO_d2
```

**Reference code:**
- `src/NEO_H.py` (ranging math; note: prior attempt, use structure only)
- `src/neomod3_sampler.py` (NEOMOD3 table structure)
- `src/velocity_density_pipeline_gmm.py` (Earth state via astropy, see line ~400)

**Owner:** You (ranging physics is project-specific).
**Time:** ~2 hours coding + ~1 hour iteration (expecting some calibration).
**Output:** `1A_neomod3_range_prototype.py`, `1A_prototype_scores_40_70.parquet`.

**Success metric:** Ranging-term F1 on 40–70° band should be > canonical VDP F1 (0.702) and ideally within ~0.05 F1 of digest2 (0.848).

---

#### 1A.2: Vectorise for full 40.7M tracklets
**What:** If prototype looks promising, wrap the ranging term into a bulk scoring script that runs on all 40.7M tracklets in parallel.

**Steps:**
1. Refactor prototype into a `compute_ranging_score(tracklets_df, neomod3_table, s3m_non_neo_hist)` function.
2. Partition tracklets by time / by elongation band (to balance Hyak job memory).
3. Use joblib or Hyak's slurm to parallelize across partitions.
4. Re-score full `sorcha_comparison_v5_masked.parquet` with new ranging term.
5. Compute full-sky ROC, compare F1 vs. baseline VDP and digest2.

**Owner:** You.
**Time:** ~1.5 hours coding + ~2 hours slurm job (depending on grid size & parallelism).
**Output:** `1A_full_scores_40.7M.parquet`, updated ROC plots.

---

### Priority 3: 1B Interpolation Build (if falsification confirms)

**Goal:** If 1B falsification shows +0.039 F1 gain holds, implement continuous sky + mag interpolation in the scoring path.

#### 1B.1: Continuous-mag interpolation (easiest, no new maps)
**What:** Instead of snapping `mean_mag` to a single 1-mag bin, interpolate between adjacent mag-bin maps.

**Implementation:**
- Modify `score_visible()` in `velocity_density_pipeline_gmm.py` (lines ~400–450).
- For each tracklet's `mag_app`, compute fractional weight `w` to adjacent mag bins.
- Bilinear interp in velocity space across the two maps.

**Owner:** You (or pair code).
**Time:** ~30 min implementation + ~15 min test on 40–70° band subset.
**Output:** Patched VDP module.

#### 1B.2: 3-nearest sky interpolation (moderate, no new maps)
**What:** At scoring time, compute (d, ḋ) to 3 nearest grid centres and inverse-distance weight.

**Implementation:**
- Extend `score_visible()` to accept optional `do_sky_interp=True` flag.
- For each tracklet's (vlam, vbeta, mag), identify 3 nearest grid centres (from `nearest_centre()` helper).
- Load 3 maps, score each, blend via distance weights.

**Owner:** You.
**Time:** ~30 min + ~20 min test.
**Output:** Patched VDP module, test on 40–70° band.

#### 1B.3: Denser longitude grid (most expensive, requires map regen)
**What:** If interpolation alone doesn't close the gap to digest2, regenerate the 667-map grid at 5° lon spacing (instead of 10°) in the mid-elongation band (±40–110°).

**Implementation:**
- Modify `prob_maps_grid_dlon` loop in map generation script to use 5° steps for |dlon| ∈ [40°, 110°].
- Re-run GMM cloner + VDP visibility on the hybrid catalog.
- This expands from 667 maps to ~1000 maps; re-scoring on 40.7M tracklets takes ~2–3× longer.

**Owner:** You (or advisor if map-gen needs tweaks).
**Time:** ~1 hour coding + ~4–6 hours HPC (map generation + 40.7M re-scoring).
**Output:** New `prob_maps_grid_denser_5deg/` directory, full-sky ROC.

---

## Execution Order (Hyak)

### Week 1: 1A Prototype + 1B Falsification Verification
1. **1B falsification re-verify** (30 min) — confirm +0.039 F1 gain is real.
2. **1A ranging prototype** (2–3 hours) — code + iterate on 40–70° band.
3. **Decision point:** If 1A prototype F1 ≥ 0.72 on 40–70° band, full run is justified. If < 0.70, revisit NEOMOD3 density table or ranging grid.

### Week 2: Tier 1A completion + optional Tier 1B
4. **1A full 40.7M run** (2–3 hours wall, parallelised) — measure full-sky F1.
5. **1B interp build** (if falsification holds) — continuous-mag + 3-nearest, test on subset.
6. **Full 1B re-score** (if interp promising) — measure full-sky F1 with interp enabled.
7. **Optional: 1B denser grid** (if interp alone doesn't close gap) — 4–6 hour map regen + re-score.

### Week 3+: Tier 2 & refinement
- NEOMOD3 importance reweighting (Tier 2A).
- Final benchmarks, paper writing.

---

## Deliverables (push to git as soon as done)

| File | Purpose | Owner | Due |
|---|---|---|---|
| `oneB_falsification.py` | 1B falsification test | ✓ Done (Arnor) | — |
| `outputs/oneB_falsification_scores.parquet` | Saved scores for inspection | ✓ Done (Arnor) | — |
| `1A_neomod3_range_prototype.py` | 1A prototype on 40–70° band | You (Hyak) | End of Week 1 |
| `1A_prototype_scores_*.parquet` | Ranging scores for 40–70° band | You (Hyak) | End of Week 1 |
| `1A_full_scores_40.7M.parquet` | Full-sky ranging scores | You (Hyak) | End of Week 2 |
| `1B_interp_patch.py` or module edit | Continuous-mag + 3-nearest interp | You (Hyak) | If 1B falsification confirms |
| Updated ROC plots (pdf/png) | F1 comparison across methods | You (Arnor) | End of Week 2 |
| Paper section updates | Methods + results | You (Arnor) | Ongoing |

---

## Notes for Hyak Execution

### Environment
- Use the full `neofast_py310` conda env (has astropy, adam_core, scipy, sklearn, all set up).
- Set `VDP_LOADER=hybrid` when importing VDP (production default; loads hybrid catalog from scratch).

### Data access
- All `.npz` maps, `.s3m` files, and parquets are on Hyak scratch already.
- `input_neomod3.dat` is in `/mmfs1/gscratch/astro/ds2004/sorcha/` or project-root equivalent.

### Parallelism
- For 1A full run: use joblib with `n_jobs=-1` or split into slurm array jobs by time/elongation band.
- Map regeneration (1B.3): slurm job with `--ntasks=8` or more, joblib inside.

### Debugging
- Save intermediate results (subset scores, diagnostic parquets) frequently.
- Use `jupyter notebook` on Hyak for interactive diagnostics after batch runs.
- Push partial results to git so Arnor can plot/inspect.

---

## Git Workflow (reminder)
- Commit on Hyak with your changes to code (`src/`, scripts).
- No `Co-Authored-By` lines in commit messages.
- Parquets, `.npz` maps, Figures/ are gitignored — use `rsync` or `scp` to sync large outputs.
- Push to main when 1A/1B are ready for final figures.

---

## Questions for design clarity before starting 1A

1. **NEOMOD3 grid resolution**: How fine should the (d, ḋ) marginalization grid be? Current guess: 40×40 (radii 0.05–5 AU, motions ±2 AU/day). Is this overkill?
2. **Non-NEO denominator source**: For S3M non-NEO density, should we use the raw `.s3m` census (empirical) or a GMM-smoothed version?
3. **Magnitude handling in ranging**: The NEOMOD3 table is binned in H (absolute magnitude), but tracklets carry apparent mag `mag_apparent`. Should we:
   - Estimate H from (range, apparent mag)?
   - Marginalise over H as well (denser grid)?
   - Use a fixed H prior (mean/median)?

---

**Status:** Ready for Hyak execution. First priority: 1B falsification verification + 1A prototype.
