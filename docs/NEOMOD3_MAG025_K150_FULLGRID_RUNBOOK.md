# NEOMOD3 0.25-mag / k=150 full-grid build — living runbook

Status: **GATE PASSED (22/22). BLOCKED ON THE NEO SOURCE, NOT ON STORAGE.** The 667-center build
is not launched: the specified HIGH NEO realization does not exist all-sky. See §13.2.

Last updated: 2026-08-16.

---

## 1. Scientific objective and frozen configuration

Build NEO/MBA/TNO/Trojan velocity-density maps and posteriors on the existing 667 sky centers at
**0.25-mag apparent-V resolution**, with a larger NEO Monte Carlo sample and a larger kNN neighbour
count for NEO, then evaluate on a new, larger held-out realization (TEST2).

| item | frozen value |
|---|---|
| sky centers | 667 existing (`outputs/geometric_density_estimator_ablation/frozen_center_list.txt`) |
| epoch / observer geometry | `2027-08-25T00:00:00`, unchanged from `MAP_BUILD_SEAL.json` |
| patch radius | 30 deg |
| velocity domain | [-5, +5] deg/day in `v_lambda`, `v_beta` |
| grid step | 0.01 deg/day → 1001 x 1001 |
| magnitude quantity | **apparent V** only, for map sources and scored tracklets |
| magnitude bins | 0.25 mag wide, `14.00 <= V < 25.00` → **44 bins** |
| bin semantics | half-open `[lo, lo+0.25)`; V=23.4 → `[23.25, 23.50)` only |
| NEO source | HIGH NEOMOD3 GEN realization (see §2, and the caveat in §13.2) |
| NEO k | **150** |
| MBA / TNO / Trojan k | **10** |
| Gaussian smoothing | **OFF for every population** |
| posterior | `P(c) = rho_c / sum_all rho`, per defined pixel |
| normalization | physical weights preserved; larger NEO MC must not raise physical NEO abundance |
| calibration | none. No Platt scaling, no density renormalization in this task |

Density integrals are recorded per (center, population, bin) for a later normalization study. They
are recorded only; nothing is rescaled by them here.

### 1.1 Apparent-V confirmation

No conversion is required for map sources. Both magnitude paths already produce apparent V with the
HG phase function (G = 0.15):

- `compute_apparent_magnitude_for_population` — docstring: *"Compute apparent V magnitude for every
  object in df at obstime ... using the HG phase function (G=0.15)"* (non-NEO path).
- `build_visible_subset_dataframe` — computes `mag_app` from the same n-body state
  (`# ---- apparent V mag from THIS SAME state (fixing_integrator.md §9.8) ----`), which is the
  column the NEOMOD3 projection cache stores (NEO path).

So `mag_app` **is** apparent V in both sources. Verified in the gate (§9).

---

## 2. Exact inputs

| role | path | hash / note |
|---|---|---|
| sealed VDP module | `neomod/src/velocity_density_pipeline_neomod_clone_only.py` | sha256 `a6de18c2197cbcb9f93014712aaad232272578150e2222f28eae723b339f923a` |
| map seal (v1) | `outputs/splits/MAP_BUILD_SEAL.json` | sha256 `c30cc29d24ead70da1111782092b37557f21d6b90f581609bf2490907a474960` |
| center list | `outputs/geometric_density_estimator_ablation/frozen_center_list.txt` | 667 centers |
| non-NEO source | `outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet` | 14,380,436 rows |
| non-NEO split manifest | `outputs/splits/nonneo_split_manifest.parquet` | GEN/CAL/TEST by parent identity |
| split provenance (1-mag) | `outputs/splits/split_provenance.json` | 8 bins only — see §13.3 |
| split provenance (0.25-mag) | `outputs/splits/split_provenance_mag025.json` | **new**, built by `build_mag025_split_provenance.py` |
| NEO base cache (1e8 draws) | `outputs/neomod3_projection_cache/` | `effective_factor_NEO` 8.746673464103498 |
| NEO HIGH extra draws | seeds 1,000,000 + 0..319, 2e6 draws each = 6.4e8 | recorded in `outputs/more_neomod_samples_knn/high_draw_manifest.json` |
| HIGH single-cell source | `outputs/more_neomod_samples_knn/source_high.parquet` | sha256 `40490b3bc4ffaec122919981396168299c1e84a384dd345c46f8a7adb20fc297`, 207,059 rows — **one center, one bin only** |

Output root (never overwrites anything previous):

    outputs/neomod3_mag025_k150_maps_v2/

---

## 3. Magnitude selection and bin-edge semantics

44 bins, `lo_i = 14.00 + 0.25 i` for `i = 0..43`, membership `lo_i <= V < lo_i + 0.25`.

- The top bin is `[24.75, 25.00)`. V = 25.00 exactly is **out of range**, not folded into the top bin.
- V = 23.4 selects **only** `[23.25, 23.50)`.
- A slice is never widened. If support is insufficient, §4 applies; the bin is never merged with a
  neighbour and never silently promoted to a 1-mag bin.

Implemented in `mag025_bins()` and asserted in `validate_neomod3_mag025_k150_maps.py`.

---

## 4. Insufficient-support rule (single, explicit, documented)

Requested k: NEO 150, MBA/TNO/Trojans 10.

For a (center, population, 0.25-mag bin) cell with `n` selected samples:

    n >= k_req + 1        -> VALID,   effective_k = k_req
    2 <= n <= k_req       -> INVALID, effective_k = None, reason = "insufficient_support"
    n < 2                 -> INVALID, effective_k = None, reason = "no_samples" / "below_min"

**The requested k is never reduced and the slice is never widened.** A cell that cannot support its
requested k produces `density = NaN` for that population, not 0. Pixels where any population density
is NaN are `P = NaN` and are counted as undefined; they are never coerced to zero.

Rationale for not falling back to a smaller effective k: the sealed path's `k_eff = min(k, n-1)`
would silently make a sparse bin use a different estimator than a dense bin, so a comparison across
bins would confound support with configuration. Refusing is auditable; degrading silently is not.

Every cell is recorded in the coverage table with `n`, `k_requested`, `k_effective`, `valid`,
`reason` — including valid ones — in `coverage_table.parquet`.

---

## 5. k by population

| population | k | source |
|---|---|---|
| NEO | 150 | HIGH NEOMOD3 GEN realization |
| MBA | 10 | S3M GEN split |
| TNO | 10 | S3M GEN split |
| Trojans | 10 | S3M GEN split |

Passed per population via `clone_sources[pop]["k_map"]`, which the sealed builder reads as
`k_eff = min(info.get("k_map", k_map), n_visible - 1)`. The builder asserts the estimator actually
received the requested k (spy wrapper), so a silent `min()` reduction cannot pass unnoticed.

---

## 6. Density estimator and the closed-form fast path

### 6.1 The identity

The sealed Bayesian kNN posterior for the local spacing scale `d0` is

    p(d0 | d_1..d_k)  ∝  prod_j (2 / (d0 Γ(j))) exp[-(d_j/d0)^2] (d_j/d0)^(2j-1)
                      =  C(d) · d0^(-k(k+1)) · exp(-S/d0^2),      S = sum_j d_j^2

because `sum_j (2j-1) = k^2` and there are `k` factors of `1/d0`. With the flat prior on `d0` that
the sealed code uses (it normalizes with `np.trapezoid(p, d0_grid)`, i.e. measure `dd0`), the
posterior mean of `n0 = 1/(pi d0^2)` has the **exact closed form**

    n0 = ( k(k+1)/2 - 1/2 ) / ( pi * sum_j d_j^2 )

This is not an approximation of the sealed estimator — it is the analytic value of the integral the
sealed estimator evaluates numerically.

### 6.2 Why it matters here

Measured against this closed form, the sealed quadrature's `n_d0_grid` default of 400 nodes is
converged at k=10 but not at larger k, because the integrand's peak narrows as `d0^-k(k+1)`:

| k | n_d0=400 | 2000 | 8000 |
|---|---|---|---|
| 10 | 7.1e-16 | 7.5e-16 | 6.6e-16 |
| 25 | 6.0e-12 | 8.5e-16 | 6.6e-16 |
| 50 | **2.6e-04** | 2.5e-15 | 1.2e-15 |
| 100 | **1.1e-02** | 6.0e-15 | 4.0e-15 |
| 150 | — | — | 8.4e-15 |
| 200 | — | — | 1.7e-14 |
| 250 | — | — | 1.8e-14 |

At NEO k=150 the production default would inject a k-dependent numerical error. The closed form has
no quadrature error at all, so it is both faster and more accurate. See §12 for the equivalence
check and the retain/reject decision.

---

## 7. Slurm commands (reproduce the build)

```bash
cd /mmfs1/gscratch/dirac/ds2004/sorcha

# 0. exact 0.25-mag split provenance (required; see §13.3)
sbatch neomod/pipeline/slurm/build_mag025_split_provenance.sbatch

# 1. single-center acceptance + timing gate
sbatch neomod/pipeline/slurm/gate_neomod3_mag025_k150.sbatch

# 2. full grid (restartable; skips completed hash-valid centers)
sbatch neomod/pipeline/slurm/build_neomod3_mag025_k150_array.sbatch

# 3. resume after preemption / partial failure — identical command, idempotent
sbatch neomod/pipeline/slurm/build_neomod3_mag025_k150_array.sbatch

# 4. full-build validation + seal
sbatch neomod/pipeline/slurm/validate_neomod3_mag025_k150.sbatch
```

---

## 8. Resources, parallelization, array design

- Partition `cpu-g2`: 68 nodes, **192 CPUs**, 1.5 TB RAM per node. Account `astro`.
- `MaxArraySize = 10001`, `MaxJobCount = 100000` → `--array=0-666` is legal in one array.
- No per-association `MaxSubmitJobs`/`MaxJobs` limit is set for this user (`sacctmgr` returns no
  association row), so concurrency is bounded by fair-share and free nodes, not policy. Recorded
  concurrency: `%MAX_CONCURRENT` default **100** (100 x 32 CPUs = 3200 CPUs ≈ 17 nodes).
- `n_jobs` is taken from `SLURM_CPUS_PER_TASK`, never hard-coded.
- Oversubscription guard, set in every sbatch:
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`.
- All work (including hashing, validation, notebook execution and "small" Python checks) runs under
  `sbatch`/`srun`. Nothing is computed on a login node.

### 8.1 GPU decision

**No GPU backend is used, because none exists in this environment.** Checked on a compute node:

| module | status |
|---|---|
| `cupy` | not installed |
| `cuml` | not installed |
| `torch` | not installed |
| `faiss` | not installed |
| `numba.cuda` | module imports (numba 0.65.1) but there is no kNN implementation |

Per the time-box rule, no new GPU kernel was written. Requesting a GPU while running NumPy/SciPy on
the CPU is forbidden and was not done. The optimized CPU path (§6, §12) is used instead.

---

## 9. Pre-build acceptance gate (one center)

`pipeline/gate_neomod3_mag025_k150.py`, run on `dlon+000_lat+00`:

1. 0.25-mag membership is exact; V=23.4 selects only `[23.25, 23.50)`, and the selected rows'
   V range is contained in that interval.
2. NEO estimator received k=150; MBA/TNO/Trojans received k=10 (spy assertion, not inspection).
3. `smooth_density_maps` is False and no smoothing pass ran.
4. NEO source hash matches the HIGH realization.
5. No GEN/CAL/TEST identity overlap for non-NEO parents.
6. Physical weights unchanged: `sum(w)` per population/bin equals `n / effective_factor`.
7. All defined posteriors finite and summing to 1 within 1e-12.
8. Closed-form vs sealed-quadrature agreement within documented tolerance.
9. Measured wall time and **measured on-disk size** for one center → storage projection for 667.

---

## 10. Validation and acceptance criteria (full build)

- exactly 667 valid centers, no missing / duplicate / malformed / mixed-configuration maps;
- every expected (population, 0.25-mag bin) cell present in the coverage table;
- sample-count and effective-k tables by center/population/bin;
- configuration identical across centers (grid, epoch, k by population, smoothing off, bin edges);
- posterior sums to 1 wherever defined; NaN preserved, never zeroed;
- every final map and config hashed → `MAP_BUILD_SEAL_V2.json`.

**TEST2 is not inspected or scored before `MAP_BUILD_SEAL_V2.json` exists.**

---

## 11. Errors to avoid

1. **`/tmp` is node-local.** Scripts staged in `/tmp` are invisible to compute nodes. Stage under
   GPFS. (Hit twice in this project.)
2. **Do not trust `n_d0_grid=400` at large k.** It is converged at k=10 only; at k>=50 it injects a
   k-dependent bias that looks like structure. §6.
3. **`_nonneo_split_fraction` raises on unknown bins by design.** Do not "fix" it by falling back to
   1.0 or to the population-level fraction — that inflates P(NEO) everywhere, uniformly, in a way no
   ROC curve reveals. Build exact 0.25-mag fractions instead. §13.3.
4. **`source_high.parquet` is one center, one bin.** It is not an all-sky NEO source. §13.2.
5. **Slurm COMPLETED does not mean the file is readable.** A corrupt HDF5 passed a file-count check
   earlier in this project. Validation must open and read every artifact.
6. **`np.array_split` is uneven.** Never hard-code shard row counts.
7. **Never widen a magnitude slice** to reach k. §4.
8. **Do not coerce NaN to 0** anywhere — density, posterior, or table.
9. Submit before cancelling when repartitioning a job; `astro` has no access to `compute`.
10. Do not write to `prob_maps_grid_s3m_nbody/` or any v1 map directory.

---

## 12. Optimizations attempted

### 12.1 Closed-form Bayesian kNN posterior mean — **RETAINED**

- **Old:** `evaluate_density_map_full_posterior_2d` — per grid point, build an `(n_d0, k)` array and
  trapezoid-integrate a normalized posterior over an 8000-node `d0` grid.
- **New:** `n0 = (k(k+1)/2 - 1/2) / (pi * sum_j d_j^2)` — one `cKDTree.query` plus a row sum.
- **Commands / resources:** `cpu-g2`, `--cpus-per-task=64`, `n_jobs=SLURM_CPUS_PER_TASK`.
- **Before / after wall time (measured, gate job 38578570, `cpu-g2`, 64 workers, k=150,
  3,000 evaluation points):** `sealed_quadrature` 6.38 s, `closed_form` 0.01 s -> **623x**.
  Whole-map effect: a (population, bin) density cell evaluates in **0.2-0.5 s** at 1001x1001
  (gate job 38578478), versus ~10 min extrapolated for the 8000-node quadrature at k=150.
- **Numerical equivalence:** max relative deviation vs the sealed quadrature at n_d0=8000 is
  8.4e-15 (k=150) / 1.8e-14 (k=250), i.e. floating-point level. The closed form is the analytic
  limit of the quadrature, so the residual is quadrature error, not fast-path error.
- **Retained:** yes — faster *and* strictly more accurate.

### 12.2 GPU kNN — **NOT ATTEMPTED (time-boxed out)**

No `cupy`/`cuml`/`torch`/`faiss` in the environment (§8.1). Writing a CUDA kNN from scratch was
outside the time box, and requesting a GPU to run CPU code is forbidden. Revisit only if a GPU
stack is installed.

### 12.3 Storage layout — see §13.1

---

## 13. Open blockers

### 13.1 Storage — **RESOLVED, fits**

Superseded by measurement. The projection in the table below assumed the production array layout;
the v2 layout stores density as float32, omits the redundant `density_raw`/`nearest_dist`, and does
not store invalid cells as all-NaN arrays. Measured in the gate:

    57.2 MB for 4 bins  ->  14.31 MB/bin  ->  630 MB/center  ->  0.420 TB for 667 centers
    free: 1.266 TB      ->  33.2% of free space, passes the >20%-headroom check

Adding the all-sky HIGH cache (§13.2, ~116 GB) and TEST2 (tens of GB) keeps the total near
0.55-0.60 TB, still under half the free space. The original worst case is retained below only to
show what the production layout would have cost.

#### Original (superseded) worst-case estimate

`df` on the compute node:

    gpfs  260T size  259T used  1.2T avail  100% /mmfs1

Free space: **1.27 TB**, on a filesystem reported 100% full.

Uncompressed projection for 44 bins x 4 populations x 1001^2:

| layout | per center | 667 centers |
|---|---|---|
| production layout (density f64 + raw f32 + support f32 + nearest f32) | 3.53 GB | **2.35 TB** |
| lean (density f32 + support i32) | 1.41 GB | **0.94 TB** |
| density f32 only | 0.71 GB | **0.47 TB** |

The production layout does not fit. Even the leanest layout is a large fraction of the remaining
free space on a shared filesystem that is already at 100%. Compressed sizes are measured in the
gate job (§9.9) before any decision.

Additional storage required by other parts of this task:
- all-sky HIGH NEO cache (§13.2): ~58 GB, plus ~58 GB for the by-pixel copy;
- TEST2 evaluation realization: tens of GB.

### 13.2 The HIGH NEO source does not exist all-sky (**ACTIVE BLOCKER**)

`outputs/more_neomod_samples_knn/source_high.parquet` contains 207,059 rows **selected for one
center (`dlon+000_lat+00`) and one bin (24 <= V < 25)**. The 6.4e8 additional draws behind it were
generated with recorded seeds (1,000,000 + 0..319), but only the survivors of that single patch and
magnitude slice were retained on disk.

To use the same HIGH realization across 667 centers and 44 bins, those seeds must be re-run and the
**all-sky** rows retained, then merged with the existing 1e8-draw cache:

    n_draws_total     = 7.4e8
    effective_factor  = 7.4e8 / 11,432,917.944 = 64.725384   (already used and validated)

This is reproducible and is genuinely the same realization — only the retention policy differs.
It is, however, a real additional build step with real storage cost (§13.1).

### 13.3 `_nonneo_split_fraction` cannot serve 0.25-mag bins

`_SPLIT_BIN_LABELS` contains only the eight 1-mag bins. Any other `(mag_min, mag_max)` raises
`ValueError` — deliberately, because guessing the retained fraction biases `P(NEO)` uniformly.

Exact 0.25-mag fractions must be computed from `nonneo_split_manifest.parquet` and the epoch-state
apparent-V magnitudes: `f(pop, bin) = N_GEN(pop, bin) / N_all(pop, bin)`. This is exact, not an
interpolation of the 1-mag values. Sparse bins are expected: TNO already has **no** `14_16` entry at
1-mag resolution, and at 0.25 mag many (population, bin) cells will be empty — those cells are
INVALID under §4, never f=1.0.

---

## 14. Final artifact manifest and hashes

Populated by `validate_neomod3_mag025_k150_maps.py` into `MAP_BUILD_SEAL_V2.json`. Empty until the
full build runs.


---

## 15. Gate results (job 38578478 build, 38578570 validate)

Center `dlon+000_lat+00`, the four 0.25-mag bins the HIGH source actually covers.
**22/22 checks PASS**, `gate_results.json` in `outputs/neomod3_mag025_k150_maps_v2_gate/`.

| gate item | result |
|---|---|
| 44 bins, 0.25 mag, half-open | PASS |
| V=23.4 -> only `[23.25, 23.50)` | PASS (`V023.25_023.50`) |
| V=25.00 not folded into top bin | PASS |
| selected V contained in slice | PASS (24.2500..24.5000) |
| NEO k=150, others k=10 | PASS |
| no silent k reduction | PASS |
| Gaussian smoothing off | PASS |
| NEO source hash = HIGH | PASS (`40490b3bc4ffaec1`) |
| NEO effective_factor 64.725384 | PASS |
| GEN/CAL/TEST disjoint | PASS (8,469,839 / 2,822,512 / 2,819,574) |
| physical weights = n / effective_factor | PASS |
| posterior sums to 1 | PASS (max 2.220e-16) |
| P in [0,1], finite | PASS |
| invalid cells absent, not zeroed | PASS |
| closed_form vs sealed_quadrature | PASS (max 7.358e-15) |
| storage fits with headroom | PASS (0.420 TB of 1.266 TB) |

Measured build: **120.0 s** for 4 bins at one center (12/16 cells valid), of which density
evaluation was ~2.5 s; the remainder is source loading and per-bin selection.

Extrapolated to 44 bins: roughly **6-8 min/center** dominated by fixed source loading, i.e.
~45-70 min wall for 667 centers at `%100` concurrency. The ~1 min/center target is not reached and
is not reachable while each task reloads the 2.6 GB epoch-state cache; the density estimator itself
is no longer the bottleneck (§12.1).

### 15.1 Trojans coverage note

All four Trojan cells at this center/bin were INVALID (`insufficient_support` / `no_samples`) — the
same emptiness the v1 `mag24+` map showed. This is recorded in the coverage table with `n`,
`k_requested`, `k_effective`, `valid`, `reason`; it is not zero-filled.
