# E0 PILOT — PREREGISTRATION (frozen 2026-08-03, before any pilot result is viewed)

Governs the GEN map pilot that gates the full 667-map rebuild
(`EVALUATION_PROTOCOL.md` v1.1 §0.3). Centers, candidate settings and acceptance criteria are fixed
**here, in advance**. Changing any of them after seeing results requires a dated amendment saying so
explicitly.

**All E0 decisions use GEN (maps) and CAL (evaluation). TEST is not drawn, not scored, not viewed.**

---

## §1. Pilot centers (16 of 667) — fixed

Chosen to span the regimes where map construction is most likely to fail, including the two centers
with the **worst** measured ΔAUC (`new_neomod_cloning.md` §9.4) so the pilot cannot pass by only
sampling easy sky.

| # | center | regime |
|---|---|---|
| 1 | `dlon+000_lat+00` | antisun, densest NEO, fast movers |
| 2 | `dlon+020_lat-12` | busiest center overall |
| 3 | `dlon+010_lat+02` | near-ecliptic dense |
| 4 | `dlon-010_lat-02` | near-ecliptic dense, mirrored |
| 5 | `dlon+000_lat+08` | ecliptic shoulder |
| 6 | `dlon+050_lat-01` | **worst-ΔAUC cluster** (§9.4) |
| 7 | `dlon+050_lat+01` | **worst-ΔAUC cluster** (§9.4) |
| 8 | `dlon-050_lat+01` | mirror of the problem cluster |
| 9 | `dlon+090_lat+00` | mid-longitude |
| 10 | `dlon+140_lat+00` | sun-exclusion edge |
| 11 | `dlon-140_lat+00` | sun-exclusion edge, mirrored |
| 12 | `dlon+000_lat+25` | mid-latitude |
| 13 | `dlon+000_lat-25` | mid-latitude, mirrored |
| 14 | `dlon+000_lat+50` | sparse polar |
| 15 | `dlon+000_lat-50` | sparse polar, mirrored |
| 16 | `dlon-050_lat-50` | **worst single ΔAUC center** (−0.171, §9.4) |

## §2. Candidate smoothing settings — the only knobs E0 may tune

Evaluated on **CAL**. Support threshold is in **raw clones** (`scale_by_clone_factor = False`).

| knob | candidates | default |
|---|---|---|
| `support_threshold` | **2, 3, 5, 10** raw clones | 2 (behaviour-preserving) |
| `sigma_pixels` | current default only | unchanged |
| `presmoothing_passes` | current default only | unchanged |

Threshold 2 exactly reproduces the historical scaled rule (`10.0 / 8.746673 → count ≥ 2`), so it is
the reference point, not a new choice. Nothing else is tuned in E0: `support_mask_min`,
nearest-dist masking and interpolation mode are **not** E0 knobs — interpolation belongs to
E1-Resolution and is decided on CAL there.

## §3. Acceptance checks — pass/fail fixed in advance

The pilot **passes only if all of A–E pass.** Failure blocks the 667-map rebuild.

### A. Density normalisation (the §0.2 trap)
For every (population, magnitude bin) with **≥ 1,000 true objects** in the patch:

```
R = ∫ρ dA  /  (true objects in that patch+bin)
```

- **PASS: 1.05 ≤ R ≤ 1.20**, and the spread of R across bins ≤ 0.10.
- Expected ≈ **1.11**: the known kNN `k/(k−1)` bias (`new_neomod_cloning.md` §11.2, **not yet
  applied**). R ≈ 1.00 would mean the split correction is *missing* — the failure this check exists
  to catch; R ≈ 1.11/f ≈ 1.85 would mean it was applied twice.
- Bins with < 1,000 objects are reported but do not gate (Poisson-dominated).

### B. Support is statistical, not inflated
- `support_count__POP__BIN` must equal a raw integer histogram: **all values integral**, and
  `max ≤ n_clones_in_bin`.
- GEN maps must show support **≈ f_GEN × (unsplit support)** in dense bins — support *falls* with
  the split (fewer real samples) even as density is corrected upward. If support rose, the split
  fraction leaked into the support path.

### C. Non-NEO leakage is gone
- Zero ObjID overlap between the GEN objects that built the maps and the CAL rows scored against
  them (recomputed from the manifest, not assumed).

### D. Abstention and coverage
- Report per classifier: `n_scored`, `n_abstain`, and abstention rate by velocity band and
  magnitude bin. **Reported, not gated** — no preregistered threshold, because no prior measurement
  exists to justify one.

### E. CAL sanity (not a performance claim)
- CAL ROC AUC for the NEOMOD3 VDP must be **finite and > 0.5** at every pilot center with ≥ 30 NEOs
  and ≥ 30 non-NEOs.
- The chosen `support_threshold` is whichever maximises **CAL** partial AUC at FPR ≤ 0.01, ties
  broken toward the **larger** threshold (more conservative smoothing).
- **No CAL number is reported as a result.** CAL exists to choose settings.

## §4. What happens after the pilot

1. If A–E pass → write `MODEL_SEAL.json` recording every frozen setting → full 667-map GEN rebuild.
2. If any fail → fix, re-run the pilot. TEST stays sealed either way.
3. TEST is drawn and scored **once**, with all frozen variants evaluated in a **single evaluation
   event**, thresholds taken from CAL and applied unchanged.

## §5. Operational note — module imports must run on compute nodes

Importing `velocity_density_pipeline_neomod_clone_only` on a **login node** aborts:

```
F....  env.cc:93] Check failed: ret == 0 (11 vs. 0)
Thread tf_XLAPjRtCpuClient creation via pthread_create() failed.
```

JAX (pulled in transitively) tries to build its CPU thread pool and hits the login node's per-user
thread/process limit, raising **SIGABRT**. This is an environment limit, **not** a code fault — the
module compiles and imports fine on a compute node.

**Therefore: every module import, unit test and map build runs under `srun`/`sbatch`.** Only
`py_compile`, `grep` and file edits are safe on the login node. A SIGABRT with this stack trace means
"you ran it in the wrong place", not "the source is broken" — do not start reverting code.

---

# AMENDMENT 1 — 2026-08-03

**The original preregistration above is preserved verbatim** (sha256[:16] `2608c472c04f2d81`). This
amendment corrects an arithmetic error, records a verified implementation defect, and freezes
semantics that were underspecified. **Written before any E0 result was viewed:** the checks job
(38078320) was cancelled while still PENDING and produced **no output**.

The 16 threshold-2 GEN maps from job 38077241 are retained as **pilot artifacts only**, valid for
reuse if their metadata and input hashes pass A1.3 and A1.7.

## A1.1 ❌ Check A failure-mode arithmetic was WRONG — corrected

The original §3.A said "R ≈ 1.00 would mean the split correction is missing". **That is incorrect.**

Let *b* ≈ 1.11 be the kNN `k/(k−1)` estimator bias (§11.2, unapplied) and *f* = the exact
`f_GEN(pop, bin)`. The kNN estimator over the GEN sample integrates to `b · N_GEN = b · f · N_full`,
so:

| state | ∫ρ dA | **R = ∫ρ / N_full** |
|---|---|---:|
| **corrected once** (intended) | `b · N_full` | **≈ b ≈ 1.11** |
| **correction missing** | `b · f · N_full` | **≈ b·f ≈ 0.67** |
| **corrected twice** | `b · N_full / f` | **≈ b/f ≈ 1.85** |

The PASS window `1.05 ≤ R ≤ 1.20` is **unchanged and still correct** — it accepts only the
corrected-once case and rejects both failure modes. Only the diagnostic annotation was wrong.

**`N_full` is defined explicitly** as GEN + CAL + TEST objects of that population inside the *same*
sky patch (≤ 30° of the center), the *same* magnitude bin, **and** the ±5 deg/day velocity domain.
Objects outside the velocity grid contribute no density and must not appear in the denominator.

## A1.2 ⚠️ VERIFIED DEFECT — `density_raw` is stored POST-smoothing

`density_downweighted_map` is **rebound inside the smoothing loop**
(`density_downweighted_map = smooth_density_map_by_support(density_downweighted_map, …)`), and the
array stored at line 2083 as `density_maps_downweighted_raw[pop]` is that smoothed result. The name
`_raw` is misleading.

Consequences:
- Check A would integrate a **smoothed** array for any smoothed population (NEO by default).
  MBA/TNO/Trojans are not in `smooth_population_names`, so their stored arrays *are* unsmoothed —
  which is why the original Check A (non-NEO only) was not actually invalidated, but it was correct
  by accident.
- Thresholds 2/3/5/10 **cannot** be derived from the stored archives, because smoothing is already
  baked in.

**Required:** store a genuinely pre-smoothing, physically normalised array
`density_unsmoothed__POP__BIN` alongside the existing keys. Check A integrates **that** array.
It is also the common starting point for every threshold candidate.

## A1.3 Threshold candidates — separate deterministic products

Thresholds **2, 3, 5, 10** (raw clones) are produced as **separate map archives** built from the
*same* GEN rows, or derived independently from the stored `density_unsmoothed` + `support_count`
arrays. Every archive **must record its own threshold in metadata**, and E0 **asserts** the recorded
value matches the one requested. A candidate whose metadata does not match is discarded, not
reinterpreted.

## A1.4 Check D replaced — CAL-row coverage, not pixel coverage

The original Check D counted *map pixels*, which is not a coverage measurement of the evaluation.
It is replaced by a **CAL-row** test reporting, per classifier:

- `n_scored`, `n_abstain`, abstention rate
- broken down by **population**, **velocity band** (§6 strata), and **magnitude bin**

Still **reported, not gated** — but on rows, not pixels.

## A1.5 CAL threshold-selection semantics — frozen

- **Metric: macro-averaged pAUC at FPR ≤ 0.01**, averaged over pilot centers with ≥ 30 NEOs and
  ≥ 30 non-NEOs. Macro, not pooled: pooling would let the few dense ecliptic centers decide a
  setting applied to all 667.
- **All four candidates are evaluated on IDENTICAL CAL rows** — the common scorable subset across
  all candidates.
- **A candidate may not win by abstaining on hard rows.** Selection is on the common subset, and
  total coverage is reported beside every candidate. Additionally, a candidate whose coverage is
  more than **2 percentage points** below the best-coverage candidate is **disqualified**
  regardless of pAUC.
- Ties broken toward the **larger** threshold (more conservative smoothing).

## A1.6 Seal sequence

1. E0 passes → write **`MAP_BUILD_SEAL.json`** (GEN rows, split hashes, smoothing threshold, code
   commit, per-archive metadata) → run the 667-map build.
2. **`MODEL_SEAL.json` is written only after ALL remaining CAL-only choices are complete** —
   interpolation mode (E1-Resolution), probability calibration, and operating thresholds.
3. TEST is drawn and scored once, after MODEL_SEAL exists.

## A1.7 NEO provenance assertion — required, not assumed

The map build **must assert at runtime** and record in archive metadata that:

- the NEO density sample comes **exclusively** from the independent NEOMOD3 **GEN** realisation —
  record the cache path, `n_draws`, seed, and `cache_metadata.json` hash;
- the **268,511 legacy S3M NEO rows** retained in the Stage-0 cache by `VDP_LOADER=s3m`
  **never enter the NEO density tree**.

Code reading confirms this holds (`df` is replaced by the NEOMOD3 cache before the sky cut, and
there is deliberately no K|M fallback for NEO), but it is currently guaranteed only by a comment.
E0 requires a runtime assertion and stored provenance.

---

**Status:** E0 is **NOT cleared**. Relaunch only after A1.2, A1.3 and A1.7 are implemented and this
amendment is committed. TEST remains untouched.
