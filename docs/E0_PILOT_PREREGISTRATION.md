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
