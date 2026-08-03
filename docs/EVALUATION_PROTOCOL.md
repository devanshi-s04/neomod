# EVALUATION PROTOCOL v1.2 — FROZEN 2026-08-03

**Status: FROZEN before reading any further results.** Amendments require a version bump (v1.1, …)
with a dated changelog entry. No metric, cut, or dataset in this document may be changed to
accommodate a result.

**Purpose.** Stop the churn. Prior sections of `new_neomod_cloning.md` reversed conclusions five
times — not because the classifiers changed, but because the *evaluation* changed underneath them
(prior, threshold placement, benchmark population, score semantics). This protocol fixes the
evaluation so that a result means the same thing every time it is quoted.

**Standing rule: no "overall winner" statement.** Report TPR/FPR and prevalence-dependent quantities
separately, always with the prior stated. A classifier must never appear to change because someone
chose a different prior.

---

## §0. SPLIT DESIGN — GEN / CAL / TEST  (v1.1, supersedes v1.0 §1)

**S3M NEOs do not appear in the primary evaluation.** NEOMOD3 is used throughout; circularity is
prevented by *independent realisations and strict dataset separation*, not by substituting an
obsolete population as a competing truth.

### §0.1 The three frozen splits

| split | role | may be used for |
|---|---|---|
| **GEN** | builds the density maps | map construction only |
| **CAL** | all model/calibration decisions | smoothing, support rules, interpolation, thresholds, calibration |
| **TEST** | locked final evaluation | **one** deterministic comparison; never a model decision |

**No object or clone lineage may cross splits.**

- **NEO** — three *independent NEOMOD3 draws* with distinct seeds. Independent draws (not a
  partition) because NEOMOD3 is a distribution; each draw is a full-sky realisation.
- **non-NEO (MBA/TNO/Trojans)** — the un-cloned Stage-0 objects **partitioned by parent ObjID**,
  disjointly, into the same three roles.

Multiple seeds measure **finite-sampling variation only**. They are not physical model uncertainty
unless NEOMOD3 parameter uncertainties are varied too, which v1.1 does not do.

### §0.2 ⚠️ NORMALISATION TRAP — partitioning changes the density scale

`density_downweighted = density_clone / effective_factor`
(`velocity_density_pipeline_neomod_clone_only.py:1991`).

Non-NEO densities are built from **real object counts**, so a GEN split holding fraction *f* of the
objects produces a density **f× too low**. Uncorrected, ρ_MBA falls while ρ_NEO does not, and
`P(NEO) = ρ_NEO/Σρ` is **inflated everywhere** — a silent, uniform bias that no ROC curve reveals.

**Correction:** `effective_factor = clone_factor × f_GEN` for each non-NEO population, so the stored
density represents the *full* population. Fold it into `effective_factor`, **not** `clone_factor` —
`clone_factor` is consumed by `np.repeat()` (lines 969–987) and must stay an integer.

NEO needs no such factor: `w_abs = total_weight / n_draws` already normalises any draw size, provided
the GEN cache's `cache_metadata.json` records **its own** `n_draws`.

**This correction is verified by E0, not assumed** (§0.3).

### §0.3 E0 — population and map fidelity (runs BEFORE any classifier comparison)

Establish that the map represents its generating population:

1. GEN vs TEST orbital element distributions (KS per element).
2. GEN vs TEST (v_λ, v_β), sky-position and magnitude distributions.
3. Genuine local support for high-velocity NEOs — real clones, not smoothing fill.
4. Effective sample counts and nearest-sample distances per region.
5. **Density normalisation**: ∫ρ dA vs true object counts per population per magnitude bin — this is
   what catches a wrong §0.2 factor. Expect ~1.0 after also applying the `(k−1)/k` correction
   (`new_neomod_cloning.md` §11.2); ~1.11 without it.
6. Explicit population priors recorded.
7. Unsupported queries → **abstention**, never automatic 0 or 1 (§3).
8. Convergence curves vs number of NEOMOD3 samples.

E0 tests whether map construction distorts NEOMOD3 *before* digest2 enters.

### §0.4 Classifiers

| classifier | status |
|---|---|
| **NEOMOD3 VDP** (GEN-built) | primary |
| **deterministic digest2** | primary comparator (§1.1) |
| legacy production VDP (±2, S3M) | **historical reference only — not a truth model** |

### §0.5 Downstream experiments

| id | what | rows |
|---|---|---|
| **E1** | controlled idealised benchmark | locked TEST |
| **E1-Resolution** | VDP ablation: (1) nearest sky+mag bin, (2) continuous mag, (3) continuous sky, (4) continuous both, (5) denser grid — interpolate densities and priors *before* forming the posterior | identical TEST rows |
| **E1-Information** | Complete-X: add Δt, per-detection times/mags/filters, astrometric+photometric errors, observer geometry, then the full two-detection tracklet; calibrated GBM density-ratio as an information diagnostic | identical TEST rows |
| **E2** | forward-simulate an independent NEOMOD3 TEST population through Sorcha (cadence, footprint, efficiency, filter photometry, noise, real tracklet building); paired clean vs observed rows for the *same* objects | new |
| **E3** | archived NEOCP/MPC candidates, original pre-confirmation tracklets, real times/obscodes/astrometry | real data |

E2 **replaces** the old S3M-NEO Sorcha test. E1-Information decides whether what the VDP discards
matters: if complete X does not improve the locked metrics, the extra complexity is unnecessary; if
it does, feature ablation says what the VDP should retain.

### §0.6 Final claim structure — kept separate, never merged

- **E1** — performance under the best-current debiased NEOMOD3 population
- **E1 ablations** — information lost to map construction and discretisation
- **E2** — robustness under realistic LSST observation
- **E3** — external behaviour on real candidates

---

## §1. Datasets

Two evaluations, reported separately and never merged.

### E1 — Geometric benchmark (populations matched to the maps)

| field | value |
|---|---|
| file | `outputs/benchmark_tracklets_neomod3/tracklets_benchmark_neomod3.parquet` |
| md5 (first 16) | `8a3d978beef124f4` |
| size / rows | 437 MB / **3,220,714** |
| NEO / MBA / TNO / Trojans | 24,998 / 3,075,543 / 37,818 / 82,355 |
| class fraction (NEO) | **0.776%** |
| epoch | 2027-08-25T00:00:00 |
| NEO source | independent NEOMOD3 draw, 30M orbits, **seed 20270825** (map cache uses 42) |
| non-NEO source | Stage-0 n-body epoch cache |
| sampling | **none** — true absolute counts, scale = 1.000000 |
| astrometry | exact (no noise) |

**Scope:** answers *"does the classifier separate these populations?"* Not an operational estimate —
no footprint, detection efficiency, trailing loss, or linking (see `new_neomod_cloning.md` §11.5).

### E2 — Operational simulation (Sorcha, detected + linked)

| field | value |
|---|---|
| file | `outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet` |
| md5 (first 16) | `06026197993b00fe` |
| rows | 648,769 (NEO 148,773 = **22.93%**) |
| ⚠️ defect | non-NEOs **capped at 500k** by pipeline stage-4 |

**E2 as it stands is NOT admissible for any prevalence-dependent metric.** Before use, either
(a) regenerate uncapped, or (b) recover exact per-row inclusion weights from the pre-stage-4 shards.
Until one of those exists, E2 supplies **TPR/FPR and AUC only**.

### Maps / classifiers scored on identical rows

| classifier | source | commit |
|---|---|---|
| NEOMOD3 VDP (±5) | `prob_maps_grid_neomod3_full/` (667 maps, 108 GB) | per §7 |
| production VDP (±2) | `prob_maps_grid_s3m_nbody/` (667 maps, 18 GB) | per §7 |
| digest2 | v0.19, **`repeatable` config** (deterministic — never `random`) | — |

### §1.1 digest2 determinism — REQUIRED SETUP

digest2 uses a **Monte Carlo** method and, by default, **seeds its PRNG randomly**, so the same
tracklet can score differently between runs. Every digest2 result in this project before
2026-08-01 was produced this way — including the §9.3 scores in `new_neomod_cloning.md` and the
discrepancies found by the night-61642 audit.

**The keyword is `repeatable`.** Source: `digest2/d2cli.c:238` parses it; `digest2/digest2.c:325`
then sets `tk->rand64 = 3`, a constant reseed per tracklet. `digest2/OPERATION.md` documents
`repeatable` / `random` as the controlling pair.

#### Where to set it

`neomod/pipeline/sorcha_phase2.py`, in the digest2 runner (~line 354):

```python
cfg_text = "noheadings\nnorms\nrepeatable\nNEO\n"
```

Keyword order does not matter; `repeatable` must simply be present.

> ⚠️ **Editing `digest2/d2_roc_comparison.config` does NOT work.** `sorcha_phase2.py` writes its
> **own temporary config** to a `NamedTemporaryFile` and passes that with `-c`, ignoring the
> checked-in `.config` file entirely. Changing the `.config` file looks like a fix and changes
> nothing. The only effective place is `cfg_text` in the source.

#### How to verify it is actually on

1. **Static** — confirm the keyword reaches the binary:
   ```bash
   grep -n "cfg_text" neomod/pipeline/sorcha_phase2.py     # must contain repeatable
   ```
2. **Empirical (authoritative)** — score the same rows twice and require identical output:
   ```bash
   for i in 1 2; do
     conda_prep/bin/python neomod/pipeline/sorcha_phase2.py run-digest2 \
       --work-dir /tmp/d2check_$i --digest2-dir digest2 \
       --d2-row-start 0 --d2-row-stop 1000 --overwrite
   done
   # scores must match EXACTLY; any difference means repeatable is not in effect
   ```

A run whose determinism has not been verified this way is **not admissible** under this protocol.
Record the verification alongside the results table (§7).

#### Why this matters beyond reproducibility

A stochastic classifier cannot be compared at a fixed threshold: part of any measured difference is
resampling noise rather than discrimination. It also makes bootstrap intervals (§5) invalid, because
they would conflate Monte Carlo noise with sampling variance.

> The commit is **not embedded here** — a commit cannot contain its own hash. Every results table
> records the commit it was produced under (§7), which is the binding reference. The working tree was
> clean at freeze; the protocol and all pipeline code were committed together.

---

## §2. Magnitude sets

Report both. Neither is "the" answer.

| set | cut | role |
|---|---|---|
| **M1 (primary)** | `mag_app < 24.5` | LSST single-visit depth |
| **M2 (sensitivity)** | `mag_app < 25.0` | the maps' own bin span |

Both use the same lower bound `mag_app >= 14.0`. A conclusion that holds in M1 but not M2 (or vice
versa) must be reported as magnitude-dependent, not stated unqualified.

---

## §3. Invalid-row policy

**An invalid or unavailable score is `NaN` — an abstention. It is NEVER silently zero.**

| situation | value | counted as |
|---|---|---|
| tracklet outside the map's velocity grid | `NaN` | abstention |
| magnitude outside all map bins | `NaN` | abstention |
| digest2 failure / timeout / hang | `NaN` | abstention |
| digest2 returns a parsed score | that score | scored |
| population density genuinely zero at an in-grid point | `0.0` | **scored** — a real prediction |

Every table reports **`n_scored` and `n_abstain` per classifier**. Metrics are computed on rows where
*all* classifiers produced a score (the common-support set), and the abstention counts are published
alongside so that coverage differences are visible rather than hidden. Silently mapping an abstention
to 0 turns "no opinion" into "confidently not a NEO" and inflates any classifier whose failures
correlate with the negative class — the exact bug found in the digest2 audit (silent-zero).

---

## §4. Metrics

### §4.1 Primary — prevalence-independent

Report these first, always:

1. **ROC AUC**
2. **Partial AUC at low FPR** — FPR ∈ [0, 0.01] and [0, 0.001], standardised (McClish). This is the
   regime a survey actually operates in.
3. **TPR at fixed FPR** — FPR = 10⁻⁴, 10⁻³, 10⁻², 10⁻¹
4. **FPR at fixed TPR** — TPR = 0.60, 0.70, 0.80, 0.85, 0.90, 0.95

**TPR and FPR are the primary reported quantities.** They are properties of the classifier and do not
move when the population mix changes.

### §4.2 Derived — prevalence-dependent, prior always stated

Contamination is *computed*, never measured directly:

```
C(π) = (1−π)·FPR / ( π·TPR + (1−π)·FPR )
```

Report `C` at each of these priors, labelled:

| π | label | basis |
|---|---|---|
| 0.00776 | **geometric synthetic prior** | E1's own composition — a property of the simulated sky, **NOT an LSST candidate-stream prior** |
| 0.05, 0.10 | illustrative | shows prior sensitivity |
| measured | **candidate-stream prior** | only when E2 is uncapped/reweighted; unavailable at v1.0 |

> ⚠️ **The 0.776% figure is a geometric synthetic prior.** It is what fraction of *objects present in
> the sky above the magnitude cut* are NEOs. The prior that matters operationally is the fraction of
> the *candidate stream reaching the classifier*, which is shaped by detection, linking and any
> upstream filtering — and is not currently measured. Do not describe 0.776% as "the LSST prior".

### §4.3 Secondary — clearly labelled prior-dependent

**Best F1**, and the threshold achieving it. Must always be printed with the prior it assumes.
**Never** used as the headline or to declare a winner: at π = 0.00776 the F1 optimum migrates to low
completeness, which is a statement about the prior, not about classifier quality.

### §4.4 Comparison rule (this caused a documented false conclusion)

**Classifiers are compared at MATCHED completeness or MATCHED FPR — never at each classifier's own
best-F1 threshold.** §9.3 compared operating points at 79.2% vs 73.6% completeness and concluded
NEOMOD3 had worse contamination; at matched completeness NEOMOD3 is *better* from 70–95%. Comparing
at self-selected thresholds measures threshold placement, not discrimination.

---

## §5. Uncertainty

**Bootstrap over sky cells, not rows.** Resample the 667 map cells with replacement (B = 1,000),
recomputing every metric per replicate; report 95% percentile intervals.

Rows within a cell share a map, a sky direction and a local population mix, so they are **not
independent**. Row-level bootstrap understates uncertainty — and sky direction is the strongest
confounder measured (ΔAUC spans +0.03 to −0.11 across directions, `new_neomod_cloning.md` §9.4).

A difference between classifiers is reported as **significant only if the paired per-cell difference
interval excludes zero.**

---

## §6. Mandatory subsets

Every primary metric is reported overall **and** for:

| subset | strata |
|---|---|
| velocity | \|v\| ≤ 0.25, 0.25–0.5, 0.5–1, 1–2, **>2 (fast movers)** |
| magnitude | the 8 map bins |
| sky direction | \|Δlon\| 0–20, 20–50, 50–90, 90–130, 130–180 |
| ecliptic latitude | \|β\| 0–2, 2–8, 8–18, 18–35, 35–50 |
| population (FP source) | MBA / TNO / Trojans / other |

Rationale: the \|v\|>2 stratum contains **zero non-NEOs** (7,727/7,727 were NEO), so it can raise TPR
without any possibility of raising FPR — pooled metrics hide this, and it must not be reported as a
completeness/contamination "trade-off". Sky direction and velocity band both flip conclusions when
pooled.

---

## §7. Versioning — required beside every table

```
protocol      : v1.0
dataset       : E1 | E2
file + md5    : <path> <md5[:16]>
rows          : n_total, n_NEO, class_fraction
magnitude set : M1 (<24.5) | M2 (<25.0)
filters       : in-grid, sun-exclusion 140 deg, mag bounds
weights       : none | per-row inclusion weights <source>
scored/abstain: per classifier
prior(s)      : listed with each C(π)
thresholds    : matched-completeness | matched-FPR (state which)
code commit   : <sha>, clean|dirty
digest2       : v0.19 repeatable
generated     : <ISO date>, job id
```

A table without this block is not quotable.

---

## §8. Deliverables at v1.0

1. **One** final comparison table on **E1** (full 3.22M rows, no subsampling) — TPR/FPR primary,
   C(π) derived at the stated priors, subsets per §6, cell-bootstrap intervals.
2. **A separate** evaluation on **E2**, uncapped or inclusion-weighted. If neither is available,
   E2 reports AUC/TPR/FPR only and says so explicitly.
3. No combined "winner" statement across E1 and E2. They answer different questions.

## §9. Open items deliberately NOT in scope at v1.0

| item | why deferred |
|---|---|
| astrometric-noise study | needs one shared perturbation applied to positions, both classifiers reading the same perturbed rows — see below |
| candidate-stream prior | requires uncapped E2 |
| resolution experiment (finer sky cells / mag bins) | tests whether the VDP's observation-space advantage is being lost to binning |
| `(k−1)/k` density correction | `new_neomod_cloning.md` §11.2; changes no score, batch with next rebuild |

**Noise study design note (for v1.1):** perturb the *observations* once — RA/Dec per detection — and
let both classifiers derive their inputs from the same perturbed rows. digest2 consumes positions
directly while the VDP consumes velocities derived from them, so noise propagates differently;
perturbing each classifier's inputs separately would confound the comparison.

**Resolution experiment note:** the VDP compares in observation space and so should lose less
information than digest2's observation → 6 orbital elements → 4-D model path. But the VDP as built
discards resolution digest2 keeps: sky position binned to 667 cells (~10°×5°), magnitude to 8 bins,
with only (vlam, vbeta) continuous. That predicts finer binning narrows the gap — a falsifiable test
of where the loss actually is.

---

## §10. Operational note — module imports must run on COMPUTE nodes

Importing the VDP module on a **login node** aborts with SIGABRT:
`Check failed: ret == 0 (11 vs. 0) Thread tf_XLAPjRtCpuClient creation via pthread_create() failed`.
JAX (pulled in transitively) cannot build its CPU thread pool under the login node's per-user thread
limit. **This is an environment limit, not a code fault.** Every module import, unit test and map
build runs under `srun`/`sbatch`; only `py_compile`, `grep` and file edits are safe on the login
node. A SIGABRT with that stack trace means the job ran in the wrong place — do not revert source.

Related: **`/tmp` is node-local.** A script written to `/tmp` on the login node is invisible to
compute nodes and the job dies with `No such file or directory`. Job scripts live on GPFS.

---

## Changelog

| version | date | change |
|---|---|---|
| v1.0 | 2026-08-01 | initial freeze |
| v1.2 | 2026-08-03 | Support threshold expressed as integer 2.0 raw clones (exactly equivalent; counts are integers). TEST seal split into `TEST_DATA_SEAL.json` (frozen now) and `MODEL_SEAL.json` (written after CAL). E0 pilot preregistered in `E0_PILOT_PREREGISTRATION.md`. §10 added: login-node JAX SIGABRT and node-local `/tmp`. |
| v1.1 | 2026-08-03 | **§0 added: GEN/CAL/TEST split design, normalisation trap, E0 fidelity, E1-Resolution/E1-Information/E2/E3.** S3M NEOs excluded from the primary evaluation; legacy VDP demoted to historical reference. E1a (built 2026-08-01) is superseded: its non-NEOs were 100% shared with the maps. |
| v1.0.1 | 2026-08-03 | added §1.1 (digest2 determinism procedure + the `.config` trap). **No specification change** — §1 already required `repeatable`; this documents how to achieve and verify it. |
