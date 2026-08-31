# Gorica, 24 Aug 2026 — Legacy vs New VDP: Advisor Requests and Living Plan

**Status:** Active living document  
**Created:** 2026-08-24  
**Owner:** Devanshi  
**Scope:** Current VDP map/scorer state, comparison with the legacy classifier and digest2, map-completeness work, and the frozen Rubin nightly scoring interface.

## How to use this document

This is the single working record for the discussion with the advisor and the coordination with Hyak. Update it whenever a result, decision, job state, file path, map seal, or priority changes.

Do not silently rewrite an old decision. Mark it **superseded**, add the replacement and date, and append a short entry to the changelog. Keep measured facts, current decisions, hypotheses, and open questions visibly separate.

## One-paragraph recap

The current NEOMOD3 VDP maps and TEST2 scoring run are complete on Hyak. Globally, the new classifier and the legacy classifier have nearly tied ROC AUC, but the legacy classifier is cleaner in the operationally interesting top-ranked range, and the new classifier substantially lowers the median NEO score while barely changing the non-NEO median. The likely explanation is excessive effective smoothing caused by increasing NEO kNN from 10 to 150 while also narrowing the magnitude bins from 1.0 to 0.25 mag, but that is still a hypothesis that must be demonstrated in a dedicated notebook. Separately, only about 44% of TEST2 rows have all four population densities, mainly because the sparse TNO and Trojan maps are incomplete. Returning `NaN` is scientifically correct when the probability denominator is incomplete, but it is not an acceptable operational endpoint. Our immediate priority is to freeze a stable, versioned Rubin tracklet-to-score interface around the exact current Hyak scorer and maps. After that, we will make the legacy-versus-new difference visible, test population-specific physical cloning as the first coverage remedy, and run selected-center mechanism ablations using existing data only.

## What Devanshi and the advisor want

1. A frozen, explicit path from a Rubin-like tracklet text/CSV file to an output CSV containing follow-up probability scores.
2. Plots that first **show exactly where** the legacy and new VDP classifiers differ. We should not label the discrepancy only as a “low-FPR problem” before the notebook proves the score/ROC regime in which it occurs.
3. Direction-specific ROC comparisons for small regions of the sky, reusing the older `mag245` n-body implementation that already made 5° × 5° single-night sky-direction ROC tiles.
4. A contained mechanism notebook that determines why the old classifier did better where it did, using controlled ablations rather than guesses.
5. A careful explanation of digest2’s apparent quantization. digest2 must remain an unmodified black box: tracklets go in and its official scores come out.
6. Operationally complete maps. A valid, in-domain Rubin tracklet should normally receive a score rather than becoming `NaN` because one simulated population was undersampled.
7. Use the existing/older evaluation data for diagnosis. Do **not** draw or generate a new final evaluation sample until the map design and current tests are finished.

## Current decisions

These are agreed unless explicitly superseded later.

- Freeze the exact latest Hyak VDP implementation, map set, and scoring semantics as Rubin scorer v1 before changing map methodology.
- “Frozen v1” means a stable engineering interface and regression target. It does **not** mean the current maps are operationally ready; current complete four-density coverage is only about 44%.
- The scorer must group tracklets by the map/sky center they fall in, load each map once, and then score the relevant magnitude-bin groups. This is the intended nightly execution pattern.
- A two- or three-population denominator must never be reported as `P(NEO)`. Missing population densities remain explicit invalid/abstention reasons in the frozen diagnostic behavior.
- The future operational design should reserve `NaN` mainly for malformed or genuinely out-of-domain inputs, not finite-simulation holes inside the declared domain.
- The first map-completeness remedy to test is clean, population-specific cloning with explicit normalization weights.
- Physical cloning is performed on the global source population **before** selecting map centers and magnitude slices: clones occupy hypothetical orbits, are propagated/projected to the target epoch, have their geometry and apparent magnitude computed, and only then enter the applicable sky and magnitude cells.
- MBA, TNO, Trojan, and—if needed—NEO may have different clone factors. Their physical weights must preserve the intended population abundance.
- Mechanism diagnosis may reuse already inspected TEST/CAL/older rows because it is a paired causal comparison, not a new out-of-sample performance claim. This limitation must be written on every such result.
- No new evaluation population will be generated yet. A fresh final evaluation is a later step, after the map policy is frozen.
- digest2 will not be modified. We will audit its inputs, official returned values, and our plotting/parsing only.
- Full-grid ablation builds are not currently appropriate. Mechanism studies will use selected centers spanning sky and coverage regimes.
- The frozen v1 input is honestly V-only. This is temporary: the eventual Rubin interface must accept the actual LSST filter and photometric uncertainty for each detection and either use validated per-band maps or an explicit, uncertainty-aware band-to-V model. It must never relabel an LSST-band magnitude as V.
- A missing TNO/Trojan density estimate is not evidence that the physical density is zero. The next scientific task is to make the dominant NEO-versus-MBA decision usable while repairing sparse-population coverage.

## Current map and scorer configuration to freeze

The latest completed Hyak configuration is currently understood to be:

- 667 sky centers.
- 44 apparent-V magnitude bins, 0.25 mag wide, covering `14 <= V < 25`.
- NEOMOD3 HIGH source for NEOs.
- `k_NEO = 150`.
- `k_nonNEO = 10` for MBA, TNO, and Trojan maps.
- Velocity grid from -5 to +5 deg/day with 0.01 deg/day spacing.
- No extra Gaussian smoothing in the new maps.
- No support masking or new calibration layer.
- Local four-density normalization for NEO, MBA, TNO, and Trojan.

Before the v1 seal is written, Hyak must confirm the exact commit, scorer source hash, map manifest/seal, map roots, configuration, dependency environment, and any runtime flags. The completed Hyak work was reported at commit:

`87f6bd82d190e946f99ab128ff8ffff380d09a7a`

That commit was reported as not yet pushed/transferred, so the Mac checkout must not be assumed to contain the authoritative implementation.

## What is already measured

### Rubin scorer v1 feasibility findings

Hyak confirmed that commit `87f6bd82d190e946f99ab128ff8ffff380d09a7a` reconstructs the TEST2 scorer: `pipeline/score_test2.py` at that commit is byte-identical to the scoring working tree. The subsequently dirty validator and runbook are not imported by the scorer.

Two related coordinate facts must not be conflated:

- Ecliptic sky coordinates `lambda, beta` and rates `v_lambda, v_beta` **can be derived from an observed two-detection RA/Dec tracklet**. This is the correct operational input path.
- TEST2’s stored `lambda, beta` were cached model-frame values from propagation, rather than values transformed back from the apparent observed RA/Dec. A real tracklet does not contain those particular cached model values.

From the two TEST2 detections, equatorial rates were recovered to at most about `2.6e-08 deg/day`, and ecliptic rates to at most about `3.4e-08 deg/day`. Thus the velocity coordinates are operationally recoverable. Transforming apparent RA/Dec back to ecliptic position differed from the stored model-frame position by about 0.0018 deg (6.5 arcsec) at the median. That changed the nearest-center assignment for 295 of 688,688 rows (0.043%), all near center boundaries.

The correct interpretation is not “lambda/beta are unobservable.” It is: **the old cached model-frame definition is unavailable to Rubin, so v1 must define and seal an observable RA/Dec-to-ecliptic transformation.** Regression should retain two paths: a precomputed-geometry path to prove exact reproduction of the old TEST2 scorer, and the operational derived-geometry path to quantify the deliberate assignment change.

Hyak also measured an MPC-80 round trip on 5,000 TEST2 tracklets. This is relevant only to the **separate digest2 comparison path**: the Rubin VDP scorer consumes the native full-precision tracklet CSV and must never round-trip its inputs through MPC-80. The MPC representation retained less precision than the canonical source values:

- RA maximum difference: about `2.1e-05 deg`.
- Dec maximum difference: about `2.7e-05 deg`.
- Magnitude maximum difference: about `0.05 mag` because these lines recorded magnitude to 0.1 mag.
- Median derived-rate changes: about 0.059 and 0.042 of a 0.01-deg/day velocity-grid cell; the maximum `v_lambda` change was about 1.29 grid cells.
- Magnitude-bin assignment matched the full-precision source for 4,372/5,000 rows (87.44%); 628 rows (12.56%) changed bins.

This does **not** mean MPC-80 has more significant figures than required by a 0.25-mag bin, nor that a 0.25-mag bin is too wide. It means the reverse: serializing the full-precision photometry to the tested MPC-80 representation rounded each magnitude to 0.1 mag. The mean of two rounded magnitudes therefore moves in 0.05-mag increments, and values near a hard 0.25-mag boundary can cross that boundary. This loss does **not** affect VDP scoring because VDP uses the original Rubin CSV. MPC-80 is generated only as digest2’s required black-box input, and the resulting digest2 score is merged back to the original tracklet ID. Its representation loss is therefore a limitation/provenance note for digest2 comparison, not a reason to alter the VDP bins or Rubin interface.

Two implementation details remain to be sealed before the observable geometry definition is final:

1. Which reference position represents a tracklet for center assignment—first detection, midpoint direction, or position propagated to the map reference time?
2. Whether the coordinate transform and antisolar longitude use the tracklet observation time, a map epoch read from the seal, or both under an explicitly declared validity window. A hard-coded 2027 epoch must not be presented as a general Rubin-nightly rule.

### TEST2 composition

TEST2 contains 688,688 rows:

| Population | Rows | Note |
|---|---:|---|
| NEO | 49,887 | Fresh independent NEOMOD3 NEO rows; physical weight about 0.0999994694 |
| MBA | 614,710 | Reused frozen prior TEST non-NEO rows |
| TNO | 7,575 | Reused frozen prior TEST non-NEO rows |
| Trojan | 16,516 | Reused frozen prior TEST non-NEO rows |

All three scoring modes completed. The common-scorable comparison set contains 302,666 rows (43.95%), including 14,026 NEOs.

### Global comparison on common-scorable rows

| Classifier | ROC AUC | historical normalized pAUC (raw/L, FPR<=0.01) | Completeness at 5% contamination | Completeness at 10% contamination |
|---|---:|---:|---:|---:|
| New VDP | 0.93814 | 0.71476 | 0.4252 | 0.5391 |
| Legacy VDP, raw | 0.93695 | 0.72007 | 0.5506 | 0.5987 |
| digest2 | 0.92474 | 0.63491 | 0.5369 | 0.5369 |

The global AUCs are nearly tied. That summary alone hides the difference that matters for ranking follow-up candidates.

### Current evidence for the affected score regime

| NEO completeness | New precision | Legacy precision |
|---:|---:|---:|
| 20% | 0.9929 | 0.9964 |
| 30% | 0.9836 | 0.9976 |
| 40% | 0.9606 | 0.9894 |
| 50% | 0.9261 | 0.9723 |

Legacy is cleaner throughout approximately 20–50% completeness, not merely at one arbitrarily named low-FPR threshold. The median NEO score is 0.2804 for the new maps versus 0.5078 for legacy, while the non-NEO medians are almost unchanged at 0.0005 and 0.0006.

These measurements motivate the phrase **“deterioration in the top-ranked/operational region.”** The evidence notebook still needs to make this visible and quantify uncertainty before we claim a mechanism.

### Four-density availability

- NEO-cell availability: 99.33%.
- Complete four-density coverage: 43.96% unweighted; approximately 45.06% weighted.
- Four-density coverage by true population:
  - NEO: 28.12%.
  - MBA: 43.71%.
  - TNO: 51.88%.
  - Trojan: 97.46%.

The dominant failure is not unstable NEO maps. It is missing rare-population density estimates—especially TNO and Trojan—in narrow 0.25-mag slices. A row with an incomplete denominator is correctly marked invalid instead of being assigned an inflated pseudo-probability, but this diagnostic policy exposes that the maps are not operationally complete.

### Strict four-population probability versus an operational NEO–MBA score

Three states must be kept distinct:

1. **Valid, physically very small density:** the population model supports a finite estimate near zero. This density can participate honestly in the four-population denominator.
2. **Valid physical zero:** the model establishes that the population has zero support in the declared cell/domain. Zero may be used only when this is a model conclusion, not because no finite simulation point happened to land there.
3. **Unavailable estimate:** the current finite map cannot estimate the population density. This is what `missing_population_density:TNO,Trojans` means today. It is not evidence for either of the first two states.

The frozen v1 behavior remains scientifically correct for the quantity it reports: `p_NEO`, `p_MBA`, `p_TNO`, and `p_Trojan` are a normalized four-class posterior and are emitted only when all four density estimates are valid. Omitting unavailable TNO/Trojan terms and still calling the result `P(NEO)` would overstate what is known.

However, that strict posterior is not sufficient as the only operational ranking output. In large regions the main decision is NEO versus MBA, and returning only `NaN` plus a missing-population reason discards useful information. The next scorer/map design should therefore expose two explicitly different quantities:

```text
p_NEO_4pop     = rho_NEO / (rho_NEO + rho_MBA + rho_TNO + rho_Trojan)
                 valid only when all four estimates are valid

score_NEO_MBA  = rho_NEO / (rho_NEO + rho_MBA)
                 valid when NEO and MBA estimates are valid
                 explicitly a conditional/binary ranking score, not full P(NEO)
```

There is a third—and potentially better operational—formulation if `OTHER` is constructed as one physically weighted aggregate of every non-NEO population:

```text
p_NEO_vs_OTHER = rho_NEO / (rho_NEO + rho_OTHER)

rho_OTHER = density of the complete non-NEO mixture
          = MBA + TNO + Trojan contributions with correct physical weights
```

If `rho_OTHER` is merely calculated by summing the three current component maps, it is algebraically identical to `p_NEO_4pop` and retains the same missing-map failure. The useful alternative is to estimate the non-NEO mixture directly from pooled, physically weighted non-NEO samples. That matches the actual follow-up question—NEO versus everything else—and does not require separate MBA/TNO/Trojan probabilities on every row.

Direct pooling must still be validated carefully. MBA abundance must not wash out compact Trojan/TNO modes, the populations must be disjoint, and cloning/weights must preserve their physical mixture. Where all four component maps are valid, `p_NEO_vs_OTHER` should agree with `p_NEO_4pop`; that is the primary cross-check. The production target may therefore become `p_NEO_vs_OTHER`, while component densities remain diagnostic outputs where available.

Rows using the binary score must also carry `missing_populations`, `fourpop_valid`, and a fallback/coverage state. This makes the row useful for follow-up ranking without pretending that TNO/Trojan contamination was modeled. The engineering reason code remains useful for debugging, but it cannot be the only operational output.

The final objective is still to make `p_NEO_4pop` available throughout the intended domain by repairing TNO/Trojan density estimation. The binary score is an honest bridge and a useful diagnostic, not a substitute for that map work.

### Rubin photometry roadmap

`mag0_V` and `mag1_V` are correct names for the current frozen maps because those maps were built in apparent V. They are not the final Rubin input contract. Rubin observations arrive in LSST passbands, so a later interface must include, per detection:

```text
mag0, band0, magerr0
mag1, band1, magerr1
```

Two future approaches should be compared rather than assuming a universal color offset:

- build/seal separate density maps in each LSST band using the population color/SED model; or
- convert/marginalize from the observed LSST band to V using a validated population-aware color model and propagate its uncertainty across magnitude bins.

Direct per-band maps are likely cleaner because asteroid colors vary and a fixed band-to-V conversion can bias populations differently. In either approach, hard assignment from an uncertain converted magnitude should be tested against interpolation or marginalization across neighboring magnitude slices. This is deferred while the current map methodology is being diagnosed, but it is a required operational milestone before claiming Rubin readiness.

## Legacy smoothing: now verified

This is no longer an open metadata question. The legacy maps used:

```text
smooth_density_maps      = True
smooth_population_names  = ['NEO']
smooth_support_threshold = 10.0
smooth_sigma_pixels      = 4.0
```

The pixel scale makes this a Gaussian sigma of 0.04 deg/day. The NEO array read by the legacy scorer is the post-smoothed array. MBA, TNO, and Trojan arrays were not Gaussian-smoothed.

## Working bandwidth hypothesis—not yet a conclusion

The previously cited neighbor-radius change of about 0.031 to 0.045 deg/day came from a one-magnitude, one-center k-sweep and is not the production configuration. Production uses 0.25-mag slices, with roughly one quarter as many samples per slice. Since a two-dimensional kNN radius scales approximately as `sqrt(k/n)`, using `k = 150` rather than 10 while narrowing the bins may give a production NEO radius around 0.09 deg/day.

For rough context only, the legacy NEO effective scale can be approximated as:

`sqrt(0.031^2 + 0.04^2) ~= 0.05 deg/day`

This makes excessive smoothing/peak flattening a plausible leading hypothesis: the new NEO score distribution moves downward while the non-NEO distribution barely moves. It is not yet proof. We need actual per-cell neighbor-distance measurements and ablations that vary `k` and magnitude-bin width separately; otherwise those two changes remain confounded.

## Cloning policy to test

Devanshi’s intended cloning interpretation is the default plan:

1. Begin with the global physical source population, not a preselected map cell.
2. Create population-specific clones by varying physically justified orbital phase/elements under a documented model.
3. Propagate/project every original and clone to the map epoch.
4. Compute sky position, apparent velocity, and apparent V magnitude at that epoch.
5. Assign the resulting realization to sky centers and magnitude slices.
6. Carry explicit physical weights through density estimation and final normalization.

If population `c` is cloned uniformly by a factor `f_c`, the simplest weight is `w_clone = w_parent / f_c`. The implementation should ideally carry per-row weights through the estimator rather than apply an ad hoc post-map correction. Different factors may be used for MBA, TNO, Trojan, and NEO.

This procedure can fill cells that were empty only because a finite realization placed no object there: clones from elsewhere in the global population can project into that cell at the target epoch. Cloning **after** selecting a cell cannot solve an empty cell.

Important caution: a cell with no relevant parent support under the chosen generative/cloning model cannot be repaired merely by increasing a clone factor. Hyak reported 10,005 TNO cells with `no_split_fraction`, meaning no GEN objects in that slice under the current construction. The pilot must determine how many of these are finite-realization holes that global pre-selection cloning fills. Any remaining genuinely unsupported region will need a principled pooled-neighbor, adaptive-bin, or generative treatment—not zero filling.

### Cloning pilot acceptance checks

- Complete four-density coverage, overall and by true population.
- Counts and weighted mass per population, center, and magnitude bin.
- Density integrals before and after cloning.
- Preservation of the intended physical population ratios.
- Duplicate/parent independence audit.
- Effective neighbor radii and density smoothness.
- Change in `no_split_fraction` and the reason distribution for any remaining holes.
- Paired score changes on the same existing diagnostic tracklets.

## Ordered work plan

### Corrected execution order as of 2026-08-25

1. Rubin scorer v1 freeze and interface separation — complete.
2. Legacy-versus-new proof notebook: prove or reject whether the observed advantage occupies the predicted low-FPR/top-ranked regime.
3. Coverage/scorer formulation study: compare four-component, NEO–MBA, and direct physically weighted NEO–OTHER probabilities/scores.
4. Selected-center population-specific cloning/weighting pilot.
5. Contained mechanism ablations for kNN k, magnitude-bin width, interpolation, and legacy Gaussian smoothing.
6. Select/preregister map v2, build only the chosen full grid, then obtain a fresh final evaluation if needed.

The proof notebook answers **where** legacy is better. Coverage work determines how to avoid unusable scores. The mechanism notebook answers **why** legacy is better. These questions must not be collapsed into one notebook or one causal claim.

### Priority 1 — Freeze Rubin scorer v1 — complete

First inventory the exact latest Hyak implementation and seal it without altering its scientific semantics.

Proposed artifacts:

```text
neomod/src/rubin_vdp_scorer_v1.py
neomod/pipeline/score_rubin_tracklets.py
neomod/docs/RUBIN_TRACKLET_SCORER_V1.md
neomod/seals/RUBIN_TRACKLET_SCORER_V1_SEAL.json
neomod/examples/rubin_tracklets_v1.csv
neomod/examples/rubin_tracklet_scores_v1.csv
```

Proposed wide CSV input contract:

```text
tracklet_id
mjd0, ra0, dec0, mag0_V
mjd1, ra1, dec1, mag1_V
observatory_code
```

The Rubin VDP CLI accepts the native canonical tracklet CSV only. It must not accept MPC-80 or convert Rubin rows through MPC-80. A separate digest2 comparison adapter may export the same identified tracklets to MPC-80, run digest2 as an unmodified black box, and merge digest2 results back by tracklet ID.

Output should preserve input order and identity and include:

```text
tracklet_id
p_NEO, p_MBA, p_TNO, p_Trojan
rho_NEO, rho_MBA, rho_TNO, rho_Trojan
valid, reason
v_lambda, v_beta, mean_V
map_center_id, magnitude_bin_id, map_id
scorer_version, scorer_seal, map_seal, input_hash
```

Execution requirements:

- Require an explicit map/scorer seal; never silently use “whatever is latest.”
- Group rows by map/sky center so each large map is loaded once, then group by magnitude bin and score in bounded chunks.
- Preserve row order, IDs, and deterministic results.
- Do not silently clip, drop, zero-fill, or renormalize incomplete rows.
- Emit machine-readable invalid reasons.
- Keep digest2 outside v1 as a separate black-box comparator.
- Regression-test against the frozen TEST2 output—ideally all 688,688 rows—for exact center/bin/reason behavior and score equivalence within a declared numerical tolerance.

**v1 acceptance gate:** A clean example input runs end to end, produces the documented schema, reproduces the current Hyak scorer, and records enough provenance to rerun the result. The documentation must display the 43.96% complete-density limitation prominently.

#### Completion record from Hyak

The expensive Task 1 implementation and validation run is complete. It must **not** be repeated merely to correct the MPC-80 interface boundary.

**Committed Hyak state:**

- Commit: `988e82cea3a3200de630728fe214eee43a340ab9` (`988e82c`).
- Subject: `Rubin Tracklet Scorer v1: frozen CLI around the exact TEST2 classifier`.
- 18 files changed, 1,918 insertions and 24 deletions.
- Working tree reported clean after the commit.
- Commit not yet pushed/transferred.
- Baseline scorer commit `87f6bd82d190e946f99ab128ff8ffff380d09a7a` was verified sufficient and byte-identical for the relevant scorer source.

**Seals used for the completed validation:**

```text
scorer seal  0e61ec37e263e9a4f9e81ed4ef25ef8fa069f630b8d2003f2ce29ce37de9eaeb
map seal     a010bb1bf2e421a3258ff007f7bc4034369918a3877418620cd14e3f0b9d19ab
```

**Full TEST2 regression:** 688,688 native canonical CSV rows, final jobs `38808527` and `38808528`, each completed in about 1 hour 50 minutes.

- Precomputed/model geometry isolation path reproduced the frozen oracle exactly:
  - center, magnitude bin, validity, and reason: 688,688/688,688 identical;
  - `P_NEO` maximum absolute difference: `3.553e-15`;
  - NaN pattern and reason-count tables identical;
  - coverage: 43.9582% in both;
  - within sealed tolerance: true.
- Observable RA/Dec-derived operational path:
  - magnitude bin: 688,688/688,688 identical;
  - center: 688,393/688,688 identical, with 295 adjacent-center reassignments;
  - validity differed for 4 rows and reason for 5 rows as a consequence of reassignment;
  - maximum absolute `P_NEO` change: `3.816e-03`;
  - coverage: 43.9579% rather than 43.9582%;
  - correctly reported as outside the old oracle tolerance rather than passed as exact.

The 295-row difference is the already documented cached-model-versus-observed-coordinate distinction. It is not a failure to derive ecliptic coordinates from RA/Dec.

**Native CSV coverage/reasons in the operational path:**

```text
ok                                             302,733  (43.9579%)
missing_population_density:Trojans             331,207
missing_population_density:TNO,Trojans          36,373
missing_population_density:TNO                  13,255
missing_population_density:NEO,TNO,Trojans       4,423
outside_velocity_grid                               507
plus three smaller missing-density combinations
```

**Determinism/invariance:** A concentrated 11,545-row test genuinely split center groups at chunk sizes 200,000, 997, and 500. Permuted input and changed chunk sizes produced maximum numeric difference `0.000e+00`, identical NaN/categorical results, and preserved each input order. When chunks split a center, the current delegated scorer reopens its map; outputs are unchanged but I/O rises. The default full TEST2 run opened each of 667 center maps once.

**Performance:**

- 6,575.7 seconds (1.83 hours) for 688,688 rows.
- 104.7 rows/second.
- 7.77 GB peak RSS.
- 667 map opens and 265.8 GB read.
- 99.52% of wall time was map I/O/decompression plus density scoring.
- Cost depends strongly on centers/cells touched, not merely row count.

**Engineering defects found and handled without changing the frozen map science:**

- Out-of-range V rows were dropped by `pandas.groupby(None)` and incorrectly retained reason `ok`; the wrapper now emits `v_out_of_range`.
- A malformed coordinate could reach `SkyCoord` before row validation and abort the run; validation now prevents this.
- `input_hash` was incompatible with the observed pandas 3.0 behavior; repaired.
- Map-open accounting undercounted opens when chunks split a center; reporting now distinguishes centers and actual opens.

**Interface cleanup completed:** Follow-up commit `6277cff` preserved `988e82c`, removed MPC-80 ingest from the VDP CLI, moved MPC export to `pipeline/export_mpc80_for_digest2.py`, added an auditable tracklet-ID mapping, and corrected the geometry/MPC wording. The immutable validation seal was not rewritten; lightweight interface metadata references it. All scientific sources were verified byte-identical, and a 32-row smoke test reproduced every scientific output column exactly. No full regression was rerun. Task 1 is closed; both commits were reported as not pushed.

#### Detailed Hyak task: freeze the current Rubin scorer interface

##### Explanation: what we want and why

We want to preserve the exact classifier that produced the completed current-map TEST2 results behind one stable, documented command-line interface. The interface should take a file of two-detection tracklets and return a row-aligned CSV containing the current VDP scores, component densities, validity state, and enough provenance to determine exactly which scorer and maps produced every result.

This gives us a durable baseline while the scientific map design changes. Later we may replace the maps, change cloning, change kNN bandwidth, or improve coverage, but a Rubin-facing caller should not need to be redesigned each time. The input/output contract and grouped execution path should remain stable while the sealed scientific assets underneath it are explicitly versioned.

This task is **freezing and packaging**, not classifier development. Hyak should not smooth maps, fill missing populations, recalibrate probabilities, change bin assignment, rebuild densities, generate new test objects, modify digest2, or start Sorcha. The first version must preserve the current scorer’s scientific behavior, including its `NaN`/invalid outcomes and reason codes. Its approximately 43.96% complete four-density coverage is a known limitation that must be reported, not hidden.

The result we want is:

```text
native full-precision two-detection Rubin CSV
                         |
                         v
              frozen Rubin scorer v1
          (explicit scorer + map seal)
                         |
                         v
row-aligned CSV with densities, probabilities, validity,
map assignment, derived motion, and full provenance
```

The interface is successful when it can reproduce the frozen TEST2 scores from ordinary input files without relying on notebook state, implicit paths, or “latest” files.

##### Step 1 — Identify the authoritative frozen baseline

Before writing the wrapper, report:

1. The exact Git commit and whether the working tree that produced TEST2 contained uncommitted changes.
2. The exact Python entry point/functions used for current-map scoring.
3. The complete map root, manifest, filenames, and existing seal or checksums.
4. All configuration values affecting sky-center assignment, magnitude-bin assignment, velocity calculation, density lookup, four-population normalization, and invalid-reason behavior.
5. The Python/dependency environment and any environment variables or scheduler/runtime flags needed for scoring.
6. The paths to the frozen TEST2 input and current-map output used as the regression oracle.

If commit `87f6bd82d190e946f99ab128ff8ffff380d09a7a` is not a complete reconstruction of the TEST2 scorer, stop and record the additional working-tree diff or artifacts before packaging it. Do not silently recreate missing behavior from memory.

##### Step 2 — Write the v1 contract before implementation

Create `neomod/docs/RUBIN_TRACKLET_SCORER_V1.md` and define:

- accepted input formats and required units;
- required and optional columns;
- how exactly two detections are paired into one tracklet;
- coordinate, time, magnitude, and observatory-code conventions;
- output columns, units, and meanings;
- allowed `valid`/`reason` states;
- behavior for malformed, incomplete, outside-magnitude, outside-velocity, missing-map, and incomplete-four-density rows;
- ordering and determinism guarantees;
- exit codes and file-level failure behavior;
- the distinction between the stable interface version and the replaceable sealed map/scorer version.

The canonical wide CSV should contain at least:

```text
tracklet_id
mjd0, ra0, dec0, mag0_V
mjd1, ra1, dec1, mag1_V
observatory_code
```

The documentation must state that `mag*_V` is apparent V. Do not silently interpret an LSST-band magnitude as V. Any future Rubin-band-to-V conversion belongs in a separately specified preprocessing stage unless the current frozen scorer already has an authoritative conversion that can be sealed and tested.

MPC 80-column parsing/export does not belong in this VDP CLI. If an MPC parser/exporter has already been implemented, move it to a clearly separate digest2 comparison tool. The data flow must be explicit:

```text
native Rubin CSV ──> frozen VDP scorer v1 ──> VDP scores
       |
       └──────────> MPC-80 export ──> unmodified digest2 ──> digest2 scores
                                                           |
                      merge by original tracklet_id <───────┘
```

The VDP branch must always retain the original CSV precision. The digest2 branch should record the MPC-80 export hash and the original tracklet-ID mapping.

##### Step 3 — Create an immutable seal

Create `neomod/seals/RUBIN_TRACKLET_SCORER_V1_SEAL.json`. It should include at least:

- interface version;
- Git commit and dirty-tree patch hash, if applicable;
- hashes of every scoring source file;
- map-set ID, map manifest path, and hashes of the manifest/maps or a documented Merkle/content seal;
- all scientific configuration values;
- supported input schema version;
- dependency/environment lock or hash;
- creation time and creator/job provenance;
- TEST2 oracle input/output paths and hashes;
- declared numerical comparison tolerance, if bit-for-bit equality is not portable.

The CLI must require an explicit seal or use one documented v1 default that resolves to the same immutable seal. It must never search for the newest maps automatically.

##### Step 4 — Separate the reusable scorer from file orchestration

Implement two layers:

1. `neomod/src/rubin_vdp_scorer_v1.py` — reusable functions for validation, derived tracklet quantities, map/bin assignment, density lookup, probability calculation, and reason reporting, while calling the authoritative current scorer code wherever possible.
2. `neomod/pipeline/score_rubin_tracklets.py` — file parsing, chunking, grouping, map loading, output ordering, logging, and CLI behavior.

Do not duplicate or algebraically rewrite the current scientific scorer merely to make the wrapper look cleaner. Prefer a thin adapter around the exact tested implementation. If refactoring is required, first lock the oracle and prove the refactor has not changed outputs.

##### Step 5 — Implement grouped-by-map execution

The execution path should:

1. Parse and validate all input rows without changing their order.
2. Compute or retrieve each tracklet’s assigned sky center and magnitude bin using the frozen rules.
3. Mark rows that cannot be assigned with explicit reasons.
4. Group valid candidates by map/sky-center ID.
5. Load a given large map once, score all of its rows—subgrouping by magnitude bin where useful—and release it before moving to the next map.
6. Process large groups in bounded chunks without changing numerical results.
7. Restore the exact original row order before writing the output.

Add concise run statistics: input rows, valid rows, reason counts, maps loaded, maximum group/chunk size, and elapsed time. The scientific result must not depend on input ordering or chunk size.

##### Step 6 — Return a complete audit record per tracklet

The output should contain at least:

```text
tracklet_id
p_NEO, p_MBA, p_TNO, p_Trojan
rho_NEO, rho_MBA, rho_TNO, rho_Trojan
valid, reason
v_lambda, v_beta, mean_V
map_center_id, magnitude_bin_id, map_id
interface_version, scorer_seal, map_seal, input_hash
```

Use the exact current scorer semantics for probability and invalid fields. Do not synthesize a four-class probability when one population density is unavailable. Do not replace a missing density with zero or epsilon. Preserve density information that the current scorer can honestly return and attach the exact invalid reason.

##### Step 7 — Add examples and failure fixtures

Create small committed fixtures covering:

- ordinary valid tracklets;
- more than one tracklet assigned to the same map, demonstrating grouped loading;
- tracklets assigned to different maps;
- magnitude and velocity boundaries;
- incomplete four-density lookup;
- malformed coordinates/times/magnitudes;
- incomplete or duplicate detection pairs;
- unsupported input-band behavior.

Include one example command and its expected output. The example must run without a notebook and without editing source paths.

##### Step 8 — Reproduce TEST2

Run the CLI against the frozen TEST2 input and compare it row by row with the authoritative current-map TEST2 output.

Compare:

- row identity and order;
- sky-center and magnitude-bin assignment;
- all four densities;
- all four returned probabilities where valid;
- `valid` and exact reason code;
- derived velocity and mean magnitude;
- total reason/coverage counts.

The target is exact reproduction. If platform/library floating-point behavior prevents bit-for-bit equality, report the maximum absolute and relative difference, show that assignments/reasons are identical, and set a strict justified tolerance in the seal. Never weaken the tolerance simply to make the test pass.

Also prove that randomly permuting input rows and changing chunk size produce the same row-keyed outputs.

##### Step 9 — Report performance without optimizing science

Measure enough runtime behavior to plan Rubin nightly use:

- total rows and wall time;
- rows/second;
- peak memory;
- number and total volume of map loads;
- time spent parsing, assigning, loading, and scoring;
- behavior with many rows sharing one map versus rows spread over many maps.

This is an engineering measurement only. Do not alter map contents or score semantics for performance. Safe caching/chunking improvements are acceptable only if regression equivalence is retained.

##### Step 10 — Deliver a completion report

Return:

1. Commit hash and changed-file list.
2. Exact example CLI command.
3. Input/output schema links.
4. Scorer and map seal IDs.
5. TEST2 regression summary, including maximum numeric discrepancies.
6. Coverage and invalid-reason counts.
7. Runtime/memory/map-load summary.
8. Any blocker or ambiguity found in the old scorer.
9. An explicit statement that no maps, scientific scoring rules, digest2 code, Sorcha data, or evaluation populations were changed.

##### Acceptance checklist for Hyak

- [x] Authoritative TEST2 scorer and map state can be reconstructed exactly.
- [x] v1 contract was documented alongside implementation.
- [x] Completed validation runs are tied to the recorded immutable scorer/map seal.
- [x] Native canonical Rubin CSV runs end to end, and the VDP CLI rejects MPC-80 ingest.
- [x] MPC-80 export/parser code is isolated in the digest2 comparison tooling and preserves an auditable tracklet-ID mapping.
- [x] Default execution loads once per grouped center rather than once per tracklet; sub-center chunking reopens maps but is explicitly counted.
- [x] Output preserves every input `tracklet_id` and original row order.
- [x] All densities, probabilities, assignments, validity states, and reasons are returned/auditable.
- [x] No silent clipping, row dropping, zero filling, band conversion, or “latest map” selection occurs in the native CSV path.
- [x] Frozen TEST2 results are reproduced exactly under the precomputed isolation path; the 295 observable-geometry changes are quantified separately.
- [x] Results are invariant to input order and chunk size.
- [x] Documentation states that current four-density coverage is approximately 43.96% and that v1 is not yet an operational-readiness claim.
- [x] Example files, tests, validation artifacts, and a completion report were delivered for the completed run.

### Priority 2 — Evidence notebook: prove where legacy and new differ

Target notebook:

`notebooks/validation/neomod3_new_vs_legacy_score_regime.ipynb`

Use the same existing rows and make direct legacy-versus-new plots. This notebook should establish the phenomenon before the mechanism notebook tries to explain it.

Required views:

- Full ROC and carefully labeled zooms; no unsupported “low-FPR” label.
- `TPR_new - TPR_legacy` evaluated on a common FPR grid.
- Precision versus NEO completeness across at least the 20–50% operating range.
- Score histograms/ECDFs by true population for new and legacy.
- Paired per-row score difference, including NEO magnitude, velocity, and sky-direction slices.
- AUC and partial-AUC definitions written explicitly.
- Paired bootstrap intervals for the main differences.
- Curve crossing points and threshold/tie behavior.
- Coverage/abstention accounting so a common-scorable restriction cannot be mistaken for whole-sample performance.

The notebook must contain plots that visibly show the legacy-versus-new differences. Its first conclusion must prove or reject whether the performance gap occupies a low-FPR/top-ranked operating regime. “Low FPR” describes **where** legacy is better; it is not itself the mechanism explaining **why**. The notebook should end with a factual description of the affected regime and a list of mechanism hypotheses ranked for later testing.

#### Detailed Hyak work order for the proof notebook

##### Purpose and stopping point

Build one self-contained notebook that establishes, visually and statistically, where the new and legacy classifiers differ on the already completed TEST2 scores. It must not rebuild maps, rescore tracklets, generate a new population, tune thresholds, or test smoothing/k/bin-width mechanisms. It stops after stating the affected operating regime and listing—not testing—mechanism hypotheses.

##### Reuse before writing new plotting code

Use the established metric/plot conventions from:

```text
mag245_nbody_benchmark_vs_sorcha_roc.ipynb
neomod3_mag245_roc_sky_magnitude.ipynb
```

In particular, preserve the distinction already documented there:

- ROC = true-positive rate/completeness versus false-positive rate.
- Completeness–contamination = recall versus `1 - precision`; this is not an ROC curve.
- **Superseded 2026-08-26.** Existing code used standardized partial AUC at `max_fpr=0.01`. The proof notebook reports **unnormalised raw partial area over each explicitly stated FPR interval** instead. McClish standardization is dropped and must not be reintroduced — see `EVALUATION_PROTOCOL.md` v1.4 section 4.1. `raw/L` survives only as "historical normalized pAUC" in the sealed-headline reproduction table.

##### Compaction/resume read order on Hyak

If the Hyak agent loses context or compacts mid-build, it must not continue only from the temporary builder script. Before resuming, read the latest repository state and these documents in order:

1. `docs/gorica_24_aug_ legacy_vs_new.md` — current scientific priorities, notebook contract, and decisions. If it is not present on Hyak, stop and request/sync the current copy rather than reconstructing it from older notes.
2. `docs/NEOMOD3_MAG025_K150_FULLGRID_RUNBOOK.md` — current Hyak environment, map-build/job layout, paths, and operational run instructions for the map set used by TEST2.
3. `docs/EVALUATION_PROTOCOL.md` — authoritative GEN/CAL/TEST roles, normalization rules, invalid-row policy, metric definitions, and compute-node requirement.
4. `docs/RUBIN_TRACKLET_SCORER_V1.md` — frozen current scorer schema, seals, geometry behavior, regression result, and coverage limitation.
5. `docs/CALIBRATION_STAGE_MANIFEST.md` and `docs/E1_INTERPOLATION_ABLATION_MANIFEST.md` — read only to preserve prior CAL/TEST decisions and avoid reopening tuned choices.

`docs/HANDOFF.md`, `docs/NEOMplanHYAK.md`, and older dated handoffs may contain useful environment history, but they are secondary and must not override the newer runbook, frozen evaluation protocol, or this living plan.

At resume, run read-only checks first:

```text
git status --short
git log -5 --oneline
git diff --stat
git diff --cached --stat
squeue/sacct for any notebook jobs already submitted
```

Then inspect `.knnexp/build_regime_nb.py`, the generated notebook, and any existing output directory. Re-run the all-code-cell syntax check before notebook execution; the captured pre-compaction transcript ended while fixing an unterminated digest2 annotation string in generated cell 33, so current on-disk state—not the transcript—decides whether that issue remains.

##### Proposed inputs

Use the immutable existing Hyak artifacts, expected under:

```text
outputs/test2_geometric/TEST2_SCORED.parquet
outputs/test2_geometric/TEST2_TRACKLETS.parquet
outputs/test2_geometric/TEST2_V1_REGRESSION.json
neomod/seals/RUBIN_TRACKLET_SCORER_V1_SEAL.json
neomod/seals/RUBIN_TRACKLET_SCORER_V1_INTERFACE.json
```

At the top of the notebook, inventory the actual column names and write input hashes, Git commit, scorer/map seals, row counts, and weight semantics to a machine-readable receipt. Do not silently guess score, truth, weight, or validity columns.

##### Row universes

Create and name the following masks before computing metrics:

1. `ALL_TEST2`: all 688,688 TEST2 rows, used for coverage accounting only.
2. `PAIRED_NEW_LEGACY`: rows with finite/valid new and legacy scores. This is the **primary** proof population.
3. `COMMON_THREEWAY`: rows with finite/valid new, legacy, and digest2 scores. Use only for panels containing digest2.

Do not restrict the primary legacy-versus-new comparison merely because digest2 is unavailable on a row. Report raw and physical-weighted counts by true population for every mask and every exclusion reason.

Use exactly the same rows, truth labels, and physical weights for new and legacy within a comparison. If multiple tracklets share a source parent/object, identify that grouping so bootstrap resampling can occur by parent rather than treating correlated rows as independent.

##### Metric definitions to print in the notebook

```text
TPR / completeness = weighted true positives / weighted NEO total
FPR                = weighted false positives / weighted non-NEO total
precision          = weighted true positives / weighted selected total
contamination      = 1 - precision
```

Report:

- full ROC AUC;
- **unnormalised raw partial area** over each predeclared FPR interval, always quoted with its interval (max value is `L`; never compare across intervals);
- TPR at fixed FPR;
- precision/contamination at fixed completeness;
- completeness at fixed contamination;
- score medians and selected quantiles by true population;
- paired new-minus-legacy differences.

Use physical weights for headline metrics and show raw counts alongside them. If weighting does not affect a metric because weights are constant within a class, state that rather than silently dropping weights.

##### Predeclared operating ranges

Do not choose one favorable zoom after seeing the curves. Evaluate the full ROC plus a fixed grid of candidate operating regions:

```text
FPR <= 0.001
FPR <= 0.005
FPR <= 0.010   # matches the older max_fpr convention
FPR <= 0.020
FPR <= 0.050
```

For logarithmic FPR plots, the left edge must be tied to the empirical resolution `1 / weighted_nonNEO_count`; do not plot or interpolate through log(0). Show reachable empirical operating points as steps/markers. Interpolated fixed-FPR summaries must be labeled as interpolation.

The notebook should call the result “low-FPR” only if legacy’s paired advantage is concentrated and statistically supported in the predeclared low-FPR windows. Otherwise describe the actual range found.

##### Required notebook sections and figures

1. **Provenance and assertions**
   - input paths/hashes, commits, seals, software versions;
   - score/truth/weight/validity columns;
   - unique tracklet IDs and join integrity;
   - assertions that row order and paired identities match.

2. **Coverage and row-flow table**
   - all rows → new valid → legacy valid → paired → three-way;
   - raw and weighted totals by NEO/MBA/TNO/Trojan;
   - exclusion reasons;
   - prominent warning that paired performance is not whole-domain coverage.

3. **Score-distribution proof**
   - weighted histograms and ECDFs for NEO and non-NEO;
   - new and legacy overlaid on identical rows;
   - paired scatter/hexbin of `score_new` versus `score_legacy`;
   - paired `score_new - score_legacy` distribution by true population;
   - median and upper-tail quantiles, explicitly checking the previously observed NEO median shift (about 0.280 versus 0.508) and nearly unchanged non-NEO medians.

4. **Full conventional ROC**
   - new and legacy on `PAIRED_NEW_LEGACY`;
   - digest2 only in a separate `COMMON_THREEWAY` panel;
   - full AUC and paired delta.

5. **Low-FPR candidate windows**
   - linear zooms at each predeclared limit or a compact multi-panel figure;
   - one log-FPR view with empirical resolution shown;
   - raw pAUC table, one row per interval, with each interval's ceiling shown;
   - no curve smoothing that invents operating points.

6. **Difference curves**
   - `TPR_new - TPR_legacy` on a common FPR grid, with zero reference;
   - a table of TPR differences at the predeclared FPR values;
   - bootstrap confidence bands/intervals.

7. **Old operational plot style**
   - completeness versus contamination/`1 - precision`, explicitly labeled “not ROC”;
   - focus on 20%, 30%, 40%, and 50% completeness;
   - reproduce/check the existing precision values that suggested legacy is cleaner throughout this band;
   - completeness at 5% and 10% contamination;
   - show reachable thresholds and ties rather than a visually smoothed fiction.

8. **Threshold/top-ranked view**
   - selected weighted candidate count and contamination as score threshold changes;
   - top-ranked quantile or top-N comparison on identical rows;
   - identify whether the gap is due to NEO scores moving down, non-NEO upper tails moving up, or both—descriptively, without assigning a map mechanism.

9. **Descriptive stratification appendix**
   - magnitude, speed, and coarse sky-center/antisolar-direction slices using existing columns only;
   - enforce declared minimum positive/negative counts;
   - mark insufficient cells rather than computing unstable AUCs;
   - do not yet build the advisor’s full 5° × 5° direction mosaic; that remains the next localization deliverable.

10. **Conclusion with a forced separation of claims**
    - `WHERE`: exact FPR/completeness/contamination range in which legacy is better, tied, or worse;
    - `CONFIDENCE`: paired uncertainty and whether the difference excludes zero;
    - `COVERAGE`: the row universe to which the claim applies;
    - `NOT YET WHY`: list excessive effective smoothing, k, bin width, interpolation, source choice, and density-validity policy only as hypotheses for the later mechanism notebook.

##### Paired uncertainty

Use a fixed random seed and a stratified paired bootstrap, keeping each row’s new/legacy scores together. Resample source parent/object clusters if repeated parents exist; otherwise resample within truth/population strata. Produce at least 500 replicates, increasing to 1,000 if runtime is reasonable.

Bootstrap at least:

- delta full AUC;
- delta raw pAUC over the intervals [0,0.001], [0,0.005], [0,0.01], [0,0.02], [0,0.05];
- delta TPR at the same fixed-FPR points;
- delta precision at 20%, 30%, 40%, and 50% completeness;
- delta completeness at 5% and 10% contamination.

Report the point estimate, 2.5/50/97.5 percentiles, and fraction of replicates favoring legacy. Do not turn the bootstrap into a new tuning loop.

##### Outputs

Create:

```text
notebooks/validation/neomod3_new_vs_legacy_score_regime.ipynb
outputs/neomod3_new_vs_legacy_score_regime/coverage.csv
outputs/neomod3_new_vs_legacy_score_regime/global_metrics.csv
outputs/neomod3_new_vs_legacy_score_regime/fpr_operating_points.csv
outputs/neomod3_new_vs_legacy_score_regime/completeness_operating_points.csv
outputs/neomod3_new_vs_legacy_score_regime/score_quantiles.csv
outputs/neomod3_new_vs_legacy_score_regime/bootstrap_intervals.csv
outputs/neomod3_new_vs_legacy_score_regime/RUN_RECEIPT.json
outputs/neomod3_new_vs_legacy_score_regime/figures/*.png
```

Keep large row-level intermediates outside Git. Commit the notebook, compact tables, receipt, and figures according to the existing repository convention.

##### Git staging after successful verification

After the notebook executes successfully and all acceptance checks pass, stage only durable deliverables:

```text
git add notebooks/validation/neomod3_new_vs_legacy_score_regime.ipynb
git add outputs/neomod3_new_vs_legacy_score_regime/coverage.csv
git add outputs/neomod3_new_vs_legacy_score_regime/global_metrics.csv
git add outputs/neomod3_new_vs_legacy_score_regime/fpr_operating_points.csv
git add outputs/neomod3_new_vs_legacy_score_regime/completeness_operating_points.csv
git add outputs/neomod3_new_vs_legacy_score_regime/contamination_operating_points.csv
git add outputs/neomod3_new_vs_legacy_score_regime/score_quantiles.csv
git add outputs/neomod3_new_vs_legacy_score_regime/bootstrap_intervals.csv
git add outputs/neomod3_new_vs_legacy_score_regime/stratification.csv
git add outputs/neomod3_new_vs_legacy_score_regime/RUN_RECEIPT.json
git add outputs/neomod3_new_vs_legacy_score_regime/figures
```

Stage any additional compact table actually produced only after documenting it. Do not stage `.knnexp/`, Slurm logs, executed temporary copies, large row-level parquets/CSVs, maps, or unrelated dirty files. If the builder is needed for reproducibility, move a cleaned version from `.knnexp/` into a deliberate repository path before staging; never stage the scratch directory itself.

Show `git status --short` and `git diff --cached --stat` after staging. Do not commit or push unless separately authorized.

##### Acceptance criteria

- Notebook executes top to bottom from the sealed existing inputs.
- No map build, rescoring, new population, threshold tuning, digest2 modification, or mechanism ablation occurs.
- Primary new-versus-legacy metrics use identical paired rows and weights.
- Coverage/abstention is visible and never conflated with performance.
- Known headline metrics are reproduced within a stated tolerance before new plots are trusted.
- Every ROC/pAUC/completeness/contamination definition is explicit.
- The notebook proves or rejects the low-FPR/top-ranked hypothesis with paired uncertainty.
- The conclusion states where legacy differs and explicitly defers why to the mechanism notebook.

### Priority 3 — Resolve incomplete-population scoring and coverage

This follows the proof notebook and has two linked but distinct deliverables.

**A. Compare operational binary formulations without lying about probability semantics.** On the same existing tracklets, compare:

- strict `p_NEO_4pop` where all four component estimates are valid;
- `score_NEO_MBA` wherever NEO and MBA estimates are valid;
- a direct physically weighted `p_NEO_vs_OTHER` pilot, where `rho_OTHER` is estimated as the pooled non-NEO mixture rather than a sum that fails when a component map is missing.

Report:

- coverage of each formulation;
- agreement of direct `p_NEO_vs_OTHER` and component-sum `p_NEO_4pop` on the all-four valid subset;
- performance on identical matched subsets;
- whether pooled MBA abundance washes out true TNO/Trojan clumps;
- where true TNO/Trojan rows receive dangerously high NEO–MBA scores;
- sky, magnitude, and velocity locations of fallback-only rows.

This is a diagnostic/scorer-design notebook first, not an immediate change to frozen v1.

**B. Repair the maps/mixture estimate.** Run the selected-center population-specific global-cloning pilot with explicit physical weights before considering arbitrary floors. Determine which missing cells are finite-realization holes filled by proper propagation/cloning and which require adaptive magnitude widths, neighboring-center pooling, a direct pooled-OTHER estimator, or a continuous magnitude/density model.

Success means that intended in-domain tracklets receive a valid NEO-versus-all-nonNEO estimate at a declared high coverage target, while physical integrals and mixture ratios remain correct. Component-level `missing_population_density:*` may remain a useful diagnostic but must not dominate the operational follow-up output.

### Priority 4 — Recover the old direction-specific ROC machinery

Locate and reuse the older `mag245` n-body notebook/code that already generated 5° × 5° sky-direction ROC tile plots for a single night. Record its exact path, coordinate convention, tiling rule, row selection, and minimum-class-count requirement before adapting it.

There is an ambiguity to resolve with the advisor: “5×5” might mean a single 5° × 5° direction patch, or a 5-by-5 mosaic of 25 small directions. The old implementation is evidence for the former, but we should confirm the requested deliverable.

Before producing a mosaic, run per-tile feasibility counts. The common-scorable sample has 14,026 NEOs and will not populate 25 directions evenly; any tile without enough positives and negatives must be labeled insufficient rather than assigned a misleading ROC.

The final localized comparison should place new VDP, legacy VDP, and black-box digest2 on identical rows in each viable region.

### Priority 5 — Short digest2 resolution/plotting audit

This is expected to be a small audit, not a source-code project.

- Compare the raw digest2 output values used in the old and current notebooks.
- Count unique score values and verify whether the official score is the same integer 0–100 quantity.
- Check parser/version/config provenance for both runs.
- Compare plotting calls. The current strong hypothesis is that old ROC plots looked smooth because `roc_curve` operating points were joined by a line, while contamination/completeness plots correctly showed discrete dots and ties.
- Plot the actual discrete operating points and ties clearly.

No digest2 code modification, hidden-score extraction, or recalibration is authorized. It remains a black box.

### Priority 3 implementation detail — Selected-center population-cloning coverage pilot

Choose a small set of centers that spans:

- high and low four-density coverage;
- different sky directions;
- common and sparse magnitude slices;
- representative `no_split_fraction` failure modes.

Run population-specific global pre-selection cloning at those centers, with explicit weights and more than one clone-factor level where useful. Determine whether clean physical cloning makes the narrow 0.25-mag four-density maps adequately complete before designing a more complex pooling model.

This coverage work comes before expensive full mechanism builds because coverage is the current operational blocker and a successful cloning policy may also change the apparent legacy/new mechanism.

### Priority 6 — Contained mechanism-ablation notebook

Target notebook:

`notebooks/validation/neomod3_new_vs_legacy_mechanism_ablation.ipynb`

Run on selected centers only, using identical source rows, weights, and diagnostic tracklets. Change one component at a time:

1. `k_NEO = 10` versus `k_NEO = 150` at fixed bin width.
2. 0.25-mag versus 1.0-mag bins at fixed `k`.
3. Hard containing-bin lookup versus magnitude interpolation.
4. Gaussian smoothing off/on using the now-verified legacy NEO-only sigma of 0.04 deg/day.
5. Current four-density validity versus the legacy zero-fill behavior, only as a labeled semantic diagnostic—not as an acceptable probability policy.
6. BASE versus HIGH source choice where source files permit a paired comparison.
7. Current versus population-specific cloned non-NEO maps after the cloning pilot defines a defensible configuration.

For every row report density integrals, effective neighbor radii, coverage, NEO/non-NEO score distributions, and paired performance changes. The goal is to disprove alternatives and identify which component causes the observed score-regime change.

### Priority 7 — Choose map v2, then final evaluation

Only after the evidence, coverage pilot, and mechanism ablation:

1. Select and preregister the next map policy.
2. Build the necessary full 667-center production maps.
3. Freeze a new map seal under the unchanged scoring interface contract.
4. Obtain a genuinely fresh final evaluation realization/partition if a publishable out-of-sample claim is needed.
5. Run the final comparison and only then proceed to downstream Sorcha/nightly operational validation.

## Compute and data constraints

- One full 667-center build took roughly 3.3 hours and 248 GB.
- Five full-grid ablation variants would require about 1.2 TB and do not fit comfortably in the available space.
- Therefore, diagnostic ablations start with selected centers. A full build occurs only for the selected policy.
- No untouched S3M non-NEO TEST partition remains: the earlier TEST work consumed every eligible TEST-split S3M parent for the relevant populations.
- Reused exposed data are acceptable for paired mechanism diagnosis, not for a fresh performance claim.
- We will use the older/existing data now and will not draw a new population merely to continue tuning.

## What not to do yet

- Do not launch another full-grid map matrix.
- Do not begin a new Sorcha production run with the incomplete maps.
- Do not generate a fresh final evaluation population while map choices are still being inspected.
- Do not tune map choices and then describe TEST2 as untouched validation.
- Do not replace missing density populations with zero or an arbitrary epsilon and call the result `P(NEO)`.
- Do not modify digest2.
- Do not assert that the mechanism is Gaussian smoothing, high `k`, narrow bins, or a “low-FPR problem” until the evidence and ablation notebooks demonstrate it.
- Do not claim the frozen v1 interface is Rubin-ready while its declared in-domain complete-density coverage remains about 44%.

## Decision gates and open questions

### Needed from Hyak for v1

- What is the authoritative current commit and working tree state?
- What exact scorer entry point produced TEST2?
- Where are the map roots, manifest, and seal?
- What source/config/dependency hashes must enter the v1 seal?
- Are the full frozen TEST2 score outputs accessible for regression testing?
- What exact RA/Dec frame does the Rubin/canonical input contract require?
- Is center assignment based on detection 0, the tracklet midpoint, or a position at the sealed map epoch?
- How is the map epoch connected to each tracklet’s MJD, and what date-validity window does a sealed map set have?

### Needed for the evidence/localization work

- What is the exact path of the old `mag245` 5° × 5° ROC implementation?
- Did the advisor mean one 5° × 5° patch or a 5-by-5 collection of 25 directions?
- What minimum NEO/non-NEO counts will be required for a reported local ROC?

### Needed for map design

- What are the measured production per-cell NEO neighbor radii under legacy and new configurations?
- Which centers and bins span the relevant coverage regimes for the pilot?
- Which orbital elements/phases may be varied in a physically justified clone model?
- What clone-factor grid should be tried for MBA, TNO, and Trojan?
- How much of `no_split_fraction` is repaired by global pre-selection cloning?
- If holes remain, should the next test use adaptive magnitude width, neighboring-cell pooling, or a continuous generative density?

## Near-term handoff message for Hyak

> Task 1 is accepted and closed at commits `988e82c` plus interface cleanup `6277cff`; do not rerun it. The immediate next deliverable is the existing-data proof notebook `notebooks/validation/neomod3_new_vs_legacy_score_regime.ipynb`. It must show, rather than assume, whether legacy’s advantage occurs in a low-FPR/top-ranked operating regime: full and zoomed ROC, `TPR_new - TPR_legacy` on a common FPR grid, precision versus 20–50% completeness, paired score distributions, bootstrap uncertainty, curve crossings/ties, and coverage accounting on identical rows. This establishes **where** legacy is better; it must not claim smoothing/k/bin width as the cause. After that notebook, begin the coverage/scorer diagnostic comparing strict four-component `p_NEO_4pop`, conditional `score_NEO_MBA`, and a directly estimated physically weighted `p_NEO_vs_OTHER = rho_NEO/(rho_NEO+rho_OTHER)`. A direct OTHER mixture is the operational target to test; simply summing current missing component maps does not solve coverage. Then run the selected-center cloning/weighting pilot. No new evaluation population, full-grid build, zero filling, or digest2 modification.

## Deliverables checklist

- [x] Authoritative Hyak scorer/map inventory recorded.
- [ ] Rubin scorer v1 input/output contract final after CSV-only/MPC separation.
- [x] Validation seal created with source, map, config, and environment provenance.
- [x] Native-CSV CLI implemented with grouped-by-map loading and deterministic output.
- [x] Full TEST2 regression comparison passed under precomputed geometry; operational geometry difference quantified.
- [x] MPC-80 input/export code isolated from the Rubin VDP CLI without rerunning the accepted full regression (`6277cff`).
- [ ] Legacy-versus-new score-regime proof notebook complete; low-FPR/top-ranked hypothesis proved or rejected.
- [ ] NEO+MBA conditional-score coverage/safety diagnostic complete.
- [ ] Direct physically weighted NEO-versus-OTHER mixture pilot complete and cross-checked against the four-component sum.
- [ ] Selected-center TNO/Trojan cloning coverage pilot complete.
- [ ] Future LSST-band photometry strategy selected and validated before operational Rubin use.
- [ ] Old direction-specific ROC implementation located and documented.
- [ ] Per-tile feasibility counts complete.
- [ ] digest2 plotting/quantization audit complete.
- [ ] Selected-center cloning coverage pilot preregistered.
- [ ] Cloning pilot complete and coverage decision recorded.
- [ ] Selected-center mechanism ablation complete.
- [ ] Future map-v2 policy selected and preregistered.
- [ ] Fresh final evaluation reserved only after policy freeze.

## Changelog

### 2026-08-24 — Initial synthesis

- Consolidated the advisor notes, TEST2 findings, coverage problem, legacy smoothing metadata, current bandwidth hypothesis, digest2 constraint, Rubin interface priority, cloning interpretation, compute limits, and evaluation-data policy.
- Recorded that proof of the legacy-versus-new score regime belongs in its own notebook before causal claims.
- Recorded global pre-selection orbital cloning with per-population physical weights as the first coverage remedy to test.
- Recorded that existing/older data will be used for diagnosis and no new evaluation sample will be drawn yet.

### 2026-08-24 — Rubin scorer v1 feasibility update

- Recorded that observed RA/Dec can be converted to ecliptic position/rates; only TEST2’s cached model-frame coordinates are unavailable to a real tracklet.
- Recorded the 295/688,688 observable-versus-cached center-assignment differences and the need for dual regression reporting.
- Corrected the MPC-80 interpretation: its tested magnitude representation loses precision through 0.1-mag rounding, causing 12.56% of the 5,000 round-tripped examples to cross a hard 0.25-mag boundary.
- Added open decisions for the tracklet reference position, coordinate frame/time, map epoch, and map date-validity window.

### 2026-08-25 — Corrected MPC-80 scope

- Corrected the Rubin interface boundary: frozen VDP scorer v1 consumes the native full-precision Rubin tracklet CSV only.
- MPC-80 is created solely for the separate black-box digest2 comparison and is never fed back into VDP scoring or VDP magnitude-bin assignment.
- Reclassified the measured MPC-80 rounding/bin changes as digest2-input provenance, not a limitation of the native Rubin VDP path.

### 2026-08-25 — Task 1 validation completed on Hyak

- Recorded Hyak commit `988e82cea3a3200de630728fe214eee43a340ab9`, scorer seal `0e61ec37...`, and map seal `a010bb1b...`; the commit was reported clean and not pushed.
- Recorded the completed 688,688-row regressions: exact precomputed-geometry reproduction and the quantified 295-row observable-geometry reassignment.
- Recorded 100% native-CSV magnitude-bin reproduction, 43.9579% operational coverage, reason counts, performance, and exact order/chunk invariance.
- Recorded wrapper fixes for out-of-range magnitude reasons, malformed-coordinate validation, input hashing, and honest map-open accounting.
- Accepted all expensive compute as final for Task 1. Only the MPC-80 interface separation, documentation correction, and associated lightweight provenance update remain; no full regression rerun is requested.

### 2026-08-25 — Task 1 interface cleanup completed

- Recorded follow-up commit `6277cff`: native-CSV-only VDP CLI, separate digest2 MPC-80 adapter, immutable validation seal preserved, interface metadata added, and all scientific source hashes unchanged.
- Recorded exact scientific-column equivalence on the lightweight smoke test and closure of Task 1 without another full TEST2 run.

### 2026-08-25 — Coverage semantics and Rubin-band roadmap

- Distinguished valid near-zero/zero population density from an unavailable density estimate; current missing TNO/Trojan cells do not prove physical absence.
- Kept the strict four-density posterior semantics, but promoted an explicitly labeled NEO-versus-MBA conditional score as the next operational diagnostic so useful rows are not reduced to a reason code alone.
- Made sparse-population coverage repair the immediate next scientific priority, beginning with population-specific global cloning and explicit weights on selected centers.
- Recorded apparent V as a temporary frozen-map contract and added future per-LSST-band maps versus uncertainty-aware band-to-V conversion as a required Rubin-readiness decision.

### 2026-08-25 — Restored proof-notebook priority and added OTHER mixture

- Restored the legacy-versus-new score-regime proof notebook as the immediate next task. It must prove or reject the predicted low-FPR/top-ranked location of the gap before mechanism or coverage work proceeds.
- Distinguished the location question (“where is legacy better?”) from the later mechanism question (“why is it better?”).
- Added `p_NEO_vs_OTHER = rho_NEO/(rho_NEO + rho_OTHER)` as the preferred operational binary formulation to test when `rho_OTHER` is estimated directly from the complete, physically weighted non-NEO mixture.
- Recorded that summing the existing component maps is equivalent to the strict four-population denominator and therefore does not repair missing-component coverage.

### 2026-08-31 — k_OTHER sensitivity and poorly-covered-centre coverage test

Both studies are **contained selected-centre diagnostics**, not global performance claims. The
earlier pilot conclusion keeps that label.

#### Controlled-change check

`--k-other` was added to `pipeline/build_neomod3_pneo_other_pilot.py`. k = 10 rebuilt under the
modified builder reproduces the original pilot maps **bit-for-bit** (identical key sets, max
absolute difference 0.0 across every array), so the flag is behaviour-preserving and `k_OTHER` is
genuinely the only quantity varied. Asserted in code; the study aborts if it fails.

#### k_OTHER = 10 versus 30

Same three centres, same `2027-08-25T00:00:00` epoch, same 6,866 TEST2 rows, same NEO numerator
copied from the frozen maps, same per-sample physical weights.

**Row coverage is unchanged by k:** 6,852/6,866 (99.80%) at both settings. Valid *bins*
differ (40-41 at k=10 versus 33-35 at k=30, since support requires n > k), but those extra bins hold
few TEST2 rows.

**The NEO score shift reverses sign with k, and the four-population reference is bracketed:**

| | old `p_NEO_4pop` | k_OTHER = 10 | k_OTHER = 30 |
|---|---:|---:|---:|
| NEO median | **0.390248** | 0.420513 | 0.313013 |
| NEO median old − new | — | **-0.013525** | **+0.054042** |
| NEO mean old − new | — | -0.028143 | +0.052565 |

k = 10 puts NEO scores **above** the four-population reference; k = 30 puts them **below**; the
reference lies between. This **confirms the bandwidth hypothesis recorded on 2026-08-31** as the
cause of the NEO shift: pooling roughly triples the sample count, so a smaller k means a narrower
kernel, a lower `rho_OTHER`, and a higher `p_NEO`. It follows that some intermediate k would null
the shift. That value was **not** searched for — this was a bounded two-point comparison, not a
tuning loop.

**Agreement with the frozen four-population reference is much better at k = 10:**

| metric (paired rows, n = 6,484) | k = 10 | k = 30 |
|---|---:|---:|
| global Spearman vs old | **0.998707** | 0.965801 |
| global Pearson vs old | 0.997800 | 0.996510 |
| median abs diff vs old | **1.739e-07** | 5.398e-05 |
| MBA Spearman vs old | 0.999958 | 0.960241 |
| TNO Spearman vs old | 1.000000 | 0.907115 |
| Trojans Spearman vs old | 0.988983 | 0.973743 |

Spearman between the two settings themselves is 0.963875.

**Provisional pilot setting: `k_OTHER = 10`.** Retained because it agrees most closely with the
frozen four-population reference, is identical to the k every frozen non-NEO population already
uses (minimal change), yields more valid bins, and costs no coverage.

**`k_OTHER = 10` is NOT physically final and NOT calibrated.** It is the provisional setting of a
five-centre pilot, chosen from a two-point comparison against a reference that is itself a
particular bandwidth choice. Nothing here establishes a physically preferred neighbour count, and
the bracketing result above shows the NEO score value depends on k. Any future map version must
select k_OTHER on its own evidence rather than inheriting this number.

Results: `outputs/neomod3_pneo_other_pilot/k_other_sensitivity/`.

#### Does the coverage gain persist at poorly covered centres?

**A structural finding first.** `center_valid_counts_all.csv` lists 341 centres, but TEST2 touches
**667**. The missing **326 centres have zero old four-population valid rows** and therefore never
appear in that file at all. Roughly half the sky is completely unscorable under the four-population
rule, and a file keyed on valid counts cannot show it.

"Poorly covered" therefore has to mean low valid *fraction*, not low count. Selection rule: lowest
old-valid fraction among centres present in the counts file with at least 200 TEST2 rows, ties by
ascending id. Selected `dlon-060_lat+00` (10.5% old-valid) and
`dlon-060_lat-05` (10.9% old-valid), against 94-95% for the first three centres.

Same method, same epoch, same NEO numerator, same weights, `k_OTHER = 10`.

| centre | rows | old valid | new valid | delta |
|---|---:|---:|---:|---:|
| `dlon-060_lat+00` | 996 | 105 (10.54%) | 977 (98.09%) | **+87.55 pp** |
| `dlon-060_lat-05` | 1,468 | 160 (10.90%) | 1,448 (98.64%) | **+87.74 pp** |
| **total** | **2,464** | **265 (10.75%)** | **2,425 (98.42%)** | **+87.66 pp** |

**The gain persists and is far larger than at the best-covered centres: +87.7 percentage points
versus +5 pp.** 2,160 rows recovered, **0 lost** — the new formulation again strictly
dominates.

By true population (no Trojans fall at these centres):

| truth | rows | old valid | new valid | recovered |
|---|---:|---:|---:|---:|
| NEO | 76 | 10.53% | 96.05% | 65 |
| MBA | 2,339 | 10.65% | 98.46% | 2,054 |
| TNO | 49 | 16.33% | 100.00% | 41 |

The old scorer abstained here overwhelmingly on `missing_population_density:Trojans` (2,028 rows) —
the sparse Trojan map — plus `TNO,Trojans` (135). Pooling absorbs both. Under the new formulation
only 39 rows remain invalid: 36 lacking a NEO density and 3 outside the velocity grid.
**`OTHER` itself never failed.** Nothing was zero-filled and no epsilon floor was used.

Agreement holds on the 265 rows valid under both: median |old − new| 1.911e-06, Spearman
0.999162.

Results: `outputs/neomod3_pneo_other_pilot/poor_centers/`.

#### Scope, unchanged

Five centres out of 667. These diagnostics support statements about denominator agreement, about
the bandwidth origin of the NEO shift, and about coverage repair at the centres tested. They remain
**selected-centre diagnostics** and support no global performance, ROC, or survey-wide coverage
claim. No full-grid build, no new evaluation population, no digest2 change, no v1 scorer change, no
zero-filling, no epsilon floor.

### 2026-08-31 — Current maps frozen; selected-centre NEO-vs-OTHER pilot complete

**Map freeze.** `docs/NEOMOD3_CURRENT_MAPS_ARCHIVE.md` and
`outputs/neomod3_current_maps_archive/current_maps_manifest.json` record the map set that produced
the accepted TEST2 results. All 265.86 GB were re-hashed from disk and compared with the per-centre
hashes in `MAP_BUILD_SEAL_V2.json`: **667/667 match, 0 mismatched, 0 missing.** The map root was not
modified.

- map root `/mmfs1/gscratch/dirac/ds2004/sorcha/outputs/neomod3_mag025_k150_maps_v2`
- **map reference epoch `2027-08-25T00:00:00`** (built 2026-08-17T06:46:31Z, archived 2026-08-31T10:39:51Z — three distinct timestamps, recorded separately)
- 667 centres, 44 bins, k = NEO 150 / others 10, smoothing False, masking False
- velocity grid [-5.0, 5.0] at 0.01; bins 14.0 <= V < 25.0 step 0.25, half-open
- map seal `a010bb1bf2e421a3258ff007f7bc4034369918a3877418620cd14e3f0b9d19ab`
- scorer seal `0e61ec37e263e9a4f9e81ed4ef25ef8fa069f630b8d2003f2ce29ce37de9eaeb`
- map-build commit `e4780ef806d4eb4a77e3d384144f8a29ff24ea1e`; sealed module `a6de18c2197cbcb9f93014712aaad232272578150e2222f28eae723b339f923a`
- the seal chain was verified to close on itself: scorer seal -> map seal -> split provenance -> epoch cache -> split manifest

**Measured fact that shaped the pilot design.** The three non-NEO populations do **not** share a
split fraction: all 42 populated bins differ (MBA ~0.60, TNO 0.52-0.74, Trojans 0.49-0.64). A pooled
non-NEO estimate must carry a per-sample physical weight `w_i = 1/f_split(pop, bin)`; concatenating
points and dividing by one scalar would be biased.

**OTHER estimator (new, contained to the pilot).** The sealed closed-form kNN estimator is used
unchanged and multiplied by the local mean physical weight of the same k neighbours:

```text
rho_OTHER(x) = n0(x) * mean_{j in kNN(x)} w_j ,   n0 = (k(k+1)/2 - 1/2) / (pi * sum_j d_j^2)
```

With homogeneous weights this is exactly `n0/eff`, identical to the frozen single-population line,
so the pilot estimator is a strict generalisation rather than a different estimator. `rho_NEO` is
copied verbatim from the frozen maps, so the entire change lives in the denominator.
`k_OTHER = 10`, the same k every frozen non-NEO population uses.

**Centres** (rule: top 3 by old-valid count, ties by ascending id; no hand-picking; no tie occurred):

| rank | centre | old-valid | rows | frac |
|---|---|---:|---:|---:|
| 1 | `dlon-130_lat+08` | 2,212 | 2,336 | 0.9469 |
| 2 | `dlon-120_lat+08` | 2,174 | 2,288 | 0.9502 |
| 3 | `dlon+130_lat+08` | 2,098 | 2,242 | 0.9358 |

**Selection caveat, recorded so it is not forgotten:** this rule necessarily selects the
best-covered centres (94.4% old-valid here against 43.96% globally). The pilot therefore
demonstrates **agreement** and understates coverage repair. A coverage-repair demonstration needs
poorly covered centres.

**Coverage result — the new formulation strictly dominates on these rows.**

| | old four-population | new NEO-vs-OTHER |
|---|---:|---:|
| valid rows | 6,484 / 6,866 (94.44%) | 6,852 / 6,866 (99.80%) |
| old-valid but new-invalid | -- | **0** |
| new-valid but old-invalid | -- | **368** |

All 368 recovered rows previously failed with exactly
`missing_population_density:TNO` — the failure mode OTHER was designed to remove. Recovered rows by
truth: MBA 324, Trojans 26, NEO 18. Nothing was zero-filled and no epsilon floor was used; the
14 rows still invalid fail on NEO or the velocity grid, never on OTHER.

At **bin** level the gain is larger than at row level, because TEST2 rows concentrate in
well-covered bins. Scorable magnitude slices per centre rise from the four-population intersection
to NEO-and-OTHER: 16->27, 16->26, 15->27 (**+11, +10, +12**). The binding constraint moves from TNO
to NEO.

**Agreement result.** On the 6,484 paired rows: median |old - new| = 1.739e-07,
Spearman 0.998707, Pearson 0.997800. Non-NEO medians agree to six decimals. The primary
cross-check in this plan — that `p_NEO_vs_OTHER` should agree with `p_NEO_4pop` where all four
component maps are valid — **passes**.

**The agreement is not uniform, and that matters.** The difference concentrates on true NEOs:

| truth | n | median old | median new | median old-new | max abs diff |
|---|---:|---:|---:|---:|---:|
| NEO | 228 | 0.390248 | 0.420513 | -0.013525 | 0.169828 |
| MBA | 5,327 | 0.000482 | 0.000482 | -0.000000 | 0.012266 |
| TNO | 45 | 0.000002 | 0.000002 | -0.000000 | 0.000000 |
| Trojans | 884 | 0.000407 | 0.000415 | -0.000008 | 0.000421 |

True-NEO scores move **up** (median 0.3902 -> 0.4205, mean shift -0.0281, max |diff| 0.1698) while
non-NEO scores barely move. That is a real systematic shift, not numerical noise, and large enough
to matter for ranking.

**Leading hypothesis, untested here:** pooling roughly triples the sample count in a bin, so
`k = 10` on the pooled set implies a *narrower* kernel than `k = 10` on MBA alone. Where NEOs live —
the low-non-NEO-density part of velocity space — a narrower kernel lowers `rho_OTHER` and therefore
raises `p_NEO`. This is the same bandwidth mechanism suspected in the legacy-versus-new work. It is
recorded as a hypothesis; this pilot is a map/score diagnostic and does not test it. Choosing
`k_OTHER` deliberately, rather than inheriting 10, is the obvious next contained experiment.

**Artifacts.**

```text
docs/NEOMOD3_CURRENT_MAPS_ARCHIVE.md
outputs/neomod3_current_maps_archive/current_maps_manifest.json
pipeline/build_neomod3_pneo_other_pilot.py
pipeline/score_neomod3_pneo_other_pilot.py          (separate; frozen v1 scorer untouched)
notebooks/validation/neomod3_pneo_other_pilot.ipynb
outputs/neomod3_pneo_other_pilot/selected_centers.csv
outputs/neomod3_pneo_other_pilot/center_valid_counts_all.csv
outputs/neomod3_pneo_other_pilot/old_vs_new_scores.csv
outputs/neomod3_pneo_other_pilot/coverage_by_center.csv
outputs/neomod3_pneo_other_pilot/coverage_by_population.csv
outputs/neomod3_pneo_other_pilot/score_comparison.csv
outputs/neomod3_pneo_other_pilot/score_comparison_by_population.csv
outputs/neomod3_pneo_other_pilot/RUN_RECEIPT.json
outputs/neomod3_pneo_other_pilot/pilot_maps/pilot_<centre>.npz     (3 files, ~235 MB each)
```

Pilot builder sha256 `71fa750e981e87d1e4ef123144be3ddb09b36c88da0c2ed3497eb3be4e7d3d83`;
pilot scorer sha256 `6fb4706f5900e16825a443808c4cd97a906275270578ad25d271b2a3d7516fbf`.
TEST2 inputs were re-verified against the proof-notebook receipt before use. No full-grid build, no
new evaluation population, no digest2 change, no frozen-scorer change, no zero-filling, no epsilon
floor, no ROC study.

### 2026-08-26 — Partial-AUC metric policy: raw area, McClish dropped

- **Superseded** the requirement to report McClish-standardized partial AUC. The proof notebook and all future analyses report **unnormalised raw partial area over an explicitly stated FPR interval**, plus TPR at fixed FPR, `TPR_new - TPR_legacy`, and completeness/precision operating points. `EVALUATION_PROTOCOL.md` bumped to v1.4.
- Reason: the standardization rescales the raw area onto a 0.5-baseline scale that adds no information beyond the raw area over a stated interval, and it obscures how small the absolute differences are. **Do not reintroduce it.**
- A raw pAUC over `[0, L]` has maximum `L`. Every reported value must carry its interval, and raw pAUC values over different intervals must never be compared with each other.
- `raw/L` is retained under the name **historical normalized pAUC** for one purpose only: reproducing the sealed `0.71476` / `0.72007` values that `pauc_std()` in `pipeline/test2_merge_evaluate.py` produced. It must never be called "standardized" and must not drive a conclusion.
- Recorded that `pipeline/test2_merge_evaluate.py` was the sole place computing `raw/L` under the name `pauc_std`, while twelve other pipelines used `roc_auc_score(max_fpr=...)`. Pre-v1.4 recorded results in those pipelines stand as historical record and were not recomputed.
- Corrected the global comparison table above: its third column is the historical normalized pAUC, not a standardized one.

### 2026-08-25 — Proof notebook compaction-resume instructions

- Recorded the mandatory read order after compaction, led by this living plan, `NEOMOD3_MAG025_K150_FULLGRID_RUNBOOK.md`, `EVALUATION_PROTOCOL.md`, and the Rubin scorer v1 documentation.
- Recorded that the captured transcript ended while repairing a generated-notebook string syntax error; Hyak must inspect current on-disk state and rerun code-cell compilation before execution.
- Added explicit durable-output staging commands and prohibited staging scratch builders, logs, large row-level artifacts, maps, or unrelated changes.
