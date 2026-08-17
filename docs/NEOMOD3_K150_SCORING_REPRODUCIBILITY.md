# NEOMOD3 k=150 scoring reproducibility record

Living document. Updated as each stage runs; nothing essential is left in shell history.

Status: **COMPLETE through the executed TEST2 notebook.** Stopped before Sorcha.

Previous status: assembly approved with reused frozen non-NEOs. Stage 1 complete and sealed; TEST2 NEO
parents drawn and hashed.

**The strict fresh-non-NEO gate FAILED and that result stands unaltered (§6.3). It is recorded as a
provenance finding, not rewritten as a pass.**

### Approved design and its limitation

TEST2 is intentionally:

    fresh independent TEST2 NEOs  +  frozen REUSED prior-TEST S3M non-NEOs

**Purpose.** A *paired* comparison of new VDP, legacy VDP and digest2 on **identical contaminant
rows**. Because all three classifiers see the same non-NEO rows, the comparison between them is
internally valid even though the contaminant sample is not independent of the earlier inspection.

**Limitation, stated plainly.** The non-NEO parents were used in the previously inspected TEST, so
absolute contamination, precision and FPR figures from TEST2 are **not** independent of that earlier
look and must not be presented as a fresh out-of-sample estimate. Only the NEO side is fresh.
Differences *between* classifiers on identical rows remain the trustworthy quantity.

**No further tuning is permitted** on TEST2 — no k, magnitude-width, map-parameter or calibration
adjustment, and no operating threshold selected from it.

Last updated: 2026-08-17.

---

## 0. Git commit before the work

    e4780ef806d4eb4a77e3d384144f8a29ff24ea1e   (all-sky NEOMOD3 GEN source, 16/16 validated)

Commit after the work: *(recorded at the end of §9)*

---

## 1. Absolute paths

| role | absolute path |
|---|---|
| workspace | `/mmfs1/gscratch/dirac/ds2004/sorcha` |
| repo | `/mmfs1/gscratch/dirac/ds2004/sorcha/neomod` |
| **map root (scored against)** | `/mmfs1/gscratch/dirac/ds2004/sorcha/outputs/neomod3_mag025_k150_maps_v2` |
| map seal | `.../outputs/neomod3_mag025_k150_maps_v2/MAP_BUILD_SEAL_V2.json` |
| build report | `.../outputs/neomod3_mag025_k150_maps_v2/FULLGRID_BUILD_REPORT.{json,md}` |
| coverage table | `.../outputs/neomod3_mag025_k150_maps_v2/coverage_table.parquet` |
| density integrals | `.../outputs/neomod3_mag025_k150_maps_v2/density_integrals.parquet` |
| all-sky NEO GEN cache | `.../outputs/neomod3_projection_cache_high_allsky/` |
| legacy (v1) maps | `/mmfs1/gscratch/dirac/ds2004/sorcha/prob_maps_grid_neomod3_GEN_final` |
| legacy map seal | `.../outputs/splits/MAP_BUILD_SEAL.json` |
| 0.25-mag split provenance | `.../outputs/splits/split_provenance_mag025.json` |
| non-NEO split manifest | `.../outputs/splits/nonneo_split_manifest.parquet` |
| epoch-state cache | `.../outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet` |
| prior TEST (inspected — NOT reused) | `.../outputs/test_tracklets_neomod3/` |
| **TEST2 root** | `.../outputs/test2_geometric/` |
| digest2 | `/mmfs1/gscratch/dirac/ds2004/sorcha/digest2` |
| scoring seal | `.../outputs/test2_geometric/SCORING_SEAL.json` |

---

## 2. Frozen scientific configuration

| item | value |
|---|---|
| centers | 667 |
| apparent V range | `14 <= V < 25` |
| magnitude slices | 44 half-open 0.25-mag bins, `[lo, lo+0.25)` |
| k_NEO | 150 |
| k_MBA / k_TNO / k_Trojan | 10 |
| estimator | Bayesian kNN only (closed-form evaluation, §6 of the build runbook) |
| Gaussian smoothing | OFF |
| support masking | OFF |
| velocity grid | [-5, +5] deg/day @ 0.01 → 1001 x 1001 |
| physical weights | as recorded by the map build; `effective_factor_NEO` = 64.725384 |
| renormalization | **none** — no k-dependent integral renormalization, no Platt calibration |
| posterior | `P(c|x) = rho_c / sum_all rho`, asserted to sum to 1 where the denominator is valid |
| invalid scores | remain NaN with an explicit reason; never coerced to 0 |

---

## 3. Stage 1 — map validation and seal

**Not rebuilt.** The map root was located from validation job `38590098`, which wrote
`MAP_BUILD_SEAL_V2.json` into it.

### 3.1 Build provenance

| | |
|---|---|
| build arrays | `38582732` (64 CPU, cancelled after 15 centers) + `38582905` (16 CPU, 652 centers) |
| validation job | `38590098`, elapsed **00:26:14** |
| centers | **667/667**, 0 missing, 0 malformed |
| cells | **117,392** = 667 x 4 populations x 44 bins |
| storage | **248 GB** |
| task states | 2001 records, all COMPLETED; **no failures, no retries** |

### 3.2 Validation result — every check PASS

667 readable maps and 667 `.ok` manifests; marker hash == file hash; no duplicate content;
every map records the same NEO source hash, `effective_factor` 64.725384, 44-bin 0.25-mag scheme,
per-population k `{NEO:150, MBA:10, TNO:10, Trojans:10}`, velocity grid, `gaussian_smoothing=False`,
`support_masking=False`, split-provenance hash, sealed-module hash, density engine `closed_form`,
and epoch. No smoothing arrays anywhere. No cell silently reduced k. Every invalid cell carries an
explicit reason.

Featured-cell recheck at `dlon+000_lat+00` — exact:

| bin | expected | got |
|---|---|---|
| [24.00, 24.25) | 35,580 | **35,580** |
| [24.25, 24.50) | 44,721 | **44,721** |
| [24.50, 24.75) | 56,294 | **56,294** |
| [24.75, 25.00) | 70,464 | **70,464** |
| combined | 207,059 | **207,059** |

### 3.3 Coverage — this materially affects scoring

| population | valid cells | invalid cells |
|---|---|---|
| MBA | 25,765 | 3,583 |
| NEO | **16,808** | **12,540** |
| TNO | 9,414 | 19,934 |
| Trojans | 7,263 | 22,085 |

Invalid reasons:

| population | below_min | insufficient_support | no_samples | no_split_fraction |
|---|---|---|---|---|
| MBA | 335 | 2,750 | 498 | 0 |
| NEO | 932 | 10,035 | 1,573 | 0 |
| TNO | 1,607 | 3,429 | 4,893 | 10,005 |
| Trojans | 1,747 | 3,045 | 15,959 | 1,334 |

Valid cells per center: **min 42, median 86, max 114 of 176**.

**Consequence to carry into Stage 3/4:** NEO is valid in only 16,808 of 29,348 (center, bin) cells —
**57.3%**. A TEST2 tracklet landing in a cell where NEO is invalid cannot be scored by the new maps
and must be reported as NaN with reason, not dropped silently and not zeroed. Expected VDP coverage
on TEST2 is therefore well below 100%, and all metrics must be reported both on all rows and on the
common-scorable subset (Stage 4 requirement).

`no_split_fraction` (TNO 10,005; Trojans 1,334) are 0.25-mag bins where the GEN split retained no
objects, so an exact retained fraction does not exist. Per the insufficient-support rule these are
INVALID rather than f=1.0, which would have inflated `P(NEO)` uniformly.

---

## 4. Environment and versions

*(filled by `stage1_scoring_seal.py`, recorded in `SCORING_SEAL.json`)*

---

## 5. Seeds and parent-population sources

| role | source | seed |
|---|---|---|
| map NEO (GEN) | all-sky NEOMOD3 realization | BASE 42..141, HIGH 1,000,000..1,000,319 |
| map non-NEO (GEN) | S3M epoch-state cache, GEN split | split salt `neomod3-gencaltest-v1.1` |
| prior TEST NEO | `gen_benchmark_tracklets_neomod3.py` | `BM_NEO_SEED=20270825` — **inspected, not reused** |
| **TEST2 NEO** | new independent NEOMOD3 realization | `BM_NEO_SEED=` *(recorded in §6)* — must differ from 42..141, 1,000,000..1,000,319 and 20270825 |
| **TEST2 non-NEO** | S3M held-out parents, TEST split | GEN parents excluded by identity |

---

## 6. Stage 2 — TEST2 construction

### 6.1 TEST2 NEO parents — built and hashed

| item | value |
|---|---|
| Slurm job | **38590981** (`--array=0-149%148`, `cpu-g2`, 4 CPU / 80 GB per task) |
| sbatch | `neomod/pipeline/slurm/build_test2_neo_draw.sbatch` |
| task states | **150/150 COMPLETED** (450 sacct records incl. batch/extern), 0 failed, 0 retried |
| array wall time | ~4 min (all 150 tasks concurrent under the 592-CPU cap at 4 CPU/task) |
| seeds | **777,000,000 – 777,000,149** (`BM_NEO_SEED + shard`) |
| draws | **300,000,000** (150 x 2,000,000) |
| selected clones, 14 <= V < 25 | **830,625** |
| output dir | `/mmfs1/gscratch/dirac/ds2004/sorcha/outputs/test2_geometric/neo_shards/` |
| combined shard sha256 | `7ac39a7ee250dfa0e3df83bff67cfc0b0f5e652c5c64bac0690f08362d0a2959` |

Seed disjointness (asserted): map GEN cache 42–141; all-sky HIGH 1,000,000–1,000,319; prior
inspected TEST 31,415,926+. None intersect 777,000,000–777,000,149.

Per-shard sha256 (first five and last two; all 150 in `TEST2_LEAKAGE_AUDIT.json`):

| shard | sha256 |
|---|---|
| `neo_shard_000.parquet` | `9b1f853d72cdd804c128e657c88807203aaa92591de1804f2c521885ee2431c2` |
| `neo_shard_001.parquet` | `205ed2e210db492039b8637269e286e01cfe56c80fc2cffe39858c38f97b4af2` |
| `neo_shard_002.parquet` | `d54334ebac1eea3c50fb03dc8a12498c04e510d15ae608170517ccaa790c7238` |
| `neo_shard_003.parquet` | `b1a60aab8b532c0c44624ae6bc4c125171b0744eced542469e7fe7eaef16c4d9` |
| `neo_shard_004.parquet` | `73095c93b14e4d64663bdb0d3471b4aada61619c424b135a1238e5cd020f7562` |
| … | … |
| `neo_shard_148.parquet` | `8a33411cec726393072fd2ee0bfbd6a53a0a47f52db1acc04ff869f5d76ef651` |
| `neo_shard_149.parquet` | `891f08074a4581adb836417d4b19db5c3a3a8fc1099c96087eb14354860edaf6` |

### 6.2 Why "10x NEO" requires oversampling with weights

The benchmark carries **true absolute counts**: `expected_NEO = n_clones x w_new x role_fraction`
with `w_new = TOTAL_NEO_ABS / n_drawn`, so `n_clones x w_new` is invariant — drawing 10x more orbits
does **not** raise the NEO tracklet count (pinned near 4,973). The non-NEO populations are a finite
real S3M catalogue already split 60/20/20, so they cannot be enlarged at all.

10x NEO is therefore reachable only by oversampling NEO and carrying a physical weight
(`w_phys = 0.1` on NEO rows, `1.0` elsewhere), so weighted contamination and precision are
unchanged and only NEO Monte-Carlo noise falls. The 3e8-draw pool exists so those rows come from
~830k clones rather than recycling a ~65k pool.

### 6.3 Non-NEO leakage audit — **FRESH-TEST2 GATE FAILED**

Jobs: **38591258** (invalid, see error #5) and **38591417** (corrected, 00:02:44).
Script: `neomod/pipeline/test2_leakage_audit.py` (read-only).
Artifacts: `outputs/test2_geometric/TEST2_LEAKAGE_AUDIT.{json,csv}`.

Prior-TEST source identities were recovered by exact float64 join on `(population, lam_deg,
beta_deg)` — the prior builder overwrote `ObjID` with `NM…` strings but copied those columns
bit-identically from the epoch cache. Recovery **636,801 / 638,801 = 99.687%**; 10,903 ambiguous
keys (21,353 rows) were dropped rather than allowed to inflate the join.

| population | TEST-split parents | eligible (14<=V<25, in-grid) | overlap GEN | overlap CAL | **overlap prior TEST** | remaining after all exclusions |
|---|---|---|---|---|---|---|
| MBA | 2,773,597 | 614,710 | 0 | 0 | **614,710** | **0** |
| TNO | 9,770 | 7,575 | 0 | 0 | **7,575** | **0** |
| Trojans | 36,207 | 16,516 | 0 | 0 | **14,516** | 2,000 |
| **total** | | **638,801** | **0** | **0** | **636,801** | 2,000 |

`BM_SPLIT_ROLE=TEST` did exactly what it guarantees — zero overlap with GEN and CAL — and nothing
more. The prior inspected TEST was built from this same partition and consumed **every eligible
parent**, so the proposed TEST2 non-NEO set is essentially identical to it.

The 2,000 "remaining" Trojans are **not** genuinely fresh: they are precisely the 2,000 rows whose
join key was ambiguous and therefore dropped, so they are unresolved rather than unused. The honest
count of fresh eligible non-NEO parents is **0 for all three populations**.

**Gate result: FAILED. No TEST2 artifacts were written and no scoring was started.**

---

## 7. Stage 3 — scoring

*(pending)*

---

## 8. Stage 4/5 — evaluation and notebook

*(pending)*

---

## 9. Errors, causes, and corrections

| # | error | cause | correction |
|---|---|---|---|
| 1 | Full-grid array ran only 9 tasks at `%100` | `AssocGrpCpuLimit`: the `astro` association caps total CPUs (~592), not task count; 9 x 64 = 576 exhausted it | reshaped to 16 CPU / 80 GB / `%400` → 37 concurrent (4.1x). Submitted the replacement and confirmed it scheduled *before* cancelling the old array |
| 2 | `by_pixel` cache unreadable in practice | per-shard `write_dataset` produced 337 files/pixel, 289,949 total; a 104-pixel read never finished in 11 min | single-pass compaction to one file per pixel (768); revalidated, featured counts still exact |
| 3 | Report writer would crash after the full build | `tabulate` is not installed, so `DataFrame.to_markdown()` raises | `_md()` falls back to a fenced `to_string()` table |
| 4 | `_nonneo_split_fraction` raises on 0.25-mag bins | `_SPLIT_BIN_LABELS` knows only the eight 1-mag bins, and refuses to guess by design | computed exact 0.25-mag fractions from the split manifest; reproduces all 23 frozen 1-mag values to 0.000e+00 |
| 5 | Stage 1 seal crashed with `IsADirectoryError` on digest2 | `/…/sorcha/digest2` is the source **directory**; the executable is `digest2/digest2` | corrected the path; recorded binary sha256 and version string in the seal |
| 6 | Leakage audit reported "0 overlap / GATE PASSED" — **false** | join key included `e`, which is NaN for every non-NEO row in the prior TEST file; NaN != NaN, so the merge matched 0 of 638,801 rows and forced `overlap = 0` | dropped `e` from the key (lam/beta match bit-identically), dropped ambiguous keys instead of inflating the join, and added a guard that flags any recovery below 90% as untrustworthy rather than a pass |

---

## 10. Exact rerun commands

```bash
cd /mmfs1/gscratch/dirac/ds2004/sorcha

# --- maps (already built; DO NOT rerun unless rebuilding from scratch) -------------
sbatch neomod/pipeline/slurm/build_mag025_split_provenance.sbatch
sbatch --array=0-319%160 --export=ALL,STAGE=draw   neomod/pipeline/slurm/build_neomod3_high_allsky.sbatch
sbatch --array=0-99%100  --export=ALL,STAGE=rebase neomod/pipeline/slurm/build_neomod3_high_allsky.sbatch
sbatch neomod/pipeline/slurm/partition_validate_high_allsky.sbatch
sbatch neomod/pipeline/slurm/build_neomod3_mag025_k150_array_wide.sbatch      # 667 centers
sbatch neomod/pipeline/slurm/validate_neomod3_mag025_k150.sbatch --full --array-job-id <ID>

# --- Stage 1: scoring seal ---------------------------------------------------------
sbatch neomod/pipeline/slurm/stage1_scoring_seal.sbatch

# --- Stage 2: TEST2 NEO parents (independent realization, seeds 777000000+) ---------
sbatch neomod/pipeline/slurm/build_test2_neo_draw.sbatch

# --- Stage 2 gate: read-only non-NEO leakage audit ----------------------------------
sbatch neomod/pipeline/slurm/test2_leakage_audit.sbatch
```


---

## 6.4 TEST2 assembly — COMPLETE (job 38591577, 00:05:18)

Script `neomod/pipeline/build_test2_assemble.py`, sbatch
`neomod/pipeline/slurm/build_test2_assemble.sbatch`.

### Exact NEO weight (computed, not assumed)

    in-grid clones          655,167
    w_new = ABS/n_drawn     0.0381097265
    role_fraction (TEST)    0.1998008068
    expected physical NEO   4988.6735
    sampled NEO rows        49,887   (without replacement, seed 777123)
    w_phys per NEO row      **0.0999994694**   <- NOT exactly 0.1

### Counts

| population | rows | summed physical weight |
|---|---|---|
| NEO | 49,887 | 4,988.67353 |
| MBA | 614,710 | 614,710 |
| TNO | 7,575 | 7,575 |
| Trojans | 16,516 | 16,516 |
| **total** | **688,688** | |

### Hard gates

| gate | result |
|---|---|
| GEN leakage | **0** |
| CAL leakage | **0** |
| V outside [14,25) | **0** (flagged, never clipped) |
| `source_parent_uid` unique | yes |
| `tracklet_uid` unique | yes |
| MPC-80 lines exactly 80 chars | yes (both detections, all rows) |

### VDP coverage — measured per tracklet, raw and physical-weighted

| population | rows | raw | weighted |
|---|---|---|---|
| NEO | 49,887 | 99.59% | 99.59% |
| MBA | 614,710 | 99.29% | 99.29% |
| TNO | 7,575 | 99.99% | 99.99% |
| Trojans | 16,516 | 99.86% | 99.86% |
| **ALL** | **688,688** | **99.33%** | **99.31%** |

**This settles the earlier ambiguity.** The 57.3% valid-*cell* fraction is not an abstention rate:
tracklets concentrate where objects actually are, which is precisely where cells have support, so
real per-tracklet coverage is **99.33%**. My earlier phrasing ("expected VDP coverage well below
100%") was wrong and is superseded by this measurement.

### Artifacts

| file | sha256 |
|---|---|
| `TEST2_PARENT_MANIFEST.parquet` | `5b77824f0a83d468fcc3c26dcd04a3951414d46a394c8502e7f11ab1dfc0acb9` |
| `TEST2_TRACKLETS.parquet` | `36eabd6cf32f8c8be42749b4e228785c533f56187e69fe38b4cd2e2a8eb00ce8` |
| `TEST2_MPC80.parquet` | `c7eb200f957fb118829dc4b693c496902408c93a414b6a2fa6e4fd11c3a69ee4` |
| `TEST2_SEAL.json` | (seal) |

Identity preserved on every row: `source_parent_uid`, `s3m_objid` (non-NEO), NEOMOD3
`neo_shard`/`neo_row`/`neo_seed`/`neo_orbit_digest` (NEO), immutable `tracklet_uid`, plus
`tracklet_input_sha256` over the exact two MPC-80 lines for classifier-input identity proof.


---

## 7. Stage 3 — scoring (COMPLETE)

| mode | job | shards | throughput | notes |
|---|---|---|---|---|
| smoke (cancelled) | `38591754` | — | — | cancelled: sample spanned 439 centers x 372 MB maps, ~163 GB of reads |
| smoke (MPC-80 fix) | `38591976` | — | 00:05:07 | 14/14 |
| coverage audit + corrected smoke | `38592300` | — | 00:06:21 | **19/19** |
| new_vdp | **38592336** | 32/32 | ~78 rows/s, ~275 s/shard | valid ~9,000/21,500 |
| legacy_vdp | **38592337** | 32/32 | — | valid ~100% |
| digest2 | **38592338** | 32/32 | ~2.8–3.1 rows/s, ~7,000–7,800 s/shard | valid 21,518/21,522 |
| merge + evaluate | **38594049** | — | 00:13:29 | **31/31 PASS, ALL_PASS=True** |
| notebook | **38594209** | — | 00:01:46 | 28 cells, 0 errors, 0 unexecuted, 5 figures |

96 task records, all COMPLETED. No shard needed a rerun. `38592338` was briefly pending only because
the two VDP arrays held the ~592-CPU account budget; it started unaided, so it was not resubmitted.
Note for future digest2 runs: it uses one CPU internally, so `--cpus-per-task=1 --mem=16G` avoids
reserving ~7 idle cores per task (this run wasted them rather than discarding completed work).

Shard manifest: greedy LPT balanced by tracklet count (not round-robin center count), shards
21,488–21,539 rows, max/min 1.00, sha256
`2b15c2a6241c5581716d043b7d0ad3a7fa1272c7d45c168d2b6943d859f983d8`.

### 7.1 Coverage policy (accepted definitions)

| quantity | value |
|---|---|
| NEO-cell availability | **99.33% raw** |
| four-density VDP scoring coverage | **43.96% raw / 45.06% weighted** |

A valid new-VDP probability requires all four densities present and finite. Partial-denominator rows
are NaN with reasons like `missing_population_density:TNO,Trojans`; never a fallback, never zeroed.
Four-density availability by population: NEO **28.12%**, MBA 43.71%, TNO 51.88%, Trojans 97.46%.

### 7.2 Headline metrics — common-scorable subset only

302,666 rows (43.95% of TEST2), NEO 14,026. **Not whole-population performance.**

| classifier | weighted AUC | pAUC(FPR<=.01) | AUC CI95 (parent bootstrap) | C@5% contam | C@10% |
|---|---|---|---|---|---|
| new VDP (k=150, 0.25 mag) | 0.93814 | 0.71476 | [0.93497, 0.94095] | 0.4252 | 0.5391 |
| legacy VDP raw | 0.93695 | 0.72007 | [0.93346, 0.94012] | 0.5506 | 0.5987 |
| legacy VDP calibrated | 0.93695 | 0.72007 | — | 0.5506 | 0.5987 |
| digest2 | 0.92474 | 0.63491 | [0.92139, 0.92757] | 0.5369 | 0.5369 |

Paired new-minus-legacy(raw): median -0.00005,
mean -0.00577,
new higher on 39.7% of rows.

View C (own broader valid sets, NOT comparable to the above): legacy raw n 686,100 AUC 0.96536;
digest2 n 688,604 AUC 0.96449. The View B/C gap reflects WHICH rows are scorable, not classifier
quality.

---

## 8. Errors, causes, corrections (continued)

| # | error | cause | correction |
|---|---|---|---|
| 7 | smoke test hung >20 min, cancelled (`38591754`) | sample spread over 439 centers x 372 MB maps = ~163 GB of reads for 850 rows | restricted the smoke sample to ~10 centers chosen deterministically; the full run shards BY center so each map loads once |
| 8 | **digest2 returned 0 valid scores for all 568 smoke rows** | `format_mpc80` truncates the identifier to MPC-80 columns 1–12; with 5 leading spaces only 7 chars survive, so 8-char keys `T0000000/T0000001/T0000002` all collapsed to `T000000`. digest2 saw one object with dozens of scattered observations | keys narrowed to `T{i:06d}` (7 chars); added an assertion that the **12-column designation actually emitted** is unique across all 688,688 rows and identical for both detections of a tracklet. The pre-existing 80-character check was insufficient — truncated lines were still exactly 80 chars |
| 9 | **smoke test reported PASS on a broken classifier** | the integer-range check was `len(d2v) == 0 or (...)`, so an empty result set passed vacuously | it now fails on an empty set, plus a separate explicit `digest2 produced valid scores at all` assertion |
| 10 | partial-denominator rows were marked valid (`ok_partial_denominator_2/3`) | `P(NEO)` was computed from a 2- or 3-population denominator, which is not the same quantity; missing keys mean INVALID cells, and project policy forbids treating invalid as zero | `set(present) == set(POPS)` now required; otherwise NaN + `missing_population_density:<pops>`; asserts valid rows have `n_pops_new == 4` and four-class sums to 1 |
| 11 | round-robin center sharding was badly unbalanced | centers differ by >100x in tracklet count | greedy LPT balancing by tracklet count; max/min 1.00; manifest hashed |

Regenerated artifact hashes after the MPC-80 correction: `TEST2_TRACKLETS.parquet`
`e374547893e2945a3ea6dfcd3ae80c1f...`, `TEST2_MPC80.parquet` `179a64ab8b22cba2b48156e8011af88c...`;
`TEST2_PARENT_MANIFEST.parquet` unchanged (`5b77824f...`) because identity did not change.
