# Sorcha Full Pipeline — S3M Linking Experiment

**Created:** 2026-07-06. **All three cases completed:** 2026-07-07.
**Purpose:** Measure how the two Rubin SSP linking constraints — the multi-night
tracklet requirement and the inter-detection motion threshold — suppress the
recovered population (especially TNOs), by re-running the full Sorcha → VDP →
digest2 pipeline three times, changing **only** those two parameters.

This document is the single source of truth for the run: the three cases, every
file involved, the resource specs, the exact submit order, and cleanup. See §9 for
final results. A related, independently-run effort — rebuilding the benchmark with
population-proportional caps ("Option A", `benchmark_comparison_s3m_v2.parquet`) and
then also epoch-matching it to the busiest Sorcha night found via this pipeline's
case1/2/3 output (`benchmark_comparison_s3m_v3.parquet`, MJD 61642 = 2027-08-25) — is
documented in `WAGG_SORCHA_HYAK_CONTEXT.md` ("benchmark population-cap discovery" and
"Benchmark v3"), not here; it shares no jobs or scripts with this pipeline, but v3's
epoch choice was derived directly from this pipeline's per-night tracklet counts
(§9 case1/2/3 outputs).

---

## 1. Why this experiment exists

The original pure-S3M Sorcha run (`sorcha_comparison_s3m.parquet`) was produced
with `Rubin_full_footprint_detections.ini`, which has **no `[LINKINGFILTER]`
section**. Sorcha therefore emitted every detection, and `sorcha_postprocess.py`
built tracklets using only a 3–90 min time filter. **No SSP linking was ever
applied.** That run is effectively "1 tracklet, no motion threshold."

To attribute the TNO loss to specific linking constraints we now let Sorcha apply
its native SSP linking, varying two parameters:

| Parameter | Meaning |
|-----------|---------|
| `SSP_number_tracklets` | nightly tracklets required (within a 15-day window) to call an object detected |
| `SSP_separation_threshold` | minimum inter-detection motion (arcsec) for a pair to count as a tracklet |

---

## 2. The three cases

Every other pipeline input and parameter is **identical** across cases. Only the
two values below differ (verified: the three config files are byte-identical apart
from these two lines).

| Case | `SSP_number_tracklets` | `SSP_separation_threshold` | Isolates |
|------|:---:|:---:|----------|
| **case1** | 1 | 0.5 arcsec | motion threshold only (no multi-night requirement) |
| **case2** | 3 | 0.5 arcsec | full LSST linking (real SSP baseline) |
| **case3** | 1 | 0.001 arcsec | neither constraint (≈ the original run) |

Interpretation of the comparisons:
- **case2 → case1**: cost of the 3-night linking requirement.
- **case1 → case3**: cost of the 0.5 arcsec motion threshold.
- **case2** is the scientifically correct "what Rubin SSP actually recovers" baseline.

> Note: `SSP_separation_threshold` must be `> 0` — Sorcha rejects `0.0`. `0.001`
> arcsec is used as an effective zero.

Config files (`neomod/pipeline/config/`):
- `s3m_case1_tk1_sep05.ini`
- `s3m_case2_tk3_sep05.ini`
- `s3m_case3_tk1_sep0.ini`

All three are the original `Rubin_full_footprint_detections.ini` (output `hdf5`
/ `all`, so postprocess can read the rate columns) **plus** a `[LINKINGFILTER]`
section. `output_format = hdf5` and `output_columns = all` are **required** —
do not switch to the stock `Rubin_full_footprint.ini`, which outputs
`sqlite3`/`basic` and breaks `sorcha_postprocess.py`.

---

## 3. Pipeline stages, scripts, and specs

All stage scripts live in `neomod/pipeline/slurm/s3m_linking/` and share
`_case_env.sh`, which maps `$CASE` → config + all paths. The case is chosen at
submit time with `--export=ALL,CASE=caseN`; the scripts themselves never change.

Resources were tuned for max concurrency on `ckpt-all` (~10.5k idle CPUs); the
per-task work matches the original validated run. Wall-time estimates come from
the original run's sacct history.

| # | Script | Resources | Array | Est. wall | Reads | Writes |
|---|--------|-----------|-------|-----------|-------|--------|
| 1 | `1_production.sh` | 16 cpu / 96G / 4h | `0-898%450` | ~1–1.5 h | s3m orbits+phys, pointings db, AR cache | `production/inst*_part*.h5` |
| 2 | `2_postprocess.sh` | 8 cpu / 24G / 2h (parallel 8) | `0-899%450` | ~10–15 min | `production/*.h5`, `prob_maps_grid` | `tracklets/tracklets_*.parquet` |
| 3 | `3_vdp.sh` | 2 cpu / 32G / 6h | `0-112%113` | ~10 min | `tracklets/`, `prob_maps_grid_s3m` | `phase2/vdp_shards/` |
| 4 | `4_sample.sh` | 1 cpu / 16G / 1h | single | ~10 min | `phase2/vdp_shards/` | `phase2/subsample.parquet` |
| 5 | `5_digest2.sh` | 4 cpu / 16G / 1.5h (parallel 4×1250) | `0-N%130` | ~15–20 min | `phase2/subsample.parquet`, `digest2/` | `phase2/digest2_shards/` |
| 6 | `6_combine.sh` | 1 cpu / 16G / 1h | single | ~10 min | `phase2/digest2_shards/` | `sorcha_comparison_<case>.parquet` |

**Total per case ≈ 2–2.5 h wall** (was 4–6 h). Partition is `ckpt-all` (falls
back to `ckpt` if the account isn't permitted there). `ckpt` is preemptible —
stages 1/2 are idempotent (`h5_is_good` / `--skip-existing`), so re-`sbatch` any
stage to refill preempted/timed-out tasks.

**Two different map directories (by design, preserved from the original):**
- `prob_maps_grid` — used in stage 2 for **footprint / cell assignment** (geometry only).
- `prob_maps_grid_s3m` — used in stage 3 for **VDP scoring** (pure-S3M densities).

**Per-case output tree** (`outputs/s3m_linking/<case>/`):
```
production/                      inst*_part*.h5         (~197 G — DELETE after stage 2)
tracklets/                       tracklets_*.parquet    (~17 G)
phase2/vdp_shards/               vdp_*.parquet
phase2/subsample.parquet
phase2/digest2_shards/           d2_*.parquet
sorcha_comparison_<case>.parquet                        (final product)
work/s3m_linking/<case>/         Sorcha scratch         (DELETE after stage 1)
```

---

## 4. Submit order

Run stages **in order per case**; each waits on the previous. You can run the
three cases concurrently if disk allows (see §6), or sequentially to cap disk at
one production run at a time.

```bash
cd /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/slurm/s3m_linking
C=case1      # then repeat for case2, case3

sbatch --export=ALL,CASE=$C 1_production.sh
# wait for completion, then:
sbatch --export=ALL,CASE=$C 2_postprocess.sh
# wait, then:
sbatch --export=ALL,CASE=$C 3_vdp.sh
# wait, then:
sbatch --export=ALL,CASE=$C 4_sample.sh          # prints subsample row count
# read the printed row count R, set array = ceil(R/5000)-1:
sbatch --export=ALL,CASE=$C --array=0-129%130 5_digest2.sh
# wait, then:
sbatch --export=ALL,CASE=$C 6_combine.sh
```

**Array sizes that vary per case:**
- Stage 3 (VDP): `NSHARDS = ceil(#tracklet_parquets / 128)`. Default `0-112`
  fits case1/case3; case2 has fewer files, extra shards no-op harmlessly. Exact:
  ```bash
  n=$(ls outputs/s3m_linking/$C/tracklets/*.parquet | wc -l); echo $(( (n+127)/128 - 1 ))
  ```
- Stage 5 (digest2): `ceil(subsample_rows / 5000)`. Stage 4 prints the row count.
  case1/case3 ≈ 648k rows → `0-129`. case2 is much smaller → set a smaller array.

---

## 5. Verifying a completed case

```bash
python - <<'PY'
import pandas as pd
for c in ["case1","case2","case3"]:
    try:
        df = pd.read_parquet(f"outputs/s3m_linking/{c}/sorcha_comparison_{c}.parquet",
                             columns=["population","P_NEO_vdp","P_NEO_d2"])
        print(c, len(df), df["population"].value_counts().to_dict())
    except FileNotFoundError:
        print(c, "not done")
PY
```
The headline number for the experiment is the **TNO count per case** (and per
magnitude bin) — that quantifies how many TNOs each linking constraint removes.

---

## 6. Disk and cleanup

The dirac allocation was at **97% (≈8.5 T free)** when this was set up. Each case
costs roughly:

| Item | Size | Lifetime |
|------|------|----------|
| `production/*.h5` | ~197 G | **delete after stage 2** |
| `work/` scratch | ~4 G | **delete after stage 1** |
| `tracklets/` | ~17 G | keep (cheap) |
| `phase2/` | ~14 G | keep |
| `sorcha_comparison_<case>.parquet` | ~1 G | keep (the product) |

**Cleanup after stage 2 verifies for a case** (reclaims ~200 G/case):
```bash
C=case1
# confirm tracklets exist first:
ls outputs/s3m_linking/$C/tracklets/*.parquet | wc -l
rm -rf outputs/s3m_linking/$C/production work/s3m_linking/$C
```

**Disk strategy:** if running all three concurrently you need ~3×200 G ≈ 600 G of
production h5 live at once (fits in 8.5 T). To be safe, run **sequentially** and
delete each case's `production/` before starting the next — peak footprint then
stays near one production run (~215 G).

---

## 7. Compute cost (read before submitting)

Stage 1 (production) is the single most expensive job in the whole project —
899 tasks × up to 4 h × 16 cores — and this runs it **three times**. Everything
downstream is cheap by comparison.

**Optimisation available (not yet wired):** Sorcha's cost is the ASSIST+REBOUND
ephemeris generation, *not* the linking. The linking is cheap post-ephemeris. If
stage 1 were run once with the ephemerides cached/exported and reused via
`ephemerides_type = external`, the three linking variants would cost ~1× the
ephemeris generation + 3× cheap linking instead of 3× the full run. This requires
verifying Sorcha's ephemeris-reuse path against this cache; deferred unless the
3× production cost is a problem.

---

## 8. Relationship to existing files

- `sorcha_comparison_s3m.parquet` (existing, 648,908 rows) = "no linking" run.
  It is roughly what **case3** reproduces; comparing case3 to it is also a sanity
  check that the native-linking path with `sep≈0, tk=1` matches the old
  postprocess-only path.
- The linking-analysis writeup that motivated this is in
  `WAGG_SORCHA_HYAK_CONTEXT.md` → "Pure-S3M benchmark vs Sorcha: structural gap
  analysis (2026-07-06)" and "TNO contamination…".

---

## 9. Results (all three cases complete, 2026-07-07)

| Case | File | Rows | NEO | MBA | TNO | Trojan | other |
|---|---|---:|---:|---:|---:|---:|---:|
| case1 (tk=1, sep=0.5″) | `outputs/s3m_linking/case1/sorcha_comparison_case1.parquet` | 648,769 | 148,773 | 463,419 | 4,965 | 12,723 | 18,889 |
| case2 (tk=3, sep=0.5″ — real LSST) | `outputs/s3m_linking/case2/sorcha_comparison_case2.parquet` | 624,783 | 124,748 | 463,831 | 4,742 | 12,788 | 18,674 |
| case3 (tk=1, sep≈0) | `outputs/s3m_linking/case3/sorcha_comparison_case3.parquet` | 648,828 | 148,828 | 463,698 | 4,930 | 12,743 | 18,629 |

**Headline finding:** TNO barely moves across all three cases (4,965 → 4,742 → 4,930,
within ~5%) — the linking-filter parameters are *not* the dominant driver of TNO
suppression, contrary to the original hypothesis motivating this experiment (§1).
**NEO drops ~16%** under case2's real linking (124,748 vs ~148,800 in case1/case3) —
the 15-day/3-tracklet window is evidently a real constraint for NEOs specifically
(fast, brief-visibility-window objects), not TNOs. Full ROC/F1 comparison across
cases not yet run as of this writing — all data is local, no further SLURM needed.

**Pathological-orbit losses were consistent and small across all three runs**,
concentrated in the same catalog region every time (instances ~849–898, the
TNO/Trojan tail of the combined S3M input):
- case1: 2/899 instances failed (2 missing parts of 14,381)
- case2: ~40/14,381 parts missing (24 instances initially failed; retried at full
  concurrency then at low concurrency to rule out contention — 16/24 remained
  failed both times, confirming genuine per-object integrator issues, not preemption)
- case3: 2/899 instances failed — **the same two instances (882, 884)** seen in
  case1/case2, confirming this is catalog-position-specific, not run-specific.

None of these were chased further (consistent with documented precedent for
S3M/hybrid production runs — see `WAGG_SORCHA_HYAK_CONTEXT.md` "Two Recovered
Problem Files" and `inst00820_part003`). Total loss across all three cases is
<0.3% of the 14.38M-object catalog.
