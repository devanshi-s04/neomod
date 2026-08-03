# New cloning scheme: NEOMOD3-only NEO clones, no MBA/TNO/Trojan cloning

**Created:** 2026-07-29.
**Status:** comparison evidence COMPLETE; new module WRITTEN and import-verified; **no maps built yet**.
**Decision (advisors, 2026-07-29):** stop cloning MBA / TNO / Trojans entirely. Clone only NEOs, and
draw those NEOs from **NEOMOD3** instead of S3M. This fixes the orbital-distribution problem the
GMM cloner was working around.
**Related:** `cloning_gmm_neo.md` (the GMM cloner this replaces), `ARNOR_NOTE_GMM_BUG.md`
(normalisation), `fixing_integrator.md` §13.6 (the hybrid-loader trap), `mag245_tests.md`
(the VDP-zero / digest2-positive cohort that motivated looking at NEO coverage).

> **Headline:** S3M is **truncated at H = 25** by construction; NEOMOD3 covers 15 < H < 28, and ~**97%**
> of what NEOMOD3 describes is a faint population S3M has **zero** objects for. S3M also shows visible
> generation artifacts in 2-D (discrete inclination banding). NEOMOD3 is debiased and smooth.
> **But** the two models genuinely **disagree about NEO orbital shape**, and that disagreement
> **survives an H-matched cut** — so switching the clone source changes the intrinsic population the
> maps represent. That is an intended consequence to report, not a bug to fix.

---

## §1. Why change

The production NEO cloner (`_clone_neo_gmm` in `velocity_density_pipeline_gmm.py`) does:

1. take the **visible S3M NEOs** in the sky patch,
2. *augment* them with NEOMOD3 samples (`_neomod3_sampler.sample_neomod3_orbits`, line ~1440),
3. fit a **GaussianMixture** (80 components) to the combined set,
4. sample clones from that GMM,
5. attach `H` drawn from the **S3M empirical** distribution — explicitly **not** NEOMOD3's H
   (`h_source = h_source_visible`, line ~1469).

So even where NEOMOD3 was already in the loop, its distribution got **re-smoothed through a GMM** and
its `H` information was **thrown away**. That is the "GMM smooths the distribution" problem recorded in
`cloning_gmm_neo.md` (the source `q` spike near 1.1–1.3 is washed out).

Separately, MBA/TNO/Trojan were cloned with conditional K|M purely to get enough points for the kNN
density estimator — which is not a physics argument, and K|M only repeats parent orbits anyway.

**New scheme:** NEO clones come straight from NEOMOD3 (no GMM, no S3M anchor); the other three
populations are not cloned at all.

---

## §2. NEOMOD3: the files, and which one we actually need

Upstream: <https://www2.boulder.swri.edu/~davidn/NEOMOD_Simulator/> (David Nesvorny, SWRI).
Paper: NEOMOD3, arXiv:2404.18805.

| file | what it is | do we need it? |
|---|---|---|
| `input_neomod3.dat` | **the binned NEOMOD3 model** (62 MB) — the datacube | **YES — already on disk** |
| `neomod3_simulator.f` | reference Fortran sampler | validation reference only |
| `neomod3.par` | its parameter file (seed, N, diameter range) | reference only |
| `give_me_prob.f` | source probabilities for a given (a,e,i,H) | not needed |
| `neomod2*`, `magmod*`, `diamod*` | older / component models | not needed |

**We already have the datacube locally** — no download needed:

```
neomod/NEOMOD3/input_neomod3.dat      62 MB, fetched 2026-06-01
neomod/src/neomod3_sampler.py         working sampler (load_neomod3_array, sample_neomod3_orbits)
neomod/src/NEOMOD3.py                 older/alternate parser -> array4D
```

> **Cert trap:** the SWRI host serves a self-signed cert. `curl -k` (or `ssl.CERT_NONE`, which
> `neomod3_sampler.load_neomod3_array()` already does in its download fallback) is required.

### §2.1 Datacube structure

`input_neomod3.dat` is **not a catalog of objects** — it is a binned 4-D probability table. Each row is
`iH ia ie ii  weight` followed by the fractional contribution of 12 source regions.

| axis | bins | min | max | width |
|---|---:|---:|---:|---:|
| H | 52 | 15.0 | 28.0 mag | 0.25 |
| a | 42 | 0.0 | 4.2 AU | 0.10 |
| e | 25 | 0.0 | 1.0 | 0.04 |
| i | 22 | 0.0 | 88.0° | 4.00 |

- total cells 52×42×25×22 = **1,201,200**; **566,835 populated (47%)**.
- 12 sources: `nu6, 3:1, 5:2, 7:3, 8:3, 9:4, 11:5, 2:1, inner, Hungaria, Phocaea, JFC`.
- header also carries the **size distribution**: reference diameter `Dref = 1.0 km` with
  `871.05` NEOs above it, and 8 log-D spline segments with slopes — this is NEOMOD3's own
  **absolute** normalisation (relevant to §6).
- **calibrated for 15 < H < 28**; outside that range the model is extrapolation (stated in the
  simulator header).

### §2.2 How the reference sampler works vs ours

| | reference `neomod3_simulator.f` | our `neomod3_sampler.py` |
|---|---|---|
| size/H | draws log-D from the spline size distribution, then albedo (size-dependent Wright-like model), then `H = -5 log10(D √pV / 1329)` | samples the H **bin** directly from the table |
| orbit | draws `(a,e,i)` **continuously**, then **trilinear interpolation** of the model + **rejection** | picks a cell by weight, then **uniform dither** inside the cell |
| result | piecewise-**linear** density | piecewise-**constant** (histogram-like) density |
| angles | n/a (not emitted) | `node, argperi, M` uniform; `t_p` reconstructed from `M` at the map epoch |

Both are valid samples of the same binned model. Ours is slightly blockier at bin edges; since clones
go through a kNN density estimator downstream, this is a **fidelity nuance, not a correctness bug**.
Our sampler's mechanics are independently validated (§3: mean-anomaly KS = 0.002, p = 0.66).

---

## §3. S3M vs NEOMOD3 — the comparison

**Setup.** Pure S3M NEOs (`VDP_LOADER=s3m`, `S0.s3m`, `t_0` = 2008): **268,511** objects.
NEOMOD3: 500,000 draws → **485,084** valid NEO orbits (`q < 1.3`). Epoch `2027-08-25T00:00:00`
(the Stage-0 n-body cache epoch, *not* the old notebook's stale 2026-05-09).

### §3.1 Full sample

| element | S3M median | NEOMOD3 median | KS |
|---|---:|---:|---:|
| a | 1.976 | 1.942 | 0.0734 |
| q | 0.791 | 0.953 | 0.2130 |
| e | 0.592 | 0.504 | 0.2227 |
| i | 18.89° | 13.29° | 0.1766 |
| **H** | 24.14 | 27.37 | **0.9694** |
| M (control) | 180.15 | 180.00 | **0.0018** (p = 0.66) |

The `M` panel is the **control**: uniform-in-phase on both sides confirms the sampler's angle/`t_p`
reconstruction is right. The `H` panel is the headline: **S3M stops dead at H = 25**.

### §3.2 H < 25 matched cut — the differences are REAL

The full-sample KS on q/e/i mixes "orbital-shape disagreement" with "different H regimes". Cutting
**both** samples at H < 25 isolates the shape difference:

| element | KS full | KS H<25 | change |
|---|---:|---:|---:|
| a | 0.0734 | 0.0743 | +0.0009 |
| q | 0.2130 | 0.2046 | −0.0084 |
| e | 0.2227 | 0.1924 | −0.0302 |
| i | 0.1766 | 0.1531 | −0.0235 |

Counts in the cut: S3M **268,416 (100.0%)** — S3M is entirely below H = 25 — vs NEOMOD3 **14,849
(3.1%)**.

**The differences survive.** At matched H, NEOMOD3 has:
- **fewer very-low-perihelion objects** (q < 0.4),
- a **lower eccentricity peak** (~0.52 vs ~0.59) and thinner high-e tail,
- **fewer high-inclination objects** (S3M extends far more heavily to i ≈ 60–80°).

So S3M and NEOMOD3 **genuinely disagree about the NEO orbital distribution**, independent of magnitude
range. Switching the clone source **changes the intrinsic NEO population the maps are built from**.

### §3.3 H coverage — what NEOMOD3 adds

| | H range | n |
|---|---|---:|
| S3M | [10.31, **25.00**] | 268,511 |
| NEOMOD3 | [15.25, **28.00**] | 485,084 |

**470,235 NEOMOD3 samples (96.9%) have H ≥ 25 — a population S3M contains ZERO objects for.**
This is the strongest argument for the switch: any cloner anchored on S3M (GMM-smoothed or not) is
*structurally incapable* of representing the faint NEO population — precisely where detection is
hardest, and where the VDP-zero / digest2-positive cohort lives (`mag245_tests.md`, and the 28
out-of-±2 NEOs from the velocity-grid work, all mag 24–25).

### §3.4 Two-dimensional structure — an S3M artifact

In the `a–i` plane **S3M shows strong horizontal striping** (discrete inclination bands) and blocky
`a–e` structure. These are **synthetic-catalog generation artifacts, not physics**. NEOMOD3 is smooth.
Overlaid contours show NEOMOD3 is the **more concentrated** distribution; S3M has broader `e` and `i`
tails.

---

## §4. The new module

`neomod/src/velocity_density_pipeline_neomod_clone_only.py` — a copy of
`velocity_density_pipeline_gmm.py` (2,839 lines) with **exactly two behavioural changes**.

### §4.1 NEO — NEOMOD3 sampled directly

New function `_sample_neo_neomod3_visible(...)`:

- loads the NEOMOD3 4-D table, optionally **zeroes H bins outside `NEOMOD3_H_RANGE` before sampling**,
- draws in batches, projects each batch through the **same** `build_visible_subset_dataframe` sky cut
  every other population gets, accumulates until `NEOMOD3_TARGET_VISIBLE` visible orbits (or
  `NEOMOD3_MAX_DRAWS`),
- returns the visible clone frame + a diagnostic dict (`n_drawn`, `n_visible`, `visible_fraction`, …).

`H` now comes from **NEOMOD3 jointly with (a, e, i)** — the GMM path deliberately substituted S3M's
empirical H; with no GMM in the path there is no reason to, and the joint H–orbit correlation is part
of what makes the model debiased.

**There is deliberately NO K|M fallback for NEO.** A fallback would silently reintroduce S3M-anchored
clones and quietly undo the point of the module. If NEOMOD3 is unusable we want it to fail loudly.

Tunables (env-overridable):

| constant | default | meaning |
|---|---:|---|
| `NEOMOD3_TARGET_VISIBLE` | 200,000 | visible NEOMOD3 orbits wanted per (center, mag-bin) |
| `NEOMOD3_MAX_DRAWS` | 20,000,000 | cap on raw draws to get there |
| `NEOMOD3_DRAW_BATCH` | 2,000,000 | draw batch size |
| `NEOMOD3_H_RANGE` | (15.0, 28.0) | H bins to sample from — **see §6.2** |

### §4.2 MBA / TNO / Trojans — no cloning

`clone_visible_df = df_cloner_input` — the **real visible source objects are the density sample**.
`clone_population_conditional_K_from_M_with_skycut` is retained in the file for reference but is
**never called**. Likewise `_clone_neo_gmm` is retained, definition-only.

`DEFAULT_POPULATION_SETTINGS` clone factors are now **1 for every population** (were NEO 80, TNO 10,
Trojans 5, MBA 5-via-CLI).

### §4.3 effective_factor after the change

```
NEO                : effective_factor = n_visible_NEOMOD3 / n_source_visible_S3M_NEO
MBA / TNO / Trojans: n_visible_clones_gmm is None -> falls through to f (= 1) -> density unchanged
```

For NEO this means **shape from NEOMOD3, amplitude expressed "per real visible S3M NEO in-patch"** —
i.e. NEO stays on the same footing the GMM pipeline used, so this is a controlled change. Whether that
is the *right* footing is §6.1.

### §4.4 Bug caught during the port

The `effective_factor` block still referenced **`gmm_success`**, a variable the rewrite deleted — a
`NameError` waiting to fire on the first NEO map. Fixed; verified no other dangling references
(`gmm_success`, `gmm_diag`, `gmm_clone_df` all gone).

### §4.5 Verification done

```
py_compile                        OK
import (PYTHONPATH=neomod/src)    OK
clone factors                     {'MBA': 1, 'NEO': 1, 'TNO': 1, 'Trojans': 1}
_NEOMOD3_AVAILABLE                True
NEOMOD3_TARGET_VISIBLE / H range  200000 / (15.0, 28.0)
_sample_neo_neomod3_visible       present & callable
grid defaults                     (-2.0, 2.0), 0.01   (unchanged)
K|M cloner called?                no
_clone_neo_gmm called?            no
```

**No maps have been generated with this module yet.**

---

### §4.6 The NEOMOD3 projection cache (added 2026-07-30)

`neomod/pipeline/neomod3_projection_cache.py` + `slurm/neomod3_projection_cache.sbatch`.

Because NEOMOD3 is a **global** model but each map covers a 30° patch (~6.7% of sky), per-center
sampling is prohibitively wasteful (§6.5). The projection is **epoch-dependent, not center-dependent**,
so it is done once:

```
sample NEOMOD3 (full H 15-28)  ->  project ALL-SKY at the map epoch  ->  parquet cache
                                        |
                    each of the 667 centers slices the SAME cache by its own sky cut
```

This plugs into **existing, already-validated code**: `build_visible_subset_dataframe` detects
`_CACHE_POSITION_COLS = (ra_deg, dec_deg, vlam, vbeta)` and skips straight to the center-dependent sky
cut — the same §12.2 Stage-0 path the n-body map regen uses. No new slice function, no second sky-cut
convention to keep in sync.

```bash
sbatch --export=ALL,N_TOTAL=100000000,NSHARDS=100 --array=0-99%100 \
       neomod/pipeline/slurm/neomod3_projection_cache.sbatch
conda_prep/bin/python neomod/pipeline/neomod3_projection_cache.py merge
conda_prep/bin/python neomod/pipeline/neomod3_projection_cache.py stats
```

Cache columns: `a e i node argperi t_p q H | ra_deg dec_deg dra_deg_day ddec_deg_day |
lam_deg beta_deg vlam vbeta | mag_app M_obs_deg` (float32). Idempotent per shard (skip-existing), so a
preempted ckpt task only needs a re-`sbatch`.

**Sizing.** Measured yield per 30° center vs what the current production maps actually contain
(`support_count__NEO__*` summed, map `dlon+000_lat+00`):

| mag bin | production GMM clones | NEOMOD3 @40M | NEOMOD3 @100M |
|---|---:|---:|---:|
| 14_16 | **0** | 0 | 0 |
| 16_18 | 105 | 17 | 42 |
| 18_20 | 189 | 155 | 388 |
| mag20 | 238 | 160 | 400 |
| mag21 | 473 | 381 | 953 |
| mag22 | 1,114 | 762 | 1,906 |
| mag23 | 2,438 | 1,917 | 4,793 |
| mag24+ | 4,815 | 4,531 | 11,327 |
| **total** | **9,372** | 7,923 | **19,808** |
| **clones with \|v\| > 2** | **0** | ~1,287 | ~3,200 |

`14_16` is empty in production too, so a zero there is not a regression. **40M** roughly matches
production; **100M** (adopted) gives ~2x production density plus real high-|v| coverage — the whole
point of the switch. Cost: ~19 s per 250k orbits per shard, so 100M over 100 shards is ~3 min wall and
~8 GB on disk.

#### BUILT + ACCEPTANCE-TESTED (2026-07-30)

`outputs/neomod3_projection_cache/neomod3_projection_20270825T000000.parquet` —
**97,056,755 rows, 7.88 GB**, 100/100 shards, zero failures, H 15.00–28.00.

Acceptance test: sliced at the antisun center (`dlon+000_lat+00`, antisun lon 331.203°, 30° cut) via
`build_visible_subset_dataframe(..., scorer=None)`. **`scorer=None` succeeding is itself the proof the
§12.2 cache branch engaged** — the propagate path would have dereferenced the scorer.

- cache load 9 s; slice 97M → **5,357,573** clones in **18 s**.
- **in map mag range: 45,776** (vs production **9,372**, ~4.9x).
- **in-map with |v| > 2: 4,910** (production: **0**).

| mag bin | NEOMOD3 cache | production | **max \|v\| NEOMOD3** | production max \|v\| |
|---|---:|---:|---:|---:|
| 14_16 | 7 | 0 | 4.43 | 0.00 |
| 16_18 | 58 | 105 | 39.73 | 1.79 |
| 18_20 | 363 | 189 | 53.60 | 1.31 |
| mag20 | 652 | 238 | 21.69 | 1.89 |
| mag21 | 1,553 | 473 | 22.53 | 1.92 |
| mag22 | 4,275 | 1,114 | 20.35 | 2.00 |
| mag23 | 11,085 | 2,438 | 9.80 | 1.96 |
| mag24+ | 27,783 | 4,815 | 7.26 | 1.86 |

**This closes the kNN-bleed / support-mask problem at its source.** Production NEO clone support hit a
hard ceiling at |v| ≈ 2.0 in *every* bin; the NEOMOD3 cache reaches **4.4–53.6 deg/day**, i.e. real
NEO density exists throughout the region where the ±5 pilot could only produce a support-mask default.
It also covers the 28 stranded night-61642 NEOs (max |v| 4.75) with genuine clones.

*Caveat on the counts:* this is the **antisun** center, the densest one — it holds ~2.5x the all-sky
average in-map density (45,776 actual vs 18,461 predicted from a flat 6.7% sky fraction), because
near-opposition objects are brighter. Off-antisun centers will get proportionally fewer; worth
re-checking per-bin counts at a high-elongation, high-latitude center before the full grid build.

### §4.7 END-TO-END SMOKE TEST — PASSED (2026-07-30)

`neomod/pipeline/test_neomod_map_smoke.py`: one center (`dlon+020_lat-12`), one mag bin (23–24),
±5 grid (1001²), all four populations, via `build_cloned_maps_for_center_magbin`. **75 s** →
~10 min for a full 8-bin map.

| pop | magcut_count (all-sky) | visible in patch | ∫ρ | **max \|v\| with support** |
|---|---:|---:|---:|---:|
| MBA | 1,023,722 | 162,276 | 185,244 | 1.60 |
| **NEO** | 64,442 | **9,449** | **1,199.5** | **5.00** |
| TNO | 16,203 | 2,852 | 26,776 | 0.03 |
| Trojans | 30,958 | **0** | 0 | — |

**1. THE VELOCITY CEILING IS GONE — the headline result.** Production S3M+GMM NEO clone support stops
at **|v| = 2.00** in every bin; NEOMOD3 support reaches **5.00**, the grid edge. The region that
produced the ±5-pilot support-mask artifact now carries **genuine NEO density**.

**2. The absolute normalisation is independently validated.**

| | N(H < 25) |
|---|---:|
| NEOMOD3 (sum of datacube weights) | 350,160 |
| S3M catalog | 268,511 |
| **ratio** | **1.30×** |

Two independently-built models agreeing to **30%** confirms the weights are absolute object counts and
that `w = 0.11433 objects/clone` is right — a units or relative-vs-absolute error would show up as
orders of magnitude. In the mag23 bin the new map predicts **1,199 NEOs vs S3M's 335 (3.6×)**: the
1.30× model difference compounded with the H>25 faint population S3M cannot represent.

**3. Caveats, all pre-existing or benign:**
- **kNN estimator has a ~11% positive bias**: ∫ρ_NEO = 1,199 vs the 1,080 implied by
  `n_clones / effective_factor`. Matches the production `map/true` median of **1.12** — so the §6.1
  acceptance test should expect **≈1.1, not exactly 1.0**.
- **TNO ∫ρ over-estimated ~9×** (26,776 from 2,852 objects). **Pre-existing** — the production
  diagnostic showed the same 2–4× effect. Cause: TNOs occupy a minuscule velocity region
  (max |v| = 0.03 deg/day) where the k=10 kNN estimate spikes. It sits at v≈0 where P_NEO ≈ 0 anyway.
- **Trojans = 0 for a SKY-GEOMETRY reason, not a cloning one:** 30,958 pass the magnitude cut all-sky
  but **zero** fall in this patch (Trojans cluster at Jupiter's L4/L5). Confirms §6.4 is a non-issue.

### §4.8 FIRST REAL MAPS BUILT + ACCEPTANCE PASSED (2026-07-31)

Driver wiring: `sorcha_gen_maps_grid.py` gained **`--vdp-module`** (default = the production GMM
pipeline, unchanged; opt-in to `velocity_density_pipeline_neomod_clone_only`). Slurm:
`neomod/pipeline/slurm/neomod3_maps.sbatch`. Two centers at ±5, into a NEW dir
`prob_maps_grid_neomod3_vlim5/` — **4 min each**, 156 MB each, no errors. Both caches in play:
Stage-0 n-body for MBA/TNO/Trojans (uncloned), NEOMOD3 projection cache for NEO (absolute norm).

⚠️ **PARTITION:** these jobs were first submitted to `ckpt-all` and were **preempted** (512/559 nodes
allocated), dying right after the 8 GB cache load and requeueing with a `(BeginTime)` backoff. Moved to
**`cpu-g2-mem2x`** (non-preemptible, astro account, 50 CPU / 707 GB free) → ran immediately. Memory
request 150G→80G (the cache needs ~20 GB). **Any job with an expensive startup should avoid ckpt** —
same lesson as `NEOMplanHYAK.md` §1.3–1.4 and `SORCHA_V5_PIPELINE.md` "Why no --requeue".

**A1. The velocity ceiling is GONE in every magnitude bin.**

| bin | prod max\|v\| | **new max\|v\|** | NEO clones (new) |
|---|---:|---:|---:|
| 14_16 | **0.00** | **4.43** | 7 |
| 16_18 | 1.79 | 4.86 | 53 |
| 18_20 | 1.31 | 4.81 | 323 |
| mag20 | 1.89 | **5.00** | 571 |
| mag21 | 1.92 | 4.97 | 1,437 |
| mag22 | 2.00 | 4.99 | 4,142 |
| mag23 | 1.96 | 4.99 | 10,980 |
| mag24+ | 1.86 | 4.99 | 27,733 |

(`dlon+000_lat+00`; `dlon+020_lat-12` is equivalent, 3.77–5.00.) Production had **zero** NEO support in
`14_16`; the new maps have genuine clones there. The ±5-pilot support-mask artifact is fixed **at the
source** — real density now exists where previously only the mask produced P_NEO = 1.

**A2. The absolute normalisation validates itself — a stronger test than planned.**

| bin | NEO frac (new) | NEO frac (prod) | ratio |
|---|---:|---:|---:|
| 18_20 | 0.00166 | 0.00171 | **0.97** |
| mag20 | 0.00205 | 0.00207 | **0.99** |
| mag21 | 0.00266 | 0.00190 | 1.40 |
| mag22 | 0.00409 | 0.00193 | 2.12 |
| mag23 | 0.00604 | 0.00165 | 3.66 |
| mag24+ | 0.00798 | 0.00211 | 3.78 |

In the **bright** bins — where S3M and NEOMOD3 both have the objects — the two independently-normalised
maps agree to **1–3%**. **If `w = 0.11433 objects/clone` were wrong by any constant factor, EVERY bin
would be off, including those.** They are not. Divergence appears **only** in the faint bins, growing to
~3.8×, exactly where S3M truncates at H=25 and NEOMOD3 supplies the population it cannot represent.
Together with the independent 1.30× N(H<25) cross-check (§4.7), option (b) is confirmed.

⚠️ **SCALING NOTE before any 667-map build:** each map task currently loads the **whole 97M-row / 8 GB**
cache to extract ~50k clones for its own patch. Fine for 2 maps; **wasteful ×667**. Pre-slice the cache
per center (small per-center files) before a full-grid run.

### §4.9 SCORING + ROC — NEOMOD3 beats production VDP *and* digest2 at both test centers (2026-07-31)

`neomod/pipeline/nm3_score_roc.py` → `outputs/neomod3_score_roc/`
(`neomod3_roc_summary.csv`, `neomod3_scored_tracklets.parquet`, `neomod3_roc.png`).
Sorcha case1 tracklets re-scored against `prob_maps_grid_neomod3_vlim5/` with production settings
(`support_mask_min=1`, nearest-dist mask OFF, Johnson-V), compared to the STORED production
`P_NEO_vdp_Vband` (±2 S3M/GMM maps) and `P_NEO_d2`.

| center | classifier | AUC | best F1 | compl% | contam% |
|---|---|---:|---:|---:|---:|
| **antisun** `dlon+000_lat+00` (n=1,876, 212 NEO) | **NEOMOD3 (±5)** | **0.9725** | **0.8768** | 84.0 | 8.2 |
| | production VDP (±2) | 0.9367 | 0.8615 | 80.7 | 7.6 |
| | digest2 | 0.9405 | 0.8145 | 79.7 | 16.8 |
| **busiest** `dlon+020_lat-12` (n=5,853, 1,066 NEO) | **NEOMOD3 (±5)** | **0.9450** | **0.8716** | 82.5 | 7.6 |
| | production VDP (±2) | 0.9200 | 0.8583 | 81.0 | 8.7 |
| | digest2 | 0.9373 | 0.8471 | 76.9 | 5.7 |

AUC **+0.025 to +0.036** over production VDP; NEOMOD3 also edges digest2 at both centers
(+0.032, +0.008).

**THE STRANDED-NEO PROBLEM IS SOLVED.** Of the tracklets at these centers with |v| > 2:

| | |
|---|---|
| count | **33 — every one a NEO**, zero non-NEO |
| production VDP scored exactly 0 | **100 %** |
| NEOMOD3 scores exactly 0 | **12.1 %** |
| median P_NEO | **0.0000 → 1.0000** |

**88 % of previously-unscoreable NEOs recover**, and from **genuine clone density** (support out to
|v| = 5.0, §4.8) rather than the support-mask default that made the ±5 pilot's recovery unacceptable.
This is ±5-pilot acceptance criterion #2 — which **FAILED** there (`vdp_bleed_diagnostic_bundle`) —
now **PASSING**. Because |v| > 2 is empirically pure-NEO, this is free recall at no contamination cost.

⚠️ **Two caveats on how far these claims go:**
1. **Two centers, not full sky.** The Stage-2 full-sky comparison (648k tracklets,
   `mag245_digest2_validation_review.ipynb`) had digest2 clearly ahead of VDP (0.943–0.953 vs
   0.866–0.881). These are ~7.7k tracklets in NEO-rich geometry. Encouraging and real, but **not yet a
   global "VDP beats digest2" claim** — that needs the full grid.
2. **The night-61642 panels are too small to quote.** Busiest center: 43 tracklets / 7 NEO (NEOMOD3
   returns a perfect AUC 1.000 / F1 1.000 — meaningless at that size). Antisun: 20 / 2, and the
   metrics guard correctly refused to compute it. Single-night scoring *looks* sensible; no curve
   should be drawn from it.

## §5. Files

### Created
| file | purpose |
|---|---|
| `neomod/src/velocity_density_pipeline_neomod_clone_only.py` | the new pipeline module |
| `neomod/pipeline/compare_s3m_vs_neomod3.py` | loads S3M + samples NEOMOD3 → 1D plot + KS stats |
| `neomod/pipeline/compare_s3m_vs_neomod3_extended.py` | H<25 cut, 2D panels, faint-tail figure |
| `neomod/pipeline/slurm/compare_s3m_vs_neomod3.sbatch` | runs the base comparison |
| `neomod/neomod3_vs_s3m_findings.ipynb` | executed findings notebook (5 figures, 4 tables, 0 errors) |
| `neomod/docs/new_neomod_cloning.md` | this document |

### Outputs (`outputs/neomod3_vs_s3m_comparison/`)
| file | contents |
|---|---|
| `s3m_vs_neomod3_1D.png` | full-sample 6-panel 1D comparison |
| `s3m_vs_neomod3_1D_Hcut25.png` | **H<25 matched** 6-panel comparison |
| `s3m_vs_neomod3_2D_full.png` | a–e, a–i, q–e: S3M / NEOMOD3 / contour overlay |
| `s3m_vs_neomod3_2D_Hcut25.png` | same, H<25 matched |
| `neomod3_adds_faint_tail.png` | H coverage + the faint population S3M lacks |
| `s3m_vs_neomod3_stats_full_and_Hcut.csv` | KS + summary stats, both samples |
| `s3m_source_neo.parquet` | 268,511 S3M NEOs (a,q,e,i,H,M) |
| `neomod3_samples.parquet` | 485,084 NEOMOD3 samples (+node,argperi,t_p) |

### Untouched
`velocity_density_pipeline_gmm.py`, all `prob_maps_*` directories, all production parquets.

---

## §6. Open decisions

### §6.1 Cross-population normalisation — **MEASURED 2026-07-30. Production is CORRECT; the job is to PRESERVE it.**

**Do not repeat the two mistakes below.** They cost time and produced two wrong conclusions that were
retracted:

1. ❌ *"Compare ∫ρ to `magcut_count`."* `magcut_count` in the npz is `len(df_sel)` — the magnitude cut
   applied to the **whole-sky** population; the 30° sky cut happens **afterwards**
   (`build_cloned_maps_for_center_magbin`, ~line 1585). Comparing a patch quantity to an all-sky one is
   meaningless. It produced a bogus "NEO recovers 3× less than MBA, matching the old GMM bug" claim —
   **that was wrong**; MBAs are simply far more concentrated toward the ecliptic/opposition than the
   more isotropic NEOs, so a 30° antisun patch holds a bigger *fraction* of them. Physics, not a bug.
2. ❌ *"NEOs cycle through the patch over 2 yr, inflating the detected fraction."* Plausible, and
   **disproven** — restricting to a single night gives the same ratio.

**The correct test** (`neomod/pipeline/calibration_check_v2.py`): compute the ground truth *directly*
from the Stage-0 cache — count objects per population actually inside the 30° patch, in that magnitude
bin, at the map epoch — and compare three quantities:

| quantity | meaning | tests |
|---|---|---|
| TRUE visible (Stage-0 cache) | what the map *should* encode | ground truth |
| MAP implied (∫ρ dv) | what the map *does* encode | **normalisation** |
| DETECTED (Sorcha tracklets) | what the survey saw | detection efficiency (physical) |

**RESULT — production n-body maps, `dlon+020_lat-12`, full 2-yr:**

| mag bin | true NEO frac | map NEO frac | **map/true** | det NEO frac | det/true |
|---|---:|---:|---:|---:|---:|
| 16_18 | 0.00237 | 0.00266 | 1.12 | 0.125 | 52.7 |
| 18_20 | 0.00092 | 0.00164 | 1.77 | 0.184 | 199 |
| mag20 | 0.00151 | 0.00241 | 1.60 | 0.177 | 117 |
| mag21 | 0.00153 | 0.00196 | 1.28 | 0.173 | 113 |
| mag22 | 0.00177 | 0.00191 | 1.08 | 0.158 | 89 |
| mag23 | 0.00203 | 0.00211 | **1.04** | 0.207 | 102 |
| mag24+ | 0.00211 | 0.00213 | **1.01** | 0.229 | 109 |

- **NORMALISATION: map/true median 1.12×, → 1.01 in the well-populated faint bins.** The existing
  `effective_factor` mechanism encodes the population mix **correctly**. There is **no normalisation
  bug**. The 1.6–1.8 values are the sparse bright bins (19–43 true NEOs) — small-number noise.
- **DETECTION: det/true ≈ 109×** — Sorcha detects ~13% of visible NEOs but ~0.09% of visible MBAs in
  the same bin. This is the **survey selection function, not a normalisation error**: each bin is 1 mag
  wide while LSST's efficiency collapses across it, and MBAs (rising far more steeply with magnitude)
  pile up at the faint edge where they are missed.

**So the question changed.** It is not "which normalisation is right" but **"how do we preserve a
working normalisation when the NEO clone source changes?"** The map's NEO amplitude must equal the
**true visible NEO count** in that patch+bin. Today that comes from S3M's count and it verifiably
works — but S3M has **zero** H>25 objects, so once the faint population matters its count is an
undercount. Only NEOMOD3's own absolute normalisation supplies the right number → **option (b)**.

**ACCEPTANCE TEST for the new pipeline:** rerun `calibration_check_v2.py` on the NEOMOD3 maps and
require `map/true ≈ 1`.
⚠️ **Interpretation caveat:** the Stage-0 cache used as ground truth is *itself* S3M-derived and so also
lacks H>25 NEOs. When the NEOMOD3 pipeline correctly includes the faint population, `map/true` **will
legitimately exceed 1**, by roughly the faint-NEO share. That is success, not failure — predict the
expected excess in advance so it is not misread as an error.

**Also found:** Trojans contribute **exactly zero density in every magnitude bin** at these centers in
the *production* maps already (`support_count__Trojans__* = 0`). So §6.4's worry about losing Trojans
by removing cloning is moot — they were never contributing.

#### The original framing (kept for context)

### §6.1b Cross-population normalisation — the original three options

`P(NEO) = ρ_NEO / Σ_pop ρ_pop` is only meaningful if the densities share a footing. Previously every
population was cloned from the **same S3M catalog** and downweighted by `n_clones/n_source`, so the
footing was implicit and consistent. Now **NEO density comes from an external absolute model** while
MBA/TNO/Trojan come from **real object counts**. Those are not automatically on the same scale — and
the ratio is exactly what the classifier reads.

The module keeps the existing `effective_factor` behaviour (a faithful structural port) and does
**not** invent a normalisation. Options:

- **(a) survey-prior rescale** — the `ARNOR_NOTE_GMM_BUG.md` fix: normalise each population's map to
  integrate to its true population count.
- **(b) NEOMOD3-absolute** — use NEOMOD3's own `D > Dref` normalisation (871.05 NEOs with D > 1 km)
  for NEO and a matched absolute count for the others.
- **(c) empirical** — calibrate the NEO:non-NEO ratio against a known-truth sample.

**Until this is settled, P(NEO) from this module is a SHAPE result, not a calibrated probability.**

#### §6.1c ADOPTED: option (b), NEOMOD3-absolute — and the number it needs

The NEOMOD3 datacube weights are **absolute object counts**, not relative weights. Verified:

- `sum(array4D) = 11,432,918` = the true number of NEOs with **15 < H < 28**.
- Cross-check against the header's own stated normalisation (`1.0 871.05 ! reference diameter Dref,
  # of NEOs with D>D_ref`): H ≈ 17.75 corresponds to D = 1 km at a typical albedo
  (`H = -5log10(D√pV/1329)`, pV ≈ 0.14), and the weights summed up to that H give **~865** vs the
  stated **871.05**. Confirms the weights are real counts.

Therefore each cached NEOMOD3 clone represents

```
w_abs = total_weight / n_draws = 11,432,918 / 100,000,000 = 0.11433 real NEOs per clone
```

and the implementation is a one-line change to the existing machinery:

```
effective_factor(NEO) = n_draws / total_weight = 8.747      # instead of n_visible/n_source
rho_NEO = density_of_clones / effective_factor              # -> integrates to the TRUE NEO count
```

This makes the NEO amplitude independent of both the sampling budget *and* S3M, which is the whole
point. `n_draws` and `total_weight` are recorded in
`outputs/neomod3_projection_cache/cache_metadata.json` so the map builder never has to re-derive them.

### §6.2 H sampling range — **RESOLVED 2026-07-30: do NOT cut. Keep the full 15–28 range.**

The question was whether to cut NEOMOD3 at `H < 25`, matching S3M's truncation, the maps' 14–25
magnitude bins, and LSST's ~24.5 limit — since "we won't observe fainter objects anyway".

**Measured answer: no.** `H` is **absolute** magnitude; the map bins and the LSST limit are
**apparent**. For NEOs the two decouple, because NEOs come close:
`mag_app ≈ H + 5·log10(r·Δ) + phase`. At Δ ≈ 0.05 AU an H = 28 object appears at mag ≈ 21.5.

2M NEOMOD3 orbits were sampled and projected through the real pipeline machinery
(`test_neomod3_H_cut_impact.py` → `outputs/neomod3_H_cut_test/`):

- objects landing in the maps' 14–25 **apparent** range span the **full** H range **15.23 → 28.00**;
  there is no safe H cut.
- **`H<25` would delete 35.7% of all detectable NEOs and 81.2% of the detectable |v|>2 NEOs.**
  (`H<26` → 63.2%; `H<27` → 35.5%.)
- the faint population is **twice as fast**: in-map median |v| **1.55** deg/day (H≥25) vs **0.76**
  (H<25); fraction above |v|=2 is **38.3%** vs **4.9%**. In-map objects with |v|>2 have median
  **H = 26.56**.

Individually faint objects are rarely visible (0.11% of H≥25 samples land in-range), but there are so
many of them that they still supply **35.7% of everything detectable**.

**This is the fix for the kNN-bleed problem.** Measured on the production n-body maps, S3M-derived
GMM NEO clone support **stops at |v| ≈ 1.8–2.0 deg/day in every magnitude bin** — it never crosses the
±2 grid boundary at all. That is the root cause of the ±5-pilot artifact
(`vdp_bleed_diagnostic_bundle`): P_NEO = 1.0 beyond |v| = 2 came from the support mask, because there
was never any NEO density out there. Cutting NEOMOD3 at H<25 would **re-create exactly that problem**;
keeping the faint end is what populates |v| > 2 with genuine clone support.

The real cost of keeping it is **sampling efficiency**, not correctness — see §6.5.

### §6.3 Sampling fidelity (low priority)

Adopt the reference Fortran's trilinear-interpolation + rejection scheme instead of bin-pick +
uniform dither? Only affects bin-edge blockiness; the kNN estimator smooths it downstream. Probably
not worth it.

### §6.4 Does removing MBA/TNO/Trojan cloning leave enough support?

Cloning existed to give the kNN estimator enough points. With clone_factor = 1 the MBA density is
built from real visible objects only. MBA is numerous so this is likely fine; **TNO and Trojans are
not** (TNO ~48k and Trojans ~180k total, before the per-patch sky cut and mag-bin split).
`DEFAULT_MIN_CONDITIONAL_CLONER_INPUT` already returns zeros when a population has too few visible
objects — expect **more empty TNO/Trojan maps** than before. Worth checking on the first test map.

---

### §6.5 Sampling efficiency — solved by a projection cache (§4.6)

Draw cost is steeply H-dependent (measured, per object landing in the maps' 14–25 range):

| H band | draws per in-map object |
|---|---:|
| 15–25 | 16 |
| 25–26 | 189 |
| 26–27 | 631 |
| 27–28 | 2,208 |

Sampling per (center, mag-bin) is therefore hopeless — a 30° patch is only ~6.7% of the sky. **But the
projection is epoch-dependent and NOT center-dependent**, so it is done **once, globally**, and all 667
centers slice the same cache. See §4.6.

## §7. Reproduce

```bash
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python

# base comparison (loads S3M, samples NEOMOD3, 1D plot + stats)
sbatch neomod/pipeline/slurm/compare_s3m_vs_neomod3.sbatch

# extended: H<25 matched cut, 2D panels, faint-tail figure (reads the saved parquets — fast)
$PY neomod/pipeline/compare_s3m_vs_neomod3_extended.py

# import check on the new module
PYTHONPATH=neomod/src $PY -c "import velocity_density_pipeline_neomod_clone_only as v; print(v.DEFAULT_POPULATION_SETTINGS)"
```

### Traps hit while doing this (all cost real time)

1. **`load_s3m_population()` defaults to the HYBRID catalog** (`t_0` = 2023), not pure S3M
   (`t_0` = 2008). Must set `VDP_LOADER=s3m`. The first comparison run silently used 295,040 hybrid
   objects instead of 268,511 S3M ones. See `fixing_integrator.md` §13.6.
2. **`s3m_loader` searches `['.', 'S3Mdata']` relative to CWD** — a script must `os.chdir` into
   `neomod/` regardless of the Slurm `--chdir`.
3. **`/tmp` is node-local on Hyak** — a script written to `/tmp` on the login node is invisible to the
   compute node. Put scripts on GPFS.
4. **SWRI serves a self-signed cert** — `curl -k` / `ssl.CERT_NONE`.
5. Standing rule (`NEOMplanHYAK.md` §1.4): **no data-touching step is small enough for the login
   node** — everything goes through Slurm.

---

## §8. HEALPix-partitioned NEOMOD3 cache — detailed plan (PREREQUISITE for the 667-map build)

**Status:** planned, not yet built. **Nothing already produced depends on this** — the two ±5 maps
(§4.8) and all results in §4.9 were built from the monolithic cache and stand on their own. This is
purely an I/O optimisation to make a 667-map build affordable.

### §8.1 The problem, measured

Every map task currently loads the **entire 97,056,755-row / 7.88 GB** cache to extract the ~50k clones
inside its own 30° patch.

| | |
|---|---|
| I/O for a 667-map build | **667 × 7.88 GB ≈ 5.3 TB** |
| RAM per task | ~20 GB (pandas overhead on a 7.88 GB parquet) |
| consequence | jobs are slow to start and **expensive to lose to preemption** — exactly what killed the first `ckpt-all` attempt (§4.8) |

### §8.2 Why per-CENTER slicing fails but per-PIXEL partitioning works

- **Per-center slicing (rejected).** The 667 patches overlap heavily: each is 6.9% of the sky, so
  667 × 6.9% ≈ **45 sky-coverages** — every clone falls inside ~45 different patches. Materialising
  per-center files would write **~4.4 billion rows**.
- **Per-pixel partitioning (adopted).** HEALPix tiles the sky into **disjoint** cells, so each clone is
  written **exactly once**. Total stays 97M rows / 7.88 GB — **zero duplication**. The overlap moves
  from *write* time to *read* time (several centers read the same pixel file), which is free.

### §8.3 Design — nside = 8

Measured pixels touched by a 30° disc (`query_disc(..., inclusive=True)`, averaged over 5 centers
spanning the antisun, the ecliptic, high latitude and the poles):

| nside | npix | touched by 30° disc | fraction | read/task |
|---|---:|---:|---:|---:|
| 4 | 192 | ~23 | 12.1% | 0.95 GB |
| **8** | **768** | **~72** | **9.4%** | **0.74 GB** |
| 16 | 3,072 | ~245 | 8.0% | 0.63 GB |
| *(floor: a 30° disc is 6.9% of the sky)* | | | 6.9% | 0.54 GB |

**nside = 8 chosen:** 768 partitions × ~126,376 rows (~10.3 MB) each. It captures nearly all the
available gain (0.74 vs a 0.54 GB floor) without the filesystem overhead of 3,072 small files.
**Net: ~11× less I/O and RAM per task.**

`max_pixrad(8) = 7.473°` — this is the safety margin added to the query radius (§8.4 step 2).

### §8.4 Implementation

**Step 1 — build the partitioned cache** (one pass, ~10 min; new script
`neomod/pipeline/neomod3_cache_healpix.py`):
```python
pix = hp.ang2pix(8, ra_deg, dec_deg, lonlat=True)     # disjoint tiles; ICRS ra/dec, degrees
# -> outputs/neomod3_projection_cache/by_pixel/pix=NNN/part.parquet   (768 files)
```
The monolithic parquet is **kept, untouched**.

**Step 2 — read only the needed partitions** at map-build time:
```python
vec    = hp.ang2vec(center_lon_ecl_to_radec..., lonlat=True)   # SAME center the sky cut uses
radius = np.radians(max_sep_deg + np.degrees(hp.max_pixrad(8)))  # 30 + 7.473 deg
pixels = hp.query_disc(8, vec, radius, inclusive=True)           # conservative superset
```
then load only those partitions and hand the result to `build_visible_subset_dataframe`, which applies
the **exact, unchanged 30° angular cut**.

**Step 3 — module wiring:** `velocity_density_pipeline_neomod_clone_only._load_neomod3_cache()` gains
an optional partitioned path (`NEOMOD3_CACHE_DIR`), selected by env/flag. Default stays the monolithic
file, so **the existing behaviour is the fallback**.

### §8.5 Why this cannot change the science

Three independent reasons, by construction:

1. **Pixel selection only has to be a SUPERSET.** `inclusive=True` returns every pixel that overlaps
   the disc even partially, and §8.4 adds a further `max_pixrad` margin.
2. **The exact 30° cut still runs afterwards**, unchanged. Reading *extra* clones is therefore
   harmless — they are cut. The **only** possible failure is reading too *few*, which (1) prevents.
3. **The clones are byte-identical.** Only their file location changes; no value is recomputed.

The residual risk is not to correctness but to *my code* — a subtle bug in the pixel-query or
coordinate handling. §8.6 is designed so any such bug is **loud** (a failed assertion) rather than
**silent** (a slightly under-dense patch edge).

### §8.6 VALIDATION PLAN — four tiers, each gating the next

**T1 — partition integrity** (before any map is built):
- **row count**: `sum over 768 partitions == 97,056,755` exactly;
- **no duplication / no loss**: column checksums (`sum`, `min`, `max` of `vlam`, `vbeta`, `mag_app`, `H`)
  must match the monolithic cache to float precision;
- **pixel assignment**: for a random 1M sample, recompute `hp.ang2pix(8, ra, dec, lonlat=True)` and
  assert it equals the partition the row was stored in (catches any ra/dec vs lon/lat, degree/radian,
  or frame mix-up);
- **partition count**: exactly 768 directories, none empty-by-accident (an empty pixel is legitimate
  only if the sky there is genuinely unpopulated — assert against the monolithic per-pixel histogram).

**T2 — clone-set equality per center** (the decisive test of the pixel query):
for **12 centers** chosen to span the extremes — antisun, `dlon±140` (sun-exclusion edge),
`lat 0, ±25, ±50` (poles), plus the two already-built test centers:
- read via partitions → 30° cut → set **A**
- read monolithic → 30° cut → set **B**
- assert `len(A) == len(B)` **and** the sorted `(vlam, vbeta, mag_app)` arrays are **element-wise
  identical**. Not a spot check — exact equality.

**T3 — end-to-end map regression** (the strongest test; we have known-good maps to reproduce):
rebuild **`dlon+000_lat+00`** and **`dlon+020_lat-12`** from the partitioned cache and compare against
the existing `prob_maps_grid_neomod3_vlim5/` maps built from the monolithic cache:
- **T3a, `--n-jobs 1`: require BIT-IDENTICAL** arrays (`density_raw`, `support_count`, `nearest_dist`,
  for every population × magnitude bin). Since T2 proves the input clone set is identical and the
  estimator is deterministic single-threaded, *any* difference here means a real bug.
- **T3b, `--n-jobs 16`: require agreement within the documented joblib scatter.** Map generation is
  **not** bit-reproducible in parallel — `SORCHA_V5_PIPELINE.md` records "known joblib parallel
  non-determinism (~0.02 F1 between regenerations)". T3a with `n_jobs=1` exists precisely to remove
  that confound; T3b confirms the parallel path stays inside the known envelope.

**T4 — scoring regression**: re-score the tracklets at those two centers with the rebuilt maps and
require the §4.9 numbers to reproduce (AUC 0.9725 / 0.9450, F1 0.8768 / 0.8716, and the 33 |v|>2 NEOs
still recovering at 88%).

**Gate:** the 667-map build starts only after T1–T4 all pass.

### §8.6b VALIDATION RESULTS — T1–T4 all PASS (2026-07-31). Gate cleared.

| tier | result | evidence |
|---|---|---|
| **T1** partition integrity | **PASS** | 97,056,755 rows exact; column checksums rel-diff 0–3e-16; 1M-sample `ang2pix` re-derivation → **0** mismatches; 768/768 dirs populated |
| **T2** clone-set equality | **PASS** | all **12** centers element-wise identical, 98–106 pixels touched each (job 37970…) |
| **T3a** bit-identical @ `n_jobs=1` | **PASS** | job 37979790 — `density_raw`, `support_count`, `nearest_dist` × 4 populations, **max&#124;diff&#124; = 0.000e+00** on all 12 arrays |
| **T3b** parallel @ `n_jobs=16` | **PASS (exceeded)** | job 37981774 — worst relative array difference **0.000e+00**, i.e. bit-identical, *better* than the "within joblib scatter" bar §8.6 set |
| **T4** scoring regression | **PASS** | AUC 0.9725 / 0.9450, F1 0.8768 / 0.8716, night-61642 AUC 1.000 — reproduce the §4.9 CSV to every printed digit; 33 &#124;v&#124;>2 NEOs still 87.9% recovered (12.1% zero) |

**Read T3a's scope honestly.** Only **NEO** is an informative comparison: MBA/TNO never touch the
NEOMOD3 cache (their identity was guaranteed), and Trojans returned all-zero maps in *both* paths
(0 objects visible within 30° of that center at mag 23–24), so that pair is vacuous. The real
evidence is NEO — **9,449 clones**, identical to the last bit. That is sufficient, because NEO is the
only population the HEALPix change can possibly affect.

**T3b came out stronger than planned.** §8.6 only required "within the documented joblib scatter",
because `SORCHA_V5_PIPELINE.md` records ~0.02 F1 non-determinism between regenerations. In fact the
16-way parallel rebuild was *bit-identical* to the monolithic reference. So for this module + fixed
seed the parallel path is deterministic, and the historical scatter is not reproduced here. Do **not**
generalise that to the whole pipeline — treat it as a bonus, not a new invariant.

#### One real semantic difference: `magcut_count__NEO__*` in the npz

The **only** field that differs between cache paths, by design:

| path | `magcut_count__NEO__mag23` @ antisun | meaning |
|---|---:|---|
| monolithic | 64,442 | NEOMOD3 objects in the mag bin **over the whole sky** |
| HEALPix | 16,373 | …**within the ~100 pixels loaded** (~13% of sky, but ~25% of NEOs — they crowd the ecliptic) |

Both then apply the identical 30° sky cut and hand the estimator the **same 9,449 clones**, which is
why every science array is identical. `magcut_counts` is written to the npz and read back as
metadata; **grep confirms no numerical consumer** anywhere in the module or pipeline (only dev
notebooks print it). Left as-is deliberately — making it all-sky under HEALPix would require reading
the whole cache, defeating the purpose.

> ⚠️ **This field bit us once already.** In the earlier calibration check, `magcut_count__NEO` was
> read as a *patch* count and compared against ∫ρ over the patch, producing a bogus "NEO
> under-weighted 3×" claim that had to be retracted. Under HEALPix it is now closer to a patch count
> but still **not** one. Never use it as a population normalisation — use `outputs/…/cache_metadata.json`
> (`w_abs_objects_per_clone`, `effective_factor_NEO`) instead.

### §8.7 Failure modes and where each is caught

| failure mode | caught by |
|---|---|
| pixel query selects too few pixels (under-read → thin patch edges) | **T2** (clone-count/element mismatch) |
| ra/dec ↔ lon/lat swap, degrees ↔ radians, wrong frame | **T1** (pixel-assignment check) and **T2** |
| rows lost or duplicated while writing partitions | **T1** (row count + checksums) |
| a partition silently missing at read time | **T1** (768-directory assert) and **T2** |
| subtle change in the density estimate | **T3a** (bit-identical requirement) |
| parallel non-determinism mistaken for a bug | **T3a vs T3b** separation |

### §8.8 Rollback

The monolithic cache is **kept**, and the partitioned path is **opt-in** (default unchanged). If any
tier fails, drop back to the monolithic reader with zero loss — the only cost is that a 667-map build
stays expensive.

### §8.9 Cost / benefit

| | monolithic | HEALPix nside=8 |
|---|---:|---:|
| cache on disk | 7.88 GB | 10.22 GB (partitioning costs ~30%: smaller row groups compress worse) |
| read per map task | 7.88 GB | **~1.38 GB** (~100 of 768 pixels) |
| I/O for 667 maps | **5.3 TB** | **0.92 TB** |
| preemption cost | high (expensive startup lost) | low |
| build cost | — | one pass, 33 s |

**Corrected 2026-07-31.** An earlier draft of this table claimed 0.74 GB/task and an 8–11× saving.
That was wrong: it applied the pixel fraction to the *monolithic* size and ignored the fact that the
partitioned cache is 30% larger. Measured reduction is **5.7×**, not 8–11×. Still decisive — it is
what makes the full grid viable on preemptible `ckpt-all` with skip-existing resume, instead of
~15 h serialised on the 64-CPU-capped non-preemptible partition.

---

## §9. FULL 667-MAP GRID BUILD — launched 2026-07-31 (job 37983895)

Steps 1–3 of the old plan are done (§6.1, §6.2, §4.8, §4.9); the T1–T4 gate is cleared (§8.6b), so
the grid rebuild is running.

**Submit:** `sbatch neomod/pipeline/slurm/neomod3_grid_full.sbatch`

| setting | value | why |
|---|---|---|
| output | `prob_maps_grid_neomod3_full/` | **new** dir; `prob_maps_grid_s3m_nbody/` never touched |
| array | `0-666%200` | 3,567 idle CPUs on `ckpt-all` ÷ 16 = ~220 slots; %200 → ~4 waves |
| partition | `ckpt-all` | preemptible, but 3,567 idle vs 81 on `cpu-g2-mem2x` (astro capped at 64 CPU there = 4 concurrent tasks = ~16 h) |
| cpus / mem | 16 / 48 G | measured peak **25.5 GB for two concurrent builds** → ~13 GB each; 48 G is 3.7× headroom |
| `NEOMOD3_CACHE_DIR` | `…/by_pixel` | the whole point — ~1.4 GB read per task, not 7.9 GB |
| `--ref-obstime` | `2027-08-25T00:00:00` | **must be explicit** |
| `--cache` | `epoch_state_2027-08-25T000000.parquet` | Stage-0 n-body states for MBA/TNO/Trojans |
| `--mba-clone-factor 1` | no cloning | the whole premise of this module |
| `--velocity-grid-limit 5.0` | ±5 deg/day | recovers the |v|>2 NEOs |
| `--overwrite` | **omitted** | makes the array idempotent — see below |

**Two traps that would have silently corrupted the entire grid:**

1. **`DEFAULT_REF_OBSTIME` in `sorcha_gen_maps_grid.py` is `2026-01-01T00:00:00`**, but the NEOMOD3
   projection cache is projected for **2027-08-25**. Omitting `--ref-obstime` would have built 667
   maps whose sky patches were selected at the wrong epoch — with no error raised. Always pass it.
2. **Preemption.** `ckpt-all` kills tasks. `--overwrite` is deliberately omitted so a task whose
   `.npz` already exists exits immediately (`sorcha_gen_maps_grid.py:240`). **Mop-up = resubmit the
   identical script**, as many times as needed, until 667 files exist. Do not add `--requeue`.

**Disk:** 667 × ~156 MB ≈ **104 GB**. GPFS reported **1.6 T free (100% used)** at launch — it fits,
but check `df -h /mmfs1` before any *further* large output, and note this is far tighter than the
8.7 T recorded in older docs.

### §9.1 BUILD RESULT (2026-07-31) — 667/667 in ~25 min, one truncated file

All 667 tasks reported `COMPLETED`; **108 GiB** written; no requeue/mop-up pass was needed.
Acceptance run by `neomod/pipeline/nm3_fullgrid_accept.py` (job 37993375):

| check | result |
|---|---|
| **A** file count | **667 / 667**, 107.9 GiB, none undersized |
| **B** validated centers unchanged | `dlon+000_lat+00` and `dlon+020_lat-12` **bit-identical** (worst rel diff 0.000e+00) to the T3b rebuild — the production array reproduces the validated build exactly |
| **C** NEO velocity support | over 61 sampled maps spanning the grid: `mag22` median &#124;v&#124;max **4.98**, `mag23` **4.99**, `mag24+` **5.00**; **61/61 maps reach past 2.0** in every bin |
| **D** full-grid scoring | **crashed** on a corrupt file — fixed, rerun in §9.3 |

**C is the headline.** The old ±2 grid pinned NEO support at exactly 2.00 everywhere (the ceiling that
stranded the fast NEOs). The new grid carries genuine NEO density out to ~5.0 at *every* sampled
center, not just the two test centers — so the |v|>2 recovery demonstrated in §4.9 is a grid-wide
property, not a local accident.

### §9.2 ⚠️ TRAP: skip-existing is keyed on file EXISTENCE, so a truncated map survives a mop-up

Check **D** died with `BadZipFile: File is not a zip file` on
**`prob_maps_grid_dlon+050_lat-25.npz`** (array task **77**): **97.3 MiB instead of ~156 MiB** — a task
killed mid-write on preemptible `ckpt-all`, leaving a truncated `.npz`.

This is the dangerous failure mode of the §9 idempotency design, and it defeats every cheap check:
- Slurm reported the array **fully COMPLETED** — no FAILED/CANCELLED task to notice.
- `ls | wc -l` says **667 / 667**.
- My own size guard (`< 10 MiB`) passed it at 97 MiB.
- **Worst: resubmitting the script would have SKIPPED it**, because `sorcha_gen_maps_grid.py:240`
  skips when the path exists. The grid would have stayed silently corrupt behind a green file count.

**Therefore the mop-up procedure is not "resubmit until 667 files exist".** It is:

```bash
# 1. integrity scan -- actually OPEN every npz and force a payload array read
conda_prep/bin/python neomod/pipeline/nm3_fullgrid_accept.py     # check A+B+C, then D
# 2. rebuild ONLY the bad indices, WITH --overwrite (skip-existing would refuse otherwise)
#    map filename -> array index via neomod/pipeline/slurm/grid_map_manifest.csv
sbatch ... sorcha_gen_maps_grid.py --task-id <IDX> ... --overwrite
```

A file-count check is necessary but **not** sufficient. Only a real open-and-read of every map is.
Rebuilt task 77 on the non-preemptible partition (job 37993986).

**Acceptance checks when it finishes:**
- `ls prob_maps_grid_neomod3_full/*.npz | wc -l` == **667** *and* all 667 open with 151 keys each
- the two centers already validated must still match `prob_maps_grid_neomod3_vlim5_hpx/`
- NEO |v| support reaches ~5.0 away from the two test centers
- full-grid scoring vs production VDP + digest2


> 📌 **All evaluation from 2026-08-01 onward follows `neomod/docs/EVALUATION_PROTOCOL.md` (v1.0,
> FROZEN).** TPR/FPR are primary; contamination and F1 are derived and prior-dependent; classifiers
> are compared at matched completeness or matched FPR, never at self-selected thresholds; invalid
> scores are NaN abstentions, never zero. The tables in §9.3–§9.9 predate that protocol and violate
> parts of it — read them with §11.6.

### §9.3 FULL-GRID SCORING RESULT (job 37994096) — beats production VDP everywhere, but does **NOT** beat digest2 grid-wide

> ⚠️ **EVERY NUMBER IN §9.3–§9.9 IS PROVISIONAL — the benchmark is the wrong population.**
> These tracklets are Sorcha detections of **S3M** objects, but the maps' NEO prior is **NEOMOD3**,
> which contains a faint (H>25) population S3M does not have at all (§9.10). Truth and prior describe
> different universes, so the digest2 comparison, the contamination analysis and the calibration
> results below cannot be taken at face value. **Re-measure on a NEOMOD3-based benchmark — see §11.1.**

648,769 tracklets, 617 centers, all 667 maps verified intact first.

| set | classifier | AUC | best F1 | completeness | contamination |
|---|---|---:|---:|---:|---:|
| full_2yr | **NEOMOD3 full grid (±5)** | **0.9288** | 0.8319 | 79.23% | 12.43% |
| full_2yr | production VDP (±2) | 0.8814 | 0.8020 | 73.62% | 11.94% |
| full_2yr | digest2 | **0.9430** | **0.8489** | 77.79% | **6.60%** |
| night_61642 | **NEOMOD3 full grid (±5)** | 0.9225 | 0.8111 | 74.88% | 11.53% |
| night_61642 | production VDP (±2) | 0.8927 | 0.7930 | 71.71% | 11.31% |
| night_61642 | digest2 | **0.9300** | **0.8143** | 72.44% | **7.04%** |

**Two conclusions, and the second one is unwelcome.**

1. **vs production VDP — a clear, unambiguous win.** AUC 0.8814 → 0.9288, F1 0.8020 → 0.8319,
   completeness 73.6% → 79.2%, at essentially unchanged contamination (11.9% → 12.4%). The new
   cloning scheme is strictly better than what it replaces. This holds on both tracklet sets.

2. **vs digest2 — we LOSE grid-wide, and this contradicts the two-center result.** §4.9 showed
   NEOMOD3 ahead of digest2 at both test centers (AUC 0.9725/0.9450 vs 0.9405/0.9373). Across all
   617 centers digest2 is **ahead**: AUC 0.9430 vs 0.9288, F1 0.8489 vs 0.8319, and — the biggest
   gap — contamination **6.60% vs 12.43%**, i.e. we produce roughly **twice** the false positives at
   the best-F1 operating point.

   **The two test centers were not representative** — but *not* for the reason first written here.
   An earlier draft of this section claimed they were "NEO-rich, which flatters a NEO-density
   classifier". **That was wrong and is retracted**: measured NEO fraction at the two test centers is
   **0.161**, *below* the grid mean of **0.277**. The real reason is **sky geometry**, quantified in
   §9.4 — both centers sit at |Δlon| ≤ 20°, |lat| ≤ 12°, which is the one region of the grid where
   NEOMOD3 does beat digest2. Do not quote the §4.9 head-to-head against digest2 as a general
   result; §9.3–§9.4 supersede it. The honest statement is: *the new cloning fixes the VDP, and the
   VDP still trails digest2 on overall ranking quality.*

**Where we are unambiguously better: the fast movers.**

| | count | production VDP | NEOMOD3 full grid |
|---|---:|---|---|
| |v|>2 tracklets, whole grid | **7,727** | scored **0** for 100.0% | scored 0 for **14.0%** |
| of which NEO | **7,727 (100%)** | median P = 0.0000 | median P = **1.0000** |
| of which non-NEO | **0** | — | — |

Every single one of the 7,727 tracklets beyond the old ±2 boundary is a **real NEO** — zero
contaminants. Production could not score them *at all* (off-grid ⇒ 0), so all 7,727 were structural
misses. The new grid recovers **86%** of them at median probability 1.0. That is ~6,600 NEOs that the
production pipeline was guaranteed to miss for a reason that had nothing to do with the science.

**What this argues for next** (not yet done):
- The gap vs digest2 is a **contamination** gap, not a completeness gap — NEOMOD3 is *more* complete
  (79.2% vs 77.8%) but twice as contaminated. Investigate the false positives at the best-F1
  threshold: which population, which magnitude bin, which velocity regime.
- A combined score (VDP × digest2) is the obvious thing to test — they fail differently, and VDP
  uniquely covers the |v|>2 regime where digest2's own short-arc behaviour is weakest.

### §9.4 CONTAMINATION PROBE — the sky-direction question, and the real cause

Script: `neomod/pipeline/nm3_contamination_probe.py`, `neomod/pipeline/nm3_support_mask_probe.py`.
Output: `outputs/neomod3_fullgrid/per_center_auc.csv`. Three hypotheses tested.

#### H2 — is the pooled comparison unfair to VDP? **NO.**

Worth ruling out first: VDP scores come from 617 *different* maps while digest2 is one global
classifier, so pooling could destroy AUC even with perfect per-center ranking (Simpson's paradox).
It does not — digest2 wins per-center too, by *more* than it wins pooled:

| | NEOMOD3 | digest2 |
|---|---:|---:|
| pooled AUC | 0.9288 | **0.9430** |
| per-center mean | 0.8921 | **0.9011** |
| per-center median | 0.8867 | **0.9187** |
| per-center, tracklet-weighted | 0.9049 | **0.9204** |

NEOMOD3 wins at **236 / 506** centers (46.6%). The deficit is real, not an artefact of pooling.

#### H1 — is it sky direction? **YES — this is why the test centers misled us.**

| |Δlon| from antisun | centers | tracklets | mean ΔAUC (NEOMOD3 − digest2) |
|---|---:|---:|---:|
| 0–20° | 88 | 195,137 | **+0.0114** ← *both test centers live here* |
| 20–50° | 132 | 228,724 | **−0.0537** |
| 50–90° | 167 | 150,078 | −0.0074 |
| 90–130° | 65 | 11,515 | +0.0306 |
| 130–180° | 32 | 6,563 | +0.0201 |

| |ecliptic lat| | centers | tracklets | NEO frac | mean ΔAUC |
|---|---:|---:|---:|---:|
| 0–2° | 115 | 117,140 | 0.132 | −0.0049 |
| 2–8° | 198 | 267,296 | 0.173 | +0.0020 |
| 8–18° | 101 | 182,927 | 0.295 | −0.0042 |
| 18–35° | 78 | 62,559 | 0.610 | −0.0307 |
| 35–50° | 14 | 12,508 | 0.934 | **−0.1108** |

The deficit concentrates in **two** regions: the 20–50° longitude annulus (the largest tracklet
population, 228k) and **high ecliptic latitude** (|lat| > 35°, where ΔAUC = −0.11). The worst
centers are almost all |lat| = 35–50, plus a cluster at Δlon ≈ +50, lat ≈ 0.

The two §4.9 test centers score **ΔAUC +0.0320** and **+0.0077** — both positive, both inside the
one favourable longitude band. **Your instinct was right: this is a sky-direction effect**, and the
test centers were chosen from the region where the new maps happen to win.

#### H3 — where do the false positives come from? **One velocity band, and it is INHERITED, not new.**

False positives at each classifier's own best-F1 threshold:

| |v|max band | non-NEOs | production VDP | NEOMOD3 | digest2 |
|---|---:|---:|---:|---:|
| 0–0.25 | 446,064 | 1,691 (0.38%) | 3,297 (0.74%) | 6,525 (1.46%) |
| **0.25–0.5** | **52,773** | **12,673 (24.01%)** | **12,950 (24.54%)** | **1,502 (2.85%)** |
| 0.5–1.0 | 1,158 | 484 (41.8%) | 491 (42.4%) | 145 (12.5%) |
| **total** | | **14,849** | **16,739** | **8,173** |

**The single band |v| ∈ (0.25, 0.5] contributes 77% of NEOMOD3's false positives** — and production
VDP is just as bad there (24.01% vs 24.54%). So the ~2× contamination gap versus digest2 is a
**pre-existing property of the VDP that the new cloning inherited, not something it introduced**.
Our change moved that band's FP rate by **+0.5 points**.

By population, NEOMOD3 is *better* than digest2 almost everywhere — TNO FP rate 4.4% vs 37.8%,
Trojan 0.09% vs 5.1% — and worse on MBA (2.6% vs 0.9%) and on `other` (24.7% vs 7.1%).

#### Mechanism (`nm3_support_mask_probe.py`, center `dlon+020_lat-12`, mag23, annulus |v| 0.25–0.5)

| | production (MBA cf=5) | NEOMOD3 (MBA cf=1) |
|---|---:|---:|
| MBA support < 1 ⇒ **ρ_MBA zeroed** | 75.1% of cells | 82.4% |
| NEO support < 1 (**exempt** from the mask) | 91.7% | 78.6% |
| resulting P(NEO) > 0.99 | **75.1%** | **82.4%** |

The support mask (`support_mask_min=1`) zeroes a population wherever its support < 1, **but NEO is
in `_support_mask_skip`** — the documented asymmetry from `vdp_bleed_diagnostic_bundle`. In the
0.25–0.5 band real MBAs are sparse enough that MBA gets zeroed across most cells while NEO never
is, so `P(NEO) = ρ_NEO / Σρ` collapses to **1.0 over ~75–82% of the band**. Every non-NEO there is
scored a certain NEO. Removing MBA cloning thinned MBA support ~5× (median support 4→1, 7→2), which
widened the zeroed region 75.1% → 82.4% — real, but a **modest worsening of a pre-existing bug**,
not the cause of it.

#### What this means

1. **The contamination gap is not a NEOMOD3 problem.** Fixing it means fixing the asymmetric support
   mask, which predates this work and equally afflicts the production maps.
2. ~~**The obvious fix to test:** stop exempting NEO from the support mask, or raise MBA's effective
   support in the tail.~~ **TESTED AND REFUTED — see §9.5.** Neither works: removing the mask
   entirely leaves the band's FP count unchanged (12,950 → 12,994), and un-exempting NEO destroys
   the classifier (AUC 0.4326, 99.3% of fast NEOs zeroed). The masking explanation below is
   correlation, not cause.
3. **Do not conclude "NEOMOD3 loses to digest2".** It wins near antisun and at high |Δlon|, loses at
   |lat| > 35° and Δlon 20–50°, and is far better on TNO/Trojan confusion. It is the |v| 0.25–0.5
   MBA band — an inherited defect — that decides the pooled number.


---

# §10. RUNBOOK — how to build the maps from scratch

Everything needed to reproduce `prob_maps_grid_neomod3_full/` (667 maps, NEOMOD3 clone-only, ±5 deg/day)
on a clean checkout. Written so it can be followed without re-reading §1–§9.

## §10.0 What a "map" actually is

One `.npz` per sky center = **151 arrays**. For each of 4 populations × 8 magnitude bins
(`14_16, 16_18, 18_20, mag20, mag21, mag22, mag23, mag24+`) it stores a 1001×1001 grid over
(`vlam`, `vbeta`) ∈ [−5, +5] deg/day at 0.01 step:

| array | meaning |
|---|---|
| `density_raw__<POP>__<BIN>` | kNN density of that population's clones, downweighted by `effective_factor` |
| `support_count__<POP>__<BIN>` | how many clones actually support each cell (drives the support mask) |
| `nearest_dist__<POP>__<BIN>` | distance to nearest clone (optional masking, **off** in production scoring) |
| `magcut_count__<POP>__<BIN>` | scalar bookkeeping — **not** a normalisation, see §8.6b |
| `x_grid`, `y_grid`, `grid_lat_deg`, `delta_lon_from_antisun_deg` | grid + center metadata |

Scoring is then `P(NEO) = ρ_NEO / Σ_pop ρ_pop`, bilinear-interpolated at the tracklet's
(`vlam`, `vbeta`) in its magnitude bin. Off-grid ⇒ 0 — which is exactly why the domain was widened
from ±2 to ±5.

## §10.1 Pipeline shape

```
 S3M .s3m files                      NEOMOD3 datacube (Nesvorny+2024)
        |                                        |
        | STAGE 0: n-body to epoch               | STEP 2: sample + project to sky
        v                                        v
 epoch_state_<EPOCH>.parquet          neomod3_projection_<EPOCH>.parquet   (97M clones, 7.9 GB)
 (14.4M objects, MBA/TNO/Trojan)                 |
        |                                        | STEP 3: HEALPix partition (nside=8)
        |                                        v
        |                             by_pixel/pix=NNN/*.parquet  (768 dirs, 10.2 GB)
        |                                        |
        +----------------+-----------------------+
                         | STEP 5: 667 array tasks, each reads ~100 pixels
                         v
              prob_maps_grid_neomod3_full/*.npz
```

**Division of labour:** NEO density comes **entirely** from NEOMOD3 (an external absolute model);
MBA/TNO/Trojans come from real S3M objects and are **not cloned** (`clone_factor=1`). That asymmetry
is the whole point of this module — see §6.1 for why it stays calibrated.

## §10.2 Prerequisites

| thing | path | note |
|---|---|---|
| python | `conda_prep/bin/python` | **never** `conda activate` — that was a past bug |
| slurm | `--account=astro` | **not** `--account=dirac` |
| module | `neomod/src/velocity_density_pipeline_neomod_clone_only.py` | the clone-only VDP; `_gmm.py` is the untouched original |
| S3M data | `neomod/S3Mdata/*.s3m` | `s3m_loader` resolves **CWD-relative** — run from `neomod/`, and set `VDP_LOADER=s3m` (default is *hybrid*, a different `t_0`) |

## §10.3 Step 1 — Stage-0 epoch-state cache (MBA/TNO/Trojans)

Provides n-body-propagated states + magnitudes at the map epoch so no map task ever propagates.
Full spec in `fixing_integrator.md` §10.5.

```bash
sbatch neomod/pipeline/slurm/stage0_epoch_state.sbatch      # array over 435 shards
conda_prep/bin/python neomod/pipeline/stage0_merge.py       # -> outputs/epoch_state_cache/epoch_state_<EPOCH>.parquet
```
Current file: `epoch_state_2027-08-25T000000.parquet` (14,380,436 objects; MBA 13,883,361).

> **Epoch is not free.** Stage-0 propagation is epoch-dependent (19.6 yr vs 21.2 yr differ). A cache
> built for one epoch **cannot** be reused at another. Epoch must match everywhere in §10.4–§10.5.

## §10.4 Step 2 — NEOMOD3 projection cache (NEO)

Sample orbits from the NEOMOD3 4-D datacube (H 52 × a 42 × e 25 × i 22) and project them to
observable (`ra`, `dec`, `vlam`, `vbeta`, `mag_app`) **once, all-sky**, at the map epoch.

```bash
sbatch neomod/pipeline/slurm/neomod3_projection_cache.sbatch          # 40 shards x 1M orbits
conda_prep/bin/python neomod/pipeline/neomod3_projection_cache.py merge
conda_prep/bin/python neomod/pipeline/neomod3_projection_cache.py stats
```
Produces **97,056,755** clones (7.88 GB) + `cache_metadata.json`.

Two settled decisions baked in here:
- **`NEOMOD3_H_RANGE = (15.0, 28.0)` — do NOT cut at H<25** (§6.2). An H<25 cut deletes **81.2%** of
  detectable |v|>2 NEOs; the fast movers are nearby and faint.
- **Absolute normalisation** (§6.1): NEOMOD3 weights are absolute object counts, so each clone is
  worth `w_abs = 0.11433` real objects ⇒ `effective_factor_NEO = 8.7467`. This lives in
  `cache_metadata.json` and is the **only** correct normalisation source.

## §10.5 Step 3 — HEALPix partition (makes the 667-map build affordable)

```bash
conda_prep/bin/python neomod/pipeline/neomod3_cache_healpix.py build
conda_prep/bin/python neomod/pipeline/neomod3_cache_healpix.py validate    # T1
conda_prep/bin/python neomod/pipeline/neomod3_cache_healpix.py validate2   # T2
```
Writes `outputs/neomod3_projection_cache/by_pixel/pix=NNN/` — 768 pixels, 10.22 GB, **33 s**.

> **The one thing that makes this fast:** the cache is in projection order (random sky positions), so
> each record batch would scatter rows across all 768 open file handles — tiny simultaneous appends
> that crawl at ~2 MB/s on GPFS (**45+ min**). Sorting by `pix` before writing turns it into 768
> sequential bulk writes: **6 s**. Do not remove the `tbl.sort_by([("pix","ascending")])`.

Read path: `hp.query_disc(nside=8, vec, radians(max_sep + max_pixrad), inclusive=True)` — the
`max_pixrad` (7.473°) margin is what guarantees no clone near a patch edge is missed.

## §10.6 Step 4 — validate BEFORE building 667 maps

`T1–T4` (§8.6, results §8.6b). Skipping these risks 108 GiB of silently wrong maps.

```bash
sbatch neomod/pipeline/slurm/hpx_T3b_T4.sbatch    # T3b + T4
```

## §10.7 Step 5 — build the grid

```bash
sbatch neomod/pipeline/slurm/neomod3_grid_full.sbatch     # 0-666%200 on ckpt-all
```
~25 min wall for 667 maps, 108 GiB. Full option rationale in §9.

**The three settings that silently corrupt everything if wrong:**

| flag | must be | if wrong |
|---|---|---|
| `--ref-obstime 2027-08-25T00:00:00` | matches both caches | **default is 2026-01-01** — builds 667 maps at the wrong epoch, no error raised |
| `--vdp-module velocity_density_pipeline_neomod_clone_only` | the clone-only module | silently builds old GMM-cloned maps |
| `NEOMOD3_CACHE_DIR=.../by_pixel` | set | falls back to monolithic: correct maps, ~6× the I/O |

Also: `VDP_LOADER=s3m`, `--mba-clone-factor 1`, `--velocity-grid-limit 5.0 --velocity-grid-step 0.01`,
and **omit `--overwrite`** so preempted tasks can be resumed.

## §10.8 Step 6 — accept the grid

```bash
conda_prep/bin/python neomod/pipeline/nm3_fullgrid_accept.py
```
Checks A–D of §9.1. **Read §9.2 before trusting a green file count** — a truncated `.npz` passes
`ls | wc -l`, passes Slurm's `COMPLETED`, and is *skipped* by a plain resubmit. Only an
open-and-read integrity scan catches it; rebuild those indices with `--overwrite`.

## §10.9 Changing things

| want to change | change | rebuild from |
|---|---|---|
| epoch | `--ref-obstime` + both caches | §10.3 (everything) |
| velocity domain | `--velocity-grid-limit` | §10.7 only |
| grid resolution | `--velocity-grid-step` | §10.7 only |
| H range / #clones | `NEOMOD3_H_RANGE`, `N_TOTAL` | §10.4 (re-partition + re-validate) |
| sky sampling | `--lon-step`, `--lat-base`, `--sun-exclusion` | §10.7 (+ regenerate the manifest) |

## §10.10 Wall-clock budget

| step | cost |
|---|---|
| Stage 0 (§10.3) | ~200 core-h ≈ 30 min wall — **one-time per epoch** |
| NEOMOD3 projection (§10.4) | 40 shards ≈ 1 h |
| HEALPix partition (§10.5) | **33 s** |
| T1–T4 (§10.6) | ~20 min |
| 667 maps (§10.7) | **~25 min** at %200 |
| acceptance (§10.8) | ~15 min |

With both caches already built, a full grid rebuild is **under an hour**.

### §9.5 SUPPORT-MASK FIX TEST — all three candidates FAIL, and the §9.4 mechanism is **REFUTED**

Script: `neomod/pipeline/nm3_support_mask_fix.py` (job 37995564, 4.6 min). The mask is applied at
**score time**, so all variants are pure re-scores of the same 667 maps — no rebuild.

| variant | AUC | best F1 | complet. | contam. | FP total | **FP in 0.25–0.5** | fast-NEO scored 0 | fast-NEO median |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline (NEO exempt) | 0.9288 | 0.8319 | 79.23% | 12.43% | 16,739 | 12,950 | 14.03% | 1.000 |
| **symmetric** (NEO masked too) | **0.4326** | 0.3871 | 27.20% | 32.90% | 19,841 | 14,525 | **99.28%** | 0.000 |
| **nomask** | **0.9295** | **0.8334** | 78.83% | **11.60%** | **15,397** | 12,994 | 14.03% | 0.651 |
| mba_relax (MBA needs support>0) | 0.9288 | 0.8319 | 79.23% | 12.43% | 16,739 | 12,950 | 14.03% | 1.000 |
| digest2 | 0.9430 | 0.8489 | 77.79% | **6.60%** | 8,173 | **1,502** | 0.01% | 1.000 |

#### ❌ The §9.4 mechanism was WRONG. Retracted.

§9.4 concluded the 0.25–0.5 false positives were caused by the support mask zeroing ρ_MBA while NEO
stayed exempt. **`nomask` disproves it**: turning the mask *completely off* changes the band's FP
count by **+44 (12,950 → 12,994)** — essentially nothing. If masking were the cause, removing it
would have collapsed that number. It does not. The correlation in §9.4 (75.1% → 82.4% of cells
zeroed) was real but **not causal**.

What is actually happening: even with ρ_MBA fully restored, **ρ_NEO still dominates the sum in that
velocity band**, so `P(NEO) = ρ_NEO/Σρ` stays high. The bled/unmasked MBA density there is tiny
compared with the absolutely-normalised, smoothed NEO density. This is a **density/normalisation
property, not a masking bug** — which also explains why production VDP (MBA cf=5) has an almost
identical 24.01% band FP rate. **Five times more MBA clones does not fix it either.**

`mba_relax` returning results **identical to baseline, to every digit**, pins this down further:
lowering MBA's threshold from `support ≥ 1` to `support > 0` changes nothing, so MBA support in
those cells is **exactly 0**, not fractional. There are no MBA clones there at all to un-mask.

#### ❌ "Stop exempting NEO" is catastrophic — the exemption is load-bearing

`symmetric` destroys the classifier: AUC **0.9288 → 0.4326** (worse than chance in ranking terms),
and **99.28%** of the 7,727 fast NEOs go to zero — the entire ±5 grid result, erased. The reason is
structural: NEO is the **smoothed** population, and smoothing deliberately fills low-support cells.
Masking by *raw* support then deletes exactly what smoothing just created. `_support_mask_skip` is
not an oversight to be removed; it is what makes smoothing coherent. **Do not "fix" it.**

The guardrail column earned its place here — on pooled AUC alone `symmetric` merely looks bad, but
the fast-NEO column shows it annihilates the specific result this work exists to deliver.

#### ~ `nomask` is a small genuine improvement, with a trade-off

Dropping the support mask entirely: AUC 0.9288 → 0.9295, F1 0.8319 → 0.8334, contamination
**12.43% → 11.60%**, 1,342 fewer false positives — while fast-NEO recovery is **unchanged** (14.03%
scored zero either way). The cost is confidence, not coverage: fast-NEO median P falls **1.000 →
0.651**. Real but marginal; it changes production scoring behaviour, so **flagged for a decision, not
adopted**.

#### Where this leaves the digest2 gap

Unclosed, and none of these levers touch it. digest2 takes **1,502** FPs in the 0.25–0.5 band against
our ~12,950 — a ~8.6× difference that survives every masking variant *and* a 5× change in MBA
cloning. The remaining hypotheses, in the order worth testing:
1. **The NEO density is genuinely too high at 0.25–0.5 deg/day** — an absolute-normalisation
   (§6.1) question, not a masking one. Check ∫ρ_NEO restricted to that band against the expected
   NEO count there.
2. **kNN bleed from the dense NEO core** into the band (the `vdp_bleed_diagnostic_bundle` failure
   mode, but for NEO rather than MBA).
3. **digest2 uses information VDP does not** — arc curvature, not just mean rate. If so the gap is
   not closable by re-weighting densities, and the combined VDP × digest2 score becomes the
   practical answer.

### §9.6 THE 0.25–0.5 BAND EXPLAINED — a real physical overlap, not a density bug

Scripts: `nm3_band_diagnosis.py`, `nm3_band_discrimination.py`. Both remaining §9.5 hypotheses are
**refuted**; the band is genuinely ambiguous, and the practical fix is a combined score.

#### The band is ~50/50 NEO/MBA in truth

| | share of band |
|---|---:|
| NEO | **47.8%** |
| MBA | **47.4%** |
| other | 4.9% |

**32.5% of all NEOs and 10.3% of all MBAs move at 0.25–0.5 deg/day.** NEO median speed is 0.466,
MBA median 0.171 with a 90th percentile of 0.251 — so this band is exactly where the slow tail of
the NEO distribution meets the fast tail of the MBA distribution. A classifier that sees only
(vlam, vbeta, mag) is looking at a genuinely ~50/50 mixture. Most of the "24.5% FP rate" is
**irreducible ambiguity, not a defect.**

#### ❌ Hypothesis 1 (ρ_NEO too high) — REFUTED, the opposite is true

Predicted ÷ true NEO fraction, by speed:

| \|v\| | true NEO frac | VDP predicted/true | digest2 predicted/true |
|---|---:|---:|---:|
| <0.15 | 0.058 | **0.089** | **1.864** |
| 0.15–0.20 | 0.072 | 0.171 | 1.358 |
| 0.20–0.25 | 0.073 | 0.280 | 1.263 |
| 0.25–0.30 | 0.225 | 0.565 | 0.848 |
| 0.35–0.40 | 0.688 | 0.917 | 0.902 |
| 0.4–0.5 | 0.816 | 0.953 | 0.905 |
| >0.7 | 1.000 | 0.999 | 0.970 |

VDP **never over-predicts** — ratio ≤ 1 everywhere, and it is *conservative* at low speed (0.089 at
|v|<0.15). It is **digest2** that over-predicts there (1.86×). The NEO density is not too high; if
anything the VDP is under-confident. Re-normalising ρ_NEO downward would make things worse.

#### ❌ Hypothesis 2 (kNN bleed from the NEO core) — REFUTED

`density_raw` is already downweighted, so ∫ρ over a region = expected object count, and
(Σ support_count) × w_abs = the count implied by the clones actually present. Their ratio is ~1 when
density is earned and ≫1 when manufactured:

| population | core \|v\|<0.25 | band 0.25–0.5 | fast \|v\|>2 |
|---|---:|---:|---:|
| NEO | 1.09–1.13 | **1.09–1.10** | 1.14–1.17 |
| MBA | 1.11–1.12 | **1.11–1.13** | ∞ (0 clones) |

The ratio is **flat at ~1.1 across every region and both populations** — a uniform smoothing/
bookkeeping offset, not region-dependent inflation. There is no bleed into the band; the NEO density
there is supported by NEO clones that genuinely are there. (MBA at |v|>2 is the one genuine bleed case:
zero clones, tiny non-zero density — but that region is irrelevant, MBAs never move that fast.)

#### ⚠️ The "8.6× FP gap in the band" was partly a threshold artefact

Restricting to the band and letting each classifier pick its **own** operating point:

| slice | classifier | AUC | best F1 | contamination |
|---|---|---:|---:|---:|
| core \|v\|<0.25 | VDP | 0.7709 | 0.4072 | 50.9% |
| core \|v\|<0.25 | digest2 | **0.8250** | **0.5345** | **30.4%** |
| **band 0.25–0.5** | VDP | 0.9107 | 0.8353 | 17.2% |
| **band 0.25–0.5** | digest2 | **0.9286** | 0.8432 | **12.1%** |
| fast \|v\|>0.5 | VDP | 0.8851 | 0.9917 | 1.6% |
| fast \|v\|>0.5 | digest2 | **0.9836** | 0.9944 | **0.7%** |

Within the band the two are **close** (AUC 0.911 vs 0.929) — not the 8.6× disaster implied by
comparing FP counts at *globally* optimised thresholds. VDP's global best-F1 threshold is **0.0114**,
very low, so at the global operating point it sweeps in band MBAs. §9.4/§9.5 overstated this as a
"band failure"; a large part of it is **where the global threshold lands**. digest2's real remaining
edge is broader — it is ahead in *every* velocity slice, most of all in the slow core (0.825 vs
0.771), which is consistent with it using arc curvature rather than mean rate alone.

#### ✅ Combined score — beats both, and costs nothing

| score | AUC | best F1 | complet. | contam. | FP in band | fast NEOs scored 0 |
|---|---:|---:|---:|---:|---:|---:|
| NEOMOD3 VDP | 0.9288 | 0.8319 | 79.2% | 12.43% | 12,950 | 14.03% |
| digest2 | 0.9430 | 0.8489 | 77.8% | 6.60% | 1,502 | 0.01% |
| product / geometric | 0.9518 | 0.8577 | **85.1%** | 13.53% | 12,403 | 14.04% |
| **rank-average** | **0.9576** | **0.8625** | 81.1% | **7.95%** | **5,881** | **0.00%** |

**The rank-average of the two scores beats both classifiers on AUC and F1**, nearly matches digest2
on contamination (7.95% vs 6.60%), less than halves the band false positives (12,950 → 5,881), and
**recovers 100% of the 7,727 fast NEOs** — strictly better than VDP's 14.03% miss rate.

Why rank-average and not the product: a product inherits VDP's zeros (0 × anything = 0), so it
keeps the 14% fast-NEO loss. Rank-averaging never zeroes — an off-grid VDP score still has a rank,
and digest2's rank carries the object. That is why it is the only variant with **0.00%** in the last
column.

**Conclusion.** The VDP and digest2 fail in different places, so combining them is not a workaround
but the correct use of two complementary classifiers.

> ⚠️ **The "adopt the rank-average" recommendation originally written here was WRONG and is
> withdrawn — see §9.7.** Rank-averaging destroys the probabilistic meaning of the score: it sums to
> **2.19×** the true NEO count. A calibrated logit-combination gets the same AUC, a *better* F1, and
> stays a real probability. Use that instead.

### §9.7 KEEPING THE PROBABILITY — calibrated combination beats the rank-average outright

Script: `nm3_calibrated_combine.py`. Fit on **half the sky centers**, evaluated on the held-out half
(324,413 tracklets / 309 centers), so nothing below is fitted on its own test data.

The objection that prompted this: rank-averaging is no longer a probability, and a probability is
precisely what distinguishes the VDP (`P = ρ_NEO/Σρ`, a real posterior) from digest2 (a heuristic
score). Optimising AUC/F1 while discarding calibration is optimising the wrong thing.

| score | AUC | best F1 | contam. | Brier ↓ | ECE ↓ | **Σp / true count** | fast NEO = 0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| VDP alone | 0.9286 | 0.8312 | 9.00% | 0.0668 | 0.0631 | **0.763** | 13.7% |
| digest2 alone | 0.9423 | 0.8454 | 7.07% | 0.0561 | 0.0296 | 1.057 | 0% |
| rank-average (§9.6) | 0.9573 | 0.8621 | 8.12% | — | — | **2.191** | 0% |
| VDP isotonic-recal | 0.9294 | 0.8312 | 9.00% | 0.0578 | 0.0059 | 1.023 | 0% |
| **logit-combine + isotonic** | **0.9576** | **0.8735** | 7.85% | **0.0439** | **0.0015** | **1.002** | **0%** |

`Σp / true count` is the test that matters for population work: sum the scores over a sample and you
should recover the number of NEOs in it.

#### The rank-average is unusable as a probability — confirmed

**Σp = 2.19× the truth.** It is a ranking statistic whose values are uniform by construction, so its
sum carries no physical meaning at all. Ranking-only applications are fine; anything that sums,
weights, or thresholds on a physical probability is not.

#### ❌ Correction: "VDP is well calibrated" (§9.6) was WRONG

§9.6 read VDP's *ratio ≤ 1 at every speed* as calibration and called it "conservative". That was a
misreading — **systematic under-prediction is a calibration failure**, just in the safe direction.
Measured: **Σp = 0.763**, i.e. VDP under-counts NEOs by **24%**, with ECE 0.063 (20× worse than the
combination). Its reliability curve is badly off in the low bins:

| VDP predicts | actually observed |
|---:|---:|
| 0.002 | **0.068** |
| 0.189 | **0.528** |
| 0.406 | 0.605 |
| 0.998 | 0.974 |

So the VDP was never the calibrated one — digest2 (Σp = 1.057) is closer out of the box. This does
not change any ranking result in §9.3–§9.6 (AUC/F1 are threshold-free), but it does mean **raw VDP
`P_NEO` must not be summed for population estimates.**

#### ✅ The calibrated combination dominates on every axis

`logit-combine` = logistic regression on `[logit P_vdp, logit P_d2]`, then isotonic recalibration —
both fitted on the training centers only. Fitted weights **[0.233 VDP, 0.436 digest2]**: digest2 gets
about twice the weight, but both contribute.

Versus the rank-average it is **equal on AUC (0.9576), better on F1 (0.8735 vs 0.8621)**, and it is
an actual probability: ECE **0.0015**, Σp **1.002**. Its reliability curve is essentially exact:

| predicted | observed |
|---:|---:|
| 0.031 | 0.031 |
| 0.162 | 0.164 |
| 0.408 | 0.391 |
| 0.614 | 0.600 |
| 0.828 | 0.825 |
| 0.997 | 0.997 |

There is no trade-off to argue about — it wins on ranking *and* calibration, and still recovers 100%
of the fast NEOs.

#### Also worth knowing: recalibrating the VDP alone is nearly free

Isotonic recalibration of VDP by itself takes ECE **0.0631 → 0.0059** and Σp **0.763 → 1.023**,
without digest2 and without touching a map. AUC is unchanged (isotonic is monotonic). **If the VDP
must stand alone as a physical probability, this step is required.**

#### Recommendation (replaces §9.6) — ⚠️ ITSELF WITHDRAWN, see §9.8

1. **Do not adopt the rank-average.**
2. For a standalone VDP probability: **apply isotonic recalibration**, fitted on held-out centers.
3. For best overall performance: **logit-combine + isotonic**, which is both the best classifier and
   a genuine probability.
4. Fit any calibration on **held-out centers** — sky position is the strongest confounder (§9.4), so
   a random row split would leak.

### §9.8 DECISION — do NOT adopt the fitted combination. §9.6 and §9.7 recommendations are both withdrawn.

§9.6 recommended a rank-average; §9.7 withdrew that and recommended `logit-combine + isotonic`.
**Both are now withdrawn.** The error was consistent: optimising AUC/F1 without asking what the score
is *for*. Recorded here because the reasoning matters more than the numbers.

#### 1. It needs truth labels, which exist only in simulation

Logistic regression and isotonic regression are **supervised**. They are fitted against
`is_neo`, which we have only because Sorcha generated it. On real LSST data there is nothing to fit.
Adopting them means the "probability" is calibrated to **simulation truth**, so every
simulation↔reality mismatch — NEOMOD3 vs the real NEO population, Sorcha's detection model, linking
efficiency — propagates silently into a number presented as a physical probability.

The VDP's current Σp = 0.763 under-count is an **honest defect of the density model**. A fitted
correction does not fix it; it conceals it behind assumptions that cannot be checked without labels.

#### 2. It destroys the independence that makes the digest2 comparison meaningful

If `P_NEO` is a function of `P_d2`, the VDP can no longer be used as an independent check on digest2.
The whole comparison this work rests on stops meaning anything.

#### 3. Per-object harm that the aggregate metrics hide (test half, 74,037 true NEOs)

| threshold | true NEOs **lost** | true NEOs gained | net |
|---|---:|---:|---:|
| P ≥ 0.5 | 282 | 7,201 | **+6,919** |
| P ≥ 0.9 | **3,840** | 3,219 | **−621** |

- **169 true NEOs drop by >0.30** (median 0.910 → 0.432); their median |v| is 0.302 — the overlap band.
- **98.6% of non-NEOs have their scores RAISED.**
- At P ≥ 0.5 the combination **adds** 1,437 false positives and removes only 451.

It improves AUC/F1 while making things *worse* at P ≥ 0.5 and demoting hundreds of genuine NEOs. The
aggregate numbers in §9.7 concealed all of this.

#### 4. An unresolved free parameter

`logit(0) = −∞`, so VDP zeros were clipped at `1e-6` (→ −13.816 in log-odds). At `1e-4` the same
object scores ~0.25 instead of ~0.107. Since ~14% of fast NEOs have VDP exactly 0, this arbitrary
epsilon materially controls their fate. It was never tuned.

#### Where this leaves things

**Legitimate use** of a fitted combination: candidate **triage** — ranking follow-up targets, where
only ordering matters and digest2 is being run anyway.

**Illegitimate**: anything claiming a physical probability, summing scores for population estimates,
or comparing VDP against digest2.

**The one genuinely open defect** is Σp = 0.763 — the VDP under-counts NEOs by 24%. It belongs to the
density model (candidates: the absolute normalisation §6.1, the magnitude binning, the support mask)
and deserves a **physical diagnosis, not a fitted patch**.

**Everything from §9.3 onward that framed "beat digest2" as the goal was benchmark-chasing.** The
objective set for this work — remove MBA/TNO/Trojan cloning, take NEOs from NEOMOD3, stop losing fast
NEOs off the grid edge — was met and validated in §8.6b, §9.1 and §9.3.

### §9.9 DENSITY AUDIT — the under-count is NOT a normalisation constant. NEOMOD3 and S3M describe different universes.

Script: `nm3_density_audit.py`. Output: `outputs/neomod3_fullgrid/density_audit.csv`. Four centers ×
8 magnitude bins. No truth labels, no simulation-transfer assumption — just
∫ρ dA (what the map says) vs clones×weight (what the clones imply) vs objects counted in the
Stage-0 n-body cache (what is actually there).

#### The control passes — the method is sound

MBA density is built from real S3M objects, so ∫ρ must reproduce the true count:

| magbin | ∫ρ | clones×w | **true count** | ∫ρ / true |
|---|---:|---:|---:|---:|
| 16_18 | 3,067 | 2,623 | **2,622** | 1.183 |
| mag21 | 52,270 | 46,270 | **46,260** | 1.132 |
| mag23 | 163,300 | 146,300 | **146,300** | 1.122 |
| mag24+ | 291,200 | 254,600 | **254,600** | 1.144 |

`clones×w` matches `true count` **to 4 significant figures** in every bin — as it must, since MBA is
no longer cloned, so its "clones" *are* the objects. This confirms the audit method works.

It also isolates a **real, separate defect**: ∫ρ overshoots by a **uniform ~13%** in every bin and
both populations. That is the kNN density estimator (plus smoothing) integrating to ~1.13× the
count it was built from. Uniform, so it cancels in `P = ρ_NEO/Σρ` — but it means **∫ρ is not
directly usable as an object count without dividing by ~1.13.**

#### NEO fails, and the failure grows with faintness

| magbin | clones×w | **true (S3M)** | **∫ρ / true** |
|---|---:|---:|---:|
| mag20 | 49.2 | 44.0 | 1.34 |
| mag21 | 119 | 80 | **1.70** |
| mag22 | 324 | 168 | **2.12** |
| mag23 | 868 | 324 | **2.97** |
| mag24+ | 2,285 | 614 | **4.12** |
| **total** | | **5,036** | **3.29** |

**NEOMOD3 puts 3.3× more NEOs in the patch than S3M has, rising to 4.1× in the faintest bin.**
`A/B` stays ~1.12 throughout, so this is **not** the estimator and **not** the support mask — it is
the clone population itself.

#### This is §6.2 showing up as a measurement, and it is expected

S3M covers **H < 25**; NEOMOD3 covers **H 15–28**. §6.2 deliberately kept the full range (an H<25 cut
would have deleted 81.2% of detectable |v|>2 NEOs). At faint *apparent* magnitude the H>25 objects
dominate — which is exactly the 0.9 → 4.1 gradient above. The findings notebook already said ~97% of
what NEOMOD3 describes is a population S3M does not contain; this is that statement, measured in the
maps.

**So the map is not broken. The map and the test set disagree about how many faint NEOs exist.**

#### What this means for the Σp = 0.763 under-count

It is **not** a single normalisation constant, and the sign is instructive: ρ_NEO is *too high*
relative to S3M by 3.3×, yet Σp comes out *low*. Both are true because `P = ρ_NEO/Σρ` is bounded:
where NEO dominates, P saturates at 1 and the excess cannot show up; the Σp deficit is generated in
the slow, MBA-dominated regime where P is small (VDP predicts 0.002 where S3M truth is 0.068) and
where most objects live. So the 24% deficit is **not** correctable by rescaling ρ_NEO — rescaling it
up would push the already-saturated regions nowhere and worsen the faint-end excess.

#### ⚠️ Consequence: VDP calibration CANNOT be validated on S3M tracklets at faint magnitudes

The prior (NEOMOD3) and the test data (S3M) describe different faint-NEO populations. Any
calibration fitted against S3M truth — including §9.7's isotonic — would be **forcing NEOMOD3's
prior back onto S3M's population, undoing the §6.2 decision through the back door**, silently and
only at faint magnitudes. This is an independent reason not to adopt §9.7, and it does not depend on
the label-availability argument in §9.8.

#### Where that leaves the calibration question

1. **∫ρ ≈ 1.13 × count** is a genuine estimator bias, uniform, and the one cleanly fixable item.
2. **The faint-end 3–4× is a population-model disagreement, not a bug.** It can only be adjudicated
   against a test set that shares NEOMOD3's H range — S3M cannot settle it.
3. **Σp = 0.763 is therefore partly an artefact of evaluating a NEOMOD3-prior classifier on
   S3M-generated tracklets** and should not be "fixed" until (1) and (2) are separated.
4. Only 0.52% of true NEOs move faster than the ±5 grid edge, so off-grid loss is negligible and is
   **not** a contributor.

### §9.10 The ~13% bias LOCATED (classic kNN bias, one-constant fix) and the faint excess DECOMPOSED

Script: `nm3_estimator_bias.py`.

#### The ~13% overshoot is the standard kNN normalisation bias

Measured against a synthetic field of **known** uniform density:

| k | measured ρ̂/ρ_true | k/(k−1) |
|---:|---:|---:|
| 5 | 1.26–1.33 | 1.250 |
| **10 (`DEFAULT_K_MAP`)** | **1.114–1.144** | **1.111** |
| 20 | 1.031–1.046 | 1.053 |
| 40 | 1.044–1.049 | 1.026 |

**Origin.** `estimate_density_full_posterior_2d` returns the posterior *mean* of `n0 = 1/(pi d0^2)`
— a mean of an inverse square. By Jensen `E[1/d0²] > 1/E[d0]²`, so averaging in density space is
biased high by construction. This is the textbook Loftsgaarden–Quesenberry result: `k/(N·V_k)` is
biased; the unbiased estimator is `(k−1)/(N·V_k)`.

**Fix:** multiply density by **(k−1)/k = 0.9** at k=10. (Raising k also works but over-smooths, which
would blur the sharp NEO/MBA velocity structure the classifier depends on — do not.)

**Scope — this changes NO classification result.** The bias is uniform across populations and
magnitude bins, so it cancels exactly in `P = ρ_NEO/Σρ`. Every AUC/F1/ROC number in this document is
unaffected. It matters only where ∫ρ is used as an **object count**, i.e. the §9.9 audit and any
population inference. Fix it for correctness of ρ as a density, not to improve scoring.

#### The faint-end excess is ~80% the H>25 tail, ~20% a real normalisation difference

| magbin | NEOMOD3 / S3M | NEOMOD3 restricted to **H<25** / S3M |
|---|---:|---:|
| mag21 | 1.91–1.93 | 1.26–1.27 |
| mag22 | 2.26–2.53 | 1.15–1.36 |
| mag23 | 3.22–3.77 | 1.28–1.40 |
| mag24+ | **4.59–4.80** | **1.23–1.28** |

Restricting to S3M's own H range collapses 4.8× → ~1.25×. So:

1. **The dominant term is the H>25 population S3M does not contain** — the §6.2 decision, measured.
2. **A residual ~1.3× survives at matched H.** That is a genuine normalisation difference between
   NEOMOD3 and S3M for the *same* population, consistent with the 1.30× recorded when `w_abs` was
   validated (§6.1). **This part IS answerable with S3M** and is the more tractable open question.

#### Why S3M cannot settle term (1) — the structural argument

Scored tracklets are Sorcha detections of **S3M** objects, so `is_neo` means "came from an S3M NEO".
**The test set contains no H>25 NEOs at all** — not undetected, *never simulated*. At mag 24–25 the
map predicts ~2,285 NEOs in a patch; the test set holds 614. Either:

- NEOMOD3 is right, and the test set is missing ~1,671 NEOs LSST would actually detect; or
- NEOMOD3 over-predicts by ~4×.

**The measurement is identical under both.** Testing a vehicle-counter on a dataset containing only
cars will always look like over-counting, and cannot distinguish a wrong model from an incomplete
dataset. Fitting a correction to that data does not resolve the ambiguity — it silently assumes the
second answer and deletes the faint population §6.2 deliberately kept. **This is the deepest reason
§9.7 must not be adopted.**

**What would settle it:** run Sorcha on a NEOMOD3-derived synthetic catalogue spanning H to 28
through the same detection + linking pipeline, so truth contains faint NEOs. Then predicted-vs-actual
counts are meaningful. Concrete work, not a fitting exercise.

#### Ranked open items

| # | item | tractable with what we have? |
|---|---|---|
| 1 | apply `(k−1)/k` density correction | **yes** — one constant, no rebuild needed for scoring |
| 2 | residual ~1.3× NEOMOD3 vs S3M at matched H | **yes** — compare the two population models directly |
| 3 | faint-end (H>25) validity | **no** — needs a NEOMOD3-based Sorcha run |

---

# §11. NEXT SESSION — the benchmark is wrong, and one estimator fix to apply

## §11.1 ⚠️ THE BENCHMARK ITSELF IS THE PROBLEM (advisor, 2026-08-01)

**We have been scoring S3M-generated tracklets with maps whose NEO prior is NEOMOD3.** The test set
and the model describe different populations, so every evaluation from §9.3 onward rests on a
benchmark that cannot answer the question asked of it. §9.10 reached this structurally; the advisor's
instruction states it directly: **pure S3M tracklets are not a sufficient benchmark.**

**What must be rebuilt:** a benchmark whose objects come from the **same populations that build the
maps** —

| population in the benchmark catalogue | source | why |
|---|---|---|
| **NEO** | **NEOMOD3 samples, H 15–28** | this is what the maps' NEO density is built from |
| MBA / TNO / Trojans | S3M (unchanged) | this is what the maps' non-NEO densities are built from |

That composition is exactly the map composition, which is the point: truth and prior then describe
the same universe, and predicted-vs-actual counts become meaningful for the first time.

**What this does and does not invalidate:**

- **Unaffected** — the map build, T1–T4, the HEALPix work, §8.6b, §9.1, §9.2, §9.10's estimator
  result. None of it depends on the tracklet set.
- **Provisional, must be re-measured on the new benchmark** — every number in §9.3–§9.9: the
  digest2 comparison, the ROC/AUC/F1 tables, the contamination analysis, the |v| 0.25–0.5 band
  result, Σp = 0.763, and the whole calibration discussion. The faint-end NEOs missing from S3M are
  precisely the ones driving those numbers.
- **Still valid and strengthened** — the decision NOT to adopt a fitted combination (§9.8). Fitting
  a correction against a benchmark now known to be the wrong population would have baked that error
  in permanently.

**Pieces that already exist:**

| need | where |
|---|---|
| NEOMOD3 orbit sampling from the datacube | `neomod3_projection_cache.py` (already samples a/e/i/H) |
| Sorcha input format + photometry convention | `prep_s3m_sorcha_inputs.py` → `inputs/s3m_sorcha_{orbits,phys}.csv` |
| orbit CSV schema | `ObjID,FORMAT,q,e,inc,node,argPeri,t_p_MJD_TDB,epochMJD_TDB` (COM) |
| phys CSV schema | `ObjID,H_r,GS,u-r,g-r,i-r,z-r,y-r`; `H_r = H − (Johnson V − LSST r)`, `GS = 0.15`, colours from `CDS_colors.parquet` seed 42 |
| full run pipeline | `neomod/pipeline/slurm/s3m_linking/` (6 stages, CASE-parametrised) |

Ask about the file Jake had — it may already do the NEOMOD3 → Sorcha-input step and save writing it.

**Sketch:** sample NEOMOD3 NEOs with counts set by its absolute normalisation (`w_abs`), convert to
COM elements + phys params in the `prep_s3m_sorcha_inputs.py` convention, concatenate with S3M's
MBA/TNO/Trojans, then run the existing 6-stage linking pipeline as a new CASE. Re-run §9.3's scoring
against that.

## §11.2 FIX TO APPLY: the kNN density normalisation bias

Diagnosed in §9.10, not yet applied. Recorded here so it is not lost.

Measured against a synthetic field of **known** uniform density:

| k | measured ρ̂/ρ_true | k/(k−1) |
|---:|---:|---:|
| 5 | 1.26–1.33 | 1.250 |
| **10 (`DEFAULT_K_MAP`)** | **1.114–1.144** | **1.111** |
| 20 | 1.031–1.046 | 1.053 |
| 40 | 1.044–1.049 | 1.026 |

It tracks **k/(k−1)**, and at the map's `k=10` that is 1.111 — the ~1.13 overshoot seen in the maps
(§9.9), where MBA's `clones×w` matched the true object count to 4 significant figures but ∫ρ came out
13% high in every single bin.

**Where it comes from.** `estimate_density_full_posterior_2d` returns the posterior *mean* of
`n0 = 1/(pi*d0^2)` — a mean of an inverse square. By Jensen's inequality `E[1/d0²] > 1/E[d0]²`, so
averaging in density space is biased high **by construction**. This is the textbook
Loftsgaarden–Quesenberry result: `ρ̂ = k/(N·V_k)` is biased, and the unbiased form is
`(k−1)/(N·V_k)`.

**The fix is one constant** — multiply by `(k−1)/k = 0.9` at k=10, in
`velocity_density_pipeline_neomod_clone_only.py:1298`:

```python
n0_mean = np.trapezoid(n0_grid * p_d0, d0_grid)
n0_mean *= (k - 1.0) / k      # Loftsgaarden-Quesenberry: E[1/d0^2] is biased high by k/(k-1)
```

Raising `k` also reduces the bias but **over-smooths**, blurring exactly the sharp NEO/MBA velocity
structure the classifier depends on — do not do that instead.

**It changes NO classification result.** The bias is uniform across populations and magnitude bins,
so it cancels exactly in `P = ρ_NEO/Σρ`. Every AUC/F1/ROC number in this document is unaffected. It
matters only where ∫ρ is used as an **object count** — the §9.9 audit and any population inference.
Fix it for correctness of ρ as a density, not to improve scoring.

**Cost:** requires a full 667-map rebuild (~25 min, §10.7) for the stored densities to be correct.
Since it changes no score, **batch it with the next rebuild rather than doing it on its own.**

## §11.3 NEOMOD3 BENCHMARK — BUILT (2026-08-01)

`neomod/pipeline/gen_benchmark_tracklets_neomod3.py` + `slurm/benchmark_neomod3.sbatch`
→ `outputs/benchmark_tracklets_neomod3/tracklets_benchmark_neomod3.parquet` (3,220,714 rows, 437 MB).

Composition mirrors the maps: **NEO from an independent NEOMOD3 draw**, MBA/TNO/Trojans from the
Stage-0 n-body epoch cache. Schema is **column-identical to v3**, so every downstream scorer runs
unchanged.

| population | rows | share |
|---|---:|---:|
| MBA | 3,075,543 | 95.493% |
| Trojans | 82,355 | 2.557% |
| TNO | 37,818 | 1.174% |
| **NEO** | **24,998** | **0.776%** |

### Four deliberate differences from the S3M benchmark

**1. NEO comes from a fresh, INDEPENDENT NEOMOD3 draw** — 30M orbits at seed 20270825, *not* the
map cache's seed 42. NEOMOD3 is a distribution, not a catalogue, so this is a second independent
realisation of the same population. This is not fussiness: maps are evaluated on a fixed 1001×1001
velocity grid and a clone typically **dominates its own cell**, so re-using cache clones as test
objects would have the map recognising its own training data — for the one population under test.
Self-consistency bonus: the independent 10M and 30M draws predict 25,115 and 24,998 in-grid NEOs
(0.5% apart), which is an end-to-end check of the absolute normalisation.

**2. Magnitude cut 14 ≤ mag_app < 25** (v3 had none). v3 is only **32.7%** scoreable by the maps and
its **NEO just 9%** (S3M NEO median mag 28.3). Without the cut the benchmark is mostly objects the
maps structurally cannot score.

**3. TRUE absolute counts — no downscaling** (`TARGET_TOTAL = None`, scale = 1.000000). Each sky
direction carries the number of objects actually there at this epoch. Verified:

| direction | benchmark rows | true objects | ratio | NEO % |
|---|---:|---:|---:|---:|
| 0–20° | 570,083 | 570,119 | **1.000** | 1.058 |
| 20–40° | 501,192 | 501,210 | **1.000** | 0.930 |
| 40–70° | 663,807 | 663,807 | **1.000** | 0.806 |
| 70–110° | 814,110 | 814,067 | **1.000** | 0.649 |
| 110–141° | 671,522 | 671,510 | **1.000** | 0.547 |

The NEO fraction falling monotonically from antisun outward (1.058% → 0.547%) is the expected
physical signature, not something imposed.

**4. Population mix from ABSOLUTE expected counts** — NEO from NEOMOD3's own normalisation
(11,432,918 NEOs at H 15–28, scaled by survival through the mag + sun-exclusion cuts), non-NEO from
real S3M counts.

### ⚠️ Guard against the v1/v3 ratio bug

v1 capped MBA at 200k and Trojan at 100k while leaving NEO and TNO uncapped — MBA suppressed **~69×**
and TNO's share inflated. **v3 kept a milder form of it**: MBA capped at 650k of 13.88M with no
magnitude cut, which is why its NEO share reads 1.76%.

Three guards now:
1. **One common scale factor for all four populations**, never per-population.
2. **Share-integrity check**: worst error **0.0000 percentage points**.
3. The silent-clamp branch is a **fatal error**, not a quiet cap — that branch *was* the v1
   mechanism, so it now refuses to run rather than distorting the mix.

### The honest S3M vs NEOMOD3 comparison

Under **identical** cuts (in-grid, 14 ≤ mag < 25, same epoch), NEO count differs only by model:

| | NEO in-grid | share |
|---|---:|---:|
| S3M | 13,225 | 0.412% |
| **NEOMOD3** | **25,115** | **0.780%** |
| ratio | **1.90×** | |

**Do not compare either to v3's 1.76%** — v3 had no magnitude cut and a 21× MBA suppression, so its
shares are not commensurate with either number.

### Ready for scoring

667/667 map cells populated; 5,050 NEOs with |v| > 2 (the fast movers the ±5 grid exists for).
Next: score with VDP + digest2 and re-measure §9.3–§9.9 on a benchmark whose populations finally
match the maps. If digest2 cost forces a subsample, do it in the **scoring** stage with a uniform
rate — never in this file, which is the truth set.

## §11.4 RUNBOOK — rebuilding the NEOMOD3 benchmark from scratch

Everything needed to regenerate `tracklets_benchmark_neomod3.parquet`. Written to be followed without
re-reading §11.1–§11.3.

### §11.4.1 What a benchmark tracklet is

**Not** a survey simulation. It is the *geometric ground-truth* set: every object that is really
there at one instant, propagated, assigned to a map cell, and turned into a synthetic 2-detection
tracklet. No detection efficiency, no linking, no losses. Sorcha tracklets (the other test set) add
all of that; this one isolates the classifier from survey effects.

Per object: position at `t0` and at `t0 + 30 min`, extrapolated linearly from the n-body rates.

### §11.4.2 Files

| file | role |
|---|---|
| `neomod/pipeline/gen_benchmark_tracklets_neomod3.py` | the builder — `neo-shard` and `build` subcommands |
| `neomod/pipeline/slurm/benchmark_neomod3.sbatch` | array job for the NEO sampling shards |
| `neomod/src/neomod3_sampler.py` | `sample_neomod3_orbits(n, epoch, rng)` — draws from the datacube |
| `neomod/src/velocity_density_pipeline_gmm.py` | `build_visible_subset_dataframe` (propagate + observables), `load_s3m_population` (observer geometry only) |
| `outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet` | MBA/TNO/Trojans, already n-body propagated |
| `outputs/neomod3_projection_cache/cache_metadata.json` | supplies `total_weight_absolute_NEO_count` = 11,432,918 |
| `neomod/pipeline/gen_benchmark_tracklets_s3m.py` | the S3M original, for reference/comparison |

**Outputs** (`outputs/benchmark_tracklets_neomod3/`): `neo_shards/neo_shard_NNN.parquet` (19 MB
total), `tracklets_{NEO,MBA,TNO,Trojans}.parquet`, `tracklets_benchmark_neomod3.parquet` (437 MB),
`benchmark_metadata.json`.

### §11.4.3 Commands

```bash
# 1. NEO: independent NEOMOD3 draw, 30 shards x 1M orbits  (~1.5 min/shard, 4 cpu, ~1.2 GB RSS)
sbatch --array=0-9  neomod/pipeline/slurm/benchmark_neomod3.sbatch
sbatch --array=10-29 neomod/pipeline/slurm/benchmark_neomod3.sbatch

# 2. assemble  (needs the 14.4M-row epoch cache in memory)
srun --account=astro --partition=cpu-g2-mem2x --cpus-per-task=8 --mem=250G --time=00:50:00 \
  --pty bash -c "VDP_LOADER=s3m PYTHONPATH=$PWD/neomod/src \
  conda_prep/bin/python neomod/pipeline/gen_benchmark_tracklets_neomod3.py build"
```

### §11.4.4 Parameters (all at the top of the builder)

| name | value | notes |
|---|---|---|
| `EPOCH` | `2027-08-25T00:00:00` | **must** match the maps AND the Stage-0 cache |
| `MAG_MIN`, `MAG_MAX` | 14.0, 25.0 | the maps' own bin span; v3 had no cut |
| `NEO_SEED` | 20270825 | **must differ from the map cache's 42** — see traps |
| `TARGET_TOTAL` | `None` | `None` = true absolute counts (no scaling) |
| `DT_DAYS` | 30/1440 | tracklet baseline, matches v3 |
| `LON_STEP`, `SUN_EXCLUSION` | 10.0, 40.0 | ⇒ `DLON_LIMIT` 140° |
| `LAT_BASE` | 0,1,2,3,4,5,8,12,18,25,35,50 | symmetric ⇒ 23 latitudes ⇒ **667 cells** |
| `--n-orbits-total` | 10,000,000 per array | ×3 arrays = 30M |

`LON_STEP`, `SUN_EXCLUSION` and `LAT_BASE` **must** equal `sorcha_gen_maps_grid.py`'s or cell
assignment silently points at the wrong maps.

### §11.4.5 Schema (25 columns, identical to v3)

`ObjID population prob_map_file prob_map n_det_per_night | mean_ra mean_dec mean_dra mean_ddec
mean_mag | ra0 dec0 mjd0_utc mag0 ra1 dec1 mjd1_utc mag1 | lam_deg beta_deg
dlon_from_antisun_deg H vlam vbeta e`

Group 2 is what the VDP scorer reads, group 3 is what digest2 reads, group 4 is analysis only.
Column parity with v3 is checked after every build — **do not** add/rename columns without
re-checking downstream scorers.

### §11.4.6 ⚠️ Traps

**1. The NEO seed must differ from the map cache's.** Maps are evaluated on a fixed 1001×1001
velocity grid and a clone typically dominates its own cell, so a test NEO that is also a map clone
gets scored against density it created itself. Silent, and it flatters exactly the population under
test.

**2. Draw at least 11.5M orbits.** `w_new = 11,432,918 / n_drawn`, and running at true absolute
counts needs `w_new ≤ 1`, i.e. **`n_drawn ≥ 11,432,918`**. Below that the NEO pool cannot cover its
own expected count and `build` exits FATAL (by design — see trap 4). 30M gives comfortable margin.

**3. All shards must draw the SAME number of orbits.** `build` infers the total as
`n_orbits_drawn[shard 0] × n_shards`. Mixing shard sizes silently corrupts `w_new`, hence the NEO
count, hence every ratio. Keep `--nshards`/`--n-orbits-total` consistent across array submissions.

**4. Never cap one population.** The scale factor is applied to all four or to none. v1 capped MBA
200k / Trojan 100k while leaving NEO and TNO uncapped (**MBA suppressed ~69×**); v3 kept a milder
form (MBA 650k of 13.88M, no mag cut) which is why its NEO share reads 1.76%. The clamp branch in
`build` is a **fatal error** for this reason — if it fires, draw more orbits or lower
`TARGET_TOTAL`, never let it clamp.

**5. `VDP_LOADER=s3m` is required.** The default loader is *hybrid*, which carries a different
`t_0`.

**6. IERS.** 2027 epochs need `iers.conf.auto_max_age = None` (already set in the builder), or
astropy refuses to transform.

**7. Changing the epoch is not a one-line change** — it needs a matching Stage-0 cache (§10.3) and a
matching NEOMOD3 projection epoch.

### §11.4.7 Built-in validation (all must pass)

| check | criterion | last run |
|---|---|---|
| share integrity | benchmark share − true share ≈ 0 for every population | **0.0000 pp** |
| per-direction realism | (rows ÷ scale) ÷ true objects ≈ 1 in all 5 direction bins | **1.000 ×5** |
| clamp guard | must not fire | did not fire |
| column parity vs v3 | no missing columns | pass |
| cell coverage | 667/667 populated | pass |
| NEO fraction vs direction | monotonically decreasing from antisun | 1.058 → 0.547% |
| independent-draw agreement | 10M vs 30M draws predict the same NEO count | 25,115 vs 24,998 (0.5%) |

The last one is the strongest: two independent realisations agreeing to 0.5% validates the absolute
normalisation end to end.

### §11.4.8 Changing things

| want | change | rebuild |
|---|---|---|
| different epoch | `EPOCH` + Stage-0 cache + NEOMOD3 projection epoch | everything |
| LSST-limit cut instead of map span | `MAG_MAX = 24.5` | `build` only |
| smaller file | `TARGET_TOTAL = <int>` (uniform factor, all populations) | `build` only |
| more NEO statistics | more shards (raises the pool, not the expected count) | shards + `build` |
| S3M NEOs instead | use `gen_benchmark_tracklets_s3m.py` | — |

Note that "more NEO statistics" only helps up to the true count: at `scale = 1` the NEO row count is
fixed at ~25,000 by NEOMOD3's absolute normalisation. Extra draws buy margin against trap 2, not
more NEOs.

## §11.5 WHAT THE BENCHMARK IS AND IS NOT REALISTIC ABOUT

### Realistic: what is THERE

| check | value |
|---|---|
| per-direction counts vs true in-grid objects | ratio **1.000** in all 5 bins |
| mean surface density (mag<25) | **131 / sq deg** |
| within \|ecliptic lat\| < 5° | **483 / sq deg** |
| implied content of one LSST visit (9.6 sq deg) | ~1,258 average, ~4,637 on the ecliptic |

Right order of magnitude, and the population mix, velocities and positions all come from n-body
propagation of real/model populations.

### NOT realistic: what would be OBSERVED

Every survey effect is absent **by construction** — that is the point of a geometric benchmark, but
it bounds what may be concluded from it:

| missing effect | consequence |
|---|---|
| **no footprint** | LSST covers ~1,000–2,000 sq deg/night = **6%** of this sky. A real night samples ~196,000 of these 3,220,714 objects, not all of them. |
| **no detection efficiency** | no fill factor, seeing, airmass, vignetting, chip gaps |
| **hard mag<25 cut** | real single-visit depth is ~24.5 and varies with band/conditions |
| **no trailing losses** | hits fast movers hardest ⇒ the \|v\|>2 NEOs are **over-represented** vs reality — exactly the population of interest |
| **perfect astrometry** | no positional noise |
| **fixed 2 detections, 30 min apart** | real tracklets vary in cadence and multiplicity |
| **no linking losses** | no SSP tracklet-building or linking filter |

**⚠️ Consequence for the digest2 comparison.** digest2 was designed for real, noisy short arcs and its
score depends on arc quality. On noiseless positions it may perform better than it operationally
does, so **a digest2-vs-VDP result on this benchmark is not an operational comparison**. It answers
*"does the classifier separate these populations?"*, not *"what will a night look like?"*. The Sorcha
regeneration answers the second.

### ⚠️ If digest2 cost forces a subsample

3.22M rows ≈ **644 digest2 tasks** (~5,000 rows each). If trimmed:

1. Do it in the **scoring stage**, never in `tracklets_benchmark_neomod3.parquet` — the truth file
   must stay a faithful snapshot.
2. Keeping all NEOs and subsampling only non-NEOs (what the existing stage 4 does with its 500k
   non-NEO cap) is fine for **AUC** — rank-based, invariant to class balance.
3. **But it inflates contamination and precision**, because the prior has changed. Any contamination
   figure from a subsampled set **must be reweighted to the true 0.776% NEO share** before being
   compared with §9.3. Forgetting this is a subtler cousin of the v1 ratio bug.

### ⚠️ Do not compare \|v\|>2 counts across the two test sets

| set | fast NEOs | what it is |
|---|---:|---|
| NEOMOD3 benchmark | **5,050** | one instant, 667 cells, geometric, no losses |
| S3M Sorcha (full) | 7,727 | **~2 years / 606 nights**, detected + linked |
| S3M Sorcha, night 61642 | **27** | one night |

The Sorcha set accumulates over 606 nights; the benchmark is a snapshot. The like-for-like per-night
figure is **5,050 vs 27**, and that gap is survey effects (6% footprint coverage, detection
efficiency, trailing, linking) — **not** a gain from the new cloning. An earlier draft of this
document compared 5,050 against "33" (the two-center subset) as if it were an improvement; that was
wrong and is retracted here.

## §11.6 ⚠️ THE PRIOR PROBLEM — §9.3–§9.9's F1 and contamination were computed at a 22.9% NEO prior

Found while planning the new benchmark's scoring. **This affects the S3M results already recorded.**

### What happened

The linking pipeline's stage-4 `sample` step caps non-NEOs at 500k. The §9.3 scored set is therefore
648,769 rows with 148,773 NEOs — a **22.93% NEO prior**. The physical prior (new benchmark, in-grid,
mag<25) is **0.776%**: about **127 non-NEOs per NEO**, not 3.

Precision — and therefore contamination and F1 — depends on the positive:negative ratio. Drop most
negatives and you drop most false positives, so precision is flattered. **AUC and completeness are
immune** (rank/positive-only statistics).

Reweighting each non-NEO by 38.0× to restore the true prior:

| classifier | AUC (unchanged) | F1 as reported → true | contamination as reported → true |
|---|---:|---|---|
| NEOMOD3 VDP | 0.9288 | 0.8319 → **0.6126** | 12.43% → **41.86%** |
| production VDP | 0.8814 | 0.8020 → **0.6556** | 11.94% → **23.30%** |
| digest2 | 0.9430 | 0.8489 → **0.7218** | 6.60% → **14.14%** |

**Every AUC-based statement in §9.3–§9.9 stands. Every F1/contamination/threshold statement does
not.**

### ❌ But "production beats NEOMOD3" is NOT the right reading — that was a threshold artefact

§9.3 compared the two at **each classifier's own best-F1 threshold**, which sit at *different*
completeness (79.2% vs 73.6%). NEOMOD3 was buying more completeness, so of course it paid more
contamination. Comparing at **matched completeness** (true prior) reverses it:

| completeness | NEOMOD3 (±5) | production VDP (±2) | digest2 |
|---|---:|---:|---:|
| 60% | 41.82% | **32.78%** | **10.38%** |
| 70% | **62.58%** | 76.66% | **37.56%** |
| 80% | **85.49%** | 92.78% | **81.68%** |
| 85% | **92.79%** | 97.23% | 92.94% |
| 90% | **97.10%** | 98.70% | **96.57%** |
| 95% | **98.62%** | 99.26% | **98.04%** |

At 80% completeness NEOMOD3 yields **18,425 false positives vs production's 40,221** — under half, at
a higher threshold (0.0086 vs 0.0039). **NEOMOD3 dominates production from 70% to 95% completeness**;
production wins only at ≤60%, which is precisely where the F1 optimum lands at a 0.78% prior. Both
statements are true at once; the single-number F1 summary just happens to land in production's one
favourable corner.

**Lesson: never compare contamination at each classifier's own operating point.** Match completeness
(or match threshold) or the comparison measures threshold placement, not discrimination.

### The finding that outranks all of it

At the true prior **every** classifier is brutal where it matters: 85% completeness costs **>92%
contamination** for all three. With ~127 non-NEOs per NEO, no single classifier gives a usable
operating point at high completeness. This is the honest framing for the VDP-vs-digest2 question —
it is not "which one wins", it is "none of them is sufficient alone at survey priors."

### Consequence for the new benchmark

**Score all 3,220,714 rows — do not subsample.** ~806 digest2 tasks (4,000 rows each, ~20 min on
`ckpt-all`, from the v3 timings) makes the prior correct **by construction** and removes reweighting
entirely. If a future run must subsample, keep NEOs whole, subsample non-NEOs at a recorded uniform
rate, and pass `sample_weight` to every precision/F1 computation.
