# Rubin Tracklet Scorer — interface v1 (FROZEN)

A stable command-line contract around the **exact** VDP classifier that produced the accepted
TEST2 results. Future map versions are meant to be swapped underneath this contract without
redesigning the nightly pipeline.

> **v1 is an engineering baseline, not an operational-readiness claim.**
> Four-density coverage on TEST2 is **43.96%** (raw) / **45.06%** (physical-weighted). A large
> minority of in-domain tracklets return `NaN` with an explicit reason. v1 freezes that limitation
> deliberately rather than papering over it.

---

## 1. What v1 does and does not do

v1 **packages**. It does not change science. The density lookup and posterior are delegated
verbatim to `score_test2.score_new_vdp`, the function that produced TEST2's `P_NEO_new` at commit
`87f6bd82d190e946f99ab128ff8ffff380d09a7a`.

Frozen behaviour:

| item | value |
|---|---|
| sky centres | 667 |
| magnitude bins | 44 half-open 0.25-mag apparent-V bins, `14 <= V < 25` |
| k | NEO 150; MBA/TNO/Trojan 10 |
| estimator | Bayesian kNN (closed-form posterior mean) |
| Gaussian smoothing | **off** |
| support masking | **off** |
| velocity grid | [-5, +5] deg/day at 0.01, bilinear **in velocity only** |
| magnitude interpolation | **none** — the single containing 0.25-mag slice is used |
| calibration | **none** (no Platt) |
| posterior | `P(c) = rho_c / sum_all rho` |
| validity | requires **all four** densities present and finite |
| missing densities | never replaced by zero or epsilon |

A 2- or 3-population denominator is **not** `P(NEO)`. Such rows return `NaN` with reason
`missing_population_density:<pops>`.

---

## 2. Commands

```bash
python neomod/pipeline/score_rubin_tracklets.py \
    --input  nightly_tracklets.csv \
    --output nightly_tracklets_scored.csv \
    --model-seal neomod/seals/RUBIN_TRACKLET_SCORER_V1_SEAL.json
```

**Native full-precision CSV only.** MPC-80 ingest was removed from this path; see §5.

`--model-seal` is **required**. The scorer never searches for, or auto-selects, "the latest" maps.
It verifies the sealed map root, the `MAP_BUILD_SEAL_V2.json` hash, and the hash of every scoring
source file, and refuses to run on a mismatch.

Optional: `--chunk-size N` (default 200000), `--summary-json PATH`,
`--geometry-source {derived,precomputed}` (see §7).

---

## 3. Canonical CSV input

| column | units / convention |
|---|---|
| `tracklet_id` | string, **unique**, preserved verbatim in the output |
| `mjd0`, `mjd1` | Modified Julian Date, **UTC**, `mjd1 > mjd0` |
| `ra0`, `ra1` | degrees, ICRS, `[0, 360)` |
| `dec0`, `dec1` | degrees, ICRS, `[-90, +90]` |
| `mag0_V`, `mag1_V` | **apparent V magnitude** (see §4) |
| `observatory_code` | string, carried through; v1 assumes geocentric geometry |

Exactly two detections per row. Extra columns are ignored except the geometry passthrough of §7.

---

## 4. Magnitude is apparent V — this is not negotiable

`mag0_V` / `mag1_V` **must** be apparent V. v1 does **not** convert LSST bands. Feeding an
`r`-band magnitude in as V silently shifts the object into the wrong 0.25-mag slice.

Any Rubin-band → V conversion belongs in a **separate, explicitly versioned preprocessing stage**
upstream of this scorer.

`mean_V = (mag0_V + mag1_V) / 2`, and the bin is the half-open interval containing it:

```
bin(V) = [14.00 + 0.25*floor((V-14)/0.25),  +0.25)
V = 23.4   ->  [23.25, 23.50)
V = 23.25  ->  [23.25, 23.50)      (lower edge is inclusive)
V = 23.50  ->  [23.50, 23.75)
V < 14 or V >= 25  ->  INVALID, reason `v_out_of_range`; never clipped into an edge bin
```

---

## 5. digest2 comparison branch — MPC-80 lives here, not in the VDP path

The VDP scorer reads only the native CSV. MPC-80 emission is a **separate adapter** used solely to
feed the unmodified digest2 binary:

```
native Rubin CSV ──> frozen VDP scorer ──> VDP scores
       │
       └──> export_mpc80_for_digest2.py ──> unmodified digest2 ──> digest2 scores
                                                                 │
                    merge by tracklet_id  <──────────────────────┘
```

```bash
python neomod/pipeline/export_mpc80_for_digest2.py \
    --input   nightly_tracklets.csv \
    --output-mpc nightly_tracklets.mpc \
    --mapping nightly_tracklets_digest2_mapping.csv \
    --audit-roundtrip
```

`--mapping` writes an auditable `tracklet_id` ↔ MPC-80 designation correspondence with both emitted
lines and their sha256, so any digest2 result traces back to the original `tracklet_id` and to the
exact bytes digest2 received. Merge the two branches on `tracklet_id`.

**Rows that cannot be represented in 80 columns** (non-finite time/position/magnitude, out-of-range
coordinates) are **excluded and recorded with a reason** in the mapping — never silently dropped,
and never emitted with a placeholder digest2 would treat as a real measurement.

Designations are `D000000`…, 7 characters. MPC-80 truncates the identifier to 12 columns and 5 are
leading spaces, so only 7 survive; the adapter asserts that the designations actually *emitted* are
unique and identical across a tracklet's two detections. (An 8-character key silently collided once,
and digest2 returned nothing usable.)

### 5.1 MPC-80 round-trip loss is provenance about digest2's input — not a VDP limitation

Round-tripping 5,000 sealed TEST2 tracklets through MPC-80 changes the 0.25-mag bin for **12.56%**
of them and perturbs a 30-minute-baseline velocity by up to 1.3 grid cells, because the format
stores magnitudes to 0.1 and RA to 0.01 s.

That figure describes **the representation digest2 receives**. It is **not** a limitation of the VDP
magnitude-bin assignment, and it **did not affect the completed native-CSV regression**: reading the
native full-precision CSV, magnitude-bin assignment reproduced the oracle **688,688 / 688,688
exactly** (§12.1). Keeping MPC-80 out of the VDP path makes that separation structural rather than
a matter of discipline.

---

## 6. Output schema

One row per input row, in the original input order, with `tracklet_id` preserved.

```
tracklet_id
p_NEO, p_MBA, p_TNO, p_Trojan
rho_NEO, rho_MBA, rho_TNO, rho_Trojan
valid, reason
v_lambda, v_beta, mean_V
map_center_id, magnitude_bin_id, map_id
interface_version, scorer_seal, map_seal, input_hash
```

- Component densities are returned whenever they were read, including on rows that end up invalid.
- `p_*` are populated **only** when all four densities are valid; otherwise `NaN`.
- Rows are never dropped, clipped, or zero-filled.

### Reason codes

| reason | meaning |
|---|---|
| `ok` | four densities valid; probabilities returned |
| `missing_population_density:<pops>` | one or more population densities unavailable in that cell |
| `v_out_of_range` | mean V outside `[14, 25)` |
| `outside_velocity_grid` | \|vlam\| or \|vbeta\| > 5 deg/day |
| `nonfinite_density` / `zero_total_density` | degenerate density at that point |
| `missing_map` | no map file for the assigned centre |
| `nonfinite_field:<col>` / `out_of_range:<col>` | schema validation failure |
| `nonpositive_time_baseline` | `mjd1 <= mjd0` |
| `mpc80_bad_pair_count:N` | MPC-80 key had N != 2 observations |
| `unsupported_magnitude_band:<bands>` | non-V band |

---

## 7. Geometry: cached-model vs observation-derived λ, β

Ecliptic λ, β **are** derivable from observed RA/Dec, and v1 derives them that way. What Rubin
cannot supply is TEST2's *particular cached model-frame* λ, β, which were carried in the
epoch-state cache and produced during n-body propagation rather than from the apparent position.

The 295-row difference below is therefore **cached-model geometry versus observation-derived
geometry** — not "observable versus unobservable ecliptic coordinates". (An earlier draft of this
document, and the `geometry_divergence_from_TEST2` note inside the immutable validation seal, used
that looser phrasing; this section is the corrected wording. The seal is not rewritten — see
`RUBIN_TRACKLET_SCORER_V1_INTERFACE.json`.)

Measured against TEST2's stored values:

| quantity | result |
|---|---|
| λ / β | median 1.8e-03 deg (~6.5 arcsec) |
| **`center_label`** | **688,393 / 688,688 identical — 295 reassigned (0.043%)** |
| `magnitude_bin` | 688,688 / 688,688 identical |
| `vlam` / `vbeta` re-derived from the two detections | max 3.4e-08 deg/day (median 2.3e-11) |

v1 uses the observation-derived values because those are what a survey provides.
`--geometry-source precomputed` consumes stored `lam_deg`, `beta_deg`, `vlam`, `vbeta` columns and
exists **only** to isolate this effect in the regression; it is not for operational use.

---

## 8. Execution model

1. Read and validate every row; preserve original index and `tracklet_id`.
2. Derive motion, ecliptic position, sky centre and magnitude bin.
3. Mark unassignable rows with explicit reasons.
4. **Group by sky centre.** With the default `--chunk-size` (200000, above the largest centre
   group) each ~390 MB map is opened **exactly once** — verified on the full TEST2 run:
   `n_centers = 667`, `n_map_opens = 667`, `largest_center_group = 3028`.
5. Chunk *within* a centre group when `--chunk-size` is smaller than that group. This costs **one
   extra map open per additional chunk** — results are identical, only I/O rises. The run summary
   reports `n_centers`, `n_map_opens` and `one_open_per_center` so the cost is visible rather than
   assumed.
6. Restore original input order before writing.

Output is invariant to input row order and to chunk size. Both are asserted in the regression.

---

## 9. Exit codes

| code | meaning |
|---|---|
| 0 | completed (individual rows may still be invalid, with reasons) |
| 2 | input or schema error |
| 3 | seal missing, malformed, or inconsistent with on-disk maps/sources |
| 4 | internal integrity assertion failed (e.g. probabilities not summing to one) |

---

## 10. Determinism

Given the same input file and the same seal, output is byte-identical. There is no RNG in the
scoring path. Map loads are read-only.

---

## 11. Files

```
neomod/src/rubin_vdp_scorer_v1.py                    library
neomod/pipeline/score_rubin_tracklets.py             CLI
neomod/pipeline/make_rubin_scorer_v1_seal.py         seal builder
neomod/seals/RUBIN_TRACKLET_SCORER_V1_SEAL.json      immutable seal
neomod/examples/rubin_tracklets_v1.csv               example input
neomod/examples/rubin_tracklets_v1_scored.csv        example output
neomod/pipeline/export_mpc80_for_digest2.py          digest2 branch adapter (NOT the VDP path)
neomod/seals/RUBIN_TRACKLET_SCORER_V1_INTERFACE.json interface metadata -> validation seal
neomod/examples/rubin_tracklets_v1_for_digest2.mpc   example MPC-80 export
neomod/examples/rubin_tracklets_v1_digest2_mapping.csv  auditable id <-> designation mapping
```

---

## 12. TEST2 regression (Step 8) — results

Jobs `38808527` (derived) and `38808528` (precomputed), both 01:50:04, under seal
`0e61ec37e263e9a4f9e81ed4ef25ef8fa069f630b8d2003f2ce29ce37de9eaeb`.
Full artifact: `outputs/test2_geometric/TEST2_V1_REGRESSION.json`.

### 12.1 Precomputed geometry — EXACT reproduction

| check | result |
|---|---|
| rows / order / id set | 688,688, order preserved, ids identical |
| `center_label` | **688,688 / 688,688** identical |
| `magnitude_bin` | **688,688 / 688,688** identical |
| `valid` | **688,688 / 688,688** identical |
| `reason` | **688,688 / 688,688** identical |
| `P_NEO` | max abs **3.553e-15**, max rel 2.698e-12, NaN pattern identical |
| `vlam` / `vbeta` | max abs 1.42e-14 / 7.11e-15 |
| coverage | 43.9582% vs 43.9582% |
| reason-count tables | identical |
| **within sealed tolerance** | **TRUE** |

Given the same geometry, the interface reproduces the frozen oracle to floating-point noise, with
every assignment and reason matching exactly.

### 12.2 Derived geometry — the operational path

| check | result |
|---|---|
| `center_label` | 688,393 / 688,688 — **295 reassigned (0.043%)** |
| `magnitude_bin` | 688,688 / 688,688 identical |
| `valid` / `reason` | 4 / 5 mismatches (knock-on from reassignment) |
| `P_NEO` | max abs 3.816e-03 |
| coverage | 43.9579% vs 43.9582% |
| within sealed tolerance | **FALSE — correctly flagged** |

All 295 are one-cell hops between adjacent centres (`dlon+120_lat-08 → dlon+110_lat-08`, etc.),
consistent with the ~6.5 arcsec difference nudging boundary rows. This is **expected and
documented** (§7), not a defect: v1 derives λ, β from the observed RA/Dec, whereas TEST2 used the
cached model-frame values that a survey does not carry.

---

## 13. Performance (Step 9)

Full 688,688-row run, single process, 8 CPUs, `cpu-g2`:

| quantity | value |
|---|---|
| wall | **6,575.7 s (1.83 h)** |
| throughput | **104.7 rows/s** |
| CSV parse | 1.19 s |
| prepare / assignment | 13.56 s |
| scoring (map I/O + density) | 6,544.4 s — **99.52% of wall** |
| peak RSS | 7.77 GB |
| map opens | 667 (one per centre) |
| map volume read | **265.8 GB**, mean 398.6 MB/map |
| effective read rate | 40.6 MB/s |
| cost per (centre, bin) cell | **352 ms** |

### 13.1 Cost scales with centres touched, not rows

| case | rows | centres | rows/s |
|---|---|---|---|
| concentrated | 11,545 | 5 | **165** |
| spread | 688,688 | 667 | **105** |

The concentrated case amortises 2,309 rows per map load against 1,033 for the full grid. Scoring
arithmetic itself is ~0.1 microseconds per tracklet; essentially all wall time is `savez_compressed`
decompression.

**Implication for nightly Rubin use:** a night's tracklets confined to the observed footprint would
touch far fewer than 667 centres, so throughput would be well above 105 rows/s. The lever for a
future v2 is the storage format (uncompressed or memory-mapped arrays), not the estimator.
