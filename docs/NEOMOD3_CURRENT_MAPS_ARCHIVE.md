# NEOMOD3 current maps — ARCHIVE / FREEZE RECORD

**Status: frozen record of the map set that produced the accepted TEST2 results.**
Written 2026-08-31T10:39:51Z. Read-only with respect to the map root — nothing under it was
written, moved or deleted.

> **The archive date is not the map epoch.**
> `archive_created_utc = 2026-08-31T10:39:51Z` is when this record was made.
> The **map reference epoch is `2027-08-25T00:00:00`**, preserved exactly from
> `MAP_BUILD_SEAL_V2.json`. The map set was *built* at `2026-08-17T06:46:31Z`.
> These three timestamps mean different things and must not be conflated.

Machine-readable companion: `outputs/neomod3_current_maps_archive/current_maps_manifest.json`
(0.44 MB; carries the size and SHA-256 of all 2,676 files individually).

---

## 1. Map root and integrity verification

| | |
|---|---|
| map root | `/mmfs1/gscratch/dirac/ds2004/sorcha/outputs/neomod3_mag025_k150_maps_v2` |
| total size | 265.86 GB |
| map files (`.npz`) | **667** (one per sky centre) |
| all files | 2,676 (667 npz + 1,336 parquet + 667 `.ok` + 3 csv + 2 json + 1 md) |
| map seal | `/mmfs1/gscratch/dirac/ds2004/sorcha/outputs/neomod3_mag025_k150_maps_v2/MAP_BUILD_SEAL_V2.json` |
| map seal sha256 | `a010bb1bf2e421a3258ff007f7bc4034369918a3877418620cd14e3f0b9d19ab` |
| scorer seal | `/mmfs1/gscratch/dirac/ds2004/sorcha/neomod/seals/RUBIN_TRACKLET_SCORER_V1_SEAL.json` |
| scorer seal sha256 | `0e61ec37e263e9a4f9e81ed4ef25ef8fa069f630b8d2003f2ce29ce37de9eaeb` |

**Independent re-verification, 2026-08-31T10:39:51Z:** every `.npz` was re-hashed from disk and
compared with the per-centre hashes recorded in `MAP_BUILD_SEAL_V2.json` at build time.

```
sealed centres        667
matched               667
mismatched            0
missing on disk       0
ALL_MAPS_MATCH_SEAL   True
```

The maps on disk are bit-for-bit the maps that were sealed. `MAP_BUILD_SEAL_V2.json` also records
`ALL_PASS = True` for its own build-time gate checks.

The seal chain closes on itself, which is worth stating because it was checked rather than assumed:

- the scorer seal's `map_build_seal_sha256` equals this map seal's own file hash (`a010bb1b…`);
- `split_provenance_mag025.json` hashes to `951b8782…`, matching `split_provenance_mag025_sha256`
  inside the map seal's frozen config;
- `epoch_state_cache` hashes to `b3b22403…`, matching `epoch_cache_sha256_16` in the split provenance;
- `nonneo_split_manifest.parquet` hashes to `7c119b43…`, matching that file's `source_hash16`.

---

## 2. Frozen scientific configuration

| setting | value |
|---|---|
| **map reference epoch** | **`2027-08-25T00:00:00`** |
| sky centres | 667 |
| patch radius | 30.0 deg |
| centre coordinate convention | custom_ecliptic; antisolar-relative longitude dlon and ecliptic latitude lat, degrees |
| centre list | `/mmfs1/gscratch/dirac/ds2004/sorcha/outputs/geometric_density_estimator_ablation/frozen_center_list.txt` |
| velocity grid | [-5.0, 5.0] deg/day at 0.01, shape [1001, 1001] |
| velocity interpolation | bilinear in velocity only |
| magnitude quantity | apparent V (HG, G=0.15) |
| magnitude bins | 44 bins, 14.0 ≤ V < 25.0, step 0.25, half-open [lo, lo+step) |
| magnitude interpolation | **none** — the single containing slice is used |
| k | NEO 150; MBA/TNO/Trojans 10 |
| density engine | closed_form (analytic posterior mean of the sealed estimator) |
| Gaussian smoothing | **False** |
| support masking | **False** |
| calibration | none (no Platt) |
| posterior | `P(c) = rho_c / sum_all rho`, all four populations required |

Example centre from the seal: `dlon+000_lat+00` at
lon 331.203000 deg, lat 0.0 deg.

---

## 3. Population sources and physical normalisation

### NEO

| | |
|---|---|
| source | all-sky GEN realization /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/neomod3_projection_cache_high_allsky/by_pixel (104 pixels, 740,000,000 draws) |
| source sha256 | `40d7d11135348f817e50ef030c41cf7a2757a0fe62844b019c17c8ae33eb69a2` |
| by-pixel tree digest (re-verified) | `f3eaf764c265f17bee0f2ee6de44cf010d493e8189afb2d7fbb020f6188decd7` |
| shards / size | 768 parquet, 229.0 MB |
| effective factor | **64.725384** |
| physical weight per sample | 1.544989e-02 |

### Non-NEO (MBA, TNO, Trojans)

effective_factor = split_fraction f(population, magnitude_bin); physical_weight_per_sample = 1/f

Split fractions are recorded per population **and** per magnitude bin in
`split_fractions_by_population_magbin` in the manifest, and hash to
`951b8782a185fe6b0a140fdf9e5490f2b603a42cd9919a4200abb57d49412640`.

**They are not shared between populations.** Measured across the 42 populated bins, all 42 differ
between the three populations (MBA ≈ 0.60 throughout; TNO 0.52–0.74; Trojans 0.49–0.64). Any pooled
non-NEO estimator must therefore carry a per-sample weight rather than one scalar. This is the fact
that determines the design of the NEO-vs-OTHER pilot.

### Clone factors

`no cloning in the current map set; physical weights come from the NEO effective factor and the non-NEO split fractions only`

---

## 4. Code and environment

| | |
|---|---|
| map-build script | `neomod/pipeline/build_neomod3_mag025_k150_maps.py`  sha256 `c5d8107e412609592114a9525667b804070ef45027e628da14b7702e93ae679e` |
| sealed VDP module | `neomod/src/velocity_density_pipeline_neomod_clone_only.py`  sha256 `a6de18c2197cbcb9f93014712aaad232272578150e2222f28eae723b339f923a` |
| map-build commit | `e4780ef806d4eb4a77e3d384144f8a29ff24ea1e` |
| repo HEAD when archived | `31d16327f136386e730c6def95fef0b0bea3c6f8` |
| all-sky GEN job script | `neomod/pipeline/slurm/build_neomod3_high_allsky.sbatch`  sha256 `3d06e2d4de3aa893d8a1d202797a2a545fdc9572afac339777d70db7fd023637` |
| TEST2 scorer | `neomod/pipeline/score_test2.py`  sha256 `73ae6f872123441a6e7ad5e09d857cbea1613194fd0177fb8a154288ccc3433c` |
| Rubin v1 CLI | `neomod/pipeline/score_rubin_tracklets.py`  sha256 `42ca08fa18f76c1b7fcf13d0ab9ec9f83dc98a6f9a5a1801746b33c842a61649` |
| Rubin v1 library | `neomod/src/rubin_vdp_scorer_v1.py`  sha256 `d599fa7eb1393a5aec471fc255b0f1c9d22058274d30e2b07a6c80ee19f6ba03` |
| environment | `neomod/seals/FROZEN_ENV_SEAL.json`; interpreter `conda_prep/bin/python` |

**Map-build command** (one Slurm array task per centre; nothing runs on a login node):

```bash
sbatch neomod/pipeline/slurm/build_neomod3_mag025_k150_array.sbatch
# each task:
conda_prep/bin/python neomod/pipeline/build_neomod3_mag025_k150_maps.py \
    --task-id $SLURM_ARRAY_TASK_ID --neo-source allsky_v2 --density-engine closed_form
```

Build cost recorded in the seal: 184.0 s for the example centre.

---

## 5. Input artifacts (sizes and SHA-256 re-verified at archive time)

| artifact | path | size | sha256 |
|---|---|---:|---|
| epoch state cache | `outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet` | 2597.8 MB | `b3b2240376a87adfe4f9d067bf8e66b5d946acff7e33b8053bcb7c18ec63d28e` |
| non-NEO split manifest | `outputs/splits/nonneo_split_manifest.parquet` | 85.3 MB | `7c119b435de077bb78666a0d2d54b18efa74f64820214639ab926ab21f413766` |
| split provenance (0.25 mag) | `outputs/splits/split_provenance_mag025.json` | 24.8 kB | `951b8782a185fe6b0a140fdf9e5490f2b603a42cd9919a4200abb57d49412640` |
| frozen centre list | `outputs/geometric_density_estimator_ablation/frozen_center_list.txt` | 10.7 kB | `aec2c2060e566be64e411722a3ea536db61b6326817447ffa81019567489db3f` |

---

## 6. TEST2 files used for the completed scoring

These are the exact files named in
`outputs/neomod3_new_vs_legacy_score_regime/RUN_RECEIPT.json`, re-hashed here and confirmed
unchanged.

| file | size | sha256 |
|---|---:|---|
| `outputs/test2_geometric/TEST2_TRACKLETS.parquet` | 188.7 MB | `e374547893e2945a3ea6dfcd3ae80c1fa4a91ac45a6ce6ebce7be23a552cc09b` |
| `outputs/test2_geometric/TEST2_SCORED.parquet` | 269.8 MB | `8b6e1de609bf4831dc111f658f29a4535dffabaa2cf2dd3fec0c19e0e02ba503` |
| `outputs/test2_geometric/TEST2_V1_REGRESSION.json` | 4.9 kB | `a138563a43951c9d3af25a6e6792e6dcc472bd8a093104f0c6444766076931fc` |

TEST2 contains 688,688 rows. Four-density coverage is **43.96%** raw / **45.06%** physical-weighted
— the known limitation this map set carries, frozen deliberately rather than papered over.

---

## 7. What this archive does not claim

- It does not claim the maps are operationally complete. They are not; see the coverage figure above.
- It does not re-open any frozen CAL/E0/E1 decision.
- It does not modify, rebuild, move or delete the map root, and it creates no new map files.
