# mag245 n-body digest2 validation package (Stages 1–3)

Reproducible package for the controlled digest2 audit + benchmark/Sorcha validation on the mag245
n-body sets (matched on `ObjID == s3m_objid`, Sorcha night 61642). **No production parquet was
modified.** Built 2026-07-22/23.

## What to open first

**`neomod/mag245_digest2_validation_review.ipynb`** — the self-contained review notebook. It is
already **executed**: every table and figure is baked in, so it reads without re-running or needing
the data files. It also documents versions, digest2 config, and the exact regeneration commands.

## Environment / config (pinned)

- digest2 **v0.19** (compiled 2026-05-28); exec sha256 `3182930840cd5ac6`, model sha256 `426fc29e9fcb3baf`.
- config `noheadings / norms / repeatable / NEO`; obscode **X05**; parse `int(2nd token)/100`, strict
  (missing→NaN never 0; out-of-range rejected/counted). `repeatable` scores are cpu-count-independent.
- VDP: `neomod/src/velocity_density_pipeline_gmm.py` (n-body integrator); maps
  `prob_maps_grid_s3m_nbody/` (667 npz, `support_mask_min=1`, nearest-dist mask OFF, Johnson V).
- astropy 7.2.0 · numpy 2.4.4 · pandas 3.0.2 · scikit-learn 1.8.0 · assist 1.2.3 · rebound 4.6.0 ·
  sorcha 1.1.1 · adam_core 0.5.5.

## Scripts (git-tracked, `neomod/pipeline/`)

| script | stage |
|---|---|
| `audit_digest2.py`, `audit_digest2_exact_input.py` | audit + exact-input repeatability |
| `stage1_deterministic_rescore.py` (+ `slurm/stage1_deterministic_rescore.sbatch`) | Stage 1 |
| `stage2_exact_classifier.py` | Stage 2 |
| `stage3_synchronized_time.py` | Stage 3 |
| `rescore_vdp_Vband.py` | reusable VDP re-score |

## Small tables (git-tracked bundle)

`neomod/mag245_validation_deliverables/` — the summary CSVs: `same_input_delta_summary.csv`,
`parse_issues.csv`, `cross_input_summary.csv`, `STAGE1_summary.csv`,
`{bench,sorcha_rawfix,sorcha_vcorr}_metrics.csv`, `STAGE2_summary.csv`, `STAGE3_decomposition.csv`.

All figures are **embedded in the notebook** (no external files needed to read it). The source PNGs
are on Hyak under each stage's `outputs/…` folder (`.png` is gitignored, so they are not in the repo).

## Full outputs (on Hyak, large; not in git)

| dir | contents |
|---|---|
| `outputs/digest2_exact_input_audit_night61642/` | `per_object_scores.parquet` (794×3×2 runs) + summaries |
| `outputs/mag245_nbody_deterministic_rescore/` | `{bench,sorcha_rawfix,sorcha_vcorr}_deterministic.parquet` (adds `P_NEO_d2_det`), scorediffs, ROC |
| `outputs/mag245_nbody_exact_classifier_test/` | `manifest_{src}.parquet` (canonical tracklets + SHA256, hash-verified), VDP scores, metrics, ROC |
| `outputs/mag245_nbody_synchronized_time_test/` | `stage3_cases_794NEO.parquet`, `stage3_per_object_794NEO.parquet`, decomposition, figures |

Production inputs (unchanged): `outputs/phase2_benchmark_s3m_nbody_mag245/benchmark_comparison_s3m_nbody_mag245.parquet`,
`outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet`.

## Get it onto the Mac

```bash
# scripts + notebook + bundle are in the git repo — from the Mac:
cd <neomod checkout> && git pull
# (or copy just the review notebook, self-contained, from Hyak:)
scp ds2004@klone.hyak.uw.edu:/mmfs1/gscratch/dirac/ds2004/sorcha/neomod/mag245_digest2_validation_review.ipynb .
```

## Headline results

1. **digest2 wins on identical inputs** (hash-verified): ΔAUC +0.065…+0.086, ΔF1 +0.05…+0.07 on all
   three variants; robust to the magnitude band. `bench` reproduces paper Fig. 2 exactly.
2. **Epoch + band explain most of the benchmark↔Sorcha digest2 gap.** Time-sync (B−A) 8.9% of NEOs
   >0.10; band (D−C) 4.8%; residual after both (D−B) only 1.3% >0.10, mean ≈ 0. Sky sep 259″→13″.
3. **VDP less time-sensitive** (B−A 2.6% >0.10) but larger rate-scatter residual (D−B 9.7% >0.10).
4. **Audit integrity**: digest2 `repeatable` has **zero** same-input variance; production's `random`
   mode + silent-zero + un-ceilinged parse were the figure anomalies (the "1.0→0.0" is a real 1.00
   silently zeroed; the 1.17 is a parse artifact). Bug fixes move the aggregate ROC negligibly.

**Not yet run (paused for review):** full 1,762-object Stage 3 expansion (all populations, ROC-capable).
