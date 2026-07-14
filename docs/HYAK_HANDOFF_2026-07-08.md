# Hyak handoff — 2026-07-08

**Read this first, then run §4.** Full context: `docs/NEOMplan.md` (§9 Log has today's results),
design rationale: `docs/D2_detail.md`.

---

## 1. Why these jobs must run on Hyak (and the others must not)

| job | where | why |
|---|---|---|
| **A. digest2 on the Kurlander referee** | **Hyak only** | `digest2` binary is not installed on Arnor. This is the missing piece of the 4-way comparison table. |
| **B. Full-population Kurlander tracklet build** | **Hyak** | The 479 GB source lives on Hyak scratch after your rsync. Building it from Arnor reads across the *shared Epyc NFS mount* — slower and inconsiderate to other astro users. Also parallelises to ~5 min with a Slurm array vs ~88 min serial. |
| decisive L1-vs-L2 test | already done on Arnor | It was a ~5 min job on a random object subset; no cluster needed. **Result already in NEOMplan §9 (RESULT 3).** |
| ranging-engine scoring | Arnor | ~3 min per 120k tracklets. |

I could not launch these myself: `ssh ds2004@klone` from Arnor returns `Permission denied`
(needs your password/2FA).

---

## 2. What already happened today (so you know what these jobs are *for*)

Three results, all in `NEOMplan.md §9`:

1. **RESULT 1 — the ranging engine beats digest2.** On the 40–70° band (20k tracklets), using the
   *same S3M population* and the *same two detections* digest2 gets:
   L1 (our engine) **F1 0.901 / AUC 0.967** vs digest2 **F1 0.833 / AUC 0.924** → **+0.067 F1**.
   The mid-elongation gap was a digest2 *implementation* gap, not a velocity-vs-orbit-space gap.
2. **RESULT 2 — NEOMOD3 cannot win on an S3M referee** (L2−L1 = −0.0021 ± 0.0009). Expected: the
   v5 truth is S3M-drawn, so an S3M prior is optimal by construction. The referee trap, made empirical.
3. **RESULT 3 — the symmetric flip is confirmed** on your Kurlander referee (NEOMOD3-drawn NEOs,
   same S3M non-NEO ⇒ exactly one variable moves):
   L2−L1 = **+0.0028 ± 0.0004** (P>0 = 1.000) with `neomod3_norm='absolute'`, and a clean **monotone
   per-elongation ramp** −0.0026 (antisun) → **+0.0116 (sunward)**. NEOMOD3 must be run with its own
   N(H) — items 1A and 2B are not separable in practice.

**Honest effect-size hierarchy:** engine vs digest2 = **+0.067 F1**; NEOMOD3 vs S3M prior =
**+0.003 F1**. The implementation matters ~20× more than the population model — *except* sunward
and at faint magnitudes, where NEOMOD3 is significantly better.

**What's missing, and exactly why job A exists:** we have *no* digest2 or VDP-map score on the
Kurlander referee. Without them we cannot say whether the engine still beats digest2 on the *fair*
(NEOMOD3-drawn) referee. That is the last cell of the table.

---

## 3. Files to send (all live in the repo)

| file | purpose |
|---|---|
| `pipeline/kurlander/build_tracklets.py` | detections → nightly tracklets (vectorised; 1.9 s/file) |
| `pipeline/kurlander/slurm_build_tracklets.sbatch` | 40-way array for job B |
| `pipeline/kurlander/merge_shards.py` | merge shards → natural + eval referee sets |
| `pipeline/kurlander/score_digest2.py` | job A (digest2 on the referee) |
| `outputs/kurlander/referee_eval.parquet` | **the 120k-tracklet referee eval set** (11 MB) — the exact rows scored on Arnor, so digest2 lands on identical tracklets |
| `src/ranging_engine.py` | the engine (if you want to re-score on Hyak) |
| `outputs/pop_cache_wide.npz` | population tables the engine needs (~30 MB) |
| `docs/NEOMplan.md`, `docs/D2_detail.md`, this file | context |

---

## 4. Run these

```bash
# ---- 0. one-time ----
mkdir -p /mmfs1/gscratch/dirac/ds2004/logs
cd /mmfs1/gscratch/dirac/ds2004/sorcha        # wherever the repo lives on klone

# ---- A. digest2 on the Kurlander referee  (the important one) ----
# point --d2-dir at your digest2 install (the dir containing the `digest2` binary + MODEL files)
python pipeline/kurlander/score_digest2.py \
    --in  outputs/kurlander/referee_eval.parquet \
    --out outputs/kurlander/referee_eval_d2.parquet \
    --d2-dir /mmfs1/gscratch/dirac/ds2004/digest2
# ~120k tracklets. If it's slow, raise --chunk or wrap in a Slurm job.
# SANITY: digest2 nan frac should print ~0.00. If it's ~1.0, the binary/MODEL path is wrong.

# ---- B. full-population tracklet build (fixes the 1/200 MBA subset caveat) ----
sbatch pipeline/kurlander/slurm_build_tracklets.sbatch     # 40-way array, ~5 min
# when it finishes:
python pipeline/kurlander/merge_shards.py \
    --glob '/mmfs1/gscratch/dirac/ds2004/kurlander2025/tracklets/shard_*.parquet' \
    --out  /mmfs1/gscratch/dirac/ds2004/kurlander2025/tracklets/referee.parquet
```

Then `scp` back to Arnor: `outputs/kurlander/referee_eval_d2.parquet` and
`.../tracklets/referee*.parquet`, and I'll finish the 4-way table.

---

## 5. Things that will bite you (read before debugging)

1. **`slurm_build_tracklets.sbatch` assumes** `--account=dirac --partition=ckpt` and that
   `conda activate neofast_py310` works from `~/.bashrc`. Fix those two lines if not.
2. **`--root`** in the sbatch is `/mmfs1/gscratch/dirac/ds2004/kurlander2025` and expects the
   *Epyc* subdirectory naming you rsync'd: `neomod/outfiles/`, `s3m/outfiles/`, `trojanmod/outfiles/`,
   `cfeps/outfiles/`. If your rsync flattened it, pass a different `--root` or fix `POPS` in
   `build_tracklets.py`.
3. **RA-rate convention (cost us real time today).** `mean_dra` written by `build_tracklets.py` is
   the **RAW α̇**, *not* α̇·cos δ — matching the v5 parquet. The engine must therefore be called with
   `dra_cosdec=False`. `score_digest2.py` likewise propagates `ra1 = ra0 + dra*dt` with no cos δ
   divide, exactly as the v5 digest2 baseline did. Do not "fix" either one.
4. **F1 is prevalence-dependent.** `referee_eval.parquet` is deliberately subsampled to NEO
   fraction **0.293** to match v5 so F1 is comparable. Never compare F1 across referees with
   different NEO fractions; use AUC for that.
5. **Do not use** `small_neo_output.h5` (its linked count disagrees with the paper) or the `_022` /
   `.csv` dev leftovers. `build_tracklets.py` already points at `neo_output_1.h5` only.
6. **Known confound to control later:** the Kurlander referee's non-NEO population is the *same*
   S3M our denominator uses, with the same magnitude cutoff — so the faint-end NEOMOD3 gain is
   partly circular. The *elongation ramp* is unaffected and is the robust result. See NEOMplan §9.

---

## 6. What I'll do with the outputs

- `referee_eval_d2.parquet` → complete the 4-way table (VDP · digest2 · L1 · L2) on the fair referee,
  answering: **does the engine still beat digest2 when the truth is NEOMOD3-drawn?**
- full `referee.parquet` → recompute RESULT 3 at the true natural population ratio, removing the
  1/200-subset caveat.
