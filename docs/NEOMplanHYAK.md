# NEOM Plan — Hyak-side log (Kurlander referee jobs)

**Date:** 2026-07-09 → 2026-07-11
**Context:** Arnor scp'd over the Kurlander-referee handoff (`HYAK_HANDOFF_2026-07-08.md`,
`D2_detail.md`, `NEOMplan.md`) plus the four job files (`build_tracklets.py`,
`slurm_build_tracklets.sbatch`, `merge_shards.py`, `score_digest2.py`) and the 120k-tracklet
`referee_eval.parquet`. This doc logs what actually ran on Hyak, every bug fixed, and the
final deliverables — the Hyak-side counterpart to `NEOMplan.md` (owned by Arnor).

Two jobs were requested:
- **Job A** — digest2 on the Kurlander referee (the missing cell in Arnor's 4-way table;
  digest2 binary isn't installed on Arnor).
- **Job B** — full-population Kurlander tracklet build (479 GB source only reachable
  efficiently from Hyak scratch after the earlier rsync; removes the 1/200 MBA-subset
  caveat in Arnor's existing result).

---

## 1. Bugs found and fixed before anything would run

The handoff doc's own §5 ("things that will bite you") flagged two of these in advance —
both were confirmed real and fixed. A third (digest2 path) and a scaling issue (found
independently) were not anticipated.

### 1.1 `slurm_build_tracklets.sbatch` — invalid Slurm account + missing conda env
- `--account=dirac` does not exist in this user's associations (valid: `astro`,
  `astro-ckpt`, `ckpt-astro`, `ckpt-dirac`). Every job this session uses `--account=astro
  --partition=ckpt-all`; changed to match.
- `conda activate neofast_py310` — no such env on this system. Replaced with the
  `conda_prep/bin/python` interpreter used by every other script in this repo.
- Added an explicit `--chdir=/mmfs1/gscratch/dirac/ds2004/sorcha` (every other sbatch
  script in the repo sets this so relative paths resolve regardless of submit dir) and
  an explicit `--error=` log path.
- `build_tracklets.py` itself needed **no changes** — its `POPS` dict
  (`neomod/outfiles/`, `s3m/outfiles/`, `trojanmod/outfiles/`, `cfeps/outfiles/`) matched
  the rsync'd directory layout exactly.

### 1.2 `score_digest2.py` — wrong `--d2-dir` in the handoff doc's own example command
The handoff doc's §4 suggested `--d2-dir /mmfs1/gscratch/dirac/ds2004/digest2` — that path
doesn't exist. The real digest2 install (binary + `.model` file) lives inside the repo at
`/mmfs1/gscratch/dirac/ds2004/sorcha/digest2/` (used throughout this session's earlier
benchmark/case1-2-3 digest2 work). Corrected at call time; no code change needed.

### 1.3 Job A run on the shared login node, then found to scale badly with big chunks
First attempt ran `score_digest2.py` directly (backgrounded with `nohup`) rather than via
Slurm, since the handoff doc suggested this was fine for ~120k rows ("if it's slow, raise
`--chunk` or wrap in a Slurm job"). Diagnosis:
- Chunk 1 (20,000 tracklets = 40,000 MPC80 lines, one digest2 subprocess call) took
  **4,793 s = 80 minutes**, sustained at 169% CPU on the **shared login node** — both
  inappropriate there and far slower than useful (~8 h projected for all 120k rows serially).
- **Fix:** added `--row-start`/`--row-stop` args to `score_digest2.py` so the input can be
  sliced before scoring, then built `slurm_score_digest2.sbatch`: a 24-task Slurm array,
  4 GNU-parallel workers per task, 1,250 rows/worker (`OUTER_CHUNK=5000, INNER_CHUNK=1250`)
  — the exact pattern already proven for case1/2/3 and benchmark v1/v2/v3 digest2 stages
  earlier this session. 24 × 5,000 = 120,000 rows exactly.
- Killed the login-node process and its digest2 child cleanly before resubmitting.

### 1.4 `merge_shards.py` — OOM on the login node
Running the merge directly (per the handoff doc's own §4 command) OOM-killed (exit 137).
The 40 shards total ~19 GB on disk; loaded + concatenated in pandas this exceeds the
login node's usable memory (Job B's task 0 alone produced 1.85M tracklets from the NEO
file, before the ~2,700 MBA/TNO chunks are added). **Fix:** wrapped the same merge call in
`slurm_merge_shards.sbatch` (128 G mem, single task) — same lesson as 1.3, generalized:
**no data-touching step in this pipeline is small enough for the login node.**

### 1.5 An unexplained tool-connection outage (session-level, not a pipeline bug)
Mid-session, every `sbatch`/`squeue`/`sacct`/`cat`-on-recent-file call started failing
with a permission-layer "Stream closed" error, while trivial commands (`echo`, `date`,
`whoami`) kept working. Root cause: a genuine Hyak/connection drop on the user's side —
confirmed when the user re-ran `squeue` locally and saw no `kurl_d2` job at all, meaning
the submission never went through despite the sbatch script being correct. Resubmitted
cleanly once the connection was restored (job 36982174). No pipeline files were at fault.

---

## 2. What actually ran (chronological)

| # | Job | ID | Result |
|---|---|---|---|
| 1 | Job B — tracklet build, 40-way array | 36968470 | 40/40 COMPLETED |
| 2 | Job A — digest2, serial, login node (**bad, killed**) | PID 78315 | killed after 80 min/chunk, ~8h projected |
| 3 | Merge shards, login node (**bad, OOM**) | PID 3026897 | exit 137 (OOM-killed) |
| 4 | Job A — digest2, parallelized 24-task array | 36982174 (first attempt never submitted — connection drop) → **36982174** (resubmitted) | 24/24 COMPLETED |
| 5 | Merge shards, proper Slurm job (128G) | 36982208 | COMPLETED |
| 6 | Merge 96 digest2 shards → final parquet (login node, small/safe) | inline | 120,000 rows, nan frac 0.00039 |

---

## 3. Final deliverables

| File | Size | Contents |
|---|---|---|
| `outputs/kurlander/referee_eval_d2.parquet` | 11 MB | digest2 `P_NEO_d2` scored on the exact 120k tracklets Arnor already scored with VDP/L1/L2 — the missing cell in the 4-way table. **Sanity check passed:** nan frac = 0.00039 ≈ 0.00 (handoff doc's own accept criterion). Score distribution sane: mean 0.314, median 0.09, full 0–1 range. |
| `kurlander2025/tracklets/referee.parquet` | 493 MB | Full-population tracklet build, subsampled to NEO frac 0.293 (matches v5 eval set for F1 comparability) — **6,346,146 tracklets**. Removes the 1/200 MBA-subset caveat. |
| `kurlander2025/tracklets/referee_natural.parquet` | 19 GB | Full natural-ratio population — **337,029,560 tracklets** (MBA 334,667,889 / NEO 1,852,257 / TNO 509,414), true NEO frac 0.005517. For the true-ratio Result-3 recompute Arnor's handoff doc mentioned. **Not yet sent** — large transfer, held pending Arnor's ask. |

### Flag for Arnor (not yet resolved)
The merged natural population has **only `MBA`/`NEO`/`TNO`** — no `Trojan`, despite
`trojanmod` being in `build_tracklets.py --pops`. Not investigated further on this pass;
worth checking whether `trojanmod/outfiles/*.h5` has a column-name mismatch against
`build_tracklets.py`'s `NEED` list (`ObjID, fieldMJD_TAI, RA_deg, Dec_deg,
trailedSourceMag, optFilter, Linked`) before trusting any Trojan-dependent conclusion
from this referee set.

---

## 4. Handoff back to Arnor

```bash
scp /mmfs1/gscratch/dirac/ds2004/sorcha/outputs/kurlander/referee_eval_d2.parquet \
    /mmfs1/gscratch/dirac/ds2004/kurlander2025/tracklets/referee.parquet \
    ds2004@arnor.astro.washington.edu:/astro/users/ds2004/vdp/outputs/kurlander/
```

`referee_natural.parquet` (19 GB) held back pending explicit request — flagged above.

With `referee_eval_d2.parquet` in hand, Arnor can complete the 4-way table (VDP · digest2
· L1 · L2) on the fair (NEOMOD3-drawn) referee and answer the open question from
`HYAK_HANDOFF_2026-07-08.md`: **does the ranging engine still beat digest2 when the truth
is NEOMOD3-drawn, not S3M-drawn?**

---

## 5. Artifacts added/changed in this repo

- `neomod/pipeline/kurlander/slurm_build_tracklets.sbatch` — fixed (account, conda, chdir)
- `neomod/pipeline/kurlander/score_digest2.py` — added `--row-start`/`--row-stop`
- `neomod/pipeline/kurlander/slurm_score_digest2.sbatch` — **new**, parallelized digest2
- `neomod/pipeline/kurlander/slurm_merge_shards.sbatch` — **new**, memory-safe merge
- `neomod/pipeline/kurlander/merge_shards.py` — unchanged, ran as-is under Slurm
- `neomod/pipeline/kurlander/build_tracklets.py` — unchanged
