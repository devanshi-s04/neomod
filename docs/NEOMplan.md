# NEOM Master Plan v2 — one document, all tasks, all decisions

**Created:** 2026-07-08 (v1). **Rewritten as v2:** 2026-07-08 (Arnor, after 1B closure + stack ceiling
measurement + benchmark-v2 audit). **Supersedes:** NEOMplan v1, `WORK_BREAKDOWN_1A_1B_HYAK.md`.
**Owners:** Devanshi (all Hyak jobs, advisor decisions), Arnor-Claude (analysis, scripts, figures, docs).

**How to read this doc:**
- §1 = facts that are settled — do not re-litigate.
- §2 = the goal, reframed once, clearly.
- §3 = **the decision register (D1–D6)** — every open choice, its options, recommendation, who decides.
- §4 = **the task cards (T0–T5)** — every task: machine, inputs, steps, time, output, accept criterion.
- §5 = execution timeline + dependency graph.
- §6 = standardized evaluation protocol (so numbers stop moving between notebooks).
- §7 = closed items (things we will NOT do, and why).

---

## §1. Frozen facts (settled — the foundation)

**F1. Current headline numbers** (Sorcha hybrid v5, 667-map GMM grid, support mask, 707,670 eval tracklets):

| Elongation band | N | VDP F1 | digest2 F1 | VDP+d2 STACK F1 |
|---|---|---|---|---|
| full sky | 707,670 | 0.808 | 0.836 | **0.865** |
| 0–20° | 106,168 | **0.876** | 0.848 | 0.890 |
| 20–40° | 92,110 | 0.817 | 0.859 | 0.879 |
| 40–70° | 102,411 | 0.762 | 0.839 | 0.862 |
| 70–110° | 47,316 | 0.771 | 0.777 | 0.807 |
| 110–180° | 5,830 | 0.866 | 0.863 | 0.888 |

Full-sky AUC: VDP 0.880, digest2 0.930, **stack 0.944**. Pure-S3M control: 0.815 / 0.851 / **0.879** F1
(AUC 0.880 / 0.944 / **0.956**). Stack = 2-feature logistic on (logit P_vdp, logit P_d2), 50/50
train/test split, seed 0.

**F2. 1B is CLOSED.** Conditioning resolution (10° lon grid, nearest-centre snapping, 1-mag bins) is
NOT the bottleneck. On a corr=1.000 correct baseline (mask flag fixed:
`mask_radius_deg_per_day=np.inf` must be passed to `from_npz` to match production scoring), both a
one-lon-step mis-assignment and 3-nearest interpolation give **+0.006 F1 = the noise floor**. The
earlier +0.039 "interp gain" was an artifact of the nearest-clone-distance mask being wrongly enabled.
Do not build sky-interp, continuous-mag interp, or the 5° denser grid.

**F3. The 40–70° gap (−0.077 on Sorcha) is information digest2 extracts that velocity-map binning
cannot:** exact per-tracklet (d, ḋ) geometry + the V↔d↔H magnitude coupling. Same observables, AUC
0.93 — the information is in the data; VDP's parameterisation doesn't extract it.

**F4. The stack gain (+0.014 to +0.030 F1) is uniform across ALL elongation bands** — including 0–20°
where VDP wins and 70–110° where both are weak. VDP and digest2 see partially independent signal
everywhere. This is measured, banked, and costs zero new modeling.

**F5. Benchmark v2 is NOT currently a usable referee.** Audit (2026-07-08, Arnor): in
`docs/benchmark_comparison_s3m_v2.parquet`, median `P_NEO_vdp` = 0.0 **for the NEOs themselves**
(AUC 0.393, i.e. anti-correlated); the same degeneracy appears in
`outputs/phase2_benchmark_s3m/benchmark_comparison_s3m.parquet` (the v1-caps file). Neither file
reproduces the notebook numbers (v1 F1=0.787, v2 F1=0.627) when scored raw — the notebooks must apply
filters not encoded in the parquets. Until T0.2 explains this, **no strategic conclusion may rest on
benchmark v1/v2**, including v1's "single most important diagnostic" role in plan v1. The Sorcha
cadenced run is the primary referee (it reproduces exactly).

**F6. The prior NEOMOD3 null result tests nothing.** It only added NEOMOD3 orbits to the GMM
*training set* (changed cloner support, not density weights; maps still normalised to S3M counts).
The reweighting idea (2A/2B) and the ranging prior (1A) remain untested.

**F7. Referee trap.** Truth sets drawn from S3M structurally penalise NEOMOD3-trained classifiers.
Any "NEOMOD3 helps" claim requires a fair referee (see D1) — this gates all of Tier 2.

---

## §2. The goal, reframed once

Old framing: *"close the F1 gap to digest2."* — **Obsolete.** The free VDP+digest2 stack already
beats digest2 by +0.029 F1 / +0.014 AUC everywhere.

**New framing — the program has three claims, in order of value:**

1. **[C1 — banked]** VDP and digest2 are complementary; their combination sets a new short-arc
   classification bar (F1 0.865 / AUC 0.944), and VDP alone wins in the antisun discovery sweet spot
   (0–20°: 0.876 vs 0.848).
2. **[C2 — the crux]** A NEOMOD3-debiased population prior improves short-arc NEO classification over
   the S3M-2011 prior both classifiers currently inherit. Demonstrated two ways:
   ranging framework (1A: NEOMOD3-ranging vs digest2 on identical tracklets) and velocity framework
   (2A/2B reweighting), with a clean model-vs-algorithm control (2D).
3. **[C3 — endgame]** The best config recovers real NEOCP NEOs that digest2 misses at detection
   (Keys' 14%) and improves follow-up queue purity (Wagg's metric).

Every task below serves one of C1/C2/C3.

---

## §3. DECISION REGISTER — every open decision, in one place

| ID | Decision | Blocks | Who | When needed |
|---|---|---|---|---|
| D1 | Fair referee for NEOMOD3 claims | T1.1 config, all of Phase 4 | **Advisor + you** | before T1.1 launch |
| D2 | 1A ranging-term design (3 sub-questions) | T2.1 | **Advisor + you** | before T2.1 build |
| D3 | Benchmark v1/v2 fate | paper §benchmark, F5 | you (after T0.2) | end of Phase 0 |
| D4 | Paper headline: stack-first or NEOMOD3-first | paper structure | **Advisor** | after Phase 2 results |
| D5 | 1A go/no-go for production scale-up | T2.2 | you (criterion below) | after T2.1 |
| D6 | Tier 4 real-data go | T5.2 | you (after T1.2 archive check) | any time |

### D1 — the referee for NEOMOD3 claims  ⟶ RECOMMENDATION: option (a)
The evaluation-referee trap (F7): retrain to NEOMOD3, evaluate on S3M-drawn truth → F1 *drops* for a
biased-referee reason, not a real one.
- **(a) NEOMOD3-drawn Sorcha run** — NEO input sampled from `neomod3_sampler` (orbits + H), non-NEO
  stays S3M, same v5 config/linking as the production run. A ~case-run-sized Hyak job. **Only clean
  fair referee currently available.** → T1.1.
- **(b) Real NEOCP data** — the ultimate referee but the archive is thin (3 snapshots, all
  2026-05-07). Blocked until T1.2 verifies the cron. Endgame (C3), not the Phase-4 referee.
- **(c) Report both S3M-referee and NEOMOD3-referee numbers** — transparent, costs one extra
  evaluation pass, defuses the reviewer question "did you just move the target?" **Do this
  regardless** once (a) exists.

### D2 — 1A ranging-term design (the three sub-questions)  ⟶ RECOMMENDATIONS inline
1. **(d, ḋ) grid:** log-spaced d from 0.05 to 100 AU (digest2's range — do NOT truncate at 5 AU or
   the denominator loses the distant-population hypotheses that suppress TNO/Trojan false positives),
   ~60 pts; ḋ ±0.05 AU/day (~±87 km/s, generous), ~40 pts; include the phase-space Jacobian
   (digest2's `dRdM` weighting). Start 60×40, halve/double to confirm convergence.
2. **Non-NEO denominator:** S3M MBA+TNO+Trojan binned in (q, e, i, H) — same binning digest2 uses —
   built empirically from the raw `.s3m` census (no GMM smoothing; the census is millions of objects,
   bins are well-populated). Apply Wagg f=0.80 NEO-side rescale only in 2B, not here (keep 1A a pure
   prior-shape test first).
3. **Magnitude handling:** marginalise, don't point-estimate. For each (d, ḋ) hypothesis, H is
   *determined* by (V_apparent, d, phase angle) up to phase-function scatter — so H rides along the
   d-grid for free (this IS digest2's V↔d↔H coupling). No extra grid axis needed; just evaluate
   ρ(q,e,i,H(d)) at each grid node.
   **Score definition (the fix for the old failed attempt):**
   `P_NEO_rng = Σ_grid w·ρ_NEOMOD3 / (Σ_grid w·ρ_NEOMOD3 + Σ_grid w·ρ_S3M,non-NEO)` — a class ratio,
   never a raw weight.

### D3 — benchmark v1/v2 fate (after T0.2 diagnosis)
- If the parquets' `P_NEO_vdp` is genuinely broken (wrong mask/epoch/assignment at scoring time):
  fix + re-score on Hyak (cheap — single epoch, ~475k rows), keep v2 as the controlled-setting
  benchmark.
- If the notebooks used filters the parquets don't encode: document the filter, decide whether the
  filtered claim is honest, and either way **remove "benchmark v2 is the single most important
  diagnostic" from plan + paper** — the Sorcha run is the referee of record (F5).

### D4 — paper headline (advisor call, after 1A prototype result)
- If 1A ≥ digest2 on identical tracklets → headline is C2 (debiased prior wins), C1 supporting.
- If 1A < digest2 → headline is C1 (complementarity + antisun win + stack), C2 lives in 2A/2B/2D
  as the model-improvement story. Either way the paper's "physics, not a bug" framing gets softened
  per F3.

### D5 — 1A go/no-go (pre-registered accept criterion — write it down NOW so we don't rationalise later)
On the Sorcha 40–70° band (same tracklets, same eval protocol §6):
- **GO** (production scale-up + paper claim) if `P_NEO_rng` best-F1 ≥ digest2's 0.839, OR if
  stack(VDP, rng) ≥ stack(VDP, d2) = 0.862 — either means the NEOMOD3 prior matches/beats digest2's
  S3M prior in its own framework.
- **ITERATE** (one round: grid resolution, denominator binning, phase function) if F1 ∈ [0.80, 0.839].
- **NO-GO** if < 0.80 after one iteration round → C2 falls back to 2A/2B/2D only; ship C1.

### D6 — Tier 4 go
Gate: T1.2 shows the NEOCP cron has been accumulating snapshots. If dead: restart it NOW (data
accrues while everything else runs), and Tier 4 slips to whenever ~4–6 weeks of stream exist.

---

## §4. TASK CARDS

Machine rule (unchanged): *needs population table / catalog / map .npz / digest2 / Slurm → Hyak.
Operates only on existing parquet columns → Arnor.* Handoff artifact is always a scored parquet
(Hyak adds a `P_NEO_*` column → scp → Arnor evaluates per §6).

### ── PHASE 0 — Arnor, this week, zero Hyak wait ──

**T0.1 Bank the stack (C1). — Arnor. ~1 h. NO dependencies.**
- Turn the ad-hoc stack measurement into a reproducible artifact: `src/stack_vdp_d2.py`.
- Inputs: `outputs/phase2/sorcha_comparison_v5_masked.parquet`,
  `outputs/phase2_s3m/sorcha_comparison_s3m.parquet`.
- Steps: logistic on (logit P_vdp, logit P_d2), 50/50 split seed 0; report full-sky + per-band
  F1/AUC/compl/contam per §6; save `outputs/stack_scores_{hybrid,s3m}.parquet` + one ROC/F1 figure.
- Also fit the 3-feature variant (+ |elongation|) — costs nothing, tests whether band-aware blending
  adds anything beyond the two scores.
- **Output:** script + 2 parquets + figure + numbers table for the paper.
- **Accept:** reproduces F1 0.865/0.879 (±0.002).

**T0.2 Audit benchmark v1/v2 parquets (feeds D3). — Arnor (diagnosis), Hyak (re-score if needed). ~2–3 h.**
- Why is median `P_NEO_vdp` = 0 for true NEOs in BOTH benchmark parquets, and why do raw-parquet
  F1s (0.106 v2 / — v1) not reproduce the notebook numbers (0.627 / 0.787)?
- Steps: open `benchmark_v5_normalisation_s3m.ipynb` + v2 notebook; extract the exact row filters and
  scoring calls; check the scoring-time mask flags (the SAME `mask_radius` bug class as F2 — prime
  suspect), map-assignment epoch, and whether P=0 rows are "outside all mag bins".
- **Output:** short writeup appended to this doc under D3; if broken → 1-line fix + Hyak re-score
  request (~30 min job).

**T0.3 Per-elongation operating thresholds (Tier-3, operational). — Arnor. ~1 h.**
- On stacked + raw scores: per-band optimal thresholds and the resulting global
  completeness/contamination table (an operations deliverable — MPC-style users don't retune AUC).
- **Output:** threshold table in `outputs/`, feeds paper ops section.

**T0.4 Reconcile "Sec-13 vs benchmark-v2" headline conflict. — Arnor. ~1 h. After T0.2.**
- Expected resolution: v2's VDP column was degenerate → conflict dissolves. Verify Sec-13's
  direction-controlled result stands on its own; write the 1-paragraph resolution for the paper.

### ── PHASE 1 — Hyak, START IMMEDIATELY (long pole; runs while everything else proceeds) ──

**T1.1 NEOMOD3-drawn Sorcha referee run (2C-a; needs D1 confirmed). — Hyak. ~days wall (case-run-sized).**
- The single most schedule-critical job: every Phase-4 result waits on it.
- Steps:
  1. Sample a full NEO population from `neomod3_sampler` (orbits + H, matched in number to the S3M
     NEO input of the v5 run). Non-NEO input: unchanged S3M.
  2. Run Sorcha with the identical v5 config + linking case as production.
  3. Phase-1 tracklets → Phase-2 VDP scoring (current maps) → Phase-3 digest2.
- **Output:** `sorcha_comparison_neomod3ref.parquet` (same schema as v5 masked parquet).
- **Accept:** row counts within ~2× of v5 run; baseline VDP + d2 numbers on it are sane
  (this ALSO immediately tells us how much the current S3M-trained classifiers degrade on
  debiased-real-world-like NEOs — a paper number by itself).

**T1.2 Verify NEOCP cron/archive (gates D6). — Hyak. 10 min.**
- `ls neomod/neocp_data/raw/` — count snapshots + date range. If stalled: restart cron now.
- **Output:** one line in this doc at D6.

### ── PHASE 2 — Hyak, the crux (needs D2; do not start coding before D2 is settled with advisor) ──

**T2.1 1A ranging-term prototype ("VDP-R"). — Hyak build + score, Arnor evaluate. ~2–3 days total.**
- Per D2: for each tracklet (α, δ, rates, V, epoch) evaluate the (d, ḋ) grid; each node →
  heliocentric (r, ṙ) → orbit hypothesis (a, e, i) via the attributable/Gauss method (reference:
  digest2's own solver structure; reuse `NEO_H.py` machinery where sound) → H(d) from V + phase →
  look up ρ_NEOMOD3(a,e,i,H) [numerator] and ρ_S3M,non-NEO(q,e,i,H) [denominator] → Jacobian-weighted
  class-ratio score (D2.3).
- Test set: **Sorcha hybrid v5, 40–70° band** (102k test-half tracklets; NOT benchmark v2 — F5).
- Sanity gates before trusting F1: (i) score distribution not degenerate at 0/1; (ii) NEO median ≫
  non-NEO median; (iii) a handful of hand-checked tracklets have physically sensible admissible
  (d, ḋ) regions.
- **Output:** `1A_neomod3_range_prototype.py` + `P_NEO_rng` column on the band → scp to Arnor.
- **Arnor follow-up (same day):** §6 evaluation + stack(VDP, rng) + stack(VDP, d2, rng) → apply D5.
- **Accept/iterate/kill:** per D5 (pre-registered).

**T2.2 1A production scale-up (only on D5 = GO). — Hyak. ~½ day coding + Slurm array.**
- Refactor to `compute_ranging_score(df, neomod3_table, non_neo_hist)`; Slurm array over
  elongation/time partitions; score the full 707k eval parquet first, the 40.7M stream only when the
  paper needs it.
- **Output:** `P_NEO_rng` on the full eval parquet; full-sky + per-band results.

### ── PHASE 3 — Hyak, mechanical, run in parallel with Phase 2 (no decisions needed) ──

**T3.1 2D control: digest2 with a NEOMOD3-rebuilt population model. — Hyak. ~1 day.**
- Rebuild digest2's binned NEO model from NEOMOD3 samples via MUK (`digest2/archive/*/muk.c`);
  non-NEO bins unchanged. Rerun digest2 on the standard eval tracklets.
- This is the experiment that cleanly separates **model advantage from algorithm advantage**, and it
  reads on 1A: if d2+NEOMOD3 > d2+S3M, the prior matters independent of our pipeline. The control
  Mario/Željko will ask for.
- **Output:** `P_NEO_d2_neomod3` column → Arnor evaluates (incl. on the T1.1 referee when it lands).
- Interpretation matrix (pre-registered):
  d2+N3 > d2 on NEOMOD3-referee AND ≈ on S3M-referee → prior helps, referee trap confirmed → C2 core
  evidence. d2+N3 ≈ d2 everywhere → NEOMOD3's shape change doesn't matter at Rubin depths → C2
  weakens, C1 becomes the paper.

### ── PHASE 4 — Hyak, the NEOMOD3-inside-VDP package (needs D1 + T1.1 referee parquet) ──

**T4.1 2A: reweight NEO density in map generation. — Hyak (667-map regen, Slurm). ~1–2 days.**
- Each source NEO gets w = ρ_NEOMOD3(H,a,e,i)/ρ_source(H,a,e,i); w carries through GMM cloning into a
  weighted kNN density. Corrects S3M's H-independence per mag bin — the faint bins that dominate the
  Rubin stream are exactly where NEOMOD2 says S3M's (a,e,i) shape is most wrong.
- **Output:** `prob_maps_grid_neomod3w/` + re-scored eval parquets (both referees, per D1-c).

**T4.2 2B: renormalise per-mag-bin prior odds. — Hyak (bake into map lookup). ~½ day, do together with T4.1 regen.**
- NEO map integral per mag bin set from NEOMOD3 N(H) (through the pV(H) albedo model); MBA side
  rescaled by Wagg f=0.80 (or self-calibrated to MPCORB at m < 20.5 — try both, report both).
- **Output:** variant maps; evaluate 2A-only, 2B-only, 2A+2B (ablation — reviewers will want it).
- **Accept (C2, velocity framework):** 2A+2B beats baseline VDP on the NEOMOD3 referee without losing
  on the S3M referee's antisun band.

### ── PHASE 5 — opportunistic + endgame ──

**T5.1 Tier-3 map-construction wins. — Hyak, piggyback on the T4.1 regen (same Slurm budget).**
- GMM components 80→200 (one line); fitted-GMM analytic density instead of kNN (kills bleed + the
  entire mask machinery); "unknown-population" maps (subtract MPCORB-knowns; digest2-noid analog).
  Fold the first into T4.1's regen; the analytic-density change is a bigger refactor — only if the
  regen cadence makes it cheap.

**T5.2 Tier-4 real NEOCP (C3; gated by D6). — Hyak score, Arnor metrics.**
- Best config vs MPC's published D2 per object: (a) fraction of the ~14% digest2-missed NEOs
  recovered; (b) top-N follow-up queue purity (Wagg metric).

**T5.3 Paper integration. — Arnor, continuous.**
- Soften "physics, not a bug" (F3); retire/repair benchmark v2 per D3; add C1 stack section
  (T0.1 numbers); add C2 results as they land; ops thresholds (T0.3).

---

## §5. Timeline + dependency graph

```
 week 1            week 2                  week 3+
 ───────           ───────                 ────────
 T0.1 stack ───────► paper C1 section
 T0.2 v2 audit ──► D3
 T0.3 thresholds
 T0.4 reconcile
 T1.2 cron check ─► D6 ..................► T5.2 (real data, when archive suffices)
 [D1 w/ advisor] ─► T1.1 NEOMOD3 Sorcha ══════════► referee parquet ─► T4.1+T4.2 ─► C2(velocity)
 [D2 w/ advisor] ─► T2.1 1A prototype ─► D5 ─► T2.2 production ─────────────────► C2(ranging)
                    T3.1 2D control (parallel, anytime) ────────────────────────► C2(control)
                                                       [D4 headline w/ advisor]
```
Critical path: **D1 → T1.1** (days of wall time; nothing in Phase 4 can start without it) — which is
why it launches in week 1 even though its results are consumed in week 3.
Crux path: **D2 → T2.1 → D5**.
Free path: **T0.x all startable today.**

**Advisor meeting agenda (one meeting unblocks everything):** D1 (referee), D2 (three ranging design
recs — say yes/no to each), preview D4.

---

## §6. Standardized evaluation protocol (use for EVERY result from now on)

- **Primary referee:** `outputs/phase2/sorcha_comparison_v5_masked.parquet` (707,670 tracklets).
  Secondary: pure-S3M variant. NEOMOD3 referee (T1.1) for all Tier-2/NEOMOD3 claims, reported
  alongside the S3M referee (D1-c). Benchmark v1/v2: quarantined pending D3.
- **Truth label:** `population == "NEO"`.
- **Bands:** |Δλ_antisun| ∈ {0–20, 20–40, 40–70, 70–110, 110–180}°.
- **Metrics:** best-F1 (+ completeness/contamination at that threshold) and AUC; full sky + per band.
- **Any trained combiner:** 50/50 train/test split, seed 0, metrics reported on the test half only.
- **Scoring-time flags (the F2 lesson):** any re-scoring must pass
  `support_mask_min=1, mask_radius_deg_per_day=np.inf` to match production columns, and must first
  demonstrate corr = 1.000 against the stored `P_NEO_vdp` on the same rows before any variant is
  trusted.

---

## §7. Closed — do NOT rebuild (and why)

| Item | Why closed |
|---|---|
| 1B sky-interp / continuous-mag / 5° denser grid | F2: +0.006 = noise; the +0.039 was a mask artifact |
| NEOMOD3 orbits into GMM training set (old attempt) | F6: changes support, not weights — provably null |
| Ranging score as raw weight | old fatal flaw; score is a class ratio (D2.3) |
| Benchmark v2 as strategic referee | F5: degenerate VDP column, numbers don't reproduce (pending T0.2) |
| "Physics, not a bug" as the whole mid-elongation story | F3: digest2 gets AUC 0.93 from identical observables |

---

## Logistics (unchanged)
Git: code (`src/`, scripts, docs) committed on the machine that made it; **no `Co-Authored-By`**;
user pushes. Parquets / .npz / Figures gitignored → scp/rsync.
Hyak↔Arnor doc sync: `scp neomod/docs/<f>.md arnor.astro.washington.edu:/astro/users/ds2004/vdp/docs/`.
Env note (Arnor): sklearn in `neofast_py310` is broken (GLIBCXX) — irrelevant for T0.x (pure numpy),
fix only if map work ever moves local.
