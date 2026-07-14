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

**STATUS (2026-07-08 — advisor away until Jul 13; decisions taken solo, documented for review):**
- **D1 DECIDED**: Kurlander set = referee/TEST set only (never a classifier ingredient); report both
  referees. Advisor sanity-check Jul 13.
- **D2 DECIDED v1**: full spec below (systematic grid, AR-restricted, uniform-in-ρ̇ + ρ² prior,
  L0/L1/L2 ablation ladder with geometry frozen via L1↔digest2 calibration). All knobs parameterised.
- **D3 (user)**: keep all benchmark parquets; T0.2 audit to be ready before the Mon Jul 13 advisor
  meeting so the three versions (v1 caps / v2 proportional / v3) can be explained with evidence.
- **D4**: clarified — it is only the paper-headline framing (NEOMOD3-first vs stack-first), decided
  by the 1A outcome; nothing to do now, advisor call later.
- **D5**: stays pre-registered as written (GO: rng F1 ≥ 0.839 or stack(VDP,rng) ≥ 0.862 on 40–70°;
  kill < 0.80 after one iteration). Blocks nothing now; kept because criteria fixed before results
  are what make solo results credible.
- **D6 DONE**: NEOCP cron restarted 2026-07-07 by user; verify accumulation in ~1 week.

### D1 — the referee for NEOMOD3 claims  ⟶ **DECIDED (2026-07-08, solo; advisor sanity-check Jul 13)**
The published, DOI-citable Kurlander et al. 2025 catalog (T1.1 Option A, verified byte-exact against
CANFAR) is the referee. **Role clarification (user-confirmed): it is a TEST set only** — we score its
tracklets with (a) the existing VDP maps, (b) digest2, (c) the new ranging term, and compare. It is
never an ingredient of any classifier: no GMM training, no density tables, no priors derive from it.
Classifier ingredients = NEOMOD3 model table (NEO numerator) + S3M census (non-NEO denominator) +
existing hybrid/GMM maps. Every result is reported on BOTH referees (S3M-drawn v5 parquet AND the
Kurlander NEOMOD3-drawn set) per D1-c. Original options kept below for the record.
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

### D2 — 1A ranging-term design  ⟶ **DECIDED v1 (2026-07-08, solo — advisor review Jul 13)**
Papers reviewed: **Farnocchia, Chesley & Micheli 2015** (systematic ranging / Scout),
**Spoto et al. 2018** (Admissible Region + MOV, OrbFit NEOCP scanner), **Solin & Granvik 2018**
(neoranger, statistical/MCMC ranging). Every knob below is an explicit parameter for later tweaking.
**Full concept-by-concept source map (section/equation/figure per paper): `docs/D2_detail.md`.**

**Method class: systematic (ρ, ρ̇) grid, NOT statistical/MCMC ranging.** neoranger costs ~1–10 min
per object (Solin & Granvik §4.4) — impossible at 10⁵–10⁶ tracklets. Scout/OrbFit's per-node cost is
the constrained least-squares attributable fit — but **our tracklets have exactly 2 detections: the
attributable A=(α,δ,α̇,δ̇) is determined EXACTLY (4 data = 4 params, χ²≡0), so no per-node fit
exists**. The whole computation is closed-form numpy over (tracklet × node); the posterior over
(ρ,ρ̇) is driven entirely by the prior/population structure — the very thing we upgrade.
(Extension knob: for ≥3-detection tracklets adopt Scout's constrained fit; not needed for v1.)

**v1 specification (engine = `P_NEO_rng`):**
1. **ρ grid:** log₁₀-uniform, ρ ∈ [0.01, 100] au, N_ρ = 64. (Farnocchia log spacing; 100 au = AR
   a_max and digest2's outer bound — do NOT truncate at 5 au or TNO/Trojan hypotheses vanish from
   the denominator. ρ<0.01 au dropped = Earth-satellite region, Spoto cond. 2.)
2. **ρ̇ sampling:** per-ρ-column admissible interval from the bounded-orbit energy condition
   (heliocentric E < −k²/(2·100 au); Spoto §2.1 cond. 1 — quadratic in ρ̇, closed form, vectorises),
   N_ρ̇ = 32 uniform inside it; node weight carries the interval length Δρ̇(ρ). Zero wasted nodes;
   every node is a bound orbit. Shooting-star cut H(ρ) ≤ 34.5 (Milani/Spoto).
3. **Node weight (prior):** w = [Δlog₁₀ρ · ln10 · ρ] (log measure, Farnocchia fn.4) × ρ^γ (spatial
   factor; **γ=2 default** per Farnocchia §3.2 f′_prior ∝ ρ²; γ ∈ {0,2,4} switch) × Δρ̇(ρ)/N_ρ̇
   (uniform in ρ̇ — Farnocchia's tested-best choice). **Jeffreys' prior explicitly rejected**
   (Farnocchia Table 3: p-values 10⁻⁵–10⁻⁶ for TRUE MBA solutions; Solin & Granvik report the same
   uniform-vs-Jeffreys pathology).
4. **H rides the ρ-grid:** H = V − 5log₁₀(ρ·r_helio) − Φ_HG(phase, G=0.15) — exactly the Kurlander
   input convention; observed filter mag → V via the same solar-color offsets our digest2 phase-3
   uses (documented knob). This is digest2's V↔d↔H coupling; no extra grid axis.
5. **Population factors:** numerator = NEOMOD3 4D table f(a,e,i,H) (`input_neomod3.dat`,
   multilinear interp). Denominator = empirical S3M histogram in (q,e,i,H) for MBA+TNO+Trojans
   (raw `.s3m` census, bin widths matched to NEOMOD3 axes: dH=0.25, da≈0.10, de=0.04, di=4°; light
   1-bin Gaussian smoothing — both knobs). Wagg ×0.80 MBA rescale = calibration-only option
   (AUC/ranking-neutral, threshold shift only).
6. **Score = class ratio (the old NEO_H.py fatal-flaw fix):**
   `P_NEO_rng = Σ w·f_NEO / (Σ w·f_NEO + Σ w·f_nonNEO)` — never a raw weight.
7. **Ablation ladder — one engine, three population settings (validation + the C2 experiment):**
   - **L0**: f ≡ 1 within the AR — geometric-only, ≈ Spoto et al.'s NEO "score".
   - **L1**: f = S3M for all classes — a **digest2 replication**. Engine validation: L1 must
     correlate strongly with stored `P_NEO_d2`. **All geometry knobs (γ, resolution, color table)
     are frozen by maximising L1↔digest2 agreement on the 40–70° band.**
   - **L2**: NEOMOD3 numerator ÷ S3M denominator — the classifier. **L2 − L1 = the isolated NEOMOD3
     effect with every knob frozen** — no post-hoc tuning on the result of interest (the solo-work
     discipline).
8. **Convergence:** rerun a 10k subsample at 128×64; require |ΔF1| < 0.002.
9. **Truth diagnostic (never scoring):** Kurlander rows carry true Range/RangeRate — compute
   Farnocchia-style p-values of the truth under the L2 posterior on the (ρ,ρ̇) grid; flags prior
   pathologies exactly as Farnocchia Table 3 flagged Jeffreys.
10. **Cost:** 707k × ~2k nodes ≈ 1.4×10⁹ closed-form node evals → chunked numpy, tens of minutes on
    one Hyak node; the 40–70° prototype band = minutes (Arnor-viable).

Original recommendation (kept for the record):
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

**T1.1 NEOMOD3 referee set (2C-a; needs D1 confirmed). — Hyak. FAST PATH FOUND (2026-07-08).**
- **Option A (adopt Jake Kurlander's public Sorcha run — Kurlander et al. 2025, AJ 170:99):**
  `https://epyc.astro.washington.edu/~jkurla/LSST_Sorcha_predictions/`. Verified from the paper +
  pilot files (`one_day_neo.h5` inspected on Arnor):
  - NEO input = **NEOMOD3** (orbits+H+albedo from the NEOMOD3 generator); MBA = **S3M × 0.80 Wagg
    rescale** (matches our 2B plan!); Trojans = Vokrouhlický+2024; TNO = CFEPS-L7 9 subpops;
    Hildas available (`hildamod/`). Population proportions = each model's absolute calibration →
    physically fair mixing.
  - Outputs are **all 5σ pre-linking detections** (SSP outcome as a `Linked` bool), 57 cols incl.
    RA/Dec (noisy + true), **RA/Dec rates**, per-visit MJD/filter/mags(PSF+trailed+σ)/SNR, full
    orbital elements + H_r per row, and **Range_LTC_km + RangeRate** (truth for validating 1A!).
  - Layout: per-population `infiles/` + `outfiles/` (s3m/outfiles = 2778 h5; neomod/outfiles = 6
    files + small-NEO d<10m subsample, weight 4.42). ~1.1B detections total. One night ≈ 470 NEO
    objects (~265 tracklet-able) → need the bulk files, not the one_day pilots.
  - Cadence caveat: **v3.4 baseline (survey start 2025-05-01), not our v5.0** — fine as a
    self-contained referee (antisun-relative maps are epoch-free), but per-band numbers are not
    directly comparable to our v5 run.
  - Remaining work: bulk-copy outfiles → Hyak → our Phase-1 tracklet builder → Phase-2 VDP →
    Phase-3 digest2 → `sorcha_comparison_neomod3ref.parquet`. **No Sorcha run needed.**
- **Access + canonicality VERIFIED (2026-07-08, Arnor):**
  - Epyc's disk is NFS-mounted on Arnor: `/astro/users/jkurla/public_html/LSST_Sorcha_predictions/`
    — direct read, fast (200k rows / 0.4 s). Hyak needs an rsync (storage confirmed available).
  - Sizes: s3m outfiles **479 GB** (2778 h5), neomod outfiles ~3.8 GB, trojanmod 21 GB, cfeps
    8.9 GB, hildamod 4.8 GB. Schema note: NEO files carry KEP elements (a, ma), others COM (q, t_p).
  - **CANFAR cross-check COMPLETE (byte-exact, via DOI 25.0062 README + listings):**
    Epyc `neo_output_1.h5` (3,032,238,888 B) **= CANFAR `neo/outputs/large_neo_output.h5`**
    (2.82 GB exact); `small_neo_output.h5` = 887.16 MB exact; s3m `0_0.h5` = 435.51 MB exact;
    s3m outfiles = CANFAR `mba/outputs/{0-234}_{0-12}.h5` (2,777 files; ignore stray `old_tr`).
    The `_022` / `.csv` files on Epyc are dev leftovers NOT in the DOI dataset — do not use.
    Repo-root `s3m_colors.h5`/`s3m_orbits.h5` on Arnor = CANFAR `mba/inputs` (byte-exact).
    CANFAR also publishes `config.ini` + `baseline_v3.4_10yrs.db` → full reproducibility.
  - **Paper-number check:** large-NEO file: 112,855 linked ≈ paper's 1.1E5 (d≥10m) ✓; MJD
    60796.0–64448.4 = exactly 10.0 yr from 2025-05-01 (v3.4 start) ✓; 406,758 NEOs observed ≥1×,
    6.7M detections. README col spec: standard Sorcha outputs, `object_linked`→`Linked`, no
    date_linked_MJD; NEO ObjID in small file = `{ObjID1}N{ObjID2}`.
  - **One residual oddity:** canonical `small_neo_output.h5` has 42,348 Linked objects vs the
    paper text's 3,026 raw (×4.425 → paper's 1.4E4). Inconsistent with Table 5 arithmetic —
    **exclude d<10 m from our benchmark for now** (H≳28 is beyond NEOMOD3's debiased core anyway);
    optionally ask Jake casually.
- **Option B (fallback, only if Jake objects to reuse):** run it ourselves as originally planned
  (`neomod3_sampler` NEO input + S3M non-NEO, v5 config, days of wall time).
- **Accept:** baseline VDP + d2 numbers on it are sane; this ALSO immediately measures how much the
  current S3M-trained classifiers degrade on debiased NEOs — a paper number by itself.
- **Leakage rule:** Jake's NEOMOD3 draw is EVAL-ONLY. Any NEOMOD3-based training/reweighting (1A,
  2A/2B) must use our own independent draw via `neomod3_sampler`.

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

## §8. CONTINUOUS QA HARNESS (the Željko protocol)

**Mechanism:** one living notebook `notebooks/qa/1A_engine_qa.ipynb`, sections QA0–QA5 matching the
stages below. Every check runs on the **frozen QA subsample** `outputs/qa_subsample.parquet`
(~20k tracklets from the v5 eval parquet, seed 0, stratified by class × elongation band) so plots are
comparable run-to-run. Plots saved to `Figures/qa/` with stage-prefixed names. The engine exposes all
intermediates (`debug=True` → per-node elements, H, weights, per-class sums) so the notebook never
reimplements pipeline math. **Gate rule: no stage's code advances until its checks pass and the plots
are in the notebook.** The notebook doubles as the Jul-13 advisor deck.

| Stage | Check | What the plot shows | Pass looks like |
|---|---|---|---|
| **QA0 inputs** | QA0a NEOMOD3 marginals | 1D a,e,i,H marginals + (a,e) density at H=18 vs H=25 | Shapes match Nesvorný+24 figs; (a,e) shape CHANGES with H (the whole point of NEOMOD3) |
| | QA0b S3M denominator | per-class q,e,i,H marginals + bin-occupancy histogram | MBA q≈1.8–3.3, Trojan a≈5.2, TNO q≳28; occupancy justifies smoothing choice |
| | QA0c NEOMOD3/S3M-NEO ratio map | ratio in (a,e) and (H) | Shows WHERE 1A can differ from digest2; not flat |
| **QA1 geometry** | QA1a truth round-trip | recovered (a,e,i,q) at TRUE (ρ,ρ̇) vs truth columns, 1:1 line | <1% scatter — validates Earth state, light-time, unit vectors, rates |
| | QA1b Farnocchia Fig.-1 panels | e,q,i,H contours over (ρ,ρ̇) for 4 exemplar tracklets (NEO; MBA@antisun; MBA@60°; TNO) | Qualitatively matches F15 Fig. 1 (H contours ~vertical; AR boundary shape). **The eyeball check** |
| | QA1c AR sanity | admissible fraction vs ρ per elongation, H=34.5 line overlaid | AR closes at large ρ; no admissible nodes with E>0 |
| **QA2 prior (L0)** | QA2a weight maps | w(ρ,ρ̇) for γ=0/2/4, same 4 exemplars | Smooth; γ visibly shifts mass outward |
| | QA2b L0 score by class | score distributions NEO/MBA/TNO/Trojan on 40–70° band | Some separation already (S18's score works); record AUC as floor |
| | QA2c truth p-values | F15-Table-3-style p-value histogram of TRUE (ρ,ρ̇) under posterior, per γ | ~Uniform. Spike at 0 = pathological prior (how F15 convicted Jeffreys) |
| **QA3 L1=digest2 replication** | QA3a L1 vs P_NEO_d2 | scatter + corr + binned calibration curve | High corr; monotone calibration |
| | QA3b knob freeze | L1↔d2 agreement over (γ, resolution, color table) | Pick max, FREEZE, record values here |
| | QA3c disagreement anatomy | L1−d2 residual vs elongation/mag/rate | Interpretable (obs-error machinery, binning), not structured by class |
| **QA4 L2=result** | QA4a L2 vs L1 by class | scatter colored by truth | NEOs move up / MBAs down where NEOMOD3≠S3M |
| | QA4b headline table | F1/AUC per band, BOTH referees, vs pre-registered D5 bar | Judged against D5 (0.839 / 0.862), no re-litigation |
| | QA4c the money plot | (a,e,i,H) location of tracklets whose class flipped L1→L2 | Matches QA0c ratio map — result explained by the model difference |
| | QA4d reliability | P(NEO) vs empirical NEO fraction, both referees | Calibrated-ish; note Wagg-0.80 effect |
| **QA5 referee pipeline** | QA5a tracklet QA | dt distribution, nightly counts, class mix vs Kurlander paper Fig. 2 | Pairs ~30 min; ~3×10⁵ det/night |
| | QA5b **antisun density maps** | (v_λ,v_β) per-class density at antisun + 60°, side-by-side with same maps from our v5 run | Same physics, same shapes — Željko's canonical check |
| | QA5c baseline scores | P_NEO_vdp + P_NEO_d2 distributions and per-band F1 on Kurlander tracklets vs v5 run | Sane baselines before any L2 claim |

## §9. WORK LOG + DETAILED STEP PLANS (living section — update as work proceeds)

**Design reference:** `docs/D2_detail.md` (paper-by-paper source map for every D2 choice).

### Step A — ranging engine + L0/L1 prototype (Arnor; me; target: before Jul 13)
- **Files:** `src/ranging_engine.py` (module) + `notebooks/qa/1A_engine_qa.ipynb` (QA0–QA3).
- **Module layout:** `build_grid(A, obs_state)` (log-ρ columns, per-column admissible ρ̇ via the
  energy quadratic) → `elements_from_nodes()` (heliocentric state → a,e,i,q; vectorised) →
  `H_from_nodes()` (V from filter mag via phase-3 color table; HG G=0.15) → `node_weights(γ)` →
  `class_score(pop_tables, level=L0|L1|L2)`. `debug=True` returns all intermediates for QA.
- **Inputs (all on Arnor):** `input_neomod3.dat`; `.s3m` census files (repo root); v5 masked parquet
  (`ra0/dec0/mjd0_utc/ra1/dec1/mjd1_utc/mag0/filter0` per tracklet); astropy Earth state.
- **Order:** QA0 input checks → geometry + QA1 (truth round-trip on benchmark/Kurlander pilot rows,
  Fig.-1 panels) → weights + QA2 → S3M tables + L1 + QA3 (freeze knobs) → L0/L1 scores on the
  40–70° band of the v5 parquet.
- **Output:** `outputs/rng_prototype_L0L1_40_70.parquet`; QA0–QA3 pages in the notebook.
- **Gate to Step D:** QA3a corr healthy + knobs frozen and recorded in §8 table.

### Step B — T0.2 benchmark audit + T0.1 stack banking (Arnor; me; target: before Mon Jul 13)
- **T0.2 (Monday ammunition):** extract exact filters/flags from `benchmark_v5_normalisation_s3m`
  + v2 notebooks; re-score a one-map subsample with `support_mask_min=1, mask_radius=np.inf` per §6;
  identify why raw-parquet columns are degenerate (prime suspect: F2-class mask flag at scoring
  time); produce a one-page summary of the three benchmark versions (v1 caps / v2 proportional /
  v3 cases) with the honest numbers for each. → appended under D3.
- **T0.1:** `src/stack_vdp_d2.py` per §4 T0.1 (logit-logistic, 50/50 seed 0, §6 protocol);
  outputs `outputs/stack_scores_{hybrid,s3m}.parquet` + ROC/F1 figure + paper table. Accept:
  reproduces F1 0.865/0.879 ±0.002.

### Step C — Kurlander referee build (Hyak; user rsync + my adapter; target: as rsync lands)
- **C1 rsync — DONE (2026-07-08, user):** `/astro/users/jkurla/public_html/LSST_Sorcha_predictions/`
  → `/mmfs1/gscratch/dirac/ds2004/kurlander2025/` (dirac allocation). Steps C2–C4 now unblocked.
  First action on Hyak: verify byte counts/row counts match the Arnor-side numbers logged in T1.1
  (large_neo_output.h5 = 3,032,238,888 B; 112,855 linked NEOs; s3m outfiles = 2,777 files).
- **C2 schema adapter (me, can prototype on Arnor over NFS before rsync completes):** map Jake's
  57-col outputs → our phase-1 inputs: `RA_deg/Dec_deg` (noisy astrometry) → ra/dec;
  `fieldMJD_TAI` → time (convert TAI→UTC as our pipeline expects); `trailedSourceMag`+`optFilter` →
  mag/filter; `ObjID` → object id; truth: NEO ⇔ q<1.3 with q=a(1−e) (KEP files) or q direct (COM);
  `Linked` kept as a column (enables linked-only sensitivity cut). **Never read** Range/RangeRate/
  elements during scoring — truth/diagnostic columns only.
- **C3 tracklet build (Hyak):** same-night detection pairs, same rules as v5 phase-1 (min separation
  as in production config); exclude small-NEO file (paper-number oddity, see T1.1 card) and Hildas
  (v1). → `tracklets_kurlander.parquet`.
- **C4 scoring (Hyak):** Phase-2 VDP (667 maps, §6 flags) + Phase-3 digest2 →
  `sorcha_comparison_neomod3ref.parquet` (v5-compatible schema). QA5 checks in the notebook.

### Step D — L2 on both referees → D5 (Hyak run, Arnor evaluation; after A+C)
- Run matrix: {L0, L1, L2} × {v5 S3M referee, Kurlander NEOMOD3 referee}, full sets (L2 first on
  the 40–70° band, then full sky). Fixed knobs from QA3b — no retuning after this point.
- Evaluate per §6; apply the pre-registered D5 criterion; QA4 pages; assemble the Jul-13 deck
  (QA notebook + stack result + T0.2 one-pager + this plan).

### Log
- 2026-07-08 (**RESULT 3 — the decisive symmetric test, CONFIRMED**, Opus).
  `notebooks/qa/qa5_kurlander_referee.py` → `outputs/qa5_kurlander_scores.parquet`.
  Eval: 120k tracklets, NEO frac forced to 0.293 (v5-matched, for F1 comparability). Knobs frozen
  by QA3 (γ=2, 128×64). Bootstrap B=300.

  **The pre-registered symmetric flip is significant:**
  | referee (truth NEOs from) | L2 − L1 (F1) | P(dF1>0) |
  |---|---|---|
  | v5 (**S3M**-drawn) | **−0.0021 ± 0.0009** | 0.01 |
  | Kurlander (**NEOMOD3**-drawn), per_H_match | +0.0008 ± 0.0004 | 0.97 |
  | Kurlander (**NEOMOD3**-drawn), absolute | **+0.0028 ± 0.0004** | 1.000 |

  The classifier prefers whichever NEO model generated the truth. Since only the NEO population
  differs between referees (both use the same S3M non-NEO, which is also our denominator), this is a
  clean one-variable experiment.

  **Per-elongation ramp (the robust signal), L2_perH − L1:** monotone across five bins,
  −0.0026±0.0005 (0–20°, antisun) → −0.0006 (20–40°) → +0.0017±0.0007 (40–70°) →
  +0.0061±0.0016 (70–110°) → **+0.0116±0.0035 (110–180°, sunward)**. Crosses zero between 20–40 and
  40–70. Physically coherent: |Δλ|→180° looks sunward where Atens/IEOs dominate — exactly the
  low-a region QA0c flagged as NEOMOD3's biggest enrichment over S3M.

  **Mechanism identified — 1A alone is hamstrung; use 1A+2B.** `absolute` is positive in EVERY band
  (antisun +0.0027±0.0006, sunward +0.0119±0.0033). The antisun mag-split proves why:
  dF1(per_H_match) goes −0.0011 → −0.0055 with faintness while dF1(absolute) goes −0.0005 → **+0.0067**.
  This is exactly QA0c's prediction: per_H_match inherits S3M-NEO's H≈25 cutoff, zeroing the NEO
  numerator for faint tracklets. **Recommendation: NEOMOD3 must be run with its own N(H)
  (`neomod3_norm='absolute'`), i.e. items 1A and 2B are not separable in practice.**

  ⚠️ **CONFOUND (must control before publishing the faint-end claim).** The Kurlander referee's
  non-NEO population is the SAME S3M (Wagg-scaled) that our denominator uses, with the same
  magnitude cutoff. So "no faint MBAs" is true *in this referee by construction*, and the faint-end
  gain of `absolute` is **partly circular** — it may not transfer to real Rubin data, where faint
  MBAs exist. The **elongation ramp is NOT affected** (it is driven by NEO orbital structure, not
  the MBA cutoff) and is therefore the robust, quotable result. Required control: repeat with a
  magnitude cut where the S3M non-NEO census is complete, and/or an extrapolated faint MBA N(H).

  ⚠️ **Referee construction caveat:** we took ALL NEO tracklets (1.85M) but only a 1/200 object-batch
  subset of MBA → the "natural" NEO fraction (0.419) is an ARTIFACT, not the Rubin ratio. Only the
  forced-0.293 eval set is meaningful. A full-population build fixes this.

  **Honest hierarchy of effect sizes.** Engine vs digest2 = **+0.067 F1** (RESULT 1). NEOMOD3 vs S3M
  prior = **+0.003 F1** aggregate (up to +0.012 sunward). *The classifier's implementation matters
  ~20× more than the debiased population model at Rubin depths* — except in the sunward/faint regime
  where the debiased model is measurably, significantly better. That is the paper's real story.

  STILL MISSING for the full 4-way table: digest2 and VDP-map scores on the Kurlander referee.

  **VENUE DECISION (corrected).** The decisive test ran on Arnor because it was a ~5 min job on a
  random object subset — that was fine. But the *full-population* build was briefly queued on Arnor
  and **killed**: it reads 479 GB across the **shared Epyc NFS mount**, which is slower and
  inconsiderate to other astro users, and defeats the purpose of the rsync that made the data local
  to Hyak. Two jobs are genuinely Hyak-only / Hyak-appropriate and are handed off in
  **`docs/HYAK_HANDOFF_2026-07-08.md`**:
    (A) **digest2 on the Kurlander referee** — the `digest2` binary is NOT installed on Arnor; this
        is the last cell of the 4-way table and answers "does the engine still beat digest2 when the
        truth is NEOMOD3-drawn?"
    (B) **full-population tracklet build** (40-way Slurm array, ~5 min) — removes the 1/200 MBA
        subset caveat and gives the true natural population ratio.
  I cannot launch these: `ssh ds2004@klone` from Arnor → `Permission denied` (needs password/2FA).
  New artifacts: `pipeline/kurlander/score_digest2.py`, `pipeline/kurlander/merge_shards.py`,
  `pipeline/kurlander/slurm_build_tracklets.sbatch`, `outputs/kurlander/referee_eval.parquet`
  (the exact 120k rows scored on Arnor, so digest2 lands on identical tracklets).
- 2026-07-08 (**Step C: Kurlander referee built — far cheaper than planned**, Opus).
  `pipeline/kurlander/build_tracklets.py` (+ `slurm_build_tracklets.sbatch`, `merge_shards.py`),
  `notebooks/qa/qa5_kurlander_referee.py`.
  **Cost re-estimate kills the "long Hyak job" premise.** Kurlander files are partitioned by OBJECT
  BATCH, not time, so (a) any subset of chunks is a uniformly random object subsample, and (b) 69%
  of object-nights already have ≥2 detections *within one chunk* → tracklets form with no cross-file
  join. Also: the naive per-group Python loop cost 36–45 s/file (→ ~30 h for 2778 files); the
  **vectorised groupby is 1.9 s/file, bit-identical (max diff 1e-13)** → full scan ≈ 88 min
  single-threaded, ~5 min on a 40-way array. Built on Arnor over Epyc NFS in ~2 min:
  **1,852,257 NEO tracklets** (complete canonical file) + 2,565,532 MBA (random subset) + Trojan/TNO.
  **Why this referee is a clean experiment:** its non-NEO population is the SAME S3M our denominator
  uses; only the NEO population changes (NEOMOD3 ← S3M). Exactly one variable moves.
  **Methodological trap recorded:** F1 is prevalence-dependent → the eval set is subsampled to the
  v5 NEO fraction (0.293) so F1 is comparable across referees; AUC (prevalence-independent) is also
  reported on the natural mix. Never compare F1 across referees with different NEO fractions.
- 2026-07-08: D1/D2 decided solo (see §3 status block); D2_detail.md written; QA harness designed
  (§8); step plans A–D added. Rsync to Hyak started by user. NEOCP cron restarted 2026-07-07.
- 2026-07-08 (Step A build, Opus): ranging engine written + validated. Artifacts:
  `src/ranging_engine.py` (build_grid/elements/H/weights/class_score/score_tracklets, all knobs),
  `src/build_population_cache.py` → `outputs/pop_cache_wide.npz` (S3M neo/mba/trojan/tno on wide
  (H,log10a,e,i) grid; MBA 13.9M rows), `notebooks/qa/qa1a_truth_roundtrip.py`. Fixed NEOMOD3 file
  path (symlink NEOMOD3/→root). Installed jplephem.
  **Three latent bugs in old NEO_H.py found + fixed:** (1) equatorial vs ecliptic inclination
  (~23°), (2) barycentric vs heliocentric state (Sun offset ~1.5e6 km), (3) missing observer diurnal
  velocity (~0.46 km/s → ~3% in a). Also fixed the two original design flaws: H now rides the grid
  (was hardcoded 19.0), score is a class ratio (was raw weight sum).
  **QA1a truth round-trip PASS:** on 861 Kurlander NEO detections, median |rel| error a=0.007%,
  e=0.007%, i=0.003%, q=0.003% (was 2.5% before the 3 fixes). Geometry validated to ground truth.
  **Engine smoke test PASS** (800 stratified v5 tracklets, balanced NEO/MBA):
  L0 geometric AUC=0.777 (corr_d2 0.54); **L1 S3M AUC=0.972, corr_d2=0.905** (NEO medP 1.000, MBA
  0.001) vs stored digest2 AUC 0.942. L1↔digest2 replication confirmed (QA3 gate essentially met on
  first try). Speed: ~0.5 s / 800 tracklets after table load; earth-state ephemeris is the bottleneck
  (~tens of min for full 707k, matches estimate). IERS polar-motion warning for post-2025 dates =
  arcsec-level, ignorable.
  **L2 smoke test PASS** (same 800 sample): AUC=0.970 (≈L1 0.972); L2−L1 mean +0.003, corr 0.9987;
  effect asymmetric in the right direction — true NEOs +0.0057, MBAs +0.0002 (NEOMOD3 lifts genuine
  NEOs). Small on a balanced sample as expected; the mid-elongation overlap (40–70° band) is where it
  should matter — that is QA4b.
  NEXT: QA notebook then QA3 knob-freeze, QA4 L2-vs-L1 on the 40–70° band vs D5 bar.
- 2026-07-08 (QA harness built + executed, Opus): `notebooks/qa/qa_lib.py` (plotting fns),
  `notebooks/qa/build_notebook.py` → **`notebooks/qa/1A_engine_qa.ipynb` executed, 8 figures embedded,
  0 errors** (kernel `neofast_py310`; installed nbformat/nbclient/ipykernel). Figures also in
  `Figures/qa/`. Fixed: earth_observer_state now accepts numeric MJD; engine default cache path is
  module-absolute (`DEFAULT_CACHE`) so it works from any CWD.
  **QA verdicts (all pass):**
  - QA1a geometry round-trip: 0.007% (already logged).
  - QA1b (Farnocchia Fig-1 eyeball): **strong match** — H contours ~vertical (V↔d↔H coupling), AR
    closes at large ρ, and the **TNO exemplar shows Spoto's second admissible component at ρ≈40 au**
    (distant-object two-component AR, q>28). Geometry qualitatively correct across NEO/MBA/TNO.
  - QA2b: L0 geometric score ranks NEO(0.64) > MBA(0.44) > Trojan(0.29) > TNO(0.00) — correct.
  - QA2c: **no pathological spike at p=0** (unlike Jeffreys, Farnocchia Table 3); prior healthy.
    High-skew toward p~0.8 flagged as a γ-tuning candidate for QA3.
  - **QA0c — strategically important finding:** S3M-NEO's N(H) peaks at H≈24–25 and cuts off while
    NEOMOD3 rises to 28. The default L2 per-H normalisation isolates orbital *shape* (clean 1A test)
    but structurally inherits S3M's faint cutoff → **1A alone is a modest orbital-shape effect; the
    large faint-end NEOMOD3 gain lives in 2B (N(H) renormalisation).** This confirms the 1A/2B split
    is the right decomposition and tempers expectations for 1A alone (consistent with the small
    L2−L1 in the smoke test). Actionable: when 1A is evaluated (QA4/D5), also plan a 2B variant
    (`neomod3_norm='absolute'` — a knob to add) to expose the faint-end effect.
- 2026-07-08 (**QA3 + QA4 — the two headline results**, Opus). Artifacts:
  `notebooks/qa/qa3_knob_sweep.py`, `notebooks/qa/qa4_L2_eval.py`, `outputs/qa3_knob_sweep.csv`,
  `outputs/qa4_band_scores.parquet`, `Figures/qa/qa3_knob_sweep.png`, `Figures/qa/qa4_L2_vs_L1.png`.

  **BUG FOUND + FIXED (RA-rate convention).** The v5 parquet's `mean_dra` is the RAW α̇, **not**
  α̇·cos δ. Verified on 2-detection tracklets against the wrap-safe finite difference of (ra0,ra1):
  raw hypothesis corr = 0.999990 (median|rel| 2.6e-3) vs cos-δ hypothesis corr = 0.939 (median|rel|
  3.1e-2 ≈ 1−cos(13°) = median |dec|). The engine expects the proper-motion convention (α̇·cos δ),
  which Kurlander's `RARateCosDec_deg_day` supplies — which is why QA1a never caught it. Engine now
  takes an explicit `dra_cosdec` flag (declare, don't guess); v5 callers pass `dra_cosdec=False`.
  Fixing it improved L1 band F1 0.893 → 0.905.

  **KNOB FREEZE (QA3).** γ: the literal protocol ("max L1↔digest2 agreement") picks γ=0
  (Spearman 0.80), but γ=0 is a *worse classifier* (AUC 0.965 / F1 0.890) than γ=2 (AUC 0.971 /
  F1 0.905) — maximising agreement means inheriting digest2's suboptimality. **Adopted γ=2 as a
  pre-registered literature default** (Farnocchia 2015 §3.2 ρ² spatial factor, with f_pop supplied
  separately = exactly our structure), independently confirmed by AUC. The γ-vs-agreement trend is
  reported as a diagnostic (digest2's effective prior looks flatter, cf. Farnocchia §3.3
  ρ^{2−5η} ≈ flat), **not** used to select. Grid: 128×64 (64×32 within 0.001 AUC; ΔF1 non-monotone
  at ~0.004 = subsample noise). p-value KS ≈ 0.52 for all γ → not γ-discriminating; the pass
  condition (no spike at 0) stands, but this diagnostic is weaker than hoped.

  **DIGEST2 BASELINE AUDIT (fairness).** `run_digest2_comparison_gmm.py` feeds digest2 exactly **2
  observations synthesised from the same (position, rate)** we use (`ra1 = ra0 + dra·dt`, no cos δ
  divide — independently confirming the raw-α̇ convention). So digest2 has *identical information*;
  we are not advantaged. `P_NEO_d2` reaching 1.98 is a handful of parse glitches (frac>1 ≈ 5e-5, all
  NEO); AUC is identical clipped (0.9303) → the column is a sound digest2 score. Residual caveats:
  digest2's score is quantised at 0.01 (inherent handicap), and our filter→V colour correction is
  still **unapplied** (a knob that should only help us further).

  ### RESULT 1 — the ranging ENGINE beats digest2 (20k tracklets, 40–70° band, knobs frozen)
  | scorer | F1 | AUC |
  |---|---|---|
  | **L1 (our engine, S3M prior)** | **0.9007** | **0.9672** |
  | L2 (NEOMOD3, per_H_match) = 1A | 0.8985 | 0.9650 |
  | L2 (NEOMOD3, absolute) = 1A+2B | 0.8990 | 0.9653 |
  | digest2 (reference) | 0.8334 | 0.9241 |

  **L1 − digest2 = +0.0673 F1 / +0.0431 AUC**, using the *same S3M population* and the *same two
  detections*. **The mid-elongation gap was never "orbit-space beats velocity-space" — it was a
  digest2 *implementation* gap** (coarse binning, quantised score, 2011-era machinery). This closes
  and reverses the −0.080 band gap (F3) that motivated the entire program. **D5 pre-registered bar
  (band F1 ≥ 0.839): cleared by the ENGINE alone (0.9007) → GO.**

  ### RESULT 2 — NEOMOD3 cannot win on an S3M referee (the trap, demonstrated)
  L2 − L1 = **−0.0022 F1** (per_H_match) / −0.0017 (absolute). NEOMOD3 *does* lift true NEOs
  (mean +0.0081 vs +0.0001 for non-NEOs) and QA4c confirms the flips land exactly where QA0c
  predicted (near q=1.3 and low-a Atens, 3.7% of tracklets) — but it loses AUC. **This is expected
  and is the referee trap (F7/D1) made empirical:** the v5 truth is drawn from S3M, so an S3M prior
  is optimal *by construction* and no better model of reality can win. The result is therefore
  **not evidence against NEOMOD3** — it is a validation that the harness detects prior/truth
  mismatch.
  → **Sharpened, symmetric, pre-registered prediction for the Kurlander referee (rsync landed):**
  **L1 > L2 on the S3M-drawn v5 referee (shown: −0.002) AND L2 > L1 on the NEOMOD3-drawn Kurlander
  referee.** A two-sided flip is far stronger evidence than a one-sided gain, and it is exactly what
  D1's dual-referee rule was designed to deliver. If the flip does not appear, NEOMOD3's orbital
  shape genuinely does not help at Rubin depths and C2 falls back to 2B/2D.

  **Strategic consequence (feeds D4).** The paper's headline candidate is now RESULT 1: a modern
  vectorised systematic-ranging classifier beats digest2 at mid-elongation on identical inputs;
  VDP remains the antisun winner; the stack should now be VDP + ranging (not VDP + digest2).
  NEOMOD3 becomes the *prior-fidelity* demonstration (RESULT 2's symmetric flip), not the headline.
  NEXT: Step C (Kurlander referee build on Hyak) is now the critical path and the decisive
  experiment. Also queued: filter→V colour correction knob; re-run the C1 stack with L1 replacing
  digest2; full-band (221k) confirmation of RESULT 1.

## Logistics (unchanged)
Git: code (`src/`, scripts, docs) committed on the machine that made it; **no `Co-Authored-By`**;
user pushes. Parquets / .npz / Figures gitignored → scp/rsync.
Hyak↔Arnor doc sync: `scp neomod/docs/<f>.md arnor.astro.washington.edu:/astro/users/ds2004/vdp/docs/`.
Env note (Arnor): sklearn in `neofast_py310` is broken (GLIBCXX) — irrelevant for T0.x (pure numpy),
fix only if map work ever moves local.
