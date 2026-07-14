# D2 in detail — the 1A ranging-term design, with paper-by-paper source map

**Purpose:** every concept in the D2 v1 spec (NEOMplan §3), traced to the exact section / equation /
figure / table of the source papers, so it can be looked up, challenged, or explained by another
agent without re-reading everything. Papers: **F15** = Farnocchia, Chesley & Micheli 2015
(*Systematic ranging and late warning asteroid impacts*, Icarus 258, arXiv:1504.00025);
**S18** = Spoto et al. 2018 (*Short arc orbit determination and imminent impactors in the Gaia era*,
A&A, arXiv:1801.04004); **SG18** = Solin & Granvik 2018 (*Monitoring NEO discoveries for imminent
impactors*, A&A 616, arXiv:1805.01315); **K19** = Keys et al. 2019 (*digest2*, PASP 131:064501).

---

## 0. The one-paragraph idea

A 2-detection tracklet fixes the sky position and motion — the **attributable**
A = (α, δ, α̇, δ̇) — but leaves topocentric range ρ and range-rate ρ̇ unconstrained (F15 §2; S18
§2.1, their Eq. 1; K19 §2 first paragraph). Every choice of (ρ, ρ̇) plus the observer's known
heliocentric state converts to one heliocentric orbit (a, e, i, …) and, through the apparent
magnitude, one absolute magnitude H. **Systematic ranging** evaluates a grid over (ρ, ρ̇) and asks,
at each node, "what kind of object would this be?" digest2 answers with S3M population counts in
(q, e, i, H) bins (K19 §3.3). Our 1A term is the same machine with the NEO side of the population
replaced by the debiased NEOMOD3 model, and the output defined as a class ratio
P(NEO) = Σw·f_NEO / Σw·(f_NEO + f_nonNEO).

---

## 1. Concept-by-concept source map

| # | Concept in our v1 spec | Where it comes from | Exact location |
|---|---|---|---|
| 1 | Attributable A=(α,δ,α̇,δ̇); ρ,ρ̇ unconstrained | F15, S18 | F15 §2 (p.4, citing Milani & Knežević 2005); S18 §2.1 Eq. (1) |
| 2 | Systematic raster over (ρ,ρ̇) vs statistical/MCMC sampling | F15, SG18 | F15 §1–2 (Chesley 2005 heritage, fn.3: Tholen & Whiteley 2002); SG18 §3.1 (review of both families) |
| 3 | Why NOT Monte-Carlo ranging at our scale | SG18 | SG18 §4.4: ranging = 1–7 min per observation set + ~1.5 min propagation, per OBJECT. We have 10⁵–10⁶ tracklets |
| 4 | Per-node constrained least-squares fit of A (the expensive part) | F15, S18 | F15 §2 Eq. (1) (Q = νᵀWν, minimised at fixed ρ,ρ̇); S18 §2.1 "doubly constrained differential corrections" |
| 5 | **Our shortcut: 2 detections ⇒ A is exact, χ²≡0, no fit** | ours | 4 observables = 4 attributable params. The F15/S18 fit machinery exists for ≥3 obs with noise; for a 2-det tracklet residuals vanish identically at every node → whole grid is closed-form |
| 6 | Grid must contain the Admissible Region (bound orbits) | F15, S18 | F15 §2 (citing Milani et al. 2004); S18 §2.1 condition 1: E_helio < −k²/(2a_max), **a_max = 100 au** |
| 7 | Not-an-Earth-satellite exclusion | S18 | S18 §2.1 condition 2, sphere-of-influence radius R_SI ≈ 0.010044 au |
| 8 | Shooting-star cut H ≤ 34.5 | S18 | S18 §2.1 ("shooting star limit", citing Milani et al. 2004) |
| 9 | Log-ρ sampling for near-field resolution | F15, S18 | F15 §3 footnote 4; S18 Table 1 (log₁₀ρ grid when object may be near / NEO score > 50%) |
| 10 | Log-measure factor (density × ρ·ln10 on a log grid) | F15, S18 | F15 fn.4 ("probability density has to be multiplied by ρ"); S18 Appendix A.3: det Mσ = log(10)·ρ |
| 11 | AR-restricted ρ̇ interval per ρ column (no wasted nodes) | S18 (vectorised by us) | S18 §2.1: "we check the value of the heliocentric energy for each grid point, and we discard those not satisfying condition 1." We solve the quadratic-in-ρ̇ energy condition analytically per column instead of discarding |
| 12 | Grid sizes precedent | S18 | S18 Table 1: 50×50 (one AR component), 100×100 (two components); two-step densification §2.1. We use 64×32 + convergence check — coarser is justified because item 5 removes all χ² structure from the plane |
| 13 | Posterior = error term × prior | F15 | F15 §3: f_post(ρ,ρ̇) ∝ f_err(ν(ρ,ρ̇))·f_prior(ρ,ρ̇). For us f_err ≡ const (item 5) ⇒ posterior = prior structure only |
| 14 | **Jeffreys prior REJECTED** | F15, SG18 | F15 §3.1 Eq. (2), §6.6 + Table 3 (p-values 10⁻⁵–10⁻⁶ for TRUE solutions of MBAs 2015 CV, 2015 BU92), §7 ("we avoid the use of Jeffreys' prior"); SG18 §4.2 (TC3/AA first tracklet: Jeffreys 81%/60% vs uniform ~10⁻¹⁰ — prior sensitivity) |
| 15 | Uniform-in-ρ̇ prior is the tested-best simple choice | F15, SG18, K19-era Scout | F15 §3.3 + §7 ("uniform … good compromise"; ρ^{2−5η} argument, η=0.35–0.47 → exponent −0.35…0.25 ≈ flat); SG18 §3.1 ("Farnocchia et al. (2015a) … concluded that a uniform prior produces the best results") |
| 16 | **Population-model prior — the template for our classifier** | F15 | F15 §3.2: f_pop(q,e,i,H) derived from S3M (Grav et al. 2011); computed iteratively with geometric factor **f′_prior ∝ ρ²** ("ρ³ and ρ̇ are uniformly distributed"). Our γ=2 default IS this factor; our upgrade = replace the NEO part of f_pop with NEOMOD3 |
| 17 | Population-prior caveat (Poisson/finite-sample risk) | F15 | F15 §7: "a population model is … affected by Poisson statistics errors … small impactors might have peculiar orbits deemed unlikely." Fine for THEIR goal (never miss an impactor); for OUR goal (classification) the population density is not a nuisance prior — it IS the signal, same as digest2 |
| 18 | Weighted-sum output over grid nodes | F15 | F15 §4: P = Σ w_ij·p_ij with w_ij ∝ f_post. We replace the impact indicator p_ij with class densities → class ratio |
| 19 | Class-score precedent (NEO/MBO/DO/SO probabilities) | S18, SG18 | S18 §2.1 (score definition; **NEO ⇔ q < 1.3 au**, our truth label too) + §4.5 (2017 AE21: score's operational value); SG18 §4.4 + Fig. 7 (NEO-class probability per tracklet). Our **L0** rung ≈ S18's score (uniform density within AR) |
| 20 | Rigorous Jacobian alternative to ad-hoc priors | S18 | S18 §3 + Appendix A: propagate the residual PDF through the MOV with det Mμ (Eq. 8) and det Mσ (A.3) — "no a priori assumption". Kept as a documented knob; for 2-det tracklets their χ-weighting degenerates and mainly the measure terms survive, which we already carry |
| 21 | H rides the ρ grid (V↔d↔H coupling) | F15, K19 | F15 §2 ("if the observations contain photometric measurements we also compute the absolute magnitude for each grid point") + Fig. 1 bottom-right panel (H contours over the (ρ,ρ̇) plane); K19 §3 (digest2 bins population in (q,e,i,**H**)) |
| 22 | Photometry NOT used as a residual/likelihood term | F15 | F15 §3 (p.7–8): photometric uncertainties + rotation trends → "we do not use the information obtained from the photometric residuals." We follow: V enters only through H(ρ), never as a fit residual |
| 23 | digest2 = same family, binned S3M look-up, our L1 target | K19 | K19 §2 (short-arc problem, admissible region), §2.1 (PANGLOSS heritage), §3.1 endpoint synthesis, §3.3 model-population look-up, Appendix A parabolic limit. Our **L1** (S3M in numerator & denominator) is a digest2 replication used to freeze geometry knobs |
| 24 | Truth p-value diagnostic | F15 | F15 Table 3 + §6.4–6.6: p-value of the true solution under the computed posterior — the tool that convicted Jeffreys. We run the same test using Kurlander's true Range/RangeRate columns (diagnostics only, never scoring) |
| 25 | 2-tracklet vs 1-tracklet information jump (context) | F15, SG18 | F15 Table 2 (TC3/AA: IP goes ~10⁻²–10⁻³ → ~1.0 adding the 2nd tracklet); SG18 Fig. 4. Context for why single-tracklet classification is fundamentally probabilistic — the regime digest2 and we both live in |

---

## 2. What we deliberately do differently (and the defence)

1. **No per-node least squares.** F15/S18 fit the attributable at every node because they have ≥3
   noisy observations. Our tracklets are 2-detection: A is exact, residuals are identically zero →
   the fit is not simplified, it is *absent*. Consequence: the whole engine is vectorised
   (tracklet × node) numpy; per-full-eval-set cost ≈ 1.4×10⁹ closed-form node evaluations.
2. **Coarser grid than S18.** Their 50×50→100×100 two-step densification chases χ² minima on the
   (ρ,ρ̇) plane. With χ²≡0 there are no minima; the integrand is smooth population density ×
   geometry. 64×32 AR-restricted nodes + documented convergence check (128×64, require |ΔF1|<0.002).
3. **Class ratio, not impact probability.** F15 §4's weighted sum with the impact flag replaced by
   population densities per class. This is also the fix for the old failed `NEO_H.py` attempt, which
   returned a raw weight instead of a ratio.
4. **NEOMOD3 in the numerator.** F15 §3.2 used S3M (Grav 2011) as f_pop and even anticipated the
   upgrade: "new population models are expected in the future (Granvik et al. 2014)." NEOMOD3
   (Nesvorný et al. 2024) is that model: debiased, H-dependent orbital distribution. digest2 still
   carries 2011 S3M (K19); this is precisely the wedge 1A tests.
5. **The L0/L1/L2 ladder with knob-freezing.** L0 (uniform-in-AR ≈ S18 score) → L1 (S3M everywhere ≈
   digest2 replication; all geometry knobs frozen by maximising L1↔`P_NEO_d2` agreement) → L2
   (NEOMOD3 numerator). L2−L1 is then the isolated NEOMOD3 effect with zero tuned parameters —
   the solo-work credibility discipline.

## 3. Open knobs (all parameterised in the engine)

| Knob | v1 value | Source / alternative |
|---|---|---|
| ρ range | [0.01, 100] au | S18 a_max=100; R_SI cut. Alt: extend down for shooting-star studies (blocked by H≤34.5 anyway) |
| N_ρ × N_ρ̇ | 64 × 32 (log-ρ × AR-restricted) | S18 Table 1 precedent; convergence check mandatory |
| Spatial prior exponent γ | 2 | F15 §3.2 (ρ²). Alternatives: 0 (F15 §3.3 uniform), 4 (full d³x d³v volume incl. velocity-box scaling). Frozen via L1 calibration |
| Population binning (denominator) | NEOMOD3-matched: dH=0.25, da≈0.10, de=0.04, di=4° | K19 §3.3 digest2 bins differ; smoothing σ=1 bin optional |
| Filter→V colors | phase-3 digest2 conversion table | Kurlander Table 3 colors as alternative |
| Phase function | HG, G=0.15 | Matches Kurlander inputs exactly (their §2.2.1) |
| Wagg ×0.80 MBA rescale | off (calibration-only) | Ranking/AUC-neutral; on for calibrated posteriors |
| ≥3-detection tracklets | not supported in v1 | Would need F15 Eq.(1) constrained fit per node |
