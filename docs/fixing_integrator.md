# Fixing the integrator: replacing two-body propagation with ASSIST n-body

**Created:** 2026-07-16. **Updated:** 2026-07-18 (Milestone B complete for case1; §13 runbook added).
**Status:** Integrator **BUILT + VALIDATED** (§9); regen **DONE for case1** — Stage 0 cache (§10.9),
production port (§11), cache-fed 667-grid n-body maps (§12.4.2), full re-score + ROC (§12.4.3).
**Result:** n-body maps discriminate **as well as** two-body (AUC identical) and **better when
epoch-matched** (night-61642: AUC +0.011, F1 +0.005); pipeline now self-consistent (n-body↔n-body).
**Now:** Arnor deliverables shipped (§13.2); **case2/case3 deferred** (grid is case-independent —
each needs only the ~9-min V-band re-score).
**Owners:** Devanshi + advisor (decisions), Hyak-Claude (build), Arnor-Claude (validation/analysis)
**Related:** `SORCHA_MAG_BAND_FINDING.md` (the *other*, independent bug — magnitude band),
`NOTE_FOR_HYAK_benchmark_objid_mapping.md`, `scoring_consistency_night61642.ipynb` (Arnor).
**Runbook + deliverables + caveats:** §13.

> **Headline:** ASSIST n-body + an origin-frame fix took `S0000vowa` from **71.95° → 0.005086°**;
> adding light-time + topocentric `X05` reaches **0.12 arcsec** (§9.9). A **second, independent
> production bug** was found en route (barycentric vs heliocentric observer — §9.2), live in
> `fast.py` / `gmm.py` / `hybrid.py`.
>
> **ADVISOR DECISION (2026-07-17): changing the integrator means regenerating ALL maps and ALL
> scoring.** A pipeline with n-body scoring against two-body maps is not defensible. §10 is the
> plan; it is **feasible and likely *cheaper* than the current two-body build** (§10.2).

---

## §1. The bug (confirmed at source)

The benchmark and Sorcha **disagree about where objects are on the sky**. Arnor found it via the
identity-matched night-61642 test; I traced it to the propagator.

| | Benchmark / VDP pipeline | Sorcha |
|---|---|---|
| Propagator | `neoscore._propagate_chunk_newton` — Newton's method on Kepler's equation | ASSIST + REBOUND |
| Physics | **pure two-body**, Sun's μ only, **zero planetary perturbations** | **full n-body** |
| Call site | `velocity_density_pipeline_*.py::build_visible_subset_dataframe` → `nsc.elements_to_helio_ecliptic_state(..., method="newton", n_iter=10)` | `ephemerides_type = ar` in the Sorcha config |
| Span | MJD 54466 (2008-01-01, S3M epoch `t_0`) → 61642 (2027-08-25) = **19.6 years** | same |

The two-body code advances mean anomaly linearly from `t_p` and solves `M = E − e·sin(E)`:
```python
E = M.copy()
for _ in range(n_iter):
    E -= (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
```
Nothing else touches it. Over ~20 years that is a fiction for anything that meets a planet.

### The evidence (Arnor, 3,205 identity-matched objects, night 61642)

| population | median sky sep | expected from motion in the epoch gap | excess |
|---|---|---|---|
| TNO | 0.040° | 0.003° | 0.036° |
| MBA | 1.401° | 0.032° | 1.369° |
| NEO | 2.271° | 0.075° | 2.196° |

- correlation(separation, hours since epoch) = **0.026** → time-of-night explains ~2%. **Not** an epoch
  offset; Arnor explicitly retracted that earlier hypothesis.
- The ordering *is* the diagnosis: divergence scales with how far the object moves along its orbit
  and how often it meets a planet. TNOs (centuries-long periods, no encounters) agree to 0.04°;
  MBAs (Jupiter perturbations) 1.4°; NEOs (planet-crossing, chaotic) 2.3°.
- Classic **along-track** signature: displaced *along its own orbit* → nearly the same apparent rate
  (which is why |Δv| median 0.0106 deg/day survives) but the **wrong place on the sky**.
- Direct link: |Δv| median is 0.0077 for pairs with sep < 0.5°, but 0.0570 for sep > 5° (**7× worse**).
  15% of pairs (481) have sep > 5°; ~half are NEOs.

### The worked example — `S0000vowa` (the acceptance test)

S3M source record (`S0.s3m`, the NEO file), decoded:

| element | value |
|---|---|
| q | 0.698174 AU |
| e | 0.692393 |
| i | 18.896° |
| Ω | 5.175° |
| ω | 236.446° |
| t_p | MJD 54100.18 |
| H | 23.974 (Johnson V) |
| **t_0 (epoch)** | **MJD 54466.0 = 2008-01-01** |

Derived a = q/(1−e) = **2.2697 AU** — matches Sorcha's `a_au = 2.2696940` exactly; e matches to 12
decimals; H matches the benchmark's `H_V` exactly. **The identity join is definitively correct — same
object, same source record, both sides.** It is a small (H≈24), high-e Apollo: q = 0.70 AU crosses
both Earth's and Venus's orbits, so ~20 years of two-body propagation is hopeless.

| field | Sorcha | benchmark |
|---|---|---|
| RA / Dec | 18.473° / −36.583° | 197.744° / −71.465° |
| ecl lon / lat | 359.723° / −40.392° | 236.990° / −56.239° |
| vlam / vbeta | −1.632 / +1.471 | +5.181 / +0.282 |
| \|v\| | 2.197 °/day | 5.188 °/day |
| prob map | `dlon+030_lat-35` | `dlon-090_lat-50` |
| **P_NEO_vdp** | **1.000** | **0.000** |
| P_NEO_d2 | 1.000 | 1.000 |

**sky separation = 71.95°.** The benchmark's own rate says it should have moved **2.03°** in the
9.39 h between epochs. It is **72°** away — a **35× discrepancy**. Different prob map, opposite score.

`P_NEO_d2 = 1.0 on both sides` is the control: digest2 works from the *orbit*, not the sky position,
so the ephemeris divergence doesn't touch it. Only VDP — which is a function of sky position — breaks.

### What holds / what needs a caveat
- **Holds:** the kinematic validation (where positions agree, |Δv| = 0.0077 deg/day); the magnitude
  /band findings (independent); "96% of objects score identically."
- **Needs a caveat:** the ~4% that differ are **not** boundary noise — they are objects the two
  pipelines place in genuinely different sky locations. Arnor's Task-3 buckets "velocity tail" and
  "map mismatch" are both **downstream of this one root cause**.

---

## §2. Constraint discovered: there is no working n-body path in this codebase

`neoscore.py` offers `method="adam"` as an alternative to `"newton"` — but it routes to
`neomod/adam_core_stub/`, which is **10 lines across 4 files**: empty class shells for
`KeplerianCoordinates`, `Timestamp`, `Origin`. A mock, not a library. **Flipping the method flag is
not a fix.** The stub is injected by `sys.path.insert(0, ".../adam_core_stub")` in the pipeline
scripts — which also means it *shadows* the real `adam_core`.

**But the real ecosystem is already here** (verified on Hyak):

| package | status |
|---|---|
| `adam_core` | **real, 0.5.5, installed** in `conda_prep` (shadowed by the stub) |
| `assist` | 1.2.3 installed (Sorcha's own n-body engine) |
| `rebound` | 4.6.0 installed |
| `adam_assist` (`ASSISTPropagator`) | **NOT installed** → this plan installs it |
| ASSIST ephemeris data | Sorcha cache has `linux_p1550p2650.440` (planets) + `de440s.bsp`; `adam-assist` bundles `jpl_small_bodies_de441_n16` (asteroid perturbers) |

---

## §3. The fix: replace **one call**, not the whole function

`build_visible_subset_dataframe` (`fast.py:475–627`) has two separable layers:

1. **Lines 504–510 — propagation.** `nsc.elements_to_helio_ecliptic_state(...)` : elements →
   heliocentric **ecliptic** state (r km, v km/s) at obstime. **This is the only two-body physics.
   This is the entire bug.**
2. **Lines 512–627 — geometry.** Ecliptic→equatorial rotation, observer subtraction
   (`scorer._get_earth_and_observer`), analytic RA/Dec rates, sky cut, and the manual
   equatorial→ecliptic rotation producing `vlam`/`vbeta`. This math is **validated against Horizons**
   and is protected by the `PROJECT_CONTEXT.md` warning (astropy's `GeocentricTrueEcliptic` transform
   silently drops differentials → wrong rates). **Do not touch it.**

So: **swap step 1 for an n-body state provider with the identical output contract**
(`r_helio_ecl_km`, `v_helio_ecl_km_s`, same row order), leave step 2 byte-identical. That also makes
the experiment clean: **exactly one variable changes** between the two pipelines.

### Deliberate choice: do NOT use `generate_ephemeris()` as the drop-in
The advisor's snippet uses `propagator.generate_ephemeris(orbits, observers)`, which does its *own*
observer subtraction, light-time correction, and frame conversion, returning RA/Dec/rates directly.
Using it would **replace our validated geometry layer** and confound the comparison (we'd change
propagation *and* geometry *and* light-time at once, and couldn't attribute any residual).
→ Use `propagate_orbits()` (orbit propagation only) and feed our existing math.
→ Keep `generate_ephemeris()` as an **independent cross-check**, not as the pipeline.

### New requirement: `t_0` becomes load-bearing
The two-body path never reads `t_0` (the osculating epoch, MJD 54466) — only `t_p`. For n-body it is
**essential**: elements→state **at t_0** is exact by definition (that is what "osculating" means —
zero propagation error at the epoch), then ASSIST integrates t_0 → obstime with full physics.
`df` already carries `t_0` from `s3m_loader` (`DEFAULT_COLS`), so it only needs threading through.
*Using `t_p` as the epoch instead of `t_0` would be a subtle, silent, catastrophic error.*

---

## §4. Build plan — `velocity_density_pipeline_adam.py`

A copy of `fast.py` with the propagation swapped. New file, so `fast.py`/`gmm.py` stay untouched and
every existing result stays reproducible.

**Step 1 — install** `adam-assist` into `conda_prep`.
- Dry-run verified: **pure additions, zero upgrades** — `adam-assist-0.3.10`, `cffi`,
  `cryptography`, `jpl_small_bodies_de441_n16`, `timezonefinder`. `adam_core` stays 0.5.5;
  `assist`/`rebound`/`numpy`/`astropy`/`pandas`/`sorcha` untouched.
- **Gate:** after install, `import sorcha` and `import velocity_density_pipeline_fast` must both
  still work. If either breaks → uninstall, fall back to §6(b).

**Step 2 — `propagate_elements_nbody()`** (the new function, the only real new code):
```
S3M elements (a,e,i,Ω,ω,t_p) @ t_0
  -> KeplerianCoordinates (adam_core)          # exact at t_0, no propagation yet
  -> Orbits(coordinates=..., time=Timestamp(t_0, scale="tdb"))
  -> ASSISTPropagator().propagate_orbits(orbits, Timestamp(obstime))   # n-body t_0 -> obstime
  -> .coordinates.to_cartesian()               # heliocentric
  -> return (r_km, v_km_s) in the SAME frame/units the newton path returned
```
Critical bookkeeping to verify, not assume:
- **Frame:** adam_core cartesian is typically **ecliptic** with `Origin("SUN")`. The newton path
  returns heliocentric **ecliptic** (`fast.py:512-519` then rotates ecliptic→equatorial). Must
  confirm adam's frame/origin and match it exactly — a silent ecliptic/equatorial mixup would look
  like a ~23.4° error.
- **Units:** adam_core uses **AU / AU-per-day**; the newton path returns **km / km-per-s**. Convert.
- **Time scale:** `t_obs = Time(obstime_str, scale="tdb")` in the existing code → keep TDB.
- **Row order:** must be preserved 1:1 (downstream indexes by position).

**Step 3 — swap the call** in the copied `build_visible_subset_dataframe`; everything below line 512
unchanged.

**Step 4 — cost check.** ASSIST n-body is *far* more expensive than closed-form Kepler. The
benchmark-scale question (8,436 objects fine; ~460k v3 objects; millions of map clones) is measured
in §5, not assumed.

---

## §5. Acceptance test (pre-registered, before any results)

**Primary — `S0000vowa`, the 72° object.** Propagate its S3M record from t_0 = 54466 to the Sorcha
detection time with ASSIST, run it through the unchanged geometry layer, and compare to **Sorcha's
own reported position** (RA 18.473°, Dec −36.583°) — Sorcha is ground truth here, since it *is* the
ASSIST n-body ephemeris we're trying to match.

| | criterion |
|---|---|
| **PASS** | sky separation **< 0.01°** (arcsec-level; residual = light-time/observer-frame detail) |
| **Partial** | < 0.1° — right physics, some frame/timing detail off; investigate |
| **FAIL** | ≳ 1° — still diverging; force model or epoch handling wrong |
| Reference | current two-body = **71.95°** |

Compare at Sorcha's **actual visit time** (mjd0 = 61642.391188), not the benchmark's 00:00 epoch, so
the epoch gap isn't confounded.

**Secondary:**
1. **Population sweep** on the 3,205 matched objects: median sky sep must collapse for **all three**
   populations (target ≪ 0.04° TNO / 1.40° MBA / 2.27° NEO). NEO improving most is the signature.
2. **Non-regression:** TNOs must *not* get worse (they already agree at 0.04°) — catches a frame/unit
   bug that a NEO-only test would miss.
3. **|Δv| link:** the sep>5° bucket (|Δv| 0.0570) should fall toward the sep<0.5° bucket (0.0077).
4. **Scoring:** S0000vowa should land in the *same prob map* as Sorcha and score P_NEO_vdp ≈ 1.0
   (vs 0.000 today).
5. **Independent cross-check:** `generate_ephemeris()` on the same object should agree with our
   geometry-layer RA/Dec — validates that reusing our own geometry was legitimate.

---

## §6. Risks & fallbacks

| risk | mitigation |
|---|---|
| pip breaks `conda_prep` (every job depends on it) | dry-run done (additions only); post-install import gate; fallback (b) |
| `adam_core_stub` shadows real `adam_core` | the adam pipeline must **not** put the stub on `sys.path` — verify `adam_core.__file__` points at `site-packages`, not the stub |
| frame/unit/origin mismatch (ecliptic vs equatorial, AU vs km) | explicit assertions + the TNO non-regression test |
| wrong epoch (`t_p` instead of `t_0`) | called out in §3; assert `t_0` is finite and ≈54466 |
| n-body too slow at map scale | measure first (§4 step 4); maps likely stay two-body (§7) |
| ASSIST ephemeris data missing/downloading | `jpl_small_bodies_de441_n16` comes with the install; Sorcha's cache has the planetary kernels |

**Fallback (b):** drive `assist`+`rebound` directly using **Sorcha's own cache files**
(`linux_p1550p2650.440`, `de440s.bsp`). Gives *exact* force-model and data parity with Sorcha —
literally the stated goal — and needs **zero new installs**. More code (barycentric↔heliocentric and
ecliptic↔equatorial bookkeeping around the integrator). Adopt if the install threatens the env.

---

## §7. Scope

**In scope:** the **benchmark / scoring** path, where per-object sky position must match Sorcha.

**Out of scope (documented approximation):** the **maps** (`prob_maps_grid_s3m`) are built with the
same two-body propagator. Two-body conserves (a, e, i), and a map is a *statistical ensemble* — the
velocity distribution at a given geometry is set by the orbit distribution, not by which individual
object sits there — so the aggregate is approximately preserved. n-body over 19.6 yr *would*
redistribute NEO elements somewhat (a real second-order effect), but n-body integrating millions of
clones is likely infeasible. **Decide explicitly with the advisor; record as a known limitation
either way.** This fix targets where it demonstrably matters: per-object positions.

**Also unaffected:** digest2 (`P_NEO_d2` = 1.0 on both sides for S0000vowa) — it reads the orbit, not
the sky position.

---

## §8. Execution order

1. `pip install adam-assist` → import gate (`sorcha`, `velocity_density_pipeline_fast`).
2. Confirm real `adam_core` resolves (not the stub) + probe the frame/origin/unit contract of
   `propagate_orbits` on one object.
3. Write `velocity_density_pipeline_adam.py` (copy of `fast.py` + `propagate_elements_nbody`).
4. **Acceptance test on S0000vowa** (§5 primary). Stop and report — pass or fail.
5. If PASS → population sweep (§5 secondary) on the 3,205 matched objects.
6. If PASS → cost measurement, then rebuild `benchmark_night61642` with the adam pipeline and hand
   to Arnor for the full same-object scoring test.
7. Update this doc with results; commit.

---

# §9. RESULTS — built, tested, verified (2026-07-16)

Implementation lives in **`codex_integrator/`** (deliberately outside `neomod/src/`, so the
production VDP modules are untouched and every existing result stays reproducible):
- `velocity_density_pipeline_adam.py` — copy of `fast.py` with the propagation swapped.
- `acceptance_s0000vowa.py` — the pre-registered §5 test.
- `README.md` — usage.

## §9.1 ACCEPTANCE TEST: **PASS** (independently reproduced)

```
object: S0000vowa
visit MJD UTC: 61642.391188036
ASSIST RA/Dec: 18.493305543, -36.583784014 deg
Sorcha RA/Dec: 18.492736355, -36.588849838 deg
sky separation: 0.005086399 deg
ASSIST vlam/vbeta: -1.595362868, 1.440684194 deg/day
Sorcha vlam/vbeta: -1.631790951, 1.471257047 deg/day
acceptance (<0.01 deg): PASS
```

| | value |
|---|---|
| two-body (before) | **71.95°** |
| n-body alone | 0.766° |
| n-body + origin fix | **0.005086°** ✅ |
| pre-registered PASS bar (§5) | < 0.01° |
| improvement | **~14,000×** |

Run it with:
```bash
MPLCONFIGDIR=/tmp/mpl PYTHONPATH=codex_integrator:neomod/src \
  conda_prep/bin/python codex_integrator/acceptance_s0000vowa.py
```
The test correctly compares at Sorcha's **actual visit time** (`mjd0_utc`), not the benchmark's
00:00 epoch, so the epoch gap is not confounded (§5 requirement).

### Implementation notes — it took the *better* route
It did **not** use `adam-assist`'s `ASSISTPropagator` (the §4 plan); it used **fallback (b)** (§6):
driving `assist`+`rebound` directly with **Sorcha's own** `universal_cartesian()`, `gm_sun` from the
ASSIST ephemeris, the production planetary kernel, `GR_SIMPLE` forces and IAS15. That is *exact
force-model parity* with Sorcha — literally the stated goal — rather than "a different n-body code
that ought to agree." `adam-assist` was still worth installing: it supplies
`jpl_small_bodies_de441_n16` (the asteroid-perturber kernel).

Verified against the §3/§5 requirements:
- ✅ **`t_0` used as the epoch** (not `t_p`) — the load-bearing requirement in §3.
- ✅ Frame chain correct: ecliptic → equatorial → +Sun (barycentric) → integrate → −Sun →
  equatorial → ecliptic; km and km/s out, matching the old `neoscore` contract.
- ✅ Row order preserved; inputs validated (finite, `a>0`, `0≤e<1`).
- ✅ Geometry layer (the Horizons-validated math) untouched.
- ✅ Real `adam_core` resolves, not the 10-line stub.

## §9.2 NEW FINDING — a second, independent bug **in production**

The acceptance test surfaced a bug that has nothing to do with two-body vs n-body:
**`_get_earth_and_observer()` returns a *barycentric* Earth state, but the VDP pipelines subtract it
from a *heliocentric* asteroid.**

Evidence, at source:
- `neoscore.py:183` docstring: *"Get Earth's **barycentric** state and observer topocentric offset."*
- `neoscore.py:190-193`: `get_body_barycentric_posvel('earth', obstime)` → `rE` is returned
  **with no Sun subtraction**. It is barycentric, full stop.
- `velocity_density_pipeline_fast.py:527` does `r_rel = r_obj_helio - (rE + r_tele)` — a
  **heliocentric** asteroid minus a **barycentric** observer. Same pattern in `gmm.py`, `hybrid.py`.

> **CORRECTION (2026-07-16).** An earlier revision of this section claimed *"neoscore's own callers
> convert, so it isn't a problem there."* **That was wrong.** `neoscore._geo_to_helio()` (line 208)
> computes `r_helio = rE + r_tele + r_geo_km` with a **barycentric** `rE` — so it returns a
> **barycentric** vector while naming itself heliocentric. The `# Convert to heliocentric` comment
> at line 101 is aspirational, not implemented. **neoscore carries the same defect.**

### Why it is fatal in `fast.py` but invisible in `neoscore.py`
Not because neoscore converts (it doesn't) — because of **what each does with the vector**:

| | uses the mixed vector to form | denominator | resulting error |
|---|---|---|---|
| `neoscore` | heliocentric **distance** (`r_sun`, for magnitude) | ~1–3 AU | 0.01/2 → **~0.01 mag** — invisible |
| `fast.py:527` | observer→object **direction** (RA/Dec) | Δ ≈ **0.3 AU** for an NEO | 0.01/0.3 → **~1–2°** — catastrophic |

Identical defect, identical ~0.005–0.01 AU offset. Harmless divided by 2 AU; fatal divided by 0.3 AU.
That is why it survived undetected, and why it is worst exactly where we care most — nearby NEOs.

The fix (in the adam copy) converts the observer to heliocentric first:
```python
r_sun_bary, v_sun_bary = get_body_barycentric_posvel("sun", t_obs)   # de432s
r_observer_helio = rE - r_sun_bary + r_tele
v_observer_helio = vE - v_sun_bary
r_rel = r_obj_helio - r_observer_helio
```

**Scale of the error:** ≈ (Sun–barycenter offset)/Δ. The offset is ~0.005–0.01 AU (Jupiter-dominated),
so the angular error is **worst for nearby objects and negligible for distant ones** — ~1° for an NEO
at Δ≈0.3 AU, ~arcsec for a TNO. It was worth **0.766°** on `S0000vowa`.

**Blast radius:** contaminates the **benchmark** and the **maps** (both route through
`build_visible_subset_dataframe`). Sorcha-side tracklets are **safe** — their RA/Dec come from Sorcha
itself, and the VDP scoring path only does the velocity-space lookup. This is a live defect in
`fast.py` / `gmm.py` / `hybrid.py` and needs its own decision (see §9.4).

## §9.3 GAPS — must close before rebuilding the benchmark

1. **✅ CLOSED (2026-07-16) — see §9.10.** ~~🔴 The magnitude path is still two-body *and* still has the origin bug.~~
   `compute_apparent_magnitude_for_population` (adam copy, line 526) still calls
   `nsc.elements_to_helio_ecliptic_state(...)` (the newton two-body propagator), and still uses the
   uncorrected `r_obs_vec = r_eq - (rE + r_tele)`. `score_orbital_df` (line 2210) feeds it into
   `__mag_app_full` → **the mag bin**.
   → The module currently yields **n-body positions with two-body magnitudes**. For `S0000vowa` the
   two-body geometry is 72° off, so its `r_sun`/`Δ` are garbage → wrong apparent mag → **wrong mag
   bin → still-wrong VDP score despite the correct position.** Because VDP scores on
   `(vlam, vbeta, mag_bin)`, **this must be fixed before the benchmark rebuild.** (The README notes
   legacy paths stay two-body — but this one sits on the critical path.)
2. **🟡 Light-time is not applied** (zero occurrences in the module) — but **Sorcha applies it**
   (its output carries `Range_LTC_km`, `RangeRate_LTC_km_s`, `Obj_Sun_*_LTC_km`).
   At Δ = 0.296 AU light-time is 148 s; at 2.2 °/day that is **0.0038°** — which accounts for
   essentially all of the residual 0.005086°. Applying it should tighten to ~0.001°.
   Optional (we already PASS), but it explains the residual and is needed for true parity.

## §9.4 Decisions needed

| # | decision | notes |
|---|---|---|
| 1 | Fix the magnitude path (§9.3.1) before rebuilding? | **Recommend yes** — otherwise the benchmark is subtly wrong in the mag bin, which is exactly what VDP scores on. |
| 2 | Apply light-time (§9.3.2)? | Recommend yes for parity; explains the residual. Not required to pass. |
| 3 | **Fix the origin bug in production** (`fast.py`/`gmm.py`/`hybrid.py`)? | It is a live defect affecting the benchmark and the maps. Fixing changes existing numbers — needs the advisor. |
| 4 | Do the **maps** get rebuilt (n-body and/or origin fix)? | §7 said maps stay two-body as a documented approximation. The origin fix is *cheap* (no propagation change) and may be worth doing even if n-body isn't. |

## §9.5 Cost — not a blocker

Measured: **~0.05 s/object** steady-state (the 2.70 s on the first call is one-time ASSIST ephemeris
load; n=5 → 0.041 s/obj, n=20 → 0.057 s/obj).

| target | objects | single-thread |
|---|---|---|
| `benchmark_night61642` | 8,436 | **~7 min** |
| benchmark v3 | 460,083 | ~6.4 h (trivially arrayed → minutes) |

So §4-step-4's cost worry is resolved: the benchmark rebuild is cheap. Available optimization if ever
needed: every S3M object shares `t_0 = 54466`, so all particles could go into **one** simulation
instead of one integration each (current code builds a `rebound.Simulation` per object).

## §9.7 REFRAME (2026-07-16, user) — the goal is **correctness**, not Sorcha-parity

An earlier revision argued *"the goal is match Sorcha, not match JPL, so parity beats convention."*
**That was a false dichotomy and is superseded.** The goal of an integrator is to be **right**;
Sorcha-agreement is then a **result we report**, not a target we tune to.

The decisive fact: **ASSIST *is* the JPL-grade integrator** (integrates against JPL DE440/441,
designed to reproduce Horizons). Its default force model is `GR_EIH` (full Einstein–Infeld–Hoffmann).
Sorcha then *deliberately downgrades it*:
```python
forces.remove("GR_EIH"); forces.append("GR_SIMPLE")   # sorcha/ephemeris/simulation_setup.py:182-183
```
So **"match Sorcha" means inheriting Sorcha's approximation.** `ASSISTPropagator` defaults are
strictly closer to truth. If our (validated) side disagrees with Sorcha, that is a finding *about
Sorcha* — which we can only say if we are independently validated.

### Honest caveat: for the *force model*, this is numerically moot
Order-of-magnitude for `S0000vowa` (a=2.27, e=0.69, ~5.7 orbits in 19.6 yr): total GR perihelion
precession ≈ **0.2″**; `GR_SIMPLE` already captures the dominant solar term, so
`GR_SIMPLE` vs `GR_EIH` ≈ **milliarcsec** — ~4 orders of magnitude below our current 0.005° = 18″
residual. **To be measured, not trusted (§9.8).** What actually matters at our precision:

| term | size | status |
|---|---|---|
| **light-time** | ~18″ | ❌ omitted — **our entire residual** |
| topocentric parallax (Δ=0.3 AU) | ~29″ | ✅ handled (`r_tele`) |
| `GR_SIMPLE` vs `GR_EIH` | ~0.001–0.01″ | irrelevant either way |

So the win from `adam_core` is **light-time + magnitudes**, not the force model.

### The validation standard this unlocks
S3M objects are **synthetic** — `S0000vowa` has no Horizons entry, so it cannot be validated against
JPL directly. Correct structure instead:
1. Validate the **code** against Horizons on a **real** object (Apophis, `query_sbdb(["99942"])` +
   `query_horizons_ephemeris`) — external, authoritative truth.
2. Apply the validated code to the synthetic S3M objects.
3. *Then* compare to Sorcha — residuals are now attributable, because our side is independently right.

This is strictly better than "match Sorcha to <0.01°", which only ever proved the two **agree**,
never that either is **correct**.

### Resulting two-path design (adopted)
| path | purpose | force model |
|---|---|---|
| **direct-assist** (`propagate_elements_nbody`, current) | **Sorcha-parity reference** — "do we reproduce Sorcha's own choices?" | `GR_SIMPLE` + Sorcha's `universal_cartesian`/kernels |
| **adam_core / `ASSISTPropagator`** (to build) | **correctness path** — `GR_EIH` default, light-time, `Observers.from_code("X05")`, `predict_magnitudes=True` | ASSIST default |

## §9.8 DECISION — propagate **once**, derive magnitude and geometry from the same state

`score_orbital_df` currently propagates **twice**: once via
`compute_apparent_magnitude_for_population` (for `__mag_app_full` → the mag bin) and once via
`build_visible_subset_dataframe` (for position/rates). With two-body Kepler that was ~free; with
n-body at ~0.05 s/object it **doubles the cost**, and — worse — it lets the two propagations drift
apart, which is exactly the §9.3.1 bug (n-body positions, two-body magnitudes).

**Adopted:** propagate once, derive both from that single state. Cheaper *and*
**consistency-by-construction** — the magnitude and the sky position can no longer disagree about
where the object is. This supersedes the "minimal fix" of merely repointing
`compute_apparent_magnitude_for_population` at the n-body propagator.

## §9.9 CORRECTNESS-PATH RESULTS (2026-07-16) — `codex_integrator/adam_correctness_checks.py`

Built the `adam_core` + `ASSISTPropagator` correctness path per §9.7 and ran all three checks.

### A. `generate_ephemeris` vs Sorcha — **0.12 arcsec**

| configuration | separation from Sorcha |
|---|---|
| two-body (the original bug) | 71.95° = 259,020″ |
| direct-assist n-body, **no light-time** (§9.1) | 0.005086° = **18.31″** |
| **adam: n-body + light-time + `X05` topocentric + `GR_EIH`** | **0.0000335° = 0.1206″** |

`light_time = 147.535 s` — matching the ~148 s predicted in §9.7 for Δ = 0.296 AU. **Light-time was
the entire 18″ residual**, exactly as diagnosed. Applying it improves agreement **~150×**.

### Magnitude — adam predicts **Johnson V natively**
`Orbits.physical_parameters` is `PhysicalParameters(H_v, G, ...)` and the output column is
**`predicted_magnitude_v`** — i.e. adam's photometric model is **Johnson V**, *the exact band our
maps use*. Supplying `H_v` = the S3M `H` (which is Johnson V) and `G = 0.15` yields:

| quantity | value |
|---|---|
| S3M `H_v` (input) | 23.9740 |
| **adam `predicted_magnitude_v`** (at the correct n-body position) | **23.2895** |
| Sorcha `mean_mag` (r/i, LSST band) | 23.1406 |
| Sorcha `mean_mag_V` (our colour conversion, §SORCHA_MAG_BAND) | 23.7426 |
| benchmark `mean_mag` (two-body synthetic V) | 24.0085 |

> Do **not** over-read the 0.45 mag gap between adam's 23.29 and our `mean_mag_V` 23.74 for *this*
> object: it has `snr0/snr1 ≈ 4.2`, so Sorcha's observed PSFMag carries ~0.25 mag of photometric
> noise, plus trailing-loss. A single low-SNR object cannot adjudicate this — it needs the
> population sweep. The benchmark's 24.01 is meaningless regardless (two-body put the object 72°
> away, so its `r_sun`/Δ are garbage).

**Consequence for §9.3.1:** `predicted_magnitude_v` gives a clean synthetic **V** magnitude at the
**correct n-body position, from the same propagation** — solving the two-body-magnitude bug *and*
the band problem in one call, and satisfying §9.8 (propagate once) by construction.

### B. `GR_SIMPLE` vs `GR_EIH` — **measured: 0.0039 arcsec** (§9.7 prediction confirmed)

`ASSISTPropagator.__init__` exposes only `min_dt, initial_dt, adaptive_mode, epsilon` — **no force
selection** — so the A/B was done in the direct-assist path, which does expose forces (a `gr_model`
parameter was added to `propagate_elements_nbody`; default `"GR_SIMPLE"` preserves Sorcha-parity).

| | value |
|---|---|
| \|Δr\| between force models (after the full 19.6 yr) | **0.833 km** |
| angular impact at Δ = 0.296 AU | **0.00388″** |
| §9.7 prediction | ~0.001–0.01″ ✅ |

So the force-model choice is **~4,700× smaller than the light-time term** and ~30× below our final
0.12″ residual. **The §9.7 principle stands (aim for correctness) and costs nothing either way.**
*Caveat:* 0.833 km is this object's *post*-amplification result; an S3M NEO with a deeper encounter
could amplify the force difference more (see C). A population-level check would settle it.

### C. Apophis vs JPL Horizons — **PASS**, and an instructive lesson in chaos

S3M objects are synthetic (no Horizons entry), so the **code** was validated on a real object.
`query_horizons_ephemeris` does not exist in adam_core 0.5.5 — the correct API is **`query_horizons`**,
which returns **state vectors**, so this compares propagated *position* directly (no geometry confound).

| from SBDB epoch (MJD 59215) | vs the 2029 close approach (MJD ≈ 62239, +8.28 yr) | \|Δr\| |
|---|---|---|
| +1.0 yr | before | **1.2 km** |
| +5.0 yr | before | 14.6 km |
| +8.0 yr | before | **38.2 km** |
| +8.5 yr | **after** | **22,870 km** |
| +10.0 yr | after | 363,691 km |

**Verdict: PASS.** Absent encounters our propagation tracks JPL to **~5 km/yr** (38 km over 8 years) —
the integrator is correct. The **600× jump across a single encounter** (38 km → 22,870 km in half a
year) is **physics, not a bug**: a deep planetary encounter exponentially amplifies the tiny
difference between SBDB's orbit solution and Horizons' internal one. It is a live demonstration of
exactly the mechanism that makes the S3M NEO problem hard (§1).

### The nuance this forces — and why Sorcha-parity is still the right standard *for synthetic objects*
§9.7 argued "aim for correctness, not Sorcha-parity." C sharpens that:
- For validating the **code**, JPL is the standard → **PASS** (~5 km/yr, pre-encounter).
- For **synthetic** S3M objects there is **no external truth** — "truth" *is* the propagation of the
  given elements. Two good integrators started from *identical* elements track each other
  (we got 0.12″ vs Sorcha); the Apophis divergence arises only because SBDB and Horizons start from
  *slightly different* solutions and an encounter amplifies it.
- Therefore: **validate the code against JPL, then match Sorcha's setup for the benchmark.** Both
  framings are right, for different jobs — and since GR_SIMPLE vs GR_EIH is 0.004″, the choice is
  moot in practice.

### Bugs found in my own checks (recorded so they aren't re-hit)
1. **`tp` is MJD, not JD**, in `CometaryCoordinates` — same scale as `time`
   (`transform.py:1434` does `dtp = tp - t0`; verified against `query_sbdb`: Apophis `tp=59100.539`
   vs epoch `59215.0`). Passing JD → `dtp ≈ 2.4e6 d` → garbage orbit → *"Light travel time is NaN"*.
   Note the trap: the **direct-assist** path *does* want JD, because Sorcha's `universal_cartesian()`
   uses the opposite convention.
2. **`Orbits` requires CARTESIAN coordinates.** Passing `CometaryCoordinates` silently yields a null
   state (`x/y/z = None`) → the same NaN-light-time error. Use `coords.to_cartesian()`.
3. `ASSISTPropagator` exposes **no force API** → do force A/Bs in the direct-assist path.
4. `query_horizons_ephemeris` **does not exist** in 0.5.5 → use `query_horizons` (returns states).

## §9.10 BLOCKER CLOSED — magnitude path fixed via the single-propagation restructure (§9.8)

Implemented in `codex_integrator/velocity_density_pipeline_adam.py`:
- `build_visible_subset_dataframe()` now derives `mag_app` from **the same n-body state and the same
  corrected (heliocentric) observer** it uses for the sky position — same HG(G=0.15) formula as the
  legacy function, just fed the right state. Magnitude depends only on `|r|` and the phase angle
  (both frame-invariant), so the existing equatorial vectors are used directly.
- `score_orbital_df()` now **consumes** that column instead of calling
  `compute_apparent_magnitude_for_population()` and propagating a second time. It raises if
  `mag_app` is absent, so the old double-propagation path cannot silently return.

### Result — independently cross-validated

| quantity | value |
|---|---|
| **NEW** single-propagation n-body `mag_app` | **23.2896** |
| adam_core `predicted_magnitude_v` (**independent code path**) | **23.2895** |
| OLD two-body benchmark `mean_mag` | 24.0085 |

Our hand-rolled HG magnitude on our direct-assist state agrees with adam_core's own photometry on
adam_core's own ASSIST propagation to **1×10⁻⁴ mag** — two entirely separate implementations landing
on the same number. And the correction is **0.72 mag** off the two-body value: **more than half a
1-mag VDP bin**, so the bug was real and would have mis-binned this object even with a perfect position.

The acceptance test is **unchanged at 0.005086399° (PASS)** — as expected, since the geometry layer
was not touched. Benefits: n-body magnitude ✅, origin fix applied to the magnitude ✅, **half the
n-body cost** (one propagation, not two) ✅, and position/magnitude consistency **by construction** ✅.

## §9.6 Next steps (superseded by §10)

1. Close §9.3.1 (magnitude path → n-body + origin fix). **Blocker.**
2. Optionally close §9.3.2 (light-time) → expect ~0.001°.
3. Re-run the acceptance test; then the §5 secondary sweep on the 3,205 matched objects —
   **including the TNO non-regression check** (TNOs already agree at 0.04°; a frame/unit bug would
   show up there and nowhere else).
4. Rebuild `benchmark_night61642` with the adam pipeline (~7 min) → hand to Arnor for the full
   same-object scoring test (expect `S0000vowa` to land in Sorcha's prob map and score ≈ 1.0, vs
   0.000 today).
5. Take §9.4 decisions to the advisor — especially the production origin bug.

---

# §10. THE FULL REGEN PLAN (advisor decision, 2026-07-17)

**Decision:** changing the integrator means **regenerating all maps and all scoring**. n-body
scoring against two-body maps is not a defensible pipeline.

This section supersedes §7 and §9.6. **§7 said "maps stay two-body as a documented approximation" —
that is now WITHDRAWN**, along with the argument behind it (see §10.1).

## §10.1 Two errors this plan corrects

1. **"Two-body conserves (a,e,i), so the map ensemble is fine."** Backwards. Two-body conserves them
   *by construction*; **n-body does not** — planetary perturbations evolve NEO elements over 19.6 yr,
   and **that evolution is precisely the physics the map should encode**. A two-body map represents
   the **2008** element distribution; the map is supposed to represent the population **at the
   observation epoch**.
2. **"n-body map generation is infeasible (~10⁴ node-hours)."** Wrong — it rested on the assumption
   that the ~37M *clones* must be n-body propagated. **They must not.** See §10.2.

## §10.2 The unlock: clones live at the observation epoch, so they need no n-body

`_gmm_feature_matrix()` builds the cloner's feature space using **mean anomaly at obstime**:
```python
M = _gmm_mean_anomaly_rad(df, obstime_str)
return np.column_stack([log(a), q, sin/cos(inc), sin/cos(node), sin/cos(argp), sin(M), cos(M)])
```
and `_gmm_unpack_features()` inverts it as:
```python
t_p = t_obs.mjd - M / n_rad_day     # t_p DERIVED from M at obstime
```

So a clone is **osculating elements defined at the observation epoch**. Propagating it from its `t_p`
to `t_obs` is a **pure Kepler solve with zero elapsed physics** — `t_p` is a bookkeeping device, not a
real epoch 19.6 years in the past. **Two-body Kepler for clones is not an approximation; it is exact
by definition.** No perturbation can accumulate over a time gap that does not exist.

**The n-body cost therefore applies only to the ~14.4M *real* S3M objects**, whose elements are given
at `t_0 = MJD 54466` and must cross 19.6 years to the map epoch.

| population | real objects (the only n-body cost) | clone factor (clones need **no** n-body) |
|---|---:|---:|
| NEO | 268,511 | ×80 |
| MBA | 13,883,361 | ×1 |
| TNO | 48,682 | ×10 |
| Trojan | 179,882 | ×5 |
| **total** | **14,380,436** | — |

**Cost:** 14.4M × ~0.05 s/object ≈ **200 core-hours** → **~24 min/task at 500-way array**. Feasible.

**And it should be *faster* than the current two-body build.** Today `sorcha_gen_maps_grid_s3m.sh` is
`--array=0-666%48, 16 cpu, 32G, 6h`, and each task re-propagates the population **3× per mag bin ×
8 mag bins** — i.e. the same object is propagated ~5,336× across the grid, because every center shares
the same obstime. Stage 0 eliminates that redundancy entirely.

## §10.3 Architecture — propagate once, globally, cache (the §9.8 principle at map scale)

```
STAGE 0  (NEW, one-time, the ONLY n-body cost)
   14.4M real S3M objects, elements @ t_0 = MJD 54466 (2008-01-01)
        │  ASSIST n-body  (19.6 yr)   ← the fix
        ▼
   state @ map epoch  ──► osculating elements @ epoch
                      ──► ra/dec, vlam/vbeta   (corrected heliocentric observer, §9.2)
                      ──► mag_app              (same state, §9.8 — free)
        │
        ▼   epoch_state cache (parquet, ~14.4M rows)

STAGE 1  (per center × magbin — 667 tasks — NO propagation)
   read cache → filter (mag bin + sky patch) → GMM fit on EPOCH elements
              → sample clones (elements @ epoch)
              → elements→state via two-body Kepler  ← EXACT (zero gap)
              → kNN density → maps
```

**Why this is more correct than today even setting the integrator aside:** the GMM currently fits
**2008** `a/e/i` with a two-body-advanced `M`. Stage 0 makes it fit the **true epoch element
distribution** — what the map is meant to represent.

## §10.4 Scope — everything downstream is invalidated and must be rebuilt

The origin fix (§9.2) alone changes every number; the integrator change compounds it.

| artifact | action |
|---|---|
| `prob_maps_grid_s3m/` (667 maps, 20 GB) | **regenerate** (Stage 0 + 1) |
| `prob_maps_grid/` (geometry/assignment maps) | regenerate **or** confirm geometry-only ⇒ unaffected |
| benchmark v1 / v2 / v3, `benchmark_night61642` | **rebuild** (n-body + origin + single-prop mag) |
| Sorcha VDP scoring: case1/2/3, v5 (`P_NEO_vdp`) | **re-score** against the new maps |
| `P_NEO_vdp_Vband` (mag-band fix) | **re-score** — must ride on top of the new maps |
| digest2 (`P_NEO_d2`) | **unaffected** — orbit-based, not sky-position-based (proved: 1.0 on both sides for `S0000vowa`) |
| Sorcha tracklets themselves | **unaffected** — RA/Dec come from Sorcha's own ASSIST |
| ROC/F1 in `NEOMplan.md`, case1/2/3 findings | **recompute** after re-scoring |

**Open for advisor:** the map **reference epoch**. Current default is `2026-01-01`
(`DEFAULT_REF_OBSTIME`), but the maps are antisun-relative and epoch-independent *by design* — while
Stage 0's n-body propagation is **not** epoch-independent (19.6 yr vs 21.2 yr to MJD 61642 differ).
If we ever want night-61642 work, Stage 0 must run at that epoch. **Cheapest resolution: run Stage 0
once per required epoch** (each ~200 core-hours), or settle on one epoch now.

## §10.5 STAGE 0 — concrete build plan (starting now)

**Deliverable:** `outputs/epoch_state_cache/epoch_state_<EPOCH>.parquet` — one row per real S3M
object, everything Stage 1 needs, zero propagation downstream.

**Schema**
| column | meaning |
|---|---|
| `ObjID`, `population` | identity (enables §NOTE_FOR_HYAK identity joins for free) |
| `a, e, i, node, argperi, t_p, M_obs_deg` | **osculating elements AT the epoch** (GMM input space) |
| `ra_deg, dec_deg, dra_deg_day, ddec_deg_day` | corrected-observer geometry (§9.2) |
| `lam_deg, beta_deg, vlam, vbeta` | ecliptic position + rates (validated manual rotation) |
| `r_sun_au, delta_au, phase_deg` | geometry |
| `mag_app` | Johnson V, from the **same** state (§9.8) |
| `H` | Johnson V absolute mag (from `.s3m`) |

**Scripts** (in `codex_integrator/`, promoted to `neomod/pipeline/` once validated)
- `stage0_epoch_state.py --pop {neo,mba,tno,trojan} --shard i --nshards N --epoch <iso> --out <p>`
- `slurm/stage0_epoch_state.sbatch` — array over shards; `--account=astro --partition=ckpt-all`,
  `conda_prep/bin/python` (**not** `conda activate`; **not** `--account=dirac` — both were bugs
  before, see `NEOMplanHYAK.md`).
- `stage0_merge.py` — concat shards → the cache; report per-population counts + mag distribution.

**Sharding / cost**
| population | objects | shards | ≈ per shard |
|---|---:|---:|---|
| MBA | 13,883,361 | 400 | ~35k obj ≈ 29 min |
| NEO | 268,511 | 20 | ~13k obj ≈ 11 min |
| Trojan | 179,882 | 10 | ~18k obj ≈ 15 min |
| TNO | 48,682 | 5 | ~10k obj ≈ 8 min |
| **total** | **14.4M** | **435** | **~200 core-hours → ~30 min wall** |

**Two optimizations to measure before launching the full run**
1. **Batched integration.** Current code builds a `rebound.Simulation` **per object**. All S3M objects
   share `t_0 = 54466`, so many massless test particles can go in **one** sim. Plausibly 10–50×.
   *Risk:* IAS15's adaptive timestep is global — one deep encounter could throttle the batch. Measure
   at batch = 1 / 100 / 1000; keep per-object as fallback.
2. **H pre-filter.** `mag_app ≈ H + 5log₁₀(r·Δ) − 2.5log₁₀(Φ)`; mag bins stop at 25. For MBAs
   (r≈2.5, Δ≈1.5 ⇒ +2.9 mag) anything with `H ≳ 22.1` can never reach mag < 25. A **conservative**
   bound (use each population's minimum plausible `r·Δ`) could remove a large fraction of the 13.9M
   MBAs before propagating. **Must be provably conservative** — verify no object within 0.5 mag of the
   cut is dropped, or skip it. Cheap insurance: run the full 14.4M anyway if the bound is at all shaky.

**Acceptance criteria (pre-registered, before results)**
1. **Re-validate `S0000vowa`** through Stage 0 at the night-61642 epoch → must reproduce the §9.1 PASS
   (< 0.01° vs Sorcha; expect ≈ 0.0051°, or ≈ 0.12″ if light-time is enabled).
2. **TNO non-regression:** TNOs already agree at 0.04° under two-body — they must **not** get worse.
   This is the frame/unit canary; a NEO-only test would miss it.
3. **Population sweep** on the 3,205 identity-matched night-61642 objects: median sky separation must
   collapse for all three populations (from 0.040° / 1.401° / 2.271°), NEO improving most.
4. **Magnitude cross-check:** `mag_app` vs adam's `predicted_magnitude_v` — agreed to **1e-4 mag** on
   `S0000vowa` (§9.10); must hold across a sample.
5. **Sanity:** no NaNs; row count = input; `|r|` within physical bounds per population.

**Execution order**
1. Measure batching (1/100/1000) + validate the H pre-filter bound. *(decides shard sizing)*
2. Write `stage0_epoch_state.py`; run **one shard** end-to-end; check acceptance #1, #2, #5.
3. Launch the 435-task array; merge; full acceptance (#3, #4).
4. **Only then** Stage 1 (map regen) — a separate plan once Stage 0's cache is validated.

**Risks**
| risk | mitigation |
|---|---|
| batching throttled by a deep encounter | measure first; per-object fallback (known-good, 0.05 s/obj) |
| H pre-filter drops a real object | conservative bound + 0.5 mag margin check; else no filter |
| 14.4M-row cache size | ~20 cols × 14.4M ≈ 1–2 GB parquet — trivial vs the 20 GB map dir |
| ckpt preemption | shard-level idempotence: skip-existing per shard, re-`sbatch` to refill |
| wrong epoch chosen | §10.4 — settle the epoch **before** launching |

## §10.6 What is NOT in scope for Stage 0
- Stage 1 (map regen) — planned separately once the cache validates.
- Porting n-body into `gmm.py` for production — the adam module is a copy of **`fast.py`** and lacks
  the **GMM cloner** and `from_npz(support_mask_min=...)`, which production requires. Stage 0 sidesteps
  this (it only needs propagation + geometry), but the port is still owed before production scoring.
- Light-time: currently **off** in the direct-assist path (0.005° vs 0.12° with it). Decide whether
  Stage 0 emits light-time-corrected positions — **Sorcha does apply it**, so parity argues yes.

## §10.7 STAGE 0 — MEASURED RESULTS + ACCEPTANCE (2026-07-17)

Built: `codex_integrator/stage0_epoch_state.py`.

### Both planned optimizations were MEASURED and REJECTED

| optimization | prediction | **measured** | verdict |
|---|---|---|---|
| batched integration (N particles / sim) | "plausibly 10–50× faster" | **N=25 → 85 ms/obj; N=200 → 205 ms/obj** vs **50 ms/obj per-object** | ❌ **counterproductive** |
| conservative H pre-filter | "could remove a large fraction of the 13.9M MBAs" | **cuts 0 / 14,380,436 (100% kept)** | ❌ **worthless** |

- **Batching is worse, and worsens with batch size** — exactly the risk flagged in §10.5: IAS15's
  adaptive timestep is **global**, so one particle in a deep encounter throttles the entire batch and
  everyone pays the smallest step. Per-object is optimal: each object takes large steps unless *it*
  has an encounter. **Keep per-object.**
- **The H pre-filter cuts nothing** because S3M is *built* as a survey-detectable catalogue — nothing
  is provably fainter than mag 25 in its most favourable geometry, and any `q ≲ 1` object has a
  `-inf` floor by construction. Left in as `--h-prefilter` (opt-in), **off by default**.

**Net: the §10.2 cost estimate stands unchanged — 14.4M objects × 0.05 s ≈ 200 core-hours.**
Measuring rather than assuming saved implementing two useless optimizations.

### Stage 0 acceptance — **PASS**

`stage0_epoch_state.py --pop neo --objid S0000vowa --epoch 61642.391188036 --epoch-scale utc`

| criterion | result |
|---|---|
| **#1 `S0000vowa` vs Sorcha** | **0.005086399°** — matches the §9.1 reference to **9 decimals**. **PASS** (bar < 0.01°) |
| **#4 magnitude cross-check** | `mag_app` **23.2896** vs adam `predicted_magnitude_v` **23.2895** → **1e-4 mag** ✅ |
| **#5 sanity** | 0 non-finite rows; `r_sun/Δ` = 1.229118 / 0.295677 AU vs Sorcha 1.229197 / 0.295702 ✅ |
| #2 TNO non-regression | **TODO** |
| #3 population sweep (3,205 matched) | **TODO** |

Stage 0 reimplements the geometry inline (rather than calling `build_visible_subset_dataframe`) and
still reproduces the validated number **exactly** — a strong independent check of the port.

### The empirical proof that the regen is necessary

| | a (AU) | e | i (°) |
|---|---|---|---|
| S3M @ 2008 (`t_0`) | 2.269694 | 0.692393 | 18.8960 |
| **n-body @ 2027 epoch** | **2.266195** | **0.694398** | **18.7502** |
| Δ over 19.6 yr | **−0.0035** | **+0.0020** | **−0.146** |

**The elements moved.** Two-body would have frozen them at the 2008 values *by construction*. This is
the concrete refutation of the withdrawn §7 argument ("two-body conserves (a,e,i), so the map ensemble
is fine") — it is backwards: **n-body does not conserve them, and that evolution is exactly what the
map must encode**. A two-body map represents the 2008 population, not the population at the
observation epoch. Advisor's call (regen everything) is empirically justified.

### Time-scale trap (fixed, worth knowing)
Map epochs are **TDB**; Sorcha visit times (`mjd0_utc`) are **UTC** — ~69 s apart, ≈0.0018° for a
2.2°/day object (18% of the PASS bar). `stage0_epoch_state.py` therefore takes an explicit
`--epoch-scale {tdb,utc}` and accepts ISO **or** bare MJD; everything downstream derives from the
single `t_obs` object (elements use `t_obs.tdb.mjd`).

### Remaining before the 435-task launch
1. **#2 TNO non-regression** — TNOs already agree at 0.040° under two-body; they must not get worse.
   This is the frame/unit canary a NEO-only test cannot catch.
2. **#3 population sweep** on the 3,205 identity-matched night-61642 objects.
3. **Settle the epoch** (§10.4) — Stage 0 is *not* epoch-independent, so this must be fixed before launch.
4. **Light-time decision** (§10.6) — off today (0.005° vs 0.12° with it). Sorcha applies it ⇒ parity argues yes.

## §10.8 DECISIONS SETTLED + STAGE 0 LAUNCHED (2026-07-17)

### Decisions (user)
| decision | resolution |
|---|---|
| **Map epoch** (§10.4) | **MJD 61642 = 2027-08-25** — the busiest Sorcha night in *every* case (case1 4,316 / case2 4,490 / case3 4,353 tracklets). Stage 0 is not epoch-independent, so this is baked into the cache. |
| **Light-time** (§10.6) | **ON.** Sorcha applies it (`Range_LTC_km`, `Obj_Sun_*_LTC_km`), so parity requires it. |
| **Concurrency** | **Max out.** |

### Light-time — implemented WITHOUT re-propagating, and cross-validated
Iterating the light-time by re-propagating n-body would ~3× the cost. Instead the position is
corrected **linearly** with the velocity already in hand (2 iterations on `lt`):
over `lt ≈ 148 s` the object moves ~4,400 km, while orbital curvature contributes only ~66 m
(`½·a·lt²` at 1.2 AU) — **0.3% of the 0.12″ (~25 km) residual**, far below our precision floor.

| | value |
|---|---|
| `light_time_s` | **147.535 s** — *identical* to adam's `generate_ephemeris` (independent implementation) |
| `S0000vowa` separation | **0.1273″** (was **18.31″** without LTC → **144× better**) |
| adam's full-iteration LTC | 0.1206″ — the 0.007″ gap **is** the linear-vs-iterated approximation, as predicted |

### ACCEPTANCE #2 — TNO non-regression: **PASS**
The frame/unit canary a NEO-only test cannot catch (TNOs already agreed at 0.040° under two-body,
so a bad rotation or AU/km slip would make them *worse*):

| | median | max |
|---|---|---|
| **n-body + light-time (6 TNOs, night 61642)** | **0.072″** | 0.143″ |
| two-body TNO baseline | **144″** (0.040°) | — |

TNOs improved **~2000×** and land in the same 0.04–0.14″ band as the NEO (0.127″) — consistent
across populations. **No regression; frame/unit handling confirmed correct.**

### Final sizing — maxed out (measured, not guessed)
| finding | consequence |
|---|---|
| binding limit is ckpt QOS **`MaxSubmitPU = 2000` tasks**, *not* CPUs (ckpt-all had **11,636 idle**) | size the array to just under 2,000 |
| ASSIST integration is **single-threaded per object** | **`cpus-per-task=1`** — extra cores are wasted; 1 cpu/task doubles task density |
| every MBA task was loading **all 14 S1 files (13.9M rows)** to take its slice | **shard MBA by file** (`--pattern S1_XX.s3m`) → ~1M rows/task, **14× less I/O** |
| batching **counterproductive** (§10.7) | per-object |

**Layout — 1,886 tasks × 1 cpu:**
```
   0-1819  MBA     14 files x 130 sub-shards   13,883,361  ~7.6k/task
1820-1854  NEO     35 shards                      268,511  ~7.7k/task
1855-1878  Trojan  24 shards                      179,882  ~7.5k/task
1879-1885  TNO      7 shards                       48,682  ~7.0k/task
```
**200 core-hours → ~10 min wall** at full concurrency. Idempotent (skip-existing per shard), so a
re-`sbatch` refills anything ckpt preempts.

**Launched:** job **37221827**, `outputs/epoch_state_cache/2027-08-25T000000/shards/`.

### Acceptance scoreboard
| # | criterion | status |
|---|---|---|
| 1 | `S0000vowa` vs Sorcha | **PASS** — 0.127″ (LTC on) / 0.005086399° (LTC off, exact vs §9.1) |
| 2 | TNO non-regression | **PASS** — 0.072″ median vs 144″ two-body |
| 3 | population sweep | **PASS** — see §10.9 (0.0001°/0.0012°/0.0024°; 600–1200× better) |
| 4 | magnitude cross-check | **PASS** — 1e-4 mag vs adam `predicted_magnitude_v` |
| 5 | sanity (NaNs, row counts, \|r\|) | **PASS** — 0 non-finite |

---

## §10.9 — Stage 0 EXECUTED, cache built, acceptance #3 PASS

**All five acceptance criteria now pass. Stage 0 is done; Stage 1 (map regen) is unblocked.**

### Stage 0 run
Job **37221827**: 1,886/1,886 tasks COMPLETED, zero failures, ~6.5 min/task, zero ckpt
preemptions — the ~10 min wall estimate held.

Merge (job **37224668**, 128 G node — the login node OOM-killed the Kurlander merge at 19 GB,
so this was deliberately put on a node) →
`outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet`, 3.5 min:

| population | rows | expected |
|---|---:|---:|
| MBA | 13,883,361 | OK |
| NEO | 268,511 | OK |
| Trojans | 179,882 | OK |
| TNO | 48,682 | OK |
| **total** | **14,380,436** | **COMPLETE** |

Zero non-finite in `ra_deg / vlam / mag_app / a / e / light_time_s`; zero duplicate ObjID.
`mag_app` median 26.06 (28.7 % brighter than 25); `light_time_s` median 1470.8 s;
`e` ∈ [0.0000, 0.9996]; `a` ∈ [0.310, 72.511] AU. Every real S3M object now carries its
n-body osculating elements at MJD 61642 plus geometry and magnitude from the *same* state,
so Stage 1 never propagates.

### Acceptance #3 — population sweep (`codex_integrator/acceptance_population_sweep.py`)

| pop | n | raw sep | motion-corrected | two-body | improvement |
|---|---:|---:|---:|---:|---:|
| TNO | 51 | 0.0031° | **0.0001°** | 0.040° | 613× |
| MBA | 3,310 | 0.0318° | **0.0012°** | 1.401° | 1179× |
| NEO | 820 | 0.0718° | **0.0024°** | 2.271° | 965× |

Objects >5° from Sorcha: **0** (two-body: 481 = 15 %). `|Δv|` median 0.0078 °/day — the
two-body bulk was 0.0106, and 0.0078 was what two-body achieved *only on its already-good
sep<0.5° subset*. The whole population is now at the accuracy two-body reached only where it
happened to work.

**Two honest caveats on this table.**

1. **The matched set is larger than Arnor's, not the same set.** We match 4,316 objects —
   every night-61642 case1 tracklet — where Arnor's two-body sweep matched 3,205. The cache
   covers all 14.4 M S3M objects, so nothing drops out; the two-body comparison could only
   match objects its own visible-subset step retained. So this is *n-body on a superset* vs
   *two-body on a subset*, which if anything understates the improvement.
2. **"Motion-corrected" is doing real work and must be read carefully.** The cache is at the
   map epoch (MJD 61642.0 TDB); each tracklet has its own visit time up to ~1 d away (median
   gap 0.189 d). We advance the cache position with its own RA/Dec rates before comparing —
   otherwise n-body is charged for ~1° of genuine NEO motion, which is the size of the effect
   being measured. Raw separations are printed alongside for exactly this reason.
   (TDB-vs-UTC in Δt is ~69 s ≈ 0.0016° for the fastest NEO — below every threshold here.)

**The leftover residual is my linear extrapolation, not the integrator.** Residual binned by
epoch gap: 0.80″ (Δt<0.05 d) → 1.61″ → 4.16″ → 5.72″ (Δt≈0.3 d), i.e. it *grows with the gap*,
the signature of neglected curvature in a linear rate extrapolation. The zero-gap floor is
sub-arcsecond, consistent with the independent 0.127″ of acceptance #1. (A global `k·Δt²` fit
gives k≈18″/day² but correlates only 0.154 — curvature is per-object, so don't read that
intercept; the small-Δt bin is the trustworthy floor estimate.)

### Where this leaves the regen
Stage 0 is the **entire n-body cost** of the map regen (§10.2: clones are osculating elements
*at* the epoch, so elements→state is a zero-gap Kepler solve, exact by definition). That cost
is now paid and cached.

**Next, and it needs its own design pass before launching:**
1. **Stage 1** — map regen from the cache (§10.6), the 435-task launch.
2. **Port n-body into `neomod/src/velocity_density_pipeline_gmm.py`** — the production module
   all three entry points import. `velocity_density_pipeline_adam.py` is a copy of `fast.py`
   and lacks the GMM cloner and `from_npz(support_mask_min=…)`; it cannot be the production
   path. **`gmm.py` still carries the live origin bug at ~line 527.**
3. **§10.4 scope** — rebuild benchmark v1/v2/v3 + night61642, re-score case1/2/3 + v5,
   recompute ROC/F1.

---

## §11 — The new integrator in the production VDP module (`velocity_density_pipeline_gmm.py`)

**Written 2026-07-17.** §3–§10 built and validated the n-body integrator in a *scratch copy*
(`codex_integrator/velocity_density_pipeline_adam.py`, itself a copy of `..._fast.py`). That copy
lacks the GMM cloner and `ProbMapSet.from_npz(support_mask_min=…)`, so it can never be the
production path. **`velocity_density_pipeline_gmm.py` is the production module** — all three entry
points import it, and it is the file that both **generates the maps** and **scores** against them.
This section documents the port of the validated fix into `gmm.py`, function by function, so it can
be read on its own (e.g. by Arnor when building the Horizons regression notebook).

### Why this module and not `fast`/`adam`
`fast.py` and `gmm.py` are ~95 % identical; the one load-bearing difference is the **NEO GMM
cloner** (`_gmm_*`, `_clone_neo_gmm`) plus the from-npz support mask. That difference is exactly what
made the port non-mechanical (see the clone-path branch below). Everything else — the geometry
layer, the sky cut, the manual equatorial→ecliptic rotation — is byte-identical between the two, so
the *propagation* changes transplant cleanly; only the *clone* interaction needed thought.

### What changed, where it came from, and why

All four changes are ports of the `adam.py` code that §9.1/§9.9/§9.10 validated, with two
deliberate deviations noted at the end.

**1. Imports + `SORCHA_PLANETS_PATH` constant** (`gmm.py:~80`, `:134`)
- Added `from pathlib import Path`; added `get_body_barycentric_posvel` and
  `solar_system_ephemeris` to the astropy import (needed by the origin fix, which runs on *every*
  call, so these are top-level).
- `SORCHA_PLANETS_PATH` points at Sorcha's own ASSIST planetary kernel
  `sorcha_cache_2025-07-06/linux_p1550p2650.440`, resolved via `Path(__file__).parents[2]`
  (**note:** `parents[2]` from `neomod/src/`, vs `parents[1]` in adam's `codex_integrator/` — the
  only path that differs by location).

**2. `propagate_elements_nbody(...)`** (`gmm.py:530`) — **new function, copied verbatim from adam**
(verified **bit-identical output**: `max|Δr| = 0`, `max|Δv| = 0` over 200 real S3M NEOs).
- Purpose: propagate osculating heliocentric-ecliptic elements at `t_0` to `obstime` with ASSIST
  n-body, returning `(r_helio_ecl_km, v_helio_ecl_km_s)` — the **exact output contract** of the old
  `nsc.elements_to_helio_ecliptic_state`, so the geometry layer downstream is untouched.
- Mechanism (per object, per §9): Sorcha's `universal_cartesian()` COM initialization →
  `ecliptic_to_equatorial` → a REBOUND `Simulation` with `assist.Extras`, `ri_ias15.adaptive_mode=1`
  → `sim.integrate(target)` → subtract the Sun to return heliocentric → `equatorial_to_ecliptic`.
- `gr_model="GR_SIMPLE"` (default) reproduces Sorcha's own downgrade
  (`sorcha/ephemeris/simulation_setup.py:182-183`); `"GR_EIH"` keeps ASSIST's full model. Measured
  difference 0.0039″ (§9.9.B), so this is moot in practice but exposed for auditing.
- **Deviation from adam (deliberate):** the `assist`/`rebound`/`sorcha`/`jpl_small_bodies` imports
  are **lazy** (inside the function), not top-level. `gmm.py` is imported by consumers that only load
  maps and score *observed* tracklets (no propagation); lazy imports keep the module importable
  without the n-body stack installed. The astropy origin-fix imports stay top-level because the
  origin fix is unconditional.
- Versions used: `adam_core 0.5.5`, `assist 1.2.3`, `rebound 4.6.0`, asteroid perturbers
  `jpl_small_bodies_de441_n16`.

**3. `build_visible_subset_dataframe(...)`** (`gmm.py:662`) — the heart of the change. Three edits:

  a. **Propagation now branches on `t_0`** (`:696`). This is the one place the port is *not* a
     straight copy, and the reason is the cloner:
     - **Real objects** (loaded via `s3m_loader`/`hybrid_loader`, which both supply a `t_0` column)
       carry their osculating epoch — S3M's MJD 54466 (2008). They get **ASSIST n-body** over the
       ~19.6 yr gap. This is the scoring path and the map-gen *source-object* path; it is the entire
       integrator bug (§1).
     - **Clones** (`raw_clone_df`, `gmm_clone_df`) are generated *at* the map epoch — the cloner sets
       `t_p` so mean anomaly is correct at `obstime` — and carry **no `t_0`**. For them the gap is
       zero, so the two-body Kepler solve (`nsc.elements_to_helio_ecliptic_state(method="newton")`)
       is **exact**. Running ASSIST on them would be wasted work *and* would crash, because adam's
       unconditional-n-body version raises without `t_0`. So the branch is both a correctness choice
       (zero-gap two-body ≡ n-body) and what keeps map generation working.
  b. **Origin fix** (`:731`), applied on **both** branches. `NEOMODScorer._get_earth_and_observer`
     returns a **barycentric** Earth state while the propagated asteroid is **heliocentric**; we
     subtract the Sun's barycentric offset (`get_body_barycentric_posvel("sun")`, `de432s`) so
     observer and object share an origin. This is the **second, independent bug** (§9.2): fatal at
     small geocentric distance (~1–2° at Δ = 0.3 AU), invisible at 2 AU.
  c. **`mag_app` from the same state** (`:751`), computed only when `H` is present (real objects;
     clones carry no `H` and don't need per-clone magnitudes — their maps are per-mag-bin). Same
     HG(G=0.15) phase function as the legacy function, fed the corrected state — §9.8's
     "propagate once, derive both" so position and magnitude cannot disagree.

**4. `score_orbital_df(...)`** (`gmm.py:2535`) — **single propagation**. It now calls
`build_visible_subset_dataframe` once and **consumes** the `mag_app` it returns, raising if absent.
The old code called `compute_apparent_magnitude_for_population()` first and propagated a **second**
time — which doubled the n-body cost and, worse, computed the magnitude with the **legacy two-body
propagator and the uncorrected observer** (§9.3.1). It requires the input df to carry `t_0` (n-body
gap) and `H`.

### Two adam bugs I did NOT copy
1. **adam's `build_visible_subset_dataframe` is unconditionally n-body** and raises without `t_0`.
   That is fine for adam (only ever run on the scoring path) but would crash `gmm.py`'s map
   generation on the clone dataframes. → fixed by the `t_0` branch above.
2. **`NameError` in adam's `score_orbital_df`:** adam removed the code that defined `mag_app_vis`
   but still references it (`adam.py:2268`), so that method would raise `NameError` if ever called.
   → the gmm version defines `mag_app_vis = visible_df["mag_app"].to_numpy(float)`.
   *(adam is a frozen scratch reference; it can keep these — the production module is `gmm.py`.)*

### Tests run on the port
| test | result |
|---|---|
| `py_compile` + import (heavy: `neoscore`, `hybrid_loader`) | OK; kernel path resolves; lazy imports confirmed (gmm imports without `assist`/`rebound`/`sorcha`) |
| `gmm.propagate_elements_nbody` vs `adam.propagate_elements_nbody`, 200 real S3M NEOs | **bit-identical** (`max|Δr| = 0 km`, `max|Δv| = 0 km/s`) |

Because the propagator is bit-identical to the adam code, **all of §9.9 applies unchanged** to the
production module — see it for the external-truth validation, restated here for the Horizons notebook:

### External validation already on record (§9.9), for the Horizons regression
These used `astroquery`/`adam_core` against **JPL Horizons** and Sorcha:
- **A — vs Sorcha (`generate_ephemeris` on `S0000vowa`):** two-body 71.95° → n-body no-LTC 18.31″ →
  **n-body + light-time + `X05` topocentric + `GR_EIH` = 0.12″** (light-time was the whole 18″
  residual, `light_time = 147.535 s` at Δ = 0.296 AU).
- **B — `GR_SIMPLE` vs `GR_EIH`:** 0.833 km after 19.6 yr = **0.0039″** — negligible, ~4,700× below
  the light-time term. Force-model choice is moot.
- **C — Apophis (`query_sbdb(["99942"])` + `query_horizons(coordinate_type="cartesian",
  location="@sun")`):** **PASS** — tracks JPL to **~5 km/yr pre-encounter** (1.2 km at +1 yr,
  38.2 km at +8 yr). The 38 km → 22,870 km jump across the 2029 encounter is **physics, not a bug**
  (a deep encounter exponentially amplifies the tiny SBDB-vs-Horizons initial-solution difference —
  the same chaos that makes the S3M NEO problem hard). S3M objects are synthetic (no Horizons
  entry), so JPL validates the **code**; Sorcha-parity is the standard for the **synthetic
  benchmark**.
- **adam_core traps** (so the notebook doesn't re-hit them): `tp` is **MJD not JD** in
  `CometaryCoordinates` (direct-assist path wants JD — opposite convention); `Orbits` needs
  **Cartesian** (`coords.to_cartesian()`); `ASSISTPropagator` exposes **no force API** (do force A/Bs
  in the direct-assist path); **`query_horizons_ephemeris` does not exist** in 0.5.5 → use
  `query_horizons`.

### §11.1 — The Stage 0 cache: an OPTIONAL efficiency path for map generation (a decision, not a blocker)

**This is deferred and does not affect the Horizons test.** Recording the decision so it isn't
re-derived.

**Where the cost is.** In map generation, `generate_probability_maps` →
`build_cloned_maps_for_center_magbin` calls `build_visible_subset_dataframe` **on the real source
population** (`df_cloner_input`, `gmm.py:~1380`) to find which real objects fall in the sky patch,
*before* cloning. With the new integrator that call is now **n-body**, and it happens **once per
(population × mag-bin × map-center)**. Summed over mag-bins it is the whole population; times the
grid of map centers it is the whole population **× N_centers** n-body propagations. That is the
~10⁴-node-hour cost that made naïve n-body map generation look infeasible.

**Why Stage 0 exists.** Every real object's n-body state at the map epoch is **the same regardless of
which map center we are building** — the propagation does not depend on the center, only the sky cut
does. So Stage 0 propagated all 14,380,436 objects to the epoch **once** and cached the result
(`outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet`: ra/dec, vlam/vbeta, mag_app,
elements). Map generation can then **look up** each object's epoch state and apply the (cheap) sky
cut per center, instead of re-propagating. Clones are still generated per center and are zero-gap
(fast) either way.

**The design: make it an OPTIONAL argument, default = propagate inline.**
- `build_visible_subset_dataframe` (and `build_cloned_maps_for_center_magbin`) gain an optional
  `precomputed_states=<cache>` argument.
- **If omitted (default):** behave exactly as now — propagate inline. Scoring, the Horizons test,
  and any one-off keep working unchanged. **Correctness is identical** because the cache *is* the
  n-body result.
- **If supplied:** skip `propagate_elements_nbody`, take r/v (or ra/dec + rates) from the cache,
  apply only the sky cut. Used solely by the Stage 1 map-regen driver.
- This is a pure **efficiency/feasibility** change, **not correctness** — the two paths produce the
  same numbers; one is just N_centers× cheaper.

**The decision for you (Devanshi):**
1. **Before sending to Arnor / for the Horizons test — do nothing.** That test propagates a handful
   of objects inline and compares to Horizons; it never touches the cache. `gmm.py` as-is is correct
   and complete for it.
2. **For Stage 1 (map regen at scale) — the cache path is effectively required** for feasibility,
   and should be added then as the optional argument above. It is safe to defer because default
   behavior is unchanged; adding it later cannot regress scoring or the Horizons test.

**Recommendation:** ship `gmm.py` as-is now (inline n-body, validated, bit-identical to adam); add
the optional `precomputed_states` hook as the first task of Stage 1, alongside the driver that
consumes it. Not before — untested cache scaffolding in the file Arnor reads would only add noise to
the integrator he is validating.

---

## §12 — STAGE 1: map regeneration (n-body), starting with one validatable center

**Written 2026-07-17.** §11 put the n-body integrator into the production module. Stage 1 rebuilds
the probability maps with it. Per the advisor, once the integrator changes we regen **all** maps and
**all** scoring (§10) — but we do **one center first**, epoch-matched to a night we can check against
Sorcha, before spending the full grid.

### 12.1 What actually has to be regenerated (and what does not)

The s3m-linking scoring pipeline (`neomod/pipeline/slurm/s3m_linking/`, `_case_env.sh`) uses **two**
map directories:

| dir | role | regen? |
|---|---|---|
| `prob_maps_grid` | **geometry only** — center lon/lat/obstime metadata, used by stage 2 to assign each tracklet to its nearest center | **no** — centers don't move |
| `prob_maps_grid_s3m` | **the densities** — the vlam/vbeta probability maps VDP scores against (stage 3) | **YES** — these are two-body; this is the whole point |

So Stage 1 = rebuild **`prob_maps_grid_s3m`** with the n-body integrator. The 667-center grid geometry
(`sorcha_gen_maps_grid.py`, `build_grid`) is unchanged. Tracklets themselves (stages 1–2 of the case
pipeline) do **not** change — only the maps and therefore the scores. So the case re-run is stages
**3–6 only** (score-vdp → sample → digest2 → combine), not a Sorcha re-run.

### 12.2 The one code change first: wire the Stage 0 cache into map generation

Without this, each of the 667 map tasks would n-body-propagate ~14.4M objects **per center** (§11.1)
— 667× the Stage 0 cost. The fix (detailed in the chat log / §11.1) is **one branch** at the top of
`build_visible_subset_dataframe`:

```python
has_positions = {"ra_deg","dec_deg","vlam","vbeta"}.issubset(df.columns)
if has_positions:
    # cache path: Phases 1-2 (propagate + geometry) already done in Stage 0; read them off
    ra_deg, dec_deg = df["ra_deg"].to_numpy(), df["dec_deg"].to_numpy()
    # vlam, vbeta, mag_app likewise straight from df
else:
    ...existing propagate (n-body if t_0 else zero-gap two-body) + geometry...
# Phase 3 (sky-cut) runs either way -- the only center-dependent step
```

- **Source objects** (map-gen) come from the cache → carry positions → propagation skipped, only the
  per-center cut runs. Cloner needs their *elements* (in cache) and the map needs *vlam/vbeta* (in
  cache), so the returned df is column-identical to the propagated one.
- **Clones** are raw elements (no positions) → fall through to the existing zero-gap two-body path,
  unchanged.
- **Scoring / Horizons test** pass raw elements (no positions) → full n-body path, unchanged →
  **cannot regress**. This is why 12.2 is safe to land before any map is built.

Plus a small hook in `sorcha_gen_maps_grid.py`: load
`outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet` once and feed its rows as the source
population, instead of loading raw elements and letting each task propagate.

**STATUS — the library half is DONE and tested (2026-07-17).** `build_visible_subset_dataframe` now:
(a) has a `_center_skycoord` helper shared by both paths; (b) takes an early cache branch when the input
carries `ra_deg/dec_deg/vlam/vbeta`. Verified: cache path selected 42,609 visible of 268,511 NEOs in
**0.06 s with no propagation and `scorer=None`**; the normal (raw-element) path is unchanged after the
center refactor (300 clones → 300 visible, finite vlam/vbeta). The driver hook (feeding the cache into
`generate_probability_maps`) is the remaining piece, done as part of Milestone A.

**Correctness invariant:** the cache used the same propagator + origin fix as the live function, so a
cache-built map is numerically the same as an inline-propagated one — a *speed* substitution, not a
*physics* one. (Footnote: the cache also applied the ~18″ light-time correction the inline path skips;
negligible for a 30° cut and element-space cloning.)

### 12.3 MILESTONE A — one center, epoch-matched, validated against Sorcha

**Purpose:** prove the whole n-body map path end-to-end (cache → clone → KDE → map → score) on a
single center before committing the 667-map array, and do a clean **n-body-map vs two-body-map**
head-to-head on identical Sorcha tracklets.

**The center (measured from `case1` night-61642 tracklets):**
- Night 61642 tracklets span **224 centers**; the busiest is **`grid_dlon-010_lat+08`**, absolute
  ecliptic **(lon 90.2°, lat +8°)**, holding **80 tracklets = 61 MBA + 19 NEO**.
- At the epoch (MJD 61642, antisun lon ≈ 331.2°) that absolute center is **dlon ≈ +119° from
  antisun** — i.e. ~61° solar elongation, *off-opposition*, where NEO/MBA velocity separation is
  inherently weaker (WAGG §"antisun always separates"). Fine for a *relative* integrator check; not a
  place to judge VDP's absolute strength.
- **Honest caveat:** 80 tracklets (19 NEO) is a **smoke test**, not a powered ROC. Statistical weight
  comes at Milestone B (the full grid). If we want more NEOs in the first check, build the **top ~10
  busiest centers** (~500+ tracklets) — still ~65× cheaper than the full grid.

**Steps:**
1. Land 12.2 (cache branch + driver hook); unit-check that a raw-element call is byte-unchanged.
2. Build **one** n-body density map at absolute center (90.2°, +8°), **obstime = 2027-08-25T00:00:00
   (MJD 61642)**, from the cache → `prob_maps_grid_s3m_nbody/…dlon-010_lat+08.npz` (new dir, so the
   two-body maps survive for the head-to-head).
3. Score the **80 night-61642 tracklets** in that footprint with `ProbMapSet.from_npz(...,
   support_mask_min=1)` + `score_observation` (rate conversion `dra = RARateCosDec/cos(dec)`,
   WAGG §"VDP Rate Conversion"). Do the same with the **existing two-body** map at the same center.
4. Compare on identical tracklets:
   - NEO vs MBA `P_NEO` separation (median P_NEO for NEO ≫ MBA);
   - per-tracklet Δ(P_NEO) between n-body and two-body maps — where and why it moves;
   - digest2 (`P_NEO_d2`, already in the parquet) as the fixed reference.

**Success criteria (Milestone A):**
- The map builds from the cache with sane coverage (non-empty clone density, NEO overlay present).
- n-body map scores are **at least as good** at NEO/MBA separation as the two-body map on these 80,
  and any large per-tracklet Δ traces to objects the integrator relocated (§1 along-track cases).
- No pipeline breakage: `score-vdp` on the case tracklets runs against the new map dir unchanged.

### 12.4 MILESTONE B — full 667-grid regen + full re-score (the §10.4 scope)

1. Run `sorcha_gen_maps_grid.py` as a Slurm array over `--task-id 0..666`, `--ref-obstime
   2027-08-25T00:00:00`, source = cache → fill `prob_maps_grid_s3m_nbody/` (667 npz). Cost is now the
   *old two-body grid* cost (clone + KDE per center; no n-body), because the propagation was paid once
   in Stage 0.
2. Re-point the case pipeline at the new maps (`PROB_MAPS_SCORE=prob_maps_grid_s3m_nbody` in
   `_case_env.sh`) and re-run **stages 3–6** for case1/2/3 (+ benchmark v-band / night61642): score-vdp
   → sample → digest2 (unchanged, orbit-based) → combine. Exact resources/arrays/submit order are in
   **`sorcha_full_pipeline.md` §3–§4** — stages 1–2 (Sorcha run + tracklets) are **skipped** (tracklets
   are unchanged by a map swap), and stage 5 (digest2) is orbit-space so its scores are identical; only
   VDP scores move. Turn-key re-score per case:
   ```bash
   cd neomod/pipeline/slurm/s3m_linking;  C=case1   # PROB_MAPS_SCORE already repointed in _case_env.sh
   sbatch --export=ALL,CASE=$C 3_vdp.sh        # 2cpu/32G, array 0-112%113, ~10 min, reads new maps
   sbatch --export=ALL,CASE=$C 4_sample.sh     # prints subsample row count R
   sbatch --export=ALL,CASE=$C --array=0-129%130 5_digest2.sh   # or reuse prior d2 shards (orbit-based, unchanged)
   sbatch --export=ALL,CASE=$C 6_combine.sh    # -> sorcha_comparison_<case>.parquet
   ```
   (Stage 5 can be *reused* from the prior run rather than recomputed — digest2 doesn't see the maps —
   but re-running is cheap and keeps the tree self-consistent.)
3. Recompute ROC / F1 (AUC, best-F1, completeness, contamination) vs digest2, per §10.4. Compare to the
   two-body-map baseline table (WAGG §"Current Status": VDP GMM F1 ≈ 0.837).

### 12.5 Epoch handling

Build at **obstime 61642** so the single-center check is epoch-matched to the night we validate
(eliminates the epoch-mismatch confound that dominated the 5-map 2yr result, WAGG §"Epoch caveat").
The full grid's `--ref-obstime` is also 61642. Whether the *2-year* re-score needs multiple map epochs
(monthly, as the 2yr hybrid run used) is a separate question deferred to Milestone B analysis — the
case1/2/3 comparison is anchored on night 61642, so a single-epoch grid is the right first target.

### 12.6 File / step chain (from WAGG_SORCHA_HYAK_CONTEXT.md)

```
cache (Stage 0)  ──►  sorcha_gen_maps_grid.py (GMM, cache-fed)  ──►  prob_maps_grid_s3m_nbody/*.npz
                                                                          │
case tracklets (stages 1-2, UNCHANGED)                                    │
   outputs/s3m_linking/caseN/tracklets/*.parquet                          ▼
        └─►  sorcha_phase2.py score-vdp  (--prob-maps-dir prob_maps_grid_s3m_nbody,
             (s3m_linking 3_vdp.sh)        --no-nearest-dist-mask --support-mask-min 1)
                                                                          │
        └─►  4_sample.sh  ─►  5_digest2.sh (orbit-based, unchanged)  ─►  6_combine.sh
                                                                          ▼
                                            sorcha_comparison_caseN.parquet  ─►  ROC/F1 notebook
```

### 12.7 Decisions (2026-07-17)
1. **Milestone A width:** **top ~10 busiest centers** (~500 tracklets, more NEOs) — decided, for more
   statistical weight in the first check while staying ~65× cheaper than the full grid.
2. **Map output dir:** **new dir `prob_maps_grid_s3m_nbody/`** — decided; two-body `prob_maps_grid_s3m`
   is kept for the head-to-head and rollback.
3. **2-year vs single-epoch grid** for the full re-score (12.5) — still open; decide after Milestone B.

### §12.3.1 — Milestone A first build: a center-placement bug (NOT an integrator bug), found and fixed

The first 10-map build (job 37278453) produced maps that scored **worse than random** (n-body
AUC 0.469, F1 0.347 vs the two-body baseline AUC 0.893, F1 0.833) — the n-body map called *every*
tracklet a NEO. Diagnosed in full; the integrator is fine, the map placement was wrong.

**What was wrong.** The grid is **antisun-relative**: `grid_dlon-010` means −10° from the antisun
*at the map's ref-obstime*. My builder reused the parquet's `prob_map_center_lon_deg` (= 90.2°), which
is the absolute center at the **old** ref-obstime 2026-01-01. At the map epoch (MJD 61642) the antisun
has moved to 331.2°, so `dlon-010` belongs at **321.2°**, not 90.2°. The maps were built **233° away
from the tracklets they were meant to score.**

**The evidence trail (all recorded so the reasoning is reproducible):**
1. n-body map MBA density centroid was at vlam **+0.35**; the MBA tracklets sit at **−0.23**.
2. The zero-gap two-body clone path **exactly** reproduces the cache velocities (`vlam diff = 0.0000`)
   — so `build_visible` is not the culprit.
3. **Identity-matched** (same ObjID) cache-vs-tracklet velocities **agree**: MBA cache vlam −0.215 vs
   tracklet −0.228, |Δv| 0.0079, sky-sep 0.035° — **the n-body integrator is correct** (acceptance #3
   holds at the object level).
4. The cache velocity **field** resolves it: at the (wrong) map center (90.2°, 8°) the cache MBA vlam
   is **+0.33**; at the **tracklet centroid (323.4°, 8.4°)** it is **−0.22** — matching the tracklets.
   The tracklets were 233° from where the map was built.

**The fix** (`stage1_build_maps_from_cache.py`): parse `dlon`/`lat` from the grid label and place the
absolute center at **antisun_lon(map_obstime) + dlon** (exactly what `sorcha_gen_maps_grid.py` does),
instead of reusing the stored 2026-01-01 absolute center. Verified: new centers land on the tracklet
centroids (dlon-010 → 321.2° vs centroid 323.7°; dlon-020 → 311.2° vs 314.3°). Rebuilding (job 37281202).

**Lesson for Milestone B:** the full 667-grid regen must be driven through `sorcha_gen_maps_grid.py`'s
antisun-relative `--task-id`/`--ref-obstime` path (which is already correct), *not* by reusing stored
absolute centers. The cache hook feeds that driver; it does not replace its geometry.

### §12.3.2 — Milestone A RESULT: PASS (2026-07-18)

10 maps rebuilt at the corrected antisun-relative centers (job 37281202, COMPLETED 53 min).
Validation (`stage1_validate_milestoneA.py`) on the **584** night-61642 tracklets in those
footprints (**111 NEO / 473 non-NEO**, 10 centers):

| classifier | AUC | best-F1 | NEO median P | MBA median P |
|---|---|---|---|---|
| **n-body map (new, epoch-matched 61642)** | **0.899** | **0.849** | 1.000 | 0.0011 |
| two-body map (old, 2026-01-01) | 0.897 | 0.832 | 1.000 | 0.0012 |
| digest2 (orbit-based) | 0.965 | 0.865 | — | — |

**Verdict: PASS** against the §12.3 criteria: (1) maps build from the cache with sane coverage;
(2) the n-body map is **at least as good** as the two-body map on identical tracklets (F1 0.849 vs
0.832, AUC 0.899 vs 0.897 — a marginal *gain*, well within noise at n=584 but clearly **not a
regression**); (3) no pipeline breakage. Clean NEO/MBA separation (P_NEO 1.000 vs 0.0011).

**Reading it honestly.** 584 tracklets / 111 NEO is a smoke test; the ΔF1 of +0.017 is not
statistically meaningful at this scale. What it *does* establish: the full n-body path
(cache → clone → KDE → map → score) is correct end-to-end, self-consistent (n-body scoring against
n-body maps, as the advisor required), and does not degrade map quality. The powered ROC is
Milestone B on the full 667-grid + full re-score.

**End-to-end path now proven.** Stage 0 (n-body cache) → §12.2 cache branch → cache-fed map gen at
antisun-relative centers → epoch-matched scoring. Ready for Milestone B.

---

## §12.4.1 — MILESTONE B execution plan (case1 only), 2026-07-18

Decision (Devanshi): add the cache hook to the **production** driver, launch Milestone B for
**case1 only**. Recorded here as the live plan.

### Step 1 — cache hook in `sorcha_gen_maps_grid.py` (DONE)
Added `--cache <epoch_state.parquet>`. When given, the driver builds `clone_sources` from the cache
(df per population = cache rows; `_mag_app` = cached magnitude → no two-body mag re-propagation; one
**shared** scorer, since it is used only for the population-independent observer geometry, avoiding a
13.9M-row MBA S3M load per task) and passes it to `generate_probability_maps` instead of
`population_settings`. The driver's antisun-relative center placement (`antisun_lon(ref_obstime)+dlon`)
is already correct (§12.3.1). **`--ref-obstime` MUST equal the cache epoch** (2027-08-25) or the cache
positions won't match the observer geometry. Compiles; grid still 667; `--list-only` intact.

### Step 2 — full 667-grid n-body build (array)
`codex_integrator/slurm/stage1_build_grid667_nbody.sbatch`: `--array=0-666`, one center per task via
`--task-id`, `--ref-obstime 2027-08-25T00:00:00`, `--cache …`, `--prob-maps-dir prob_maps_grid_s3m_nbody`,
`--mba-clone-factor 5`, 4 cpu / 32G / 1.5h. This is the right way to "go to the max": 667 light tasks
fan out across ckpt-all rather than one serial job. Preceded by a 1-task smoke test (job 37293616).
Cost is the *old two-body grid* cost — the n-body propagation was paid once in Stage 0.

### Step 3 — re-score case1 (the comparison)
`score_vdp_frame` (stage 3) groups tracklets by their **stored** `prob_map_file` (antisun-relative
assignment from stage 2, unchanged) and loads that map from the given dir. So pointing it at
`prob_maps_grid_s3m_nbody/` scores each tracklet against the n-body map for its cell. Two readouts:

- **PRIMARY — night-61642, powered & epoch-matched.** Score all **4,316** night-61642 case1 tracklets
  (across their ~224 centers) with the n-body maps, vs the stored two-body `P_NEO_vdp_Vband` and
  `P_NEO_d2`. This is Milestone A scaled from 10 centers to the full night — a powered, epoch-matched
  ROC with no epoch-mismatch confound. `stage1_validate_milestoneA.py` already does exactly this; it
  just needs the full 667 grid present (all 224 centers).
- **SECONDARY — full 2-yr stages 3–6 (follow-on).** Re-score every case1 tracklet against the single-
  epoch (61642) antisun-relative n-body grid — methodologically identical to the original two-body
  single-epoch (2026-01-01) grid, so apples-to-apples. Needs the V-band scoring plumbing
  (`mean_mag_V`) and a separate output (`sorcha_comparison_case1_nbody.parquet`) to preserve the
  two-body baseline. Deferred until the PRIMARY night-61642 result is in.

### Step 4 — ROC/F1
AUC, best-F1, completeness, contamination for n-body vs two-body vs digest2, on the night-61642 set
(primary) and — if run — the full-2yr set. Compare to the two-body baseline (WAGG "Current Status":
GMM F1 ≈ 0.837).

### Scope guard
**case1 only** for now (case2/case3 deferred per Devanshi). The 667-grid is case-independent (built
from the S3M cache, not from tracklets), so once built it serves all three cases; only the re-score
(step 3) is per-case.

### §12.4.2 — MILESTONE B primary result: full 667-grid built, night-61642 re-score PASS (2026-07-18)

**667-grid built**: array 37294079, **667/667 COMPLETED, zero failures**, ~6 min/task,
`prob_maps_grid_s3m_nbody/` — all n-body, antisun-relative at epoch 61642, via the production driver's
cache hook. Preceded by a 1-task smoke test (job 37293616) that confirmed the hook end-to-end.

**Powered night-61642 re-score** (`stage1_validate_milestoneA.py`, all **4,316** case1 night-61642
tracklets across **224** centers, **820 NEO / 3,496 non-NEO**):

| classifier | AUC | best-F1 | NEO median P | MBA median P |
|---|---|---|---|---|
| **n-body map (new, epoch-matched 61642)** | **0.893** | **0.793** | 1.000 | 0.0013 |
| two-body map (old, 2026-01-01) | 0.883 | 0.788 | 1.000 | 0.0013 |
| digest2 (orbit-based) | 0.931 | 0.814 | — | — |

**Verdict: PASS.** On the full powered, epoch-matched night, the n-body map is **marginally better
than the two-body map** (AUC +0.010, F1 +0.005) and **not a regression** — consistent with Milestone
A (10 centers) now confirmed at 224 centers / 4,316 tracklets. Clean NEO/MBA separation (P_NEO 1.000
vs 0.0013). digest2 still leads (0.814), as in all prior comparisons.

**How to read the small margin.** The headline deliverable of the regen is **correctness /
self-consistency** — the advisor's requirement that n-body scoring run against n-body maps — not a
large F1 jump, and that is exactly what this shows: swapping two-body → n-body maps does **not**
degrade discrimination, and slightly improves it. The integrator bug's biggest effect (§1) was
along-track displacement of **fast, planet-crossing NEOs** far from opposition; the night-61642
tracklet set is dominated by moderate-elongation objects where two-body was already adequate, so the
map-level ΔF1 is small here. The tails (fast NEOs) and the full-2-yr set are where a larger separation
could appear — the SECONDARY readout (§12.4.1 step 3) tests that.

**Milestone B (case1) primary is complete.** Remaining/optional: the full-2-yr stages 3–6 re-score
(secondary), then case2/case3 (deferred).

### §12.4.3 — MILESTONE B secondary result: full-2yr case1 re-score (2026-07-18)

Full-2yr case1 re-scored in V-band against the n-body grid (`rescore_vdp_Vband.py --maps-dir
prob_maps_grid_s3m_nbody`, job 37295357, 648,769 rows, 9 min, label self-check 100.00%,
`P_NEO_vdp_Vband` 100% non-null). New file `sorcha_comparison_case1_nbody_Vband.parquet` (baseline
`..._Vband.parquet` preserved). ROC via `stage1_roc_case1.py` (NEO=positive; MBA/Trojan/TNO/other
negative):

| set | classifier | AUC | best-F1 | complete | contam |
|---|---|---|---|---|---|
| **night-61642 (epoch-matched)** | **n-body** | **0.893** | **0.793** | 71.7% | 11.3% |
| | two-body | 0.882 | 0.788 | 69.6% | 9.4% |
| | digest2 | 0.930 | 0.814 | 72.4% | 7.0% |
| **full 2-yr (mixed epochs)** | n-body | 0.881 | 0.802 | 73.6% | 11.9% |
| | two-body | 0.882 | **0.816** | 75.0% | 10.6% |
| | digest2 | 0.943 | 0.849 | 77.8% | 6.6% |

**Interpretation (a hypothesis raised and then REJECTED by the data).** The full-2yr row shows n-body
with a *lower* best-F1 (0.802 vs 0.816). I hypothesized a grid-epoch confound (n-body grid at 61642 /
Aug-2027 vs two-body grid at 2026-01-01, so the full-2yr comparison partly measures which single epoch
better represents the 2-yr average). **The data rejects that:** AUC binned by |night − 61642| is
**flat** — Δ(AUC) = +0.002 / −0.001 / +0.001 / −0.003 / −0.001 across 0–15 … 300+ day bins. There is
no near-epoch advantage for n-body, so it is not a grid-epoch effect.

**What the numbers actually say:**
- **Discrimination is identical.** AUC ≈ 0.88 for both maps *everywhere* (full-2yr 0.881 vs 0.882;
  every time bin within ±0.003). The n-body and two-body maps rank objects equally well.
- **On the controlled, epoch-matched comparison (night-61642), n-body is genuinely better** on both
  AUC (+0.011) and F1 (+0.005). This is the clean test and it favors n-body.
- **The full-2yr best-F1 gap (−0.014) is an operating-point/curve-shape artifact**, not a
  discrimination loss: two maps with equal AUC can differ slightly in ROC shape near the F1 optimum.
  Each classifier's F1 uses its OWN best threshold (n-body 0.015, two-body 0.010).
- digest2 still leads on both sets (F1 0.849 / 0.814), as in every prior comparison.

**Bottom line for the regen.** The n-body maps **discriminate as well as the two-body maps and better
when epoch-matched** — the integrator swap does not degrade VDP, and the pipeline is now
self-consistent (n-body scoring against n-body maps), which was the advisor's requirement. The regen
is validated for case1; a **large** F1 jump was never the expected outcome (§12.4.2) because VDP
discrimination is driven by the antisun-relative velocity structure, which both integrators capture at
these elongations — the integrator bug's damage was to *individual fast-NEO positions* (§1), which
maps average over.

**Milestone B (case1) COMPLETE — primary + secondary.** case2/case3 deferred (grid is
case-independent; only the per-case V-band re-score remains for them).

---

## §13 — PROCESS RUNBOOK (how the n-body regen was actually done) + Arnor deliverables

End-to-end, reproducible. Everything runs from `conda_prep/bin/python`, `PYTHONPATH=codex_integrator:neomod/src`,
`--account=astro --partition=ckpt-all`.

### 13.1 The pipeline, step by step

| # | step | script / job | output | notes |
|---|---|---|---|---|
| 0 | **n-body cache** — propagate all 14,380,436 S3M objects from t_0 (2008) → epoch MJD 61642, once | `stage0_epoch_state.py` (array 0-1885) → `stage0_merge.py` | `outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet` | the entire n-body cost; light-time ON; validated §10.9 |
| 1 | **integrator into production `gmm.py`** — `propagate_elements_nbody`, `t_0`-branched `build_visible_subset_dataframe`, origin fix, single-prop magnitude | edits to `neomod/src/velocity_density_pipeline_gmm.py` | — | §11; bit-identical to adam; lazy n-body imports |
| 2 | **cache branch** — `build_visible_subset_dataframe` skips propagation when the input already carries `ra_deg/dec_deg/vlam/vbeta` | same file (`_center_skycoord` helper + early branch) | — | §12.2; makes the grid regen cheap |
| 3 | **cache hook in the driver** — `--cache` builds `clone_sources` from the cache (shared observer-only scorer, `_mag_app` from cache) | `neomod/pipeline/sorcha_gen_maps_grid.py` | — | §12.4.1; antisun-relative centers are the driver's own (correct) |
| 4 | **667-grid build** — one center per task, `--ref-obstime 2027-08-25T00:00:00 --cache …` | `slurm/stage1_build_grid667_nbody.sbatch` (array 0-666) | `prob_maps_grid_s3m_nbody/` (667 npz) | §12.4.2; 667/667, ~6 min/task |
| 5a | **primary re-score** — night-61642, epoch-matched | `stage1_validate_milestoneA.py` | printout | §12.4.2; n-body ≥ two-body |
| 5b | **secondary re-score** — full-2yr V-band | `slurm/stage1_rescore_case1_nbody.sbatch` (`rescore_vdp_Vband.py --maps-dir prob_maps_grid_s3m_nbody`) | `sorcha_comparison_case1_nbody_Vband.parquet` | §12.4.3; 9 min |
| 6 | **ROC** — n-body vs two-body vs digest2, full-2yr + night-61642 | `stage1_roc_case1.py` | printout | §12.4.3 |
| 7 | **Arnor deliverables** — join predictions (cache) ⋈ Sorcha obs ⋈ scores | `make_arnor_deliverables.py` | `arnor_case1_nbody_{full,night61642}.parquet` | §13.2 |

**Gotchas hit and fixed along the way** (so they aren't re-hit): map-gen must run on a node
(login-node thread blowup); the map grid is **antisun-relative** — build at `antisun(ref_obstime)+dlon`,
not the stored 2026-01-01 absolute center (§12.3.1, cost us the first Milestone-A build); `score-vdp`
reuses the **stored** `prob_map_file` (antisun-relative) so re-scoring is just a `--maps-dir` swap;
the ROC F1 sweep must be vectorized (O(n log n), not O(unique·n)).

### 13.2 Arnor deliverables (`outputs/s3m_linking/case1/`)

Two parquets, built by `make_arnor_deliverables.py`. Each row = one case1 tracklet with **both**
Sorcha's observed values (`obs_*`) and our n-body prediction for that object (`pred_*`, the Stage 0
cache at MJD 61642 TDB), plus all three scores.

| file | rows | contents |
|---|---|---|
| `arnor_case1_nbody_full.parquet` | 648,769 (486,011 objects) | full 2-yr case1 subsample |
| `arnor_case1_nbody_night61642.parquet` | 4,316 (820 NEO) | night-61642 subset — **epoch-matched, the clean position test** |

Columns: `obs_ra_deg/obs_dec_deg/obs_ecl_lon/obs_ecl_lat/obs_vlam/obs_vbeta/obs_mag_V/obs_mjd_utc`
(Sorcha) · `pred_ra_deg/pred_dec_deg/pred_ecl_lon/pred_ecl_lat/pred_vlam/pred_vbeta/pred_mag_V`
(our n-body @ 61642) + `pred_a/e/i`, `pred_delta_au`, `pred_light_time_s` · `pred_epoch_mjd_tdb`
(=61642) · `dt_days` (=`obs_mjd_utc − 61642`) · `P_NEO_nbody`, `P_NEO_twobody`, `P_NEO_d2` · `population`.

Self-check on the epoch-matched night (reproduces acceptance #3): median sky-sep 120.9″ (RAW — see
below), median |Δv| **0.0078** deg/day, median |Δmag_V| **0.039**.

### 13.3 Is anything amiss for the comparison? — ONE caveat, then all-clear

**The one thing to flag: the epoch gap.** Predictions are at MJD **61642.0 TDB**; each tracklet was
observed at `obs_mjd_utc`. `dt_days` gives the gap. The object **really moves** in that gap, so the
**raw** sky-sep is dominated by orbital motion, NOT integrator error:
- **night-61642**: `|dt_days|` ≤ ~1 → raw sep ≈ 0.03° (121″) is intra-night motion. **Motion-correct**
  (subtract `pred_vlam/pred_vbeta × dt_days`, as acceptance #3 did) → residual collapses to the
  sub-arcsecond integrator floor (§10.9: 0.8″ at Δt<0.05 d, curvature-limited). Do the comparison here
  for positions.
- **full 2-yr**: `|dt_days|` up to ±365 → raw position comparison is meaningless without motion
  correction; use this file for the **score/ROC** comparison (`P_NEO_*`) and velocity (`|Δv|`, which is
  epoch-robust), not raw positions.

**Everything else is clean and matched:**
- **Velocity** `|Δv|` is epoch-robust (0.0078 deg/day, matches Sorcha) — valid on both files.
- **Magnitude**: `pred_mag_V` (S3M H_V + HG, G=0.15) and `obs_mag_V` (Sorcha LSST PSFMag → V) are both
  Johnson V; |Δmag| 0.039. Light-time is applied on both sides.
- **Convention**: `vlam = dλ/dt`, `vbeta = dβ/dt`, deg/day (NOT ×cosβ) — the Horizons-validated
  convention.
- **0 rows unmatched** to the cache; identity join is exact (same ObjID, same S3M source record).
- The full file is the ROC **subsample** (all NEOs + proportionally sampled non-NEOs), not the complete
  population — correct for ROC, not for absolute population counts.

**No blockers.** Do positions on night-61642 (motion-corrected); do scores/velocity on either.

### §13.4 — Arnor's ~20″ position offset: diagnosed (NOT aberration, NOT integrator) and FIXED

Arnor (`NOTE_FOR_HYAK_arnor_deliverable_position_offset.md`) found that `pred_ecl_lon/lat` vs
`obs_ecl_lon/lat` disagree by ~15–29″ **even at zero time gap** (|dt|<0.02 d, 63 rows), systematic in
direction (dlon −18.7″, dlat +15.6″), not scaling with light-time — and hypothesised **aberration**
(κ≈20.5″). Investigated end-to-end; the hypothesis is **wrong**, and the cause is a Sorcha *column*
artifact, not physics.

**The discriminating test — equatorial vs ecliptic, motion-corrected, binned by |dt|:**

| \|dt\| (d) | EQ (ra/dec) | ECL (Sorcha footprint lon/lat) |
|---|---|---|
| 0.00–0.02 | **0.97″** | 17.95″ |
| 0.02–0.05 | **0.78″** | 11.17″ |

The **equatorial** position agrees to **sub-arcsecond** — reproducing acceptance #3's 0.80″. The ~18″
lives **only** in the ecliptic columns.

**Aberration ruled out quantitatively.** Earth's heliocentric speed at 61642 is 29.46 km/s → v/c =
20.3″. Predicted aberration shift on the 63 near-zero-dt rows: median 14.4″. But aberration displaces
the real apparent direction, so it must appear in **both** frames — and the measured **equatorial**
residual is 1.0″, not 14″. Aberration is excluded.

**Root cause — Sorcha's `ecl_lon`/`ecl_lat` is an approximate footprint column.** Self-consistency
checks:
- cache `lam_deg` vs `_ecl_lon_from_radec(cache ra/dec)`: **0.00″** — the cache ecliptic is exact.
- Sorcha `ecl_lon` vs `_ecl_lon_from_radec(Sorcha mean_ra/dec)`: a **constant 13.858″** across ALL 4,316
  night-61642 rows (median 13.858, max 13.9 — a fixed systematic offset, not motion).

`sorcha_postprocess._ecl_lon_from_radec` is documented "sufficient for a footprint proximity test" and
is used only to assign tracklets to 30° map cells, where a constant ~14″ is negligible. It is **not** a
precision position, and it carries a fixed ~14″ offset from the exact transform of its own ra/dec. So
comparing the exact cache ecliptic to this approximate Sorcha ecliptic manufactured the ~18″.

**The fix (deliverables regenerated).** `make_arnor_deliverables.py` now computes `obs_ecl_lon/lat`
**exactly** from `obs_ra/dec` with the same 23.439291° obliquity the cache uses (the original Sorcha
column is retained as `sorcha_ecl_lon_footprint`), and adds `pred_dra_deg_day/pred_ddec_deg_day` for
equatorial motion correction. Self-check, |dt|<0.02 d, motion-corrected: **EQ 0.97″, ECL(exact) 0.97″**
(was ~18″), `|dv|` 0.0078, `|dmag_V|` 0.039. Both frames now agree to the sub-arcsecond integrator
floor.

**Conclusion.** The integrator and the absolute position are correct (0.8–1.0″ vs Sorcha); there is no
aberration issue and nothing to fix in the propagation. The apparent ~20″ was an artifact of comparing
against Sorcha's approximate footprint `ecl_lon`. §12 results stand (they never used that column). The
corrected deliverables make the ecliptic position comparison clean; equatorial ra/dec was always exact.

### §13.5 — CORRECTION to §13.4 (Arnor follow-up): the footprint offset has TWO components, PM-dominated

Arnor (`NOTE_FOR_HYAK_footprint_offset_correlates_with_pm.md`) rightly pushed back on §13.4's "constant
13.86″" characterization. Reconciled — we were each measuring a *different* difference, and both are
real:

| comparison | median | max | corr(\|pm\|) |
|---|---:|---:|---:|
| footprint vs `ecl(mean_ra/mean_dec)` **(§13.4's number)** | 13.86″ | 13.9″ | −0.25 |
| footprint vs `ecl(ra0/dec0)` = `obs_ecl` **(Arnor's `foff`)** | 20.57″ | 237.2″ | **+0.886** |

So the stored `ecl_lon/ecl_lat` footprint column = `ecl(mean_ra,mean_dec)` **plus a ~13.86″
near-constant** transform offset (its provenance differs from the exact fixed-obliquity formula — an
older/different transform in whatever produced that column; it is NOT `_ecl_lon_from_radec(mean)` to
<0.1″ on any of the 4,316 rows). **On top of that**, comparing it to `obs_ecl` adds a
**proper-motion-scaled ½-arc term**, because the footprint is referenced to the tracklet **mean** while
`obs_ecl`/`obs_ra_deg` are referenced to the **first detection `ra0`** (`corr(foff, ½·pm·span)=+0.908`).
That PM term dominates and is why §13.4's "constant" framing was wrong: TNOs (slow) sit near 14″, fast
NEOs reach 237″.

**Arnor's "detection-epoch mismatch" is the correct root cause of the dominant term**, and it's the
more useful framing for the docs. My §13.4 constant is a real but secondary component.

**None of this changes the conclusions:**
- The corrected deliverable is unaffected — `obs_ecl_lon/lat` is the exact transform of `obs_ra/dec`
  (Arnor reproduced: 0.0000″ median), and the position agrees with Sorcha to 0.97″ (both frames,
  motion-corrected). Still no aberration, still no integrator error.
- **Grid assignment is safe even for the fastest NEOs:** worst-case footprint error 237.2″ = 0.066°,
  vs a 5° half-cell (10° lon-step) and 30° cell radius — 1.3% of a half-cell, so a mis-assignment would
  require a tracklet to straddle a cell boundary within 0.066°. Negligible, and adjacent cells carry
  near-identical velocity structure anyway.
- **`obs_ecl` uses `ra0` deliberately**: `obs_mjd_utc = mjd0_utc` (first-detection time), so `obs`
  (position at `ra0`, time `mjd0_utc`) is the self-consistent reference for motion-correcting the cache
  (`dt_days = mjd0_utc − 61642`). Using the mean would desync position and time.

**One-line cleanup if we ever want the footprint column exact** (Arnor's suggestion): have
`sorcha_postprocess` feed the same mean-position input, or just derive `ecl_lon/lat` from `ra0/dec0`
consistently — but it only drives 30°-cell footprint assignment, so it is not worth a repipeline.

### §13.6 — The BENCHMARK regen: generator was doubly broken; rebuilt from the cache

Context: "benchmark" = our pipeline's OWN synthetic tracklets made from pure S3M (propagated with our
integrator), compared object-by-object against Sorcha's tracklets. NOT Sorcha scored on VDP. When the
n-body benchmark was regenerated by re-running `gen_benchmark_night61642.py`, it did **not** match
Sorcha (raw sky-sep still 1.5°, |dv| 0.011, 15% of objects >5° — i.e. two-body-like). Two independent
bugs, both found:

1. **Wrong elements (hybrid loader).** `load_s3m_population` defaults to `VDP_LOADER=hybrid`, loading
   `hybrid_elements.parquet` with `t_0 = 60065` (2023). Those osculating elements propagate (n-body)
   **0.47° median / 77° max** away from the validated pure-S3M cache (`t_0 = 54466`, 2008). The hybrid
   Keplerian conversion is not consistent with the S3M orbit the cache and Sorcha agree on.
2. **`adam_core_stub` on `sys.path`.** The benchmark scripts `sys.path.insert(0, adam_core_stub)` — a
   legacy shim so `vdp` imports without the heavy real `adam_core`. But the n-body path now needs the
   **real** `adam_core`/`assist`/`sorcha`; with the stub shadowing it, `score_orbital_df` silently ran
   two-body (6 s vs 16 s on 100 objects; benchmark MBA task 2 min vs the ~20 min n-body would take).

**Fix: build the benchmark directly from the Stage 0 CACHE** (`build_benchmark_night61642_from_cache.py`),
the validated n-body state (pure S3M, matches Sorcha to 0.03° in acceptance #3), reformatted into the
benchmark tracklet schema (antisun-relative `prob_map` assignment; `ra1 = ra0 + dra·DT_DAYS`). Result
vs Sorcha (shared `s3m_objid`, night 61642):

| metric | two-body benchmark | **n-body benchmark (cache)** |
|---|---|---|
| \|Δv\| median (deg/day) | 0.0106 | **0.0078** |
| raw sky-sep median | 1.50° | **0.034°** (45× collapse) |
| fraction > 5° | 15.0% | **0.0%** |
| target coverage | 8,436 / 11,612 | **11,612 / 11,612** |

So the n-body benchmark now agrees with Sorcha, as expected. **The generator bugs (hybrid loader +
adam_core_stub) also mean any future `score_orbital_df`/benchmark run must (a) set `VDP_LOADER=s3m`
and (b) NOT put `adam_core_stub` on the path** — but building from the cache sidesteps both and reuses
the already-paid, validated propagation. The two-body benchmark is preserved at
`outputs/benchmark_night61642_twobody/` for A/B. Map generation was unaffected (it was fed the cache
via the §12.2 branch, never `load_s3m_population`'s hybrid path).

**Full benchmark (v3-scale, scored) is still TODO** — same cache-based approach, plus VDP scoring
against the n-body maps and digest2, replicating v3's population caps.

### §13.7 — Full n-body benchmark: cache-based, self-contained overnight job (submitted 2026-07-19)

The full v3-scale benchmark, rebuilt cache-based (§13.6 established the generator is broken for
n-body). Single self-contained Slurm job `codex_integrator/slurm/benchmark_full_nbody.sbatch`
(job 37307535, 32 cpu / 96G / 12h, ckpt-all) — no dependency chain, so it survives unattended:

1. `build_benchmark_tracklets_from_cache.py` — per-population **filtered** cache reads (never loads
   all 14.4M at once; the login-node run OOM'd, a node with 96G is fine), sun-exclusion + v3
   proportional caps (NEO 12,900 / MBA 650,000 / TNO 1,600 / Trojans 6,000), synthetic 2-detection
   tracklets → `outputs/benchmark_tracklets_s3m_nbody/tracklets_benchmark_nbody.parquet` (single
   combined file = v3 tracklet-dir layout, so score-vdp reads exactly one shard).
2. `sorcha_phase2.py score-vdp` against the **n-body** maps `prob_maps_grid_s3m_nbody` (benchmark
   `mean_mag` is already Johnson V from S3M H, matching the V-band maps — no colour re-score needed).
3. `run-digest2` (orbit-based) via GNU parallel, 32 workers × 1000-row chunks, failure-tolerant.
4. `combine` → `outputs/phase2_benchmark_s3m_nbody/benchmark_comparison_s3m_nbody.parquet`.

Verified before submit: night-61642 cache benchmark matched Sorcha (|dv| 0.0078, sky-sep 0.034°,
§13.6); all score-vdp/run-digest2/combine args exist; builder byte-compiles; tracklet schema == v3.

### §13.8 — Full n-body benchmark RESULT (2026-07-19): works, beats two-body, matches Sorcha

Job 37307535 COMPLETED (3:19): `benchmark_comparison_s3m_nbody.parquet`, **670,500 rows** (MBA 650k /
NEO 12.9k / Trojans 6k / TNO 1.6k), 100% VDP + digest2 scored.

**Read the metric on the SCORABLE subset.** The full-set VDP AUC looks broken (0.074) because **68% of
benchmark objects are fainter than mag 25** — outside the maps' 14–25 range, so unscorable (`P=0`).
That is a *selection* artifact, not a scoring bug, and it is shared with v3 (67% faint): the builder
samples the whole cache, but Sorcha only *detects* bright objects. On the scorable subset (mag in map
range):

| benchmark | scorable rows | VDP AUC | VDP F1 | complete | contam |
|---|---:|---:|---:|---:|---:|
| **n-body (new)** | 215,645 | **0.878** | **0.662** | 57% | 22% |
| two-body v3 (old) | 150,250 | 0.873 | 0.589 | 51% | 31% |

**The n-body benchmark beats two-body v3** (F1 +0.073, contamination 22% vs 31%) and its **AUC 0.878
matches Sorcha's n-body VDP** (0.881 full-2yr / 0.893 night-61642, §12.4) — benchmark and Sorcha now
agree, which is the point of the whole regen. In-range NEO median `P_vdp = 1.000` vs MBA 0.001 (clean
separation). digest2 (orbit-based) ~0.9 for both, as expected.

**One caveat / refinement.** The benchmark includes undetectable faint objects (mag > 25) that Sorcha
would never see, inflating row count and the full-set metric. For a cleaner benchmark-vs-Sorcha
comparison the builder should apply a detection cut (e.g. mag < 24.5) so its object set matches
Sorcha's detected set. The parquet carries `mean_mag`/`mag_bin_label`, so this is a one-line filter in
analysis, or a small builder change if we want the on-disk set to match Sorcha. Deferred (v3 had the
same, so the comparison above is apples-to-apples).

**Benchmark regen COMPLETE** (night-61642 §13.6 + full §13.7/§13.8), both cache-based and validated
against Sorcha. Deliverables: `benchmark_night61642/benchmark_night61642.parquet` (velocity/position,
identity-shared) and `phase2_benchmark_s3m_nbody/benchmark_comparison_s3m_nbody.parquet` (full, scored).

### §13.9 — Night-61642 benchmark SCORED (Arnor request): scoring-consistency test passes

Arnor needed the night-61642 benchmark with a VDP score (the two files sent each had only half: the
night set OR the score). Built `benchmark_night61642_nbody_scored.parquet` via
`score_night_benchmark_nbody.py` = `benchmark_night61642.parquet` scored against `prob_maps_grid_s3m_nbody`
through the SAME `rescore_vdp_Vband.rescore` path Sorcha's `P_NEO_nbody` used (support_mask_min=1, mask
OFF, Johnson V) — benchmark `mean_mag` is already V, no conversion. 11,612 rows, `P_NEO_vdp` **100%
non-null** (night set is the detected union, 100% mag<25), NEO median 1.000 / MBA 0.001.

**Scoring-consistency (identity join, 4,316 matched, Sorcha `P_NEO_nbody` vs benchmark `P_NEO_vdp`):**
**|ΔP| median = 0.0002, 97.6% within 0.1** (NEO |ΔP| 0.0000, MBA 0.0002). With the integrator and maps
now n-body on BOTH sides, benchmark and Sorcha score the same object identically — the same-object
scoring-consistency test (previously two-body-benchmark vs LSST-then-V-Sorcha) now closes to ~0.
