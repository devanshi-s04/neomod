# Fix plan: put Sorcha VDP scoring in a magnitude band consistent with the maps

**Date:** 2026-07-16
**Depends on:** `SORCHA_MAG_BAND_FINDING.md` (the ~0.6 mag V-vs-LSST-band mismatch).
**Goal:** score VDP with the tracklet magnitude in the **same band as the map mag bins**,
so a mag bin means the same thing on maps, benchmark, and Sorcha.

---

## First: "can we just tell Sorcha to output V?" — No.

Checked the Sorcha docs (outputs page) and the config. Sorcha reports apparent magnitude
**only in the observed LSST filter** per visit (`PSFMag`/`trailedSourceMag` in ugrizy);
there is **no config option** to emit a fixed/main-filter or Johnson-V magnitude. LSST
observes ugrizy, so there is no "V" exposure to report. Setting `observing_filters` or the
physical-file "main filter" changes which band `H` is defined in, but the per-visit output
is still whatever LSST filter was used. **So the band conversion must be done by us, at
scoring time — not inside Sorcha.** The good news: Sorcha already gives us everything
needed (`PSFMag`, `H_filter`, `optFilter`), and our phys file carries the per-object V→LSST
color.

---

## Two directions (pick based on the science goal)

### Option A — convert Sorcha mags to Johnson V (match the existing V maps)
Keep the maps as-is; fix the *scoring input*. In `sorcha_postprocess.py`, compute a
V-band apparent magnitude per detection and use it for `mean_mag` (hence the VDP mag bin):

```
reduced   = PSFMag − H_filter            # filter-independent geometry (distance + phase)
mag_V      = reduced + H_V               # H_V = Johnson-V absolute mag (S3M .s3m H)
mean_mag_V = 0.5 * (mag_V_det0 + mag_V_det1)
```

- `H_filter` is already in the Sorcha output. `H_V` we have per object: it's the `.s3m` `H`
  (also recoverable as `H_r + (Johnson V − LSST r)` from the phys-file color, seed 42).
- **Cost:** a postprocess column change + re-run **stage 3 (`3_vdp.sh`)** on the existing
  tracklets. No Sorcha re-run, no new production. Cheap.
- **Caveat (from the finding):** the reduced-mag conversion removes the ~0.6 mag *band*
  term but leaves ~0.2 mag of synthetic-vs-observed residual (phase `G`, trailing loss,
  photometric noise). For *internal consistency with the V maps* that residual doesn't
  matter — we only need Sorcha and the maps to be in the same system, and both then use V.
- **This is what makes the current benchmark comparison valid**, and is the minimal change.

### Option B — rebuild the maps in an LSST band (e.g. r) [more correct long-term]
Real Rubin data arrives in ugrizy, never V. A production NEO classifier will *always* see
LSST-band magnitudes. So arguably the maps should be built in an LSST band, and everything
(maps, Sorcha, real NEOCP) scored in that band — no V anywhere.

- Rebuild `prob_maps_grid_s3m` with apparent **r** (from `H_r`) instead of apparent V; then
  Sorcha's r-ish mags are consistent by construction, and the pipeline matches reality.
- **Cost:** full map regeneration (the GMM cloner + VDP visibility over the hybrid/S3M
  population, per mag bin) — the expensive path, a SLURM map-build like the original grid.
- Also requires the benchmark to be re-scored on r (trivial — it has `H` and we can derive
  r via the color), so the benchmark↔Sorcha comparison stays apples-to-apples.

---

## Recommendation

1. **Now (unblocks Arnor's test):** Option A. Add `mean_mag_V` to the postprocess (or as a
   standalone correction column keyed by `ObjID`), re-run stage 3 to produce a
   `P_NEO_vdp_Vband` for the case tracklets, and re-score the benchmark the same way (it's
   already V, so it's unchanged). Then the "same object → same score" test is clean on the
   band axis and any residual ΔP is real (velocity tail, the ~0.2 photometry residual, or a
   genuine map disagreement).
2. **For the paper / production framing:** decide Option B with the advisor. If the claim is
   "a Rubin-era classifier," the maps arguably *should* be in an LSST band, and V is only an
   intermediate. That's a bigger regen but the scientifically defensible end state. Option A
   and Option B give the same *relative* answer for the benchmark test; they differ in which
   band is treated as ground truth.

## Concrete next steps if we go with Option A
- [ ] In `sorcha_postprocess.py`, join `H_V` per `ObjID` (from the `.s3m` census or
      `H_r + colorVr`) and write `mean_mag_V` alongside `mean_mag` (keep both — don't
      overwrite, so nothing silently changes).
- [ ] Add a `--mag-col mean_mag_V` option to `score-vdp` (or a one-off scoring run) so the
      VDP mag-bin lookup uses the V mag; re-run stage 3 for case1/2/3 → `P_NEO_vdp_Vband`.
- [ ] Re-score the benchmark identically (unchanged band) for the comparison.
- [ ] Quantify the F1/AUC delta vs the old LSST-band scoring — this tells us how much the
      mis-binning actually cost, which is a paper-worthy number in its own right.
- [ ] Sanity: after the fix, re-run Arnor's mag-bin-crossing check — the ~60% cross rate
      should collapse toward the noise floor for same-velocity objects.

## What NOT to do
- Don't try to make Sorcha output V (impossible, per the docs).
- Don't overwrite `mean_mag` in place — add a new column, so existing results stay
  reproducible and the change is auditable.
- Don't assume the ~0.2 mag residual is a bug — it's synthetic-vs-observed photometry and is
  expected; only the ~0.6 mag band term is the defect.
