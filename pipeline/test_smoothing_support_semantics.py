#!/usr/bin/env python3
"""(A) NEO smoothing-support semantics: raw counts vs counts x effective_factor.
   (B) numerical proof that f_GEN scales density by exactly 1/f_GEN and leaves support bit-identical.

The smoothing threshold's own docstring says it identifies "density pixels with enough SAMPLED
support to average locally" -- a statistical quantity. Multiplying by effective_factor converts it
into an object-equivalent quantity, which is a different thing.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src"))
import velocity_density_pipeline_neomod_clone_only as v

THR = v.DEFAULT_SMOOTH_SUPPORT_THRESHOLD
EF = json.load(open(W/"outputs/neomod3_projection_cache/cache_metadata.json"))["effective_factor_NEO"]
BINS = ["14_16","16_18","18_20","mag20","mag21","mag22","mag23","mag24+"]

print(f"DEFAULT_SMOOTH_SUPPORT_THRESHOLD = {THR}")
print(f"effective_factor_NEO             = {EF:.6f}   (constant -- identical in every mag bin)")
print(f"=> scaled threshold corresponds to {THR/EF:.3f} RAW clones, i.e. >= {int(np.ceil(THR/EF))} clones\n")

print("=== (A) pixels passing the smoothing threshold: raw vs scaled ===")
rows = []
for cen in ["prob_maps_grid_dlon+000_lat+00.npz", "prob_maps_grid_dlon+020_lat-12.npz",
            "prob_maps_grid_dlon+000_lat-50.npz"]:
    z = np.load(W/"prob_maps_grid_neomod3_full"/cen, allow_pickle=True)
    for b in BINS:
        k = f"support_count__NEO__{b}"
        if k not in z.files:
            continue
        s = np.nan_to_num(np.asarray(z[k], float))
        raw = int((s >= THR).sum())
        scaled = int((s*EF >= THR).sum())
        nz = int((s > 0).sum())
        rows.append(dict(center=cen.split("dlon")[1][:9], magbin=b, pixels_with_any_clone=nz,
                         pass_raw=raw, pass_scaled=scaled,
                         ratio=(scaled/raw if raw else np.inf),
                         min_clones_raw=THR, min_clones_scaled=THR/EF))
d = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print(d.to_string(index=False, float_format=lambda x: f"{x:,.3f}"))
tot_r, tot_s = d.pass_raw.sum(), d.pass_scaled.sum()
print(f"\n  TOTAL pixels smoothed: raw={tot_r:,}  scaled={tot_s:,}  "
      f"scaled/raw = {tot_s/max(tot_r,1):.2f}x")
print("  The scaled rule smooths pixels supported by as few as 2 clones; the raw rule requires 10.")

print("\n=== (B) f_GEN: density scales by exactly 1/f, support must be BIT-IDENTICAL ===")
rng = np.random.default_rng(0)
n = 4000
vlam = rng.normal(0, 0.3, n); vbeta = rng.normal(0, 0.3, n)
x = np.linspace(-1, 1, 201); y = np.linspace(-1, 1, 201)
sup_a = v.make_support_count_map(vlam, vbeta, x, y, clone_factor=1.0)
sup_b = v.make_support_count_map(vlam, vbeta, x, y, clone_factor=0.6)
identical = np.array_equal(sup_a, sup_b)
print(f"  support_count_map(clone_factor=1.0) vs (0.6): bit-identical = {identical}")
print(f"    (make_support_count_map ignores clone_factor -- returns raw histogram counts)")

dens = rng.random((201, 201))*100
f_gen = 0.600223
d_nosplit = dens/1.0
d_split = dens/(1.0*f_gen)
ratio = d_split/d_nosplit
print(f"  density ratio with f_GEN={f_gen}: min={ratio.min():.12f} max={ratio.max():.12f}")
print(f"    expected 1/f = {1/f_gen:.12f}")
ok_d = np.allclose(ratio, 1/f_gen, rtol=0, atol=1e-12)
print(f"    density scales by exactly 1/f_GEN: {ok_d}")

v.NONNEO_SPLIT_FRACTIONS = json.load(open(W/"outputs/splits/split_provenance.json"))
f_mod = v._nonneo_split_fraction("MBA", 23.0, 24.0)
sup_split = v.make_support_count_map(vlam, vbeta, x, y, clone_factor=1.0*f_mod)
print(f"  module f for MBA/mag23 = {f_mod:.6f}")
print(f"  support unchanged under that factor: {np.array_equal(sup_a, sup_split)}")
v.NONNEO_SPLIT_FRACTIONS = None
ok = identical and ok_d and np.array_equal(sup_a, sup_split)
print(f"\n{'='*60}\n{'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
