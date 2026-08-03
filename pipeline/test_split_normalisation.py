#!/usr/bin/env python3
"""Unit-test the frozen split before any map is built (protocol v1.1 §0).

1. manifest: zero parent overlap, every object assigned exactly once
2. every (population, mag bin) fraction in the provenance recomputed from the manifest + cache
3. the module resolves the SAME fraction it will use at map-build time
4. a missing/invalid fraction is a HARD ERROR, not a silent 1.0
5. density and support factors are separated
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src"))
import velocity_density_pipeline_neomod_clone_only as v

man = pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
prov = json.load(open(W/"outputs/splits/split_provenance.json"))
cache = pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet",
                        columns=["ObjID", "population", "mag_app"])
cache = cache[cache.population != "NEO"]
BINS = [(14.,16.,"14_16"),(16.,18.,"16_18"),(18.,20.,"18_20"),(20.,21.,"mag20"),
        (21.,22.,"mag21"),(22.,23.,"mag22"),(23.,24.,"mag23"),(24.,25.,"mag24+")]
fail = 0

print("=== 1. manifest integrity ===")
assert man.ObjID.duplicated().sum() == 0, "duplicate ObjID in manifest"
n_multi = man.groupby("ObjID", observed=True).split.nunique().gt(1).sum()
print(f"  objects: {len(man):,}   duplicated: 0   in >1 split: {n_multi}")
assert n_multi == 0
print(f"  split sizes: {man.split.value_counts().to_dict()}")

print("\n=== 2. recompute every fraction from the manifest ===")
j = cache.merge(man[["ObjID", "split"]], on="ObjID", how="inner")
print(f"  joined {len(j):,} of {len(cache):,} cache rows")
for pop, g in j.groupby("population", observed=True):
    for lo, hi, lab in BINS:
        m = (g.mag_app >= lo) & (g.mag_app < hi)
        if not m.sum():
            continue
        recomputed = (g.split[m] == "GEN").sum()/m.sum()
        stored = prov["f_by_population_magbin"].get(pop, {}).get(lab)
        ok = stored is not None and abs(recomputed - stored) < 1e-9
        if not ok:
            fail += 1
            print(f"  MISMATCH {pop}/{lab}: stored={stored} recomputed={recomputed:.6f}")
print(f"  fractions verified: {'ALL MATCH' if fail == 0 else f'{fail} MISMATCHES'}")

print("\n=== 3. module resolves the same fraction ===")
v.NONNEO_SPLIT_FRACTIONS = prov
for pop in ["MBA", "TNO", "Trojans"]:
    for lo, hi, lab in BINS:
        stored = prov["f_by_population_magbin"].get(pop, {}).get(lab)
        if stored is None:
            continue
        got = v._nonneo_split_fraction(pop, lo, hi)
        if abs(got - stored) > 1e-12:
            fail += 1; print(f"  MISMATCH module {pop}/{lab}: {got} vs {stored}")
print("  module lookup: OK" if fail == 0 else "  module lookup: FAILED")
assert v._nonneo_split_fraction("NEO", 23., 24.) == 1.0, "NEO must never be split-corrected"
print("  NEO exempt: OK")

print("\n=== 4. invalid provenance must RAISE, not silently return 1.0 ===")
for bad, why in [({"f_by_population_magbin": {}}, "missing population"),
                 ({"f_by_population_magbin": {"MBA": {}}}, "missing bin"),
                 ({"f_by_population_magbin": {"MBA": {"mag23": 0.0}}}, "zero fraction"),
                 ({"f_by_population_magbin": {"MBA": {"mag23": 1.7}}}, "fraction > 1")]:
    v.NONNEO_SPLIT_FRACTIONS = bad
    try:
        v._nonneo_split_fraction("MBA", 23., 24.); fail += 1
        print(f"  NO RAISE for {why}  <-- BUG")
    except ValueError:
        print(f"  raises on {why}: OK")
v.NONNEO_SPLIT_FRACTIONS = None
assert v._nonneo_split_fraction("MBA", 23., 24.) == 1.0, "no-split default must be 1.0"
print("  no-split default 1.0: OK")

print("\n=== 5. density vs support factor separation ===")
src = (W/"neomod/src/velocity_density_pipeline_neomod_clone_only.py").read_text()
assert "clone_factor=support_factor" in src, "support map still uses effective_factor"
assert "support_count_map * support_factor" in src, "smoothing still scales by effective_factor"
assert "effective_factor = f * _nonneo_split_fraction" in src, "density factor lost the split fraction"
print("  support uses support_factor (no split), density uses effective_factor (with split): OK")

print(f"\n{'='*60}\n{'ALL CHECKS PASS' if fail == 0 else f'{fail} FAILURES'}")
sys.exit(1 if fail else 0)
