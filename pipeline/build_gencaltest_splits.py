#!/usr/bin/env python3
"""Freeze the GEN / CAL / TEST splits (EVALUATION_PROTOCOL.md v1.1 §0).

NEO  : independent NEOMOD3 draws (distinct seeds) -- handled by gen_benchmark_tracklets_neomod3.py
       GEN  = seed 42        (the existing 100M-draw projection cache that built the maps)
       CAL  = seed 20270825  (the existing 30M draw, used only for the E1a smoke test)
       TEST = fresh seed, sealed until the single final scoring
non-NEO: the un-cloned Stage-0 objects PARTITIONED BY PARENT ObjID, disjointly, 60/20/20.

Assignment is by a STABLE HASH of ObjID (md5 with a recorded salt), not by row order or an RNG
stream, so the split is reproducible, order-independent, and unchanged if rows are added later.

Writes a frozen manifest + provenance recording the exact WEIGHTED retained fraction per population
(protocol §0.2: f is measured per population, never hardcoded to 0.6) and per (population, mag bin)
so E0 can check the binomial scatter is negligible.
"""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
CACHE = W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"
OUT = W/"outputs/splits"; OUT.mkdir(parents=True, exist_ok=True)
SALT = "neomod3-gencaltest-v1.1"          # recorded; changing it changes the split
FRAC = {"GEN": 60, "CAL": 20, "TEST": 20}  # percent, cumulative below
MAG_BINS = [(14,16),(16,18),(18,20),(20,21),(21,22),(22,23),(23,24),(24,25)]
BIN_LABELS = ["14_16","16_18","18_20","mag20","mag21","mag22","mag23","mag24+"]


def assign(objids: np.ndarray) -> np.ndarray:
    """Stable per-object split label from md5(salt|ObjID) -> [0,100)."""
    out = np.empty(len(objids), dtype=object)
    for i, o in enumerate(objids):
        h = hashlib.md5(f"{SALT}|{o}".encode()).digest()
        v = int.from_bytes(h[:8], "big") % 100
        out[i] = "GEN" if v < 60 else ("CAL" if v < 80 else "TEST")
    return out


def main():
    src_hash = hashlib.md5(open(CACHE, "rb").read(1 << 20)).hexdigest()[:16]  # header hash (cheap id)
    df = pd.read_parquet(CACHE, columns=["ObjID", "population", "mag_app"])
    df = df[df.population != "NEO"].reset_index(drop=True)     # NEO comes from NEOMOD3 draws
    print(f"non-NEO objects: {len(df):,}")
    df["split"] = assign(df.ObjID.to_numpy())

    # --- zero-overlap assertion (protocol §0.1: no parent may cross splits) ---
    per_obj = df.groupby("ObjID", observed=True).split.nunique()
    assert per_obj.max() == 1, f"{(per_obj>1).sum()} ObjIDs appear in more than one split"
    sets = {k: set(df.ObjID[df.split == k]) for k in FRAC}
    for a in FRAC:
        for b in FRAC:
            if a < b:
                ov = len(sets[a] & sets[b])
                assert ov == 0, f"{a}/{b} overlap: {ov}"
    print("zero parent overlap: ASSERTED")

    # --- exact WEIGHTED retained fraction, per population (never hardcoded) ---
    df["w"] = 1.0                                   # Stage-0 rows are unweighted objects
    prov = {"salt": SALT, "source": str(CACHE), "source_hash16": src_hash,
            "target_pct": FRAC, "f_gen_by_population": {}, "f_by_population_magbin": {},
            "counts": {}}
    print(f"\n{'population':>9} {'GEN':>10} {'CAL':>10} {'TEST':>10} {'f_GEN(weighted)':>16}")
    for pop, g in df.groupby("population", observed=True):
        tot_w = g.w.sum()
        f_gen = g.w[g.split == "GEN"].sum()/tot_w
        prov["f_gen_by_population"][pop] = float(f_gen)
        prov["counts"][pop] = {k: int((g.split == k).sum()) for k in FRAC}
        print(f"{pop:>9} {prov['counts'][pop]['GEN']:>10,} {prov['counts'][pop]['CAL']:>10,} "
              f"{prov['counts'][pop]['TEST']:>10,} {f_gen:>16.6f}")
        per_bin = {}
        for (lo, hi), lab in zip(MAG_BINS, BIN_LABELS):
            m = (g.mag_app >= lo) & (g.mag_app < hi)
            if m.sum():
                per_bin[lab] = float(g.w[m & (g.split == "GEN")].sum()/g.w[m].sum())
        prov["f_by_population_magbin"][pop] = per_bin

    print("\nper-magbin f_GEN (binomial scatter around the global value):")
    fb = pd.DataFrame(prov["f_by_population_magbin"]).T
    print(fb.to_string(float_format=lambda v: f"{v:.4f}"))
    print("  max |deviation| from the population-level f_GEN: "
          f"{max(abs(fb[c][p]-prov['f_gen_by_population'][p]) for p in fb.index for c in fb.columns if pd.notna(fb[c][p])):.4f}")

    df[["ObjID", "population", "split"]].to_parquet(OUT/"nonneo_split_manifest.parquet", index=False)
    json.dump(prov, open(OUT/"split_provenance.json", "w"), indent=2)
    print(f"\nwrote {OUT/'nonneo_split_manifest.parquet'} and split_provenance.json")
    print("NEO splits (independent NEOMOD3 draws, no partition needed):")
    print("  GEN  = seed 42        (existing 100M-draw projection cache -- built the maps)")
    print("  CAL  = seed 20270825  (existing 30M draw -- E1a smoke test only)")
    print("  TEST = fresh seed     (sealed; drawn separately)")


if __name__ == "__main__":
    main()
