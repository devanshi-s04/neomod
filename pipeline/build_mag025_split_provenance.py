#!/usr/bin/env python3
"""Exact GEN retained fractions at 0.25-mag apparent-V resolution.

WHY THIS EXISTS
---------------
`velocity_density_pipeline_neomod_clone_only._nonneo_split_fraction` resolves the GEN retained
fraction f from `_SPLIT_BIN_LABELS`, which contains only the eight 1-mag bins. Any other
(mag_min, mag_max) raises ValueError -- deliberately, because a map built from a GEN split holding
fraction f of the real objects carries a density f x too low, and silently falling back to f=1.0
would inflate P(NEO) uniformly everywhere in a way no ROC curve reveals.

So 0.25-mag maps need 0.25-mag fractions. This computes them EXACTLY from the same two artifacts the
1-mag values came from -- the split manifest and the epoch-state apparent-V magnitudes -- rather
than interpolating or reusing the enclosing 1-mag value.

    f(pop, bin) = N_GEN(pop, bin) / N_all(pop, bin)

Cells with no objects get f = None and are recorded as EMPTY; the builder treats them as INVALID
under the insufficient-support rule. They are never assigned f = 1.0.

Writes outputs/splits/split_provenance_mag025.json. Read-only w.r.t. everything else.
"""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
OUT = W / "outputs" / "splits" / "split_provenance_mag025.json"
NONNEO_POPS = ("MBA", "TNO", "Trojans")
MAG_LO, MAG_HI, MAG_STEP = 14.0, 25.0, 0.25


def mag025_bins():
    """44 half-open bins [lo, lo+0.25) covering 14.00 <= V < 25.00."""
    n = int(round((MAG_HI - MAG_LO) / MAG_STEP))
    out = []
    for i in range(n):
        lo = MAG_LO + i * MAG_STEP
        hi = lo + MAG_STEP
        out.append({"label": f"V{lo:06.2f}_{hi:06.2f}", "mag_min": round(lo, 2),
                    "mag_max": round(hi, 2), "index": i})
    return out


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    cache_p = W / "outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"
    man_p = W / "outputs/splits/nonneo_split_manifest.parquet"
    cache = pd.read_parquet(cache_p, columns=["ObjID", "population", "mag_app"])
    man = pd.read_parquet(man_p, columns=["ObjID", "split"])
    print(f"epoch cache {len(cache):,} rows | manifest {len(man):,} rows", flush=True)

    gen = set(man.ObjID[man.split == "GEN"])
    bins = mag025_bins()
    f_by = {p: {} for p in NONNEO_POPS}
    counts = {p: {} for p in NONNEO_POPS}
    empty = []
    for pop in NONNEO_POPS:
        sub = cache[cache.population == pop]
        v = sub.mag_app.to_numpy(float)
        isgen = sub.ObjID.isin(gen).to_numpy()
        for b in bins:
            m = np.isfinite(v) & (v >= b["mag_min"]) & (v < b["mag_max"])
            n_all = int(m.sum())
            n_gen = int((m & isgen).sum())
            counts[pop][b["label"]] = {"n_all": n_all, "n_gen": n_gen}
            if n_all == 0:
                f_by[pop][b["label"]] = None
                empty.append((pop, b["label"]))
                continue
            if n_gen == 0:
                # objects exist but none were retained by GEN: f=0 would divide by zero downstream
                f_by[pop][b["label"]] = None
                empty.append((pop, b["label"] + " (n_gen=0)"))
                continue
            f_by[pop][b["label"]] = n_gen / n_all

    # cross-check: the 0.25-mag fractions must aggregate to the frozen 1-mag values
    old = json.load(open(W / "outputs/splits/split_provenance.json"))
    check = {}
    for pop in NONNEO_POPS:
        sub = cache[cache.population == pop]
        v = sub.mag_app.to_numpy(float)
        isgen = sub.ObjID.isin(gen).to_numpy()
        for lab, (lo, hi) in {"14_16": (14, 16), "16_18": (16, 18), "18_20": (18, 20),
                              "mag20": (20, 21), "mag21": (21, 22), "mag22": (22, 23),
                              "mag23": (23, 24), "mag24+": (24, 25)}.items():
            m = np.isfinite(v) & (v >= lo) & (v < hi)
            if not m.sum():
                continue
            f_new = float((m & isgen).sum() / m.sum())
            f_old = old["f_by_population_magbin"].get(pop, {}).get(lab)
            if f_old is not None:
                check[f"{pop}/{lab}"] = {"recomputed": f_new, "frozen": float(f_old),
                                         "abs_diff": abs(f_new - float(f_old))}
    worst = max((c["abs_diff"] for c in check.values()), default=0.0)
    print(f"\n1-mag reproduction check: {len(check)} bins, max |recomputed - frozen| = {worst:.3e}",
          flush=True)
    for k, c in sorted(check.items(), key=lambda kv: -kv[1]["abs_diff"])[:5]:
        print(f"    {k:16s} recomputed {c['recomputed']:.6f}  frozen {c['frozen']:.6f}  "
              f"diff {c['abs_diff']:.3e}", flush=True)

    n_ok = sum(1 for p in NONNEO_POPS for b in bins if f_by[p][b["label"]] is not None)
    print(f"\n0.25-mag cells: {len(bins)*len(NONNEO_POPS)} total, {n_ok} with a valid fraction, "
          f"{len(empty)} empty/unretained", flush=True)

    out = {
        "salt": old.get("salt"), "source": str(man_p), "source_hash16": sha256_file(man_p)[:16],
        "epoch_cache": str(cache_p), "epoch_cache_sha256_16": sha256_file(cache_p)[:16],
        "magnitude_quantity": "apparent V (HG, G=0.15)",
        "bin_scheme": {"lo": MAG_LO, "hi": MAG_HI, "step": MAG_STEP, "n_bins": len(bins),
                       "semantics": "half-open [lo, lo+step)"},
        "bins": bins,
        "f_by_population_magbin": f_by,
        "counts_by_population_magbin": counts,
        "empty_or_unretained_cells": [f"{p}/{b}" for p, b in empty],
        "one_mag_reproduction_check": check,
        "one_mag_reproduction_max_abs_diff": worst,
        "note": ("f is EXACT per 0.25-mag bin, not interpolated from the 1-mag values. Cells with "
                 "no objects or no GEN objects carry f=None and are INVALID for map building; they "
                 "are never assigned f=1.0."),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
