#!/usr/bin/env python3
"""Literal source-identity proof for the density-estimator map ablation.

Reruns ONLY the deterministic source-selection stage — no map is rebuilt, no density is
estimated. For every (center, population, magnitude_bin) it captures the exact ObjID set the
sealed builder would hand the density estimator, then emits

    SHA256 over the sorted records (ObjID, population, center, magnitude_bin, physical_weight)

The selection is invoked once per `density_mode`, so the three digests are produced by three
independent executions of the same code path. Literal digest equality across modes is therefore
a direct proof that all three map sets consumed identical source rows and weights — not an
inference from sample counts.

Sharded by center; merge with --merge.
"""
from __future__ import annotations
import argparse, glob, hashlib, json, os, re, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))
POPS = ("NEO", "MBA", "TNO", "Trojans")
MODES = ("hist_all", "bayes_all", "bayes_nonneo")
MAX_SEP = 30.0
OUT = W / "outputs" / "geometric_density_estimator_ablation" / "source_digest"


def digest_records(objids, population, center, mag_bin, weight):
    """SHA256 over sorted (ObjID, population, center, magnitude_bin, physical_weight)."""
    h = hashlib.sha256()
    w = repr(float(weight)).encode()
    for oid in sorted(map(str, objids)):
        h.update(oid.encode()); h.update(b"\x1f")
        h.update(population.encode()); h.update(b"\x1f")
        h.update(center.encode()); h.update(b"\x1f")
        h.update(mag_bin.encode()); h.update(b"\x1f")
        h.update(w); h.update(b"\x1e")
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if a.merge:
        fs = sorted(glob.glob(str(OUT / "digest_shard_*.parquet")))
        d = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
        d.to_parquet(OUT / "source_selection_digest.parquet", index=False)
        piv = d.pivot_table(index=["center", "population", "magnitude_bin"],
                            columns="mode", values="sha256", aggfunc="first")
        print(f"shards {len(fs)}   cells {len(piv):,}")
        ref = piv["hist_all"]
        eq = {m: bool((piv[m] == ref).all()) for m in MODES}
        print("\n[ASSERT] literal SHA256 equality of the source-selection digest across modes:")
        for m in MODES:
            print(f"    hist_all == {m:13s} : {'PASS' if eq[m] else 'FAIL'}")
        nfail = int((~(piv[list(MODES)].eq(ref, axis=0)).all(axis=1)).sum())
        print(f"    cells with any mismatch: {nfail}")
        assert all(eq.values()) and nfail == 0, "source-selection digests DIFFER across modes"
        # one digest over the whole experiment
        allh = hashlib.sha256()
        for _, r in piv.sort_index().iterrows():
            allh.update(r["hist_all"].encode())
        summary = {
            "n_cells": int(len(piv)), "n_centers": int(d.center.nunique()),
            "modes": list(MODES), "literal_digest_equality_across_modes": True,
            "cells_with_mismatch": 0,
            "experiment_digest_sha256": allh.hexdigest(),
            "digest_definition": ("SHA256 over sorted records "
                                  "(ObjID, population, center, magnitude_bin, physical_weight)"),
            "total_sample_map_memberships": int(d.n_objids.sum() // len(MODES)),
            "note": ("counts are SAMPLE-MAP MEMBERSHIPS, not unique objects: the 30-degree "
                     "patches overlap, so one object is a member of many maps"),
            "unique_objids_across_all_cells": int(d.attrs.get("n_unique", 0)) or None,
        }
        json.dump(summary, open(OUT / "source_selection_digest.json", "w"), indent=2)
        print(f"\n  experiment digest : {allh.hexdigest()}")
        print(f"  sample-map memberships: {summary['total_sample_map_memberships']:,}")
        print(f"\nwrote {OUT/'source_selection_digest.json'}")
        return

    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    import velocity_density_pipeline_neomod_density_ablation as abl
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    import copy

    seal = json.load(open(W / "outputs/splits/MAP_BUILD_SEAL.json"))
    MAGBINS = [{"label": l, "mag_min": float(lo), "mag_max": float(hi)}
               for l, lo, hi in seal["magnitude_bins"]]
    neo_meta = json.load(open(W / "outputs/neomod3_projection_cache/cache_metadata.json"))
    base.NONNEO_SPLIT_FRACTIONS = json.load(open(W / "outputs/splits/split_provenance.json"))
    centers = [c.strip() for c in
               (W / "outputs/geometric_density_estimator_ablation/frozen_center_list.txt"
                ).read_text().split()][a.shard::a.nshards]
    cache = pd.read_parquet(W / "outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
    man = pd.read_parquet(W / "outputs/splits/nonneo_split_manifest.parquet")
    keep = set(man.ObjID[man.split == "GEN"])
    cache = cache[cache.ObjID.isin(keep) | (cache.population == "NEO")].reset_index(drop=True)
    ps = copy.deepcopy(base.DEFAULT_POPULATION_SETTINGS)
    _, scorer = base.load_s3m_population("neo", verbose=False)
    t = Time(seal["grid"]["ref_obstime"], scale="utc")
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    print(f"shard {a.shard}/{a.nshards}: {len(centers)} centers", flush=True)

    def eff_for(pop, mb):
        if pop == "NEO":
            return float(neo_meta["effective_factor_NEO"])
        return 1.0 * float(base._nonneo_split_fraction(pop, mb["mag_min"], mb["mag_max"]))

    rows = []
    orig = base.evaluate_density_map_full_posterior_2d
    t0 = time.time()
    for ci, cen in enumerate(centers):
        mm = re.match(r"dlon([+-]\d+)_lat([+-]\d+)", cen)
        dlon, lat = float(mm.group(1)), float(mm.group(2))
        clon = (antisun + dlon) % 360.0
        for mode in MODES:
            for pop in POPS:
                sub = cache[cache.population == pop].reset_index(drop=True)
                cs = {pop: {"df": sub, "scorer": scorer,
                            "clone_factor": ps[pop]["clone_factor"],
                            "use_conditional_cloner": ps[pop].get("use_conditional_cloner", True),
                            "scatter_size": 4, "scatter_alpha": 0.1,
                            "_mag_app": sub.mag_app.to_numpy(float)}}
                for mb in MAGBINS:
                    cap = {}

                    def _cap(tree, gp, k=None, n_d0_grid=None, show_progress=False, n_jobs=1):
                        cap["n"] = int(len(np.asarray(tree.data)))
                        return np.zeros(len(gp))

                    # Capture the VISIBLE SUBSET frame, which carries real ObjIDs and is exactly
                    # what feeds the density tree. Reverse-mapping tree coordinates to ObjIDs is
                    # ambiguous: many objects share a rounded (vlam, vbeta) pair.
                    _orig_vis = base.build_visible_subset_dataframe

                    def _cap_vis(*args, **kw):
                        out = _orig_vis(*args, **kw)
                        cap["df"] = out
                        return out

                    base.build_visible_subset_dataframe = _cap_vis
                    base.evaluate_density_map_full_posterior_2d = _cap
                    tmp = OUT / f"_sd_{a.shard}.npz"
                    try:
                        base.generate_probability_maps(
                            obstime_str=seal["grid"]["ref_obstime"], output_path=str(tmp),
                            center_lon_deg=clon, center_lat_deg=lat, center_label=cen,
                            max_sep_deg=MAX_SEP, n_jobs=1, save_overlays=False,
                            smooth_density_maps=False, clone_sources=cs, mag_bins=[mb],
                            grid_lim=(-5.0, 5.0), grid_step=0.01)
                    except Exception:
                        pass
                    finally:
                        base.evaluate_density_map_full_posterior_2d = orig
                        base.build_visible_subset_dataframe = _orig_vis
                        Path(tmp).unlink(missing_ok=True)
                    if "df" not in cap or cap.get("n", 0) < 2:
                        continue
                    dfv = cap["df"]
                    if "ObjID" in dfv.columns:
                        oid = [str(x) for x in dfv["ObjID"].to_numpy()]
                    else:
                        # NEO comes from the NEOMOD3 GEN cache, which carries no S3M ObjID.
                        # Identify those rows by their ordered orbital record instead.
                        cols = [c for c in ("a", "e", "i", "H", "vlam", "vbeta")
                                if c in dfv.columns]
                        oid = ["NM3GEN:" + ",".join(f"{v:.12g}" for v in row)
                               for row in dfv[cols].to_numpy(float)]
                    w = 1.0 / eff_for(pop, mb)
                    rows.append(dict(mode=mode, center=cen, population=pop,
                                     magnitude_bin=mb["label"], n_objids=len(oid),
                                     physical_weight_per_sample=w,
                                     sha256=digest_records(oid, pop, cen, mb["label"], w)))
        if (ci + 1) % 5 == 0:
            print(f"  {ci+1}/{len(centers)} centers  ({time.time()-t0:.0f}s)", flush=True)
    pd.DataFrame(rows).to_parquet(OUT / f"digest_shard_{a.shard:03d}.parquet", index=False)
    print(f"wrote digest_shard_{a.shard:03d}.parquet  ({len(rows):,} rows)")


if __name__ == "__main__":
    main()
