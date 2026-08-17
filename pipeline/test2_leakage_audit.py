#!/usr/bin/env python3
"""READ-ONLY non-NEO leakage audit for TEST2, plus TEST2 NEO shard hashing.

Writes nothing except its own report. Builds no tracklets, scores nothing.

WHY THIS IS NOT A FORMALITY
---------------------------
`BM_SPLIT_ROLE=TEST` guarantees separation from map-building GEN and from CAL, because those are
disjoint partitions of the same manifest. It does NOT guarantee separation from the PRIOR inspected
geometric TEST, which was itself built with `BM_SPLIT_ROLE=TEST` -- i.e. from the same partition.
Any overlap with the prior TEST is a FAILED fresh-TEST2 gate.

RECOVERING THE PRIOR TEST'S SOURCE IDENTITIES
---------------------------------------------
`gen_benchmark_tracklets_neomod3._tracklet_frame` overwrites `ObjID` with generated `NM...` strings,
so the prior tracklet file does not carry its source S3M ObjID. It does, however, copy
`lam_deg`, `beta_deg` and `e` bit-identically out of the epoch-state cache
(`lam, beta = sub.lam_deg.to_numpy(float), sub.beta_deg.to_numpy(float)`), so the original ObjID is
recoverable by an exact-value join on those float64 columns. Join uniqueness is asserted, not
assumed; unmatched rows are reported rather than silently dropped.
"""
from __future__ import annotations
import glob, hashlib, json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
OUT = W / "outputs" / "test2_geometric"
SHARDS = OUT / "neo_shards"
EPOCH_CACHE = W / "outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"
MANIFEST = W / "outputs/splits/nonneo_split_manifest.parquet"
PRIOR_TEST = W / "outputs/test_tracklets_neomod3/tracklets_benchmark_neomod3.parquet"
NONNEO = ("MBA", "TNO", "Trojans")
V_LO, V_HI = 14.0, 25.0
DLON_LIMIT = 140.0
EPOCH = "2027-08-25T00:00:00"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    t0 = time.time()
    print("=" * 74); print("TEST2 NEO SHARD HASHES"); print("=" * 74)
    fs = sorted(glob.glob(str(SHARDS / "neo_shard_*.parquet")))
    shard_h, n_clones, drawn = {}, 0, 0
    for f in fs:
        d = pd.read_parquet(f, columns=["mag_app", "n_orbits_drawn"])
        n_clones += len(d)
        drawn += int(d.n_orbits_drawn.iloc[0]) if len(d) else 0
        shard_h[Path(f).name] = sha256_file(f)
    combined = hashlib.sha256(b"".join(bytes.fromhex(shard_h[k]) for k in sorted(shard_h))).hexdigest()
    print(f"  shards {len(fs)}   clones in [{V_LO},{V_HI}) {n_clones:,}   draws {drawn:,}")
    print(f"  combined shard sha256: {combined}")
    print(f"  ({time.time()-t0:.0f}s)", flush=True)

    print("\n" + "=" * 74); print("NON-NEO LEAKAGE AUDIT (read-only)"); print("=" * 74)
    man = pd.read_parquet(MANIFEST, columns=["ObjID", "population", "split"])
    gen = {p: set(man.ObjID[(man.split == "GEN") & (man.population == p)]) for p in NONNEO}
    cal = {p: set(man.ObjID[(man.split == "CAL") & (man.population == p)]) for p in NONNEO}
    tst = {p: set(man.ObjID[(man.split == "TEST") & (man.population == p)]) for p in NONNEO}

    cache = pd.read_parquet(EPOCH_CACHE,
                            columns=["ObjID", "population", "mag_app", "lam_deg", "beta_deg", "e"])
    # in-domain and in-grid, exactly as the benchmark builder selects
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    t = Time(EPOCH, scale="utc")
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    dlon = ((cache.lam_deg.to_numpy(float) - antisun + 180.0) % 360.0) - 180.0
    v = cache.mag_app.to_numpy(float)
    eligible = np.isfinite(v) & (v >= V_LO) & (v < V_HI) & (np.abs(dlon) <= DLON_LIMIT + 0.5)
    cache["eligible"] = eligible
    print(f"  epoch cache {len(cache):,} rows; eligible (14<=V<25, |dlon|<=140): "
          f"{int(eligible.sum()):,}", flush=True)

    # ---- recover the prior TEST's source ObjIDs by exact-value join ----------
    prior = pd.read_parquet(PRIOR_TEST, columns=["ObjID", "population", "lam_deg", "beta_deg", "e"])
    prior_nonneo = prior[prior.population.isin(NONNEO)]
    print(f"  prior TEST rows: {len(prior):,}  (non-NEO {len(prior_nonneo):,})", flush=True)
    # `e` is NaN for every non-NEO row in the prior TEST file, and NaN != NaN in a merge, so
    # including it silently produced ZERO matches and a meaningless "no overlap" result.
    # lam_deg/beta_deg are copied bit-identically from the cache and match 614,709/614,709.
    key = ["population", "lam_deg", "beta_deg"]
    ckey = cache[cache.population.isin(NONNEO)][key + ["ObjID"]].copy()
    dup = int(ckey.duplicated(subset=key).sum())
    # Ambiguous keys would inflate the join; drop them and report, never silently expand.
    ckey_u = ckey.drop_duplicates(subset=key, keep=False)
    print(f"  cache join-key duplicates dropped: {dup} "
          f"({len(ckey)-len(ckey_u):,} rows removed as ambiguous)", flush=True)
    ckey = ckey_u
    joined = prior_nonneo.merge(ckey, on=key, how="left", suffixes=("_gen", "_src"))
    matched = joined.ObjID_src.notna()
    print(f"  prior TEST non-NEO rows matched back to a source ObjID: "
          f"{int(matched.sum()):,}/{len(joined):,} ({100*matched.mean():.3f}%)", flush=True)
    prior_ids = {p: set(joined.ObjID_src[matched & (joined.population == p)]) for p in NONNEO}

    # ---- the exact set TEST2 would propose under BM_SPLIT_ROLE=TEST ---------
    rows, remaining = [], {}
    for p in NONNEO:
        sub = cache[(cache.population == p) & cache.eligible]
        proposed = set(sub.ObjID) & tst[p]        # what the builder would offer TEST2
        o_gen, o_cal, o_prior = proposed & gen[p], proposed & cal[p], proposed & prior_ids[p]
        rest = proposed - gen[p] - cal[p] - prior_ids[p]
        remaining[p] = len(rest)
        rows.append(dict(population=p,
                         test_split_parents=len(tst[p]),
                         eligible_test_parents=len(proposed),
                         overlap_GEN=len(o_gen), overlap_CAL=len(o_cal),
                         overlap_prior_TEST=len(o_prior),
                         prior_TEST_used=len(prior_ids[p]),
                         remaining_after_exclusions=len(rest)))
    tab = pd.DataFrame(rows)
    print("\n  OVERLAP OF THE PROPOSED TEST2 NON-NEO PARENTS")
    print(tab.to_string(index=False), flush=True)

    fail = tab.overlap_prior_TEST.sum() > 0
    print(f"\n  overlap with GEN       : {int(tab.overlap_GEN.sum()):,}")
    print(f"  overlap with CAL       : {int(tab.overlap_CAL.sum()):,}")
    print(f"  overlap with prior TEST: {int(tab.overlap_prior_TEST.sum()):,}")
    print(f"\n  FRESH-TEST2 GATE: {'FAILED' if fail else 'PASSED'}", flush=True)
    if fail:
        print("  Any overlap with the prior inspected TEST fails the gate. BM_SPLIT_ROLE=TEST")
        print("  separates from GEN/CAL only; the prior TEST drew from this same partition.")
    # A near-zero match fraction means the recovery FAILED, not that there is no overlap.
    if matched.mean() < 0.90:
        print(f"\n  WARNING: only {100*matched.mean():.3f}% of prior TEST rows were recovered; the")
        print("  overlap numbers above are NOT trustworthy and must not be read as a pass.")

    rep = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "read_only": True,
           "neo_shards": {"dir": str(SHARDS), "n_shards": len(fs), "n_clones": int(n_clones),
                          "n_draws": int(drawn), "seed_base": 777000000,
                          "seeds": [777000000 + i for i in range(len(fs))],
                          "combined_sha256": combined, "per_shard_sha256": shard_h},
           "eligibility": {"v_range": [V_LO, V_HI], "dlon_limit": DLON_LIMIT,
                           "n_eligible_cache_rows": int(eligible.sum())},
           "prior_test_recovery": {"file": str(PRIOR_TEST),
                                   "rows_nonneo": int(len(prior_nonneo)),
                                   "matched": int(matched.sum()),
                                   "match_fraction": float(matched.mean()),
                                   "cache_key_duplicates": dup,
                                   "method": "exact float64 join on (population, lam_deg, beta_deg, e)"},
           "overlap_table": tab.to_dict("records"),
           "remaining_after_all_exclusions": remaining,
           "fresh_test2_gate_passed": (not fail)}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "TEST2_LEAKAGE_AUDIT.json").write_text(json.dumps(rep, indent=2, default=str))
    tab.to_csv(OUT / "TEST2_LEAKAGE_AUDIT.csv", index=False)
    print(f"\nwrote {OUT/'TEST2_LEAKAGE_AUDIT.json'} and .csv  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
