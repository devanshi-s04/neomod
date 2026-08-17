#!/usr/bin/env python3
"""Stage 1: verify the validated k=150 map set and write a READ-ONLY scoring seal.

Does NOT rebuild maps. Confirms the map root associated with validation job 38590098 and
MAP_BUILD_SEAL_V2.json, re-checks the frozen configuration independently of the build-time
validator, and records absolute paths, hashes, environment and versions for the scoring stage.

Stops with a non-zero exit if any check fails.
"""
from __future__ import annotations
import hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
MAPS = W / "outputs" / "neomod3_mag025_k150_maps_v2"
LEGACY_MAPS = W / "prob_maps_grid_neomod3_GEN_final"
ALLSKY = W / "outputs" / "neomod3_projection_cache_high_allsky"
OUT = W / "outputs" / "test2_geometric"
POPS = ("NEO", "MBA", "TNO", "Trojans")
K_BY_POP = {"NEO": 150, "MBA": 10, "TNO": 10, "Trojans": 10}
EXPECT_FEATURED = {"V024.00_024.25": 35580, "V024.25_024.50": 44721,
                   "V024.50_024.75": 56294, "V024.75_025.00": 70464}

R = []


def chk(name, cond, detail=""):
    ok = bool(cond)
    R.append({"check": name, "pass": ok, "detail": str(detail)})
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""), flush=True)
    return ok


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    print("=" * 72); print("STAGE 1 — MAP VERIFICATION AND SCORING SEAL"); print("=" * 72)
    seal_p = MAPS / "MAP_BUILD_SEAL_V2.json"
    chk("map root exists", MAPS.is_dir(), str(MAPS))
    chk("MAP_BUILD_SEAL_V2.json exists", seal_p.exists(), str(seal_p))
    if not seal_p.exists():
        sys.exit(1)
    seal = json.load(open(seal_p))
    rep = json.load(open(MAPS / "FULLGRID_BUILD_REPORT.json"))

    chk("build seal reports ALL_PASS", seal.get("ALL_PASS") is True)
    chk("build seal covers 667 centers", seal.get("n_centers") == 667, seal.get("n_centers"))
    chk("build seal records 44 bins", seal.get("n_bins") == 44, seal.get("n_bins"))
    chk("build seal k by population is frozen", seal.get("k_by_population") == K_BY_POP,
        json.dumps(seal.get("k_by_population")))
    chk("build seal records smoothing OFF", seal.get("gaussian_smoothing") is False)
    chk("build seal records support masking OFF", seal.get("support_masking") is False)
    cfg = seal.get("frozen_config") or {}
    chk("frozen grid is [-5,5] @ 0.01", cfg.get("grid_lim") == [-5.0, 5.0]
        and cfg.get("grid_step") == 0.01, f"{cfg.get('grid_lim')} @ {cfg.get('grid_step')}")
    chk("frozen effective_factor is 64.725384",
        abs(float(cfg.get("neo_effective_factor", 0)) - 64.725384) < 1e-6)
    chk("frozen density engine is the Bayesian kNN closed form",
        cfg.get("density_engine") == "closed_form", cfg.get("density_engine"))
    chk("frozen magnitude quantity is apparent V",
        "apparent V" in str(cfg.get("magnitude_quantity")), cfg.get("magnitude_quantity"))
    bs = cfg.get("bin_scheme", {})
    chk("bin scheme is 44 x 0.25 mag, half-open, 14..25",
        bs.get("n_bins") == 44 and bs.get("step") == 0.25 and bs.get("lo") == 14.0
        and bs.get("hi") == 25.0 and "half-open" in str(bs.get("semantics")), json.dumps(bs))

    # independent re-verification against the files themselves
    n_npz = len(list(MAPS.glob("mag025_k150_*.npz")))
    n_ok = len(list(MAPS.glob("mag025_k150_*.ok")))
    chk("667 map files on disk", n_npz == 667, n_npz)
    chk("667 manifests on disk", n_ok == 667, n_ok)
    chk("seal records a hash for every center", len(seal.get("map_sha256", {})) == 667,
        len(seal.get("map_sha256", {})))

    # spot-verify three map hashes against the seal (full 248 GB rehash is the build validator's job)
    import random
    rnd = random.Random(0)
    labs = rnd.sample(sorted(seal["map_sha256"]), 3)
    ok_h = True
    for lab in labs:
        f = MAPS / f"mag025_k150_{lab}.npz"
        h = sha256_file(f)
        if h != seal["map_sha256"][lab]:
            ok_h = False
        print(f"    {lab}: {h[:16]} vs seal {seal['map_sha256'][lab][:16]}", flush=True)
    chk("spot-checked map hashes match the seal", ok_h, ",".join(labs))

    # featured-cell counts and coverage
    cov = pd.read_parquet(MAPS / "coverage_table.parquet")
    chk("coverage table has 117,392 cells", len(cov) == 117392, f"{len(cov):,}")
    fc = cov[(cov.center == "dlon+000_lat+00") & (cov.population == "NEO")
             & (cov.magnitude_bin.isin(EXPECT_FEATURED))]
    got = {r.magnitude_bin: int(r.n_visible) for r in fc.itertuples()}
    chk("featured-cell counts match the completed build", got == EXPECT_FEATURED, json.dumps(got))
    chk("featured four-bin total is 207,059", sum(got.values()) == 207059, sum(got.values()))
    chk("no cell silently reduced k",
        ((cov.k_effective.isna()) | (cov.k_effective == cov.k_requested)).all())
    chk("requested k matches frozen values per population",
        all((cov[cov.population == p].k_requested == K_BY_POP[p]).all() for p in POPS))

    # all-sky input hashes match what the build consumed
    asrep = json.load(open(ALLSKY / "validation_report.json"))
    asman = json.load(open(ALLSKY / "manifest.json"))
    chk("all-sky GEN source validation passed", asrep.get("ALL_PASS") is True)
    chk("all-sky source records 7.4e8 draws", asman.get("total_draws") == 740_000_000)
    chk("all-sky source effective_factor is 64.725384",
        abs(float(asman["effective_factor_NEO"]) - 64.725384) < 1e-6)
    chk("maps record the all-sky validation-report hash as their NEO source",
        cfg.get("neo_source_sha256") == sha256_file(ALLSKY / "validation_report.json"),
        str(cfg.get("neo_source_sha256"))[:16])
    chk("all-sky featured counts match the map featured counts",
        asrep["featured_center_counts"] == EXPECT_FEATURED)

    coverage_summary = (cov.groupby(["population", "valid"]).size().unstack(fill_value=0)
                        .rename(columns={True: "valid", False: "invalid"}))
    neo_valid = int(coverage_summary.loc["NEO", "valid"])
    neo_total = int(coverage_summary.loc["NEO"].sum())
    print(f"\n  NEO valid cells: {neo_valid:,}/{neo_total:,} = {100*neo_valid/neo_total:.1f}%",
          flush=True)
    print("  (a TEST2 tracklet landing in an invalid cell scores NaN with reason, never 0)",
          flush=True)

    # environment
    def ver(mod):
        try:
            m = __import__(mod)
            return getattr(m, "__version__", "?")
        except Exception:
            return None
    # NOTE: W/"digest2" is the source DIRECTORY; the executable is W/"digest2"/"digest2",
    # which is how sorcha_test_score_shard.py has always invoked it.
    d2bin = W / "digest2" / "digest2"
    d2ver = ""
    try:
        d2 = subprocess.run([str(d2bin), "--version"], capture_output=True, text=True,
                            timeout=60)
        d2ver = (d2.stdout + d2.stderr).strip().splitlines()[0] if (d2.stdout or d2.stderr) else ""
    except Exception as e:
        d2ver = f"(unavailable: {e})"
    env = {"python": sys.version.split()[0], "platform": platform.platform(),
           "numpy": ver("numpy"), "pandas": ver("pandas"), "scipy": ver("scipy"),
           "pyarrow": ver("pyarrow"), "healpy": ver("healpy"), "astropy": ver("astropy"),
           "sklearn": ver("sklearn"), "matplotlib": ver("matplotlib"),
           "digest2_binary": str(d2bin),
           "digest2_model_dir": str(W / "digest2"),
           "digest2_sha256": sha256_file(d2bin) if d2bin.is_file() else None,
           "digest2_version_string": d2ver}
    print("\n  environment:", json.dumps(env, indent=2), flush=True)

    allpass = all(r["pass"] for r in R)
    OUT.mkdir(parents=True, exist_ok=True)
    out = {
        "stage": "1 — map verification and scoring seal",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "read_only": True,
        "git_commit_before": subprocess.run(
            ["git", "-C", str(W / "neomod"), "rev-parse", "HEAD"],
            capture_output=True, text=True).stdout.strip(),
        "map_root": str(MAPS),
        "map_build_seal": str(seal_p),
        "map_build_seal_sha256": sha256_file(seal_p),
        "map_build_validation_job": "38590098",
        "legacy_map_root": str(LEGACY_MAPS),
        "legacy_map_seal": str(W / "outputs/splits/MAP_BUILD_SEAL.json"),
        "allsky_source_root": str(ALLSKY),
        "allsky_validation_report_sha256": sha256_file(ALLSKY / "validation_report.json"),
        "allsky_manifest_sha256": sha256_file(ALLSKY / "manifest.json"),
        "frozen_config": cfg,
        "k_by_population": K_BY_POP,
        "n_centers": 667, "n_bins": 44, "n_cells": int(len(cov)),
        "coverage_by_population": coverage_summary.reset_index().to_dict("records"),
        "neo_valid_cell_fraction": neo_valid / neo_total,
        "featured_cell_counts": got,
        "map_sha256": seal["map_sha256"],
        "storage_bytes": int(rep.get("storage_bytes", 0)),
        "environment": env,
        "checks": R, "ALL_PASS": allpass,
    }
    p = OUT / "SCORING_SEAL.json"
    p.write_text(json.dumps(out, indent=2, default=str))
    os.chmod(p, 0o444)                     # read-only seal
    print(f"\nSTAGE 1 {'PASSED' if allpass else 'FAILED'} -> {p}", flush=True)
    if not allpass:
        sys.exit(2)


if __name__ == "__main__":
    main()
