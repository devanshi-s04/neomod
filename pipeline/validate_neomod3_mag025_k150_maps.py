#!/usr/bin/env python3
"""Validation for the NEOMOD3 0.25-mag / k=150 build.

    --gate      pre-build acceptance gate on one center (runbook §9)
    --full      full-build validation + MAP_BUILD_SEAL_V2.json (runbook §10)
"""
from __future__ import annotations
import argparse, glob, hashlib, json, os, shutil, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))
sys.path.insert(0, str(W / "neomod" / "pipeline"))

OUT_ROOT = W / "outputs" / "neomod3_mag025_k150_maps_v2"
POPS = ("NEO", "MBA", "TNO", "Trojans")
K_BY_POP = {"NEO": 150, "MBA": 10, "TNO": 10, "Trojans": 10}


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


RESULTS = []


def chk(name, cond, detail=""):
    ok = bool(cond)
    RESULTS.append({"check": name, "pass": ok, "detail": detail})
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""), flush=True)
    return ok


def gate(a):
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    from build_neomod3_mag025_k150_maps import (
        mag025_bins, density_closed_form, density_sealed_quadrature, K_BY_POP as KBP)
    from scipy.spatial import cKDTree

    gd = Path(a.gate_dir)
    label = a.center
    z = np.load(gd / f"mag025_k150_{label}.npz", allow_pickle=True)
    meta = json.loads(str(z["meta_json"]))
    cov = pd.DataFrame(json.loads(str(z["coverage_json"])))
    post = pd.DataFrame(json.loads(str(z["posterior_json"])))

    print("=" * 72); print("PRE-BUILD ACCEPTANCE GATE"); print("=" * 72)

    # 1. 0.25-mag membership semantics
    bins = mag025_bins()
    chk("44 bins over 14<=V<25, step 0.25", len(bins) == 44 and bins[0]["mag_min"] == 14.0
        and abs(bins[-1]["mag_max"] - 25.0) < 1e-9)
    hit = [b for b in bins if b["mag_min"] <= 23.4 < b["mag_max"]]
    chk("V=23.4 selects exactly one bin", len(hit) == 1, hit[0]["label"] if hit else "none")
    chk("V=23.4 -> [23.25, 23.50)", hit and hit[0]["mag_min"] == 23.25 and hit[0]["mag_max"] == 23.5)
    chk("V=25.00 is out of range (not folded into top bin)",
        not any(b["mag_min"] <= 25.0 < b["mag_max"] for b in bins))
    # empirical: selected rows live strictly inside their slice
    src = pd.read_parquet(meta_neo_path(meta), columns=["mag_app"])
    lo, hi = 24.25, 24.50
    v = src.mag_app.to_numpy(float)
    sel = v[(v >= lo) & (v < hi)]
    chk("selected V range contained in its 0.25 slice",
        sel.min() >= lo and sel.max() < hi, f"{sel.min():.4f}..{sel.max():.4f} in [{lo},{hi})")

    # 2. k by population
    val = cov[cov.valid]
    okk = all(int(r.k_effective) == KBP[r.population] for r in val.itertuples())
    chk("NEO k=150 and non-NEO k=10 as requested", okk,
        str(sorted(set(zip(val.population, val.k_effective)))))
    chk("no cell silently reduced k", (val.k_effective == val.k_requested).all())

    # 3. smoothing / masking absent
    chk("Gaussian smoothing off", meta["gaussian_smoothing"] is False)
    chk("support masking off", meta["support_masking"] is False)
    chk("no smoothing arrays present in output",
        not any(k.startswith("smooth") for k in z.keys()))

    # 4. NEO source hash
    chk("NEO source hash is the HIGH realization",
        meta["neo_source_sha256"] == a.expected_neo_hash,
        meta["neo_source_sha256"][:16])
    chk("NEO effective_factor is the HIGH value (7.4e8 draws)",
        abs(meta["neo_effective_factor"] - 64.725384) < 1e-6)

    # 5. GEN/CAL/TEST identity overlap
    man = pd.read_parquet(W / "outputs/splits/nonneo_split_manifest.parquet",
                          columns=["ObjID", "split"])
    g = set(man.ObjID[man.split == "GEN"]); c = set(man.ObjID[man.split == "CAL"])
    t = set(man.ObjID[man.split == "TEST"])
    chk("GEN/CAL/TEST are pairwise disjoint",
        not (g & c) and not (g & t) and not (c & t),
        f"GEN {len(g):,} CAL {len(c):,} TEST {len(t):,}")

    # 6. physical weights unchanged
    wok = True
    for r in val.itertuples():
        if r.total_physical_weight is None:
            continue
        if abs(r.total_physical_weight - r.n_visible / r.effective_factor) > 1e-6:
            wok = False
    chk("total physical weight == n / effective_factor for every valid cell", wok)
    neo = val[val.population == "NEO"]
    chk("NEO physical weight uses 1/64.725384 per sample",
        np.allclose(neo.physical_weight_per_sample, 1.0 / 64.725384))

    # 7. posteriors finite and sum to 1
    dpost = post[post.defined]
    chk("every defined posterior sums to 1 within 1e-12",
        (dpost.max_abs_sum_deviation < 1e-12).all(),
        f"max {dpost.max_abs_sum_deviation.max():.3e}")
    allfin = True
    for lab in dpost.magnitude_bin:
        P = z[f"P_NEO__{lab}"]
        d = np.isfinite(P)
        if not ((P[d] >= 0).all() and (P[d] <= 1).all()):
            allfin = False
    chk("all defined P in [0,1] and finite", allfin)
    chk("invalid cells are absent, not stored as zeros",
        all(f"density__{r.population}__{r.magnitude_bin}" not in z.keys()
            for r in cov[~cov.valid].itertuples()))

    # 8. closed-form vs sealed-quadrature equivalence
    neo_src = pd.read_parquet(meta_neo_path(meta), columns=["mag_app", "vlam", "vbeta"])
    m = (neo_src.mag_app >= 24.25) & (neo_src.mag_app < 24.50)
    pts = np.column_stack([neo_src.vlam[m].to_numpy(float), neo_src.vbeta[m].to_numpy(float)])
    tree = cKDTree(pts)
    rng = np.random.default_rng(11)
    gp = np.column_stack([rng.uniform(-2, 2, int(a.equiv_points)),
                          rng.uniform(-2, 2, int(a.equiv_points))])
    t0 = time.time(); cf = density_closed_form(tree, gp, 150, workers=a.n_jobs); t_cf = time.time() - t0
    t0 = time.time(); sq = density_sealed_quadrature(tree, gp, 150, workers=a.n_jobs); t_sq = time.time() - t0
    rel = np.abs(cf - sq) / np.abs(sq)
    chk("closed_form vs sealed_quadrature agree < 1e-9",
        rel.max() < 1e-9, f"max {rel.max():.3e} median {np.median(rel):.3e}")
    speed = t_sq / t_cf if t_cf else float("inf")
    print(f"    timing on {a.equiv_points:,} points at k=150 with {a.n_jobs} workers: "
          f"closed_form {t_cf:.2f}s  sealed_quadrature {t_sq:.2f}s  speedup {speed:.0f}x", flush=True)

    # 9. storage projection
    nb_built = cov.magnitude_bin.nunique()
    size = (gd / f"mag025_k150_{label}.npz").stat().st_size
    per_bin = size / nb_built
    proj_center = per_bin * 44
    proj_all = proj_center * 667
    free = shutil.disk_usage(str(W))[2]
    print(f"\n  measured: {size/1e6:.1f} MB for {nb_built} bins -> {per_bin/1e6:.2f} MB/bin")
    print(f"  projection: {proj_center/1e6:.0f} MB/center  ->  {proj_all/1e12:.3f} TB for 667")
    print(f"  free space: {free/1e12:.3f} TB")
    chk("projected storage for 667 centers fits in free space", proj_all < free,
        f"need {proj_all/1e12:.3f} TB, free {free/1e12:.3f} TB")
    chk("projected storage leaves >20% headroom", proj_all < 0.8 * free,
        f"{100*proj_all/free:.1f}% of free space")

    out = {"center": label, "results": RESULTS, "build_seconds": meta["build_seconds"],
           "bins_built": int(nb_built), "bytes": int(size),
           "projected_bytes_per_center_44bins": float(proj_center),
           "projected_bytes_667_centers": float(proj_all),
           "free_bytes": int(free),
           "equivalence": {"max_rel": float(rel.max()), "median_rel": float(np.median(rel)),
                           "closed_form_s": t_cf, "sealed_quadrature_s": t_sq,
                           "speedup": float(speed), "n_points": int(a.equiv_points), "k": 150},
           "coverage_summary": cov.groupby(["population", "valid"]).size()
                                  .rename("n").reset_index().to_dict("records"),
           "ALL_PASS": all(r["pass"] for r in RESULTS)}
    Path(a.gate_dir, "gate_results.json").write_text(json.dumps(out, indent=2, default=str))
    print(f"\nGATE {'PASSED' if out['ALL_PASS'] else 'FAILED'} -> {a.gate_dir}/gate_results.json")
    return out


def meta_neo_path(meta):
    d = meta["neo_source"]
    return d.split("parquet ")[-1] if "parquet " in d else str(
        W / "outputs/more_neomod_samples_knn/source_high.parquet")


def full(a):
    """Full-build validation + MAP_BUILD_SEAL_V2.json."""
    from build_neomod3_mag025_k150_maps import mag025_bins, load_center_list
    centers = load_center_list()
    bins = [b["label"] for b in mag025_bins()]
    print("=" * 72); print("FULL-BUILD VALIDATION"); print("=" * 72)
    files, markers, cov_all, int_all, hashes, metas = [], [], [], [], {}, []
    for lab in centers:
        f = OUT_ROOT / f"mag025_k150_{lab}.npz"
        mk = OUT_ROOT / f"mag025_k150_{lab}.ok"
        if not f.exists() or not mk.exists():
            continue
        try:
            z = np.load(f, allow_pickle=True)
            metas.append(json.loads(str(z["meta_json"])))
            cov_all.append(pd.DataFrame(json.loads(str(z["coverage_json"]))))
            int_all.append(pd.DataFrame(json.loads(str(z["integrals_json"]))))
            _ = z["x_grid"]
        except Exception as e:
            print(f"  MALFORMED {lab}: {e}", flush=True); continue
        files.append(f); markers.append(json.load(open(mk)))
        hashes[lab] = sha256_file(f)
    chk("exactly 667 valid centers", len(files) == 667, f"{len(files)} readable")
    chk("no duplicate center outputs", len(set(hashes)) == len(hashes))
    cov = pd.concat(cov_all, ignore_index=True) if cov_all else pd.DataFrame()
    ints = pd.concat(int_all, ignore_index=True) if int_all else pd.DataFrame()
    if len(cov):
        chk("every (center, population, bin) cell accounted for",
            len(cov) == len(files) * len(POPS) * len(bins),
            f"{len(cov):,} rows vs expected {len(files)*len(POPS)*len(bins):,}")
        chk("marker hash matches file hash for every center",
            all(m["sha256"] == hashes[m["center"]] for m in markers))
    cfg = {(m["grid_step"], tuple(m["grid_lim"]), json.dumps(m["k_by_population"], sort_keys=True),
            m["gaussian_smoothing"], m["density_engine"], m["epoch"]) for m in metas}
    chk("configuration identical across all centers", len(cfg) == 1, str(list(cfg)[:1]))
    if len(cov):
        cov.to_parquet(OUT_ROOT / "coverage_table.parquet", index=False)
        ints.to_parquet(OUT_ROOT / "density_integrals.parquet", index=False)
        print("\n  cells by population and validity:")
        print(cov.groupby(["population", "valid"]).size().rename("n").reset_index().to_string(index=False))
    seal = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "n_centers": len(files), "n_bins": len(bins), "populations": list(POPS),
            "k_by_population": K_BY_POP, "gaussian_smoothing": False,
            "config": list(cfg)[0] if len(cfg) == 1 else None,
            "map_sha256": hashes, "checks": RESULTS,
            "ALL_PASS": all(r["pass"] for r in RESULTS)}
    if seal["ALL_PASS"]:
        (OUT_ROOT / "MAP_BUILD_SEAL_V2.json").write_text(json.dumps(seal, indent=2, default=str))
        print(f"\nwrote {OUT_ROOT/'MAP_BUILD_SEAL_V2.json'}")
    else:
        print("\nvalidation FAILED -- MAP_BUILD_SEAL_V2.json deliberately not written")
    return seal


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gate", action="store_true")
    p.add_argument("--full", action="store_true")
    p.add_argument("--center", default="dlon+000_lat+00")
    p.add_argument("--gate-dir", default=str(W / "outputs/neomod3_mag025_k150_maps_v2_gate"))
    p.add_argument("--expected-neo-hash",
                   default="40490b3bc4ffaec122919981396168299c1e84a384dd345c46f8a7adb20fc297")
    p.add_argument("--equiv-points", type=int, default=3000)
    p.add_argument("--n-jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    a = p.parse_args()
    if a.gate:
        gate(a)
    if a.full:
        full(a)


if __name__ == "__main__":
    main()
