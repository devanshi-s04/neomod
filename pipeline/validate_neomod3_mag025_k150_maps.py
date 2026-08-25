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
    """Full-build validation, coverage reporting, seal, and FULLGRID_BUILD_REPORT."""
    import subprocess, shutil
    from build_neomod3_mag025_k150_maps import mag025_bins, load_center_list, K_BY_POP as KBP
    centers = load_center_list()
    binlabels = [b["label"] for b in mag025_bins()]
    print("=" * 72); print("FULL-BUILD VALIDATION"); print("=" * 72)

    files, markers, cov_all, int_all, hashes, metas, missing = [], [], [], [], {}, [], []
    malformed = []
    for lab in centers:
        f = OUT_ROOT / f"mag025_k150_{lab}.npz"
        mk = OUT_ROOT / f"mag025_k150_{lab}.ok"
        if not f.exists() or not mk.exists():
            missing.append(lab); continue
        try:
            z = np.load(f, allow_pickle=True)
            m = json.loads(str(z["meta_json"]))
            cov_all.append(pd.DataFrame(json.loads(str(z["coverage_json"]))))
            int_all.append(pd.DataFrame(json.loads(str(z["integrals_json"]))))
            _ = z["x_grid"], z["y_grid"]
            nkeys = len([k for k in z.keys() if k.startswith("density__")])
        except Exception as e:
            malformed.append((lab, str(e))); continue
        metas.append(m); files.append(f); markers.append(json.load(open(mk)))
        hashes[lab] = sha256_file(f)

    chk("exactly 667 completed, readable map files", len(files) == 667,
        f"{len(files)} readable, {len(missing)} missing, {len(malformed)} malformed")
    chk("exactly 667 manifests (.ok markers)", len(markers) == 667, str(len(markers)))
    chk("no malformed maps", not malformed, str(malformed[:3]))
    chk("no duplicate map content across centers",
        len(set(hashes.values())) == len(hashes),
        f"{len(set(hashes.values()))} distinct of {len(hashes)}")
    chk("marker hash matches file hash for every center",
        all(m["sha256"] == hashes.get(m["center"]) for m in markers))

    # ---- frozen configuration recorded in EVERY map -------------------------
    def uni(key):
        return {json.dumps(m.get(key), sort_keys=True) for m in metas}
    chk("every map records the frozen NEO source hash", len(uni("neo_source_sha256")) == 1,
        list(uni("neo_source_sha256"))[0][:34] if metas else "")
    chk("every map records the same NEO effective_factor (physical weights)",
        len(uni("neo_effective_factor")) == 1 and
        abs(metas[0]["neo_effective_factor"] - 64.725384) < 1e-6,
        str(metas[0]["neo_effective_factor"]) if metas else "")
    chk("every map records the 44-bin 0.25-mag scheme", len(uni("bin_scheme")) == 1 and
        metas[0]["bin_scheme"]["n_bins"] == 44 and metas[0]["bin_scheme"]["step"] == 0.25,
        json.dumps(metas[0]["bin_scheme"]) if metas else "")
    chk("every map records population-specific k (NEO 150, others 10)",
        len(uni("k_by_population")) == 1 and metas[0]["k_by_population"] == KBP,
        json.dumps(metas[0]["k_by_population"]) if metas else "")
    chk("every map records the same velocity grid", len(uni("grid_lim")) == 1 and
        len(uni("grid_step")) == 1 and metas[0]["grid_lim"] == [-5.0, 5.0] and
        metas[0]["grid_step"] == 0.01)
    chk("every map records gaussian_smoothing = False",
        all(m["gaussian_smoothing"] is False for m in metas))
    chk("every map records support_masking = False",
        all(m["support_masking"] is False for m in metas))
    chk("every map records the same 0.25-mag split provenance hash",
        len(uni("split_provenance_mag025_sha256")) == 1)
    chk("every map records the same sealed-module hash", len(uni("sealed_module_sha256")) == 1)
    chk("every map records the same density engine", len(uni("density_engine")) == 1,
        metas[0]["density_engine"] if metas else "")
    chk("every map records the same epoch", len(uni("epoch")) == 1)
    chk("no smoothing arrays present anywhere",
        not any(k.startswith("smooth") for f in files[:1]
                for k in np.load(f, allow_pickle=True).keys()))

    cov = pd.concat(cov_all, ignore_index=True) if cov_all else pd.DataFrame()
    ints = pd.concat(int_all, ignore_index=True) if int_all else pd.DataFrame()
    exp_cells = len(files) * len(POPS) * len(binlabels)
    chk("every (center, population, bin) cell accounted for", len(cov) == exp_cells,
        f"{len(cov):,} vs expected {exp_cells:,}")
    if len(cov):
        chk("no cell silently reduced k",
            ((cov.k_effective.isna()) | (cov.k_effective == cov.k_requested)).all())
        chk("every invalid cell carries an explicit reason",
            (cov.loc[~cov.valid, "reason"] != "ok").all(),
            str(sorted(cov.loc[~cov.valid, "reason"].unique().tolist())))
        chk("requested k matches the frozen per-population values",
            all((cov[cov.population == p].k_requested == KBP[p]).all() for p in POPS))

    # ---- featured-center recheck -------------------------------------------
    EXPECT = {"V024.00_024.25": 35580, "V024.25_024.50": 44721,
              "V024.50_024.75": 56294, "V024.75_025.00": 70464}
    fc = cov[(cov.center == "dlon+000_lat+00") & (cov.population == "NEO") &
             (cov.magnitude_bin.isin(EXPECT))] if len(cov) else pd.DataFrame()
    got = {r.magnitude_bin: int(r.n_visible) for r in fc.itertuples()} if len(fc) else {}
    for k2 in EXPECT:
        print(f"    featured {k2}: got {got.get(k2, 0):,}  expected {EXPECT[k2]:,}  "
              f"{'MATCH' if got.get(k2) == EXPECT[k2] else 'MISMATCH'}", flush=True)
    chk("dlon+000_lat+00 NEO counts in the four 24-25 bins match validated values",
        got == EXPECT, json.dumps(got))
    chk("featured-center four-bin total is 207,059", sum(got.values()) == 207059,
        str(sum(got.values())))

    # ---- coverage reporting -------------------------------------------------
    by_pop = by_bin = by_center = pd.DataFrame()
    if len(cov):
        by_pop = (cov.groupby(["population", "valid"]).size().unstack(fill_value=0)
                  .rename(columns={True: "valid", False: "invalid"}).reset_index())
        by_bin = (cov.groupby(["magnitude_bin", "population"]).valid.sum().unstack(fill_value=0)
                  .reset_index())
        by_center = (cov.groupby("center").valid.agg(["sum", "count"])
                     .rename(columns={"sum": "valid_cells", "count": "cells"}).reset_index())
        cov.to_parquet(OUT_ROOT / "coverage_table.parquet", index=False)
        ints.to_parquet(OUT_ROOT / "density_integrals.parquet", index=False)
        by_pop.to_csv(OUT_ROOT / "coverage_by_population.csv", index=False)
        by_bin.to_csv(OUT_ROOT / "coverage_by_magnitude_bin.csv", index=False)
        by_center.to_csv(OUT_ROOT / "coverage_by_center.csv", index=False)
        print("\n  VALID / INVALID CELLS BY POPULATION")
        print(by_pop.to_string(index=False))
        print("\n  REASONS FOR INVALID CELLS")
        print(cov.loc[~cov.valid].groupby(["population", "reason"]).size()
              .rename("n").reset_index().to_string(index=False))
        print("\n  VALID BINS PER CENTER (min/median/max)")
        print(f"    {by_center.valid_cells.min()} / {int(by_center.valid_cells.median())} "
              f"/ {by_center.valid_cells.max()} of {len(POPS)*len(binlabels)}")

    # ---- job accounting -----------------------------------------------------
    jobinfo = {}
    if a.array_job_id:
        r = subprocess.run(["sacct", "-j", str(a.array_job_id), "-n", "-P",
                            "--format=JobID,State,Elapsed,MaxRSS"],
                           capture_output=True, text=True).stdout.strip().splitlines()
        rows = [x.split("|") for x in r if x and "." not in x.split("|")[0]]
        states = {}
        for jid, st, el, _mx in rows:
            states[st] = states.get(st, 0) + 1
        secs = []
        for jid, st, el, _mx in rows:
            try:
                parts = [float(x) for x in el.replace("-", ":").split(":")]
                while len(parts) < 3:
                    parts.insert(0, 0.0)
                secs.append(parts[-3] * 3600 + parts[-2] * 60 + parts[-1])
            except Exception:
                pass
        jobinfo = {"array_job_id": str(a.array_job_id), "task_states": states,
                   "n_tasks": len(rows),
                   "elapsed_seconds": {"min": min(secs) if secs else None,
                                       "median": float(np.median(secs)) if secs else None,
                                       "max": max(secs) if secs else None,
                                       "sum": float(np.sum(secs)) if secs else None}}
        chk("no failed array tasks",
            not any(k for k in states if k not in ("COMPLETED",)), json.dumps(states))

    size = sum(f.stat().st_size for f in OUT_ROOT.glob("*.npz"))
    free = shutil.disk_usage(str(W))[2]

    seal = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "n_centers": len(files), "n_bins": len(binlabels), "populations": list(POPS),
            "k_by_population": K_BY_POP, "gaussian_smoothing": False, "support_masking": False,
            "frozen_config": metas[0] if metas else None,
            "provenance_commit": a.provenance_commit,
            "map_sha256": hashes, "checks": RESULTS,
            "ALL_PASS": all(r["pass"] for r in RESULTS)}
    if seal["ALL_PASS"]:
        (OUT_ROOT / "MAP_BUILD_SEAL_V2.json").write_text(json.dumps(seal, indent=2, default=str))
        print(f"\nwrote {OUT_ROOT/'MAP_BUILD_SEAL_V2.json'}")
    else:
        print("\nvalidation FAILED -- MAP_BUILD_SEAL_V2.json deliberately not written")

    rep = {"created_utc": seal["created_utc"], "n_centers": len(files),
           "missing_centers": missing, "malformed_centers": malformed,
           "n_cells": int(len(cov)), "n_valid_cells": int(cov.valid.sum()) if len(cov) else 0,
           "coverage_by_population": by_pop.to_dict("records") if len(by_pop) else [],
           "invalid_reasons": (cov.loc[~cov.valid].groupby(["population", "reason"]).size()
                               .rename("n").reset_index().to_dict("records")) if len(cov) else [],
           "featured_center_recheck": {"got": got, "expected": EXPECT,
                                       "total": int(sum(got.values()))},
           "storage_bytes": int(size), "free_bytes": int(free),
           "job": jobinfo, "checks": RESULTS,
           "ALL_PASS": seal["ALL_PASS"]}
    (OUT_ROOT / "FULLGRID_BUILD_REPORT.json").write_text(json.dumps(rep, indent=2, default=str))
    _write_report_md(OUT_ROOT / "FULLGRID_BUILD_REPORT.md", rep, by_pop, by_bin, by_center,
                     metas[0] if metas else {}, a)
    print(f"wrote {OUT_ROOT/'FULLGRID_BUILD_REPORT.json'} and .md")
    return seal


def _md(df):
    """Markdown table, falling back to plain text if `tabulate` is unavailable."""
    if df is None or not len(df):
        return "(none)"
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def _write_report_md(path, rep, by_pop, by_bin, by_center, meta, a):
    L = []
    L.append("# NEOMOD3 0.25-mag / k=150 full-grid build report\n")
    L.append(f"Generated {rep['created_utc']}. "
             f"Overall: **{'PASS' if rep['ALL_PASS'] else 'FAIL'}**\n")
    L.append("## Frozen configuration\n")
    L.append("| item | value |\n|---|---|")
    for k in ("epoch", "grid_lim", "grid_step", "k_by_population", "gaussian_smoothing",
              "support_masking", "density_engine", "magnitude_quantity", "neo_source",
              "neo_source_sha256", "neo_effective_factor", "split_provenance_mag025_sha256",
              "sealed_module_sha256", "code_commit"):
        if k in meta:
            L.append(f"| `{k}` | `{meta[k]}` |")
    L.append(f"| provenance commit | `{a.provenance_commit}` |")
    L.append("")
    L.append("## Scale\n")
    L.append(f"- centers: **{rep['n_centers']}**")
    L.append(f"- cells: **{rep['n_cells']:,}**, valid **{rep['n_valid_cells']:,}**")
    L.append(f"- storage: **{rep['storage_bytes']/1e9:.1f} GB**, "
             f"free after: **{rep['free_bytes']/1e12:.3f} TB**")
    if rep["job"]:
        j = rep["job"]; e = j.get("elapsed_seconds", {})
        L.append(f"- array job: **{j['array_job_id']}**, tasks {j['n_tasks']}, "
                 f"states {j['task_states']}")
        if e.get("median"):
            L.append(f"- per-task elapsed: min {e['min']:.0f}s, median {e['median']:.0f}s, "
                     f"max {e['max']:.0f}s, total CPU-wall {e['sum']/3600:.1f} h")
    L.append(f"- missing centers: {rep['missing_centers'] or 'none'}")
    L.append(f"- malformed centers: {rep['malformed_centers'] or 'none'}")
    L.append("")
    L.append("## Valid / invalid cells by population\n")
    L.append(_md(by_pop))
    L.append("")
    L.append("## Reasons for invalid cells\n")
    if rep["invalid_reasons"]:
        L.append(_md(pd.DataFrame(rep["invalid_reasons"])))
    L.append("")
    L.append("## Valid cells by magnitude bin and population\n")
    L.append(_md(by_bin))
    L.append("")
    L.append("## Featured-center recheck (dlon+000_lat+00, NEO)\n")
    fc = rep["featured_center_recheck"]
    L.append("| bin | expected | got | |\n|---|---|---|---|")
    for k2, v in fc["expected"].items():
        g = fc["got"].get(k2, 0)
        L.append(f"| {k2} | {v:,} | {g:,} | {'MATCH' if g == v else 'MISMATCH'} |")
    L.append(f"| **combined** | 207,059 | {fc['total']:,} | "
             f"{'MATCH' if fc['total'] == 207059 else 'MISMATCH'} |")
    L.append("")
    L.append("## Validation checks\n")
    L.append("| check | result | detail |\n|---|---|---|")
    for c in rep["checks"]:
        L.append(f"| {c['check']} | {'PASS' if c['pass'] else 'FAIL'} | {c['detail']} |")
    L.append("")
    L.append("## Per-center valid-cell distribution\n")
    if len(by_center):
        L.append(f"- min {by_center.valid_cells.min()}, "
                 f"median {int(by_center.valid_cells.median())}, "
                 f"max {by_center.valid_cells.max()} of {by_center.cells.iloc[0]} cells")
    path.write_text("\n".join(L))


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
    p.add_argument("--array-job-id", default=None)
    p.add_argument("--provenance-commit", default="e4780ef806d4eb4a77e3d384144f8a29ff24ea1e")
    a = p.parse_args()
    if a.gate:
        gate(a)
    if a.full:
        full(a)


if __name__ == "__main__":
    main()
