#!/usr/bin/env python3
"""TEST2 assembly: fresh NEOMOD3 NEOs + frozen reused prior-TEST S3M non-NEOs.

PROVENANCE, STATED NOT HIDDEN
-----------------------------
The strict fresh-non-NEO gate FAILED (TEST2_LEAKAGE_AUDIT.json): the prior inspected TEST consumed
every eligible TEST-split S3M parent, so the non-NEO rows here are REUSED, not independent. This
assembly is approved on that basis. Its purpose is a PAIRED comparison of new VDP, legacy VDP and
digest2 on IDENTICAL contaminant rows; absolute contamination/precision/FPR are therefore not a
fresh out-of-sample estimate. Only the NEO side is fresh.

IDENTITY IS PRESERVED
---------------------
The legacy builder overwrote source IDs with generated `NM...` strings, which is what made the
leakage audit require a float-join reconstruction. Every row here keeps:
    source_parent_uid     stable, population-prefixed source identity
    s3m_objid             original S3M ObjID for non-NEO (None for NEO)
    neo_shard/neo_row/neo_seed/neo_orbit_digest   for NEO
    tracklet_uid          immutable, independent of row order

WEIGHTS
-------
NEO is oversampled ~10x its physical expectation, so each NEO row carries
    w_phys = expected_physical_NEO / n_sampled_NEO
computed EXACTLY from the realised counts (never assumed to be 0.1). Non-NEO rows carry 1.0.
"""
from __future__ import annotations
import argparse, glob, hashlib, json, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))
sys.path.insert(0, str(W / "neomod" / "pipeline"))

OUT = W / "outputs" / "test2_geometric"
SHARDS = OUT / "neo_shards"
EPOCH_CACHE = W / "outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"
MANIFEST = W / "outputs/splits/nonneo_split_manifest.parquet"
NM3_META = W / "outputs/neomod3_projection_cache/cache_metadata.json"
SPLIT_PROV = W / "outputs/splits/split_provenance.json"
MAPS_V2 = W / "outputs/neomod3_mag025_k150_maps_v2"

EPOCH = "2027-08-25T00:00:00"
DT_DAYS = 30.0 / 1440.0
V_LO, V_HI = 14.0, 25.0
DLON_LIMIT = 140.0
NEO_OVERSAMPLE = 10
SAMPLE_SEED = 777123
NONNEO = ("MBA", "TNO", "Trojans")
MAG_STEP = 0.25


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def mag_bin_label(v):
    """Half-open [lo, lo+0.25); NaN outside [14,25) -- never clipped."""
    v = np.asarray(v, float)
    idx = np.floor((v - V_LO) / MAG_STEP).astype("float64")
    bad = ~np.isfinite(v) | (v < V_LO) | (v >= V_HI)
    idx[bad] = np.nan
    lab = np.array([None] * len(v), dtype=object)
    ok = ~np.isnan(idx)
    lo = V_LO + idx[ok] * MAG_STEP
    lab[ok] = [f"V{a:06.2f}_{a+MAG_STEP:06.2f}" for a in lo]
    return lab, bad


def center_grid():
    lats = sorted({float(v) for v in [0, 1, 2, 3, 4, 5, 8, 12, 18, 25, 35, 50]} |
                  {float(-v) for v in [0, 1, 2, 3, 4, 5, 8, 12, 18, 25, 35, 50]})
    dl = [round(d, 6) for d in np.arange(-180.0, 180.0, 10.0) if abs(d) <= DLON_LIMIT + 1e-9]
    return np.array([(d, l) for l in lats for d in dl], dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    a = ap.parse_args()
    t0 = time.time()
    from astropy.time import Time
    from astropy.coordinates import GeocentricTrueEcliptic, SkyCoord, get_sun
    from astropy.utils import iers; iers.conf.auto_max_age = None
    import astropy.units as u
    from sorcha_phase2 import format_mpc80

    t = Time(EPOCH, scale="utc"); mjd_ref = float(t.mjd)
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    grid = center_grid()
    print(f"epoch {EPOCH}  antisun {antisun:.4f}  centers {len(grid)}", flush=True)

    # ---------------- NEO: fresh independent realization ---------------------
    fs = sorted(glob.glob(str(SHARDS / "neo_shard_*.parquet")))
    frames = []
    for i, f in enumerate(fs):
        d = pd.read_parquet(f)
        d["neo_shard"] = i
        d["neo_row"] = np.arange(len(d), dtype=np.int64)
        d["neo_seed"] = np.int64(777000000 + i)
        frames.append(d)
    neo = pd.concat(frames, ignore_index=True); del frames
    n_drawn = int(neo.n_orbits_drawn.iloc[0]) * len(fs)
    print(f"NEO clones {len(neo):,} from {n_drawn:,} draws ({time.time()-t0:.0f}s)", flush=True)

    c = SkyCoord(ra=neo.ra_deg.to_numpy() * u.deg, dec=neo.dec_deg.to_numpy() * u.deg, frame="icrs")
    e = c.transform_to(GeocentricTrueEcliptic(obstime=t))
    neo["lam_deg"], neo["beta_deg"] = e.lon.deg, e.lat.deg
    neo["dlon"] = ((neo.lam_deg - antisun + 180.0) % 360.0) - 180.0
    neo = neo[np.abs(neo.dlon) <= DLON_LIMIT + 0.5].reset_index(drop=True)

    total_neo_abs = float(json.load(open(NM3_META))["total_weight_absolute_NEO_count"])
    w_new = total_neo_abs / n_drawn
    cnt = json.load(open(SPLIT_PROV))["counts"]
    n_test = sum(cnt[p]["TEST"] for p in NONNEO)
    n_all = sum(sum(cnt[p].values()) for p in NONNEO)
    role_fraction = n_test / n_all
    expected_neo = len(neo) * w_new * role_fraction
    n_target = int(round(NEO_OVERSAMPLE * expected_neo))
    if n_target > len(neo):
        raise RuntimeError(f"cannot sample {n_target:,} NEO rows without replacement "
                           f"from {len(neo):,} clones")
    rng = np.random.default_rng(SAMPLE_SEED)
    pick = rng.choice(len(neo), size=n_target, replace=False)
    assert len(set(pick.tolist())) == n_target, "NEO sampling drew a duplicate"
    neo_s = neo.iloc[np.sort(pick)].reset_index(drop=True)
    w_neo = expected_neo / len(neo_s)
    print(f"NEO in-grid {len(neo):,}  w_new {w_new:.8f}  role_fraction {role_fraction:.10f}",
          flush=True)
    print(f"expected physical NEO {expected_neo:,.4f}  sampled {len(neo_s):,}  "
          f"w_phys {w_neo:.10f}  (10x check: {expected_neo/len(neo_s)*NEO_OVERSAMPLE:.6f})",
          flush=True)

    dig = [hashlib.sha256(("|".join(f"{v:.12g}" for v in r)).encode()).hexdigest()[:24]
           for r in neo_s[["a", "e", "i", "node", "argperi", "t_p", "H"]].to_numpy(float)]
    neo_s["neo_orbit_digest"] = dig
    neo_s["source_parent_uid"] = [f"NEO2:{s}:{r}:{d}" for s, r, d in
                                  zip(neo_s.neo_shard, neo_s.neo_row, dig)]
    neo_s["s3m_objid"] = None
    neo_s["w_phys"] = w_neo
    neo_s["population"] = "NEO"

    # ---------------- non-NEO: frozen REUSED prior-TEST parents --------------
    cache = pd.read_parquet(EPOCH_CACHE)
    man = pd.read_parquet(MANIFEST, columns=["ObjID", "split"])
    test_ids = set(man.ObjID[man.split == "TEST"])
    gen_ids = set(man.ObjID[man.split == "GEN"])
    cal_ids = set(man.ObjID[man.split == "CAL"])
    parts = []
    for p in NONNEO:
        sub = cache[cache.population == p]
        v = sub.mag_app.to_numpy(float)
        dl = ((sub.lam_deg.to_numpy(float) - antisun + 180.0) % 360.0) - 180.0
        ok = (np.isfinite(v) & (v >= V_LO) & (v < V_HI) & (np.abs(dl) <= DLON_LIMIT + 0.5)
              & sub.ObjID.isin(test_ids).to_numpy())
        s = sub[ok].reset_index(drop=True)
        s["dlon"] = dl[ok]
        s["source_parent_uid"] = "S3M:" + s.ObjID.astype(str)
        s["s3m_objid"] = s.ObjID.astype(str)
        s["w_phys"] = 1.0
        for col in ("neo_shard", "neo_row", "neo_seed", "neo_orbit_digest"):
            s[col] = None
        parts.append(s)
        print(f"  {p}: {len(s):,} eligible TEST parents", flush=True)
    nn = pd.concat(parts, ignore_index=True); del parts, cache

    leak_gen = int(nn.s3m_objid.isin(gen_ids).sum())
    leak_cal = int(nn.s3m_objid.isin(cal_ids).sum())
    print(f"\nHARD GATE  GEN leakage {leak_gen}   CAL leakage {leak_cal}", flush=True)
    if leak_gen or leak_cal:
        raise RuntimeError(f"GEN/CAL leakage detected: GEN {leak_gen}, CAL {leak_cal}")

    # ---------------- unified frame + tracklets ------------------------------
    keep = ["source_parent_uid", "s3m_objid", "population", "neo_shard", "neo_row", "neo_seed",
            "neo_orbit_digest", "w_phys", "ra_deg", "dec_deg", "dra_deg_day", "ddec_deg_day",
            "mag_app", "lam_deg", "beta_deg", "dlon", "H", "vlam", "vbeta"]
    for c_ in keep:
        if c_ not in neo_s.columns:
            neo_s[c_] = np.nan
        if c_ not in nn.columns:
            nn[c_] = np.nan
    df = pd.concat([neo_s[keep], nn[keep]], ignore_index=True)
    assert df.source_parent_uid.is_unique, "source_parent_uid is not unique"

    df["tracklet_uid"] = ["T2:" + hashlib.sha256(u_.encode()).hexdigest()[:32]
                          for u_ in df.source_parent_uid]
    assert df.tracklet_uid.is_unique, "tracklet_uid is not unique"

    ci = np.argmin((df.dlon.to_numpy(float)[:, None] - grid[None, :, 0]) ** 2 +
                   (df.beta_deg.to_numpy(float)[:, None] - grid[None, :, 1]) ** 2, axis=1)
    cd, cl = grid[ci, 0], grid[ci, 1]
    df["center_label"] = [f"dlon{int(round(x)):+04d}_lat{int(round(y)):+03d}" for x, y in zip(cd, cl)]
    df["map_file"] = "mag025_k150_" + df.center_label + ".npz"
    df["legacy_map_file"] = [f"prob_maps_grid_dlon{int(round(x)):+04d}_lat{int(round(y)):+03d}.npz"
                             for x, y in zip(cd, cl)]
    lab, bad = mag_bin_label(df.mag_app.to_numpy(float))
    df["magnitude_bin"] = lab
    df["v_out_of_range"] = bad
    print(f"\nV outside [14,25): {int(bad.sum())} (flagged, not clipped)", flush=True)

    ra0 = df.ra_deg.to_numpy(float); dec0 = df.dec_deg.to_numpy(float)
    dra = df.dra_deg_day.to_numpy(float); ddec = df.ddec_deg_day.to_numpy(float)
    df["mjd0_utc"], df["mjd1_utc"] = mjd_ref, mjd_ref + DT_DAYS
    df["ra0"], df["dec0"] = ra0, dec0
    df["ra1"] = (ra0 + dra * DT_DAYS) % 360.0
    df["dec1"] = np.clip(dec0 + ddec * DT_DAYS, -89.99, 89.99)
    df["mag0_V"] = df["mag1_V"] = df["mean_mag_V"] = df.mag_app.to_numpy(float)
    df["n_det_per_night"] = 2

    # ---------------- MPC-80 -------------------------------------------------
    # MPC-80 columns 1-12 hold the identifier and format_mpc80 TRUNCATES to 12 chars. With 5
    # leading spaces only 7 characters survive, so an 8-char key silently collided: T0000000,
    # T0000001, T0000002 all became "T000000" and digest2 saw one object with many observations,
    # returning nothing that matched. 7 chars (T000000..T688687) is unique for 688,688 rows.
    keys = [f"T{i:06d}" for i in range(len(df))]
    if len(set(keys)) != len(keys):
        raise RuntimeError("d2_key collision before formatting")
    df["d2_key"] = keys
    l0 = [format_mpc80(f"     {k}", m, r, d, g) for k, m, r, d, g in
          zip(keys, df.mjd0_utc, df.ra0, df.dec0, df.mag0_V)]
    l1 = [format_mpc80(f"     {k}", m, r, d, g) for k, m, r, d, g in
          zip(keys, df.mjd1_utc, df.ra1, df.dec1, df.mag1_V)]
    bad80 = [i for i, (x, y) in enumerate(zip(l0, l1)) if len(x) != 80 or len(y) != 80]
    if bad80:
        raise RuntimeError(f"{len(bad80)} MPC-80 lines are not 80 characters")
    # THE check that catches truncation: the identifier actually present in the emitted line must be
    # unique, and both detections of a tracklet must carry the same one.
    des0 = [x[:12] for x in l0]; des1 = [y[:12] for y in l1]
    if len(set(des0)) != len(des0):
        raise RuntimeError(f"MPC-80 designations collide: {len(des0)-len(set(des0))} duplicates "
                           "(identifier truncated to 12 columns)")
    if des0 != des1:
        raise RuntimeError("the two detections of a tracklet carry different MPC-80 designations")
    df["mpc80_designation"] = des0
    df["mpc80_line0"], df["mpc80_line1"] = l0, l1
    df["tracklet_input_sha256"] = [hashlib.sha256((x + "\n" + y).encode()).hexdigest()
                                   for x, y in zip(l0, l1)]

    # ---------------- coverage against the v2 coverage table -----------------
    cov = pd.read_parquet(MAPS_V2 / "coverage_table.parquet")
    valid = cov[cov.valid][["center", "population", "magnitude_bin"]]
    neo_ok = set(map(tuple, valid[valid.population == "NEO"][["center", "magnitude_bin"]].to_numpy()))
    pairs = list(zip(df.center_label, df.magnitude_bin))
    df["vdp_cell_valid"] = [pr in neo_ok for pr in pairs]
    df.loc[df.v_out_of_range, "vdp_cell_valid"] = False
    reason = np.where(df.v_out_of_range, "v_out_of_range",
                      np.where(df.vdp_cell_valid, "ok", "neo_cell_invalid"))
    df["vdp_expected_reason"] = reason

    print("\nVDP COVERAGE (raw rows and physical-weighted)", flush=True)
    for p in ("NEO",) + NONNEO:
        g = df[df.population == p]
        raw = g.vdp_cell_valid.mean() if len(g) else np.nan
        wt = (g.w_phys[g.vdp_cell_valid].sum() / g.w_phys.sum()) if len(g) else np.nan
        print(f"  {p:8s} rows {len(g):>8,}  raw {100*raw:6.2f}%  weighted {100*wt:6.2f}%", flush=True)
    tot_raw = df.vdp_cell_valid.mean()
    tot_wt = df.w_phys[df.vdp_cell_valid].sum() / df.w_phys.sum()
    print(f"  {'ALL':8s} rows {len(df):>8,}  raw {100*tot_raw:6.2f}%  weighted {100*tot_wt:6.2f}%",
          flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    pm = df[["source_parent_uid", "s3m_objid", "population", "neo_shard", "neo_row", "neo_seed",
             "neo_orbit_digest", "w_phys", "tracklet_uid"]]
    pm.to_parquet(OUT / "TEST2_PARENT_MANIFEST.parquet", index=False)
    tk = df.drop(columns=["mpc80_line0", "mpc80_line1"])
    tk.to_parquet(OUT / "TEST2_TRACKLETS.parquet", index=False)
    df[["tracklet_uid", "d2_key", "mpc80_line0", "mpc80_line1",
        "tracklet_input_sha256"]].to_parquet(OUT / "TEST2_MPC80.parquet", index=False)

    counts = df.groupby("population").agg(rows=("tracklet_uid", "size"),
                                          weight=("w_phys", "sum")).reset_index()
    print("\nCOUNTS BY POPULATION"); print(counts.to_string(index=False), flush=True)

    seal = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "design": "fresh independent TEST2 NEOs + frozen REUSED prior-TEST S3M non-NEOs",
            "fresh_nonneo_gate": "FAILED (see TEST2_LEAKAGE_AUDIT.json) — reuse is intentional",
            "limitation": ("non-NEO parents were used in the previously inspected TEST, so absolute "
                           "contamination/precision/FPR are not independent of that look; the paired "
                           "classifier comparison on identical rows is the valid quantity"),
            "epoch": EPOCH, "dt_days": DT_DAYS, "v_range": [V_LO, V_HI],
            "dlon_limit": DLON_LIMIT, "mag_step": MAG_STEP,
            "neo": {"seed_base": 777000000, "n_shards": len(fs), "n_draws": n_drawn,
                    "n_clones_allsky": int(len(neo)) , "w_new": w_new,
                    "role_fraction": role_fraction,
                    "expected_physical_neo": expected_neo,
                    "oversample_factor": NEO_OVERSAMPLE,
                    "n_sampled": int(len(neo_s)), "w_phys": w_neo,
                    "sample_seed": SAMPLE_SEED, "sampling": "without replacement"},
            "nonneo": {"source": "epoch_state_cache, split=TEST", "w_phys": 1.0,
                       "gen_leakage": leak_gen, "cal_leakage": leak_cal},
            "counts": counts.to_dict("records"),
            "coverage": {"raw": float(tot_raw), "weighted": float(tot_wt)},
            "v_out_of_range_rows": int(bad.sum()),
            "artifacts": {n: {"path": str(OUT / n), "sha256": sha256_file(OUT / n)}
                          for n in ("TEST2_PARENT_MANIFEST.parquet", "TEST2_TRACKLETS.parquet",
                                    "TEST2_MPC80.parquet")},
            "map_root": str(MAPS_V2),
            "scoring_seal": str(OUT / "SCORING_SEAL.json")}
    (OUT / "TEST2_SEAL.json").write_text(json.dumps(seal, indent=2, default=str))
    print(f"\nwrote TEST2 artifacts + seal in {time.time()-t0:.0f}s", flush=True)
    for n, v in seal["artifacts"].items():
        print(f"  {n}  {v['sha256'][:32]}", flush=True)


if __name__ == "__main__":
    main()
