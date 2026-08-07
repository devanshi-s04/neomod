#!/usr/bin/env python3
"""CANONICAL Sorcha TEST tracklet builder -- used for BOTH pilot and production.

Supersedes build_pilot_tracklets.py, correct_pilot_tracklets.py and every generated sed copy.

Pair construction: every eligible distinct-visit pair within each (ObjID, night) is examined,
but only TRANSIENTLY. Candidates are generated one lag at a time and immediately reduced to the
running earliest per object-night, so the full Cartesian pair table is NEVER materialized
(~274M rows at production scale; up to 35,437 eligible pairs for one deep-drilling object-night).
Gates 0 < dt <= MAX_DT and separation >= MIN_SEP are applied BEFORE the earliest-pair choice.

Canonical products:
  TRACKLETS        one canonical earliest qualifying pair per (ObjID, night)
  FIRST_TRACKLET   one earliest tracklet per parent, taken from those nightly tracklets
  PAIR_DIAGNOSTICS compact per-object-night alternative-pair summary (counts, dt range only)

Geometry: spherical midpoint of the two measured positions, evaluated at that row's own
mjd_mid = (mjd0+mjd1)/2 -- never a nightly average.

Magnitudes: mag_V = mag - colour(filter) + (V-r), with (V-r) read from the sidecar manifest.

Linking: the EXACT sorcha PPLinkingFilter (drop_unlinked=False). Failure is a HARD ERROR.
"""
from __future__ import annotations
import argparse, glob, hashlib, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
C1 = json.load(open(W / "neomod/pipeline/config/neomod3_test_case1_linking.json"))["params"]
MAX_DT = float(C1["SSP_maximum_time"])                 # 0.0625 d = 90 min
MIN_SEP = float(C1["SSP_separation_threshold"])        # 0.5 arcsec
NIGHT_H = float(C1["SSP_night_start_utc"])
N_OBS = int(C1["SSP_number_observations"])
N_TRK = int(C1["SSP_number_tracklets"])
TRK_WIN = int(float(C1["SSP_track_window"]))           # MUST be int: PPMiniDifi numba kernel
SSP_EFF = float(C1["SSP_detection_efficiency"])
COLOR_COLS = {"u": "u-r", "g": "g-r", "i": "i-r", "z": "z-r", "y": "y-r"}   # r -> 0
LON_STEP, SUN_EXCL = 10.0, 40.0
DLON_LIMIT = 180.0 - SUN_EXCL
LAT_BASE = [0, 1, 2, 3, 4, 5, 8, 12, 18, 25, 35, 50]
DOM_MAIN = (14.0, 24.5)
DOM_FULL = (14.0, 25.0)
VMANIFEST = W / "outputs" / "sorcha_test_inputs" / "sorcha_test_vminusr.parquet"


def sha256(p, ch=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def build_grid():
    dl = [round(d, 6) for d in np.arange(-180.0, 180.0, LON_STEP) if abs(d) <= DLON_LIMIT + 1e-9]
    lats = sorted({float(v) for v in LAT_BASE} | {float(-v) for v in LAT_BASE})
    return np.array([(d, l) for l in lats for d in dl], dtype=float)


def map_filename(d, l):
    return f"prob_maps_grid_dlon{int(round(d)):+04d}_lat{int(round(l)):+03d}.npz"


def earliest_pair_per_object_night(det):
    """The canonical EARLIEST qualifying pair per (ObjID, night), plus compact diagnostics.

    All eligible pairs are enumerated TRANSIENTLY, one lag at a time, and immediately reduced
    to the running earliest per (ObjID, night). The complete Cartesian pair table is never
    materialized: deep-drilling fields reach 35,437 eligible pairs for a single object-night,
    which projects to ~274M rows at production scale.

    Ordering: earliest = smallest mjd0, ties broken by smallest mjd1.
    Returns (i_idx, j_idx, diagnostics) where diagnostics is one compact row per object-night.
    """
    mjd = det.mjd_utc.to_numpy(); fid = det.FieldID.to_numpy()
    ra = det.RA_deg.to_numpy(); dec = det.Dec_deg.to_numpy()
    n = len(det)
    # integer group id per detection; det is pre-sorted by (ObjID, night, mjd)
    gid = (det.ObjID.to_numpy() != np.roll(det.ObjID.to_numpy(), 1))
    gid |= (det.night.to_numpy() != np.roll(det.night.to_numpy(), 1))
    gid[0] = True
    gid = np.cumsum(gid) - 1

    best = None            # bounded: at most one row per object-night
    diag = None
    n_examined, k = 0, 1
    while k < n:
        i = np.arange(0, n - k); j = i + k
        ok = (gid[i] == gid[j]) & (fid[i] != fid[j])
        dt = mjd[j] - mjd[i]
        ok &= (dt > 0) & (dt <= MAX_DT)
        if not ok.any():
            break          # dt is monotonic in lag within a group -> no larger lag can qualify
        ii, jj = i[ok], j[ok]
        cosd = np.cos(np.deg2rad(0.5*(dec[ii]+dec[jj])))
        dra = ((ra[jj]-ra[ii]+180.0) % 360.0)-180.0
        sep = 3600.0*np.hypot(dra*cosd, dec[jj]-dec[ii])
        good = sep >= MIN_SEP                       # BOTH gates before any selection
        ii, jj = ii[good], jj[good]
        if len(ii):
            n_examined += len(ii)
            cand = pd.DataFrame({"gid": gid[ii], "mjd0": mjd[ii], "mjd1": mjd[jj],
                                 "i": ii, "j": jj})
            # reduce IMMEDIATELY: keep only the running earliest per object-night
            best = cand if best is None else pd.concat([best, cand], ignore_index=True)
            best = (best.sort_values(["gid", "mjd0", "mjd1"], kind="mergesort")
                        .drop_duplicates("gid", keep="first").reset_index(drop=True))
            d_ = pd.DataFrame({"gid": gid[ii], "dt": mjd[jj]-mjd[ii]})
            agg = d_.groupby("gid").dt.agg(n_eligible_pairs="size", dt_min_days="min",
                                           dt_max_days="max").reset_index()
            if diag is None:
                diag = agg
            else:
                diag = (pd.concat([diag, agg], ignore_index=True).groupby("gid")
                        .agg(n_eligible_pairs=("n_eligible_pairs", "sum"),
                             dt_min_days=("dt_min_days", "min"),
                             dt_max_days=("dt_max_days", "max")).reset_index())
        k += 1
    if best is None or not len(best):
        return np.array([], int), np.array([], int), pd.DataFrame()
    diag = diag.merge(best[["gid", "i"]], on="gid", how="left")
    diag["ObjID"] = det.ObjID.to_numpy()[diag.i.to_numpy()]
    diag["night"] = det.night.to_numpy()[diag.i.to_numpy()]
    diag = diag.drop(columns=["gid", "i"])
    print(f"    lag expansion stopped at k={k}; object-nights with a qualifying pair "
          f"{len(best):,}; eligible pairs examined {n_examined:,} (NEVER materialized)", flush=True)
    return best.i.to_numpy(), best.j.to_numpy(), diag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prod-dir", required=True)
    ap.add_argument("--objects", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tag", required=True)
    # Sorcha writes disjoint object sets per chunk, so tracklets AND the exact PPLinkingFilter
    # are correct per shard: no object's detections are ever split across two HDF5 files.
    # Needed because "all distinct-visit pairs" is ~274M rows at production scale (deep-drilling
    # fields give up to 35,437 pairs for a single object-night) -- far too large for one frame.
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    a_ = ap.parse_args()
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, GCRS, GeocentricTrueEcliptic, get_sun
    import astropy.units as u
    t0 = time.time()
    OUT = Path(a_.out_dir); OUT.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(Path(a_.prod_dir) / "*.h5")))
    if a_.nshards > 1:
        files = files[a_.shard::a_.nshards]
    print(f"reading {len(files)} HDF5 chunks (shard {a_.shard}/{a_.nshards})...", flush=True)
    assert files, f"no HDF5 files for shard {a_.shard}"
    det = pd.concat([pd.read_hdf(f) for f in files], ignore_index=True)
    obj = pd.read_parquet(a_.objects)
    det = det.merge(obj, on="ObjID", how="left")
    assert det.population.notna().all(), "unlabelled detection"
    det["mjd_utc"] = Time(det.fieldMJD_TAI.to_numpy(), format="mjd", scale="tai").utc.mjd
    det["night"] = np.floor(det.mjd_utc - NIGHT_H / 24.0).astype(np.int64)
    print(f"raw detections {len(det):,}  objects {det.ObjID.nunique():,}", flush=True)

    # ---- V magnitudes from the sidecar manifest (no CDS re-derivation) -----------------
    vm = pd.read_parquet(VMANIFEST)
    phys = pd.read_csv(W/"outputs/sorcha_test_inputs/sorcha_test_phys.csv",
                       usecols=["ObjID"]+list(COLOR_COLS.values()), dtype={"ObjID": str})
    det = det.merge(vm[["ObjID", "V_minus_r", "H_V"]], on="ObjID", how="left")
    assert det.V_minus_r.notna().all(), "object missing from the V manifest"
    colmap = {b: dict(zip(phys.ObjID, phys[c])) for b, c in COLOR_COLS.items()}
    col = np.zeros(len(det)); f = det.optFilter.to_numpy(); oid = det.ObjID.to_numpy()
    for b in COLOR_COLS:
        mk = f == b
        if mk.any():
            col[mk] = pd.Series(oid[mk]).map(colmap[b]).to_numpy(float)
    det["mag_V"] = det.trailedSourceMag.to_numpy(float) - col + det.V_minus_r.to_numpy(float)
    assert np.isfinite(det.mag_V).all(), "non-finite mag_V"
    print(f"mag_V built (manifest V_minus_r); shift median "
          f"{np.median(det.mag_V-det.trailedSourceMag):+.4f}", flush=True)

    # ---- ALL distinct-visit pairs ------------------------------------------------------
    det = det.sort_values(["ObjID", "night", "mjd_utc"], kind="mergesort").reset_index(drop=True)
    ii, jj, diag = earliest_pair_per_object_night(det)
    a, b = det.iloc[ii].reset_index(drop=True), det.iloc[jj].reset_index(drop=True)
    dt = b.mjd_utc.values - a.mjd_utc.values
    ra0, dec0, ra1, dec1 = a.RA_deg.values, a.Dec_deg.values, b.RA_deg.values, b.Dec_deg.values
    cosd = np.cos(np.deg2rad(0.5*(dec0+dec1)))
    dra = ((ra1-ra0+180.0) % 360.0) - 180.0
    sep = 3600.0*np.hypot(dra*cosd, dec1-dec0)
    print(f"canonical nightly tracklets (earliest qualifying pair per ObjID/night): {len(a):,}",
          flush=True)

    def unit(r, d):
        r, d = np.deg2rad(r), np.deg2rad(d)
        return np.stack([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)], axis=1)
    v = unit(ra0, dec0) + unit(ra1, dec1)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    T = pd.DataFrame({
        "ObjID": a.ObjID.values, "population": a.population.values, "night": a.night.values,
        "ra0": ra0, "dec0": dec0, "mjd0_utc": a.mjd_utc.values,
        "ra1": ra1, "dec1": dec1, "mjd1_utc": b.mjd_utc.values,
        "mag0": a.trailedSourceMag.values, "mag1": b.trailedSourceMag.values,
        "mag0_V": a.mag_V.values, "mag1_V": b.mag_V.values,
        "filter0": a.optFilter.values, "filter1": b.optFilter.values,
        "fieldid0": a.FieldID.values, "fieldid1": b.FieldID.values,
        "dt_days": dt, "sep_arcsec": sep,
        "ra_mid": np.rad2deg(np.arctan2(v[:, 1], v[:, 0])) % 360.0,
        "dec_mid": np.rad2deg(np.arcsin(np.clip(v[:, 2], -1, 1))),
        "mjd_mid": 0.5*(a.mjd_utc.values+b.mjd_utc.values),
        "mean_dra": dra*cosd/dt, "mean_ddec": (dec1-dec0)/dt,
        "snr0": a.SNR.values, "snr1": b.SNR.values,
        "astrom_sigma_deg": 0.5*(a.astrometricSigma_deg.values+b.astrometricSigma_deg.values),
        "trailing_loss_mag": 0.5*((a.trailedSourceMagTrue-a.PSFMagTrue).values
                                  + (b.trailedSourceMagTrue-b.PSFMagTrue).values),
        "H_r": a.H_r.values, "H_V": a.H_V.values, "V_minus_r": a.V_minus_r.values})
    T["mean_mag"] = 0.5*(T.mag0+T.mag1)
    T["mean_mag_V"] = 0.5*(T.mag0_V+T.mag1_V)

    print("transforming at exact per-row midpoint times...", flush=True)
    tt = Time(T.mjd_mid.to_numpy(), format="mjd", scale="utc")
    sc = SkyCoord(ra=T.ra_mid.to_numpy()*u.deg, dec=T.dec_mid.to_numpy()*u.deg,
                  pm_ra_cosdec=T.mean_dra.to_numpy()*u.deg/u.day,
                  pm_dec=T.mean_ddec.to_numpy()*u.deg/u.day, frame=GCRS(obstime=tt))
    e = sc.transform_to(GeocentricTrueEcliptic(obstime=tt))
    T["lam_deg"], T["beta_deg"] = e.lon.deg, e.lat.deg
    T["vlam"] = e.pm_lon_coslat.to(u.deg/u.day).value
    T["vbeta"] = e.pm_lat.to(u.deg/u.day).value
    antisun = (get_sun(tt).transform_to(GeocentricTrueEcliptic(obstime=tt)).lon.deg+180.0) % 360.0
    T["dlon_from_antisun_deg"] = ((T.lam_deg-antisun+180.0) % 360.0)-180.0
    grid = build_grid()
    d2 = ((T.dlon_from_antisun_deg.to_numpy()[:, None]-grid[None, :, 0])**2
          + (T.beta_deg.to_numpy()[:, None]-grid[None, :, 1])**2)
    T["prob_map_file"] = [map_filename(*grid[i]) for i in np.argmin(d2, axis=1)]
    T["n_det_per_night"] = N_OBS
    T["mean_ra"], T["mean_dec"] = T.ra_mid, T.dec_mid

    # ---- earliest qualifying pair per parent -------------------------------------------
    T = T.sort_values(["ObjID", "mjd0_utc", "mjd1_utc"], kind="mergesort").reset_index(drop=True)
    T["tracklet_index"] = T.groupby("ObjID").cumcount()
    T["is_first_tracklet"] = T.tracklet_index == 0
    T["in_domain_mag245"] = (T.mean_mag_V >= DOM_MAIN[0]) & (T.mean_mag_V < DOM_MAIN[1])
    T["in_domain_map"] = (T.mean_mag_V >= DOM_FULL[0]) & (T.mean_mag_V < DOM_FULL[1])

    # ---- EXACT PPLinkingFilter -- HARD ERROR on failure ---------------------------------
    from sorcha.modules.PPLinkingFilter import PPLinkingFilter
    obs = det[["ObjID", "FieldID", "fieldMJD_TAI", "RA_deg", "Dec_deg"]].copy()
    linked = PPLinkingFilter(obs, SSP_EFF, N_OBS, N_TRK, TRK_WIN, MIN_SEP, MAX_DT,
                             NIGHT_H, drop_unlinked=False)   # no try/except: must not fall back
    lk = linked.groupby("ObjID").object_linked.max()
    T["sorcha_object_linked"] = T.ObjID.map(lk).fillna(False).astype(bool)
    print(f"EXACT PPLinkingFilter: {int(lk.sum()):,}/{lk.size:,} objects linked", flush=True)

    # stable identifier for merging sharded scores back
    T["tracklet_uid"] = T.ObjID.astype(str) + "__" + T.night.astype(str)
    assert T.tracklet_uid.is_unique, "tracklet_uid must be one row per ObjID/night"

    tag = a_.tag
    sfx = "" if a_.nshards == 1 else f"_shard{a_.shard:03d}"
    det.to_parquet(OUT/f"SORCHA_TEST_DETECTIONS_{tag}{sfx}.parquet", index=False)
    T.to_parquet(OUT/f"SORCHA_TEST_TRACKLETS_{tag}{sfx}.parquet", index=False)
    if len(diag):
        diag.to_parquet(OUT/f"SORCHA_TEST_PAIR_DIAGNOSTICS_{tag}{sfx}.parquet", index=False)
    F = T[T.is_first_tracklet].reset_index(drop=True)
    F.to_parquet(OUT/f"SORCHA_TEST_FIRST_TRACKLET_{tag}{sfx}.parquet", index=False)
    print(f"\nDETECTIONS {len(det):,}  TRACKLETS {len(T):,}  FIRST {len(F):,}")
    print(f"  in_domain_mag245 {int(F.in_domain_mag245.sum()):,}  "
          f"in_domain_map {int(F.in_domain_map.sum()):,}  out {int((~F.in_domain_map).sum()):,}")
    print(f"  linked objects {int(F.sorcha_object_linked.sum()):,}")
    print(f"elapsed {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
