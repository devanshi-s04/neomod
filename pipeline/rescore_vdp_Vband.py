#!/usr/bin/env python3
"""Re-score Sorcha VDP in Johnson-V magnitude (Option A of SORCHA_MAG_BAND_FIX_PLAN.md).

The VDP maps' mag bins are Johnson V, but Sorcha's stored P_NEO_vdp was scored with the
observed LSST-band `mean_mag` (~0.6 mag brighter) -> objects mis-binned in magnitude. This
converts each Sorcha tracklet's magnitude to Johnson V and re-scores VDP with it, adding a
`P_NEO_vdp_Vband` column WITHOUT touching the existing `mean_mag`/`P_NEO_vdp`.

Magnitude conversion (see the fix-plan doc for the derivation):
    reduced mag (distance+phase) is filter-independent, so
      m_V = m_f - H_f + H_V = m_f - (f-r color) + (H_V - H_r)
    per tracklet (two detections f0,f1):
      mean_mag_V = mean_mag - 0.5*(color_f0 + color_f1) + (H_V - H_r)
where
    color_f  = (filter - r) LSST color, from inputs/s3m_sorcha_phys.csv  (r -> 0)
    H_r      = LSST-r absolute mag, already in the comparison parquet
    H_V      = Johnson-V absolute mag, from the S3M .s3m census (join by ObjID)

Kinematics (vlam,vbeta) are band-independent and reused as-is; only the mag bin changes,
so we re-score with ProbMapSet.score_visible(vlam, vbeta, mean_mag_V) against
prob_maps_grid_s3m using the SAME flags as Sorcha stage 3
(mask_radius_deg_per_day=inf == --no-nearest-dist-mask, support_mask_min=1).

Usage:
    python rescore_vdp_Vband.py --comparison <in.parquet> --out <out.parquet>
"""
from __future__ import annotations
import argparse, glob, os, sys, time
import numpy as np
import pandas as pd

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
NEOMD = os.path.join(WORKDIR, "neomod")
sys.path.insert(0, os.path.join(NEOMD, "src"))
sys.path.insert(0, os.path.join(NEOMD, "adam_core_stub"))

PHYS_CSV = os.path.join(WORKDIR, "inputs", "s3m_sorcha_phys.csv")
S3M_GLOB = os.path.join(NEOMD, "S3Mdata", "S*.s3m")
MAPDIR = os.path.join(WORKDIR, "prob_maps_grid_s3m")   # absolute: rescore() chdir's to NEOMD
COLOR_COLS = {"u": "u-r", "g": "g-r", "i": "i-r", "z": "z-r", "y": "y-r"}  # r -> 0

# Mag-bin layout, must match vdp.DEFAULT_MAG_BINS exactly (left-closed [lo,hi);
# outside [14,25) -> None). Verified to reproduce the stored `mag_bin_label` from
# `mean_mag` at 100.00% on all three cases.
MAG_EDGES = [14, 16, 18, 20, 21, 22, 23, 24, 25]
MAG_LABELS = ["14_16", "16_18", "18_20", "mag20", "mag21", "mag22", "mag23", "mag24+"]


def mag_bin_label(mag) -> pd.Series:
    """Apparent mag -> VDP mag-bin label, matching the pipeline's DEFAULT_MAG_BINS."""
    s = pd.cut(mag, bins=MAG_EDGES, labels=MAG_LABELS, right=False).astype(object)
    return s.where(s.notna(), None)

_HV_CACHE = None


def build_hv_lookup(objids: set[str]) -> pd.Series:
    """ObjID -> Johnson-V absolute mag H_V, read from the .s3m census (col0=OID, col8=H_V).
    Restricted to `objids` to keep it light."""
    global _HV_CACHE
    if _HV_CACHE is not None:
        return _HV_CACHE
    parts = []
    for f in sorted(glob.glob(S3M_GLOB)):
        df = pd.read_csv(f, sep=r"\s+", comment="!", header=None,
                         usecols=[0, 8], names=["OID", "H_V"], dtype={0: str})
        df["OID"] = df["OID"].str.strip()
        df = df[df["OID"].isin(objids)]
        if len(df):
            parts.append(df)
        print(f"  [.s3m] {os.path.basename(f):14s} -> {len(df):>7,} matched", flush=True)
    hv = pd.concat(parts, ignore_index=True).drop_duplicates("OID").set_index("OID")["H_V"]
    _HV_CACHE = hv
    return hv


def compute_mean_mag_V(df: pd.DataFrame) -> np.ndarray:
    """mean_mag_V = mean_mag - 0.5*(color_f0+color_f1) + (H_V - H_r)."""
    phys = pd.read_csv(PHYS_CSV, dtype={"ObjID": str})
    phys["ObjID"] = phys["ObjID"].str.strip()
    phys = phys.set_index("ObjID")

    oid = df["ObjID"].astype(str).str.strip()
    hv = build_hv_lookup(set(oid.unique()))

    def color_for(filters):
        c = np.zeros(len(df))
        f = filters.astype(str).str.strip().values
        for band, col in COLOR_COLS.items():
            m = f == band
            if m.any():
                c[m] = phys.loc[oid.values[m], col].to_numpy(float)
        # r -> 0 (already), unknown filters -> 0 (rare); flag if any
        return c

    color0 = color_for(df["filter0"])
    color1 = color_for(df["filter1"])
    H_V = hv.reindex(oid.values).to_numpy(float)
    H_r = df["H_r"].to_numpy(float)
    cVr = H_V - H_r                                   # (V - r) color, per object
    mean_mag_V = df["mean_mag"].to_numpy(float) - 0.5 * (color0 + color1) + cVr
    n_missing = int(np.isnan(mean_mag_V).sum())
    if n_missing:
        print(f"  WARNING: {n_missing:,} rows have no H_V (not in .s3m) -> mean_mag_V NaN", flush=True)
    return mean_mag_V


def rescore(df: pd.DataFrame, mean_mag_V: np.ndarray) -> np.ndarray:
    """VDP re-score with V mag: group by prob_map_file, load S3M map, score_visible.
    (score_visible returns only the per-population probabilities -- no bin labels --
    so the V mag-bin label is derived separately via mag_bin_label().)"""
    import velocity_density_pipeline_gmm as vdp
    P = np.full(len(df), np.nan)
    vl = df["vlam"].to_numpy(float); vb = df["vbeta"].to_numpy(float)
    maps = df["prob_map_file"].astype(str).values
    uniq = pd.unique(maps)
    print(f"  re-scoring {len(df):,} rows over {len(uniq)} maps", flush=True)
    for k, name in enumerate(uniq, 1):
        if "grid" not in name and "antisun" not in name:
            continue
        path = os.path.join(MAPDIR, name)
        if not os.path.exists(path):
            print(f"  MISSING map {path}", flush=True); continue
        m = maps == name
        pm = vdp.ProbMapSet.from_npz(path, support_mask_min=1,
                                     mask_radius_deg_per_day=np.inf)
        out = pm.score_visible(vl[m], vb[m], mean_mag_V[m])
        P[m] = out["NEO"]
        del pm
        if k % 50 == 0:
            print(f"    {k}/{len(uniq)} maps", flush=True)
    return P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comparison", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.chdir(WORKDIR)

    t0 = time.time()
    df = pd.read_parquet(a.comparison)
    print(f"loaded {len(df):,} rows from {a.comparison}", flush=True)

    mean_mag_V = compute_mean_mag_V(df)
    dmed = np.nanmedian(mean_mag_V - df["mean_mag"].to_numpy(float))
    print(f"mean_mag_V computed: median (V - LSST) = {dmed:+.3f} mag "
          f"(expect ~ +0.6, i.e. V fainter)", flush=True)

    os.chdir(NEOMD)
    P_V = rescore(df, mean_mag_V)
    os.chdir(WORKDIR)

    df["mean_mag_V"] = mean_mag_V
    df["mag_bin_label_Vband"] = mag_bin_label(pd.Series(mean_mag_V, index=df.index))
    df["P_NEO_vdp_Vband"] = P_V

    # self-check: the same labeller must reproduce the stored LSST-band mag_bin_label
    if "mag_bin_label" in df.columns:
        repro = mag_bin_label(df["mean_mag"])
        orig = df["mag_bin_label"].astype(object)
        agree = ((repro == orig) | (repro.isna() & orig.isna())).mean()
        print(f"  label self-check (reproduces stored mag_bin_label from mean_mag): "
              f"{100*agree:.2f}%  (expect 100.00)", flush=True)
    df.to_parquet(a.out, index=False)

    ok = np.isfinite(P_V)
    print(f"\nwrote {a.out}: {len(df):,} rows in {time.time()-t0:.0f}s", flush=True)
    print(f"  P_NEO_vdp_Vband non-null: {ok.mean()*100:.1f}%", flush=True)
    print(f"  median |P_Vband - P_orig| where both present: "
          f"{np.nanmedian(np.abs(df.P_NEO_vdp_Vband - df.P_NEO_vdp)):.4f}", flush=True)


if __name__ == "__main__":
    main()
