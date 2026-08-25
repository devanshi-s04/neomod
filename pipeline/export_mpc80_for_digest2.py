#!/usr/bin/env python3
"""MPC-80 export adapter for the digest2 comparison branch — NOT part of the VDP scorer.

    native Rubin CSV ──> frozen VDP scorer ──> VDP scores
           │
           └──> MPC-80 export (this file) ──> unmodified digest2 ──> digest2 scores
                                                                   │
                          merge by tracklet_id  <───────────────────┘

WHY THIS IS A SEPARATE FILE
---------------------------
MPC-80 is a lossy transport: it stores magnitudes to 0.1 and RA to 0.01 s. Measured on 5,000
sealed TEST2 tracklets, round-tripping through it changes the 0.25-mag bin for **12.56%** of
tracklets and perturbs a 30-minute-baseline velocity by up to 1.3 grid cells.

That is provenance about **the representation digest2 receives**. It is NOT a VDP magnitude-bin
limitation: the frozen VDP scorer reads the native full-precision CSV, whose magnitude-bin
assignment reproduced the oracle 688,688/688,688 exactly. Keeping the export here makes it
structurally impossible for the VDP branch to ingest the rounded representation.

AUDITABILITY
------------
`--mapping` writes the tracklet_id <-> MPC-80 designation correspondence together with both
emitted lines and their sha256, so a digest2 result can always be traced back to the original
tracklet_id and to the exact bytes digest2 was given.
"""
from __future__ import annotations

import argparse, hashlib, sys
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
for _p in (W / "neomod" / "src", W / "neomod" / "adam_core_stub", W / "neomod" / "pipeline"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

# MPC-80 columns 1-12 carry the identifier and are TRUNCATED to 12 characters. With 5 leading
# spaces only 7 survive, so an 8-character key silently collides (this bit us once: T0000000,
# T0000001 and T0000002 all became "T000000" and digest2 returned nothing usable).
DESIG_PREFIX = "D"
DESIG_WIDTH = 6
MAX_TRACKLETS = 10 ** DESIG_WIDTH


def designation(i: int) -> str:
    return f"     {DESIG_PREFIX}{i:0{DESIG_WIDTH}d}"


def _sx(tok, sign=False):
    parts = tok.split()
    if len(parts) != 3:
        raise ValueError(f"bad sexagesimal field {tok!r}")
    a, b, c = float(parts[0]), float(parts[1]), float(parts[2])
    s = -1.0 if (sign and parts[0].strip().startswith("-")) else 1.0
    return s * (abs(a) + b / 60.0 + c / 3600.0)


def parse_mpc80(path: Path) -> pd.DataFrame:
    """Read back 80-column records. Used for auditing the export, never by the VDP scorer.

    Pairing rule: the key is columns 1-12 verbatim; observations sharing a key are sorted by time;
    exactly two form a tracklet. A key with any other count is reported, never truncated.
    """
    from astropy.time import Time
    rows = []
    for ln_no, raw in enumerate(Path(path).read_text().splitlines(), 1):
        if not raw.strip():
            continue
        ln = raw.ljust(80)
        key = ln[0:12]
        try:
            d = ln[15:32].strip().split()
            frac = float(d[2])
            mjd = Time(f"{int(d[0]):04d}-{int(d[1]):02d}-{int(frac):02d}T00:00:00",
                       scale="utc").mjd + (frac - int(frac))
            ra = _sx(ln[32:44].strip()) * 15.0
            dec = _sx(ln[44:56].strip(), sign=True)
            magtok = ln[65:70].strip()
            mag = float(magtok) if magtok else np.nan
            band = ln[70].strip() or ""
            obs = ln[77:80].strip()
        except Exception as e:
            raise ValueError(f"line {ln_no}: unparseable MPC-80 record ({e})") from None
        rows.append(dict(key=key, mjd=mjd, ra=ra, dec=dec, mag=mag, band=band, obs=obs,
                         line_no=ln_no))
    obs_df = pd.DataFrame(rows)
    out = []
    for key, g in obs_df.groupby("key", sort=False):
        g = g.sort_values(["mjd", "line_no"], kind="mergesort")
        out.append(dict(mpc80_designation=key.strip(), n_obs=len(g),
                        mjd0=g.mjd.iloc[0], ra0=g.ra.iloc[0], dec0=g.dec.iloc[0],
                        mag0=g.mag.iloc[0],
                        mjd1=g.mjd.iloc[-1] if len(g) > 1 else np.nan,
                        ra1=g.ra.iloc[-1] if len(g) > 1 else np.nan,
                        dec1=g.dec.iloc[-1] if len(g) > 1 else np.nan,
                        mag1=g.mag.iloc[-1] if len(g) > 1 else np.nan,
                        band=g.band.iloc[0], obs=g.obs.iloc[0],
                        first_line=int(g.line_no.iloc[0])))
    return pd.DataFrame(out).sort_values("first_line", kind="mergesort").reset_index(drop=True)


def export(csv_path: Path, mpc_path: Path, mapping_path: Path | None) -> pd.DataFrame:
    from sorcha_phase2 import format_mpc80
    df = pd.read_csv(csv_path)
    need = ["tracklet_id", "mjd0", "ra0", "dec0", "mag0_V", "mjd1", "ra1", "dec1", "mag1_V"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"input is missing required columns: {missing}")
    if df.tracklet_id.duplicated().any():
        raise ValueError("tracklet_id must be unique")
    if len(df) > MAX_TRACKLETS:
        raise ValueError(f"{len(df):,} tracklets exceeds the {MAX_TRACKLETS:,} distinct MPC-80 "
                         f"designations available at width {DESIG_WIDTH}")

    # A row with a non-finite time/position/magnitude cannot be represented in 80 columns.
    # format_mpc80 raises on NaN, which would abort the whole export because of one bad row.
    # Such rows are EXCLUDED and recorded with a reason -- never silently dropped, and never
    # emitted with a placeholder value that digest2 would treat as real.
    num = ["mjd0", "ra0", "dec0", "mag0_V", "mjd1", "ra1", "dec1", "mag1_V"]
    finite = np.ones(len(df), dtype=bool)
    why = np.array(["ok"] * len(df), dtype=object)
    for c in num:
        bad = ~np.isfinite(pd.to_numeric(df[c], errors="coerce").to_numpy(float))
        why = np.where(bad & finite, f"not_representable:{c}", why)
        finite &= ~bad
    for c, lo, hi in (("ra0", 0.0, 360.0), ("ra1", 0.0, 360.0),
                      ("dec0", -90.0, 90.0), ("dec1", -90.0, 90.0)):
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
        bad = np.isfinite(v) & ((v < lo) | (v > hi))
        why = np.where(bad & finite, f"out_of_range:{c}", why)
        finite &= ~bad
    excluded = pd.DataFrame({"tracklet_id": df.tracklet_id.to_numpy()[~finite],
                             "reason": why[~finite]})
    df = df[finite].reset_index(drop=True)
    if len(excluded):
        print(f"  excluded {len(excluded):,} row(s) not representable in MPC-80: "
              f"{excluded.reason.value_counts().to_dict()}")

    des = [designation(i) for i in range(len(df))]
    l0 = [format_mpc80(k, m, r, d, g) for k, m, r, d, g in
          zip(des, df.mjd0, df.ra0, df.dec0, df.mag0_V)]
    l1 = [format_mpc80(k, m, r, d, g) for k, m, r, d, g in
          zip(des, df.mjd1, df.ra1, df.dec1, df.mag1_V)]

    bad = [i for i, (x, y) in enumerate(zip(l0, l1)) if len(x) != 80 or len(y) != 80]
    if bad:
        raise ValueError(f"{len(bad)} emitted lines are not 80 characters")
    d0 = [x[:12] for x in l0]
    d1 = [y[:12] for y in l1]
    if len(set(d0)) != len(d0):
        raise ValueError(f"emitted MPC-80 designations collide "
                         f"({len(d0) - len(set(d0))} duplicates) -- identifier truncated to 12 columns")
    if d0 != d1:
        raise ValueError("the two detections of a tracklet carry different designations")

    Path(mpc_path).parent.mkdir(parents=True, exist_ok=True)
    Path(mpc_path).write_text("\n".join([l for pair in zip(l0, l1) for l in pair]) + "\n")

    mp = pd.DataFrame({"tracklet_id": df.tracklet_id.to_numpy(),
                       "mpc80_designation": [x.strip() for x in d0],
                       "mpc80_line0": l0, "mpc80_line1": l1,
                       "mpc80_sha256": [hashlib.sha256((x + "\n" + y).encode()).hexdigest()
                                        for x, y in zip(l0, l1)],
                       "exported": True})
    if len(excluded):
        ex = excluded.copy()
        for c in ("mpc80_designation", "mpc80_line0", "mpc80_line1", "mpc80_sha256"):
            ex[c] = None
        ex["exported"] = False
        mp = pd.concat([mp, ex.rename(columns={"reason": "exclusion_reason"})],
                       ignore_index=True)
    if mapping_path:
        Path(mapping_path).parent.mkdir(parents=True, exist_ok=True)
        mp.to_csv(mapping_path, index=False)
    return mp


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="native full-precision Rubin CSV")
    ap.add_argument("--output-mpc", required=True)
    ap.add_argument("--mapping", default=None,
                    help="CSV of tracklet_id <-> designation + emitted lines + sha256")
    ap.add_argument("--audit-roundtrip", action="store_true",
                    help="re-parse the emitted file and report the precision actually delivered "
                         "to digest2 (provenance only; does not affect VDP)")
    a = ap.parse_args()

    mp = export(Path(a.input), Path(a.output_mpc), Path(a.mapping) if a.mapping else None)
    n_exp = int(mp["exported"].sum())
    print(f"wrote {a.output_mpc}  ({2*n_exp:,} observations, {n_exp:,} tracklets exported"
          + (f", {len(mp)-n_exp:,} excluded)" if len(mp) > n_exp else ")"))
    if a.mapping:
        print(f"wrote {a.mapping}  (auditable tracklet_id <-> designation mapping)")

    if a.audit_roundtrip:
        src = pd.read_csv(a.input)
        back = parse_mpc80(Path(a.output_mpc))
        j = mp[mp.exported].merge(back, on="mpc80_designation", how="inner").merge(
            src, on="tracklet_id", how="inner", suffixes=("_mpc", "_src"))
        print("\nround-trip precision DELIVERED TO digest2 (provenance about the representation,")
        print("not a property of the VDP scorer, which reads the native CSV):")
        for nm, a_, b_, u in (("mjd0", j.mjd0_mpc, j.mjd0_src, "d"),
                              ("ra0", j.ra0_mpc, j.ra0_src, "deg"),
                              ("dec0", j.dec0_mpc, j.dec0_src, "deg"),
                              ("mag0", j.mag0, j.mag0_V, "mag")):
            d = np.abs(pd.to_numeric(a_, errors="coerce").to_numpy(float)
                       - pd.to_numeric(b_, errors="coerce").to_numpy(float))
            f = np.isfinite(d)
            if f.any():
                print(f"  {nm:5s} max {np.nanmax(d[f]):.3e} {u}   median {np.nanmedian(d[f]):.3e}")


if __name__ == "__main__":
    main()
