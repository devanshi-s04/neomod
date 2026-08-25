#!/usr/bin/env python3
"""Rubin tracklet scoring CLI -- interface v1 (FROZEN scientific behaviour).

    python neomod/pipeline/score_rubin_tracklets.py \
        --input  nightly_tracklets.csv \
        --output nightly_tracklets_scored.csv \
        --model-seal neomod/seals/RUBIN_TRACKLET_SCORER_V1_SEAL.json

NATIVE FULL-PRECISION CSV ONLY.

MPC-80 ingest was removed from this path deliberately. MPC-80 is a lossy transport (magnitudes to
0.1, RA to 0.01 s), and letting the VDP scorer read it would silently degrade the magnitude-bin
assignment that the native-CSV regression showed to be exact (688,688/688,688). The digest2
comparison branch gets its own adapter:

    native Rubin CSV ──> this scorer ──> VDP scores
           │
           └──> neomod/pipeline/export_mpc80_for_digest2.py ──> unmodified digest2
                                                              │
                          merge by tracklet_id  <─────────────┘

The seal is REQUIRED. The scorer never searches for or auto-selects "the latest" maps.

Exit codes
    0  completed (rows may still be individually invalid, with reasons)
    2  input or schema error
    3  seal missing, malformed, or inconsistent with the on-disk maps/science sources
    4  internal integrity assertion failed
"""
from __future__ import annotations

import argparse, hashlib, json, os, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
for _p in (W / "neomod" / "src", W / "neomod" / "adam_core_stub", W / "neomod" / "pipeline"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import rubin_vdp_scorer_v1 as v1

INTERFACE_SOURCE = "pipeline/score_rubin_tracklets.py"

OUT_COLUMNS = ["tracklet_id",
               "p_NEO", "p_MBA", "p_TNO", "p_Trojan",
               "rho_NEO", "rho_MBA", "rho_TNO", "rho_Trojan",
               "valid", "reason",
               "v_lambda", "v_beta", "mean_V",
               "map_center_id", "magnitude_bin_id", "map_id",
               "interface_version", "scorer_seal", "map_seal", "input_hash"]


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# --------------------------------------------------------------------------- seal
def load_seal(p: Path) -> dict:
    if not Path(p).exists():
        print(f"ERROR: model seal not found: {p}", file=sys.stderr); sys.exit(3)
    seal = json.load(open(p))
    for k in ("interface_version", "map_root", "map_build_seal_sha256", "scorer_sources"):
        if k not in seal:
            print(f"ERROR: seal is missing required key {k!r}", file=sys.stderr); sys.exit(3)
    if seal["interface_version"] != v1.INTERFACE_VERSION:
        print(f"ERROR: seal interface_version {seal['interface_version']!r} != "
              f"{v1.INTERFACE_VERSION!r}", file=sys.stderr); sys.exit(3)
    root = Path(seal["map_root"])
    if not root.is_dir():
        print(f"ERROR: sealed map_root does not exist: {root}", file=sys.stderr); sys.exit(3)
    mbs = root / "MAP_BUILD_SEAL_V2.json"
    if not mbs.exists():
        print(f"ERROR: sealed map root has no MAP_BUILD_SEAL_V2.json", file=sys.stderr); sys.exit(3)
    got = sha256_file(mbs)
    if got != seal["map_build_seal_sha256"]:
        print(f"ERROR: map build seal hash mismatch\n  seal says {seal['map_build_seal_sha256']}"
              f"\n  on disk  {got}", file=sys.stderr); sys.exit(3)
    # The validation seal pins BOTH the scientific path and the CLI wrapper as they stood at
    # validation. The scientific sources stay hard-pinned -- a change there invalidates the
    # completed regression. This CLI file is the interface layer; it was edited after validation
    # to remove MPC-80 ingest, which cannot alter any probability. Its current hash is recorded in
    # the interface metadata alongside a reference to the accepted validation seal, rather than
    # silently overwriting that seal.
    for rel, want in seal["scorer_sources"].items():
        f = W / "neomod" / rel
        if not f.exists():
            print(f"ERROR: scorer source missing: {rel}", file=sys.stderr); sys.exit(3)
        got = sha256_file(f)
        if got == want:
            continue
        if rel == INTERFACE_SOURCE:
            meta = Path(p).parent / "RUBIN_TRACKLET_SCORER_V1_INTERFACE.json"
            if not meta.exists():
                print(f"ERROR: {rel} differs from the validation seal and no interface metadata "
                      f"is present at {meta}", file=sys.stderr); sys.exit(3)
            im = json.load(open(meta))
            if im.get("validation_seal_sha256") != sha256_file(p):
                print("ERROR: interface metadata does not reference this validation seal",
                      file=sys.stderr); sys.exit(3)
            if im.get("interface_source_sha256") != got:
                print(f"ERROR: {rel} matches neither the validation seal nor the interface "
                      f"metadata", file=sys.stderr); sys.exit(3)
            continue
        print(f"ERROR: scientific source changed since sealing: {rel}", file=sys.stderr)
        sys.exit(3)
    return seal


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--model-seal", required=True)
    ap.add_argument("--chunk-size", type=int, default=200_000)
    ap.add_argument("--geometry-source", choices=("derived", "precomputed"), default="derived",
                    help="'derived' (default, operational) computes ecliptic coordinates from the "
                         "observed RA/Dec. 'precomputed' consumes lam_deg/beta_deg/vlam/vbeta "
                         "columns and exists ONLY to isolate the TEST2 geometry divergence.")
    ap.add_argument("--summary-json", default=None)
    a = ap.parse_args()

    t0 = time.perf_counter()
    seal = load_seal(Path(a.model_seal))
    scorer_seal_id = sha256_file(Path(a.model_seal))

    t_parse0 = time.perf_counter()
    try:
        df = pd.read_csv(a.input)
        df["_parse_reason"] = "ok"
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr); sys.exit(2)
    t_parse = time.perf_counter() - t_parse0

    if "tracklet_id" not in df.columns:
        print("ERROR: input has no tracklet_id column", file=sys.stderr); sys.exit(2)
    df = df.reset_index(drop=True)
    df["_row_order"] = np.arange(len(df))

    ih = v1.input_hash(df) if all(c in df.columns for c in v1.REQUIRED_COLUMNS) else "n/a"

    stats: dict = {}
    try:
        scored = v1.score(df, map_root=Path(seal["map_root"]),
                          geometry_source=a.geometry_source,
                          chunk_size=a.chunk_size, stats=stats)
    except AssertionError as e:
        print(f"INTEGRITY FAILURE: {e}", file=sys.stderr); sys.exit(4)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr); sys.exit(2)

    # parse-stage reasons win over scoring reasons (the row never reached the maps)
    pr = df["_parse_reason"].to_numpy()
    bad = pr != "ok"
    if bad.any():
        scored.loc[bad, "valid"] = False
        scored.loc[bad, "reason"] = pr[bad]
        for c in v1.POPS:
            scored.loc[bad, f"p_{c}"] = np.nan
            scored.loc[bad, f"rho_{c}"] = np.nan

    out = pd.DataFrame(index=scored.index)
    out["tracklet_id"] = scored["tracklet_id"].to_numpy()
    for c, name in zip(v1.POPS, ("NEO", "MBA", "TNO", "Trojan")):
        out[f"p_{name}"] = scored[f"p_{c}"].to_numpy()
    for c, name in zip(v1.POPS, ("NEO", "MBA", "TNO", "Trojan")):
        out[f"rho_{name}"] = scored[f"rho_{c}"].to_numpy()
    out["valid"] = scored["valid"].to_numpy()
    out["reason"] = scored["reason"].to_numpy()
    out["v_lambda"] = scored["vlam"].to_numpy()
    out["v_beta"] = scored["vbeta"].to_numpy()
    out["mean_V"] = scored["mean_V"].to_numpy()
    out["map_center_id"] = scored["map_center_id"].to_numpy()
    out["magnitude_bin_id"] = scored["magnitude_bin_id"].to_numpy()
    out["map_id"] = [None if c is None or (isinstance(c, float) and np.isnan(c))
                     else f"mag025_k150_{c}.npz" for c in scored["map_center_id"]]
    out["interface_version"] = v1.INTERFACE_VERSION
    out["scorer_seal"] = scorer_seal_id
    out["map_seal"] = seal["map_build_seal_sha256"]
    out["input_hash"] = ih

    out = out.loc[scored["_row_order"].sort_values().index]          # restore input order
    out = out[OUT_COLUMNS]
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.output, index=False)

    wall = time.perf_counter() - t0
    rc = pd.Series(out["reason"]).value_counts()
    summary = {"input": str(a.input), "output": str(a.output),
               "input_format": "csv", "geometry_source": a.geometry_source,
               "chunk_size": a.chunk_size,
               "interface_version": v1.INTERFACE_VERSION,
               "scorer_seal": scorer_seal_id, "map_seal": seal["map_build_seal_sha256"],
               "map_root": seal["map_root"], "input_hash": ih,
               "rows_in": int(len(df)), "rows_out": int(len(out)),
               "rows_valid": int(out["valid"].sum()),
               "coverage_raw": float(out["valid"].mean()) if len(out) else float("nan"),
               "reason_counts": {str(k): int(v) for k, v in rc.items()},
               "wall_seconds": round(wall, 3),
               "rows_per_second": round(len(out) / max(wall, 1e-9), 1),
               "parse_seconds": round(t_parse, 3),
               **stats}
    try:
        import resource
        summary["peak_rss_gb"] = round(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 3)
    except Exception:
        pass
    print(json.dumps(summary, indent=2))
    if a.summary_json:
        Path(a.summary_json).write_text(json.dumps(summary, indent=2))
    sys.exit(0)


if __name__ == "__main__":
    main()
