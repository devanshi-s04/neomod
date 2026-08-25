#!/usr/bin/env python3
"""Step 8: row-by-row regression of Rubin Scorer v1 against the frozen TEST2 oracle.

Compares the CLI output with the authoritative TEST2_SCORED.parquet columns produced by
score_test2.score_new_vdp. Reports both geometry paths:

  derived      operational: ecliptic coordinates from observed RA/Dec (the only observable path)
  precomputed  isolation:   TEST2's stored model-frame lam/beta/vlam/vbeta
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
T2 = W/"outputs"/"test2_geometric"
SEAL = json.load(open(W/"neomod"/"seals"/"RUBIN_TRACKLET_SCORER_V1_SEAL.json"))
TOL = SEAL["numerical_tolerance"]


def compare(path: Path, tag: str, oracle: pd.DataFrame) -> dict:
    got = pd.read_csv(path)
    r = {"tag": tag, "path": str(path), "rows": int(len(got))}
    print("=" * 74); print(f"{tag}  ({len(got):,} rows)"); print("=" * 74)

    r["row_count_ok"] = len(got) == len(oracle)
    r["order_preserved"] = bool((got.tracklet_id.to_numpy() == oracle.tracklet_uid.to_numpy()).all())
    r["ids_identical"] = bool(set(got.tracklet_id) == set(oracle.tracklet_uid))
    print(f"  rows {len(got):,} vs oracle {len(oracle):,}   order preserved: {r['order_preserved']}"
          f"   id set identical: {r['ids_identical']}")

    g = got.set_index("tracklet_id").loc[oracle.tracklet_uid]

    # ---- categorical assignments must match EXACTLY -------------------------
    for name, a, b in (("center_label", g.map_center_id, oracle.center_label),
                       ("magnitude_bin", g.magnitude_bin_id, oracle.magnitude_bin),
                       ("valid", g.valid.astype(bool), oracle.new_vdp_valid.astype(bool)),
                       ("reason", g.reason, oracle.new_vdp_reason)):
        aa = a.astype(str).to_numpy(); bb = np.asarray(b).astype(str)
        same = aa == bb
        r[f"{name}_identical"] = int(same.sum()); r[f"{name}_mismatch"] = int((~same).sum())
        print(f"  {name:14s} identical {int(same.sum()):>7,}/{len(same):,}"
              f"   MISMATCH {int((~same).sum()):,}")
        if (~same).any() and name == "center_label":
            ex = pd.DataFrame({"oracle": bb[~same], "v1": aa[~same]}).value_counts().head(4)
            print("      example reassignments:"); print("      " + ex.to_string().replace("\n", "\n      "))

    # ---- numeric ------------------------------------------------------------
    for name, a, b, key in (("P_NEO", g.p_NEO, oracle.P_NEO_new, "probability"),
                            ("vlam", g.v_lambda, oracle.vlam, "velocity"),
                            ("vbeta", g.v_beta, oracle.vbeta, "velocity"),
                            ("mean_V", g.mean_V, oracle.mean_mag_V, "probability")):
        x = pd.to_numeric(a, errors="coerce").to_numpy(float)
        y = pd.to_numeric(b, errors="coerce").to_numpy(float)
        both = np.isfinite(x) & np.isfinite(y)
        nan_agree = int((np.isnan(x) == np.isnan(y)).sum())
        d = np.abs(x[both] - y[both])
        rel = d/np.maximum(np.abs(y[both]), 1e-30)
        r[f"{name}_max_abs"] = float(d.max()) if both.any() else 0.0
        r[f"{name}_max_rel"] = float(rel.max()) if both.any() else 0.0
        r[f"{name}_nan_pattern_identical"] = int(nan_agree) == len(x)
        print(f"  {name:14s} both-finite {int(both.sum()):>7,}   max|abs| {r[f'{name}_max_abs']:.3e}"
              f"   max|rel| {r[f'{name}_max_rel']:.3e}   NaN pattern identical: "
              f"{r[f'{name}_nan_pattern_identical']}")

    # ---- coverage / reasons -------------------------------------------------
    cov_v1 = float(g.valid.astype(bool).mean()); cov_or = float(oracle.new_vdp_valid.mean())
    r["coverage_v1"] = cov_v1; r["coverage_oracle"] = cov_or
    print(f"  coverage  v1 {100*cov_v1:.4f}%   oracle {100*cov_or:.4f}%")
    rc_v1 = g.reason.value_counts().to_dict(); rc_or = oracle.new_vdp_reason.value_counts().to_dict()
    r["reason_counts_v1"] = {str(k): int(v) for k, v in rc_v1.items()}
    r["reason_counts_oracle"] = {str(k): int(v) for k, v in rc_or.items()}
    r["reason_counts_identical"] = r["reason_counts_v1"] == r["reason_counts_oracle"]
    print(f"  reason-count tables identical: {r['reason_counts_identical']}")

    exact = (r["row_count_ok"] and r["order_preserved"] and r["center_label_mismatch"] == 0
             and r["magnitude_bin_mismatch"] == 0 and r["valid_mismatch"] == 0
             and r["reason_mismatch"] == 0
             and r["P_NEO_max_abs"] <= TOL["probability_abs"]
             and r["vlam_max_abs"] <= TOL["velocity_abs_deg_day"])
    r["within_sealed_tolerance"] = bool(exact)
    print(f"  WITHIN SEALED TOLERANCE: {exact}")
    return r


def main():
    o = pd.read_parquet(T2/"TEST2_SCORED.parquet",
                        columns=["tracklet_uid", "P_NEO_new", "center_label", "magnitude_bin",
                                 "new_vdp_valid", "new_vdp_reason", "vlam", "vbeta", "mean_mag_V"])
    out = []
    for tag, p in (("precomputed geometry (isolation)", T2/"TEST2_V1_SCORED_precomputed.csv"),
                   ("derived geometry (OPERATIONAL)", T2/"TEST2_V1_SCORED_derived.csv")):
        if Path(p).exists():
            out.append(compare(Path(p), tag, o))
        else:
            print(f"[skip] {p} not present")
    rep = {"seal": SEAL["interface_version"], "tolerance": TOL, "comparisons": out}
    (T2/"TEST2_V1_REGRESSION.json").write_text(json.dumps(rep, indent=2, default=str))
    print(f"\nwrote {T2/'TEST2_V1_REGRESSION.json'}")


if __name__ == "__main__":
    main()
