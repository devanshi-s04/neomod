#!/usr/bin/env python3
"""Production physical-properties manifest: carry V_minus_r and H_V per object DIRECTLY.

The sealed phys CSV is hash-locked by SORCHA_TEST_INPUT_SEAL.json and must not be edited, so
(V-r) is carried in a SIDECAR manifest instead. It is recovered by matching each object's stored
LSST colour vector back to its CDS_colors row.

That match is asserted to be BOTH complete (every object matched) AND unique (no two CDS rows
share a colour vector, so the match is one-to-one and cannot be ambiguous).
"""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
SEALED = W / "outputs" / "sorcha_test_inputs"
KEY = ["u-r", "g-r", "i-r", "z-r", "y-r"]
OUT = SEALED / "sorcha_test_vminusr.parquet"


def sha256(p, ch=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main():
    seal = json.load(open(W / "outputs" / "splits" / "SORCHA_TEST_INPUT_SEAL.json"))
    phys_p = SEALED / "sorcha_test_phys.csv"
    assert sha256(phys_p) == seal["outputs"]["phys"]["sha256"], "sealed phys CSV modified"
    print("sealed phys CSV hash VERIFIED")

    phys = pd.read_csv(phys_p, dtype={"ObjID": str})
    cds = pd.read_parquet(W / "CDS_colors.parquet").rename(columns={
        "LSST (u-r)": "u-r", "LSST (g-r)": "g-r", "LSST (i-r)": "i-r",
        "LSST (z-r)": "z-r", "LSST (y-r)": "y-r"})[KEY + ["johnson V - LSST r"]]
    for d in (phys, cds):
        for c in KEY:
            d[c] = d[c].round(9)

    # --- UNIQUENESS: the colour vector must identify exactly one CDS row -------------
    dup = int(cds.duplicated(KEY).sum())
    print(f"CDS rows {len(cds)};  duplicated colour vectors: {dup}")
    assert dup == 0, ("CDS colour vectors are NOT unique -- the colour->(V-r) match would be "
                      f"ambiguous for {dup} rows. Refusing to build the manifest.")
    grp = cds.groupby(KEY)["johnson V - LSST r"].nunique()
    assert int(grp.max()) == 1, "a colour vector maps to more than one (V-r)"
    print("UNIQUENESS asserted: colour vector -> (V-r) is one-to-one")

    m = phys.merge(cds, on=KEY, how="left", validate="many_to_one")
    n_missing = int(m["johnson V - LSST r"].isna().sum())
    print(f"COMPLETENESS: matched {len(m)-n_missing:,}/{len(m):,};  unmatched {n_missing}")
    assert n_missing == 0, f"{n_missing} objects have no CDS colour match"

    out = pd.DataFrame({"ObjID": m.ObjID, "V_minus_r": m["johnson V - LSST r"].to_numpy(float),
                        "H_r": m.H_r.to_numpy(float)})
    out["H_V"] = out.H_r + out.V_minus_r          # H_r = H_V - (V-r)
    assert out.ObjID.is_unique and np.isfinite(out.V_minus_r).all() and np.isfinite(out.H_V).all()
    out.to_parquet(OUT, index=False)
    print(f"\nwrote {OUT}  rows {len(out):,}  sha256 {sha256(OUT)}")
    print(f"  V_minus_r: {out.V_minus_r.min():.4f} .. {out.V_minus_r.max():.4f}  "
          f"mean {out.V_minus_r.mean():.4f}")
    print(f"  H_V      : {out.H_V.min():.3f} .. {out.H_V.max():.3f}")
    neo = out.ObjID.str.startswith("NM3T")
    print(f"  NEO H_V  : {out.H_V[neo].min():.3f} .. {out.H_V[neo].max():.3f}  "
          f"(NEOMOD3 model span is 15-28)")


if __name__ == "__main__":
    main()
