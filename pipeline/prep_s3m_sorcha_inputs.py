"""
Build pure-S3M Sorcha input CSVs from the original .s3m files (no MPCORB).

The .s3m files (neomod/S3Mdata/) are cometary (COM) Keplerian elements:
    S3MID FORMAT q e i Omega argperi t_p H t_0 INDEX N_PAR MOID COMPCODE
(two `!!` header lines per file). Sorcha accepts COM orbits directly
(OrbitAuxReader cometary columns: ObjID,q,e,inc,node,argPeri,t_p_MJD_TDB,epochMJD_TDB).

Output (matches the hybrid prep convention so the run is directly comparable):
    inputs/s3m_sorcha_orbits.csv   COM format
    inputs/s3m_sorcha_phys.csv     ObjID,H_r,GS,u-r,g-r,i-r,z-r,y-r

Photometry follows scripts/prep_sorcha_inputs.py exactly: H_r = H - (Johnson V - LSST r),
colours sampled from CDS_colors.parquet with rng seed 42. The .s3m files carry no phase
slope, so GS = 0.15 (standard HG slope) is assigned for all objects.
"""
import numpy as np
import pandas as pd

BASE = "/mmfs1/gscratch/dirac/ds2004/sorcha"
S3MDIR = f"{BASE}/neomod/S3Mdata"
COLORS = pd.read_parquet(f"{BASE}/CDS_colors.parquet")
RNG = np.random.default_rng(42)

ORBIT_OUT = f"{BASE}/inputs/s3m_sorcha_orbits.csv"
PHYS_OUT = f"{BASE}/inputs/s3m_sorcha_phys.csv"

# the 17 matched S3M files, in a fixed order
FILES = (["S0.s3m"]
         + [f"S1_{i:02d}.s3m" for i in range(14)]
         + ["ST.s3m", "St5.s3m"])

COLS = ["ObjID", "FMT", "q", "e", "inc", "node", "argPeri",
        "t_p", "H", "t_0", "INDEX", "N_PAR", "MOID", "COMPCODE"]

first = True
total = 0
for fn in FILES:
    d = pd.read_csv(f"{S3MDIR}/{fn}", sep=r"\s+", comment="!",
                    header=None, names=COLS)
    n = len(d)
    orbits = pd.DataFrame({
        "ObjID": d.ObjID.astype(str),
        "FORMAT": "COM",
        "q": d.q.to_numpy(),
        "e": d.e.to_numpy(),
        "inc": d.inc.to_numpy(),
        "node": d.node.to_numpy(),
        "argPeri": d.argPeri.to_numpy(),
        "t_p_MJD_TDB": d.t_p.to_numpy(),
        "epochMJD_TDB": d.t_0.to_numpy(),
    })
    c = COLORS.iloc[RNG.integers(0, len(COLORS), size=n)].reset_index(drop=True)
    phys = pd.DataFrame({
        "ObjID": d.ObjID.astype(str),
        "H_r": d.H.to_numpy() - c["johnson V - LSST r"].to_numpy(),
        "GS": 0.15,
        "u-r": c["LSST (u-r)"].to_numpy(),
        "g-r": c["LSST (g-r)"].to_numpy(),
        "i-r": c["LSST (i-r)"].to_numpy(),
        "z-r": c["LSST (z-r)"].to_numpy(),
        "y-r": c["LSST (y-r)"].to_numpy(),
    })
    orbits.to_csv(ORBIT_OUT, mode="w" if first else "a", index=False, header=first)
    phys.to_csv(PHYS_OUT, mode="w" if first else "a", index=False, header=first)
    total += n
    print(f"{fn:12s} {n:>9,}  (cumulative {total:,})", flush=True)
    first = False

print(f"\nDONE. total S3M objects written: {total:,}")
print(f"  orbits -> {ORBIT_OUT}")
print(f"  phys   -> {PHYS_OUT}")
