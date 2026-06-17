"""
hybrid_catalog_prep.py
One-time conversion: hybrid.h5 (14.4M Cartesian state vectors)
to hybrid_elements.parquet (Keplerian elements + population label).

Run on login node. No Slurm needed.
  /mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python hybrid_catalog_prep.py

Reads block0_values in 500k-row chunks to keep peak RAM under 2 GB.
"""

import os
import numpy as np
import pandas as pd
import tables

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
INPUT   = os.path.join(WORKDIR, "hybrid.h5")
OUTPUT  = os.path.join(WORKDIR, "hybrid_elements.parquet")
CHUNK   = 500_000

MU = 4 * np.pi**2 / 365.25**2   # AU^3/day^2  heliocentric


def _keplerian(r_vec, v_vec, t_0):
    """Vectorised Cartesian -> Keplerian. Returns (a, e, q, i_deg, node_deg, argperi_deg, t_p)."""
    r  = np.linalg.norm(r_vec, axis=1)
    v2 = np.sum(v_vec**2, axis=1)

    a = -MU / (2.0 * (v2 / 2.0 - MU / r))

    h_vec = np.cross(r_vec, v_vec)
    h     = np.linalg.norm(h_vec, axis=1)

    i_deg = np.degrees(np.arccos(np.clip(h_vec[:, 2] / h, -1.0, 1.0)))

    e_vec = np.cross(v_vec, h_vec) / MU - r_vec / r[:, np.newaxis]
    e     = np.linalg.norm(e_vec, axis=1)
    q     = a * (1.0 - e)

    # Longitude of ascending node: N = z_hat x h = (-hy, hx, 0)
    node_deg = np.degrees(np.arctan2(h_vec[:, 0], -h_vec[:, 1])) % 360.0

    # Argument of perihelion
    n_vec  = np.column_stack([-h_vec[:, 1], h_vec[:, 0], np.zeros(len(h_vec))])
    n_norm = np.linalg.norm(n_vec, axis=1)
    safe   = (n_norm > 1e-12) & (e > 1e-12)
    cos_w  = np.where(safe,
                      np.clip(np.sum(n_vec * e_vec, axis=1) / (n_norm * e + 1e-300), -1.0, 1.0),
                      1.0)
    arg_rad = np.arccos(cos_w)
    arg_rad = np.where(e_vec[:, 2] < 0, 2.0 * np.pi - arg_rad, arg_rad)
    argperi_deg = np.degrees(arg_rad)

    # True anomaly (quadrant from sign of r.v)
    rdotv = np.sum(r_vec * v_vec, axis=1)
    cos_f = np.clip(np.where(e > 1e-12,
                             np.sum(e_vec * r_vec, axis=1) / (e * r + 1e-300),
                             1.0), -1.0, 1.0)
    f_rad = np.arctan2(np.sign(rdotv) * np.sqrt(np.maximum(1.0 - cos_f**2, 0.0)), cos_f)

    # Eccentric anomaly -> mean anomaly -> t_p
    E_rad  = 2.0 * np.arctan2(np.sqrt(np.maximum(1.0 - e, 0.0)) * np.sin(f_rad / 2.0),
                               np.sqrt(np.maximum(1.0 + e, 0.0)) * np.cos(f_rad / 2.0))
    M_rad  = E_rad - e * np.sin(E_rad)
    n_mean = np.sqrt(MU / np.maximum(np.abs(a)**3, 1e-300))
    t_p    = t_0 - M_rad / n_mean

    return a, e, q, i_deg, node_deg, argperi_deg, t_p


def _classify(a, e, q):
    pop = np.full(len(a), "other", dtype=object)
    pop[q < 1.3] = "NEO"
    pop[(a >= 1.7) & (a < 4.1) & (q >= 1.3)] = "MBA"
    pop[(a > 4.7) & (a < 5.9) & (e < 0.3) & (q >= 1.3)] = "Trojan"
    pop[a > 30.0] = "TNO"
    return pop


# ── Main ─────────────────────────────────────────────────────────────────────

print(f"Opening {INPUT} ...", flush=True)
with tables.open_file(INPUT, "r") as f:
    # axis0 = column names (10), axis1 = row index (14.4M object IDs)
    # block0_values shape = (14444912, 10) float64
    cols    = [c.decode() if isinstance(c, bytes) else c for c in f.root.df.axis0[:]]
    n_total = f.root.df.axis1.shape[0]
    print(f"  {n_total:,} objects  columns: {cols}", flush=True)

    print("  Reading OIDs ...", flush=True)
    oids = np.array([x.decode() if isinstance(x, bytes) else str(x)
                     for x in f.root.df.axis1[:]])

    ci       = {c: i for i, c in enumerate(cols)}
    n_chunks = (n_total + CHUNK - 1) // CHUNK
    parts    = []

    for k in range(n_chunks):
        s   = k * CHUNK
        end = min(s + CHUNK, n_total)
        print(f"  chunk {k+1}/{n_chunks}  rows {s:,}-{end:,}", flush=True)

        blk   = f.root.df.block0_values[s:end, :]
        r_vec = blk[:, [ci['x'], ci['y'], ci['z']]]
        v_vec = blk[:, [ci['vx'], ci['vy'], ci['vz']]]
        t_0c  = blk[:, ci['t_0']]
        H_c   = blk[:, ci['H']]

        a, e, q, i_d, nd, apd, t_p = _keplerian(r_vec, v_vec, t_0c)
        pop = _classify(a, e, q)

        parts.append(pd.DataFrame({
            "OID":        oids[s:end],
            "q":          q,
            "e":          e,
            "i":          i_d,
            "node":       nd,
            "argperi":    apd,
            "t_p":        t_p,
            "H":          H_c,
            "t_0":        t_0c,
            "a":          a,
            "population": pop,
        }))

print("Concatenating ...", flush=True)
out = pd.concat(parts, ignore_index=True)

print(f"Writing {OUTPUT} ...", flush=True)
out.to_parquet(OUTPUT, index=False)
print(f"Done. {len(out):,} rows  {os.path.getsize(OUTPUT)/1e9:.2f} GB on disk")

print("\nPopulation breakdown:")
print(out.groupby("population").size().sort_values(ascending=False).to_string())

mask = out["OID"] == "S0000001a"
if mask.any():
    row = out[mask].iloc[0]
    print(f"\nSanity check S0000001a:  q={row.q:.5f}  e={row.e:.5f}  "
          f"i={row.i:.5f}  a={row.a:.5f}")
else:
    print("\nWARNING: S0000001a not found in OIDs")
