#!/usr/bin/env python3
"""Full-TEST2 density-availability audit + immutable balanced center->shard manifest.

Distinguishes two DIFFERENT quantities that were previously conflated:

  NEO-cell availability        the NEO 0.25-mag cell exists (what the 99.33% figure measured)
  four-density scoring coverage ALL FOUR population densities exist -> the only case in which
                               P(NEO) = rho_NEO / sum(rho) is the quantity it claims to be

Reads only. Writes the audit JSON/CSV and the shard manifest.
"""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
T2 = W / "outputs" / "test2_geometric"
MAPS = W / "outputs" / "neomod3_mag025_k150_maps_v2"
POPS = ("NEO", "MBA", "TNO", "Trojans")
NSHARDS = 32


def main():
    te = pd.read_parquet(T2 / "TEST2_TRACKLETS.parquet",
                         columns=["tracklet_uid", "population", "center_label", "magnitude_bin",
                                  "w_phys", "vlam", "vbeta", "v_out_of_range"])
    cov = pd.read_parquet(MAPS / "coverage_table.parquet")
    valid = cov[cov.valid]
    avail = {p: set(map(tuple, valid[valid.population == p][["center", "magnitude_bin"]].to_numpy()))
             for p in POPS}
    pairs = list(zip(te.center_label, te.magnitude_bin))
    have = {p: np.array([pr in avail[p] for pr in pairs]) for p in POPS}
    n_have = sum(have[p].astype(int) for p in POPS)
    inb = (te.vlam.abs() <= 5.0) & (te.vbeta.abs() <= 5.0) & (~te.v_out_of_range)
    te["n_densities"] = n_have
    te["neo_available"] = have["NEO"]
    te["all_four"] = (n_have == 4) & inb.to_numpy()
    te["in_grid"] = inb.to_numpy()

    cats = {
        "all_four_densities_available": te.all_four,
        "exactly_three_available": (n_have == 3),
        "exactly_two_available": (n_have == 2),
        "exactly_one_available": (n_have == 1),
        "neo_density_unavailable": ~te.neo_available,
        "outside_velocity_grid_or_v_range": ~te.in_grid,
    }
    wtot = te.w_phys.sum()
    rows = []
    for name, m in cats.items():
        m = np.asarray(m)
        rows.append(dict(category=name, rows=int(m.sum()),
                         raw_frac=float(m.mean()),
                         weight=float(te.w_phys[m].sum()),
                         weighted_frac=float(te.w_phys[m].sum() / wtot)))
    tab = pd.DataFrame(rows)
    print("FULL-TEST2 DENSITY AVAILABILITY"); print(tab.to_string(index=False), flush=True)

    per_pop = []
    for p in POPS:
        g = te[te.population == p]
        per_pop.append(dict(population=p, rows=len(g),
                            neo_cell_availability_raw=float(g.neo_available.mean()),
                            four_density_raw=float(g.all_four.mean()),
                            four_density_weighted=float(g.w_phys[g.all_four].sum() / g.w_phys.sum()),
                            missing_TNO=int((~np.array([pr in avail["TNO"] for pr in
                                                        zip(g.center_label, g.magnitude_bin)])).sum()),
                            missing_Trojans=int((~np.array([pr in avail["Trojans"] for pr in
                                                            zip(g.center_label, g.magnitude_bin)])).sum())))
    pp = pd.DataFrame(per_pop)
    print("\nBY POPULATION"); print(pp.to_string(index=False), flush=True)

    neo_raw = float(te.neo_available.mean())
    four_raw = float(te.all_four.mean())
    four_wt = float(te.w_phys[te.all_four].sum() / wtot)
    print(f"\n  NEO-cell availability (the earlier '99.33%' figure): raw {100*neo_raw:.2f}%")
    print(f"  FOUR-DENSITY scoring coverage                      : raw {100*four_raw:.2f}%  "
          f"weighted {100*four_wt:.2f}%", flush=True)

    # ---- immutable balanced center -> shard manifest -------------------------
    cc = te.groupby("center_label").size().sort_values(ascending=False)
    loads = np.zeros(NSHARDS, dtype=np.int64)
    assign = {}
    for cen, n in cc.items():                      # greedy longest-processing-time first
        j = int(np.argmin(loads))
        assign[cen] = j
        loads[j] += int(n)
    man = {"nshards": NSHARDS, "policy": "greedy LPT balanced by tracklet count, not center count",
           "n_centers": int(len(cc)), "n_rows": int(len(te)),
           "shard_rows": loads.tolist(),
           "balance": {"min": int(loads.min()), "max": int(loads.max()),
                       "max_over_min": float(loads.max() / max(loads.min(), 1))},
           "assignment": assign}
    body = json.dumps(man, indent=2, sort_keys=True)
    (T2 / "CENTER_SHARD_MANIFEST.json").write_text(body)
    mh = hashlib.sha256(body.encode()).hexdigest()
    (T2 / "CENTER_SHARD_MANIFEST.sha256").write_text(mh + "\n")
    print(f"\n  shard rows min {loads.min():,} max {loads.max():,} "
          f"(max/min {loads.max()/max(loads.min(),1):.2f})")
    print(f"  manifest sha256 {mh}", flush=True)

    out = {"neo_cell_availability_raw": neo_raw,
           "four_density_coverage_raw": four_raw,
           "four_density_coverage_weighted": four_wt,
           "note": ("the 99.33% figure is NEO-cell availability, NOT four-density scoring "
                    "coverage; only all-four rows yield a valid P(NEO)"),
           "categories": tab.to_dict("records"),
           "by_population": pp.to_dict("records"),
           "center_shard_manifest_sha256": mh}
    (T2 / "TEST2_COVERAGE_AUDIT.json").write_text(json.dumps(out, indent=2))
    tab.to_csv(T2 / "TEST2_COVERAGE_AUDIT.csv", index=False)
    pp.to_csv(T2 / "TEST2_COVERAGE_BY_POPULATION.csv", index=False)
    print(f"\nwrote {T2/'TEST2_COVERAGE_AUDIT.json'}", flush=True)


if __name__ == "__main__":
    main()
