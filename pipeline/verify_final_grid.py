#!/usr/bin/env python3
"""Acceptance for the sealed 667-map GEN grid.

A file count is NOT sufficient: a task killed mid-write leaves a truncated .npz that Slurm reports
COMPLETED, that `ls | wc -l` counts, and that a plain resubmit SKIPS (docs new_neomod_cloning.md
§9.2). Every archive is opened and a payload array is read. Policy conformance is checked against
MAP_BUILD_SEAL.json rather than assumed.
"""
import glob, hashlib, json, os, sys
from pathlib import Path
import numpy as np
from joblib import Parallel, delayed
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
D = W/"prob_maps_grid_neomod3_GEN_final"
seal = json.loads((W/"outputs/splits/MAP_BUILD_SEAL.json").read_text())
FP = seal["FROZEN_POLICY"]
files = sorted(glob.glob(str(D/"*.npz")))
print(f"archives: {len(files)} (expect 667)   total {sum(os.path.getsize(f) for f in files)/2**30:.1f} GiB")

def check(f):
    try:
        z = np.load(f, allow_pickle=True)
        nk = len(z.files)
        _ = np.asarray(z["density_raw__NEO__mag23"]).sum()          # force payload read
        thr = float(z["smooth_support_threshold"])
        scale = bool(z["smooth_support_scale_by_clone_factor"])
        step = float(z["x_grid"][1] - z["x_grid"][0])
        lim = float(np.abs(z["x_grid"]).max())
        nu = sum(1 for k in z.files if k.startswith("density_unsmoothed__"))
        prov = json.loads(str(z["neo_provenance_json"])) if "neo_provenance_json" in z.files else {}
        return dict(f=os.path.basename(f), ok=True, nk=nk, thr=thr, scale=scale,
                    step=round(step, 6), lim=round(lim, 4), n_unsm=nu,
                    role=prov.get("source_role"), seed=prov.get("seed"),
                    s3m_dropped=prov.get("n_s3m_neo_rows_discarded"))
    except Exception as e:
        return dict(f=os.path.basename(f), ok=False, err=f"{type(e).__name__}: {e}")

res = Parallel(n_jobs=32)(delayed(check)(f) for f in files)
bad = [r for r in res if not r["ok"]]
good = [r for r in res if r["ok"]]
print(f"\n=== integrity (open + payload read) ===")
print(f"  readable {len(good)}   CORRUPT {len(bad)}")
for r in bad[:10]:
    print(f"    !! {r['f']}  {r['err'][:70]}")

def uniq(k): return sorted({r[k] for r in good})
print(f"\n=== policy conformance vs MAP_BUILD_SEAL ({seal['seal']}) ===")
checks = [
    ("smooth_support_threshold", uniq("thr"), [float(FP["smooth_support_threshold_raw_clones"])]),
    ("scale_by_clone_factor",    uniq("scale"), [FP["smooth_support_scale_by_clone_factor"]]),
    ("velocity_grid_step",       uniq("step"), [FP["velocity_grid_step"]]),
    ("velocity_grid_limit",      uniq("lim"),  [FP["velocity_grid_limit_deg_per_day"]]),
    ("neo source_role",          uniq("role"), ["NEOMOD3_GEN"]),
    ("neo seed",                 uniq("seed"), [seal["gen_provenance"]["seed"]]),
    ("density_unsmoothed arrays", uniq("n_unsm"), [32]),
    ("npz key count",            uniq("nk"), None),
]
allok = not bad
for name, got, want in checks:
    ok = (want is None) or (got == want)
    allok &= (want is None) or ok
    print(f"  {name:28s} {str(got)[:48]:50s} {'OK' if want is None or ok else 'MISMATCH want '+str(want)}")
sd = uniq("s3m_dropped")
print(f"  {'s3m NEO rows discarded':28s} {str(sd)[:48]:50s} (recorded per archive)")
print(f"\n{'='*70}\nFINAL GRID: {'PASS' if allok else 'FAIL'}")
json.dump(dict(n_archives=len(files), readable=len(good), corrupt=len(bad),
               corrupt_files=[r["f"] for r in bad],
               policy_ok=bool(allok), seal_sha256=hashlib.sha256(
                   (W/"outputs/splits/MAP_BUILD_SEAL.json").read_bytes()).hexdigest()),
          open(W/"outputs/splits/FINAL_GRID_ACCEPTANCE.json","w"), indent=2)
sys.exit(0 if allok else 1)
