#!/usr/bin/env python3
"""Build the immutable Rubin Scorer v1 seal."""
from __future__ import annotations
import hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); sys.path.insert(0, str(W/"neomod"/"pipeline"))
import rubin_vdp_scorer_v1 as v1

MAPS = W/"outputs"/"neomod3_mag025_k150_maps_v2"
T2 = W/"outputs"/"test2_geometric"
SEAL = W/"neomod"/"seals"/"RUBIN_TRACKLET_SCORER_V1_SEAL.json"
SOURCES = ["src/rubin_vdp_scorer_v1.py", "pipeline/score_rubin_tracklets.py",
           "pipeline/score_test2.py", "src/velocity_density_pipeline_neomod_clone_only.py"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def git(*a):
    return subprocess.run(["git", "-C", str(W/"neomod"), *a],
                          capture_output=True, text=True).stdout.strip()


def ver(m):
    try:
        return getattr(__import__(m), "__version__", "?")
    except Exception:
        return None


mbs = json.load(open(MAPS/"MAP_BUILD_SEAL_V2.json"))
dirty = git("status", "--porcelain", "-uno")
patch = git("diff", "HEAD")
seal = {
    "interface_version": v1.INTERFACE_VERSION,
    "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "purpose": ("frozen engineering baseline of the exact classifier that produced the accepted "
                "TEST2 P_NEO_new column; NOT an operational-readiness claim"),
    "git_commit": git("rev-parse", "HEAD"),
    "git_dirty_tracked": dirty,
    "git_dirty_patch_sha256": hashlib.sha256(patch.encode()).hexdigest() if patch else None,
    "baseline_commit_for_scorer": "87f6bd82d190e946f99ab128ff8ffff380d09a7a",
    "scorer_sources": {rel: sha(W/"neomod"/rel) for rel in SOURCES},
    "map_set_id": "neomod3_mag025_k150_maps_v2",
    "map_root": str(MAPS),
    "map_build_seal": str(MAPS/"MAP_BUILD_SEAL_V2.json"),
    "map_build_seal_sha256": sha(MAPS/"MAP_BUILD_SEAL_V2.json"),
    "n_map_files": len(mbs.get("map_sha256", {})),
    "map_content_hashes_location": ("per-centre sha256 for all 667 maps are carried inside "
                                    "MAP_BUILD_SEAL_V2.json, which is itself hashed above"),
    "scientific_configuration": {
        "epoch": v1.EPOCH,
        "n_centers": 667,
        "center_grid": {"lon_step_deg": v1.LON_STEP, "dlon_limit_deg": v1.DLON_LIMIT,
                        "lat_base_deg": v1.LAT_BASE,
                        "assignment": "nearest centre in (dlon_from_antisun, ecliptic_latitude)"},
        "magnitude": {"quantity": "apparent V (HG, G=0.15)", "lo": v1.V_LO, "hi": v1.V_HI,
                      "step": v1.MAG_STEP, "n_bins": v1.N_BINS,
                      "semantics": "half-open [lo, lo+step); outside [14,25) -> invalid, never clipped",
                      "cross_bin_interpolation": False},
        "k_by_population": mbs.get("k_by_population"),
        "estimator": "Bayesian kNN, closed-form posterior mean",
        "gaussian_smoothing": False, "support_masking": False,
        "velocity_grid": {"limit_deg_day": v1.VEL_LIMIT, "step_deg_day": 0.01, "n": 1001,
                          "interpolation": "bilinear in velocity only"},
        "normalization": {"neo_effective_factor":
                          (mbs.get("frozen_config") or {}).get("neo_effective_factor"),
                          "posterior": "P(c) = rho_c / sum_all rho"},
        "validity_policy": ("a probability requires ALL FOUR population densities present and "
                            "finite; partial denominators return NaN with reason "
                            "missing_population_density:<pops>; missing densities are never "
                            "replaced by zero or epsilon"),
        "platt_calibration": False,
        "supported_magnitude_bands": list(v1.SUPPORTED_BANDS),
        "geometry_source_operational": "derived",
        "geometry_divergence_from_TEST2": (
            "TEST2 assigned centres using epoch-state-cache ecliptic coordinates (model frame). "
            "v1 derives them from observed RA/Dec, the only quantity a survey supplies. Median "
            "difference 1.8e-03 deg (~6.5 arcsec); 295 of 688,688 TEST2 rows (0.043%) change "
            "centre as a result."),
    },
    "coverage_note": {
        "four_density_coverage_raw_TEST2": 0.4396,
        "four_density_coverage_weighted_TEST2": 0.4506,
        "neo_cell_availability_TEST2": 0.9933,
        "warning": ("~43.96% of in-domain TEST2 tracklets receive a probability; the remainder "
                    "return NaN with an explicit reason. v1 freezes this limitation rather than "
                    "papering over it."),
    },
    "regression_oracle": {
        "tracklets": str(T2/"TEST2_TRACKLETS.parquet"),
        "tracklets_sha256": sha(T2/"TEST2_TRACKLETS.parquet"),
        "scored": str(T2/"TEST2_SCORED.parquet"),
        "scored_sha256": sha(T2/"TEST2_SCORED.parquet"),
        "mpc80": str(T2/"TEST2_MPC80.parquet"),
        "mpc80_sha256": sha(T2/"TEST2_MPC80.parquet"),
        "oracle_columns": ["P_NEO_new", "center_label", "magnitude_bin",
                           "new_vdp_valid", "new_vdp_reason", "vlam", "vbeta", "mean_mag_V"],
    },
    "numerical_tolerance": {
        "assignments_and_reasons": "must match EXACTLY",
        "probability_abs": 1e-12,
        "probability_rel": 1e-9,
        "velocity_abs_deg_day": 1e-7,
        "justification": ("velocity is re-derived from the two detections rather than read from "
                          "the model state; measured max |dv| = 3.4e-08 deg/day, which is 3.4e-06 "
                          "of one 0.01 deg/day grid cell"),
    },
    "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                    "numpy": ver("numpy"), "pandas": ver("pandas"), "scipy": ver("scipy"),
                    "pyarrow": ver("pyarrow"), "astropy": ver("astropy")},
}
SEAL.parent.mkdir(parents=True, exist_ok=True)
SEAL.write_text(json.dumps(seal, indent=2))
os.chmod(SEAL, 0o444)
print(f"wrote {SEAL}")
print(f"  git_commit {seal['git_commit']}")
print(f"  map_build_seal_sha256 {seal['map_build_seal_sha256'][:32]}...")
print(f"  dirty tracked files: {dirty or '(none)'}")
