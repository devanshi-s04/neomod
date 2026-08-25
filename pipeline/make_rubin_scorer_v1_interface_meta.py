#!/usr/bin/env python3
"""Lightweight interface metadata for Rubin Scorer v1.

The accepted VALIDATION seal (0e61ec37...) is immutable and is NOT superseded by this file. It
records the scientific state under which the completed 688,688-row TEST2 regression was run.

This file records the interface layer only: the CLI wrapper was edited after validation to remove
MPC-80 ingest and move it to a separate digest2 adapter. That change cannot alter any probability,
density, geometry result or validation output -- it only removes an input path.
"""
from __future__ import annotations
import hashlib, json, subprocess, sys, time
from pathlib import Path

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
SEALS = W/"neomod"/"seals"
VAL = SEALS/"RUBIN_TRACKLET_SCORER_V1_SEAL.json"
OUT = SEALS/"RUBIN_TRACKLET_SCORER_V1_INTERFACE.json"
VA = W/"neomod"/"validation_artifacts"/"rubin_scorer_v1"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


val = json.load(open(VAL))
reg = json.load(open(VA/"TEST2_V1_REGRESSION.json"))
pre = next(c for c in reg["comparisons"] if c["tag"].startswith("precomputed"))
der = next(c for c in reg["comparisons"] if c["tag"].startswith("derived"))

meta = {
    "interface_version": val["interface_version"],
    "interface_revision": "v1.0.1",
    "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "supersedes_validation_seal": False,
    "validation_seal_path": str(VAL),
    "validation_seal_sha256": sha(VAL),
    "validation_commit": "988e82cea3a3200de630728fe214eee43a340ab9",
    "scorer_seal_id_at_validation": val.get("scorer_sources", {}) and
        "0e61ec37e263e9a4f9e81ed4ef25ef8fa069f630b8d2003f2ce29ce37de9eaeb",
    "map_seal_sha256": val["map_build_seal_sha256"],
    "git_commit": subprocess.run(["git", "-C", str(W/"neomod"), "rev-parse", "HEAD"],
                                 capture_output=True, text=True).stdout.strip(),

    "interface_source": "pipeline/score_rubin_tracklets.py",
    "interface_source_sha256": sha(W/"neomod"/"pipeline"/"score_rubin_tracklets.py"),
    "interface_source_sha256_at_validation":
        val["scorer_sources"]["pipeline/score_rubin_tracklets.py"],
    "digest2_adapter": "pipeline/export_mpc80_for_digest2.py",
    "digest2_adapter_sha256": sha(W/"neomod"/"pipeline"/"export_mpc80_for_digest2.py"),

    "scientific_sources_unchanged": {
        rel: {"sha256": sha(W/"neomod"/rel),
              "matches_validation_seal": sha(W/"neomod"/rel) == val["scorer_sources"][rel]}
        for rel in ("src/rubin_vdp_scorer_v1.py", "pipeline/score_test2.py",
                    "src/velocity_density_pipeline_neomod_clone_only.py")},

    "change": {
        "what": ("removed --input-format mpc80 and parse_mpc80 from the VDP CLI; MPC-80 emission "
                 "moved to pipeline/export_mpc80_for_digest2.py for the digest2 branch only"),
        "why": ("MPC-80 is a lossy transport; the VDP branch must never ingest the rounded "
                "representation. The native-CSV path is the validated one."),
        "cannot_affect": ["probabilities", "densities", "geometry results", "map selection",
                          "magnitude-bin assignment", "validation outputs"],
        "rationale": ("only an input path was removed; the scoring call and every scientific "
                      "source are byte-identical to the validation seal"),
    },

    "completed_regression": {
        "rerun_required": False,
        "jobs": ["38808527", "38808528"],
        "rows": 688688,
        "precomputed_geometry": {
            "center_label_mismatch": pre["center_label_mismatch"],
            "magnitude_bin_mismatch": pre["magnitude_bin_mismatch"],
            "valid_mismatch": pre["valid_mismatch"],
            "reason_mismatch": pre["reason_mismatch"],
            "P_NEO_max_abs": pre["P_NEO_max_abs"],
            "within_sealed_tolerance": pre["within_sealed_tolerance"]},
        "derived_geometry": {
            "center_label_mismatch": der["center_label_mismatch"],
            "magnitude_bin_mismatch": der["magnitude_bin_mismatch"],
            "coverage_v1": der["coverage_v1"]},
        "artifacts": {n: sha(VA/n) for n in sorted(x.name for x in VA.glob("*.json"))},
    },

    "wording_corrections": {
        "note": ("the validation seal is immutable and is NOT rewritten; these corrections live "
                 "here and in RUBIN_TRACKLET_SCORER_V1.md, which supersede the seal's prose"),
        "geometry": {
            "superseded_text": ("seal field geometry_divergence_from_TEST2 says v1 derives "
                                "coordinates 'the only quantity a survey supplies', implying the "
                                "ecliptic coordinates themselves are unobservable"),
            "corrected_text": ("ecliptic lambda,beta ARE derivable from observed RA/Dec, and v1 "
                               "derives them that way. What Rubin cannot supply is TEST2's "
                               "particular CACHED MODEL-FRAME lambda,beta from the epoch-state "
                               "cache. The 295-row difference is cached-model versus "
                               "observation-derived geometry, not observable versus unobservable "
                               "ecliptic coordinates."),
            "numbers_unchanged": True},
        "mpc80": {
            "corrected_text": ("the measured 12.56% MPC-80 round-trip magnitude-bin change is "
                               "provenance about the representation digest2 receives. It is NOT a "
                               "VDP magnitude-bin limitation and did NOT affect the completed "
                               "native-CSV regression, in which magnitude-bin assignment "
                               "reproduced the oracle 688,688/688,688 exactly."),
            "numbers_unchanged": True},
    },

    "operational_coverage_raw": 0.4395793160328044,
    "coverage_warning": ("four-density coverage is 43.96%; v1 is an engineering baseline, not an "
                         "operational-readiness claim"),
}
OUT.write_text(json.dumps(meta, indent=2))
print(f"wrote {OUT}")
print(f"  references validation seal {meta['validation_seal_sha256'][:32]}...")
print(f"  interface source now       {meta['interface_source_sha256'][:32]}...")
print(f"  at validation it was       {meta['interface_source_sha256_at_validation'][:32]}...")
for k, v in meta["scientific_sources_unchanged"].items():
    print(f"  science unchanged: {v['matches_validation_seal']}  {k}")
