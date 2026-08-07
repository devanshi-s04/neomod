#!/usr/bin/env python3
"""Chunk the sealed Sorcha TEST inputs for the full production array, and write the
preproduction configuration/provenance seal.

The sealed inputs are hash-verified before anything is written. Nothing is modified.
"""
from __future__ import annotations
import hashlib, json, time
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
SEALED = W/"outputs"/"sorcha_test_inputs"
SPLITS = W/"outputs"/"splits"
PROD = W/"outputs"/"sorcha_production"
WORK = PROD/"work"
N_CHUNKS = 320
BASE_SEED = 20260806          # production base seed; distinct from the pilot's 20260805


def sha256(p, ch=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    (PROD/"production").mkdir(parents=True, exist_ok=True)
    seal = json.load(open(SPLITS/"SORCHA_TEST_INPUT_SEAL.json"))
    orb_p, phy_p = SEALED/"sorcha_test_orbits.csv", SEALED/"sorcha_test_phys.csv"
    for k, p in (("orbits", orb_p), ("phys", phy_p)):
        got = sha256(p)
        assert got == seal["outputs"][k]["sha256"], f"sealed {k} modified"
        print(f"sealed {k} VERIFIED {got[:32]}...")
    vm_p = SEALED/"sorcha_test_vminusr.parquet"
    assert vm_p.exists(), "V manifest missing -- run sorcha_test_vmanifest.py first"

    orb = pd.read_csv(orb_p, dtype={"ObjID": str})
    phy = pd.read_csv(phy_p, dtype={"ObjID": str})
    assert (orb.ObjID.values == phy.ObjID.values).all()
    man = pd.read_parquet(SPLITS/"nonneo_split_manifest.parquet")[["ObjID", "population"]]
    lab = orb[["ObjID"]].merge(man, on="ObjID", how="left")
    lab.loc[orb.ObjID.str.startswith("NM3T").values, "population"] = "NEO"
    assert lab.population.notna().all()
    lab.to_parquet(PROD/"SORCHA_TEST_OBJECTS_production.parquet", index=False)
    print(f"\nobjects {len(lab):,}: {lab.population.value_counts().to_dict()}")

    idx = np.array_split(np.arange(len(orb)), N_CHUNKS)
    chunks = []
    for i, ix in enumerate(idx):
        d = WORK/f"chunk_{i:04d}"; d.mkdir(parents=True, exist_ok=True)
        of, pf = d/"orbits.csv", d/"physical.csv"
        if not (of.exists() and pf.exists()):
            orb.iloc[ix].to_csv(of, index=False)
            phy.iloc[ix].to_csv(pf, index=False)
        chunks.append({"chunk": i, "n": int(len(ix)), "seed": BASE_SEED+i})
    print(f"wrote {N_CHUNKS} chunks of ~{len(idx[0]):,} objects")

    prereg = {
        "seal": "SORCHA_TEST_PREPRODUCTION_SEAL", "version": "1.0",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "status": "FROZEN BEFORE THE PRODUCTION RUN.",
        "inputs": {"orbits_sha256": seal["outputs"]["orbits"]["sha256"],
                   "phys_sha256": seal["outputs"]["phys"]["sha256"],
                   "vminusr_sha256": sha256(vm_p),
                   "input_seal_sha256": sha256(SPLITS/"SORCHA_TEST_INPUT_SEAL.json"),
                   "model_seal_sha256": sha256(SPLITS/"MODEL_SEAL.json")},
        "seeds": {"SORCHA_SEED_base": BASE_SEED,
                  "SORCHA_SEED_rule": "SORCHA_SEED = base + array_task_id",
                  "PPLinkingFilter": "seed=0, hardcoded inside sorcha PPMiniDifi (not tunable)",
                  "digest2": "repeatable"},
        "config": {"run": "neomod/pipeline/config/neomod3_test_detections.ini",
                   "linking": "case1 params applied downstream via the EXACT PPLinkingFilter",
                   "pointing_db": "baseline_v5.0.0_2yrs.db",
                   "native_LINKINGFILTER": "NOT run -- raw detections preserved"},
        "canonical_code": {
            "builder": "neomod/pipeline/sorcha_test_build_tracklets.py",
            "scorer": "neomod/pipeline/sorcha_test_score_tracklets.py",
            "builder_sha256": sha256(W/"neomod/pipeline/sorcha_test_build_tracklets.py"),
            "scorer_sha256": sha256(W/"neomod/pipeline/sorcha_test_score_tracklets.py"),
            "vmanifest_sha256": sha256(W/"neomod/pipeline/sorcha_test_vmanifest.py")},
        "frozen_samples": {
            "primary_all_tracklets": "all truth-identity same-night distinct-visit pairs, 0<dt<=90min, sep>=0.5\"",
            "primary_object_level": "is_first_tracklet -- earliest qualifying pair per parent",
            "secondary_case1": "sorcha_object_linked -- exact PPLinkingFilter",
            "main_analysis": "sample_main_mag245: 14<=mean_mag_V<24.5 AND both classifiers valid",
            "full_map_diagnostic": "sample_full_map: 14<=mean_mag_V<25 AND both classifiers valid"},
        "invariants": {
            "model": "MODEL_SEAL applied unchanged; no calibration/threshold/sample change after results",
            "invalid_scores": "explicit NaN, never zero",
            "out_of_domain": "rows preserved with explicit flags, never dropped"},
        "n_chunks": N_CHUNKS, "n_objects": int(len(lab)), "chunks": chunks,
    }
    p = SPLITS/"SORCHA_TEST_PREPRODUCTION_SEAL.json"
    p.write_text(json.dumps(prereg, indent=2))
    print(f"\nwrote {p}\n  sha256 {sha256(p)}")


if __name__ == "__main__":
    main()
