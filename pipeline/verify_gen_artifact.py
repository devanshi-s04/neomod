#!/usr/bin/env python3
"""Verify every frozen GEN artifact byte and issue a receipt tied to the manifest SHA.

Run ONCE before E0. Map jobs then require the receipt (they check the receipt's manifest_sha256
matches the live manifest) instead of re-hashing 24 GiB per task.
"""
import hashlib, json, sys, time
from pathlib import Path
from joblib import Parallel, delayed
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
MAN = W/"outputs/splits/GEN_MANIFEST.json"
ART = ("datacube", "cache_meta", "monolithic", "shards", "by_pixel")

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()

man = json.loads(MAN.read_text())
man_sha = sha(MAN)
print(f"manifest       : {MAN}")
print(f"manifest sha256: {man_sha}")
print(f"artifact sha256: {man.get('frozen_artifact_sha256')}")
print(f"git commit     : {man.get('git_commit','')[:12]}  dirty={man.get('git_dirty')}")
if man.get("git_dirty"):
    sys.exit("REFUSING: manifest was built from a dirty worktree")

t0 = time.time(); bad = []; n = 0
for g in ART:
    files = man["groups"][g]["files"]
    res = Parallel(n_jobs=32)(delayed(sha)(W/rel) for rel in files)
    for (rel, want), got in zip(files.items(), res):
        n += 1
        if got != want:
            bad.append((rel, want, got))
    print(f"  {g:12s} {len(files):>5} files verified  ({time.time()-t0:.0f}s)", flush=True)

# recompute the artifact-only digest independently
art = {k: v for g in ART for k, v in man["groups"][g]["files"].items()}
art_sha = hashlib.sha256(json.dumps(art, sort_keys=True).encode()).hexdigest()
art_ok = art_sha == man.get("frozen_artifact_sha256")
print(f"\nfiles verified          : {n:,}")
print(f"hash mismatches         : {len(bad)}")
print(f"artifact digest matches : {art_ok}")
for rel, want, got in bad[:10]:
    print(f"   MISMATCH {rel}\n     want {want[:32]}\n     got  {got[:32]}")

ok = (not bad) and art_ok
receipt = {
    "receipt": "GEN_ARTIFACT_VERIFICATION",
    "result": "PASS" if ok else "FAIL",
    "manifest_path": str(MAN.relative_to(W)),
    "manifest_sha256": man_sha,
    "frozen_artifact_sha256": man.get("frozen_artifact_sha256"),
    "recomputed_artifact_sha256": art_sha,
    "files_verified": n,
    "by_pixel_partitions_verified": man["groups"]["by_pixel"]["n_files"],
    "mismatches": len(bad),
    "git_commit": man.get("git_commit"), "git_dirty": man.get("git_dirty"),
    "verified_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "note": "Map jobs must require this receipt and assert receipt.manifest_sha256 == sha256(live "
            "GEN_MANIFEST.json). Any artifact edit changes the manifest and invalidates the receipt.",
}
rp = W/"outputs/splits/GEN_VERIFICATION_RECEIPT.json"
rp.write_text(json.dumps(receipt, indent=2, sort_keys=True))
print(f"\nreceipt -> {rp}")
print(f"  receipt sha256 = {sha(rp)}")
print(f"  RESULT: {receipt['result']}")
sys.exit(0 if ok else 1)
