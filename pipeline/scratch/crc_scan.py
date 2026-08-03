import glob, os, zipfile
from joblib import Parallel, delayed
def chk(f):
    try:
        with zipfile.ZipFile(f) as z:
            bad = z.testzip()
        return (f, bad)
    except Exception as e:
        return (f, f"{type(e).__name__}: {e}")
fs = sorted(glob.glob("/mmfs1/gscratch/dirac/ds2004/sorcha/prob_maps_grid_neomod3_full/*.npz"))
r = Parallel(n_jobs=32)(delayed(chk)(f) for f in fs)
bad = [(f, b) for f, b in r if b]
print(f"CRC scan: {len(fs)} archives, {len(bad)} CORRUPT")
for f, b in bad:
    print(f"  {os.path.basename(f)}  {str(b)[:60]}")
