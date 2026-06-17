import sys, os, numpy as np, pandas as pd

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
sys.path.insert(0, os.path.join(WORKDIR, "neomod/src"))
os.chdir(os.path.join(WORKDIR, "neomod"))

import s3m_loader

neo = s3m_loader.define_s3m("neo", verbose=False)
obj = neo.loc[neo["OID"] == "S0000001a"].iloc[0]
print("S3M elements for S0000001a:")
print(f"  q={obj.q:.8f} e={obj.e:.8f} i={obj.i:.8f}")
print(f"  node={obj.node:.8f} argperi={obj.argperi:.8f} t_p={obj.t_p:.4f}")
print(f"  t_0={obj.t_0:.1f} a={obj.a:.8f}")

os.chdir(WORKDIR)
import tables
with tables.open_file("hybrid.h5", "r") as f:
    # axis0 = columns (10), axis1 = row index (14.4M object IDs)
    cols = [c.decode() if isinstance(c, bytes) else c for c in f.root.df.axis0[:]]
    print(f"\nColumns: {cols}")
    idx = f.root.df.axis1[:]
    idx_decoded = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in idx])
    row = int(np.where(idx_decoded == "S0000001a")[0][0])
    print(f"Row position: {row}")
    # block0_values is (14444912, 10) — row-major
    data = f.root.df.block0_values[row, :]
    hobj = dict(zip(cols, data))

print("\nhybrid.h5 for S0000001a:")
print(f"  x={hobj['x']:.8f} y={hobj['y']:.8f} z={hobj['z']:.8f}")
print(f"  vx={hobj['vx']:.10f} vy={hobj['vy']:.10f} vz={hobj['vz']:.10f}")
print(f"  t_0={hobj['t_0']:.1f}")

mu = 4*np.pi**2/365.25**2
a, e = obj.a, obj.e
i_r, node_r, wp_r = np.deg2rad(obj.i), np.deg2rad(obj.node), np.deg2rad(obj.argperi)
n = np.sqrt(mu/a**3)
M = (n*(hobj['t_0'] - obj.t_p)) % (2*np.pi)
E = M
for _ in range(100):
    E -= (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
f = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
r = a*(1 - e*np.cos(E))
xo, yo = r*np.cos(f), r*np.sin(f)
vs = np.sqrt(mu/(a*(1-e**2)))
vxo, vyo = -vs*np.sin(f), vs*(e+np.cos(f))
co, so = np.cos(node_r), np.sin(node_r)
cw, sw = np.cos(wp_r),   np.sin(wp_r)
ci, si = np.cos(i_r),    np.sin(i_r)
Px, Qx = co*cw - so*sw*ci, -co*sw - so*cw*ci
Py, Qy = so*cw + co*sw*ci, -so*sw + co*cw*ci
Pz, Qz = sw*si, cw*si
xe  = Px*xo  + Qx*yo;  ye  = Py*xo  + Qy*yo;  ze  = Pz*xo  + Qz*yo
vxe = Px*vxo + Qx*vyo; vye = Py*vxo + Qy*vyo; vze = Pz*vxo + Qz*vyo

print("\nPropagated (heliocentric ecliptic J2000):")
print(f"  x={xe:.8f} y={ye:.8f} z={ze:.8f}")
print(f"  vx={vxe:.10f} vy={vye:.10f} vz={vze:.10f}")

dr = np.sqrt((xe-hobj['x'])**2 + (ye-hobj['y'])**2 + (ze-hobj['z'])**2)
dv = np.sqrt((vxe-hobj['vx'])**2 + (vye-hobj['vy'])**2 + (vze-hobj['vz'])**2)
print(f"\nResiduals: |dr|={dr:.3e} AU  ({dr*1.496e8:.0f} km)")
print(f"           |dv|={dv:.3e} AU/day  ({dv*1.496e8/86400*1000:.2f} m/s)")
print("\n[PASS] heliocentric ecliptic J2000 confirmed." if dr < 0.01 else "\n[FAIL] large residual — check frame")
