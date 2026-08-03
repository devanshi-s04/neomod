#!/usr/bin/env python3
"""NEOMOD3 map acceptance: velocity reach + normalisation, vs production."""
import numpy as np, sys
NEW="prob_maps_grid_neomod3_vlim5"; PROD="prob_maps_grid_s3m_nbody"
labs=['14_16','16_18','18_20','mag20','mag21','mag22','mag23','mag24+']
for cen in ["prob_maps_grid_dlon+000_lat+00.npz","prob_maps_grid_dlon+020_lat-12.npz"]:
    zn=np.load(f"{NEW}/{cen}",allow_pickle=True); zp=np.load(f"{PROD}/{cen}",allow_pickle=True)
    xn,yn=zn['x_grid'],zn['y_grid']; Xn,Yn=np.meshgrid(xn,yn); Vn=np.maximum(np.abs(Xn),np.abs(Yn))
    xp,yp=zp['x_grid'],zp['y_grid']; Xp,Yp=np.meshgrid(xp,yp); Vp=np.maximum(np.abs(Xp),np.abs(Yp))
    celln=(xn[1]-xn[0])*(yn[1]-yn[0]); cellp=(xp[1]-xp[0])*(yp[1]-yp[0])
    print(f"\n{'='*104}\n{cen}")
    print(f"{'bin':>8} | {'NEO max|v|':>10} {'(prod)':>7} | {'NEO clones':>10} | {'int rho_NEO':>11} {'(prod)':>9} | {'NEOfrac new':>11} {'(prod)':>8}")
    for L in labs:
        scn=zn[f'support_count__NEO__{L}']; hn=scn>=1
        scp=zp[f'support_count__NEO__{L}']; hp=scp>=1
        vn=Vn[hn].max() if hn.any() else 0; vp=Vp[hp].max() if hp.any() else 0
        pops=[str(p) for p in zn['population_names']]
        In={p:float(zn[f'density_raw__{p}__{L}'].sum())*celln for p in pops}
        Ip={p:float(zp[f'density_raw__{p}__{L}'].sum())*cellp for p in pops}
        fn=In['NEO']/sum(In.values()) if sum(In.values())>0 else np.nan
        fp=Ip['NEO']/sum(Ip.values()) if sum(Ip.values())>0 else np.nan
        print(f"{L:>8} | {vn:>10.2f} {vp:>7.2f} | {scn.sum():>10,.0f} | {In['NEO']:>11,.0f} {Ip['NEO']:>9,.0f} | {fn:>11.5f} {fp:>8.5f}")
