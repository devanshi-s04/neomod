#!/usr/bin/env python3
"""Why are non-NEOs at |v| 0.25-0.5 scored as NEO? (docs §9.5 open question)

TEST A -- truth. From the Sorcha tracklets: what fraction of objects at each |v| REALLY are NEO,
          versus what the VDP predicts. If VDP says P=1 where truth says 0.2, rho_NEO is too high.

TEST B -- is rho_NEO in that band EARNED or BLED? density_raw is already downweighted, so
          integral(rho) over a region = expected object count, and (sum of support_count) x w_abs
          = the count implied by the clones actually there. Their ratio is ~1 where density comes
          from real clones, and >>1 where it was manufactured by kNN bleed + smoothing.
          Comparing the ratio between the NEO core and the 0.25-0.5 band is robust to any
          constant normalisation error, which is what makes it decisive.
"""
import json, glob, os
import numpy as np, pandas as pd
W = "/mmfs1/gscratch/dirac/ds2004/sorcha"
OUT = f"{W}/outputs/neomod3_fullgrid"
pd.set_option("display.width", 250)

d = pd.read_parquet(f"{OUT}/fullgrid_scored_tracklets.parquet")
d["nm3"] = d.P_NEO_neomod3.fillna(0.0); d["d2"] = d.P_NEO_d2.fillna(0.0)
EDGES = [0, .15, .2, .25, .3, .35, .4, .5, .7, 1.0, 5.1]
d["vb"] = pd.cut(d.maxabs_v, EDGES)

print("="*104)
print("TEST A -- truth vs prediction as a function of speed")
g = d.groupby("vb", observed=True)
tab = pd.DataFrame({
    "n":              g.size(),
    "true_NEO":       g.is_neo.sum(),
    "TRUE_NEO_frac":  g.is_neo.mean(),
    "VDP_mean_P":     g.nm3.mean(),
    "digest2_mean_P": g.d2.mean(),
    "VDP_P>0.99_%":   g.apply(lambda x: 100*(x.nm3 > 0.99).mean(), include_groups=False),
})
tab["VDP_over"] = tab.VDP_mean_P / tab.TRUE_NEO_frac.replace(0, np.nan)
tab["d2_over"] = tab.digest2_mean_P / tab.TRUE_NEO_frac.replace(0, np.nan)
print(tab.to_string(float_format=lambda v: f"{v:,.3f}"))
print("\n  VDP_over = predicted / true. 1.0 = calibrated, >1 = over-predicting NEO.")

print("\n" + "="*104)
print("TEST A2 -- of the non-NEOs in 0.25-0.5, what ARE they, and how fast do NEOs actually move?")
band = (d.maxabs_v > 0.25) & (d.maxabs_v <= 0.5)
print(f"  tracklets in band: {band.sum():,}   of which true NEO: {d[band].is_neo.sum():,} "
      f"({100*d[band].is_neo.mean():.1f}%)")
print("  population mix in band:")
print((d[band].population.value_counts(normalize=True)*100).round(2).to_string())
for pop in ["NEO", "MBA"]:
    s = d[d.population == pop].maxabs_v
    print(f"  {pop:>4} speed: median {s.median():.3f}  90th pct {s.quantile(.9):.3f}  "
          f"frac in 0.25-0.5 {100*((s > .25) & (s <= .5)).mean():5.1f}%")

print("\n" + "="*104)
print("TEST B -- is rho_NEO in the band earned by clones, or bled into it?")
CENTERS = ["prob_maps_grid_dlon+000_lat+00.npz", "prob_maps_grid_dlon+020_lat-12.npz",
           "prob_maps_grid_dlon+050_lat-01.npz"]
w_abs = json.load(open(f"{W}/outputs/neomod3_projection_cache/cache_metadata.json"))["w_abs_objects_per_clone"]
rows = []
for cen in CENTERS:
    z = np.load(f"{W}/prob_maps_grid_neomod3_full/{cen}", allow_pickle=True)
    x, y = z["x_grid"], z["y_grid"]
    dA = float((x[1]-x[0])*(y[1]-y[0]))
    X, Y = np.meshgrid(x, y); R = np.maximum(np.abs(X), np.abs(Y))
    for b in ["mag22", "mag23"]:
        for pop in ["NEO", "MBA"]:
            rho = np.nan_to_num(np.asarray(z[f"density_raw__{pop}__{b}"], float))
            sup = np.nan_to_num(np.asarray(z[f"support_count__{pop}__{b}"], float))
            for lab, m in [("core |v|<0.25", R <= .25), ("band 0.25-0.5", (R > .25) & (R <= .5)),
                           ("fast |v|>2", R > 2.0)]:
                integral = rho[m].sum()*dA
                clones = sup[m].sum()
                implied = clones*(w_abs if pop == "NEO" else 1.0)
                rows.append(dict(center=cen.split("dlon")[1][:9], magbin=b, pop=pop, region=lab,
                                 int_rho=integral, clones=clones, implied=implied,
                                 ratio=integral/implied if implied > 0 else np.inf))
r = pd.DataFrame(rows)
piv = r.pivot_table(index=["magbin", "pop", "region"], values=["int_rho", "clones", "ratio"],
                    aggfunc="mean").reindex(
    [(m, p, g) for m in ["mag22", "mag23"] for p in ["NEO", "MBA"]
     for g in ["core |v|<0.25", "band 0.25-0.5", "fast |v|>2"]])
print(piv.to_string(float_format=lambda v: f"{v:,.4g}"))
print("\n  ratio = integral(rho) / (clones present x weight).")
print("  ~1  -> the density there is supported by clones that are actually there.")
print("  >>1 -> the density was manufactured (kNN bleed / smoothing filling empty cells).")
r.to_csv(f"{OUT}/band_diagnosis.csv", index=False)
