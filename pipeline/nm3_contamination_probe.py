#!/usr/bin/env python3
"""Why does NEOMOD3 have ~2x digest2's contamination grid-wide, when it won at the 2 test centers?

H1 (user's): it's a SKY-DIRECTION effect -- some centers are bad and drag the pooled number down.
H2: it's a POOLING/calibration artefact -- VDP scores come from 617 DIFFERENT maps, digest2 is one
    global classifier. If per-center discrimination is fine but the score scale differs per center,
    pooling destroys AUC even with perfect per-center ranking. This is a Simpson's-paradox trap and
    it would make the pooled comparison in §9.3 unfair to VDP.
H3: it's a real population/magnitude/velocity failure mode.
"""
import re
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, precision_recall_curve
pd.set_option("display.width", 250)

d = pd.read_parquet("outputs/neomod3_fullgrid/fullgrid_scored_tracklets.parquet")
d["nm3"] = d.P_NEO_neomod3.fillna(0.0); d["d2"] = d.P_NEO_d2.fillna(0.0)
cen = d.prob_map_file.str.extract(r"dlon([+-]\d+)_lat([+-]\d+)")
d["dlon"] = cen[0].astype(int); d["lat"] = cen[1].astype(int)

def best_thr(y, s):
    p, r, t = precision_recall_curve(y, s)
    f1 = np.divide(2*p*r, p+r, out=np.zeros_like(p), where=(p+r) > 0)
    i = int(np.argmax(f1[:-1])); return float(t[i]) if i < len(t) else 1.0

y = d.is_neo.to_numpy()
T_nm3, T_d2 = best_thr(y, d.nm3.to_numpy()), best_thr(y, d.d2.to_numpy())
print(f"global best-F1 thresholds:  NEOMOD3 {T_nm3:.4f}   digest2 {T_d2:.4f}")

print("\n" + "="*78); print("H2 -- POOLED vs PER-CENTER AUC (the Simpson's-paradox check)")
rows = []
for c, g in d.groupby("prob_map_file"):
    yy = g.is_neo.to_numpy()
    if yy.sum() < 10 or (1-yy).sum() < 10: continue
    rows.append(dict(center=c, dlon=g.dlon.iloc[0], lat=g.lat.iloc[0], n=len(g), n_neo=int(yy.sum()),
                     neo_frac=yy.mean(),
                     auc_nm3=roc_auc_score(yy, g.nm3), auc_d2=roc_auc_score(yy, g.d2),
                     fp_nm3=int(((g.nm3 >= T_nm3) & (yy == 0)).sum()),
                     fp_d2=int(((g.d2 >= T_d2) & (yy == 0)).sum()),
                     tp_nm3=int(((g.nm3 >= T_nm3) & (yy == 1)).sum()),
                     tp_d2=int(((g.d2 >= T_d2) & (yy == 1)).sum())))
pc = pd.DataFrame(rows); pc["dauc"] = pc.auc_nm3 - pc.auc_d2
print(f"  centers with enough both-class data: {len(pc)}")
print(f"  POOLED     AUC:  NEOMOD3 {roc_auc_score(y, d.nm3):.4f}   digest2 {roc_auc_score(y, d.d2):.4f}")
print(f"  per-center AUC, unweighted mean:   NEOMOD3 {pc.auc_nm3.mean():.4f}   digest2 {pc.auc_d2.mean():.4f}")
print(f"  per-center AUC, median:            NEOMOD3 {pc.auc_nm3.median():.4f}   digest2 {pc.auc_d2.median():.4f}")
print(f"  per-center AUC, tracklet-weighted: "
      f"NEOMOD3 {np.average(pc.auc_nm3, weights=pc.n):.4f}   digest2 {np.average(pc.auc_d2, weights=pc.n):.4f}")
print(f"  centers where NEOMOD3 beats digest2: {(pc.dauc > 0).sum()} / {len(pc)}  "
      f"({100*(pc.dauc > 0).mean():.1f}%)")

print("\n" + "="*78); print("H1 -- is the deficit a SKY-DIRECTION effect?")
pc["absdlon"] = pc.dlon.abs(); pc["abslat"] = pc.lat.abs()
for col, bins in [("absdlon", [0,20,50,90,130,180]), ("abslat", [-1,2,8,18,35,50])]:
    g = pc.groupby(pd.cut(pc[col], bins))
    t = g.agg(centers=("center","size"), tracklets=("n","sum"), neo_frac=("neo_frac","mean"),
              auc_nm3=("auc_nm3","mean"), auc_d2=("auc_d2","mean"), dauc=("dauc","mean"))
    print(f"\n  by {col}:"); print(t.to_string(float_format=lambda v: f"{v:,.4f}"))
print("\n  10 WORST centers for NEOMOD3 vs digest2:")
print(pc.nsmallest(10, "dauc")[["center","dlon","lat","n","neo_frac","auc_nm3","auc_d2","dauc"]]
      .to_string(index=False, float_format=lambda v: f"{v:,.4f}"))
print("\n  the two §4.9 TEST centers:")
print(pc[pc.center.str.contains(r"dlon\+000_lat\+00|dlon\+020_lat-12")]
      [["center","dlon","lat","n","neo_frac","auc_nm3","auc_d2","dauc"]]
      .to_string(index=False, float_format=lambda v: f"{v:,.4f}"))
print(f"\n  neo_frac: test centers {pc[pc.center.str.contains(chr(92)+'+000_lat.00|.020_lat-12')].neo_frac.mean():.3f}"
      f"   grid mean {pc.neo_frac.mean():.3f}")

print("\n" + "="*78); print("H3 -- what ARE the false positives? (global best-F1 threshold)")
fp = d[(d.nm3 >= T_nm3) & (d.is_neo == 0)]
fpd = d[(d.d2 >= T_d2) & (d.is_neo == 0)]
print(f"  NEOMOD3 FP: {len(fp):,}    digest2 FP: {len(fpd):,}")
comp = pd.DataFrame({
    "all_nonNEO": d[d.is_neo == 0].population.value_counts(),
    "FP_neomod3": fp.population.value_counts(),
    "FP_digest2": fpd.population.value_counts()}).fillna(0).astype(int)
comp["nm3_FPrate_%"] = 100*comp.FP_neomod3/comp.all_nonNEO
comp["d2_FPrate_%"] = 100*comp.FP_digest2/comp.all_nonNEO
print(comp.to_string(float_format=lambda v: f"{v:,.2f}"))
mb = pd.cut(d.mean_mag_V, [0,20,21,22,23,24,25,99])
print("\n  by magnitude (V):")
print(pd.DataFrame({"nonNEO": d[d.is_neo==0].groupby(mb[d.is_neo==0],observed=True).size(),
                    "FP_nm3": fp.groupby(mb[fp.index],observed=True).size(),
                    "FP_d2": fpd.groupby(mb[fpd.index],observed=True).size()}).fillna(0).astype(int).to_string())
vb = pd.cut(d.maxabs_v, [0,0.25,0.5,1.0,2.0,5.1])
print("\n  by |v| max:")
print(pd.DataFrame({"nonNEO": d[d.is_neo==0].groupby(vb[d.is_neo==0],observed=True).size(),
                    "FP_nm3": fp.groupby(vb[fp.index],observed=True).size(),
                    "FP_d2": fpd.groupby(vb[fpd.index],observed=True).size()}).fillna(0).astype(int).to_string())
pc.to_csv("outputs/neomod3_fullgrid/per_center_auc.csv", index=False)
print("\nwrote outputs/neomod3_fullgrid/per_center_auc.csv")
