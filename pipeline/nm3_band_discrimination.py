#!/usr/bin/env python3
"""The band is genuinely ~50/50 NEO/MBA and VDP is well calibrated there -- so the remaining
question is DISCRIMINATION: inside the band, can digest2 tell NEO from MBA when VDP cannot?
If yes, digest2 uses information VDP structurally lacks (arc curvature, not just mean rate), the
gap is not closable by re-weighting densities, and a combined score is the practical answer."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, precision_recall_curve
W = "/mmfs1/gscratch/dirac/ds2004/sorcha"; OUT = f"{W}/outputs/neomod3_fullgrid"
pd.set_option("display.width", 250)
d = pd.read_parquet(f"{OUT}/fullgrid_scored_tracklets.parquet")
d["nm3"] = d.P_NEO_neomod3.fillna(0.0); d["d2"] = d.P_NEO_d2.fillna(0.0)

def mets(y, s):
    y = np.asarray(y, int); s = np.nan_to_num(np.asarray(s, float))
    if y.sum() < 5 or (1-y).sum() < 5: return None
    p, r, t = precision_recall_curve(y, s)
    f1 = np.divide(2*p*r, p+r, out=np.zeros_like(p), where=(p+r) > 0); i = int(np.argmax(f1[:-1]))
    return dict(AUC=roc_auc_score(y, s), bestF1=f1[i], completeness=r[i]*100,
                contamination=(1-p[i])*100, thresh=float(t[i]) if i < len(t) else 1.0)

print("="*94); print("DISCRIMINATION WITHIN velocity slices (can it separate NEO from MBA at the SAME speed?)")
rows = []
for lab, m in [("core  |v|<0.25", d.maxabs_v <= .25), ("BAND  0.25-0.5", (d.maxabs_v > .25) & (d.maxabs_v <= .5)),
               ("fast  |v|>0.5", d.maxabs_v > .5), ("ALL", d.maxabs_v > -1)]:
    s = d[m]; y = s.is_neo.to_numpy()
    for name, col in [("NEOMOD3 VDP", "nm3"), ("digest2", "d2")]:
        r = mets(y, s[col].to_numpy())
        if r: rows.append(dict(slice=lab, n=len(s), neo_frac=y.mean(), classifier=name, **r))
print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:,.4g}"))

print("\n" + "="*94); print("COMBINED SCORES (they fail differently -> combine)")
r_nm3 = d.nm3.rank(pct=True); r_d2 = d.d2.rank(pct=True)
d["comb_prod"] = d.nm3 * d.d2
d["comb_rank"] = (r_nm3 + r_d2)/2
d["comb_geom"] = np.sqrt(d.nm3.clip(0) * d.d2.clip(0))
y = d.is_neo.to_numpy()
out = []
for name in ["nm3", "d2", "comb_prod", "comb_rank", "comb_geom"]:
    r = mets(y, d[name].to_numpy())
    band = (d.maxabs_v > .25) & (d.maxabs_v <= .5)
    fp_band = int(((d[name].to_numpy() >= r["thresh"]) & (y == 0) & band.to_numpy()).sum())
    fast = (d.maxabs_v > 2.0).to_numpy()
    out.append(dict(score=name, **r, FP_band=fp_band,
                    fastNEO_zero_pct=100*float((d[name].to_numpy()[fast] == 0).mean())))
print(pd.DataFrame(out).to_string(index=False, float_format=lambda v: f"{v:,.4g}"))
print("\n  fastNEO_zero_pct: a combined score inherits digest2's coverage of the |v|>2 NEOs --")
print("  check it does not re-break the thing the ±5 grid was built to fix.")
