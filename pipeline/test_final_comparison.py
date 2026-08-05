#!/usr/bin/env python3
"""Final TEST comparison: frozen VDP (variant C + Platt) vs deterministic digest2, identical rows.
Nothing is tuned. Invalid digest2 scores are ABSTENTIONS (NaN), never zero.
"""
import glob, hashlib, json, os, re, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, precision_recall_curve, roc_curve, brier_score_loss, log_loss
W=Path("/mmfs1/gscratch/dirac/ds2004/sorcha"); OUT=W/"outputs/test_final"
SEAL=json.loads((W/"outputs/splits/MODEL_SEAL.json").read_text())
THR=SEAL["operating_thresholds_from_CAL"]
te=pd.read_parquet(OUT/"TEST_SCORED.parquet")
d2=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(
    str(W/"outputs/phase2_test_neomod3/digest2_shards/*.parquet")))],ignore_index=True)
print(f"TEST rows {len(te):,}   digest2 rows {len(d2):,}")
# EXPLICIT column. Auto-detection previously matched P_NEO_vdp (the VDP score carried in the
# digest2 INPUT file) and compared VDP against itself -- identical metrics, zero paired differences.
sc="P_NEO_d2"
assert sc in d2.columns, f"{sc} missing; have {list(d2.columns)}"
_v=d2[sc].dropna()
print(f"digest2 column {sc}: n={len(_v):,} min={_v.min():.4f} max={_v.max():.4f}")
assert "P_NEO_vdp" != sc
d2=d2[["ObjID",sc]].rename(columns={sc:"P_d2_raw"}).drop_duplicates("ObjID")
m=te.merge(d2,on="ObjID",how="left")
# digest2 scores are 0-100; invalid -> NaN abstention, NEVER zero
# digest2 scores here are already on [0,1]; validate rather than assume a scale
_r=m.P_d2_raw.dropna()
_scale=100.0 if _r.max()>1.5 else 1.0
print(f"  digest2 scale detected: max={_r.max():.4f} -> dividing by {_scale}")
m["P_d2"]=np.where(np.isfinite(m.P_d2_raw)&(m.P_d2_raw>=0)&(m.P_d2_raw<=_scale),m.P_d2_raw/_scale,np.nan)
valid_v=m.technical_valid.to_numpy(); valid_d=np.isfinite(m.P_d2.to_numpy())
common=valid_v&valid_d
y=m.is_neo.to_numpy(); vmax=np.maximum(m.vlam.abs(),m.vbeta.abs()).to_numpy()
print(f"\n=== COMMON SCORABLE SET ===")
print(f"  VDP valid {int(valid_v.sum()):,}   digest2 valid {int(valid_d.sum()):,}   COMMON {int(common.sum()):,}")
ex=pd.DataFrame({"truth":m.population,"vdp":valid_v,"d2":valid_d,"common":common})
print("\n  excluded from the common set, by truth:")
print(ex.groupby("truth").apply(lambda t:pd.Series({
    "n":len(t),"vdp_only_invalid":int((~t.vdp&t.d2).sum()),"d2_only_invalid":int((t.vdp&~t.d2).sum()),
    "both_invalid":int((~t.vdp&~t.d2).sum()),"common":int(t.common.sum())}),include_groups=False).to_string())
yv=y[common]; pv=m.P_cal.to_numpy()[common]; pd2=m.P_d2.to_numpy()[common]; vv=vmax[common]
print(f"\n  common rows {len(yv):,}  NEO {int(yv.sum()):,}  prior {100*yv.mean():.4f}%")
def mets(s,nm):
    pr,rc,th=precision_recall_curve(yv,s)
    f1=np.divide(2*pr*rc,pr+rc,out=np.zeros_like(pr),where=(pr+rc)>0); i=int(np.argmax(f1[:-1]))
    return dict(classifier=nm,ROC_AUC=roc_auc_score(yv,s),
        stdpAUC_fpr01=roc_auc_score(yv,s,max_fpr=0.01),maxF1=f1[i],
        Brier=brier_score_loss(yv,np.clip(s,0,1)),logloss=log_loss(yv,np.clip(s,1e-15,1-1e-15)))
T=pd.DataFrame([mets(pv,"VDP (variant C + Platt)"),mets(pd2,"digest2 (repeatable)")])
pd.set_option("display.width",240)
print("\n=== RANKING + PROBABILITY METRICS, identical rows ===")
print(T.to_string(index=False,float_format=lambda z:f"{z:,.6g}"))
print("\n=== PAIRED DIFFERENCES (VDP - digest2), truth-stratified bootstrap B=500 ===")
rng=np.random.default_rng(0); i1=np.where(yv==1)[0]; i0=np.where(yv==0)[0]
D={k:[] for k in ["ROC_AUC","stdpAUC_fpr01","maxF1","Brier"]}
for _ in range(500):
    b=np.concatenate([rng.choice(i1,len(i1),True),rng.choice(i0,len(i0),True)])
    yy,a_,c_=yv[b],pv[b],pd2[b]
    if yy.sum()<5: continue
    D["ROC_AUC"].append(roc_auc_score(yy,a_)-roc_auc_score(yy,c_))
    D["stdpAUC_fpr01"].append(roc_auc_score(yy,a_,max_fpr=0.01)-roc_auc_score(yy,c_,max_fpr=0.01))
    def f1o(t):
        pr,rc,_=precision_recall_curve(yy,t)
        f=np.divide(2*pr*rc,pr+rc,out=np.zeros_like(pr),where=(pr+rc)>0); return f[:-1].max()
    D["maxF1"].append(f1o(a_)-f1o(c_))
    D["Brier"].append(brier_score_loss(yy,np.clip(a_,0,1))-brier_score_loss(yy,np.clip(c_,0,1)))
B=pd.DataFrame([dict(metric=k,point=T.iloc[0][k]-T.iloc[1][k],lo=np.percentile(v,2.5),
    hi=np.percentile(v,97.5),excludes_zero=(np.percentile(v,2.5)>0)or(np.percentile(v,97.5)<0))
    for k,v in D.items()])
print(B.to_string(index=False,float_format=lambda z:f"{z:,.3e}"))
print("\n=== FROZEN CAL THRESHOLDS on the common set (VDP) ===")
rows=[]
for nm,t_ in [("best-F1 (diag)",THR["best_F1_DIAGNOSTIC_ONLY"]),
              ("CAL <=5%",THR["contamination_le_5pct"]),("CAL <=10%",THR["contamination_le_10pct"])]:
    s_=pv>=t_; tp=int((s_&(yv==1)).sum()); fp=int((s_&(yv==0)).sum())
    rows.append(dict(rule=nm,threshold=t_,completeness=100*tp/max(int(yv.sum()),1),
                     contamination=100*fp/max(tp+fp,1),TP=tp,FP=fp))
print(pd.DataFrame(rows).to_string(index=False,float_format=lambda z:f"{z:,.6g}"))
print("\n=== fast NEOs (|v|>2) on the common set ===")
f_=(vv>2)&(yv==1); print(f"  n={int(f_.sum())}  VDP median {np.median(pv[f_]):.4f} "
      f"frac>=0.5 {100*np.mean(pv[f_]>=0.5):.1f}%   digest2 median {np.median(pd2[f_]):.4f} "
      f"frac>=0.5 {100*np.mean(pd2[f_]>=0.5):.1f}%")
fig,ax=plt.subplots(1,2,figsize=(13,5.2))
for s_,nm,c in [(pv,"VDP (C+Platt)","tab:green"),(pd2,"digest2","tab:orange")]:
    fpr,tpr,_=roc_curve(yv,s_); ax[0].plot(fpr,tpr,lw=1.8,color=c,label=f"{nm}  AUC={roc_auc_score(yv,s_):.4f}")
    pr,rc,_=precision_recall_curve(yv,s_); ax[1].plot(rc*100,(1-pr)*100,lw=1.8,color=c,label=nm)
ax[0].plot([0,1],[0,1],"k--",lw=.8,alpha=.5); ax[0].set_xlabel("FPR"); ax[0].set_ylabel("TPR")
ax[0].set_title(f"TEST ROC — {len(yv):,} common rows, {int(yv.sum()):,} NEOs"); ax[0].legend(); ax[0].grid(alpha=.3)
ax[1].set_xlabel("NEO completeness (%)"); ax[1].set_ylabel("contamination (%)")
ax[1].set_xlim(0,100); ax[1].set_ylim(0,100); ax[1].set_title("TEST completeness vs contamination")
ax[1].legend(); ax[1].grid(alpha=.3); fig.tight_layout()
fig.savefig(OUT/"TEST_final_comparison.png",dpi=150); plt.close(fig)
m.to_parquet(OUT/"TEST_SCORED_WITH_DIGEST2.parquet",index=False)
T.to_csv(OUT/"TEST_FINAL_METRICS.csv",index=False)
json.dump(dict(common_rows=int(common.sum()),neo=int(yv.sum()),prior=float(yv.mean()),
    metrics=T.to_dict("records"),paired=B.to_dict("records"),operating=rows),
    open(OUT/"TEST_FINAL_COMPARISON.json","w"),indent=2)
print(f"\nwrote TEST_final_comparison.png, TEST_SCORED_WITH_DIGEST2.parquet, TEST_FINAL_METRICS.csv")
