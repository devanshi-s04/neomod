#!/usr/bin/env python3
"""Select Platt (P), fit the final global transform on all CAL rows, emit receipt + operating table.
Bounded reporting corrections applied. No new candidates. TEST not opened. No MODEL_SEAL written.
"""
import hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logit
from sklearn.metrics import roc_auc_score, precision_recall_curve, brier_score_loss, log_loss
W=Path("/mmfs1/gscratch/dirac/ds2004/sorcha"); C=W/"outputs/calibration"; S=W/"outputs/splits"
EPS=1e-12; NB=15
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(1<<22), b""): h.update(b)
    return h.hexdigest()
z=np.load(C/"CAL_OOF_SCORES.npz"); yv=z["y"]; P_raw=z["R"]; P_oof=z["P"]; fold=z["fold"]
cmp_=json.loads((C/"CALIBRATION_COMPARISON.json").read_text())
pi_cal=cmp_["pi_cal"]; pi_ref=cmp_["pi_ref"]
s=logit(np.clip(P_raw,EPS,1-EPS))
print(f"rows {len(yv):,}  NEO {int(yv.sum()):,}  pi_CAL {pi_cal:.6e}")

# ---- correction 1: Wilson intervals ----
def wilson(k,n,zc=1.959963985):
    if n==0: return (np.nan,np.nan)
    p=k/n; d=1+zc*zc/n; c=(p+zc*zc/(2*n))/d
    h=zc/d*np.sqrt(p*(1-p)/n+zc*zc/(4*n*n)); return (max(0,c-h),min(1,c+h))
def rel(pp,nb=NB):
    e=np.linspace(0,1,nb+1); idx=np.clip(np.digitize(pp,e)-1,0,nb-1); out=[]
    for b in range(nb):
        m=idx==b
        if not m.any(): continue
        k=int(yv[m].sum()); n=int(m.sum()); lo,hi=wilson(k,n)
        out.append(dict(bin=b,n=n,pred=float(pp[m].mean()),obs=k/n,wilson_lo=lo,wilson_hi=hi))
    return pd.DataFrame(out)

# ---- final global Platt on ALL technical-valid CAL rows ----
def nll(th):
    p=np.clip(expit(np.exp(th[0])*s+th[1]),EPS,1-EPS)
    return -(yv*np.log(p)+(1-yv)*np.log(1-p)).mean()
r=minimize(nll,[0.0,0.0],method="Nelder-Mead",options=dict(maxiter=4000,xatol=1e-10,fatol=1e-12))
a=float(np.exp(r.x[0])); b=float(r.x[1])
P_fin=expit(a*s+b)
print(f"\n=== FINAL GLOBAL PLATT (all {len(yv):,} rows) ===")
print(f"  a = {a:.8f}  (= exp(alpha), constrained > 0)    b = {b:.8f}")
fp=pd.read_csv(C/"CALIBRATION_FOLD_PARAMS.csv")
print(f"  per-fold a range [{fp.a.min():.5f}, {fp.a.max():.5f}]   b range [{fp.b.min():.5f}, {fp.b.max():.5f}]")

# ---- correction 2: ranking invariance of the FINAL transform ----
rk_raw=roc_auc_score(yv,P_raw); rk_fin=roc_auc_score(yv,P_fin)
order_same=bool(np.array_equal(np.argsort(np.argsort(P_raw,kind="stable")),
                               np.argsort(np.argsort(P_fin,kind="stable"))))
print(f"\n=== ranking invariance (FINAL global transform) ===")
print(f"  ROC raw {rk_raw:.12f}  ->  final {rk_fin:.12f}   |d| = {abs(rk_fin-rk_raw):.3e}")
print(f"  rank order identical: {order_same}   (a>0 => strictly monotonic)")
print(f"  NOTE: the OOF table showed small ROC/pAUC/F1 movement because EACH FOLD used a DIFFERENT")
print(f"        fitted transform, so OOF scores are a piecewise mixture. The single global")
print(f"        transform below preserves the raw ranking EXACTLY.")

print("\n=== reliability with WILSON intervals (raw vs OOF-Platt) ===")
for nm,pp in [("raw R",P_raw),("OOF Platt",P_oof)]:
    print(f"\n  {nm}:"); print(rel(np.clip(pp,EPS,1-EPS)).to_string(index=False,float_format=lambda v:f"{v:,.4f}"))

# ---- operating points ----
def ops(pp,tag):
    pr,rc,th=precision_recall_curve(yv,pp)
    f1=np.divide(2*pr*rc,pr+rc,out=np.zeros_like(pr),where=(pr+rc)>0); i=int(np.argmax(f1[:-1]))
    rows=[dict(source=tag,rule="best-F1 (DIAGNOSTIC ONLY)",threshold=float(th[i]) if i<len(th) else 1.0,
               completeness=rc[i]*100,contamination=(1-pr[i])*100,F1=f1[i])]
    for c in (0.05,0.10):
        ok=(1-pr[:-1])<=c
        if ok.any():
            j=int(np.argmax(np.where(ok,rc[:-1],-1)))
            rows.append(dict(source=tag,rule=f"contamination <= {int(c*100)}%",
                threshold=float(th[j]),completeness=rc[j]*100,contamination=(1-pr[j])*100,
                F1=f1[j]))
        else:
            rows.append(dict(source=tag,rule=f"contamination <= {int(c*100)}%",threshold=np.nan,
                completeness=0.0,contamination=np.nan,F1=np.nan))
    return rows
OP=pd.DataFrame(ops(P_oof,"OOF Platt (honest estimate)")+ops(P_fin,"final global Platt (in-sample)"))
print("\n=== CAL-ONLY OPERATING POINTS ===")
print(OP.to_string(index=False,float_format=lambda v:f"{v:,.6g}"))
print(f"\n  prior on these rows: pi_CAL = {100*pi_cal:.4f}%  -- all contamination figures assume it")

rec={"receipt":"CALIBRATION_SELECTION","selected":"P (Platt)",
 "created_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
 "rationale":"Platt and isotonic are effectively tied on Brier (2.15015e-03 vs 2.14967e-03); Platt "
   "is marginally better on log loss; parameters are stable across folds; it retains continuous "
   "score resolution (isotonic collapses 644,013 distinct scores to 408); simpler and more portable "
   "under prior adjustment.",
 "transform":{"form":"P_cal = sigmoid(a * logit(P_raw) + b)","a":a,"b":b,
   "a_parametrisation":"a = exp(alpha), constrained > 0","fitted_on":"all technical-valid CAL rows",
   "n_rows":int(len(yv)),"n_neo":int(yv.sum()),
   "clipping":{"input":"P_raw clipped to [1e-12, 1-1e-12] before logit",
               "output":"none applied; sigmoid is already in (0,1)"}},
 "priors":{"pi_CAL":pi_cal,
   "pi_ref_overlap_weighted_density_integral":pi_ref,
   "pi_ref_caveat":"pi_ref = 0.533% is an OVERLAP-WEIGHTED DENSITY-INTEGRAL ESTIMATE: the 667 maps "
     "are 30-deg patches at 10-deg spacing, so sky is multiply counted. The 27% difference from "
     "pi_CAL must NOT be read as a definitive physical prior mismatch without a non-overlapping, "
     "area-weighted recomputation.",
   "ratio_pi_ref_over_pi_CAL":pi_ref/pi_cal},
 "prior_transfer":{"logLR_cal":"logit(P_cal) - logit(pi_CAL)",
   "P_target":"sigmoid(logLR_cal + logit(pi_target))",
   "assumption":"LABEL SHIFT: p(X|y) unchanged between CAL and target; only p(y) differs. Invalid "
     "under covariate shift, a different detection/linking pipeline, or a non-NEOMOD3 target."},
 "ranking_invariance":{"ROC_raw":rk_raw,"ROC_final":rk_fin,"abs_diff":abs(rk_fin-rk_raw),
   "rank_order_identical":order_same,
   "note":"OOF ROC/pAUC/F1 movement arises because each fold applies a DIFFERENT fitted transform; "
          "the single global transform is strictly monotonic and preserves raw ranking exactly."},
 "oof_comparison":cmp_["results"],"fold_params":cmp_["fold_params"],
 "ece_rule":"15 equal-width bins on [0,1]; per-bin WILSON 95% intervals; 500-replicate "
            "truth-stratified bootstrap CI on ECE",
 "operating_points":OP.to_dict("records"),
 "inputs":{"cal_seal_sha256":sha(S/"CAL_DATASET_SEAL.json"),
   "variant_C_receipt_sha256":sha(S/"E1_INTERPOLATION_SELECTION_RECEIPT.json"),
   "calibration_manifest_sha256":sha(W/"neomod/docs/CALIBRATION_STAGE_MANIFEST.md"),
   "map_build_seal_sha256":sha(S/"MAP_BUILD_SEAL.json"),
   "oof_scores_sha256":sha(C/"CAL_OOF_SCORES.npz")},
 "code_commit":subprocess.run(["git","-C",str(W/"neomod"),"rev-parse","HEAD"],
                              capture_output=True,text=True).stdout.strip(),
 "test_status":"SEALED -- not opened"}
p=S/"CALIBRATION_SELECTION_RECEIPT.json"; p.write_text(json.dumps(rec,indent=2,default=str))
OP.to_csv(C/"CAL_OPERATING_POINTS.csv",index=False)
print(f"\nCALIBRATION_SELECTION_RECEIPT.json sha256 = {sha(p)}")
