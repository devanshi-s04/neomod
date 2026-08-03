#!/usr/bin/env python3
"""Can we get digest2's discrimination WITHOUT throwing away VDP's calibration?

The §9.6 rank-average wins on AUC/F1 but is no longer a probability -- which discards the one thing
that makes the VDP scientifically different from digest2 (P = rho_NEO/sum rho is a real posterior;
digest2's score is a heuristic). This tests probabilistic combinations that stay probabilities.

Metrics reported for each: ranking quality (AUC, F1) AND calibration quality --
  Brier score   : mean (P - y)^2, lower better
  ECE           : expected calibration error, |predicted - observed| averaged over probability bins
  count ratio   : sum(P) / true NEO count  -- 1.000 means the score can be summed to estimate a
                  population, which is what calibration is FOR. Rank scores fail this by design.
Honest split: fit on half the sky centers, evaluate on the other half.
"""
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, precision_recall_curve, brier_score_loss
W = "/mmfs1/gscratch/dirac/ds2004/sorcha"; OUT = f"{W}/outputs/neomod3_fullgrid"
pd.set_option("display.width", 260)

d = pd.read_parquet(f"{OUT}/fullgrid_scored_tracklets.parquet")
d["nm3"] = d.P_NEO_neomod3.fillna(0.0).clip(0, 1)
d["d2"] = d.P_NEO_d2.fillna(0.0).clip(0, 1)
y = d.is_neo.to_numpy()

rng = np.random.default_rng(7)
cens = d.prob_map_file.unique(); rng.shuffle(cens)
train_c = set(cens[:len(cens)//2])
tr = d.prob_map_file.isin(train_c).to_numpy(); te = ~tr
print(f"train {tr.sum():,} tracklets / {len(train_c)} centers   test {te.sum():,} / {len(cens)-len(train_c)}")

def lg(p, eps=1e-6): 
    p = np.clip(p, eps, 1-eps); return np.log(p/(1-p))

X = np.column_stack([lg(d.nm3), lg(d.d2)])

def ece(yy, pp, nb=15):
    b = np.clip((pp*nb).astype(int), 0, nb-1); e = 0.0
    for k in range(nb):
        m = b == k
        if m.sum(): e += m.mean()*abs(pp[m].mean() - yy[m].mean())
    return e

def ev(name, s_te, is_prob):
    yy = y[te]; s = np.nan_to_num(s_te)
    p_, r_, t_ = precision_recall_curve(yy, s)
    f1 = np.divide(2*p_*r_, p_+r_, out=np.zeros_like(p_), where=(p_+r_) > 0); i = int(np.argmax(f1[:-1]))
    row = dict(score=name, AUC=roc_auc_score(yy, s), bestF1=f1[i], contam=(1-p_[i])*100,
               fastNEO_zero=100*float((s[(d.maxabs_v.to_numpy()[te] > 2)] == 0).mean()))
    if is_prob:
        row.update(Brier=brier_score_loss(yy, s), ECE=ece(yy, s), count_ratio=s.sum()/yy.sum())
    else:
        row.update(Brier=np.nan, ECE=np.nan, count_ratio=s.sum()/yy.sum())
    return row

rows = [ev("VDP alone (a real posterior)", d.nm3.to_numpy()[te], True),
        ev("digest2 alone (heuristic)", d.d2.to_numpy()[te], True)]

# rank-average, computed on the test set (transductive, as in §9.6)
rk = (pd.Series(d.nm3[te]).rank(pct=True).to_numpy() + pd.Series(d.d2[te]).rank(pct=True).to_numpy())/2
rows.append(ev("rank-average (§9.6)", rk, False))

# isotonic recalibration of each alone, fit on train
for nm, col in [("VDP isotonic-recal", "nm3"), ("digest2 isotonic-recal", "d2")]:
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1).fit(d[col].to_numpy()[tr], y[tr])
    rows.append(ev(nm, iso.predict(d[col].to_numpy()[te]), True))

# probabilistic combination: logistic regression on the two log-odds, fit on train
lr = LogisticRegression(max_iter=2000).fit(X[tr], y[tr])
p_lr = lr.predict_proba(X[te])[:, 1]
rows.append(ev("logit-combine (LR)", p_lr, True))

# LR then isotonic-recalibrated on train predictions
iso2 = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1).fit(lr.predict_proba(X[tr])[:, 1], y[tr])
rows.append(ev("logit-combine + isotonic", iso2.predict(p_lr), True))

r = pd.DataFrame(rows)
print("\n" + "="*128)
print(r.to_string(index=False, float_format=lambda v: f"{v:,.4g}"))
print("\n  count_ratio = sum(P)/true NEO count. 1.000 = can be summed for population estimates.")
print(f"  LR coefficients on [logit P_vdp, logit P_d2]: {lr.coef_[0].round(3)}  intercept {lr.intercept_[0]:.3f}")
r.to_csv(f"{OUT}/calibrated_combine.csv", index=False)

print("\n" + "="*128); print("RELIABILITY -- predicted vs observed NEO fraction (test half)")
for nm, s in [("VDP alone", d.nm3.to_numpy()[te]), ("digest2 alone", d.d2.to_numpy()[te]),
              ("logit-combine+iso", iso2.predict(p_lr))]:
    b = pd.cut(s, [0, .1, .3, .5, .7, .9, 1.001], include_lowest=True)
    t = pd.DataFrame({"pred": s, "obs": y[te]}).groupby(b, observed=True).agg(
        n=("obs", "size"), predicted=("pred", "mean"), observed=("obs", "mean"))
    print(f"\n  {nm}"); print(t.to_string(float_format=lambda v: f"{v:,.3f}"))
