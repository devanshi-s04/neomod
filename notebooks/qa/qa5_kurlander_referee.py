#!/usr/bin/env python3
"""QA5 / Step C — the DECISIVE experiment: L1 vs L2 on the NEOMOD3-drawn Kurlander referee.

Controlled design: the Kurlander non-NEO population is the SAME S3M our denominator uses;
only the NEO population differs (NEOMOD3 instead of S3M). So a single variable moves.

Pre-registered symmetric prediction (NEOMplan §9, before this ran):
    L1 > L2 on the S3M-drawn v5 referee      [OBSERVED: dF1 = -0.0022]
    L2 > L1 on the NEOMOD3-drawn Kurlander referee   [THIS TEST]

Knobs frozen by QA3: gamma=2, 128x64, dra_cosdec=False (our builder writes RAW alpha_dot).
F1 is prevalence-dependent -> eval set is subsampled to the v5 NEO fraction (0.293).
AUC is prevalence-independent -> also reported on the natural mix.
"""
import os, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src")); sys.path.insert(0, HERE)
import ranging_engine as re
import qa_lib as q
from qa3_knob_sweep import auc, best_f1

GAMMA, N_RHO, N_RHODOT = 2.0, 128, 64
NEO_FRAC, N_EVAL, SEED = 0.293, 120000, 0
K = os.path.join(ROOT, "outputs/kurlander")


def load_referee():
    parts = []
    for f, pop in [("tr_neo.parquet","NEO"),("tr_mba.parquet","MBA"),
                   ("tr_trojan.parquet","Trojan"),("tr_tno.parquet","TNO")]:
        p = os.path.join(K, f)
        if os.path.exists(p):
            parts.append(pd.read_parquet(p))
    df = pd.concat(parts, ignore_index=True)
    print("raw tracklets:", df.population.value_counts().to_dict(),
          f"| natural NEO frac {df.is_neo.mean():.5f}")
    rng = np.random.default_rng(SEED)
    neo = df[df.is_neo]; non = df[~df.is_neo]
    n_neo = int(N_EVAL*NEO_FRAC); n_non = N_EVAL-n_neo
    neo = neo.iloc[rng.choice(len(neo), min(n_neo,len(neo)), replace=False)]
    non = non.iloc[rng.choice(len(non), min(n_non,len(non)), replace=False)]
    ev = pd.concat([neo,non], ignore_index=True)
    # antisun elongation
    c = SkyCoord(ra=ev.mean_ra.values*u.deg, dec=ev.mean_dec.values*u.deg).barycentrictrueecliptic
    ev["ecl_lon"] = c.lon.deg
    ev["absdlon"] = np.abs(q.antisun_dlon(ev.ecl_lon.values, ev.mjd_mean.values))
    return ev.reset_index(drop=True)


def score(df, level, norm="per_H_match"):
    return re.score_tracklets(df.mean_ra.values, df.mean_dec.values, df.mean_dra.values,
                              df.mean_ddec.values, df.mean_mag.values, df.mjd_mean.values,
                              level=level, n_rho=N_RHO, n_rhodot=N_RHODOT, gamma=GAMMA,
                              chunk=250, dra_cosdec=False, neomod3_norm=norm)


if __name__ == "__main__":
    ev = load_referee()
    y = ev.is_neo.values.astype(int)
    print(f"\neval set: {len(ev):,} tracklets, NEO frac {y.mean():.4f} "
          f"(matched to v5 for F1 comparability)")

    res = {}
    for name, lvl, nrm in [("L1 (S3M prior)","L1",None),
                           ("L2 (NEOMOD3 per_H_match)","L2","per_H_match"),
                           ("L2 (NEOMOD3 absolute)","L2","absolute")]:
        t=time.time(); P = score(ev,lvl,nrm) if nrm else score(ev,lvl); res[name]=P
        print(f"  {name:28s} F1={best_f1(y,P):.4f} AUC={auc(y,P):.4f}  ({time.time()-t:.0f}s)", flush=True)

    L1,L2p,L2a = res["L1 (S3M prior)"], res["L2 (NEOMOD3 per_H_match)"], res["L2 (NEOMOD3 absolute)"]
    print("\n" + "="*72)
    print("DECISIVE TEST — NEOMOD3-drawn referee (only the NEO population differs)")
    print(f"  L2_perH - L1 : dF1={best_f1(y,L2p)-best_f1(y,L1):+.4f}  dAUC={auc(y,L2p)-auc(y,L1):+.4f}")
    print(f"  L2_abs  - L1 : dF1={best_f1(y,L2a)-best_f1(y,L1):+.4f}  dAUC={auc(y,L2a)-auc(y,L1):+.4f}")
    print(f"  [v5 S3M referee gave  dF1=-0.0022  dAUC=-0.0022]")
    flip = (best_f1(y,L2p)-best_f1(y,L1)) > 0
    print(f"\n  SYMMETRIC PREDICTION {'CONFIRMED' if flip else 'NOT confirmed'}: "
          f"L2 {'>' if flip else '<='} L1 on NEOMOD3 truth")
    print("="*72)

    # per-band
    print("\nper elongation band (dF1 = L2_perH - L1):")
    for lo,hi in [(0,20),(20,40),(40,70),(70,110),(110,180)]:
        m=(ev.absdlon>=lo)&(ev.absdlon<hi)
        if m.sum()<500: continue
        print(f"  {lo:>3}-{hi:<3} N={m.sum():>6,} NEO={y[m].mean()*100:4.1f}%  "
              f"F1: L1={best_f1(y[m],L1[m]):.4f} L2={best_f1(y[m],L2p[m]):.4f} "
              f"dF1={best_f1(y[m],L2p[m])-best_f1(y[m],L1[m]):+.4f}")

    out=ev[["population","q_au","absdlon","mean_mag","is_neo"]].copy()
    out["P_L1"],out["P_L2_perH"],out["P_L2_abs"]=L1,L2p,L2a
    out.to_parquet(os.path.join(ROOT,"outputs/qa5_kurlander_scores.parquet"),index=False)
    print("\nsaved outputs/qa5_kurlander_scores.parquet")
