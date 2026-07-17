#!/usr/bin/env python3
"""QA4 — L2 (NEOMOD3) vs L1 (S3M) vs digest2 on the 40-70 deg band, knobs FROZEN by QA3.

Frozen (pre-registered, not tuned on L2): gamma=2.0 (Farnocchia 2015 §3.2 rho^2 spatial
factor, confirmed by AUC), grid 128x64, rho in [0.01,100] au, uniform-in-rho_dot inside the
Admissible Region, class-ratio score, dra_cosdec=False (v5 parquet convention).

Reports both NEOMOD3 normalisations:
  per_H_match -> item 1A: isolates the debiased ORBITAL SHAPE (inherits S3M's H~25 cutoff)
  absolute    -> item 2B territory: NEOMOD3's own N(H) (see faint-H caveat in load_tables)

QA4c: do the tracklets whose score flips L1->L2 live where QA0c predicted (the (a,e)
regions where NEOMOD3 differs from S3M-NEO)?
"""
import os, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src")); sys.path.insert(0, HERE)
import ranging_engine as re
import qa_lib as q
from qa3_knob_sweep import auc, best_f1, spearman

# ---- FROZEN KNOBS (QA3) ----
GAMMA, N_RHO, N_RHODOT = 2.0, 128, 64
N_BAND, BAND = 20000, (40.0, 70.0)
D5_BAR_F1 = 0.839          # pre-registered: digest2 band F1 from the corrected-baseline analysis


def score(df, level, norm="per_H_match", chunk=250):
    return re.score_tracklets(df.mean_ra.values, df.mean_dec.values, df.mean_dra.values,
                              df.mean_ddec.values, df.mean_mag.values, df.mjd0_utc.values,
                              level=level, n_rho=N_RHO, n_rhodot=N_RHODOT, gamma=GAMMA,
                              chunk=chunk, dra_cosdec=False, neomod3_norm=norm)


if __name__ == "__main__":
    df = pd.read_parquet(q.V5, columns=["mean_ra","mean_dec","mean_dra","mean_ddec","mean_mag",
                                        "mjd0_utc","population","ecl_lon","P_NEO_d2","q_au","a_au","e"])
    df["absdlon"] = np.abs(q.antisun_dlon(df.ecl_lon.values, df.mjd0_utc.values))
    band = df[(df.absdlon >= BAND[0]) & (df.absdlon < BAND[1])].sample(N_BAND, random_state=1).reset_index(drop=True)
    y = (band.population == "NEO").values.astype(int)
    d2 = band.P_NEO_d2.values
    print(f"QA4 on {len(band)} tracklets, {BAND[0]:.0f}-{BAND[1]:.0f} deg, NEO frac {y.mean():.3f}")
    print(f"frozen: gamma={GAMMA} grid={N_RHO}x{N_RHODOT}\n")

    res = {}
    for name, lvl, nrm in [("L1 (S3M)", "L1", None),
                           ("L2 (NEOMOD3, per_H_match) = 1A", "L2", "per_H_match"),
                           ("L2 (NEOMOD3, absolute) = 1A+2B", "L2", "absolute")]:
        t = time.time()
        P = score(band, lvl, nrm) if nrm else score(band, lvl)
        res[name] = P
        print(f"{name:34s} F1={best_f1(y,P):.4f}  AUC={auc(y,P):.4f}  ({time.time()-t:.0f}s)")
    print(f"{'digest2 (reference)':34s} F1={best_f1(y,d2):.4f}  AUC={auc(y,d2):.4f}")

    L1, L2p, L2a = res["L1 (S3M)"], res["L2 (NEOMOD3, per_H_match) = 1A"], res["L2 (NEOMOD3, absolute) = 1A+2B"]
    print(f"\n--- D5 pre-registered bar (band F1 >= {D5_BAR_F1}) ---")
    for name, P in res.items():
        print(f"  {name:34s} F1={best_f1(y,P):.4f}  {'GO' if best_f1(y,P)>=D5_BAR_F1 else 'below bar'}")

    print(f"\n--- the actual NEOMOD3 effect (isolated, knobs frozen) ---")
    print(f"  1A   (L2_perH - L1): dF1={best_f1(y,L2p)-best_f1(y,L1):+.4f}  dAUC={auc(y,L2p)-auc(y,L1):+.4f}")
    print(f"  1A+2B(L2_abs  - L1): dF1={best_f1(y,L2a)-best_f1(y,L1):+.4f}  dAUC={auc(y,L2a)-auc(y,L1):+.4f}")
    print(f"  engine alone (L1 - digest2): dF1={best_f1(y,L1)-best_f1(y,d2):+.4f}  dAUC={auc(y,L1)-auc(y,d2):+.4f}")

    # ---- QA4c: do flips land where QA0c predicted? ----
    d = L2p - L1; ok = np.isfinite(d)
    flip = ok & (np.abs(d) > 0.05)
    print(f"\nQA4c: |L2-L1|>0.05 for {flip.sum()} tracklets ({flip.mean()*100:.1f}%)")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    ax[0].scatter(np.log10(band.a_au[ok]), band.e[ok], s=2, c="lightgrey", label="all")
    sc = ax[0].scatter(np.log10(band.a_au[flip]), band.e[flip], s=6, c=d[flip], cmap="RdBu_r",
                       vmin=-0.3, vmax=0.3, label="flipped")
    aa = np.logspace(-0.1, 0.8, 100); ax[0].plot(np.log10(aa), 1-1.3/aa, "k--", lw=1)
    ax[0].set_xlabel("log10 a (TRUE)"); ax[0].set_ylabel("e (TRUE)"); ax[0].set_xlim(-0.05, 0.75); ax[0].set_ylim(0, 1)
    ax[0].set_title("QA4c — where L2 differs from L1 (compare to QA0c ratio map)")
    fig.colorbar(sc, ax=ax[0], label="L2-L1")
    ax[1].hist(d[ok], bins=60, color="C0"); ax[1].set_yscale("log")
    ax[1].set_xlabel("L2 - L1 (per_H_match)"); ax[1].set_ylabel("# tracklets")
    ax[1].set_title(f"NEOMOD3 shift: mean {np.nanmean(d):+.4f}, NEOs {np.nanmean(d[y==1]):+.4f}, "
                    f"non-NEOs {np.nanmean(d[y==0]):+.4f}")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "Figures/qa/qa4_L2_vs_L1.png"), dpi=90, bbox_inches="tight")

    out = band[["population","q_au","a_au","e","absdlon","P_NEO_d2"]].copy()
    out["P_L1"], out["P_L2_perH"], out["P_L2_abs"] = L1, L2p, L2a
    out.to_parquet(os.path.join(ROOT, "outputs/qa4_band_scores.parquet"), index=False)
    print("\nsaved outputs/qa4_band_scores.parquet + Figures/qa/qa4_L2_vs_L1.png")
