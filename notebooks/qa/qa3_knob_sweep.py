#!/usr/bin/env python3
"""QA3 — knob-freeze sweep (NEOMplan §8).

Freeze the engine's GEOMETRY knobs (prior exponent gamma, grid resolution) by maximising
L1<->digest2 agreement on the 40-70 deg band, BEFORE ever looking at L2. This is the
solo-work discipline: the L2-L1 comparison then has zero tuned parameters.

Efficiency: observer state + grid + elements + H depend only on resolution; gamma only
re-weights. So loop resolution (outer, expensive) x gamma (inner, cheap), chunked over
tracklets to bound memory.

Also: QA2c p-value uniformity vs gamma (KS distance from uniform) as a prior-health tiebreak.
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

N_BAND = 4000           # tracklets sampled from the 40-70 deg band
GAMMAS = [0.0, 2.0, 4.0]
RESOLUTIONS = [(32, 16), (64, 32), (128, 64)]
BAND = (40.0, 70.0)


# ---------- metrics ----------
def _rank(x):
    o = np.argsort(x, kind="mergesort"); r = np.empty(len(x), float); r[o] = np.arange(len(x))
    return r

def spearman(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    ra, rb = _rank(a[m]), _rank(b[m])
    return np.corrcoef(ra, rb)[0, 1]

def auc(y, s):
    m = np.isfinite(s); y = y[m]; s = s[m]
    r = _rank(s); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1*(n1-1)/2) / (n1*n0)

def best_f1(y, s):
    m = np.isfinite(s); y = y[m].astype(float); s = s[m]
    o = np.argsort(-s, kind="mergesort"); ys = y[o]
    tp = np.cumsum(ys); fp = np.cumsum(1-ys)
    P = tp/np.maximum(tp+fp, 1e-12); R = tp/max(ys.sum(), 1)
    return np.max(2*P*R/np.maximum(P+R, 1e-12))


# ---------- sweep ----------
def run_sweep(df_band, chunk=250):
    y = (df_band.population == "NEO").values.astype(int)
    d2 = df_band.P_NEO_d2.values
    neo_tab, nonneo_tab, edges = re.load_tables("L1")
    rows = []
    for (nr, nrd) in RESOLUTIONS:
        t0 = time.time()
        P = {g: np.full(len(df_band), np.nan) for g in GAMMAS}
        for s in range(0, len(df_band), chunk):
            e = min(s+chunk, len(df_band))
            sub = df_band.iloc[s:e]
            r_obs, v_obs = re.earth_observer_state(sub.mjd0_utc.values)
            l_hat, v_ang = re.tracklet_vectors(sub.mean_ra.values, sub.mean_dec.values,
                                               sub.mean_dra.values, sub.mean_ddec.values, dra_cosdec=False)
            g_ = re.build_grid(l_hat, v_ang, r_obs, v_obs, n_rho=nr, n_rhodot=nrd)
            el = re.elements_from_nodes(g_, l_hat, v_ang)
            Hn = re.H_from_nodes(g_, sub.mean_mag.values)
            for gam in GAMMAS:                      # cheap: only re-weight
                w = re.node_weights(g_, gamma=gam)
                P[gam][s:e] = re.class_score(el, Hn, w, neo_tab, nonneo_tab, edges, "L1")
            del g_, el, Hn
        for gam in GAMMAS:
            rows.append(dict(n_rho=nr, n_rhodot=nrd, gamma=gam,
                             spearman_d2=spearman(P[gam], d2),
                             auc_L1=auc(y, P[gam]), f1_L1=best_f1(y, P[gam]),
                             nan_frac=np.mean(~np.isfinite(P[gam]))))
        print(f"  res {nr}x{nrd} done in {time.time()-t0:.0f}s", flush=True)
    res = pd.DataFrame(rows)
    res["auc_d2"] = auc(y, d2); res["f1_d2"] = best_f1(y, d2)
    return res


def pvalue_uniformity(gammas=GAMMAS, n=400):
    """KS distance of the truth p-value distribution from uniform, per gamma (prior health)."""
    out = {}
    for gam in gammas:
        _, pv = q.qa2c_pvalue(n=n, gamma=gam)
        pv = np.sort(pv); cdf = np.arange(1, len(pv)+1)/len(pv)
        out[gam] = float(np.max(np.abs(cdf - pv)))     # KS vs uniform
        plt.close("all")
    return out


if __name__ == "__main__":
    df = pd.read_parquet(q.V5, columns=["mean_ra","mean_dec","mean_dra","mean_ddec","mean_mag",
                                        "mjd0_utc","population","ecl_lon","P_NEO_d2","q_au"])
    df["absdlon"] = np.abs(q.antisun_dlon(df.ecl_lon.values, df.mjd0_utc.values))
    band = df[(df.absdlon >= BAND[0]) & (df.absdlon < BAND[1])]
    band = band.sample(N_BAND, random_state=0).reset_index(drop=True)
    print(f"QA3 sweep on {len(band)} tracklets in {BAND[0]:.0f}-{BAND[1]:.0f} deg "
          f"(NEO frac {np.mean(band.population=='NEO'):.3f})")

    res = run_sweep(band)
    print("\n=== L1 <-> digest2 agreement (freeze on max spearman) ===")
    print(res.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    best = res.loc[res.spearman_d2.idxmax()]
    print(f"\nBEST: n_rho={int(best.n_rho)} n_rhodot={int(best.n_rhodot)} gamma={best.gamma} "
          f"spearman={best.spearman_d2:.4f} auc_L1={best.auc_L1:.4f} (auc_d2={res.auc_d2.iloc[0]:.4f})")

    print("\n=== prior health: KS(truth p-value, uniform) per gamma (lower=healthier) ===")
    ks = pvalue_uniformity()
    for g, v in ks.items(): print(f"  gamma={g}: KS={v:.3f}")

    # convergence check at best gamma: 64x32 vs 128x64
    sub = res[res.gamma == best.gamma].sort_values("n_rho")
    if len(sub) >= 2:
        d = abs(sub.f1_L1.iloc[-1] - sub.f1_L1.iloc[-2])
        print(f"\nconvergence |dF1| between {int(sub.n_rho.iloc[-2])}x{int(sub.n_rhodot.iloc[-2])} "
              f"and {int(sub.n_rho.iloc[-1])}x{int(sub.n_rhodot.iloc[-1])} = {d:.4f} "
              f"({'PASS <0.002' if d < 0.002 else 'above 0.002 threshold'})")

    os.makedirs(os.path.join(ROOT, "outputs"), exist_ok=True)
    res.to_csv(os.path.join(ROOT, "outputs/qa3_knob_sweep.csv"), index=False)
    # heatmap
    piv = res.pivot(index="gamma", columns="n_rho", values="spearman_d2")
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(piv.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(piv.columns)), [f"{c}x{c//2}" for c in piv.columns])
    ax.set_yticks(range(len(piv.index)), [f"γ={g:g}" for g in piv.index])
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{piv.values[i,j]:.3f}", ha="center", va="center", color="w", fontsize=9)
    ax.set_title("QA3 — Spearman(L1, digest2) on 40-70° band\n(knobs frozen at the max)")
    fig.colorbar(im); fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "Figures/qa/qa3_knob_sweep.png"), dpi=90, bbox_inches="tight")
    print("\nsaved outputs/qa3_knob_sweep.csv + Figures/qa/qa3_knob_sweep.png")
