#!/usr/bin/env python3
"""qa_lib.py — plotting functions for the 1A ranging-engine QA harness (NEOMplan §8).
Each qaXX() returns a matplotlib Figure. Runnable as a script to regenerate all
figures into Figures/qa/. The notebook 1A_engine_qa.ipynb imports these so the heavy
logic is testable outside Jupyter."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
import ranging_engine as re
import neomod3_sampler as nm

V5 = os.path.join(ROOT, "outputs/phase2/sorcha_comparison_v5_masked.parquet")
CACHE = os.path.join(ROOT, "outputs/pop_cache_wide.npz")
FIGDIR = os.path.join(ROOT, "Figures/qa")
AU_DAY_TO_KMS = re.AU_KM / 86400.0     # 1 AU/day in km/s


# ------------------------------------------------------------------ shared context
def load_ctx():
    arr, cen, wid = nm.load_neomod3_array()
    z = np.load(CACHE)
    edges = {k: z[k] for k in ("H_edges", "loga_edges", "e_edges", "i_edges")}
    hists = {k: z[k] for k in ("hist_neo", "hist_mba", "hist_trojan", "hist_tno")}
    nm_wide = re.neomod3_to_wide(edges)
    return dict(nm_arr=arr, nm_cen=cen, edges=edges, hists=hists, nm_wide=nm_wide)


def _cc(edges_key, edges):
    e = edges[edges_key]; return 0.5 * (e[1:] + e[:-1])


def antisun_dlon(ecl_lon, mjd):
    T = (mjd + 2400000.5) - 2451545.0
    anti = ((280.46 + 0.9856474 * T) % 360 + 180) % 360
    return ((ecl_lon - anti + 180) % 360) - 180


# ================================================================== QA0 inputs
def qa0a_neomod3(ctx):
    """NEOMOD3 marginals + (a,e) density at two H slices (shape must change with H)."""
    arr, cen = ctx["nm_arr"], ctx["nm_cen"]
    fig, ax = plt.subplots(2, 3, figsize=(14, 7))
    for k, name in zip(range(4), ["H", "a", "e", "i"]):
        axis = tuple(j for j in range(4) if j != k)
        m = arr.sum(axis=axis)
        r, c = divmod(k, 3)
        ax[r, c].plot(cen[name], m, lw=2)
        ax[r, c].set_xlabel(name); ax[r, c].set_ylabel("count")
        ax[r, c].set_title(f"NEOMOD3 marginal {name}")
        if name == "H":
            ax[r, c].set_yscale("log")
    # (a,e) at H=18 and H=25
    for j, Hval in enumerate([18.0, 25.0]):
        iH = int(np.argmin(np.abs(cen["H"] - Hval)))
        slab = arr[iH].sum(axis=2)          # (a,e)
        axx = ax[1, 1] if j == 0 else ax[1, 2]
        im = axx.pcolormesh(cen["a"], cen["e"], slab.T, shading="auto", cmap="viridis")
        axx.set_xlabel("a [au]"); axx.set_ylabel("e"); axx.set_title(f"NEOMOD3 (a,e) @ H={Hval:.0f}")
        # q=1.3 boundary a(1-e)=1.3
        aa = np.linspace(cen["a"][0], cen["a"][-1], 100); axx.plot(aa, 1 - 1.3/aa, "r--", lw=1)
        axx.set_ylim(0, 1)
        fig.colorbar(im, ax=axx, fraction=0.046)
    fig.suptitle("QA0a — NEOMOD3 numerator: marginals + (a,e) shape change with H "
                 "(faint bins should differ from bright)", fontsize=11)
    fig.tight_layout(); return fig


def qa0b_s3m(ctx):
    """S3M denominator per-class marginals (log10a, e, i, H) + bin occupancy."""
    edges, hists = ctx["edges"], ctx["hists"]
    loga = _cc("loga_edges", edges); ec = _cc("e_edges", edges)
    ic = _cc("i_edges", edges); Hc = _cc("H_edges", edges)
    # histogram axes are (H=0, loga=1, e=2, i=3)
    axes_specs = [(1, loga, "log10 a"), (2, ec, "e"), (3, ic, "i [deg]"), (0, Hc, "H")]
    fig, ax = plt.subplots(2, 3, figsize=(14, 7))
    colors = {"hist_mba": "C0", "hist_trojan": "C1", "hist_tno": "C2", "hist_neo": "k"}
    for k, (haxis, cc, lbl) in enumerate(axes_specs):
        r, c = divmod(k, 3)
        for hk, col in colors.items():
            axsum = tuple(j for j in range(4) if j != haxis)
            ax[r, c].plot(cc, hists[hk].sum(axis=axsum), col, lw=1.5, label=hk.replace("hist_", ""))
        ax[r, c].set_xlabel(lbl); ax[r, c].set_yscale("log"); ax[r, c].set_title(f"S3M marginal {lbl}")
        if k == 0:
            ax[r, c].legend(fontsize=8)
    # occupancy: fraction of non-empty bins per class
    nonneo = hists["hist_mba"] + hists["hist_trojan"] + hists["hist_tno"]
    ax[1, 1].hist(np.log10(nonneo[nonneo > 0].ravel()), bins=50, color="C3")
    ax[1, 1].set_xlabel("log10 count in occupied non-NEO bins"); ax[1, 1].set_ylabel("# bins")
    ax[1, 1].set_title(f"denominator occupancy: {np.mean(nonneo>0)*100:.1f}% bins filled")
    # (log10a, e) of denominator
    im = ax[1, 2].pcolormesh(loga, ec, np.log10(nonneo.sum(axis=(0, 3)).T + 1),
                             shading="auto", cmap="magma")
    ax[1, 2].set_xlabel("log10 a"); ax[1, 2].set_ylabel("e"); ax[1, 2].set_title("non-NEO log count (a,e)")
    fig.colorbar(im, ax=ax[1, 2], fraction=0.046)
    fig.suptitle("QA0b — S3M denominator: MBA a~1.8-3.3, Trojan a~5.2, TNO a~30-70", fontsize=11)
    fig.tight_layout(); return fig


def qa0c_ratio(ctx):
    """NEOMOD3 / S3M-NEO conditional-shape ratio in (a,e) and per-H. Shows WHERE L2!=L1."""
    edges = ctx["edges"]; nm_wide = ctx["nm_wide"]; s3m_neo = ctx["hists"]["hist_neo"]
    loga = _cc("loga_edges", edges); ec = _cc("e_edges", edges); Hc = _cc("H_edges", edges)
    # per-H normalise each to unit sum over (loga,e,i), then ratio in (loga,e)
    def norm_ae(h):
        s = h.sum(axis=(1, 2, 3), keepdims=True)
        hh = np.divide(h, np.maximum(s, 1e-30))
        return hh.sum(axis=(0, 3))  # collapse H,i -> (loga,e) of per-H-normalised
    nmn = norm_ae(nm_wide); s3n = norm_ae(s3m_neo)
    ratio = np.log2(np.divide(nmn, s3n, out=np.full_like(nmn, np.nan), where=s3n > 0))
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    im = ax[0].pcolormesh(loga, ec, ratio.T, shading="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    aa = 10**loga; ax[0].plot(loga, 1 - 1.3/aa, "k--", lw=1)
    ax[0].set_xlabel("log10 a"); ax[0].set_ylabel("e")
    ax[0].set_title("log2(NEOMOD3 / S3M-NEO) in (a,e)  [red=NEOMOD3 richer]")
    ax[0].set_xlim(loga[0], np.log10(4.2)); ax[0].set_ylim(0, 1)
    fig.colorbar(im, ax=ax[0], fraction=0.046)
    # per-H total (N(H) shape, normalised)
    nH = nm_wide.sum(axis=(1, 2, 3)); sH = s3m_neo.sum(axis=(1, 2, 3))
    ax[1].plot(Hc, nH/nH.sum(), label="NEOMOD3 N(H)"); ax[1].plot(Hc, sH/sH.sum(), label="S3M-NEO N(H)")
    ax[1].set_xlabel("H"); ax[1].set_ylabel("normalised N(H)"); ax[1].legend()
    ax[1].set_title("N(H) shape (L2 rescales NEOMOD3 to S3M per-H; this is the 2B lever)")
    fig.suptitle("QA0c — where NEOMOD3 disagrees with S3M-NEO (predicts L2-L1 flip regions)", fontsize=11)
    fig.tight_layout(); return fig


# ================================================================== QA1 geometry
def _mesh_elements(row, n_rho=90, n_rhodot=90, rho=(0.01, 100), rhodot_auday=0.05):
    """Elements on a RECTANGULAR (rho,rhodot) mesh for one tracklet (Farnocchia Fig-1 style)."""
    r_obs, v_obs = re.earth_observer_state([row.mjd0_utc])
    l_hat, v_ang = re.tracklet_vectors(np.array([row.mean_ra]), np.array([row.mean_dec]),
                                       np.array([row.mean_dra]), np.array([row.mean_ddec]), dra_cosdec=False)
    r_obs, v_obs, l_hat, v_ang = r_obs[0], v_obs[0], l_hat[0], v_ang[0]
    rho_au = np.logspace(np.log10(rho[0]), np.log10(rho[1]), n_rho)
    rhodot = np.linspace(-rhodot_auday, rhodot_auday, n_rhodot) * AU_DAY_TO_KMS   # km/s
    RHO, RD = np.meshgrid(rho_au, rhodot, indexing="ij")                          # (nr,nrd)
    rho_km = RHO * re.AU_KM
    r_helio = r_obs + rho_km[..., None] * l_hat
    v = v_obs + RD[..., None] * l_hat + rho_km[..., None] * v_ang
    r_ecl = r_helio @ re._R_EQ2ECL.T; v_ecl = v @ re._R_EQ2ECL.T
    r = np.linalg.norm(r_ecl, axis=-1); v2 = np.einsum('...k,...k->...', v_ecl, v_ecl)
    energy = 0.5*v2 - re.MU_SUN/r; a = -re.MU_SUN/(2*energy)/re.AU_KM
    h = np.cross(r_ecl, v_ecl); hn = np.linalg.norm(h, axis=-1)
    evec = np.cross(v_ecl, h)/re.MU_SUN - r_ecl/r[..., None]; e = np.linalg.norm(evec, axis=-1)
    inc = np.degrees(np.arccos(np.clip(h[..., 2]/hn, -1, 1)))
    q = a*(1-e)
    # H(rho)
    r_sun_obj = np.linalg.norm(r_helio, axis=-1)/re.AU_KM
    obj_sun = -r_helio; obj_obs = r_obs - r_helio
    ca = np.einsum('...k,...k->...', obj_sun, obj_obs)/(np.linalg.norm(obj_sun,axis=-1)*np.linalg.norm(obj_obs,axis=-1))
    alpha = np.arccos(np.clip(ca, -1, 1)); ta = np.tan(np.clip(alpha,0,np.pi-1e-6)/2)
    phase = -2.5*np.log10(0.85*np.exp(-3.33*ta**0.63)+0.15*np.exp(-1.87*ta**1.22))
    H = row.mean_mag - 5*np.log10(r_sun_obj*RHO) - phase
    bound = (energy < 0) & (e < 1)
    return dict(rho=rho_au, rhodot_auday=rhodot/AU_DAY_TO_KMS, a=a, e=e, q=q, inc=inc, H=H, bound=bound)


def qa1b_farnocchia(exemplars):
    """Element contours over (rho,rhodot) for 4 exemplar tracklets — Farnocchia Fig 1."""
    fig, axes = plt.subplots(len(exemplars), 4, figsize=(16, 3.4*len(exemplars)))
    for ri, (label, row) in enumerate(exemplars):
        m = _mesh_elements(row)
        X, Y = np.meshgrid(m["rho"], m["rhodot_auday"], indexing="ij")
        for ci, (fld, name, levels) in enumerate([
                ("e", "eccentricity", [0.2,0.5,0.8,1.0]),
                ("q", "perihelion q [au]", [0.1,0.5,1.0,1.3,3.0]),
                ("inc", "inclination [deg]", [2,5,10,30]),
                ("H", "H", np.arange(15,34,3))]):
            ax = axes[ri, ci] if len(exemplars) > 1 else axes[ci]
            Z = np.where(m["bound"], m[fld], np.nan)
            cs = ax.contour(X, Y, Z, levels=levels, colors="k", linewidths=0.8)
            ax.clabel(cs, fontsize=6, fmt="%g")
            # AR boundary (bound edge) shaded
            ax.contourf(X, Y, m["bound"].astype(float), levels=[0.5,1.5], colors=["#cde"], alpha=0.4)
            ax.set_xscale("log"); ax.set_xlabel("rho [au]"); ax.set_ylabel("rhodot [au/d]")
            if ri == 0: ax.set_title(name, fontsize=9)
            if ci == 0: ax.text(-0.35, 0.5, label, transform=ax.transAxes, rotation=90,
                                va="center", fontsize=9, weight="bold")
    fig.suptitle("QA1b — element maps over (rho,rhodot) per exemplar (cf. Farnocchia 2015 Fig 1; "
                 "shaded=bound orbits). H contours ~vertical is the V-d-H coupling.", fontsize=11)
    fig.tight_layout(); return fig


def qa1c_ar(df_sample):
    """Admissible fraction vs rho per elongation band."""
    bands = [(0,20),(40,70),(70,110),(110,180)]
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    for lo, hi in bands:
        sub = df_sample[(df_sample.absdlon>=lo)&(df_sample.absdlon<hi)]
        if len(sub) < 20: continue
        sub = sub.sample(min(200, len(sub)), random_state=0)
        r_obs, v_obs = re.earth_observer_state(sub.mjd0_utc.values)
        l_hat, v_ang = re.tracklet_vectors(sub.mean_ra.values, sub.mean_dec.values,
                                           sub.mean_dra.values, sub.mean_ddec.values, dra_cosdec=False)
        g = re.build_grid(l_hat, v_ang, r_obs, v_obs, n_rho=64, n_rhodot=32)
        frac = g["admissible"].any(axis=2).mean(axis=0)   # frac of tracklets admissible at each rho
        rho = g["rho_au"][0, :, 0]
        ax.plot(rho, frac, label=f"{lo}-{hi} deg", lw=2)
    ax.set_xscale("log"); ax.set_xlabel("rho [au]"); ax.set_ylabel("fraction of tracklets with bound orbit")
    ax.legend(); ax.set_title("QA1c — Admissible Region closes at large rho (bounded-orbit condition)")
    fig.tight_layout(); return fig


# ================================================================== QA2 prior / L0
def qa2a_weights(row):
    """Node-weight maps for gamma=0,2,4 on one tracklet's admissible grid."""
    r_obs, v_obs = re.earth_observer_state([row.mjd0_utc])
    l_hat, v_ang = re.tracklet_vectors(np.array([row.mean_ra]), np.array([row.mean_dec]),
                                       np.array([row.mean_dra]), np.array([row.mean_ddec]), dra_cosdec=False)
    g = re.build_grid(l_hat, v_ang, r_obs, v_obs, n_rho=64, n_rhodot=32)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for j, gam in enumerate([0, 2, 4]):
        w = re.node_weights(g, gamma=gam)[0]        # (nr,nrd)
        rho = g["rho_au"][0, :, 0]; rd = g["rho_dot"][0]     # rd (nr,nrd) km/s
        wm = np.where(w > 0, np.log10(w + 1e-30), np.nan)
        im = ax[j].pcolormesh(np.log10(rho), np.arange(w.shape[1]), wm.T, shading="auto", cmap="viridis")
        ax[j].set_xlabel("log10 rho [au]"); ax[j].set_ylabel("rhodot node index")
        ax[j].set_title(f"log10 weight, gamma={gam}"); fig.colorbar(im, ax=ax[j], fraction=0.046)
    fig.suptitle("QA2a — geometric prior weight maps (gamma shifts mass outward in rho)", fontsize=11)
    fig.tight_layout(); return fig


def qa2b_L0(df_sample):
    """L0 (geometric) score distributions by class."""
    P = re.score_tracklets(df_sample.mean_ra.values, df_sample.mean_dec.values,
                           df_sample.mean_dra.values, df_sample.mean_ddec.values,
                           df_sample.mean_mag.values, df_sample.mjd0_utc.values,
                           level="L0", chunk=1000, dra_cosdec=False)
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    for pop, col in [("NEO","C3"),("MBA","C0"),("TNO","C2"),("Trojan","C1")]:
        m = (df_sample.population == pop).values & np.isfinite(P)
        if m.sum() > 5:
            ax.hist(P[m], bins=30, histtype="step", density=True, lw=2, color=col,
                    label=f"{pop} (med {np.median(P[m]):.2f})")
    ax.set_xlabel("L0 geometric P(NEO)"); ax.set_ylabel("density"); ax.legend()
    ax.set_title("QA2b — L0 geometric score by class (already separates; = Spoto score)")
    fig.tight_layout(); return fig, P


def qa2c_pvalue(n=400, gamma=2.0):
    """Truth p-value histogram: where the TRUE (rho,rhodot) falls in the posterior weight
    ranking. Uniform => healthy prior; spike at 0 => pathological (Farnocchia's Jeffreys test)."""
    D = "/astro/users/jkurla/public_html/LSST_Sorcha_predictions"
    df = pd.read_hdf(f"{D}/one_day_neo.h5").reset_index(drop=True)
    df = df.sample(min(n, len(df)), random_state=0).reset_index(drop=True)
    r_obs, v_obs = re.earth_observer_state(df.fieldMJD_TAI.values)
    l_hat, v_ang = re.tracklet_vectors(df.RATrue_deg.values, df.DecTrue_deg.values,
                                       df.RARateCosDec_deg_day.values, df.DecRate_deg_day.values)
    g = re.build_grid(l_hat, v_ang, r_obs, v_obs, n_rho=64, n_rhodot=32)
    w = re.node_weights(g, gamma=gamma)                         # (N,nr,nrd)
    # posterior weight at the node nearest the truth (rho,rhodot)
    rho_true = df.Range_LTC_km.values / re.AU_KM
    rd_true = df.RangeRate_LTC_km_s.values
    rho_nodes = g["rho_au"]; rd_nodes = g["rho_dot"]
    d2 = (np.log10(rho_nodes) - np.log10(rho_true)[:, None, None])**2 + \
         ((rd_nodes - rd_true[:, None, None]) / 30.0)**2
    flat = d2.reshape(len(df), -1); nn = np.argmin(flat, axis=1)
    wflat = w.reshape(len(df), -1)
    w_at_truth = wflat[np.arange(len(df)), nn]
    # p-value = fraction of admissible weight LESS likely than truth's node
    pval = np.array([(wflat[k][wflat[k] > 0] < w_at_truth[k]).mean() if w_at_truth[k] > 0 else 0.0
                     for k in range(len(df))])
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.hist(pval, bins=20, color="C4", edgecolor="k")
    ax.set_xlabel("p-value of TRUE (rho,rhodot) under prior weight"); ax.set_ylabel("# objects")
    ax.set_title(f"QA2c — truth p-value (gamma={gamma}); uniform=healthy, "
                 f"spike@0=pathological (Farnocchia Table 3)")
    fig.tight_layout(); return fig, pval


# ------------------------------------------------------------------ exemplar picker
def pick_exemplars(df):
    df = df.copy(); df["absdlon"] = np.abs(antisun_dlon(df.ecl_lon.values, df.mjd0_utc.values))
    def one(mask):
        s = df[mask]; return s.iloc[len(s)//2] if len(s) else None
    ex = []
    ex.append(("NEO", one((df.population=="NEO") & (df.absdlon<30))))
    ex.append(("MBA antisun", one((df.population=="MBA") & (df.absdlon<20))))
    ex.append(("MBA 40-70", one((df.population=="MBA") & (df.absdlon>=40) & (df.absdlon<70))))
    ex.append(("TNO", one(df.population=="TNO")))
    return [(l, r) for l, r in ex if r is not None]


if __name__ == "__main__":
    os.makedirs(FIGDIR, exist_ok=True)
    print("loading context + sample ...")
    ctx = load_ctx()
    df = pd.read_parquet(V5, columns=["mean_ra","mean_dec","mean_dra","mean_ddec","mean_mag",
                                      "mjd0_utc","population","ecl_lon","q_au","P_NEO_d2"])
    df["absdlon"] = np.abs(antisun_dlon(df.ecl_lon.values, df.mjd0_utc.values))
    samp = df.groupby("population", group_keys=False).apply(lambda d: d.sample(min(400,len(d)), random_state=0))
    ex = pick_exemplars(df)
    figs = {
        "qa0a_neomod3": qa0a_neomod3(ctx),
        "qa0b_s3m": qa0b_s3m(ctx),
        "qa0c_ratio": qa0c_ratio(ctx),
        "qa1b_farnocchia": qa1b_farnocchia(ex),
        "qa1c_ar": qa1c_ar(samp),
        "qa2a_weights": qa2a_weights(ex[0][1]),
    }
    figs["qa2b_L0"] = qa2b_L0(samp)[0]
    figs["qa2c_pvalue"] = qa2c_pvalue()[0]
    for name, fig in figs.items():
        p = os.path.join(FIGDIR, name + ".png"); fig.savefig(p, dpi=90, bbox_inches="tight")
        print("saved", p)
    print("done")
