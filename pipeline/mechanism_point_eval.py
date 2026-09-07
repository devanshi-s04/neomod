#!/usr/bin/env python3
"""Point-evaluation harness for existing-data mechanism ablations.

Instead of building 1001x1001 maps for every ablation variant, this evaluates the SAME sealed
kNN density estimator directly at the TEST2 tracklet velocities. That is ~9,330 query points per
variant instead of 1,002,001 grid points, so a mechanism can be varied and rescored in seconds
without writing a single map file.

Two facts make this exact rather than approximate:

  1. `build_visible_subset_dataframe` applies a purely GEOMETRIC cut (separation from the centre),
     so it is independent of magnitude. The visible subset is therefore built ONCE per
     (centre, population) and filtered by magnitude afterwards. Verified against the frozen build's
     own per-bin `n_visible` records: exact match.

  2. The density estimator is the sealed closed form
         n0 = ( k(k+1)/2 - 1/2 ) / ( pi * sum_j d_j^2 )
     evaluated at whatever points it is given. The frozen pipeline evaluates it on a grid and then
     bilinearly interpolates to the tracklet; this harness evaluates it at the tracklet directly.
     The only difference is the grid+bilinear step, which `validate_against_frozen` measures.

NOTHING here writes a map, modifies the frozen map root, or touches the frozen v1 scorer.
"""
from __future__ import annotations
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
for p in ("neomod/pipeline", "neomod/src", "neomod/adam_core_stub"):
    sys.path.insert(0, str(W/p))

import build_neomod3_mag025_k150_maps as FZ

POPS    = ("NEO", "MBA", "TNO", "Trojans")
NONNEO  = ("MBA", "TNO", "Trojans")
FROZEN_K = {"NEO": 150, "MBA": 10, "TNO": 10, "Trojans": 10}


def closed_form(d, k):
    """Sealed closed-form posterior mean of n0 = 1/(pi d0^2). d: (n_query, k) distances."""
    S = np.einsum("ij,ij->i", d, d)
    return (k*(k+1)/2.0 - 0.5) / (np.pi * S)


def load_nonneo_gen(verbose=True):
    """Read the epoch-state cache once and return the GEN parents per non-NEO population."""
    cache = pd.read_parquet(FZ.EPOCH_CACHE)
    man   = pd.read_parquet(FZ.SPLIT_MANIFEST, columns=["ObjID","split"])
    gen   = set(man.ObjID[man.split == "GEN"])
    out = {}
    for pop in NONNEO:
        out[pop] = cache[(cache.population==pop) & (cache.ObjID.isin(gen))].reset_index(drop=True)
        if verbose: print(f"  {pop}: {len(out[pop]):,} GEN parents", flush=True)
    del cache
    return out


class CentreSamples:
    """Visible GEN samples for one centre, built once, reusable for every mechanism variant."""

    def __init__(self, centre, epoch=None, verbose=True):
        import velocity_density_pipeline_neomod_clone_only as base
        self.centre = centre
        seal = json.load(open(W/"outputs/splits/MAP_BUILD_SEAL.json"))
        self.epoch = epoch or seal["grid"]["ref_obstime"]
        dlon, lat = FZ.parse_center(centre)
        self.clon, self.clat = FZ.center_lonlat(self.epoch, dlon, lat)
        self.f_by = json.load(open(FZ.SPLIT_MAG025))["f_by_population_magbin"]
        self.eff_neo = None
        self.vis = {}

    def load(self, neo_source_df, eff_neo, nonneo_gen=None, verbose=True):
        """neo_source_df: the same all-sky NEO GEN slice the frozen builder used for this centre.
        nonneo_gen: optional {pop: DataFrame} from load_nonneo_gen(), so the 2.6 GB epoch cache is
        read once for a whole study rather than once per centre."""
        import velocity_density_pipeline_neomod_clone_only as base
        _, scorer = base.load_s3m_population("neo", verbose=False)
        self.eff_neo = float(eff_neo)
        src = nonneo_gen if nonneo_gen is not None else load_nonneo_gen()
        for pop in NONNEO:
            sub = src[pop]
            self.vis[pop] = base.build_visible_subset_dataframe(
                sub, obstime_str=self.epoch, scorer=scorer, max_sep_deg=FZ.MAX_SEP_DEG,
                chunk=100_000, show_progress=False, center_mode="custom_ecliptic",
                center_lon_deg=self.clon, center_lat_deg=self.clat)
            if verbose: print(f"  [{self.centre}] {pop}: {len(self.vis[pop]):,} visible", flush=True)
        self.vis["NEO"] = base.build_visible_subset_dataframe(
            neo_source_df, obstime_str=self.epoch, scorer=scorer, max_sep_deg=FZ.MAX_SEP_DEG,
            chunk=100_000, show_progress=False, center_mode="custom_ecliptic",
            center_lon_deg=self.clon, center_lat_deg=self.clat)
        if verbose: print(f"  [{self.centre}] NEO: {len(self.vis['NEO']):,} visible", flush=True)
        return self

    def eff(self, pop, bin_label):
        """Physical normalisation, identical rule to the frozen builder."""
        if pop == "NEO":
            return self.eff_neo
        f = self.f_by.get(pop, {}).get(bin_label)
        return None if f is None else float(f)

    def _sample_weights(self, pop, mags):
        """Per-sample physical weight 1/f_split(pop, the sample's own 0.25-mag bin).

        The frozen builder never needs this because its windows ARE the 0.25-mag bins, so one
        scalar suffices. A wider ablation window spans several bins whose split fractions differ,
        so the weight must be per sample. NEO's effective factor carries no bin dependence, so for
        NEO this is a constant and the wider-window ablation is exact.
        Returns (weights, n_unweightable) -- samples in a bin with no split fraction get NaN and
        are excluded rather than zero-filled."""
        if pop == "NEO":
            return np.full(len(mags), 1.0/self.eff_neo), 0
        edges = 14.0 + 0.25*np.floor((mags - 14.0)/0.25)
        w = np.full(len(mags), np.nan)
        for lo in np.unique(edges[np.isfinite(edges)]):
            f = self.f_by.get(pop, {}).get(f"V{lo:06.2f}_{lo+0.25:06.2f}")
            if f is not None:
                w[edges == lo] = 1.0/float(f)
        return w, int(np.isnan(w).sum())

    def density_at(self, pop, lo, hi, k, xs, ys, bin_label_for_eff=None, workers=8):
        """Sealed kNN density at the query points, over samples with lo <= mag_app < hi.

        rho(x) = n0(x) * mean_{j in kNN(x)} w_j, with n0 the sealed closed form and w_j the
        per-sample physical weight. For a window equal to one 0.25-mag bin every w_j is identical
        and this is exactly the frozen n0/eff. Returns (density, n_samples, reason); NaN when
        invalid, never zero-filled."""
        from scipy.spatial import cKDTree
        s = self.vis[pop]
        m = s["mag_app"].to_numpy(float)
        sel = np.isfinite(m) & (m >= lo) & (m < hi)
        if int(sel.sum()) == 0:
            return np.full(len(xs), np.nan), 0, "no_samples"
        sub = s[sel]
        mags = sub["mag_app"].to_numpy(float)
        w, n_bad = self._sample_weights(pop, mags)
        keep = np.isfinite(w)
        if not keep.any():
            return np.full(len(xs), np.nan), 0, "no_split_fraction"
        sub = sub[keep]; w = w[keep]
        n = len(sub)
        if n < 2:  return np.full(len(xs), np.nan), n, "below_min"
        if n <= k: return np.full(len(xs), np.nan), n, "insufficient_support"
        tree = cKDTree(np.column_stack([sub["vlam"].to_numpy(float),
                                        sub["vbeta"].to_numpy(float)]))
        d, idx = tree.query(np.column_stack([xs, ys]), k=int(k), workers=int(workers))
        if d.ndim == 1: d = d[:, None]; idx = idx[:, None]
        return closed_form(d, k) * w[idx].mean(axis=1), n, "ok"

    def neighbour_radius(self, pop, lo, hi, k, xs, ys, workers=8):
        """Distance to the k-th neighbour at each query point -- the effective bandwidth."""
        from scipy.spatial import cKDTree
        s = self.vis[pop]
        m = s["mag_app"].to_numpy(float)
        sel = np.isfinite(m) & (m >= lo) & (m < hi)
        if int(sel.sum()) <= k: return None
        sub = s[sel]
        tree = cKDTree(np.column_stack([sub["vlam"].to_numpy(float),
                                        sub["vbeta"].to_numpy(float)]))
        d, _ = tree.query(np.column_stack([xs, ys]), k=int(k), workers=int(workers))
        if d.ndim == 1: d = d[:, None]
        return d[:, -1]


def load_neo_allsky_slice(clon, clat, epoch,
                          allsky_dir=W/"outputs"/"neomod3_projection_cache_high_allsky"):
    """The same conservative HEALPix superset slice the frozen builder takes for a centre,
    copied from build_neomod3_mag025_k150_maps.py so the NEO sample is identical."""
    import healpy as hp
    import pyarrow.dataset as pads, pyarrow.compute as pc
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, GCRS, GeocentricTrueEcliptic
    import astropy.units as u
    bypix = Path(allsky_dir)/"by_pixel"
    t_obs = Time(epoch, scale="tdb")
    c = SkyCoord(lon=clon*u.deg, lat=clat*u.deg, distance=1.0*u.AU,
                 frame=GeocentricTrueEcliptic(obstime=t_obs)).transform_to(GCRS(obstime=t_obs))
    vec = hp.ang2vec(float(c.ra.deg), float(c.dec.deg), lonlat=True)
    rad = np.radians(FZ.MAX_SEP_DEG + np.degrees(hp.max_pixrad(8)))
    want = hp.query_disc(8, vec, rad, inclusive=True).tolist()
    dset = pads.dataset(str(bypix), format="parquet", partitioning="hive")
    neo_df = dset.to_table(filter=pc.field("pix").isin(want), use_threads=True).to_pandas()
    man = json.load(open(Path(allsky_dir)/"manifest.json"))
    return neo_df, float(man["effective_factor_NEO"])


def four_population_score(cs, rows, k_by_pop=None, bin_width=0.25, workers=8):
    """P(NEO) = rho_NEO / sum_c rho_c at the tracklet points. All four densities required, exactly
    as the frozen scorer requires them. Returns a DataFrame aligned to `rows`."""
    k_by_pop = k_by_pop or FROZEN_K
    n = len(rows)
    out = {f"rho_{p}": np.full(n, np.nan) for p in POPS}
    P = np.full(n, np.nan); reason = np.array(["ok"]*n, dtype=object)
    xs = rows["vlam"].to_numpy(float); ys = rows["vbeta"].to_numpy(float)
    mv = rows["mean_mag_V"].to_numpy(float)
    lo_edges = 14.0 + bin_width*np.floor((mv - 14.0)/bin_width)
    inb = (np.abs(xs) <= 5.0) & (np.abs(ys) <= 5.0)
    for lo in np.unique(lo_edges[np.isfinite(lo_edges)]):
        ii = np.where(lo_edges == lo)[0]
        hi = lo + bin_width
        if lo < 14.0 or hi > 25.0 + 1e-9:
            reason[ii] = "v_out_of_range"; continue
        dens, miss = {}, []
        for p in POPS:
            dd, ns, why = cs.density_at(p, lo, hi, k_by_pop[p], xs[ii], ys[ii], workers=workers)
            if why != "ok": miss.append(p)
            else: dens[p] = dd; out[f"rho_{p}"][ii] = dd
        if miss:
            reason[ii] = "missing_population_density:" + ",".join(miss); continue
        tt = sum(dens[p] for p in POPS)
        good = inb[ii] & np.isfinite(tt) & (tt > 0)
        pn = np.full(len(ii), np.nan); pn[good] = dens["NEO"][good]/tt[good]
        P[ii] = pn
        reason[ii] = np.where(~inb[ii], "outside_velocity_grid",
                      np.where(~np.isfinite(tt), "nonfinite_density",
                      np.where(tt <= 0, "zero_total_density", "ok")))
    df = pd.DataFrame({"P_NEO": P, "reason": reason}, index=rows.index)
    for p in POPS: df[f"rho_{p}"] = out[f"rho_{p}"]
    df["valid"] = np.isfinite(P)
    return df
