# s3m_loader.py

import os
import glob
import pandas as pd
import numpy as np


# Map population -> filename glob pattern
# MBA:      "S1_*.s3m"
# NEO only: "S0.s3m"
# TNO:      "ST.s3m"
# Trojans:  "St5.s3m"
# Full:     "S*.s3m"
POP_TO_GLOB = {
    "mba": "S1_*.s3m",
    "neo": "S0.s3m",
    "tno": "ST.s3m",
    "trojan": "St5.s3m",
    "all": "S*.s3m",
}

DEFAULT_COLS = [
    "OID", "FORMAT", "q", "e", "i", "node", "argperi",
    "t_p", "H", "t_0", "INDEX", "N_PAR", "MOID", "COMPCODE"
]


def define_s3m(pop="mba", pattern=None, cols=None, verbose=True, search_dirs=None):
    """
    Load one or more S3M files into a single dataframe.

    Parameters
    ----------
    pop : str, optional
        Population name from:
        "mba", "neo", "tno", "trojan", "all".
        Ignored if pattern is provided.
    pattern : str or None, optional
        Explicit glob pattern e.g. "S1_*.s3m" or "/path/to/files/S1_*.s3m".
        If None, uses POP_TO_GLOB[pop].
    cols : list of str or None, optional
        Column names for the S3M files.
    verbose : bool, optional
        If True, print search information and load summary.
    search_dirs : list of str or None, optional
        Directories to search when pattern is relative.
        Default is [".", "S3Mdata"].

    Returns
    -------
    s3m : pandas DataFrame
        Concatenated S3M dataframe with an added semi-major axis column "a".
    """
    if cols is None:
        cols = DEFAULT_COLS

    if pattern is None:
        try:
            pattern = POP_TO_GLOB[pop.lower()]
        except KeyError as e:
            raise ValueError(
                f"Unknown pop='{pop}'. Choose from {list(POP_TO_GLOB.keys())} "
                f"or pass pattern=..."
            ) from e

    if search_dirs is None:
        search_dirs = [".", "S3Mdata"]

    files = []

    
    if os.path.isabs(pattern) or os.path.dirname(pattern):
        files = sorted(glob.glob(pattern))
    else:
        for d in search_dirs:
            test_pattern = os.path.join(d, pattern)
            matched = sorted(glob.glob(test_pattern))
            if matched:
                files.extend(matched)

    
    files = list(dict.fromkeys(files))

    if verbose:
        print("Search pattern:", pattern)
        print("Search dirs:", search_dirs)
        print("Matched files:", files)

    if len(files) == 0:
        raise FileNotFoundError(
            f"No files matched pattern '{pattern}' in directories {search_dirs}"
        )

    dfs = []
    for f in files:
        df = pd.read_csv(
            f,
            sep=r"\s+",
            comment="!",
            names=cols,
            engine="python"
        )
        dfs.append(df)

    s3m = pd.concat(dfs, ignore_index=True)

    # Add semi-major axis column
    s3m["a"] = s3m["q"] / (1 - s3m["e"])

    if verbose:
        print(f"Loaded {len(files)} file(s).")
        print(f"Final dataframe shape: {s3m.shape}")
        print("The final df is done.")

    return s3m


def s3m_array(s3m, n_H=52, n_a=42, n_e=25, n_i=22):
    """
    Build a 4D histogram in (H, a, e, i) from S3M dataframe.

    Parameters
    ----------
    s3m : pandas DataFrame
        S3M dataframe containing columns H, a, e, i.
    n_H, n_a, n_e, n_i : int, optional
        Number of bins in each dimension.

    Returns
    -------
    s3m : pandas DataFrame
        Original dataframe.
    array4D : ndarray
        4D histogram counts in (H, a, e, i).
    H_center, a_center, e_center, i_center : ndarray
        Bin centers in each dimension.
    """
    a_min, a_max = s3m["a"].min(), s3m["a"].max()
    e_min, e_max = s3m["e"].min(), s3m["e"].max()
    i_min, i_max = s3m["i"].min(), s3m["i"].max()
    H_min, H_max = s3m["H"].min(), s3m["H"].max()

    points = np.vstack([s3m["H"], s3m["a"], s3m["e"], s3m["i"]]).T

    edges = [
        np.linspace(H_min, H_max, n_H + 1),
        np.linspace(a_min, a_max, n_a + 1),
        np.linspace(e_min, e_max, n_e + 1),
        np.linspace(i_min, i_max, n_i + 1),
    ]

    array4D, edges_out = np.histogramdd(points, bins=edges)

    H_edges, a_edges, e_edges, i_edges = edges_out
    H_center = 0.5 * (H_edges[1:] + H_edges[:-1])
    a_center = 0.5 * (a_edges[1:] + a_edges[:-1])
    e_center = 0.5 * (e_edges[1:] + e_edges[:-1])
    i_center = 0.5 * (i_edges[1:] + i_edges[:-1])

    return s3m, array4D, H_center, a_center, e_center, i_center


def build_scorer(nsc_module, pop="mba", pattern=None, **define_kwargs):
    """
    Convenience function to build an NEOMOD scorer from S3M files.

    Parameters
    ----------
    nsc_module : module
        Module containing NEOMODScorer.
    pop : str, optional
        Population key passed to define_s3m().
    pattern : str or None, optional
        Explicit glob pattern for S3M files.
    **define_kwargs
        Extra keyword arguments passed to define_s3m().

    Returns
    -------
    scorer, s3m, array4D, Hc, ac, ec, ic
    """
    s3m = define_s3m(pop=pop, pattern=pattern, **define_kwargs)
    s3m, array4D, Hc, ac, ec, ic = s3m_array(s3m)
    scorer = nsc_module.NEOMODScorer(array4D, Hc, ac, ec, ic)
    return scorer, s3m, array4D, Hc, ac, ec, ic