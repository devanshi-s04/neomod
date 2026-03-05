# s3m_loader.py
import glob
import pandas as pd
import numpy as np

# Map population -> filename glob pattern 
# MBA:      "S1_*.s3m"  :contentReference[oaicite:0]{index=0}
# NEO only: "S0.s3m"    :contentReference[oaicite:1]{index=1}
# TNO:      "ST.s3m"    :contentReference[oaicite:2]{index=2}
# Trojans:  "St5.s3m"   :contentReference[oaicite:3]{index=3}
# Full:     "S*.s3m"    :contentReference[oaicite:4]{index=4}
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


def define_s3m(pop="mba", pattern=None, cols=None, verbose=True):
#    Load S3M files into a single dataframe.

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

    files = sorted(glob.glob(pattern))
    if verbose:
        print(files)

    if len(files) == 0:
        raise FileNotFoundError(
            f"No files matched pattern '{pattern}'. "
            f"Are you in the directory with the .s3m files?"
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
    # add semi major axis column
    s3m["a"] = s3m["q"] / (1 - s3m["e"])

    if verbose:
        print("The final df is done.")
    return s3m


def s3m_array(s3m, n_H=52, n_a=42, n_e=25, n_i=22):

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
    Returns (scorer, s3m, array4D, Hc, ac, ec, ic)

    Use:
        scorer, s3m, array4D, Hc, ac, ec, ic = sm.build_scorer(nsc, pop="neo")
    """
    s3m = define_s3m(pop=pop, pattern=pattern, **define_kwargs)
    s3m, array4D, Hc, ac, ec, ic = s3m_array(s3m)
    scorer = nsc_module.NEOMODScorer(array4D, Hc, ac, ec, ic)
    return scorer, s3m, array4D, Hc, ac, ec, ic