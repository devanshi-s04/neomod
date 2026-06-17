# hybrid_loader.py
# Drop-in replacement for s3m_loader.py.
# Used by velocity_density_pipeline_hybrid.py via:
#   import hybrid_loader as load_s3m
# The pipeline calls: load_s3m.define_s3m(pop=...) and load_s3m.s3m_array(df)

import pandas as pd

HYBRID_PARQUET = "/mmfs1/gscratch/dirac/ds2004/sorcha/hybrid_elements.parquet"

# Maps the same pop keys s3m_loader accepts -> population column label in parquet
POP_MAP = {
    "neo":    "NEO",
    "mba":    "MBA",
    "tno":    "TNO",
    "trojan": "Trojan",
    "all":    None,   # no filter, return all populations
}


def define_s3m(pop="neo", verbose=True, **kwargs):
    """Load hybrid catalog for one population.

    Drop-in for s3m_loader.define_s3m. Accepts the same pop keys:
    "neo", "mba", "tno", "trojan", "all".

    Returns a DataFrame with columns:
        OID, q, e, i, node, argperi, t_p, H, t_0, a, population
    (superset of what s3m_array needs: H, a, e, i)
    """
    pop_key = pop.lower()
    if pop_key not in POP_MAP:
        raise ValueError(
            f"Unknown pop='{pop}'. Choose from {list(POP_MAP.keys())}"
        )

    df = pd.read_parquet(HYBRID_PARQUET)

    label = POP_MAP[pop_key]
    if label is not None:
        df = df[df["population"] == label].reset_index(drop=True)

    if verbose:
        print(f"[hybrid_loader] pop={pop!r}: {len(df):,} objects loaded from {HYBRID_PARQUET}")

    return df


def s3m_array(df, n_H=52, n_a=42, n_e=25, n_i=22):
    """Build the 4D (H, a, e, i) histogram used by NEOMODScorer.

    Identical interface to s3m_loader.s3m_array — delegates to it since
    the DataFrame column format is the same.
    """
    import s3m_loader
    return s3m_loader.s3m_array(df, n_H, n_a, n_e, n_i)
