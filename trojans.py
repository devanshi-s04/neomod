import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import importlib
from functools import reduce

# to download the s3m files, access the epyc filesystem /epyc/projects/jpl_survey_sim/S3M_v09.05.15/


def define_s3m():
    cols = [
        "OID", "FORMAT", "q", "e", "i", "node", "argperi",
        "t_p", "H", "t_0", "INDEX", "N_PAR", "MOID", "COMPCODE"
    ]
    
    
    files = sorted(glob.glob("St5.s3m"))
    print(files)
    
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
        
    
    
    trojans_s3m = pd.concat(dfs, ignore_index=True)
    
    trojans_s3m["a"] = trojans_s3m["q"] / (1 - trojans_s3m["e"])
    print("The final df is done.")
    return trojans_s3m

def s3m_array(trojans_s3m):
# defining the min/max values for a,e,i,H
    a_min = trojans_s3m["a"].min()
    a_max = trojans_s3m["a"].max()
    e_min = trojans_s3m["e"].min()
    e_max = trojans_s3m["e"].max()
    i_min = trojans_s3m["i"].min()
    i_max = trojans_s3m["i"].max()
    H_min = trojans_s3m["H"].min()
    H_max = trojans_s3m["H"].max()

    n_H, n_a, n_e, n_i = 52,42,25,22
#n_H, n_a, n_e, n_i = 152,142,125,122
#n_H, n_a, n_e, n_i = 352,342,325,322
    points = np.vstack([trojans_s3m["H"], trojans_s3m["a"], trojans_s3m["e"], trojans_s3m["i"]]).T
    edges = [
        np.linspace(H_min, H_max, n_H + 1),
        np.linspace(a_min, a_max, n_a + 1),
        np.linspace(e_min, e_max, n_e + 1),
        np.linspace(i_min, i_max, n_i + 1)
    ]
    
    
    array4D, edges_out = np.histogramdd(points, bins=edges)
    
    H_edges = edges_out[0]
    a_edges = edges_out[1]
    e_edges = edges_out[2]
    i_edges = edges_out[3]
    
    
    H_center = 0.5 * (H_edges[1:] + H_edges[:-1])
    a_center = 0.5 * (a_edges[1:] + a_edges[:-1])
    e_center = 0.5 * (e_edges[1:] + e_edges[:-1])
    i_center = 0.5 * (i_edges[1:] + i_edges[:-1])

    return trojans_s3m, array4D, H_center, a_center, e_center, i_center


