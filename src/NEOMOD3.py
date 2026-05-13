import re
import pandas as pd
from urllib.request import urlopen
from io import StringIO
import numpy as np


def getNEOMOD3orbits():
    
    #url = "https://www2.boulder.swri.edu/~davidn/NEOMOD_Simulator/input_neomod3.dat"
    #r = requests.get(url, timeout=30, verify=certifi.where())
    #r.raise_for_status()
    #text = r.text

    #text = urlopen(url).read().decode("utf-8", errors="ignore")
    #lines = text.splitlines()
    with open("input_neomod3.dat") as fp:
        lines = fp.readlines()
    # to find the header row
    header_idx = next(
        i for i, ln in enumerate(lines)
        if re.match(r"^\s*iH\s+ia\s+ie\s+ii\s+binned\s+mod\.\s+nu6\b", ln)
    )

    # merging binned and mod. into one, otherwise it becomes weird :\
    raw_tokens = re.split(r"\s+", lines[header_idx].strip())
    names = []
    i = 0
    while i < len(raw_tokens):
        if raw_tokens[i] == "binned" and i + 1 < len(raw_tokens) and raw_tokens[i+1].startswith("mod."):
            names.append("binned mod.")
            i += 2
        else:
            names.append(raw_tokens[i])
            i += 1

    #skips blank or dashed lines
    data_lines = []
    for ln in lines[header_idx + 1:]:
        s = ln.strip()
        if not s or set(s) <= set("- "):    # separators/blank
            continue
        if re.match(r"^\d", s):             
            data_lines.append(s)


    df = pd.read_csv(StringIO("\n".join(data_lines)), sep=r"\s+", header=None, names=names)
    # df.head()

    array4D, H_center, a_center, e_center, i_center = getNEOMOD3bins(df)
    return df, array4D, H_center, a_center, e_center, i_center

    
def getNEOMOD3bins(df):
    
    #to calc center of bins
    n_H, H_min, H_max = 52, 15.0, 28.0
    n_a, a_min, a_max = 42, 0.0, 4.2
    n_e, e_min, e_max = 25, 0.0, 1.0
    n_i, i_min, i_max = 22, 0.0, 88.0

    da = (a_max - a_min)/(n_a) #bin width
    de = (e_max - e_min)/(n_e)
    di = (i_max - i_min)/(n_i)
    dH = (H_max - H_min)/n_H

    #bin center coords
    a_center = a_min + (np.arange(n_a) + 0.5) * da
    e_center = e_min + (np.arange(n_e) + 0.5) * de
    i_center = i_min + (np.arange(n_i) + 0.5) * di
    H_center = H_min + (np.arange(n_H) + 0.5) * dH

    array4D = np.zeros((n_H, n_a, n_e, n_i))

    idx = df[['iH', 'ia', 'ie', 'ii']].astype(int).values - 1  # 0-based indices
    vals = df['binned mod.'].astype(float).values
    array4D[idx[:,0], idx[:,1], idx[:,2], idx[:,3]] = vals

    df["a_AU"] = a_min + (df["ia"] - 0.5) * da 
    df["e"] = e_min + (df["ie"] - 0.5) * de
    return array4D, H_center, a_center, e_center, i_center

