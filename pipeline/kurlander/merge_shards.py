#!/usr/bin/env python3
"""Merge Slurm shard parquets into one tracklet table and build the eval referee.

Two referee variants are written:
  *_natural.parquet — every tracklet, true Rubin population mix (NEO frac ~1e-3).
                      Use AUC here (prevalence-independent).
  *_eval.parquet    — non-NEO subsampled so NEO fraction matches the v5 eval set
                      (default 0.293), so F1 is directly comparable to the v5 numbers.
                      *** F1 is prevalence-dependent: never compare F1 across referees
                      with different NEO fractions. ***
"""
import argparse, glob
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--glob", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--neo-frac", type=float, default=0.293, help="target NEO fraction for the eval set")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

files = sorted(glob.glob(a.glob))
print(f"merging {len(files)} shards")
df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
print(f"total tracklets: {len(df):,}")
print(df.population.value_counts().to_dict())
print(f"natural NEO frac: {df.is_neo.mean():.6f}")

nat = a.out.replace(".parquet", "_natural.parquet")
df.to_parquet(nat, index=False); print("wrote", nat)

rng = np.random.default_rng(a.seed)
neo = df[df.is_neo]; non = df[~df.is_neo]
n_non = int(len(neo) * (1 - a.neo_frac) / a.neo_frac)
if n_non < len(non):
    non = non.iloc[rng.choice(len(non), n_non, replace=False)]
ev = pd.concat([neo, non], ignore_index=True).sample(frac=1, random_state=a.seed)
ev.to_parquet(a.out, index=False)
print(f"wrote {a.out}: {len(ev):,} tracklets, NEO frac {ev.is_neo.mean():.4f}")
