# Unvalidated prototypes — DO NOT USE without a parity test

Written during the production scoring run as a faster replacement layout. **Never executed
against production data. No parity test has been run.**

- `sorcha_test_digest2_shard.py` — digest2-only shard. Preserves `repeatable` + `--cpu 1`
  exactly. Adds batch progress logging and immediate per-shard output.
- `sorcha_test_vdp_maplocal.py` — **custom VDP implementation**, shards by map cell.
  *Must not be used until it is shown to reproduce the frozen scorer
  (`sorcha_test_score_shard.py`) bit-for-bit on a common set of rows.*

The production result was produced entirely by the FROZEN scorer
`sorcha_test_score_shard.py` (jobs 38203287 and 38263818).

## Measured facts worth carrying forward
- digest2 `--cpu 1`, uncontended: **3.20 rows/s** (312.2s per 1,000 tracklets)
- digest2 `--cpu 1`, 80-way concurrent: **~1.25 rows/s** — a **2.6x contention penalty**
- VDP cost is dominated by `.npz` map loading, not row count: 14,891 rows over ~667 cells
  took 1009s, while 930 rows over far fewer cells took 30-47s
- `ckpt-all` is a scavenger partition: shard 79 was preempted 3+ times, losing all progress
  each time because shards write only at the end. Use `cpu-g2`, and split work so a
  preemption costs minutes rather than hours.
