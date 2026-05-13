# NEOCP Data Collector

This folder archives live MPC NEO Confirmation Page data for later analysis.

What gets collected:

- `raw/`: raw `neocp.txt` and NEOCP HTML snapshots.
- `snapshots/`: one parsed object-list CSV per run.
- `ephemerides/`: raw ephemeris HTML batches and one parsed ephemeris CSV per run.
- `tables/neocp_objects_history.csv`: append-only object-list history.
- `tables/neocp_ephemerides_history.csv`: append-only parsed ephemeris rows.
- `tables/neocp_runs.csv`: one row per collector run with counts/status.
- `logs/neocp_collect.log`: cron/runtime log.

The collector uses MPC's stable text list:

```text
https://minorplanetcenter.net/iau/NEO/neocp.txt
```

For ephemerides, it posts the active checkbox designations from the NEOCP HTML form
to:

```text
https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi
```

Default ephemerides are geocentric, 1-hour interval, starting now, full output,
full sexagesimal RA/Dec, total motion in arcsec/minute.

Run manually:

```bash
cd /Users/devanshisingh/Downloads/research/NEO_probability/neomod/neocp_data
./run_neocp_collect.sh
```

Install cron manually:

```bash
crontab /Users/devanshisingh/Downloads/research/NEO_probability/neomod/neocp_data/neocp_crontab.txt
```

Cron schedule is every 4 hours:

```text
0 */4 * * *
```

Important: this archives live MPC-published list data and MPC-generated
ephemerides. It does not VDP-score these objects. VDP scoring would require
matching probability maps for the appropriate epoch/field and conversion of
ephemeris positions/motions into the VDP feature frame.

MPC caution: NEOCP observations/orbits/ephemerides are preliminary. Keep this as
a local research archive unless the objects/data have appeared in published MPC
products such as MPECs/MPS.
