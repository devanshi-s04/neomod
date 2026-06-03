# VDP / digest2 context, 2026-05-06

## Restored previous context: Digest2 vs VDP Context Handoff

Date: 2026-05-05

Project directory:

```text
/Users/devanshisingh/Downloads/research/NEO_probability
```

Main files discussed:

```text
neomod/digest2_comparison.ipynb
neomod/run_digest2_comparison.py
neomod/velocity_density_pipeline.py
neomod/ROCcurve.ipynb
neomod/paper/NEOrocks.tex
neomod/paper/NEOrocks.pdf
```

### Big Picture

We are comparing two ways of assigning Near-Earth Object probability to synthetic moving-object tracklets:

1. **VDP**, the velocity-density probability method in this project.
2. **digest2**, the MPC-style classifier used for NEOCP-style candidate scoring.

The goal is to compare them on the same synthetic objects, under the same observing geometry, and understand why the plots and metrics differ.

The corrected result is:

```text
VDP recovers substantially more true NEOs at roughly the same contamination level.
VDP best F1 ~= 0.845
digest2 best F1 ~= 0.663
```

This means VDP is not just 3-4% better. In the corrected `digest2_comparison` run, VDP is much more complete at similar contamination.

### Corrected Geometry And Assumptions

Earlier notebook text incorrectly said the comparison used:

```text
|beta| < 3 deg
observatory 568 / Mauna Kea
```

That was wrong for this comparison.

The corrected setup is:

```text
Date: 2025-03-21
Sky region: circular patch centered on geocentric ecliptic lon=180 deg, lat=0 deg
Sky radius: 30 deg
Magnitude range: 14 <= mag_app <= 26
Observatory: Rubin Observatory in Chile
MPC observatory code for digest2 tracklets: X05
Tracklet: two detections separated by 30 minutes
```

Important distinction:

```text
The old |beta| < 3 deg filter was for ecliptic-strip plots only.
The digest2 comparison should use the same 30-degree opposition-centered patch as VDP.
```

Rubin observatory geometry comes from `neomod/neoscore.py`, which hardcodes Rubin's approximate geodetic location:

```text
longitude = -70.7366 deg
latitude  = -30.2407 deg
height    = 2647 m
```

### Code Changes Made

File changed:

```text
neomod/run_digest2_comparison.py
```

Main fixes:

1. Set digest2 observatory code to Rubin/X05:

```python
OBSCODE = "X05"
```

2. Removed the incorrect `|beta| < 3 deg` filtering logic.

3. Removed custom propagation/scoring code that duplicated VDP behavior.

4. Reused the canonical VDP scoring path from `velocity_density_pipeline.py`:

```python
_, pop_df = prob_map_set.score_orbital_df(
    df=df,
    scorer=scorer,
    obstime_str=OBSTIME_STR,
    max_sep_deg=prob_map_set.max_sep_deg,
    chunk=50_000,
    show_progress=False,
    return_visible=True,
)
```

5. Used `data["P_NEO"]` returned by `score_orbital_df` as `data["P_NEO_vdp"]`.

6. Removed the ad hoc `v_sky >= 0.2` cut, because in VDP the `0.2` parameter refers to nearest-clone support in velocity-map space, not a sky-speed selection cut.

7. Added a checkpoint after VDP scoring:

```text
neomod/s3m_digest2_comparison_vdp_input.parquet
```

8. Fixed digest2 config permission issues by using temporary config files instead of writing into the external digest2 checkout.

9. Fixed digest2 timeout issues by chunking digest2 calls:

```python
DIGEST2_CHUNK_TRACKLETS = 5_000
DIGEST2_TIMEOUT_SEC = 1_800
```

10. Added parser warnings for duplicate, unparsed, or missing digest2 IDs.

### Commands Run

Dependency installed because Astropy's `de432s` ephemeris needed it:

```bash
python -m pip install jplephem
```

The corrected comparison was run from:

```bash
cd /Users/devanshisingh/Downloads/research/NEO_probability/neomod
MPLCONFIGDIR=/private/tmp/mplconfig python run_digest2_comparison.py
```

The run took about an hour, mostly due to digest2 chunking.

Digest2 was run in 8 chunks:

```text
7 chunks of 5000 tracklets
1 final chunk of 3045 tracklets
```

One tiny issue:

```text
38,045 digest2 output lines
38,044 unique parsed IDs
1 object got the default digest2 score 0.0
```

This is negligible for the metrics, and the script now warns about cases like this.

The notebook was later executed successfully with:

```bash
jupyter nbconvert --to notebook --execute digest2_comparison.ipynb --inplace --ExecutePreprocessor.timeout=300
```

This needed escalation because the Jupyter kernel binds local ports.

### Fresh Generated Files

Fresh corrected output files:

```text
neomod/s3m_digest2_comparison.parquet
neomod/s3m_digest2_comparison_vdp_input.parquet
neomod/roc_comparison_vdp_digest2.png
```

Do not rely on:

```text
neomod/digest2_comparison_log.txt
```

That log may still contain stale output from the earlier incorrect beta-strip/Mauna-Kea-style run.

### Corrected Dataset Size

The corrected comparison has:

```text
Total rows: 38,045
N_NEO: 2,233
N_non-NEO: 35,812
```

Population counts:

```text
MBA        15,350
Trojans    13,035
TNO         7,427
NEO         2,233
```

### Corrected Metrics

Notebook output:

```text
N_NEO=2233  N_non-NEO=35812
VDP     best: t=0.028  completeness=76.7%  contamination=6.0%  F1=0.845
Digest2 best: t=0.970  completeness=51.3%  contamination=6.4%  F1=0.663
```

Script output may show VDP threshold as about `0.030` depending on the threshold grid, with essentially the same interpretation.

Core interpretation:

```text
VDP recovers about three quarters of the true NEOs while keeping contamination near 6%.
digest2 must use a very high threshold near 0.970 to keep contamination similarly low, but then it recovers only about half of the true NEOs.
```

### F1 Meaning

Definitions used in the notebook:

```text
completeness = fraction of true NEOs recovered
contamination = fraction of selected candidates that are false positives
purity = 1 - contamination
```

F1 combines purity and completeness:

```text
F1 = 2 * purity * completeness / (purity + completeness)
```

A high F1 means the selected candidate list is both complete and clean.

For this project:

```text
VDP F1 ~= 0.845
digest2 F1 ~= 0.663
```

VDP has the better balance.

## 2026-05-07 Updates: Notebook Diagnostics And NEOCP Archive

### `digest2_comparison.ipynb` Plot Updates

The notebook was updated and re-executed successfully with:

```bash
cd /Users/devanshisingh/Downloads/research/NEO_probability/neomod
jupyter nbconvert --to notebook --execute digest2_comparison.ipynb --inplace --ExecutePreprocessor.timeout=300
```

As before, notebook execution needs permission to bind local Jupyter kernel ports.

The score-distribution figure in section 3 now has three panels:

```text
left:   VDP NEO probability, P_NEO, by true population
middle: direct object-by-object scatter of VDP P_NEO vs digest2 NEO score
right:  digest2 NEO score by true population
```

The middle direct-score plot showed the important qualitative behavior:

```text
Most non-NEOs form a vertical wall near VDP P_NEO = 0.
digest2 spreads some of those same low-VDP objects up to high digest2 scores.
True NEOs cover a wider part of the plane, including many high-score regions.
The two classifiers agree on many objects but are not redundant.
```

Interpretation:

```text
VDP is much more selective about assigning nonzero NEO probability to MBAs,
TNOs, and Trojans.

digest2 is sensitive to many true NEOs, but also gives high scores to some
objects VDP thinks are strongly non-NEO. This helps explain why digest2 needs
a very high threshold near 0.970 to keep contamination low.
```

The ROC/threshold figure labels were cleaned up:

```text
Legends now use short labels only.
Best-F1 threshold values are annotated on the plot instead of making the legend huge.
The threshold-axis label was simplified to "Classifier threshold".
The threshold lines are visible and labeled by nearby text annotations.
```

### New Diagnostic Sections Added To `digest2_comparison.ipynb`

Several new sections were added after the ROC section.

#### 5. Score-threshold quadrants

This plot reuses the direct VDP-vs-digest2 score scatter and draws:

```text
vertical line:   VDP best-F1 threshold
horizontal line: digest2 best-F1 threshold
```

It divides the score plane into four regions:

```text
both select      = high VDP, high digest2
VDP only         = high VDP, low digest2
digest2 only     = low VDP, high digest2
neither selects  = low VDP, low digest2
```

The notebook also displays a population summary table for these quadrants.

Interpretation:

```text
The "VDP only" quadrant shows candidates VDP recovers but digest2 misses at its
best-F1 threshold.

The "digest2 only" quadrant is a contamination diagnostic: these are objects
digest2 selects even though VDP assigns low NEO probability.
```

#### 6. 2D score density

Added hexbin density plots for:

```text
all populations
true NEOs
non-NEOs
```

Purpose:

```text
The scatter plot is heavily overplotted near P_NEO = 0. The hexbin views show
the dense low-VDP wall and separate true-NEO vs non-NEO score structure.
```

#### 7. Rank comparison

Added rank-rank comparison:

```text
x-axis: log10 VDP descending score rank
y-axis: log10 digest2 descending score rank
```

The plot includes guide lines for top-N candidate-list cutoffs.

Purpose:

```text
Thresholds answer "what happens at a fixed score cutoff?"
Rank comparison answers "are the top follow-up candidates the same objects?"
```

This is useful because operational follow-up usually works from a ranked queue,
not only a calibrated score threshold.

#### 8. Top-N follow-up yield

Added observer-facing curves:

```text
x-axis: number of highest-scoring candidates followed up
y-axis, panel 1: NEO completeness
y-axis, panel 2: contamination among followed objects
```

The notebook also displays a summary table at selected top-N values.

Purpose:

```text
This answers: if an observer can follow up only the top 50, 100, 500, 1000,
... candidates, which classifier recovers more true NEOs and how contaminated
is the follow-up list?
```

#### 9. Live NEOCP published-score context

An optional live NEOCP section was added to the notebook.

Important clarification:

```text
This section does NOT VDP-score NEOCP objects.
It does NOT propagate live NEOCP RA/Dec to 2025-03-21.
It does NOT use the VDP probability maps.
```

What it does:

```text
Pulls the current MPC NEOCP tabular page.
Reads the MPC-published Score column.
Normalizes the Score column from 0-100 to 0-1.
Compares that live MPC-published score distribution to the synthetic digest2
score distribution.
```

The notebook now prints this warning explicitly:

```text
these are MPC-published NEOCP scores only; no VDP scoring or propagation to
2025-03-21 is performed here.
```

Why the distinction matters:

```text
The current VDP maps are for 2025-03-21 and a specific opposition-centered
field geometry. Live NEOCP objects are May 2026 objects, in different sky
locations and observing circumstances. A direct VDP-vs-NEOCP comparison
requires either new VDP maps for matching epochs/fields, or a carefully defined
mapping from live ephemerides into the appropriate VDP feature frame.
```

### NEOCP Longitudinal Archive Plan

The current plan is:

```text
We cannot easily get the live NEOCP state from a year ago.
Instead, start archiving the live NEOCP list now.
Later, adapt/recompute VDP maps for the dates/fields of archived NEOCP objects
and compare VDP to MPC/digest2-style NEOCP scores on real candidates.
```

Possible historical data sources may exist, but they are not guaranteed to
reconstruct the same live hourly NEOCP page state:

```text
MPC previous-designation/outcome pages
MPECs
JPL Scout
NEODyS / NEOScan
third-party NEOCP mirrors
```

These may be useful later, but the robust path is to build our own time-series
archive from now forward.

### New `neomod/neocp_data/` Collector

Created a new folder:

```text
neomod/neocp_data/
```

Main files:

```text
neomod/neocp_data/collect_neocp.py
neomod/neocp_data/run_neocp_collect.sh
neomod/neocp_data/neocp_crontab.txt
neomod/neocp_data/README.md
```

The collector is self-contained and standard-library-only so cron can run it
without requiring the notebook environment.

It uses MPC's stable text list:

```text
https://minorplanetcenter.net/iau/NEO/neocp.txt
```

This file matches the NEOCP table layout, but with:

```text
R.A. converted to decimal hours
Decl. converted to decimal degrees
```

For ephemerides, the collector replays the form used by the NEOCP page:

```text
POST https://cgi.minorplanetcenter.net/cgi-bin/confirmeph2.cgi
```

It selects active checkbox objects from:

```text
https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html
```

Default ephemeris settings:

```text
geocentric ephemerides
1-hour interval
start now + 0 hours
full output
full sexagesimal RA/Dec
total motion in arcsec/minute and position angle
no sun/altitude suppression
```

Important implementation note:

```text
MPC's form rejected "0.0" for numeric fields such as start time.
The collector now formats these compactly as "0" using Python's :g format.
```

### NEOCP Data Products

The collector writes both raw and parsed data.

Raw/per-run files:

```text
neomod/neocp_data/raw/
neomod/neocp_data/snapshots/
neomod/neocp_data/ephemerides/
```

Append-only tables:

```text
neomod/neocp_data/tables/neocp_objects_history.csv
neomod/neocp_data/tables/neocp_ephemerides_history.csv
neomod/neocp_data/tables/neocp_runs.csv
```

Runtime logs:

```text
neomod/neocp_data/logs/neocp_collect.log
```

Fields in `neocp_objects_history.csv` include:

```text
snapshot_id
fetched_at_utc
designation
score
score_0_1
discovery_date_mpc
ra_hours
ra_deg
dec_deg
v_mag
status
updated_mpc
note
nobs
arc_days
h_mag
not_seen_days
source_url
```

Fields in `neocp_ephemerides_history.csv` include:

```text
snapshot_id
fetched_at_utc
batch_id
designation
ephem_time_utc
ra_hms
ra_deg
dec_dms
dec_deg
elong_deg
v_mag
motion_arcsec_per_min
pa_deg
uncertainty
source_url
```

Each cron run automatically appends to the history tables if it succeeds. It
also writes per-run snapshots and raw MPC responses. If parsing logic changes
later, the raw HTML/text files can be re-parsed.

### Cron Job Installed

Installed user crontab:

```cron
# NEOCP list + ephemeris collector: every 4 hours.
SHELL=/bin/bash
PATH=/opt/anaconda3/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
PYTHON_BIN=/opt/anaconda3/bin/python

0 */4 * * * /Users/devanshisingh/Downloads/research/NEO_probability/neomod/neocp_data/run_neocp_collect.sh
```

This runs every 4 hours local time.

Manual run command:

```bash
cd /Users/devanshisingh/Downloads/research/NEO_probability/neomod/neocp_data
./run_neocp_collect.sh
```

A limited test run can be done with:

```bash
./run_neocp_collect.sh --max-objects 1 --request-pause 0
```

### NEOCP Collector Test Results

Initial sandboxed runs failed because network access was restricted:

```text
URLError: nodename nor servname provided, or not known
```

After allowing network access, one limited ephemeris test succeeded:

```text
objects=20
active ephemeris objects=17
ephemeris rows=25
```

Then one full seed run succeeded:

```text
objects=20
active ephemeris objects=17
ephemeris rows=425
```

The successful full seed run means the archive already contains:

```text
current NEOCP object list snapshot
current MPC-generated geocentric ephemerides for active NEOCP objects
parsed append-only table rows
raw MPC response files
```

### MPC Data Caution

The NEOCP page itself warns that observations and orbits available through the
service are preliminary. The collector README notes:

```text
Keep this as a local research archive unless the objects/data have appeared in
published MPC products such as MPECs/MPS.
```

The NEOCP archive is therefore intended for local research and later VDP method
development, not public redistribution of preliminary MPC NEOCP material.

### Per-Population Breakdown

Fresh per-population summary:

```text
MBA      n=15350  VDP median=0.000 above_thr=0.5% | d2 median=0.010 above_thr=0.0%
NEO      n=2233   VDP median=0.860 above_thr=76.7% | d2 median=0.980 above_thr=51.3%
TNO      n=7427   VDP median=0.002 above_thr=0.1% | d2 median=0.670 above_thr=1.1%
Trojans  n=13035  VDP median=0.000 above_thr=0.2% | d2 median=0.050 above_thr=0.0%
```

Interpretation:

```text
VDP strongly separates NEOs from non-NEOs.
Most MBAs, TNOs, and Trojans receive near-zero VDP probabilities.
NEOs have a high VDP median and high recovery fraction.
```

For digest2:

```text
NEOs often score very high, with median around 0.980.
However, TNOs also receive many moderately/high digest2 scores, with median around 0.670.
The high threshold of 0.970 keeps TNO false positives mostly under control, but at the cost of losing many true NEOs.
```

Important TNO caveat:

```text
This comparison is conditional on synthetic TNOs being visible in the selected sky patch and being fed to digest2 as Rubin/X05 30-minute tracklets.
In the real NEOCP submission stream, far TNOs are less likely to be detected/submitted at all.
So the TNO behavior here is a classifier stress test, not a statement that real NEOCP would be flooded with TNOs.
```

### Score Distribution Plots

VDP score histogram:

```text
True NEOs cluster toward high VDP P_NEO.
Non-NEOs pile up near P_NEO = 0.
This shows VDP is strongly population-selective under the corrected geometry.
```

digest2 score histogram:

```text
Many NEOs also score near 1.0.
But digest2 scores are messier across populations, especially for TNOs.
The high density of digest2 NEO scores near 1.0 does not automatically mean better classification, because non-NEOs can also receive substantial digest2 scores.
```

Why digest2 can show more NEO density near 1.0 than VDP:

```text
digest2 scores are more saturated/quantized near the top.
VDP is a continuous probability from a velocity-density map and spreads NEO scores more gradually.
The ROC/F1 metrics matter more than the visual height of the near-1.0 histogram bin.
```

### ROC Plot Interpretation

ROC-style plot axes:

```text
x-axis = completeness (%)
y-axis = contamination (%)
```

Ideal location:

```text
bottom-right = high completeness, low contamination
```

Observed behavior:

```text
VDP stays at lower contamination while reaching higher completeness.
digest2 contamination rises more quickly as completeness increases.
```

Best-F1 points:

```text
VDP: about 76.7% completeness at about 6.0% contamination
digest2: about 51.3% completeness at about 6.4% contamination
```

Conclusion:

```text
At similar contamination, VDP recovers many more true NEOs.
```

### Threshold vs Performance Plot

VDP:

```text
The best threshold is low, around P_NEO > 0.028.
Even a low VDP threshold is meaningful because most non-NEOs have VDP scores near zero.
```

digest2:

```text
The best threshold is high, around digest2 score > 0.970.
This is needed because digest2 gives moderately high scores to more non-NEOs.
The high threshold keeps contamination down but reduces NEO completeness.
```

### VDP Observable Plots

The notebook includes three VDP observable panels:

```text
VDP P_NEO vs v_lambda
VDP P_NEO vs v_beta
VDP P_NEO vs mag_app
```

Definitions:

```text
v_lambda = apparent motion in ecliptic longitude, deg/day
v_beta = apparent motion in ecliptic latitude, deg/day
mag_app = apparent magnitude
```

Interpretation:

```text
v_lambda is the clearest separator in this geometry.
True NEOs occupy a broad faster-motion region, especially around roughly -0.8 to -0.3 deg/day.
VDP assigns high NEO probability to that region.
Most MBAs, TNOs, and Trojans sit near low VDP scores.
```

For `v_beta`:

```text
v_beta adds information about out-of-ecliptic motion and inclination-like behavior.
It is useful, but less cleanly separating by itself than v_lambda.
```

For `mag_app`:

```text
Magnitude is not acting alone as a NEO classifier.
In VDP, magnitude selects the velocity-density map/bin being used.
The classification comes from the local density of NEO and non-NEO populations in velocity-magnitude space.
```

Threshold line:

```text
VDP threshold ~= 0.028
Objects above this line are selected as VDP NEO candidates.
```

### Digest2 Observable Plots

The notebook also includes digest2 score plotted against:

```text
v_lambda
v_beta
mag_app
```

Important caveat:

```text
These plots are diagnostic only.
They show how digest2 scores correlate with these observables in our sample.
They do not prove digest2 uses ecliptic v_lambda or v_beta internally.
```

digest2 receives:

```text
MPC-style tracklets: two sky positions, two times, magnitude, and observatory code X05.
```

So digest2 may be using information related to apparent motion, but not necessarily in the same explicit ecliptic-velocity-map way as VDP.

Interpretation:

```text
NEOs often score near digest2 = 1.0.
TNOs can also receive high digest2 scores, especially in the faint magnitude range.
Trojans and MBAs mostly remain lower, but some overlap appears.
The strict digest2 threshold near 0.970 is what prevents these high-scoring non-NEO populations from becoming large contamination.
```

Threshold line:

```text
digest2 threshold ~= 0.970
Objects above this line are selected as digest2 NEO candidates.
```

### Research Poster Caption Drafts

Short captions for the three figure groups:

1. Score distributions:

```text
VDP and digest2 scores by true population. VDP concentrates non-NEOs near zero while assigning high probabilities to most NEOs; digest2 also scores many NEOs highly but shows broader overlap with distant populations.
```

2. ROC / threshold performance:

```text
Completeness-contamination tradeoff for VDP and digest2. At comparable contamination, VDP recovers a larger fraction of true NEOs, giving the higher best-F1 operating point.
```

3. Observable-space diagnostic:

```text
VDP probability across apparent motion and magnitude. High VDP scores align with the NEO-like velocity region, while most MBAs, TNOs, and Trojans remain below the selected threshold.
```

### Notebook Written Segments Added

The notebook `neomod/digest2_comparison.ipynb` was updated with explanatory markdown segments covering:

```text
dataset/sample context
score distributions
threshold/F1 meaning
ROC interpretation
per-population breakdown
VDP score vs v_lambda, v_beta, mag_app
digest2 score vs v_lambda, v_beta, mag_app
```

The notebook was executed successfully after the edits.

Fresh notebook outputs should show:

```text
38,045 objects
N_NEO=2233 N_non-NEO=35812
Digest2 best: t=0.970 completeness=51.3% contamination=6.4% F1=0.663
```

### Digest2 Comparison Plots Chat Addendum

Main comparison figure:

```text
neomod/roc_comparison_vdp_digest2.png
```

The left panel is the completeness-contamination tradeoff. The right panel is threshold versus completeness/contamination. Read the left panel as the main scientific result and the right panel as the operating-point explanation.

Useful verbal summary:

```text
VDP's curve reaches the lower-right region more effectively: high NEO completeness while keeping contamination low.
digest2 can also be clean, but only at a very high score threshold; that threshold rejects many true NEOs.
```

Concrete selected-candidate counts at the notebook best thresholds:

```text
VDP threshold > 0.028:
  selected = 1,821
  true NEO selected = 1,713
  false positives = 108
  missed NEOs = 520
  completeness = 76.7%
  contamination = 5.9% / about 6.0%

digest2 threshold > 0.970:
  selected = 1,224
  true NEO selected = 1,146
  false positives = 78
  missed NEOs = 1,087
  completeness = 51.3%
  contamination = 6.4%
```

False-positive composition at these thresholds:

```text
VDP false positives:
  MBA = 77
  Trojans = 28
  TNO = 3

digest2 false positives:
  TNO = 78
```

This makes the plot interpretation sharper:

```text
VDP admits slightly more false positives in absolute count at its best-F1 threshold, but it recovers 567 more true NEOs than digest2 while staying at nearly the same contamination.
digest2's strict threshold suppresses MBA and Trojan false positives, but the surviving contamination is entirely from high-scoring TNOs in this synthetic stress-test sample.
```

Score-distribution plot interpretation:

```text
VDP: non-NEO populations are compressed near P_NEO = 0; true NEOs are broadly shifted high.
digest2: true NEOs often saturate near 1.0, but TNOs also receive many elevated scores, which forces the clean threshold close to 1.0.
```

Observable-space plot interpretation:

```text
VDP vs v_lambda is the cleanest diagnostic plot. It shows that the high-VDP region follows the NEO-like apparent longitudinal velocity space.
VDP vs mag_app is useful mainly to show that magnitude alone is not driving the result.
VDP vs v_beta adds context, but it is less visually decisive than v_lambda.
```

Digest2 observable-space caveat:

```text
Do not describe digest2 as using v_lambda/v_beta maps.
The digest2 plots are correlations against the same diagnostic coordinates used for VDP.
digest2 actually receives two MPC-formatted observations, a magnitude, timing, and observatory code X05.
```

Best plot set for a poster or slide:

```text
1. ROC / threshold figure: neomod/roc_comparison_vdp_digest2.png
2. Score distributions by true population from digest2_comparison.ipynb
3. Poster-friendly VDP observable plot: P_NEO vs v_lambda and mag_app
```

Suggested one-sentence caption for the main ROC figure:

```text
On the corrected Rubin/X05 opposition-patch synthetic sample, VDP recovers 76.7% of true NEOs at about 6% contamination, whereas digest2 recovers 51.3% at a similar contamination level.
```

### Final Takeaway

The corrected comparison uses Rubin/X05 and the same 30-degree opposition-centered sky patch as the VDP map. Under that setup, VDP is cleaner and more complete than digest2 for this synthetic population test.

The clearest one-sentence result is:

```text
VDP recovers about 77% of true NEOs at about 6% contamination, while digest2 recovers about 51% at a similar contamination level.
```

### How `run_digest2_comparison.py` Works

`run_digest2_comparison.py` is a self-contained pipeline that scores the same set of synthetic S3M objects with both VDP and digest2, then plots a side-by-side ROC comparison. It is a standalone script designed to be run once, producing a parquet file of scored objects and a PNG figure.

Step-by-step:

1. Loads precomputed VDP probability maps:

```python
prob_map_set = vdp.ProbMapSet.from_npz(PROB_MAPS_PATH)
```

2. Loads and propagates each S3M population, optionally subsampling MBA and Trojans for speed, then calls `prob_map_set.score_orbital_df()` to propagate to `2025-03-21`, compute apparent magnitude, apply the 30 degree sky patch, convert to `(v_lambda, v_beta)`, and score from the correct magnitude-bin map.

3. Filters to the VDP operating region: within 30 degrees of the map centre, magnitude 14-26, and velocity-grid bounds `[-0.8, 0.8]` deg/day.

4. Builds MPC 80-column tracklets with two synthetic observations separated by 30 minutes and observatory code X05.

5. Runs digest2 in 5,000-tracklet chunks using config:

```text
noheadings
norms
NEO
```

6. Computes threshold sweeps for completeness, contamination, purity, and F1.

7. Saves:

```text
s3m_digest2_comparison_vdp_input.parquet
s3m_digest2_comparison.parquet
roc_comparison_vdp_digest2.png
```

Key design choices:

| Choice | Reason |
|---|---|
| `OBSCODE = "X05"` | Rubin observatory, matching the VDP map geometry |
| 30-min tracklet baseline | Mimics realistic Rubin nightly pair cadence |
| Chunked digest2 calls (5k each) | Prevents subprocess timeout on large inputs |
| Temporary config files | Avoids writing into the external digest2 checkout |
| Reuses `score_orbital_df` for VDP | Guarantees identical sky cut and propagation logic as the map-build step |
| `v_sky >= 0.2` cut removed | That cut was a VDP support-mask concept, not a valid pre-filter for digest2 |

---

## Current task

Advisor suggested smoothing the granular texture in velocity density maps before deriving probability maps. We chose Option B: build a per-pixel support/count map and use that as the mask for Gaussian smoothing, instead of thresholding directly on probability or density value.

## Decision

- Smooth density maps, not probability maps.
- Probability maps continue to be computed as `density_pop / sum(density_all_pops)`, so density smoothing automatically propagates into VDP scoring and plotting.
- Use raw cloned-point support counts: histogram the cloned visible points onto the velocity grid and smooth density pixels where the raw cloned support is at least `support_threshold`. The first attempt divided by `clone_factor`, but that made the threshold far too conservative for NEO maps; many NEO bins had zero pixels passing `>=3`, so they were effectively unsmoothed.
- Default smoothing config added:
  - `DEFAULT_SMOOTH_DENSITY_MAPS = True`
  - `DEFAULT_SMOOTH_SUPPORT_THRESHOLD = 3.0`
  - `DEFAULT_SMOOTH_SIGMA_PIXELS = 1.5`
  - `DEFAULT_SMOOTH_TRUNCATE_SIGMA = 5.0`

## Implementation notes

- Main file changed: `neomod/velocity_density_pipeline.py`.
- New helpers:
  - `_grid_center_edges`
  - `make_support_count_map`
  - `smooth_density_map_by_support`
- Smoothing is applied inside `build_cloned_maps_for_center_magbin` after the kNN cloned density is downweighted by `clone_factor`.
- Smoothing uses a normalized Gaussian convolution over the local nonzero-support footprint, then replaces only pixels where `support_count_map >= support_threshold`. This preserves unsupported/sharp boundaries better than convolving the full map with zeros, while still letting high-support pixels average with nearby sampled neighbors.
- Second correction after preview: the smoothing mask now uses Gaussian-weighted local support around each pixel, not exact-bin support only. Exact-bin support left too much kNN texture untouched. For classifier testing, use the modest linear-density setting `sigma=1.5`, `truncate=5`, `support_threshold=3`; stronger log-space smoothing is considered figure/display-only and should not be fed into ROC comparisons unless deliberately testing that separately.
- `.npz` save/load now includes:
  - `support_count__POP__BIN`
  - smoothing metadata: `smooth_density_maps`, `smooth_support_threshold`, `smooth_sigma_pixels`, `smooth_truncate_sigma`
- `ProbMapSet` stores this metadata on `prob_mapset.smoothing`.
- Existing older `.npz` files remain loadable; missing support maps default to zero and smoothing metadata defaults to disabled/unknown.

## Figure notebook

- `neomod/paper_figures_ecliptic_strip.ipynb` already plots densities from `prob_mapset.results[label]["density_maps_downweighted_raw"]` and probabilities from `prob_mapset.get_probability_map(...)`.
- Therefore, once `prob_maps_2025-03-21.npz` is regenerated, the density and probability figures in that notebook will reflect the smoothed maps automatically.
- Added a print line for `prob_mapset.smoothing` so the notebook shows whether the loaded map artifact contains smoothing metadata.
- Added `neomod/smoothing_preview.ipynb` so smoothing parameters can be previewed without rerunning the expensive map-generation notebook. It loads the existing `.npz`, applies the current convolution helper to copies of the loaded density maps, and plots before/after density and probability panels. For the interim May 6 `.npz`, it can rescale saved support maps by clone factor because that file used the first original-equivalent support convention.

## Verification so far

- Ran `python -m py_compile neomod/velocity_density_pipeline.py`; it passed.

## Follow-up reruns needed

After regenerating the `.npz`, rerun downstream artifacts that depend on the probability maps:

- `vdp.generate_probability_maps(...)` / `download_prob_maps.ipynb` to recreate `prob_maps_2025-03-21.npz`.
- `paper_figures_ecliptic_strip.ipynb` to refresh density/probability figures.
- ROC-related notebooks/scripts, including `run_digest2_comparison.py`, because VDP probabilities will change.
- Any parquet checkpoints derived from VDP scoring, e.g. `s3m_digest2_comparison_vdp_input.parquet` and later comparison parquet outputs.

## Caveat

I cannot literally detect when the conversation reaches 250k tokens, so this file was created proactively to preserve the context.

---

## 2026-05-07 Addendum

### Why digest2 did not really "change"

- The added convolution only changes the VDP density maps and therefore the derived VDP probabilities.
- digest2 does not read the VDP maps and does not use this convolution.
- The small digest2 difference seen after the rerun is best interpreted as a rerun/evaluation-level shift, not a classifier-method change.

Relevant numbers:

```text
Earlier corrected digest2 best:
  threshold ~ 0.970
  completeness ~ 51.3%
  contamination ~ 6.4%
  F1 ~ 0.663

2026-05-07 rerun digest2 best:
  threshold = 0.970
  completeness = 51.5%
  contamination = 6.2%
  F1 = 0.665
```

This is only a handful of objects:

```text
51.3% of 2233 NEOs  ≈ 1146 selected NEOs
51.5% of 2233 NEOs  ≈ 1150 selected NEOs
```

So the digest2 shift is on the order of about 4 true NEOs, with a similarly tiny false-positive change. That is effectively stable compared with the much larger VDP-side movement.

Most likely interpretation:

- digest2 scores themselves are effectively unchanged in method.
- the displayed optimum moved slightly because the evaluation was recomputed from a fresh rerun, with the usual tiny bookkeeping sensitivity from threshold sweeps and the known one duplicate / one missing digest2 ID warning.

### New VDP ROC result after smoothing

`run_digest2_comparison.py` completed successfully on the smoothed-map `.npz`.

```text
VDP optimal:
  threshold = 0.005
  completeness = 78.3%
  contamination = 5.7%
  F1 = 0.856

Digest2 optimal:
  threshold = 0.970
  completeness = 51.5%
  contamination = 6.2%
  F1 = 0.665
```

Compared with the earlier corrected VDP result:

```text
old VDP F1 ~ 0.845
new VDP F1 ~ 0.856
```

### Current rerun state

Artifacts successfully refreshed so far:

```text
neomod/prob_maps_2025-03-21.npz
neomod/s3m_digest2_comparison_vdp_input.parquet
neomod/s3m_digest2_comparison.parquet
neomod/roc_comparison_vdp_digest2.png
```

Still should rerun for final figure consistency:

```text
neomod/paper_figures_ecliptic_strip.ipynb
neomod/ROCcurve.ipynb
```

`neomod/digest2_comparison.ipynb` has since been rerun successfully on
2026-05-07 after the new diagnostic sections and NEOCP published-score context
were added. It no longer belongs on the outstanding-rerun list.

## 2026-05-09 VDP Production-Setting Updates

File changed:

```text
neomod/velocity_density_pipeline.py
```

Current user-side changes that are now part of the working VDP configuration:

```python
DEFAULT_K_MAP = 10
```

The final apparent-magnitude bin was tightened:

```python
{"label": "mag24+", "mag_min": 24.0, "mag_max": 25.0}
```

So the default generated-map magnitude coverage is now:

```text
14 <= mag_app < 25
```

instead of the previous:

```text
14 <= mag_app < 26
```

New Gaussian smoothing decision from `smoothing_preview.ipynb`:

```python
DEFAULT_SMOOTH_DENSITY_MAPS = True
DEFAULT_SMOOTH_POPULATION_NAMES = ("NEO",)
DEFAULT_SMOOTH_SUPPORT_THRESHOLD = 10.0
DEFAULT_SMOOTH_SIGMA_PIXELS = 4.0
DEFAULT_SMOOTH_TRUNCATE_SIGMA = 5.0
```

Important behavioral change:

```text
Gaussian smoothing is now applied only to NEO density maps.
MBA, TNO, and Trojans keep their unsmoothed downweighted density maps.
```

The probability maps are therefore computed from:

```text
smoothed NEO density + unsmoothed MBA/TNO/Trojan densities
```

The smoothing implementation itself is still the support-masked helper:

```python
smooth_density_map_by_support(...)
```

Only the default parameters and population-selection gate changed.

The production smoothing gate is in `build_cloned_maps_for_center_magbin`:

```python
if smooth_density_maps and pop_name in smooth_population_names:
    density_downweighted_map = smooth_density_map_by_support(...)
```

Metadata update:

```text
Newly generated .npz map files now save smooth_population_names.
Older .npz files without that field still load, with an empty population-name
metadata list.
```

Verification performed:

```bash
cd /Users/devanshisingh/Downloads/research/NEO_probability/neomod
python -m py_compile velocity_density_pipeline.py
```

The compile passed. The shell still prints the unrelated local conda/profile
warning:

```text
CondaError: Run 'conda init' before 'conda activate'
```

### Regenerating Maps For A New Observation Time

To generate maps for a different epoch, edit the call in:

```text
neomod/download_prob_maps.ipynb
```

and change both `obstime_str` and `output_path`, e.g.:

```python
vdp.generate_probability_maps(
    obstime_str="2026-05-09T00:00:00",
    output_path="prob_maps_2026-05-09.npz",
)
```

The output path should change so the old 2025-03-21 maps are not overwritten
accidentally.

Any downstream notebook or script must then load the new `.npz` file. The VDP
maps are tied to:

```text
obstime_str
sky center
sky radius
magnitude bins
density-estimator settings
smoothing settings
```

### Rerun Implication

The previously recorded smoothed-map ROC numbers in this context file came from
an earlier smoothing configuration. After the 2026-05-09 changes, the maps
should be regenerated before treating any VDP-vs-digest2 metric as final for the
new configuration.

---

## 2026-05-09 NEOCP Collector Migration To Arnor

### Mac Cron Job Disabled

The Mac cron job was already disabled (line commented out in crontab). No action
needed on the Mac side.

### Arnor Setup

The project folder on Arnor was renamed from `digest3` to `vdp`. The home
directory on Arnor is `/astro/users/ds2004/`, accessible via JupyterHub at
`arnor.astro.washington.edu/jupyter/user/ds2004/lab`.

The NEOCP collector now lives at:

```text
~/vdp/neocp_data/collect_neocp.py
~/vdp/neocp_data/run_neocp_loop_screen.sh
```

Both files use dynamic path resolution (`Path(__file__).resolve().parent` and
`$(dirname "${BASH_SOURCE[0]}")`) so the folder name does not matter.

### Screen Session

The collector runs in a persistent screen session on Arnor, started from the
JupyterHub terminal:

```bash
chmod +x ~/vdp/neocp_data/run_neocp_loop_screen.sh
screen -S neocp
~/vdp/neocp_data/run_neocp_loop_screen.sh
# Ctrl+A, D to detach
```

Session name: `neocp`. Check with `screen -ls`. Reattach with `screen -r neocp`.

The loop script logs to:

```text
~/vdp/neocp_data/logs/neocp_collect.log
```

Default interval is 4 hours (14400s). Can be overridden without editing the
file by setting `NEOCP_INTERVAL_SEC` before running:

```bash
NEOCP_INTERVAL_SEC=7200 ~/vdp/neocp_data/run_neocp_loop_screen.sh  # 2 hours
NEOCP_INTERVAL_SEC=3600 ~/vdp/neocp_data/run_neocp_loop_screen.sh  # 1 hour
```

To change the interval on a running session: `screen -r neocp`, Ctrl+C, then
restart with the new env var.

### First Successful Run On Arnor

Test run output:

```text
20260509T230955Z: status=ok, objects=26, active=24, ephemeris_rows=25
```

First full loop run output:

```text
20260509T231155Z: status=ok, objects=26, active=24, ephemeris_rows=600
```

Data directories created automatically under `~/vdp/neocp_data/`:

```text
ephemerides/   logs/   raw/   snapshots/   tables/
```

### Data Model: Append-Only

`neocp_objects_history.csv` is append-only. Each collection run adds new rows
with a new `snapshot_id` timestamp. The same designation appears once per run.

This means object score, nobs, arc_days, and not_seen_days can be tracked over
time:

```python
objs[objs["designation"] == "ST26E14"][
    ["snapshot_id", "fetched_at_utc", "score", "nobs", "arc_days"]
]
```

### Monitoring Notebook

Created `~/vdp/check_neocp.ipynb` (lives in `~/vdp/`, one level above
`neocp_data/`). Re-run any time to see latest data. Four cells:

```text
1. Recent runs (neocp_runs.csv) — status and counts per run
2. Latest snapshot objects sorted by MPC score
3. Collection stats — snapshots, total rows, unique designations
4. Score histogram + sky positions scatter plot
```

### What The Collector Captures For VDP Comparison

The collector saves enough to VDP-score real NEOCP objects using
`score_observation` (not `score_orbital_df`, since we have no orbital elements
for real objects).

Fields available:

```text
ra_deg, dec_deg          — sky position
motion_arcsec_per_min    — total angular speed on sky
pa_deg                   — position angle of motion (N through E)
v_mag                    — apparent magnitude
ephem_time_utc           — epoch (determines which VDP map to use)
score / score_0_1        — MPC-published NEOCP score (the digest2-style score)
```

### Correct Decomposition Of Motion Into dRA/dDec

`score_observation` expects coordinate rates (dα/dt, dδ/dt), not angular rates.
MPC's `motion * sin(PA)` is the angular rate in the RA direction, which equals
`dα/dt * cos(δ)`. So the cos(Dec) factor is required for dRA:

```python
arcsec_per_min_to_deg_per_day = 60 * 24 / 3600  # = 0.4

ddec_deg_day = (motion_arcsec_per_min
                * np.cos(np.deg2rad(pa_deg))
                * arcsec_per_min_to_deg_per_day)

dra_deg_day = (motion_arcsec_per_min
               * np.sin(np.deg2rad(pa_deg))
               / np.cos(np.deg2rad(dec_deg))
               * arcsec_per_min_to_deg_per_day)
```

Omitting the `/ cos(dec)` factor is wrong and gives incorrect ecliptic
velocities, especially at high declination.

### No Need To Synthesize Digest2 Tracklets

The MPC score is already collected as `score_0_1`. There is no need to build
MPC 80-column tracklets or run digest2 on these objects. The comparison is:

```text
VDP P_NEO  (computed via score_observation)
vs.
MPC score_0_1  (already in the collected CSV)
```

### Downstream Steps For Live Comparison

1. Generate VDP maps for the relevant epoch — maps must match the epoch of the
   collected ephemerides (e.g. `prob_maps_2026-05-09.npz`). Use
   `download_prob_maps.ipynb` with the matching `obstime_str`.

2. For each collected object, apply the motion decomposition above to get
   `dra_deg_day` and `ddec_deg_day`.

3. Call `pms.score_observation(ra_deg, dec_deg, dra_deg_day, ddec_deg_day,
   mag_app)` to get `P_NEO_vdp`.

4. Compare `P_NEO_vdp` against `score_0_1` from the CSV.

Note: the VDP map sky center (opposition patch) must overlap with where the
NEOCP objects actually are. Objects outside the map's 30-degree sky patch will
return P=0 and should be excluded from the comparison.

---

## 2026-05-10 VDP Two-Stage NEO Smoothing Update

File changed:

```text
neomod/velocity_density_pipeline.py
```

Motivation:

```text
The sigma=4 NEO panel in smoothing_preview.ipynb looked smoother than the first
freshly regenerated production maps. The reason is that smoothing_preview.ipynb
was loading an older .npz that already had light smoothing baked in
(threshold=3, sigma=1.5), then applying the preview sigma=4, threshold=10 pass
on top of that loaded density.
```

Decision:

```text
Production VDP should intentionally reproduce that level of NEO smoothing.
```

Current default production smoothing is now:

```python
DEFAULT_SMOOTH_DENSITY_MAPS = True
DEFAULT_SMOOTH_POPULATION_NAMES = ("NEO",)
DEFAULT_SMOOTH_PRESMOOTHING_PASSES = (
    {"support_threshold": 3.0, "sigma_pixels": 1.5, "truncate_sigma": 5.0},
)
DEFAULT_SMOOTH_SUPPORT_THRESHOLD = 10.0
DEFAULT_SMOOTH_SIGMA_PIXELS = 4.0
DEFAULT_SMOOTH_TRUNCATE_SIGMA = 5.0
```

So every new `vdp.generate_probability_maps(...)` run applies:

```text
NEO only:
  pass 1: support_threshold=3.0,  sigma_pixels=1.5
  pass 2: support_threshold=10.0, sigma_pixels=4.0

MBA/TNO/Trojans:
  no Gaussian smoothing
```

The actual smoothing helper is still:

```python
smooth_density_map_by_support(...)
```

The code now normalizes smoothing passes through:

```python
_normalize_smoothing_passes(...)
```

and the production loop applies each pass in order:

```python
if smooth_density_maps and pop_name in smooth_population_names:
    for pass_cfg in smooth_passes:
        density_downweighted_map = smooth_density_map_by_support(...)
```

New `.npz` metadata now records the full ordered pass list:

```text
smooth_pass_support_thresholds
smooth_pass_sigma_pixels
smooth_pass_truncate_sigmas
```

and `ProbMapSet.from_npz(...).smoothing` exposes a `passes` list, e.g.:

```python
[
    {"support_threshold": 3.0, "sigma_pixels": 1.5, "truncate_sigma": 5.0},
    {"support_threshold": 10.0, "sigma_pixels": 4.0, "truncate_sigma": 5.0},
]
```

Important:

```text
Existing .npz files are not modified. Any map file must be regenerated after
this code change to get the two-stage NEO smoothing baked in.
```

Example regeneration call for the NEOCP May 9 map:

```python
import velocity_density_pipeline as vdp

vdp.generate_probability_maps(
    obstime_str="2026-05-09T22:00:00",
    output_path="prob_maps_2026-05-09T22_neocp.npz",
    center_lon_deg=229.0,
    center_lat_deg=0.0,
)
```

Because these are now defaults, no explicit smoothing arguments are needed in
`download_prob_maps.ipynb` unless deliberately overriding the production
settings.

Verification performed after code edit:

```bash
cd /Users/devanshisingh/Downloads/research/NEO_probability/neomod
python -m py_compile velocity_density_pipeline.py
```

The compile passed. The recurring local shell warning is unrelated:

```text
CondaError: Run 'conda init' before 'conda activate'
```

Also created a comparison notebook:

```text
neomod/sigma4_map_comparison.ipynb
```

Purpose:

```text
Compare the old smoothing-preview sigma=4 target against a newly generated .npz
using paper_figures_ecliptic_strip-style density/probability grids before
pointing paper figures at the new map file.
```

---

## 2026-05-11 Latest VDP / Notebook Context

### Support-Scaled Smoothing Fix

File changed:

```text
neomod/velocity_density_pipeline.py
```

Problem found:

```text
smoothing_preview.ipynb multiplied the NEO support map by the NEO clone factor
(300) before applying the support threshold. Production VDP was using raw
support counts directly. Therefore only the tiny high-support core was smoothed
in generated maps, while the preview smoothed a much broader NEO footprint.
```

Current production behavior now includes:

```python
DEFAULT_SMOOTH_SUPPORT_SCALE_BY_CLONE_FACTOR = True
```

During smoothing, VDP now uses:

```python
support_for_smoothing = support_count_map * clone_factor
```

For NEO maps this means:

```text
support_for_smoothing = support_count_map * 300
```

The final production smoothing behavior is:

```text
NEO only:
  pass 1: support_threshold=3.0,  sigma_pixels=1.5
  pass 2: support_threshold=10.0, sigma_pixels=4.0
  support thresholds evaluated on support_count_map * clone_factor

MBA/TNO/Trojans:
  no Gaussian smoothing
```

New `.npz` metadata now includes:

```text
smooth_support_scale_by_clone_factor
```

and `ProbMapSet.from_npz(...).smoothing` should show:

```python
"support_scale_by_clone_factor": True
```

Verification after the code edit:

```bash
cd /Users/devanshisingh/Downloads/research/NEO_probability/neomod
python -m py_compile velocity_density_pipeline.py
```

Compile passed. The recurring shell warning remains unrelated:

```text
CondaError: Run 'conda init' before 'conda activate'
```

Important:

```text
Any .npz generated before the support-scaling fix must be regenerated. The code
change does not modify existing map files.
```

### Poster/Figure Map Regeneration

The user is regenerating `download_prob_maps.ipynb` with:

```text
center_lon_deg=180.0
center_lat_deg=0.0
```

Reason:

```text
Poster/paper figures should keep the original center=(180,0) geometry but use
the newer support-scaled, two-stage NEO smoothing.
```

This is distinct from the NEOCP test map:

```text
prob_maps_2026-05-09T22_neocp.npz
center=(229,0)
```

For poster replacement figures, use a regenerated center=(180,0) map rather
than the NEOCP center=(229,0) map.

### Paper Figures Notebook State

File touched earlier:

```text
neomod/paper_figures_ecliptic_strip.ipynb
```

It was temporarily updated to load:

```python
PROB_MAPS_PATH = "prob_maps_2026-05-09T22_neocp.npz"
```

For final poster/paper replacement plots, this should be changed to the newly
regenerated center=(180,0) `.npz` file once it exists.

The plotted density maps use:

```python
prob_mapset.results[label]["density_maps_downweighted_raw"][pop]
```

Despite the historical name, for newly generated files this array contains the
saved downweighted density after any production smoothing has been applied.

### Cloning Test Notebook

Created:

```text
neomod/cloning_test_ZI.ipynb
```

Purpose:

```text
Trace S3M orbital elements through conditional K|M cloning, sky/velocity
projection, v_lambda/v_beta scatter, density-map construction, and smoothing
diagnostics.
```

Key notebook objects:

```text
source_df_sample       sampled source orbital-element table
source_df_with_mag     sample with mag_app and mag_bin_label
source_df              actual input to cloning
raw_clone_df           cloned orbital elements before projection
clone_visible_df       cloned orbital elements plus ra/dec/lam/beta/vlam/vbeta
```

The notebook includes:

```text
- loading all four S3M DataFrames
- accessing orbital-element columns
- apparent-magnitude bin split diagnostics
- optional one-bin cloning via USE_MAG_BIN_FOR_CLONING
- direct clone_population_conditional_K_from_M usage
- build_visible_subset_dataframe projection
- v_lambda/v_beta scatter and binned support
- density map construction
- production-style smoothing and smoothed-raw delta plot
- examples for selecting objects by velocity-space region and inspecting their
  orbital elements
```

Important default:

```python
USE_MAG_BIN_FOR_CLONING = False
```

Reason:

```text
If USE_MAG_BIN_FOR_CLONING=True with a small MAX_SOURCE_OBJECTS sample, the
selected mag bin can be too sparse and the velocity plots look nearly empty.
The mag-bin section is therefore an inspection diagnostic by default, while the
main cloning test uses the full sampled source table.
```

To test a single mag bin, set:

```python
USE_MAG_BIN_FOR_CLONING = True
```

and usually increase:

```python
MAX_SOURCE_OBJECTS
CLONE_FACTOR_TEST
```

### Desired Repository Reorganization

The user wants to reorganize `neomod/` so active code, notebooks, generated
outputs, data products, old notebooks, and paper/poster assets are clearly
separated. This will require moving files and then updating import paths,
relative file paths, notebook load/save paths, and script paths throughout the
project.
