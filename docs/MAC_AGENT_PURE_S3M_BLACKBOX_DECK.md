# Mac Agent Handoff: Pure-S3M Blackbox Notebook And Deck

## Mission

Build the final advisor-facing notebook and deck from the Hyak product bundle.

The Mac is the synthesis/presentation machine. Hyak produces or stages the full
grid products; the Mac turns them into a coherent story, tables, figures, and a
Google-Slides-targeted deck.

## Story To Preserve

The deck answers:

> When the same VDP blackbox is given pure-S3M benchmark tracklets versus
> pure-S3M Sorcha tracklets, why does performance change, and what do we need
> to fix?

The deck must both defend the benchmark result and ask for concrete help with
the Sorcha gap.

Do not turn this into a polished conference talk. It is an advisor-debug deck.

## Expected Inputs

Local known files already on the Mac:

```text
S3Mdata/*.s3m
s3m_digest2_comparison_s3m_gmm.parquet
s3m_digest2_comparison_vdp_input_s3m_gmm.parquet
outputs/phase2_s3m/sorcha_comparison_s3m.parquet
s3m_digest2_comparison_2026-05-09T22_neocp.parquet
s3m_digest2_comparison_2026-05-09_antisun_minus_45.parquet
s3m_digest2_comparison_2026-05-09_antisun_minus_90.parquet
s3m_digest2_comparison_2026-05-09_antisun_minus_120.parquet
```

Expected Hyak bundle:

```text
outputs/presentation_pure_s3m_blackbox/
  README.txt
  pure_s3m_blackbox_hyak_validation.txt
  grid_manifest_gmm.csv
  benchmark_metrics_by_direction_bin.csv
  benchmark_metrics_by_mag_bin.csv
  sorcha_metrics_by_direction_bin.csv
  sorcha_metrics_by_mag_bin.csv
  figures/
```

## Source Notebook

Create one notebook:

```text
presentation_pure_s3m_blackbox.ipynb
```

The notebook is the source of truth for all deck calculations and figures.

Required notebook sections:

1. **Experiment definition**
   - define input worlds: benchmark pure-S3M and Sorcha pure-S3M
   - define map branch: GMM 667-map grid (kNN branch dropped)
   - define direction bins
   - define metrics

2. **Load products**
   - load local benchmark/Sorcha parquets
   - load Hyak metrics/figure manifests
   - print row counts and known caveats

3. **Metric functions**
   - best-F1 threshold
   - completeness
   - contamination
   - F1
   - optional AUC

4. **Known local checkpoints**
   - Benchmark pure-S3M (GMM 667-map grid, 475k tracklets) — generated on Hyak,
     use outputs/phase2_benchmark_s3m/benchmark_comparison_s3m.parquet
   - Sorcha pure-S3M combined result should reproduce approximately:
     - VDP F1 about 0.867 (antisun bin)
     - Digest2 F1 about 0.859 (antisun bin)

5. **Tables**
   - full matrix table
   - direction-bin table
   - mag-bin table
   - benchmark-vs-Sorcha delta table

6. **Figures**
   - blackbox schematic
   - tracklet construction schematic
   - Digest2 scoring schematic
   - ROC panels
   - density/P(NEO) map panels
   - mask ON/OFF comparison panels
   - selected object rate sanity check

7. **Deck export manifest**
   - list every figure/table filename consumed by the deck

## Deck Structure

Build a `.pptx` first. Import to Google Slides only if the connector/tool is
available.

Suggested deck title:

```text
Pure-S3M VDP Blackbox Test:
Benchmark vs Sorcha Inputs
```

Slides:

1. **Title / meeting purpose**
   - "Controlled experiment to isolate why Sorcha performance changes."

2. **The one question**
   - same blackbox, different pure-S3M input worlds.

3. **Blackbox contract**
   - `tracklets + maps -> VDP -> P(NEO)`
   - `same tracklets -> Digest2 -> P(NEO)_d2`

4. **Experimental matrix**
   - benchmark vs Sorcha, both GMM 667-map grid.

5. **Tracklet construction**
   - benchmark synthetic 2-det X05 tracklets
   - Sorcha observed cadence tracklets
   - same scoring/evaluation contract after tracklets exist

6. **Map branch**
   - GMM 667 maps, ref epoch 2026-01-01
   - same sky-grid bins for both input worlds

7. **Masks**
   - nearest-distance mask
   - support-count mask
   - what mask OFF/ON means

8. **Digest2 comparison**
   - MPC80 formatting
   - same tracklets
   - Digest2 score parsing

9. **Headline result matrix**
   - full-sky metrics for benchmark vs Sorcha.

10. **Direction-bin metrics**
   - full sky plus 0-20, 20-40, 40-70, 70-110, 110-141.

11. **Magnitude-bin metrics**
   - show where the method fails or holds.

12. **Benchmark vs Sorcha ROC**
   - side-by-side, GMM branch.

13. **Density maps**
   - benchmark and Sorcha side-by-side if available.

14. **P(NEO) maps mask OFF**
   - GMM, representative direction and mag bins.

15. **P(NEO) maps mask ON**
   - GMM, representative direction and mag bins.

16. **Where the gap appears**
   - direction/mag region where Sorcha loses to Digest2.

18. **Rate sanity check**
   - selected objects versus JPL Horizons/trusted reference.

19. **Interpretation candidates**
   - tracklet/rate bug
   - cadence/selection effect
   - label/catalog mismatch
   - map-generation issue
   - real algorithmic limitation

20. **Advisor ask**
   - "The benchmark says it can work. Sorcha still does not meet expectations.
     I need help deciding what is wrong or what must change."

Use appendix slides for extra direction-bin, mag-bin, and mask panels.

## Design Rules

- Prefer dense-but-readable advisor debugging slides.
- Keep labels explicit: benchmark vs Sorcha, mask ON vs OFF.
- Use consistent colors:
  - VDP: blue
  - Digest2: orange
  - GMM branch: blue/teal accents

- Every result slide must state the input world and map branch.
- Do not hide the Sorcha underperformance.

## Acceptance Checks

Before delivery:

- Notebook runs top-to-bottom.
- Known local metrics reproduce.
- Every deck figure exists and is nonblank.
- Direction-bin labels are consistent across slides.
- Deck has no overlap/clipping.
- Final advisor ask is explicit and direct.

## Out Of Scope

- Hybrid catalog results.
- New algorithm changes.
- Paper prose rewrite.
- Claims that VDP already beats Digest2 on Sorcha if the current pure-S3M result
  does not show that.

