# Pure-S3M Blackbox Advisor Deck

## Purpose

This presentation is for Mario and Zeljko after the June advisor meeting where the
main critique was that the benchmark/Sorcha story skipped too many middle steps.

The deck should answer one question:

> When the same VDP blackbox is given pure-S3M benchmark tracklets versus
> pure-S3M Sorcha tracklets, why does performance change, and what do we need
> to fix?

The tone is both:

- **Results defense:** the controlled benchmark shows the classifier can work.
- **Debugging ask:** the Sorcha pure-S3M result still does not meet the expected
  standard of outperforming Digest2, so advisors need to help identify the bug,
  mismatch, or algorithmic change.

This is an advisor-debug deck, not a polished conference talk.

## Core Experiment

There are two input worlds:

| Input world | Meaning |
| --- | --- |
| Benchmark pure-S3M | Pure S3M objects, propagated to a map epoch, made into synthetic two-detection Rubin/X05 tracklets. No Sorcha cadence. |
| Sorcha pure-S3M | Pure S3M objects observed through Sorcha/Rubin cadence, then postprocessed into tracklets. |

Both go into the same evaluation contract:

```text
tracklets + probability maps -> VDP blackbox -> P(NEO)_vdp
same tracklets              -> Digest2       -> P(NEO)_d2
```

The key point is that benchmark-vs-Sorcha is a controlled input comparison,
not a hand-wavy "Tom Wagg did it" comparison.

## Map Branch

All scoring uses the GMM cloning pipeline (667 antisun-relative maps, ref epoch
2026-01-01). kNN/K|M maps are not included in this deck.

The controlled matrix is:

| Input type | Map branch | Comparator |
| --- | --- | --- |
| Benchmark pure-S3M | GMM 667 maps | Digest2 |
| Sorcha pure-S3M | GMM 667 maps | Digest2 |

Benchmark versus Sorcha is the single controlled axis.

## Standard Bins

Every result table and ROC panel should use the same directional structure:

- full sky
- 0-20 deg from antisun
- 20-40 deg from antisun
- 40-70 deg from antisun
- 70-110 deg from antisun
- 110-141 deg from antisun

Also include magnitude-bin versions where available.

Do not mix exact old offsets like -45 deg with the new binned story unless they
are clearly labelled as historical diagnostics.

## Required Figures

For each major comparison, show benchmark and Sorcha side by side where possible:

- log-density maps in ecliptic rate space
- P(NEO) maps
- mask OFF panels
- mask ON panels
- direction-bin ROC curves
- magnitude-bin ROC curves
- summary results tables with completeness, contamination, F1, and threshold
- tracklet-building schematic for benchmark versus Sorcha
- Digest2 input/score construction slide
- object-level sanity check comparing ecliptic rates to JPL Horizons for selected objects

The main deck can be dense because it is for advisor debugging. Put extra panels
after the main story so Devanshi can choose live what to show.

## Key Slide Story

1. **Question.** Why does performance change when the input becomes Sorcha?
2. **Blackbox.** Same VDP scoring interface; same Digest2 comparator.
3. **Inputs.** Benchmark synthetic tracklets versus Sorcha cadence tracklets.
4. **Map branch.** GMM 667-map grid, ref epoch 2026-01-01.
5. **Masks.** Explain nearest-distance mask versus support-count mask, and why mask ON/OFF panels matter.
6. **Results matrix.** Full-sky and direction-bin metrics for each input/map branch.
7. **Diagnostics.** Show rate-space density/P(NEO) panels and ROC curves.
8. **Sanity checks.** Tracklet/rate validation, including selected JPL Horizons comparisons.
9. **Advisor ask.** The benchmark says the classifier can work; Sorcha still does not satisfy expectations. Please help decide what is wrong or what must change.

## Final Advisor Ask

Say this explicitly:

```text
The controlled benchmark says the classifier can work.
The Sorcha pure-S3M result says the current model still does not meet the
expectation of outperforming Digest2.

I need your help deciding whether the remaining gap is:
1. a bug in tracklet construction or rate conversion,
2. a benchmark-vs-Sorcha cadence/selection effect,
3. a label/catalog mismatch,
4. a map-generation issue,
5. or a real algorithmic limitation requiring changes to VDP.
```

## Scope Boundary

This deck is pure-S3M only.

Hybrid catalog tests come after the pure-S3M controlled experiment is understood.

