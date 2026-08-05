# Report feature-group EIC mirror plots as PDFs

For each polarity with MS1 feature groups, generate
[`plot_xcms_feature_group_EIC_comparasion`](https://drruili.github.io/MSdev/reference/plot_xcms_feature_group_EIC_comparasion.md)
for every group with at least two features and append pages into
`Feature_group_EIC_Positive.pdf` / `Feature_group_EIC_Negative.pdf`
under `projectDir`. Reuses stored chromatograms when available
(`Positive_Chromatograms` / `Negative_Chromatograms`).

## Usage

``` r
Report_MSdev_feature_group_EIC(
  object,
  expandRt = 2,
  min_width = 20,
  max_features = NULL,
  sample_index = 1L,
  min_features = 2L,
  width = NULL,
  height = NULL
)
```

## Arguments

- object:

  An `MSdev` object with `featureGroups` on MS1 xcms objects (e.g. after
  [`MSdev_group_feature_EIC`](https://drruili.github.io/MSdev/reference/MSdev_group_feature_EIC.md)).

- expandRt:

  numeric(1). Seconds added on each side of each feature's peak RT
  window, passed to `plot_xcms_feature_group_EIC_comparasion` (default
  `2`).

- min_width:

  numeric(1). Minimum RT window width (seconds), passed to
  `plot_xcms_feature_group_EIC_comparasion` (default `20`).

- max_features:

  Maximum features per group plot; `NULL` keeps all (default).

- sample_index:

  Sample index into stored `EIC_Similarity` (default `1L`).

- min_features:

  Minimum group size to plot (default `2L`).

- width, height:

  PDF page size in inches. If `NULL`, scaled from the number of features
  in each group.

## Value

Invisible character vector of written PDF path(s).
