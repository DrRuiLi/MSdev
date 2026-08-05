# Mirror EIC comparison grid for feature groups

Pairwise chromatogram comparison for features in one or more xcms
feature groups on an `XcmsExperiment` (or compatible) object. Panel
\\(i,j)\\ shows the column feature EIC above \\y = 0\\ and the row
feature EIC below (each cropped to its peak RT window via
[`XChromatograms_subset_feature`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md),
exactly as the EIC grouping does, then normalized to its own max). Cell
labels are EIC similarities from
`MsExperiment::otherData(xcms)$EIC_Similarity` when available (e.g.
after
[`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md)).
Strip backgrounds are colored by feature group when ggh4x is installed.

## Usage

``` r
plot_xcms_feature_group_EIC_comparasion(
  xcms,
  feature_group,
  chroms = NULL,
  expandRt = 2,
  min_width = 20,
  max_features = NULL,
  sample_index = 1L,
  title = NULL
)
```

## Arguments

- xcms:

  An `XcmsExperiment` / `MsExperiment` (or `XCMSnExp`) with feature
  definitions and `featureGroups`.

- feature_group:

  Character vector of feature group id(s) to compare.

- chroms:

  Optional `XChromatograms` for the features (rownames = feature ids)
  with a `featureDefinitions` slot carrying `peakRtMin`/`peakRtMax`. If
  `NULL`, extracted via
  [`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md).

- expandRt:

  numeric(1). Seconds added on each side of each feature's
  `[peakRtMin, peakRtMax]` window when cropping via
  [`XChromatograms_subset_feature`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
  (default `2`).

- min_width:

  numeric(1). Minimum RT window width (seconds) after `expandRt`;
  shorter windows are padded equally on both sides (default `20`).

- max_features:

  Maximum features kept per group (ordered by `rtmed`). `NULL` keeps all
  members (default).

- sample_index:

  Sample index into `otherData(xcms)$EIC_Similarity` when several
  samples are stored (default `1L`).

- title:

  Optional plot title. Default is
  `"Feature group: <id> (<rt center>), ..."`, where the RT center is the
  median `rtmed` of each group's member features.

## Value

A `ggplot` object.
