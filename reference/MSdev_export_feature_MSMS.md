# Export MS/MS spectrum and chromatogram for all features

Loop
[`export_MSdev_feature_MSMS`](https://drruili.github.io/MSdev/reference/export_MSdev_feature_MSMS.md)
over features. Feature IDs are taken from `advancedAna$featureRaw`, else
`feature.se`, else xcms `featureDefinitions` (`_pos`/ `_neg` suffix).
Writes `{feature_id}.MSMS.png` (experimental vs reference mirror when a
reference spectrum is present) and `{feature_id}.Chrom.png`. Failures on
individual features are warned and skipped.

## Usage

``` r
MSdev_export_feature_MSMS(
  object,
  out.dir = file.path(object@projectInfo$projectDir, "MSMS"),
  feature_id = NULL,
  cpdb_path = object@projectInfo$CompoundDB_path
)
```

## Arguments

- object:

  MSdev object

- out.dir:

  Output directory. Default `object@projectInfo$projectDir/MSMS`.

- feature_id:

  Optional character vector of feature IDs. Default all features
  discovered as above.

- cpdb_path:

  Optional path to a CompoundDb SQLite file. Default
  `object@projectInfo$CompoundDB_path` (set by
  [`MSdev_annotation`](https://drruili.github.io/MSdev/reference/MSdev_workflow.md)).

## Value

Invisible character vector of feature IDs attempted.

## See also

[`export_MSdev_feature_MSMS`](https://drruili.github.io/MSdev/reference/export_MSdev_feature_MSMS.md),
[`plot_MSdev_feature_spectrum`](https://drruili.github.io/MSdev/reference/plot_MSdev_feature_spectrum.md)
