# Export feature EICs for all features

Write `{feature_id}.EIC.png` for each feature. Prefers stored
chromatograms from
[`get_MSdev_Chromatogram`](https://drruili.github.io/MSdev/reference/get_MSdev_Chromatogram.md)
(after
[`MSdev_get_feature_chrom`](https://drruili.github.io/MSdev/reference/MSdev_get_feature_chrom.md));
otherwise extracts with
[`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
(not
[`xcms::chromatogram`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)).
Feature IDs are resolved as in
[`MSdev_export_feature_MSMS`](https://drruili.github.io/MSdev/reference/MSdev_export_feature_MSMS.md).
Failures on individual features are warned and skipped.

## Usage

``` r
MSdev_export_feature_Chromatographs(
  object,
  out.dir = file.path(object@projectInfo$projectDir, "EIC"),
  feature_id = NULL
)
```

## Arguments

- object:

  MSdev object

- out.dir:

  Output directory. Default `object@projectInfo$projectDir/EIC`.

- feature_id:

  Optional character vector of feature IDs. Default all features
  discovered as above.

## Value

Invisible character vector of feature IDs attempted.

## See also

[`MSdev_export_feature_MSMS`](https://drruili.github.io/MSdev/reference/MSdev_export_feature_MSMS.md),
[`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md),
[`get_MSdev_Chromatogram`](https://drruili.github.io/MSdev/reference/get_MSdev_Chromatogram.md)
