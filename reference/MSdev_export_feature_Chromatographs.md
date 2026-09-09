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
Set `re_extract = TRUE` to ignore stored chromatograms and extract with
the arguments below. Feature IDs are resolved as in
[`MSdev_export_feature_MSMS`](https://drruili.github.io/MSdev/reference/MSdev_export_feature_MSMS.md).
Failures on individual features are warned and skipped.

## Usage

``` r
MSdev_export_feature_Chromatographs(
  object,
  out.dir = file.path(object@projectInfo$projectDir, "EIC"),
  feature_id = NULL,
  re_extract = FALSE,
  selected_sample = NULL,
  rt = c("expand", "identity", "all"),
  expandRt = 15,
  mz.expand = 0,
  expandMzppm = 2,
  aggregationFun = "max",
  attachPeaks = FALSE,
  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
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

- re_extract:

  logical(1). If `TRUE`, skip stored chromatograms and extract with
  [`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md).
  Default `FALSE`. Extract arguments (`rt`, `expandRt`, `mz.expand`,
  `expandMzppm`, `aggregationFun`) apply only when extracting; they have
  no effect on already-stored chromatograms unless `re_extract = TRUE`.

- selected_sample:

  Sample selection, passed to
  [`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
  and used to subset stored chromatogram columns. `NULL` (default)
  overlays one sample per group (or up to five samples). `"maxo"` uses
  the highest-intensity sample; `"all"` uses all samples; integer
  indices or sample name(s) select those samples.

- rt:

  one of `c("expand", "identity", "all")`. Passed to
  [`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
  when extracting. Default `"expand"`.

- expandRt:

  seconds added each side when `rt = "expand"`. Default `15`.

- mz.expand:

  fraction of mz width to expand on each side. Default `0`.

- expandMzppm:

  extra m/z pad in ppm applied after `mz.expand`. Default `2`.

- aggregationFun:

  `"max"` or `"sum"`, passed to
  [`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md).
  Default `"max"`.

- attachPeaks:

  logical; attach feature chromPeaks when extracting. Default `FALSE`.

- BPPARAM:

  BiocParallel backend for extraction. Default
  `SerialParam(progressbar = FALSE)`.

## Value

Invisible character vector of feature IDs attempted.

## See also

[`MSdev_export_feature_MSMS`](https://drruili.github.io/MSdev/reference/MSdev_export_feature_MSMS.md),
[`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md),
[`get_MSdev_Chromatogram`](https://drruili.github.io/MSdev/reference/get_MSdev_Chromatogram.md)
