# Extract chromatograms for features

Extract chromatograms for specified features from xcms data via
[`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
(full RT, all samples), storing them as on-disk data.

## Usage

``` r
MSdev_get_feature_chrom(
  object,
  BPPARAM = SnowParam(workers = max(1L, floor(parallel::detectCores()/2)), progressbar =
    T),
  feature.list = NULL,
  expandMzppm = 2
)
```

## Arguments

- object:

  MSdev object

- BPPARAM:

  BiocParallel backend for parallel processing

- feature.list:

  optional list of feature IDs with names "Positive" and "Negative"

- expandMzppm:

  numeric(1). Passed to
  [`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md):
  extra m/z pad in ppm (`mzmin * (1 - ppm/1e6)`,
  `mzmax * (1 + ppm/1e6)`). Default `2`.

## Value

MSdev object with chromatograms stored
