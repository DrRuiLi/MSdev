# Compute MS1 purity matrix for xcms features

Calculate a feature-by-sample MS1 purity matrix by extracting, for each
feature in `xcms.xcms`, the closest MS1 scan (by retention time) from
`xcms.ms1.sp` in each sample file (matched by
[`Spectra::dataOrigin`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)).
Purity is calculated within an isolation window around the feature m/z.

If `xcms.ms1.sp` is not provided, MS1 spectra are taken via
`ProtGenerics::spectra(xcms.xcms)` and filtered to MS level 1.

Optimized path: closest-scan lookup via `findInterval`, a single
`peaksData()` extraction, and purity computed on peak matrices (no
per-feature `Spectra` subsetting).

## Usage

``` r
get_xcms_feature_purity_matrix(
  xcms.xcms,
  xcms.ms1.sp = NULL,
  ppm = 5,
  isolation_half_window = 0.2
)
```

## Arguments

- xcms.xcms:

  `XCMSnExp` with grouped features (must have `featureDefinitions`).

- xcms.ms1.sp:

  Optional MS1 `Spectra` covering the same files as `xcms.xcms`. If
  missing or `NULL`, taken from
  [`ProtGenerics::spectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html).

- ppm:

  numeric, ppm tolerance for m/z window.

- isolation_half_window:

  numeric, half isolation window (m/z).

## Value

numeric matrix with rows = `feature_id`, columns = sample files (by
`dataOrigin`).
