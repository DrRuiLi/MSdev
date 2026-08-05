# Compute MS1 feature purity for both polarities

Per-polarity wrapper around
[`xcms_get_feature_purity`](https://drruili.github.io/MSdev/reference/xcms_get_feature_purity.md).
For each available MS1 `XcmsExperiment` (`PositiveMS1` / `NegativeMS1`),
computes the feature-by-sample purity matrix and stores it as assay
`"purity_matrix"` in `qdata`.

Skips a polarity when `purity_matrix` already exists unless
`force = TRUE`. Failures are caught per polarity so one missing xcms
slot does not abort the other.

## Usage

``` r
MSdev_get_feature_purity(
  object,
  ppm = 10,
  isolation_half_window = 0.2,
  force = FALSE,
  drop_blank = TRUE,
  ...
)
```

## Arguments

- object:

  MSdev object

- ppm:

  numeric, ppm tolerance for m/z window. Default 10.

- isolation_half_window:

  numeric, half isolation window (m/z). Default 0.2.

- force:

  logical; if `FALSE` (default), skip a polarity when `purity_matrix`
  already exists in qdata. If `TRUE`, recompute.

- drop_blank:

  logical; retained for API compatibility (unused; blanks are not
  stripped from qdata).

- ...:

  ignored (retained for API compatibility).

## Value

MSdev object with updated MS1 xcms `qdata` assays.

## See also

[`xcms_get_feature_purity`](https://drruili.github.io/MSdev/reference/xcms_get_feature_purity.md)
