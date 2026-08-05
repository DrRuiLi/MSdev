# Match Ions by M/Z and Retention Time

Matches two lists of ions based on mass-to-charge ratio (m/z) and
retention time (RT). Returns all matched pairs including multiple
matches, with error values for each match.

## Usage

``` r
match_mz_rt(
  mz1,
  rt1 = rep(NA, length(mz1)),
  mz2,
  rt2 = rep(NA, length(mz2)),
  mz.ppm = 10,
  rt.tol = Inf
)
```

## Arguments

- mz1:

  Numeric vector of m/z values for the first ion list.

- rt1:

  Numeric vector of retention times for the first ion list. Default is
  NA values of same length as mz1.

- mz2:

  Numeric vector of m/z values for the second ion list.

- rt2:

  Numeric vector of retention times for the second ion list. Default is
  NA values of same length as mz2.

- mz.ppm:

  Numeric tolerance for m/z matching in parts per million (ppm). Default
  is 10.

- rt.tol:

  Numeric tolerance for retention time matching. Default is Inf (no RT
  filtering).

## Value

A `data.table` with columns: `ion1` (index in first list), `ion2` (index
in second list), `mz.error` (relative m/z error), and `rt.error`
(absolute RT difference).

## Details

Matching is done in two steps for efficiency:

1.  [`match_mz_foverlaps`](https://drruili.github.io/MSdev/reference/match_mz_foverlaps.md)
    finds m/z candidate pairs within `mz.ppm` using
    [`data.table::foverlaps()`](https://rdrr.io/pkg/data.table/man/foverlaps.html)
    on ppm-expanded intervals (denominator `mz1`).

2.  In `data.table`, absolute RT differences are computed and rows with
    `rt.error < rt.tol` (or `NA` RT error) are kept.

`mz.error` is the relative m/z error `abs(mz1 - mz2) / mz1` (i.e.
`abs(mz.ppm) / 1e6` from `match_mz_foverlaps`), not ppm units. This
function returns all surviving pairs; callers decide how to pick one
match per ion (e.g. closest RT then m/z in
`get_Spectra_ms2_feature_id`).

## See also

[`match_mz_foverlaps`](https://drruili.github.io/MSdev/reference/match_mz_foverlaps.md),
[`match_mz`](https://drruili.github.io/MSdev/reference/match_mz.md),
`get_Spectra_ms2_feature_id`
