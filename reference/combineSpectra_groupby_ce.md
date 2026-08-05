# Combine Spectra by Collision Energy

Combines spectra that share the same collision energy using a
TIC-weighted peak combination method. This is useful for merging
technical replicates or averaged spectra at each collision energy level.

## Usage

``` r
combineSpectra_groupby_ce(sp, minProp = 0.5, ppm = 10, plot = F, ...)
```

## Arguments

- sp:

  A `Spectra` object to combine.

- minProp:

  Numeric between 0 and 1. Minimum proportion of spectra that must
  contain a peak for it to be retained in the combined spectrum. Default
  is `0.5`.

- ppm:

  Numeric. Parts per million tolerance for grouping peaks. Default is
  `10`.

- plot:

  Logical. If `TRUE`, generates plots of the combination results.
  Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`Spectra::combineSpectra()`](https://rdrr.io/pkg/Spectra/man/combineSpectra.html).

## Value

A `Spectra` object with combined spectra for each unique collision
energy.

## Details

combineSpectra_groupby_ce
