# Build MS2 precursor peak table from spectra

Construct an xcms `featureDefinitions`-style peak table from MS2
precursor m/z and retention time. Continuous MS2 spectra
(`|delta rt| < rt_tol` and m/z within `ppm` after m/z grouping) are
collapsed to one peak. Result is stored in
`object@advancedAna$MS2_Precursor`.

## Usage

``` r
MSdev_get_peak_table_from_spectra(object, rt_tol = 5, ppm = 5)
```

## Arguments

- object:

  MSdev object with `MS2_Spectra`

- rt_tol:

  numeric; max RT gap (seconds) between consecutive MS2 in the same peak
  (default `10`)

- ppm:

  numeric; m/z tolerance in ppm for grouping precursors (default `10`)

## Value

MSdev object with `advancedAna$MS2_Precursor` set
