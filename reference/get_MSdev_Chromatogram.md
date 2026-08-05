# Retrieve feature chromatograms from MSdev object

Retrieve stored feature chromatograms from
`object@xcmsData$Positive_Chromatograms` / `$Negative_Chromatograms`
(on-disk when applicable), filtered by polarity.

## Usage

``` r
get_MSdev_Chromatogram(object, polarity = c(0, 1))
```

## Arguments

- object:

  MSdev object

- polarity:

  integer; polarity/polarities to retrieve (`0` negative, `1` positive).
  Default `c(0, 1)`. A single polarity returns the chromatogram object;
  multiple polarities return a named list (`Negative` / `Positive`).

## Value

Chromatograms object, or a named list when `length(polarity) > 1`
