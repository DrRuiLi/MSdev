# Retrieve spectra from MSdev object

Retrieve stored Spectra from `object@spectra` filtered by MS level and
polarity. MS1 / MS2 are read from `MS1_Spectra` / `MS2_Spectra` (on-disk
when applicable) and combined when both levels are requested.

## Usage

``` r
get_MSdev_Spectra(object, msLevel = c(1, 2), polarity = c(0, 1))
```

## Arguments

- object:

  MSdev object

- msLevel:

  integer; MS level(s) to retrieve. Default `c(1, 2)`.

- polarity:

  integer; polarity/polarities to keep (`0` negative, `1` positive).
  Default `c(0, 1)`.

## Value

Spectra object
