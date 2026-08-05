# Get Spectra Data

Extracts data from a Spectra object into a data.frame with columns for
spectrum ID, mz, intensity, and requested spectrum variables.

## Usage

``` r
get_Spectra_data(sp, var = c("precursorMz", "collisionEnergy"))
```

## Arguments

- sp:

  A `Spectra` object containing mass spectrometry data.

- var:

  Character vector of spectrum variables to include in the output.
  Default is `c("precursorMz", "collisionEnergy")`. Should be valid
  names from
  [`Spectra::spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html).

## Value

A data.frame with columns: `sp.id` (spectrum identifier), `mz`
(mass-to-charge ratio), `intensity` (peak intensity), and any additional
variables specified in `var`.

## Details

get_Spectra_data
