# Filter Peaks Below Noise Level

Filters out mass spectral peaks that are below the noise level. Requires
that a `noise` spectrum variable has been previously calculated (e.g.,
using
[Spectra_get_noise](https://drruili.github.io/MSdev/reference/Spectra_get_noise.md)).
Peaks with intensity less than or equal to the noise level are removed.

## Usage

``` r
Spectra_filter_noise(sp)
```

## Arguments

- sp:

  A `Spectra` object with a `noise` spectrum variable.

## Value

A filtered `Spectra` object with peaks below the noise level removed.

## Details

Spectra_filter_noise
