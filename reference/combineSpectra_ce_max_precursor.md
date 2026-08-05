# Select Spectra with Maximum Precursor Intensity per Collision Energy

For each unique collision energy in the Spectra object, selects the
spectrum with the highest precursor intensity. Useful for selecting the
most representative spectrum when multiple spectra exist at each
collision energy.

## Usage

``` r
combineSpectra_ce_max_precursor(sp)
```

## Arguments

- sp:

  A `Spectra` object containing multiple spectra.

## Value

A `Spectra` object with one spectrum per collision energy (the one with
maximum precursor intensity).

## Details

combineSpectra_ce_max_precursor
