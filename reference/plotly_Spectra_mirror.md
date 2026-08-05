# Create Interactive Mirror Plot of Two Spectra

Creates an interactive plotly mirror visualization comparing two mass
spectra. The first spectrum is displayed on top, the second on bottom
(negative intensity). Matched peaks are highlighted with markers.

## Usage

``` r
plotly_Spectra_mirror(sp1, sp2)
```

## Arguments

- sp1:

  A `Spectra` object for the top spectrum.

- sp2:

  A `Spectra` object for the bottom spectrum.

## Value

A `plotly` object displaying a mirror plot of the two spectra.

## Details

plotly_Spectra_mirror
