# Plot Mass Spectra

Creates a ggplot2 bar plot visualization of mass spectra data. Peaks
above 10% of maximum intensity are highlighted in blue.

## Usage

``` r
plot_Spectra(sp, label.top = 10)
```

## Arguments

- sp:

  A `Spectra` object to plot.

- label.top:

  Integer. Number of highest intensity peaks to label with m/z values.
  Default is `10`.

## Value

A `ggplot` object displaying the mass spectrum.

## Details

plot_Spectra
