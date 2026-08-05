# Create Interactive Isotopic Mirror Plot

Creates an interactive plotly mirror visualization comparing a spectrum
with its isotopically labeled counterpart. Peaks are matched based on
the expected isotopic mass difference and highlighted accordingly.

## Usage

``` r
plotly_Spectra_iso_mirror(
  sp,
  sp.iso,
  ppm = 10,
  iso_mass_diff = 1.003355,
  iso_count = 1
)
```

## Arguments

- sp:

  A `Spectra` object for the natural abundance spectrum (top).

- sp.iso:

  A `Spectra` object for the isotopically labeled spectrum (bottom).

- ppm:

  Numeric. Parts per million tolerance for matching peaks. Default is
  `10`.

- iso_mass_diff:

  Numeric. Expected mass difference between isotopes. Default is
  `1.003355` (carbon-13 difference).

- iso_count:

  Integer. Number of isotope labels to account for. Default is `1`.

## Value

A `plotly` object displaying an isotopic mirror plot with matched peaks
highlighted.

## Details

plotly_Spectra_iso_mirror
