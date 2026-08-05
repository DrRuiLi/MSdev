# Create Interactive Plotly Mass Spectrum

Creates an interactive plotly visualization of a single mass spectrum.
Hovering over peaks shows m/z and intensity values. The top N peaks are
highlighted.

## Usage

``` r
plotly_Spectra(
  sp,
  label.top = 10,
  show.info = F,
  transform = c("identity", "log10")
)
```

## Arguments

- sp:

  A `Spectra` object. If multiple spectra are provided, only the first
  is used.

- label.top:

  Integer. Number of highest intensity peaks to highlight. Default is
  `10`.

- show.info:

  Logical. If `TRUE`, displays precursor information (m/z, intensity,
  collision energy) in the plot title. Default is `FALSE`.

- transform:

  Character. Y-axis transform: `"identity"` (default) or `"log10"`.

## Value

A `plotly` object displaying an interactive mass spectrum.

## Details

plotly_Spectra
