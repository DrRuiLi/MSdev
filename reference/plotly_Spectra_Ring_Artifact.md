# Interactive spectrum with click-locked reference peak

Extends
[`plotly_Spectra()`](https://drruili.github.io/MSdev/reference/plotly_Spectra.md)
with a plotly.js click handler (no Shiny). Click a peak to lock it as
the reference; hover labels then show \\\Delta\\m/z and relative
intensity versus that peak. Click the same peak again to clear the
selection. Useful for inspecting ring artifacts and related
constant-spacing patterns.

## Usage

``` r
plotly_Spectra_Ring_Artifact(
  sp,
  label.top = 10,
  show.info = FALSE,
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

A `plotly`/`htmlwidget` object with client-side click/hover behavior.

## Details

plotly_Spectra_Ring_Artifact
