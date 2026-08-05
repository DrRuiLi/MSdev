# Export MRM chromatograms as PDF

Export MRM chromatograms (XChromatograms or MChromatograms) to PDF file.

## Usage

``` r
export_MChromatograms_Metabolites(
  mchroms,
  file = tempfile(fileext = ".pdf"),
  patched = F
)
```

## Arguments

- mchroms:

  XChromatograms or MChromatograms object

- file:

  output file path, default tempfile(fileext = ".pdf")

- patched:

  if TRUE, combine all plots using patchwork; if FALSE, append each plot

## Value

file path of exported PDF

## Functions

- `export_MChromatograms_Metabolites()`: MRM data analysis
