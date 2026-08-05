# Export Graph2pdf

Graph2pdf.

## Usage

``` r
export_graph2pdf(p, file_path, append = F, ...)
```

## Arguments

- p:

  ggplot

- file_path:

  file path

- append:

  logic

- ...:

  additional arguments passed to export functions

## Value

null

## Details

export plot to pdf file pdf() and export::graph2pdf() not support
`append` arg using qpdf::pdf_combine() to realize that function
