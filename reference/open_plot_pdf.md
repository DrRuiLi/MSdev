# open_plot_pdf

Write a plot to a temporary PDF and open it (Windows default viewer).
Supports `ggplot`, ComplexHeatmap (`Heatmap` / `HeatmapList`), and other
objects handled by
[`export::graph2pdf()`](https://rdrr.io/pkg/export/man/graph2vector.html).

## Usage

``` r
open_plot_pdf(p, width = 5, height = 4)
```

## Arguments

- p:

  A `ggplot`, ComplexHeatmap, or other plot object.

- width:

  PDF width in inches.

- height:

  PDF height in inches.

## Value

`NULL`, invisibly. Side effect: opens the temp PDF via `open_file()`.

## See also

[`open_plot_win`](https://drruili.github.io/MSdev/reference/open_plot_win.md),
[`open_plot_ppt`](https://drruili.github.io/MSdev/reference/open_plot_ppt.md),
[`export_graph2pdf`](https://drruili.github.io/MSdev/reference/export_graph2pdf.md)
