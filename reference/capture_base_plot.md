# capture_base_plot

Evaluate base graphics in a null device and return a `recordedplot`
object for replay or export (e.g. with
[`open_plot_win`](https://drruili.github.io/MSdev/reference/open_plot_win.md)).

## Usage

``` r
capture_base_plot(expr, envir = parent.frame())
```

## Arguments

- expr:

  Base graphics commands, typically passed as a braced expression.

- envir:

  Environment in which `expr` is evaluated.

## Value

A `recordedplot` object (from
[`recordPlot`](https://rdrr.io/r/grDevices/recordplot.html)).
