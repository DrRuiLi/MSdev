# Mirror plot of two chromatograms

Overlay two chromatograms on a shared RT axis: `chrom1` above \\y = 0\\
and `chrom2` mirrored below. Intensities are optionally normalized to
each chromatogram's own maximum.

## Usage

``` r
plot_Chromatograph_mirror(
  chrom1,
  chrom2,
  labels = c("chrom1", "chrom2"),
  normalize = TRUE,
  colors = NULL,
  title = NULL
)
```

## Arguments

- chrom1, chrom2:

  A `Chromatogram` / `XChromatogram`, a one-cell `MChromatograms` /
  `XChromatograms`, or a `data.frame` with columns `rt` and `intensity`.

- labels:

  Character length-2 labels for the legend (default
  `c("chrom1", "chrom2")`).

- normalize:

  Logical; if `TRUE` (default), scale each trace to its own max
  intensity.

- colors:

  Optional length-2 color vector for the two traces.

- title:

  Optional plot title.

## Value

A `ggplot` object.
