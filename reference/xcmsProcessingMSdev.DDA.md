# Process DDA data using xcms

Perform peak detection on DDA data using xcms, grouping features across
samples.

## Usage

``` r
xcmsProcessingMSdev.DDA(
  object,
  BPPARAM = BiocParallel::SnowParam(workers = max(1L, floor(parallel::detectCores()/3)),
    progressbar = TRUE),
  ...
)
```

## Arguments

- object:

  MSdev

- BPPARAM:

  BiocParallel backend passed to xcms peak processing.

- ...:

  additional arguments passed to
  [`xcmsProcessingMS1`](https://drruili.github.io/MSdev/reference/xcmsProcessingMS1.md).

## Value

MSdev
