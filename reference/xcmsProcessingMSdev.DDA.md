# Process DDA data using xcms

Perform peak detection on DDA data using xcms, grouping features across
samples.

## Usage

``` r
xcmsProcessingMSdev.DDA(
  object,
  BPPARAM = BiocParallel::SnowParam(workers = 4, progressbar = TRUE),
  ...
)
```

## Arguments

- object:

  MSdev

- BPPARAM:

  BiocParallel backend passed to xcms peak processing.

- ...:

  additional arguments passed to xcms functions

## Value

MSdev
