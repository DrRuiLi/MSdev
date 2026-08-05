# Group features across samples using xcms

Group detected features across samples based on retention time and
intensity correlation.

## Usage

``` r
MSdev_xcms_group_features(object, diffRt = 5, intCor = 0.5, eicCor = 0.3, ...)
```

## Arguments

- object:

  MSdev object

- diffRt:

  maximum retention time difference for grouping (seconds)

- intCor:

  minimum intensity correlation threshold

- eicCor:

  minimum EIC correlation threshold

- ...:

  additional arguments passed to xcms grouping functions

## Value

MSdev object with updated feature groups
