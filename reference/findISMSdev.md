# Find internal standard features in MSdev

Find features of internal standards listed in
`object@experimentInfo@Internal_Standard` by `Exact_mass` and
`Retention_time` (if provided). Only `\[M+H\]` and `\[M-H\]` adducts are
considered. Correlation and intensity will be plotted based on
`object@advancedAna[["featureRaw"]]`. A column "internal_standard" will
be added to `object@advancedAna[["featureRaw"]]`

## Usage

``` r
findISMSdev(object, to.adjust = "featureRaw", corr.thred = 0.6)
```

## Arguments

- object:

  MSdev

- corr.thred:

  cor

## Value

MSdev
