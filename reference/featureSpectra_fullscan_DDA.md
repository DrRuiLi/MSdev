# featureSpectra_fullscan_DDA

extrat spectra from `MSdev@spectra` according to mz and rt of feature,
extracted spectra store in `object@spectra$positiveFeatureMS2` and
`object@spectra$negativeFeatureMS2`, a list contain `Spectra` object of
each feature, empty `Spectra` with precursorMz and rtime

## Usage

``` r
featureSpectra_fullscan_DDA(object)
```

## Arguments

- object:

  a `MSdev` object

## Value

a `MSdev` object
