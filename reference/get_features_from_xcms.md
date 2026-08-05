# get_features_from_xcms

extract feature data from xcms::XCMSnExp, calculate RSD of QC and Sample
( note this rely on character "QC" and "Sample" in
`sampleNames(xcms.xcms)` )

## Usage

``` r
get_features_from_xcms(xcms.xcms, missing = NA)
```

## Arguments

- xcms.xcms:

  XCMSnExp object

## Value

xcms a SummarizedExperiment subject
