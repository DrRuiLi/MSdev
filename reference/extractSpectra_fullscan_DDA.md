# extractSpectra_fullscan_DDA

extract all MS2 Spectra from `object@sampleInfo$msDataFile` which
`sampleInfo$xcmsProcessing` %in% % c("Both","MS2"), return store in
`object@spectra$positiveMS2`

## Usage

``` r
extractSpectra_fullscan_DDA(object)
```

## Arguments

- object:

  a `MSdev` object

## Value

a `MSdev` object
