# Fix overly wide xcms chromPeaks mz window

For chromPeaks with abnormal m/z width larger than `ppm`, this function
does not remove peaks. Instead, it recalculates and replaces `mzmin` and
`mzmax` around peak center `mz` so final width equals the target ppm
window.

## Usage

``` r
fix_xcms_chromPeaks_mz_width(xcms.xcms, ppm = 20, verbose = TRUE)
```

## Arguments

- xcms.xcms:

  XCMSnExp object.

- ppm:

  numeric ppm threshold/target width, default `20`.

- verbose:

  logical, print summary message.

## Value

XCMSnExp object with fixed `chromPeaks` m/z windows.
