# Update feature mz/rt using peak-intensity weighted means

Recomputes `mzmed` and `rtmed` in `xcms::featureDefinitions(xcms.xcms)`
using peak-level `mz`/`rt` weighted by peak `maxo` (maximum intensity)
from `xcms::chromPeaks(xcms.xcms)`.

## Usage

``` r
xcms_get_feature_wmean(xcms.xcms)
```

## Arguments

- xcms.xcms:

  An
  [`xcms::XCMSnExp`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
  object with feature definitions and chromPeaks.

## Value

An updated
[`xcms::XCMSnExp`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
object where `featureDefinitions(object)$mzmed` and
`featureDefinitions(object)$rtmed` are replaced by the
intensity-weighted means across each feature's constituent peaks.
