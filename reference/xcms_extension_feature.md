# XCMS feature helpers

Xcms feature se.

Extracts and adds median retention time, signal-to-noise ratio, and
maximum intensity for each feature. While xcms::featureDefinitions()
provides median mz and rt, this function calculates median values across
all peaks within a feature: peakRtMin, peakRtMax, peakWidth, peakMzMin,
peakMzMax, peakSN, peakMaxo, and polarity.

Convert xcms features to SummarizedExperiment and compute feature-level
statistics.

## Usage

``` r
get_xcms_feature_se(xcms.xcms, ...)

xcms_get_feature_def_stat(xcms.xcms)
```

## Arguments

- xcms.xcms:

  XCMSnExp object with feature definitions and chromPeaks.

- missing:

  how missing values should be reported. Allowed values are NA (the
  default), a numeric or missing = "rowmin_half". The latter replaces
  any NA with half of the row's minimal (non-missing) value.

## Value

SummarizedExperiment

XCMSnExp object with updated featureDefinitions containing additional
statistics.

## Functions

- `get_xcms_feature_se()`: extract feature data from xcms, convert to
  SummarizedExperiment

- `xcms_get_feature_def_stat()`: calculate feature statistics
