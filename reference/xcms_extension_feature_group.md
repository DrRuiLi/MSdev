# Group xcms Features

Groups features from an XCMSnExp object using multiple criteria:
similarity in retention time, abundance (intensity) correlation, and EIC
(extracted ion chromatogram) correlation. RT/abundance grouping uses
MsFeatures; EIC similarity uses xcms::EicSimilarityParam.

Group xcms features by retention time, abundance correlation, and EIC
similarity.

## Usage

``` r
xcms_get_feature_group(xcms.xcms, diffRt = 5, intCor = 0.5, eicCor = 0.5)
```

## Arguments

- xcms.xcms:

  XCMSnExp object containing feature definitions.

- diffRt:

  numeric. Maximum allowed retention time difference for grouping by
  SimilarRtimeParam. If NULL, retention time grouping is skipped.
  Default is 5.

- intCor:

  numeric. Threshold for abundance similarity (correlation) grouping
  using AbundanceSimilarityParam. If NULL, intensity correlation
  grouping is skipped. Default is 0.5.

- eicCor:

  numeric. Threshold for EIC similarity grouping using
  xcms::EicSimilarityParam. If NULL, EIC correlation grouping is
  skipped. Default is 0.5.

## Value

XCMSnExp object with featureGroups column added or updated.

## Functions

- `xcms_get_feature_group()`: group features
