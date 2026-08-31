# xcmsProcessingMS1

Import `msDataFiles`, filter `ion_mode`, find peaks using
`centWaveParam`, correct RT, group peaks using `peaksGroup`, fill peaks
by xcms at MS1 Level

## Usage

``` r
xcmsProcessingMS1(
  xcms.xcms,
  ion_mode = NA,
  xcms_param = list(findChromPeaks = xcms::CentWaveParam(), groupChromPeaks =
    xcms::PeakDensityParam(sampleGroups = "A")),
  adjustRT = T,
  chromPeaks_fix_mz_ppm = NULL,
  chromPeaks_max_mz_ppm = NULL,
  beta_cor_thresh = NULL,
  BPPARAM = BiocParallel::SnowParam(workers = 4, progressbar = T),
  ...
)
```

## Arguments

- ion_mode:

  to filter ion_mode, 1: positive, 0: negative, import when scans with
  both pos and neg

- beta_cor_thresh:

  optional numeric; if set, drop chromPeaks with `beta_cor` below this
  value after peak picking (NA scores are kept). Requires CentWave
  `verboseBetaColumns`. Default `NULL` skips filtering.

- msDataFiles:

  `char` ms file (full) paths

- peaksGroup:

  `vector` to xcms::PeakGroupsParam(sampleGroups), should contain "QC"

- centWaveParam:

  xcms::CentWaveParam()

## Value

xcms
