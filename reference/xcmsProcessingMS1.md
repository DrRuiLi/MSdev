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
  BPPARAM = BiocParallel::SnowParam(workers = max(1L, floor(parallel::detectCores()/3)),
    progressbar = T),
  ...
)
```

## Arguments

- ion_mode:

  to filter ion_mode, 1: positive, 0: negative, import when scans with
  both pos and neg

- beta_cor_thresh:

  optional numeric; if set, fill NA `beta_cor` after merge with
  `chromPeakSummary()` (finite CentWave scores are kept), then drop
  chromPeaks below this value (remaining NA scores are dropped). Default
  `NULL` skips fill and filtering.

- msDataFiles:

  `char` ms file (full) paths

- peaksGroup:

  `vector` to xcms::PeakGroupsParam(sampleGroups), should contain "QC"

- centWaveParam:

  xcms::CentWaveParam()

## Value

xcms
