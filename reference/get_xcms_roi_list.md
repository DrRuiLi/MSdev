# Build xcms centWave roiList from mz/rt targets

Construct a `roiList` accepted by `xcms::CentWaveParam(roiList = ...)`.
Input a matrix/data.frame with columns `mz` and `rt` (seconds). For each
target, `mzmin/mzmax` are calculated using ppm tolerance and the RT
window `rtmin/rtmax` is mapped to scan indices.

Scan indices are computed per file from MS1 spectra, then expanded to
the union range across files and **clamped** to `[1, min(n_MS1)]` so one
shared `roiList` stays valid for every file in
[`findChromPeaks()`](https://rdrr.io/pkg/xcms/man/findChromPeaks.html).
Unclamped unions can overflow shorter runs and raise
`Error in scanrange` inside `.centWave_orig`.

Note: `scmin`/`scmax` must be finite integer scan indices (not
`c(0, Inf)`). `centWave` uses them in arithmetic
(`N <- scmax - scmin + 1`) before clipping to each file's
`length(scantime)`.

## Usage

``` r
get_xcms_roi_list(mzrt, xcms.xcms, ppm = 10, rt_tol = 30, ion_mode = NULL)
```

## Arguments

- mzrt:

  matrix/data.frame with columns `mz` and `rt`.

- xcms.xcms:

  `XCMSnExp`, `MsExperiment`, or `XcmsExperiment` used to map RT to scan
  indices.

- ppm:

  numeric, ppm tolerance for mz window.

- rt_tol:

  numeric, RT tolerance in seconds.

- ion_mode:

  optional integer 1 (positive) or 0 (negative). If NULL, inferred from
  scan polarity; must be unique.

## Value

list of ROI objects (each ROI is a list with
`scmin, scmax, mzmin, mzmax, length, intensity`).
