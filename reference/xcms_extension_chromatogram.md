# XCMS chromatogram helpers

Fast per-file EIC extractor. Loads MS1 peaks once per sample, fills an
intensity matrix for all mz-rt boxes, then wraps
[`MSnbase::Chromatogram`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
objects. Drop-in style replacement for
[`xcms::chromatogram()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
for rectangular mz/rt region extraction.

Extract EICs for chromatographic peaks (xcms `chromPeakChromatograms`
analogue). Uses `get_xcms_chromatogram` as the engine.

Extract EICs for features (xcms `featureChromatograms` analogue) via
`get_xcms_chromatogram`. One shared mz-rt box per feature is applied to
each selected sample.

Changes the retention time units of XChromatograms objects. In some
situations (e.g., SRM data from Thermo), retention times are recorded in
minutes, which can cause errors during peak detection. This function
converts between seconds and minutes.

Crop each feature row of an `XChromatograms` object to its
`peakRtMin`/`peakRtMax` window read from the `featureDefinitions` slot.
Optional `expandRt` (seconds each side) widens the window; if the
resulting width is still below `min_width`, both sides are padded
equally until the width reaches `min_width`. Requires `peakRtMin` and
`peakRtMax` in `obj@featureDefinitions` (e.g. after
`get_xcms_feature_chromatogram` on an object that has run
[`xcms_get_feature_def_stat`](https://drruili.github.io/MSdev/reference/xcms_extension_feature.md)).
`NA` intensities inside the window are set to `0` so every feature
carries a full baseline over its window; because all features share the
sample's master RT grid, this keeps RTs aligned across features for
pairwise comparison.

When using xcms::findChromPeaks, chromatograms with fewer than two data
points cause errors. This function identifies such chromatograms and
adds a duplicate point (time +1, intensity 0) to ensure at least two
points exist.

Plots XChromatograms data as line plots, with options to normalize
intensities to 0-1 range, offset chromatograms for clarity, and
customize colors. Returns a ggplot object.

Extract, adjust, and plot chromatograms from xcms objects.

## Usage

``` r
get_xcms_chromatogram(
  object,
  mz,
  rt,
  aggregationFun = "max",
  BPPARAM = SerialParam(),
  msLevel = 1L,
  ...
)

get_xcms_peaks_chromatogram(
  xcms.xcms,
  peaks.id,
  selected_sample = NULL,
  rt.range = c("expand", "identity", "all"),
  expandRt = 15,
  aggregationFun = "max",
  BPPARAM = SerialParam()
)

get_xcms_feature_chromatogram(
  xcms.xcms,
  feature.id = NULL,
  selected_sample = "maxo",
  rt = c("expand", "identity", "all"),
  expandRt = 15,
  mz.expand = 0,
  aggregationFun = "max",
  attachPeaks = TRUE,
  BPPARAM = SerialParam(progressbar = TRUE)
)

XChromatograms_rt_unit(
  xchroms,
  unit_to = "s",
  BPPARAM = BatchtoolsParam(progressbar = T, log = F, registryargs =
    batchtoolsRegistryargs(packages = c("MSnbase")))
)

XChromatograms_subset_feature(xchroms, expandRt = 0, min_width = 0)

XChromatograms_fill_2point(xchroms)

plot_XChromatograms(
  xchroms,
  norm = T,
  move = T,
  color_by = c("column", "row"),
  color_f = NULL,
  label_df = NULL
)
```

## Arguments

- object:

  XCMSnExp / XcmsExperiment (or single-file subset).

- mz:

  numeric matrix with columns mzmin, mzmax (one row per EIC).

- rt:

  one of `c("all","expand","identity")`.

- aggregationFun:

  passed to `get_xcms_chromatogram`.

- BPPARAM:

  BiocParallel backend for parallel processing. Default is
  BatchtoolsParam.

- msLevel:

  integer; kept for API compatibility (MS1 extraction).

- ...:

  ignored (compatibility with older callers).

- xcms.xcms:

  XCMSnExp / XcmsExperiment with featureDefinitions.

- peaks.id:

  character or numeric peak IDs / indices.

- selected_sample:

  Sample selection. `"maxo"` (default) uses the sample with highest mean
  feature value; `"all"` uses all samples; integer indices or sample
  name(s) select those samples.

- rt.range:

  one of `c("all","identity","expand")`.

- expandRt:

  numeric(1). Seconds added on each side of `[peakRtMin, peakRtMax]`.
  Default `0` (no expansion).

- feature.id:

  character/numeric feature IDs (default all).

- mz.expand:

  fraction of mz width to expand on each side.

- attachPeaks:

  logical; attach feature chromPeaks into `XChromatograms` (needed for
  `removeIntensity(..., "outside_chromPeak")`).

- xchroms:

  XChromatograms object to plot.

- unit_to:

  Target unit: "s" (seconds) multiplies by 60, "m" (minutes) divides
  by 60. Default is "s".

- min_width:

  numeric(1). Minimum RT window width (seconds) after `expandRt`. If the
  window is shorter, pad both sides equally. Default `0` (no minimum).

- norm:

  logical. If TRUE, normalize intensities to 0-1 range (default TRUE).

- move:

  logical. If TRUE, offset chromatograms by index for better visibility
  (default TRUE).

- color_by:

  Character indicating grouping for coloring: "column" (by sample) or
  "row" (by feature). Default is "column".

- color_f:

  Optional character vector of colors for groups. If NULL, uses
  distinctColorPalette.

- label_df:

  Optional data frame with columns x, y, label for adding text labels
  via ggrepel.

## Value

`MChromatograms` with rows = regions, columns = samples.

MChromatograms / column-bound chromatograms.

XChromatograms. The parent `featureDefinitions` rows for the extracted
features are stored in the `featureDefinitions` slot (access with
`obj@featureDefinitions`), so the object is self-describing (carries
`feature_id`, `mzmed`, `rtmed`, and `peakRtMin`/`peakRtMax` when
present).

XChromatograms object with converted retention times.

XChromatograms with each feature's EIC cropped to its RT window and `NA`
intensities filled with `0`. The `phenoData`, `featureData`, and
`featureDefinitions` slots are carried over from the input unchanged.

XChromatograms object with chromatograms having at least two data
points.

ggplot object.

## Functions

- `get_xcms_chromatogram()`: extract chromatograms

- `get_xcms_peaks_chromatogram()`: extract chromatogram for a peak

- `get_xcms_feature_chromatogram()`: extract chromatograms for features

- `XChromatograms_rt_unit()`: convert retention time units

- `XChromatograms_subset_feature()`: subset feature chromatograms by
  peak RT

- `XChromatograms_fill_2point()`: fill chromatograms with fewer than two
  data points

- `plot_XChromatograms()`: plot chromatograms
