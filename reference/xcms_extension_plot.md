# XCMS diagnostic plots

export peaks data by xcms::chromPeaks and plot by ggplot2

Visualizes the distribution of detected features in a 2D space of
retention time (x-axis) vs m/z (y-axis). Point size represents peak
width, color represents log10 intensity. Includes peak detection
parameters in subtitle.

extract Chromatogram from xcms according to feature's mz range and plot

plot scans number of MS1 levels in each peak, note that to many peaks
will lead to stuck, apply `filterFile` to decrease peaks count

Dot-plot MS1 scan frequency along retention time. Scans are counted in
successive `rt_window`-wide RT bins (per file); frequency is
`scan_count / rt_window`.

Visualizes the number of MS2 scans that overlap each chromatographic
peak based on retention time and m/z ranges. Produces a scatter plot
with jitter, violin distribution, and counts of peaks with 0-5 MS2
scans.

extract EIC according to peaks' mzrange and rtrange, note that if
multiple sample in xcms object, only first sample will be extracted

plot feature's intensity, ordered by
`Biobase::pData(xcms.xcms)$analysis.time.positive` or
`Biobase::pData(xcms.xcms)$analysis.time.negative`

Plot MS1 total ion chromatograms (TIC) for an XCMSnExp object, colored
by sample group from `Biobase::pData(xcms.xcms)$group`.

Diagnostic and overview plots for xcms peaks, features, chromatograms,
and TIC/XIC.

## Usage

``` r
plot_xcms_peaks_distribution(
  xcms.xcms,
  plot.title = "Peaks distribution",
  type = "o"
)

plot_xcms_features_distribution(
  xcms.xcms,
  plot.title = "Features distribution"
)

plot_xcms_feature_chromatogram(xcms.xcms, feature.id, sampleNames = NULL)

plot_xcms_peaks_ms1_scans(xcms.xcms, plot.title = "Peaks Sans of MS1")

plot_xcms_ms1_scan_freq(xcms, rt_window = 5, plot.title = "MS1 Scan Frequency")

plot_xcms_peaks_ms2_scans(xcms.xcms, plot.title = "Peaks Sans of MS2")

plot_xcms_peaks_Chromatogram(xcms.xcms, peak_id, rt = "expand")

plot_xcms_feature_intensity(xcms.xcms, feature_id_to_show)

plot_xcms_TIC(xcms.xcms, col.group = NULL, title = "TIC")

plot_xcms_xic(
  xcms.filt,
  mzr = NULL,
  rtr = NULL,
  title = NULL,
  subtitle = NULL,
  base_size = 6,
  return.data = FALSE
)
```

## Arguments

- xcms.xcms:

  XCMSnExp object

- plot.title:

  Character title for the plot (default "Peaks Sans of MS2").

- type:

  `"o"`, for geom_point, `"l"`, for geom_segment

- feature.id:

  feature id

- sampleNames:

  sample names to include

- xcms:

  XCMSnExp / XcmsExperiment object

- rt_window:

  positive numeric; RT window width (same unit as retention time,
  typically seconds)

- peak_id:

  peak id

- rt:

  expansion range for rt

- feature_id_to_show:

  feature id to plot

- col.group:

  named character vector of colors for groups. If `NULL`, Blank/QC use
  fixed colors and remaining groups use `ggsci::pal_aaas()` (or an
  interpolated palette when there are more than 10 groups)

- title:

  Optional plot title.

- xcms.filt:

  `XCMSnExp` after `filterRt()` and `filterMz()`.

- mzr:

  Optional m/z range used for extraction (for axis limits).

- rtr:

  Optional RT range used for extraction (for axis limits).

- subtitle:

  Optional subtitle.

- base_size:

  Base font size.

- return.data:

  If `TRUE`, return a list with `plot`, `p_chr`, `p_mz`, `chrom`, and
  `points`.

## Value

ggplot object

ggplot object.

ggplot object

ggplot object

ggplot object

ggplot object.

ggplot object

ggplot object

ggplot object

A patchwork object with two panels, or a list when `return.data = TRUE`.

## Functions

- `plot_xcms_peaks_distribution()`: plot peaks distribution

- `plot_xcms_features_distribution()`: plot features distribution

- `plot_xcms_feature_chromatogram()`: plot feature chromatogram

- `plot_xcms_peaks_ms1_scans()`: plot MS1 scan counts for peaks

- `plot_xcms_ms1_scan_freq()`: plot MS1 scan frequency vs RT

- `plot_xcms_peaks_ms2_scans()`: plot MS2 scan counts for peaks

- `plot_xcms_peaks_Chromatogram()`: plot chromatogram for a peak

- `plot_xcms_feature_intensity()`: plot feature intensity

- `plot_xcms_TIC()`: plot TIC

- `plot_xcms_xic()`: ggplot2 XIC plot matching xcms
  `plot(type = \"XIC\")`

  Upper panel: extracted ion chromatogram (intensity vs retention time).
  Lower panel: m/z vs retention time with points coloured by intensity.
