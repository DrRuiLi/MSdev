# Visualizing xcms data

These helpers take an **`XcmsExperiment`** or **`XCMSnExp`** (after peak
picking and, for feature plots, correspondence) and return ggplot2,
patchwork, or ComplexHeatmap objects. Chromatogram overlays take
**`XChromatograms`** / **`XChromatogram`**.

Examples below assume `xcms` is already processed. They are **not
evaluated** when pkgdown builds this page (no dataset in CI).

``` r

library(ggplot2)
library(patchwork)

# xcms: XcmsExperiment or XCMSnExp after findChromPeaks / groupChromPeaks
```

``` mermaid
flowchart TD
  xcms["XcmsExperiment or XCMSnExp"]
  chrom["XChromatograms or XChromatogram"]
  qc["Run QC: TIC, MS1 scan frequency"]
  peaks["Peaks: RT-mz map, scan counts, peak EIC"]
  feats["Features: RT-mz map, feature EIC, injection intensity"]
  xic["XIC two-panel"]
  overlay["plot_XChromatograms or mirror"]
  fg["Feature-group EIC or similarity"]
  xcms --> qc
  xcms --> peaks
  xcms --> feats
  xcms --> xic
  xcms --> chrom
  chrom --> overlay
  xcms --> fg
  chrom --> fg
```

------------------------------------------------------------------------

## 1. Choosing a plotter

| Need | Function | Input class / slots | Returns |
|----|----|----|----|
| MS1 TIC overlay by `pData$group` | [`plot_xcms_TIC()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `XcmsExperiment` / `XCMSnExp` | ggplot |
| MS1 scan frequency vs RT | [`plot_xcms_ms1_scan_freq()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | spectra on the xcms object | ggplot |
| All chromPeaks on the RT–m/z plane | [`plot_xcms_peaks_distribution()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `chromPeaks` | ggplot |
| MS1 scans overlapping each peak | [`plot_xcms_peaks_ms1_scans()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `chromPeaks` + spectra | ggplot |
| MS2 scans overlapping each peak | [`plot_xcms_peaks_ms2_scans()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `chromPeaks` + MS2 spectra | ggplot |
| EIC for one chromPeak | [`plot_xcms_peaks_Chromatogram()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `chromPeaks` | ggplot |
| Features on the RT–m/z plane | [`plot_xcms_features_distribution()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `featureDefinitions` | ggplot |
| EIC for one feature (few samples) | [`plot_xcms_feature_chromatogram()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `featureDefinitions` | ggplot |
| Feature intensity vs injection order | [`plot_xcms_feature_intensity()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | `featureValues` + `pData` | ggplot |
| Two-panel XIC (EIC + m/z–RT points) | [`plot_xcms_xic()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md) | mz/RT-filtered xcms object | patchwork |
| Overlay chromatogram matrix | [`plot_XChromatograms()`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md) | `XChromatograms` / `MChromatograms` | ggplot |
| Tidy chromatogram table | [`get_chroms_data()`](https://drruili.github.io/MSdev/reference/get_chroms_data.md) | `XChromatogram` / `XChromatograms` | data.frame |
| Mirror two traces | [`plot_Chromatograph_mirror()`](https://drruili.github.io/MSdev/reference/plot_Chromatograph_mirror.md) | `XChromatogram` / chromatogram cell | ggplot |
| Pairwise EIC overlay in a feature group | [`plot_xcms_feature_group_EIC_comparasion()`](https://drruili.github.io/MSdev/reference/plot_xcms_feature_group_EIC_comparasion.md) | `featureGroups` + EICs | ggplot |
| EIC similarity heatmap | [`plot_xcms_feature_group_similarity()`](https://drruili.github.io/MSdev/reference/plot_xcms_feature_group_similarity.md) | `otherData$EIC_Similarity` | Heatmap |

------------------------------------------------------------------------

## 2. Run-level QC

### 2.1 Total ion chromatogram

[`plot_xcms_TIC()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
plots MS1 TIC traces, one line per file, colored by
`Biobase::pData(xcms)$group`. Blank and QC get fixed colors (`grey` /
teal); remaining groups use `ggsci::pal_aaas()`. Line alpha shrinks when
many files are present.

``` r

plot_xcms_TIC(xcms, title = "TIC")
```

### 2.2 MS1 scan frequency

[`plot_xcms_ms1_scan_freq()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
bins MS1 scans in successive `rt_window`-wide RT windows (default 5 s)
**per file** and plots `scan_count / rt_window`. Use it to check whether
the MS1 cycle is stable across the gradient (drops often mark
DDA/MS2-heavy regions or source instability).

``` r

plot_xcms_ms1_scan_freq(xcms, rt_window = 5)
```

------------------------------------------------------------------------

## 3. Chromatographic peaks

These plotters need peak picking
([`xcms::chromPeaks`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)).
Subtitles pull CentWave parameters from
[`processHistory()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
(`ppm`, `snthresh`, `prefilter`).

### 3.1 Peak map

[`plot_xcms_peaks_distribution()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
draws every (non-merged, width \< 60 s) peak on RT vs m/z. Color is
`log10(maxo)`; `type = "o"` uses point size for peak width, `type = "l"`
draws a horizontal segment from `rtmin` to `rtmax`.

``` r

plot_xcms_peaks_distribution(xcms, type = "o")
plot_xcms_peaks_distribution(xcms, type = "l")
```

### 3.2 Scans inside each peak

[`plot_xcms_peaks_ms1_scans()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
counts MS1 scans whose RT falls in `[rtmin, rtmax]`. A horizontal line
at 7 scans is a useful CentWave rule of thumb. With many peaks this can
be slow; subset files first.

[`plot_xcms_peaks_ms2_scans()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
counts MS2 scans whose precursor m/z **and** RT fall in the peak box,
and annotates how many peaks have 0–5 MS2 events. The xcms object must
still carry MS2 spectra; an MS1-only `XcmsExperiment` cannot produce
that overlay.

``` r

plot_xcms_peaks_ms1_scans(xcms)
plot_xcms_peaks_ms2_scans(xcms)
```

### 3.3 Peak EIC

[`plot_xcms_peaks_Chromatogram()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
extracts one peak’s EIC via
[`get_xcms_peaks_chromatogram()`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
and fills the `[rtmin, rtmax]` window.

``` r

pids <- rownames(xcms::chromPeaks(xcms))
plot_xcms_peaks_Chromatogram(xcms, peak_id = pids[[1]], rt = "expand")
```

`rt` is passed through as `rt.range`: `"expand"` (default, ±15 s),
`"identity"` (peak box only), or `"all"` (full run). Extraction details:
[Fast chromatogram
extraction](https://drruili.github.io/MSdev/articles/xcms_chromatogram_extraction.md).

------------------------------------------------------------------------

## 4. Features

These plotters need correspondence
([`xcms::featureDefinitions`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)).

### 4.1 Feature map

[`plot_xcms_features_distribution()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
is the feature analogue of the peak map: `rtmed` vs `mzmed`, color =
median `maxo`, size = `peakWidth`.

``` r

plot_xcms_features_distribution(xcms)
```

### 4.2 Feature EIC (inspection)

[`plot_xcms_feature_chromatogram()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
builds a shared mz–RT box from that feature’s `chromPeaks`, extracts
chromatograms, and overlays selected samples. If more than five samples
are present, it keeps **one sample per `pData$group`** (or the first
five if `group` is missing).

``` r

fids <- rownames(xcms::featureDefinitions(xcms))
plot_xcms_feature_chromatogram(xcms, feature.id = fids[[1]])
```

For many features or full-run traces, extract once with
[`get_xcms_feature_chromatogram()`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md),
then plot with
[`plot_XChromatograms()`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
(section 6). The inspection helper still calls
[`xcms::chromatogram()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
internally and is meant for a handful of features.

### 4.3 Intensity along the injection sequence

[`plot_xcms_feature_intensity()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
plots `featureValues` against injection order. Positive polarity sorts
`pData$analysis.time.positive`; negative uses `analysis.time.negative`.
Points are colored by `sample.type` (`Blank` / `QC` / `Sample`).

``` r

plot_xcms_feature_intensity(xcms, feature_id_to_show = fids[[1]])
```

------------------------------------------------------------------------

## 5. Two-panel XIC

[`plot_xcms_xic()`](https://drruili.github.io/MSdev/reference/xcms_extension_plot.md)
is a ggplot2 stand-in for xcms `plot(..., type = "XIC")`:

- **Upper panel:** extracted-ion chromatogram (intensity vs RT)
- **Lower panel:** centroid points (m/z vs RT, colored by intensity)

Supply an **already filtered** `XcmsExperiment` / `XCMSnExp` plus the
same mz/RT windows used for filtering (axis limits).

``` r

mzr <- c(760.58, 760.60)
rtr <- c(400, 460)
xcms.filt <- xcms::filterRt(xcms::filterMz(xcms, mz = mzr), rt = rtr)

plot_xcms_xic(
  xcms.filt,
  mzr = mzr,
  rtr = rtr,
  title = "XIC",
  subtitle = sprintf("mz %.4f–%.4f; rt %.0f–%.0f s", mzr[1], mzr[2], rtr[1], rtr[2])
)
```

Set `return.data = TRUE` to get the patchwork plus the chromatogram and
point `data.frame`s for custom ggplot layers.

------------------------------------------------------------------------

## 6. `XChromatogram` / `XChromatograms`

Extraction and plotting are separate. Convert `XChromatograms` /
`MChromatograms` / `XChromatogram` to a tidy table with
[`get_chroms_data()`](https://drruili.github.io/MSdev/reference/get_chroms_data.md)
(columns `rt`, `intensity`, `row`, `col`), or pass the S4 object to
[`plot_XChromatograms()`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md).

``` r

chroms <- get_xcms_feature_chromatogram(
  xcms, feature.id = fids[1:8], sample = "all", rt = "expand"
)

plot_XChromatograms(chroms, norm = TRUE, move = TRUE, color_by = "column")
plot_XChromatograms(chroms, norm = FALSE, move = FALSE, color_by = "row")
```

| Argument | Effect |
|----|----|
| `norm = TRUE` | [`MSnbase::normalise()`](https://lgatto.github.io/MSnbase/reference/normalise-methods.html) then scale intensities to 0–100 |
| `move = TRUE` | Offset each trace in RT and intensity so overlays do not sit on top of each other |
| `color_by` | `"column"` = sample; `"row"` = feature / region |
| `label_df` | Optional `data.frame(x, y, label)` drawn with `ggrepel` |

[`plot_Chromatograph_mirror()`](https://drruili.github.io/MSdev/reference/plot_Chromatograph_mirror.md)
compares two traces (`XChromatogram`, `Chromatogram`, an
`XChromatograms` cell, or `data.frame(rt, intensity)`), optionally
normalizing each to its own max and flipping the second below zero.

``` r

plot_Chromatograph_mirror(
  chroms[1, 1],
  chroms[2, 1],
  labels = c(fids[[1]], fids[[2]])
)
```

------------------------------------------------------------------------

## 7. Feature-group EIC views

After feature compounding, `xcms::featureGroups()` labels groups.
Pairwise EIC scores (when stored) sit on
`MsExperiment::otherData(xcms)$EIC_Similarity`.

``` r

fg <- unique(stats::na.omit(xcms::featureGroups(xcms)))

plot_xcms_feature_group_EIC_comparasion(
  xcms,
  feature_group = fg[[1]],
  expandRt = 2,
  min_width = 20
)

plot_xcms_feature_group_similarity(xcms, order_by = "feature_group")
```

How groups and the similarity matrix are built: [Feature grouping with
`EicSimilarityParam`](https://drruili.github.io/MSdev/articles/xcms-feature-group-EicSimilarityParam.md).

------------------------------------------------------------------------

## 8. Saving figures

Plotters return ggplot / patchwork / Heatmap objects. Save with ggplot2
or the package export helpers:

``` r

p <- plot_xcms_TIC(xcms)
ggplot2::ggsave("TIC.pdf", p, width = 8, height = 4)
open_plot_pdf(p, width = 8, height = 4)
```

------------------------------------------------------------------------

## 9. Related reading

- [Fast chromatogram
  extraction](https://drruili.github.io/MSdev/articles/xcms_chromatogram_extraction.md)
  — `get_xcms_chromatogram`, `get_xcms_peaks_chromatogram`,
  `get_xcms_feature_chromatogram`
- [Feature grouping with
  `EicSimilarityParam`](https://drruili.github.io/MSdev/articles/xcms-feature-group-EicSimilarityParam.md)
- Code: `R/dev_xcms.R`, `R/dev_plot.R`
