# Heatmap of xcms feature-group EIC similarity

Draws the EIC similarity matrix produced by
[`MSdev_group_feature_EIC`](https://drruili.github.io/MSdev/reference/MSdev_group_feature_EIC.md)
/
[`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md)
(stored on the xcms object at
`MsExperiment::otherData(xcms)$EIC_Similarity`) as a ComplexHeatmap.
Rows/columns are ordered either by `rtmed` (`order_by = "rtmed"`) or by
feature group (`order_by = "feature_group"`: groups sorted by each
group's median `rtmed`, features within a group by `rtmed`). The stored
matrix is sparse: only compared pairs (within `rt_tol`) carry a score,
so absent (uncompared) entries are densified to `NA` and shown in
`na_col` (grey), while a genuine zero correlation stays on the ramp
(white). A top color bar encodes `rtmed` (plus a feature-group bar when
groups exist).

## Usage

``` r
plot_xcms_feature_group_similarity(
  xcms,
  order_by = c("rtmed", "feature_group"),
  sample = 1L,
  rt_window = NULL,
  na_col = "#BDBDBD",
  box_top_n = 10L
)
```

## Arguments

- xcms:

  An xcms object (`XcmsExperiment` / `XCMSnExp`) that carries feature
  definitions, `featureGroups`, and an `EIC_Similarity` entry (e.g.
  `ms@xcmsData$PositiveMS1`).

- order_by:

  `"rtmed"` (default) or `"feature_group"`.

- sample:

  Sample name or index selecting the per-sample similarity matrix
  (default first sample).

- rt_window:

  Optional half-width (seconds) of a fixed RT window centered on the
  densest RT region, keeping the diagonal band visible. `NULL` (default)
  keeps all features.

- na_col:

  Color for uncompared (`NA`) cells (default `"#BDBDBD"`, grey).

- box_top_n:

  Number of largest feature groups to outline with a box on the
  diagonal. Only applies to `order_by = "feature_group"` (where each
  group is a contiguous block). Set `0` / `NULL` to disable (default
  `10`).

## Value

(Invisibly) a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object.
