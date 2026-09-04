# Group features by EIC similarity within RT tolerance

Per-polarity wrapper around
[`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md).
Reuses chromatograms from
[`MSdev_get_feature_chrom`](https://drruili.github.io/MSdev/reference/MSdev_get_feature_chrom.md),
updates MS1 `feature_group` labels, and optionally stores per-sample
sparse similarity matrices in `otherData(xcms)$EIC_Similarity`.

## Usage

``` r
MSdev_group_feature_EIC(
  object,
  rt_tol = 5,
  threshold = 0.5,
  expandRt = 2,
  min_width = 20,
  selected_sample = NULL,
  forceExtractChrom = FALSE,
  keep_Similarity_Matrix = TRUE,
  absent_sim = 0,
  method = c("complete_linkage", "hclust_average"),
  n_chunks = NULL,
  BPPARAM = SnowParam(workers = max(1L, floor(parallel::detectCores()/2)), progressbar =
    TRUE)
)
```

## Arguments

- object:

  MSdev object

- rt_tol:

  numeric(1). Maximum absolute RT difference (seconds) for which EIC
  similarity is computed. Default 5.

- threshold:

  numeric(1). Similarity cut-off for the chosen `method` (see
  [`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md)).
  Default 0.5.

- expandRt:

  numeric(1). Seconds added on each side of each feature's
  `peakRtMin`/`peakRtMax` window when cropping EICs via
  [`XChromatograms_subset_feature`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
  before correlation. Default `2`.

- min_width:

  numeric(1). Minimum RT window width (seconds) after `expandRt`;
  shorter windows are padded equally on both sides. Default `20`.

- selected_sample:

  NULL, integer index/indices, or sample name(s) (`sample.name` /
  chromatogram colnames). NULL uses all samples.

- forceExtractChrom:

  logical(1). If TRUE, (re)extract chromatograms via
  `MSdev_get_feature_chrom` even if already stored.

- keep_Similarity_Matrix:

  logical(1). Passed to `xcms_group_feature_EIC`. Default TRUE.

- absent_sim:

  numeric(1). Passed to `xcms_group_feature_EIC`. Default `0`; set `NA`
  for unknown absents.

- method:

  character(1). Grouping method passed to
  [`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md):
  `"complete_linkage"` or `"hclust_average"`. Default
  `"complete_linkage"`.

- n_chunks:

  NULL or positive integer. Passed to
  [`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md)
  /
  [`get_xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/get_xcms_feature_EIC_similarity.md).
  Use on large batch projects to avoid oversized Snow payloads (e.g.
  `200L`).

- BPPARAM:

  BiocParallel backend for chromatogram extraction and EIC similarity
  scoring.

## Value

MSdev object with updated MS1 `feature_group` labels; when
`keep_Similarity_Matrix` is TRUE, also `otherData(xcms)$EIC_Similarity`
on each polarity's MS1 object.
