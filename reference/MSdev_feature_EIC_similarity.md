# EIC similarity per polarity without grouping

Per-polarity wrapper around
[`xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/xcms_feature_EIC_similarity.md).
Shares chromatogram handling with
[`MSdev_group_feature_EIC`](https://drruili.github.io/MSdev/reference/MSdev_group_feature_EIC.md)
but does not assign `feature_group` labels.

## Usage

``` r
MSdev_feature_EIC_similarity(
  object,
  rt_tol = 5,
  expandRt = 2,
  min_width = 20,
  expandMzppm = 2,
  selected_sample = NULL,
  forceExtractChrom = FALSE,
  prob = 0.75,
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

- expandRt:

  numeric(1). Seconds added on each side of each feature's
  `peakRtMin`/`peakRtMax` window when cropping EICs via
  [`XChromatograms_subset_feature`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
  before correlation. Default `2`.

- min_width:

  numeric(1). Minimum RT window width (seconds) after `expandRt`;
  shorter windows are padded equally on both sides. Default `20`.

- expandMzppm:

  numeric(1). Extra m/z pad in ppm used when extracting chromatograms
  via
  [`MSdev_get_feature_chrom`](https://drruili.github.io/MSdev/reference/MSdev_get_feature_chrom.md)
  (`mzmin * (1 - ppm/1e6)`, `mzmax * (1 + ppm/1e6)`). Default `2`. Has
  no effect on already-stored chromatograms unless
  `forceExtractChrom = TRUE`.

- selected_sample:

  NULL, integer index/indices, or sample name(s) (`sample.name` /
  chromatogram colnames). NULL uses all samples.

- forceExtractChrom:

  logical(1). If TRUE, (re)extract chromatograms via
  `MSdev_get_feature_chrom` even if already stored.

- prob:

  numeric(1). Quantile used to aggregate samples. Default 0.75.

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

MSdev object; each polarity's MS1 xcms object carries
`otherData(xcms)$EIC_Similarity` and
`otherData(xcms)$EIC_Similarity_agg`.
