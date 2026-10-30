# EIC similarity for xcms features without grouping

Runs
[`get_xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/get_xcms_feature_EIC_similarity.md),
aggregates the per-sample matrices (`prob` quantile per edge, over
samples where the edge exists), and stores both in `otherData(xcms)`:
`EIC_Similarity` (named list of per-sample `dgCMatrix`) and
`EIC_Similarity_agg` (one symmetric `dgCMatrix`, diagonal 1; pairs
outside `rt_tol` are not stored). `featureGroups` is left untouched, so
callers can build their own (e.g. overlapping) groups with
[`groupSimilarityMatrix_neighborhood`](https://drruili.github.io/MSdev/reference/groupSimilarityMatrix_neighborhood.md).

## Usage

``` r
xcms_feature_EIC_similarity(
  xcms.xcms,
  chroms,
  rt_tol = 5,
  expandRt = 2,
  min_width = 20,
  selected_sample = NULL,
  prob = 0.75,
  n_chunks = NULL,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- xcms.xcms:

  XCMSnExp / XcmsExperiment with feature definitions.

- chroms:

  Chromatograms (e.g. from
  [`get_xcms_feature_chromatogram`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md))
  with feature rownames matching `featureDefinitions`.

- rt_tol:

  numeric(1). Maximum absolute RT difference (seconds) for which EIC
  similarity is computed. Default 5.

- expandRt:

  numeric(1). Seconds added on each side of each feature's
  `peakRtMin`/`peakRtMax` window when cropping EICs via
  [`XChromatograms_subset_feature`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
  before correlation (requires `peakRtMin`/`peakRtMax` in
  `chroms@featureDefinitions`). Default `2`.

- min_width:

  numeric(1). Minimum RT window width (seconds) after `expandRt`;
  shorter windows are padded equally on both sides. Default `20`.

- selected_sample:

  NULL, integer index/indices, or sample name(s) (`sample.name` /
  chromatogram colnames). NULL uses all samples.

- prob:

  numeric(1). Quantile used to aggregate samples. Default 0.75.

- n_chunks:

  NULL or positive integer. Number of BiocParallel chunks into which
  `n_samples * n_rt_pairs` similarity jobs are split. When `NULL`
  (default), uses `min(bpnworkers(BPPARAM), n_jobs)` for backward
  compatibility. Set explicitly (e.g. `200L`) on large batch projects to
  keep Snow payloads small; `BPPARAM` workers then control parallelism
  only.

- BPPARAM:

  BiocParallel backend. Default `SerialParam()`.

## Value

Updated `xcms.xcms`.
