# Pairwise EIC similarity for xcms features

Compute per-sample feature-by-feature EIC similarity matrices for pairs
with `|rtmed_i - rtmed_j| < rt_tol`. Used by
[`MSdev_group_feature_EIC`](https://drruili.github.io/MSdev/reference/MSdev_group_feature_EIC.md).
Results are stored as symmetric `Matrix` sparse matrices; absent entries
are non-neighbor pairs (filled when densified for grouping; see
`absent_sim` in
[`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md)).
Pairwise scores are pooled across selected samples and split into
`n_chunks` BiocParallel jobs (each chunk ships only the chromatogram
columns it needs). EICs are expected to be zero-filled over their
windows on the shared sample RT grid (via
[`XChromatograms_subset_feature`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)),
so correlation reflects the full peak window rather than only the
non-`NA` overlap.

## Usage

``` r
get_xcms_feature_EIC_similarity(
  xcms.xcms,
  chroms,
  rt_tol = 5,
  expandRt = 2,
  min_width = 20,
  selected_sample = NULL,
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

Named list of feature-by-feature `dgCMatrix` similarity matrices (one
per selected sample).
