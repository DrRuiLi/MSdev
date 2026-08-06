# Group xcms features by EIC similarity within RT tolerance

Compare extracted-ion chromatogram shapes for feature pairs with
`|rtmed_i - rtmed_j| < rt_tol` via
[`get_xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/get_xcms_feature_EIC_similarity.md),
aggregate across selected samples (75\\ labels with
[`groupSimilarityMatrix_completeLinkage`](https://drruili.github.io/MSdev/reference/groupSimilarityMatrix_completeLinkage.md),
and optionally store per-sample sparse similarity matrices in
`otherData(xcms)$EIC_Similarity`.

## Usage

``` r
xcms_group_feature_EIC(
  xcms.xcms,
  chroms,
  rt_tol = 5,
  threshold = 0.5,
  expandRt = 2,
  min_width = 20,
  selected_sample = NULL,
  keep_Similarity_Matrix = TRUE,
  absent_sim = 0,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- xcms.xcms:

  XcmsExperiment / MsExperiment with feature definitions and an
  `otherData` slot.

- chroms:

  Chromatograms with feature rownames matching `featureDefinitions`.

- rt_tol:

  numeric(1). Maximum absolute RT difference (seconds) for which EIC
  similarity is computed. Default 5.

- threshold:

  numeric(1). Similarity cut-off for
  [`groupSimilarityMatrix_completeLinkage`](https://drruili.github.io/MSdev/reference/groupSimilarityMatrix_completeLinkage.md).
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

  NULL, integer index/indices, or sample name(s). NULL uses all samples.

- keep_Similarity_Matrix:

  logical(1). If TRUE (default), store the named list of per-sample
  `dgCMatrix` similarity matrices in `otherData(xcms)$EIC_Similarity`.
  If FALSE, matrices are discarded after grouping.

- absent_sim:

  numeric(1). Fill value for pairs outside `rt_tol` when densifying the
  sparse similarity matrix for
  [`groupSimilarityMatrix_completeLinkage`](https://drruili.github.io/MSdev/reference/groupSimilarityMatrix_completeLinkage.md).
  Default `0` (non-overlap treated as dissimilar). Set `NA` / `NA_real_`
  to leave absents unknown.

- BPPARAM:

  BiocParallel backend passed to
  [`get_xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/get_xcms_feature_EIC_similarity.md).
  Default `SerialParam()`.

## Value

Updated `xcms.xcms` with `featureGroups` set and a `feature_group_rt`
column added to `featureDefinitions` (the median `rtmed` of each group's
member features); when `keep_Similarity_Matrix` is TRUE, also
`otherData(xcms)$EIC_Similarity`.
