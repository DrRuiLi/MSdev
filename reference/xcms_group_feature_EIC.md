# Group xcms features by EIC similarity within RT tolerance

Compare extracted-ion chromatogram shapes for feature pairs with
`|rtmed_i - rtmed_j| < rt_tol` via
[`get_xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/get_xcms_feature_EIC_similarity.md),
aggregate across selected samples (75\\ labels with `method`, and
optionally store per-sample sparse similarity matrices in
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
  method = c("complete_linkage", "hclust_average"),
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- xcms.xcms:

  XcmsExperiment / MsExperiment with feature definitions and an
  `otherData` slot. (all methods)

- chroms:

  Chromatograms with feature rownames matching `featureDefinitions`.
  (all methods)

- rt_tol:

  numeric(1). Maximum absolute RT difference (seconds) for which EIC
  similarity is computed. Default 5. (all methods)

- threshold:

  numeric(1). Similarity cut-off for the chosen `method` (see Details).
  Default 0.5. (all methods)

- expandRt:

  numeric(1). Seconds added on each side of each feature's
  `peakRtMin`/`peakRtMax` window when cropping EICs via
  [`XChromatograms_subset_feature`](https://drruili.github.io/MSdev/reference/xcms_extension_chromatogram.md)
  before correlation. Default `2`. (all methods)

- min_width:

  numeric(1). Minimum RT window width (seconds) after `expandRt`;
  shorter windows are padded equally on both sides. Default `20`. (all
  methods)

- selected_sample:

  NULL, integer index/indices, or sample name(s). NULL uses all samples.
  (all methods)

- keep_Similarity_Matrix:

  logical(1). If TRUE (default), store the named list of per-sample
  `dgCMatrix` similarity matrices in `otherData(xcms)$EIC_Similarity`.
  If FALSE, matrices are discarded after grouping. (all methods)

- absent_sim:

  numeric(1). Fill value for pairs outside `rt_tol` when densifying the
  sparse similarity matrix before grouping. Default `0` (non-overlap
  treated as dissimilar). Interpretation depends on `method` (see
  Details). (all methods)

- method:

  character(1). Grouping algorithm on the dense similarity matrix. One
  of `"complete_linkage"`, `"hclust_average"`. Default
  `"complete_linkage"`. See Details for arguments used by each method.

- BPPARAM:

  BiocParallel backend passed to
  [`get_xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/get_xcms_feature_EIC_similarity.md).
  Default `SerialParam()`. (all methods)

## Value

Updated `xcms.xcms` with `featureGroups` set and a `feature_group_rt`
column added to `featureDefinitions` (the median `rtmed` of each group's
member features); when `keep_Similarity_Matrix` is TRUE, also
`otherData(xcms)$EIC_Similarity`.

## Details

**Shared arguments** (used before grouping for every `method`):
`xcms.xcms`, `chroms`, `rt_tol`, `expandRt`, `min_width`,
`selected_sample`, `keep_Similarity_Matrix`, `absent_sim`, `BPPARAM`.

**Arguments by `method`** (how grouping uses densified similarity):

- `complete_linkage`:

  Calls
  [`groupSimilarityMatrix_completeLinkage`](https://drruili.github.io/MSdev/reference/groupSimilarityMatrix_completeLinkage.md)`(x, threshold)`.
  Available args: `threshold` (minimum finite similarity to every
  current group member to join); `absent_sim` (`NA` leaves pairs unknown
  and blocks joins that require those pairs; numeric fill is used as the
  similarity value). Helper-only args not exposed here: `full`, `...`.

- `hclust_average`:

  Calls
  [`groupSimilarityMatrix_hclustAverage`](https://drruili.github.io/MSdev/reference/groupSimilarityMatrix_hclustAverage.md)`(x, threshold)`.
  Available args: `threshold` (tree cut at height `1 - threshold`, i.e.
  average similarity \\\ge\\ threshold); `absent_sim` (non-finite
  entries are set to similarity 0 before converting to distance
  `1 - sim`).
