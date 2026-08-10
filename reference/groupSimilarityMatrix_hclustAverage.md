# Group features by average-linkage hierarchical clustering on similarity

Converts a similarity matrix to distance `1 - similarity`, runs
[`stats::hclust`](https://rdrr.io/r/stats/hclust.html) with
`method = "average"` (UPGMA), and cuts the tree at height
`1 - threshold` so clusters keep average similarity at least
`threshold`. Non-finite similarities are set to 0 before clustering.

## Usage

``` r
groupSimilarityMatrix_hclustAverage(x, threshold = 0.9)
```

## Arguments

- x:

  numeric square similarity matrix.

- threshold:

  numeric(1). Minimum average similarity within a cluster (cut height
  `1 - threshold`). Default `0.9`.

## Value

Integer vector of group IDs, length `nrow(x)`.

## Details

Available arguments: `x`, `threshold`.

## See also

[`groupSimilarityMatrix_completeLinkage`](https://drruili.github.io/MSdev/reference/groupSimilarityMatrix_completeLinkage.md),
[`xcms_group_feature_EIC`](https://drruili.github.io/MSdev/reference/xcms_group_feature_EIC.md)
