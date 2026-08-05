# Group features by complete-linkage similarity threshold

Reimplementation of
[`MsFeatures::groupSimilarityMatrix`](https://rdrr.io/pkg/MsFeatures/man/groupSimilarityMatrix.html)
with a fixed group-ID lookup when joining an existing group (named keys
instead of integer positional indexing). Join rule is complete-linkage:
a feature may join a group only if similarity to *all* current members
is finite and `>= threshold`.

## Usage

``` r
groupSimilarityMatrix_completeLinkage(x, threshold = 0.9, full = TRUE, ...)
```

## Arguments

- x:

  numeric square similarity matrix (symmetric when `full = TRUE`).

- threshold:

  numeric(1). Minimum similarity to join / stay linked. Default `0.9`.

- full:

  logical(1). If `FALSE`, only the upper triangle is used.

- ...:

  Ignored; kept for call compatibility with MsFeatures.

## Value

Integer vector of group IDs, length `nrow(x)`.
