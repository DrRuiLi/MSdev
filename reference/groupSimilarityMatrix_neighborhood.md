# Seed-centered overlapping groups from a similarity matrix

For each seed `s`, the group is `s` plus every feature `x` with
`sim(s, x) > threshold` (closed neighborhood). A feature may belong to
several groups. Pairs with no stored value in a sparse matrix (outside
`rt_tol`) and non-finite values are never neighbors.

## Usage

``` r
groupSimilarityMatrix_neighborhood(x, seeds, threshold = 0.9, prefix = "FG.")
```

## Arguments

- x:

  Symmetric similarity matrix (`dgCMatrix` or dense) with feature ids as
  dimnames.

- seeds:

  Character feature ids (matching `colnames(x)`) or integer column
  indices.

- threshold:

  numeric(1). Strict lower bound on similarity to the seed. Default 0.9.

- prefix:

  character(1). Group label prefix; label is `paste0(prefix, seed_id)`.
  Default `"FG."`.

## Value

`data.frame(feature_group, seed_id, feature_id, sim)`, seed row first in
each group with `sim = 1`, members by decreasing `sim`.
