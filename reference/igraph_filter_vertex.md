# Subset an igraph to selected vertices

Keeps only vertices in `v` and removes all others.

## Usage

``` r
igraph_filter_vertex(ig, v)
```

## Arguments

- ig:

  An `igraph` object.

- v:

  Vertex indices, names, or a logical vector (passed to
  [`igraph::V`](https://r.igraph.org/reference/V.html)).

## Value

Subsetted `igraph` object.
