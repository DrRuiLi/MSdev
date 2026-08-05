# Connected-component membership for igraph vertices

Returns the component id for each vertex (same as
`igraph::components(ig)$membership`).

## Usage

``` r
get_igraph_membership(ig)
```

## Arguments

- ig:

  An `igraph` object.

## Value

Named integer vector of component memberships (names are vertex names).
