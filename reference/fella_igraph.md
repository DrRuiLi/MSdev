# Generate FELLA Enrichment igraph

Creates an igraph object from FELLA enrichment results, adding
enrichment scores and visual attributes.

## Usage

``` r
fella_igraph(fella.fella, p = 0.05, node = 1000)
```

## Arguments

- fella.fella:

  A FELLA enrichment object.

- p:

  P-value threshold for including nodes (default 0.05).

- node:

  Maximum number of nodes to include (default 1000).

## Value

An igraph object with vertex attributes for enrichment p-values, colors,
and types.
