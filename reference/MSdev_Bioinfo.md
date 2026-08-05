# Bioinformatic analysis

PlotPathwayEnrichment.

## Usage

``` r
plotPathwayEnrichment(
  pathway.table,
  top = 20,
  method = c("set1", "set2", "bubble", "class"),
  title = NULL
)
```

## Arguments

- pathway.table:

  table
  from[`analyzePathwayGlobalTest()`](https://drruili.github.io/MSdev/reference/stat-pathway.md)
  or
  [`analyzePathwayHyperTest()`](https://drruili.github.io/MSdev/reference/stat-pathway.md)

- top:

  top n

- method:

  plot style, option: "set1","set2","bubble","path_classify1"

- title:

  title of ggplot

## Value

ggplot

## Functions

- `plotPathwayEnrichment()`: plotPathwayEnrichment

## Examples

``` r
### pathway analysis
path <- analyzePathwayHyperTest()
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Error in pathway.hyper.test$pathway.id[i] <- pathway$ENTRY: replacement has length zero
p <- plotPathwayEnrichment(path)
#> Error: object 'path' not found
```
