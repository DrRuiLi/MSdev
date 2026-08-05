# analyzePathwayGlobalTest

pathway enrichment by global test

AnalyzePathwayHyperTest.

## Usage

``` r
analyzePathwayGlobalTest(pathway.matrix, pathway.group, filter_Metabolism = F)

analyzePathwayHyperTest(kegg.id = "C00024", filter_Metabolism = F)
```

## Arguments

- pathway.matrix:

  should be a matrix, sample as rowname, kegg id as colname

- pathway.group:

  should be a vector with length same as nrow(pathway.matrix), indicate
  group of sample

- filter_Metabolism:

  only output pathway of Metabolism

- kegg.id:

  kegg.id

## Value

global.test.result

DF

## Functions

- `analyzePathwayGlobalTest()`: analyzePathwayGlobalTest

- `analyzePathwayHyperTest()`: analyzePathwayHypertest
