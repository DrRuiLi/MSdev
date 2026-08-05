# Msdev Export

Export `advancedAna` SummarizedExperiment objects to Excel under
`projectDir/Statistic`. Skips metabolite and candidate exports when
those slots are absent (e.g. after
[`MSdev_get_Se`](https://drruili.github.io/MSdev/reference/MSdev_workflow.md)
only).

## Usage

``` r
MSdev_export(object, candi = FALSE)
```

## Arguments

- object:

  MSdev

- candi:

  logical; export `candidate.se` when present

## Details

MSdev_export
