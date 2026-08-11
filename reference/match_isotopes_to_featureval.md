# Append feature intensity matrix to isotope matches

Append feature intensity matrix to isotope matches

## Usage

``` r
match_isotopes_to_featureval(matched_table, featureval)
```

## Arguments

- matched_table:

  Output from
  [`match_isotopes_to_featuredef()`](https://drruili.github.io/MSdev/reference/match_isotopes_to_featuredef.md).

- featureval:

  Numeric matrix of feature intensities (e.g. from xcms
  `featureValues()`), with feature IDs as row names.

## Value

Data frame with matched rows and joined intensity columns.
