# Final annotation pick on a featureDefinitions-style data.frame

Final annotation pick on a featureDefinitions-style data.frame

## Usage

``` r
fdf_get_feature_annotation(
  fdf,
  cpdb,
  cpdb.keys = c("name", "formula", "smiles"),
  weight_mz = 0.1,
  weight_ms2 = 0.7,
  weight_isopattern = 0.2,
  ...
)
```

## Arguments

- fdf:

  data.frame with candidate.\* and score.ms2

- cpdb:

  CompoundDb object

- cpdb.keys:

  CompDb columns to attach

- weight_mz, weight_ms2, weight_isopattern:

  score weights

- ...:

  unused

## Value

annotated `fdf`
