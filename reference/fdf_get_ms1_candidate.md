# MS1 candidate matching on a featureDefinitions-style data.frame

MS1 candidate matching on a featureDefinitions-style data.frame

## Usage

``` r
fdf_get_ms1_candidate(
  fdf,
  cpdb,
  polarity,
  mz.ppm = 10,
  mz_range = NULL,
  selected_adduct = MSCC::adduct.table$Adduct,
  ...
)
```

## Arguments

- fdf:

  data.frame with at least `mzmed` (and preferably `feature_id`)

- cpdb:

  CompoundDb object

- polarity:

  integer polarity/polarities to keep for adducts

- mz.ppm:

  m/z tolerance in ppm

- mz_range:

  numeric length-2 m/z range for filtering DB adducts; default
  `range(fdf$mzmed)`

- selected_adduct:

  adduct names

- ...:

  unused

## Value

`fdf` with candidate.\* list columns
