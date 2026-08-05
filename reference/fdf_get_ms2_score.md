# MS2 spectral similarity scores on a featureDefinitions-style data.frame

MS2 spectral similarity scores on a featureDefinitions-style data.frame

## Usage

``` r
fdf_get_ms2_score(fdf, cpdb, sp.ms2, polarity, ...)
```

## Arguments

- fdf:

  data.frame with `feature_id`, `candidate.*`, `ms2_id`

- cpdb:

  CompoundDb object

- sp.ms2:

  Spectra object (experimental MS2)

- polarity:

  integer polarity for filtering CompDb Spectra

- ...:

  unused

## Value

`fdf` with `score.ms2` list column
