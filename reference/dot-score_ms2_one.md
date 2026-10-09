# Score one feature MS2 job

Score one feature MS2 job

## Usage

``` r
.score_ms2_one(job)
```

## Arguments

- job:

  List with `exp` (experimental Spectra), `ref` (reference Spectra), and
  `keys` (candidate `compound_id_adduct` labels). Defined in the package
  namespace so Snow workers serialize this function, not the caller
  frame.

## Value

Unnamed numeric scores aligned to `keys`, or `NULL`.
