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

Unnamed numeric scores aligned to `keys`, or `NULL`. Each value is the
highest library-spectrum score for that `compound_id` + adduct.
Experimental spectra are still reduced with `max` per library spectrum
before that.
