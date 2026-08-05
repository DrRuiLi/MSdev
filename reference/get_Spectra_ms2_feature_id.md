# Assign each MS2 spectrum to an xcms feature

Match MS2 spectra to feature definitions by precursor m/z and retention
time using
[`match_mz_rt`](https://drruili.github.io/MSdev/reference/match_mz_rt.md),
then pick one feature per spectrum.

## Usage

``` r
get_Spectra_ms2_feature_id(sp, featuredef, ppm = 5, rt.tol = 5)
```

## Arguments

- sp:

  `Spectra` object (MS2; filtered to `msLevel == 2` internally)

- featuredef:

  Feature definitions data.frame (must contain `mzmed`, `rtmed`,
  `feature_id`)

- ppm:

  m/z tolerance in ppm (default 5)

- rt.tol:

  retention time tolerance in seconds (default 5)

## Value

`sp.data` data.frame (MS2 spectra metadata) with a `feature_id` column

## Details

Candidates are those with feature `mzmed`/`rtmed` within `ppm` and
`rt.tol` of the spectrum `precursorMz`/`rtime`. When several features
remain for one MS2 (`sp_id`), the closest RT to `rtmed` is preferred
(MS2 is typically triggered near the chromatographic apex); ties are
broken by smaller relative m/z error. Unmatched spectra get
`feature_id = NA`.

## See also

[`match_mz_rt`](https://drruili.github.io/MSdev/reference/match_mz_rt.md),
[`MSdev_assign_MS2`](https://drruili.github.io/MSdev/reference/MSdev_workflow.md)
