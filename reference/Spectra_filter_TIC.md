# Filter Spectra by TIC (Top N per group)

Split a `Spectra` object by metadata variables and keep the top-N
spectra by total ion current (TIC) in each split.

## Usage

``` r
Spectra_filter_TIC(
  sp,
  topN = 10,
  split_var = c("sample.source", "CE", "polarity")
)
```

## Arguments

- sp:

  A `Spectra` object.

- topN:

  Integer. Number of spectra to keep per split (default 10). If `NULL`
  or `Inf`, returns `sp` unchanged.

- split_var:

  Character vector of column names (default
  `c("sample.source","CE","polarity")`). Variables not present in
  `Spectra::spectraData(sp)` are ignored. `"CE"` is treated as an alias
  for `"collisionEnergy"`.

## Value

A filtered `Spectra` object.
