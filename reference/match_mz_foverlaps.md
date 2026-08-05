# Match m/z values using interval overlaps

Match two numeric m/z vectors within a ppm tolerance using
[`data.table::foverlaps()`](https://rdrr.io/pkg/data.table/man/foverlaps.html)
on ppm-expanded intervals.

The function returns index pairs (`ion1`, `ion2`) and the ppm error
(`mz.ppm`). When `length(mz1) > length(mz2)`, it builds intervals around
`mz2` using a global window derived from `max(ppm.base)`, then filters
to the final `ppm`.

## Usage

``` r
match_mz_foverlaps(mz1, mz2, ppm.base = mz1, ppm = 10)
```

## Arguments

- mz1:

  Numeric vector of m/z values (query set 1).

- mz2:

  Numeric vector of m/z values (query set 2).

- ppm.base:

  Numeric vector used as ppm denominator. Defaults to `mz1`. For typical
  use, `ppm.base` should have the same length as the set referenced by
  `ion1` in the returned table (usually `mz1`).

- ppm:

  Numeric scalar ppm tolerance (default 10).

## Value

A `data.table` with columns:

- `ion1`: index in `mz1`

- `ion2`: index in `mz2`

- `mz.ppm`: ppm error (signed when `length(mz1) > length(mz2)`,
  otherwise absolute)

## Examples

``` r
mz1 <- c(100.0000, 200.0000)
mz2 <- c(100.0008, 199.9990)
match_mz_foverlaps(mz1, mz2, ppm = 10)
#>     ion1  ion2 mz.ppm
#>    <int> <int>  <num>
#> 1:     1     1      8
#> 2:     2     2      5
```
