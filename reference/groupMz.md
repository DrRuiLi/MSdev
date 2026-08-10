# Group m/z Values by ppm Tolerance

Groups numeric m/z values into clusters within a given ppm tolerance
using
[`MsCoreUtils::group()`](https://rdrr.io/pkg/MsCoreUtils/man/group.html).
`NA` inputs are preserved in the output (as `NA` group ids).

## Usage

``` r
groupMz(x, ppm = 10, return.type = c("vector", "data.frame"))
```

## Arguments

- x:

  Numeric vector of m/z values.

- ppm:

  Numeric scalar ppm tolerance passed to
  [`MsCoreUtils::group()`](https://rdrr.io/pkg/MsCoreUtils/man/group.html).
  Default is 10.

- return.type:

  Character; either `"vector"` (default) or `"data.frame"`.

## Value

If `return.type = "vector"`, an integer vector of group ids the same
length as `x` (`NA` where `x` is `NA`). If `return.type = "data.frame"`,
a data frame with columns `mz`, `mz.group`, `mz.center` (median m/z per
group), `mz.diff`, `mz.ppm`, `mz.width`, and `mz.width.ppm`.

## See also

[`MsCoreUtils::group`](https://rdrr.io/pkg/MsCoreUtils/man/group.html)

## Examples

``` r
mz <- c(100.0000, 100.0008, 200.0000, NA, 200.0015)
groupMz(mz, ppm = 10)
#> [1]  1  1  2 NA  2
groupMz(mz, ppm = 10, return.type = "data.frame")
#> # A tibble: 5 × 7
#>      mz mz.group mz.center   mz.diff mz.ppm  mz.width mz.width.ppm
#>   <dbl>    <int>     <dbl>     <dbl>  <dbl>     <dbl>        <dbl>
#> 1  100         1      100.  0.000400   4.00  0.000800         8.00
#> 2  100.        1      100.  0.000400   4.00  0.000800         8.00
#> 3  200         2      200.  0.000750   3.75  0.00150          7.50
#> 4   NA        NA       NA  NA         NA    NA               NA   
#> 5  200.        2      200.  0.000750   3.75  0.00150          7.50
```
