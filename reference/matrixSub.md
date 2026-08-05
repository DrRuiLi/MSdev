# Matrix Subtraction of Two Vectors

Expands two vectors into matrices and computes element-wise subtraction,
creating a full difference matrix where each element `(i,j)` equals
`v1[i] - v2[j]`.

## Usage

``` r
matrixSub(v1, v2)
```

## Arguments

- v1:

  Numeric vector. Optionally can have names which will be used as row
  names.

- v2:

  Numeric vector. Optionally can have names which will be used as column
  names.

## Value

A matrix of dimensions length(v1) x length(v2) containing all pairwise
differences. Row and column names are preserved if present in the input
vectors.

## Examples

``` r
a <- 3:8
b <- 1:2
matrixSub(a,b)
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    3    2
#> [3,]    4    3
#> [4,]    5    4
#> [5,]    6    5
#> [6,]    7    6
```
