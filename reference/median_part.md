# Extract Middle Portion of Vector

Similar to head() and tail(), but returns elements from the middle
portion of a vector centered around the median position. Useful for
examining central tendency of ordered data.

## Usage

``` r
median_part(x, n = 10)
```

## Arguments

- x:

  A numeric or character vector from which to extract elements.

- n:

  Integer specifying the number of elements to return. Default is 10.

## Value

A vector containing up to n elements from the middle portion of the
input vector.
