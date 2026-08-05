# Normalize Values Using Min-Max Scaling

Applies min-max normalization to scale values to the range `[0, 1]`.
Works on vectors, matrices, and data frames (normalizes each row for 2D
structures).

## Usage

``` r
normalize_max_min(x)
```

## Arguments

- x:

  A numeric vector, matrix, or data frame to be normalized.

## Value

The normalized data with the same structure as the input. For vectors,
returns a numeric vector. For matrices/data frames, returns a matrix
with values normalized row-wise.
