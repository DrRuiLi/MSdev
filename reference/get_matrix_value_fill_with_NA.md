# Subset a matrix with NA fill for missing rows/columns

Returns values from `mat` for each row/column name pair in
`rownames_vec` and `colnames_vec`. Missing row or column names are
filled with `NA`.

## Usage

``` r
get_matrix_value_fill_with_NA(
  mat,
  rownames_vec = rownames(mat),
  colnames_vec = colnames(mat),
  drop = T
)
```

## Arguments

- mat:

  Matrix with dimnames.

- rownames_vec:

  Character vector of row names to extract.

- colnames_vec:

  Character vector of column names to extract.

- drop:

  If `TRUE` and the result is 1x1, return a scalar (default `TRUE`).

## Value

Matrix (or vector when `drop` applies) of extracted values.
