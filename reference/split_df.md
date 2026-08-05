# Split Data Frame Randomly

Randomly splits a data frame into multiple subsets of approximately
equal size. Useful for creating training/test splits or cross-validation
folds.

## Usage

``` r
split_df(df, n = 2)
```

## Arguments

- df:

  A data frame to be split.

- n:

  Integer specifying the number of subsets to split the data frame into.
  Default is 2.

## Value

A list of data frames, with each element containing one subset.
