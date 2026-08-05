# Wrapper for t.test

A convenience wrapper that returns the p-value from t.test, or 1 on
error.

## Usage

``` r
get_p_t_test(...)
```

## Arguments

- ...:

  arguments passed to
  [`stats::t.test`](https://rdrr.io/r/stats/t.test.html)

## Value

numeric p-value (or 1 if the test fails)
