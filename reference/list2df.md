# list2df

convert list to data.frame. In each list, sub-list will be replaced with
"large list" and missing value will be fill by NA. See
[`data.table::rbindlist`](https://rdrr.io/pkg/data.table/man/rbindlist.html)

## Usage

``` r
list2df(x)
```

## Arguments

- x:

  a list to convert to data.frame

## Value

data.frame
