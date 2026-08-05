# edit_df_in_excel

Edit df in excel.

## Usage

``` r
edit_df_in_excel(df = data.frame(), rowname = T, read_line = T)
```

## Arguments

- df:

  data.frame

- rowname:

  Logical; write row names to the worksheet.

- read_line:

  Logical; if `TRUE` (default), wait for user input and return the
  edited data frame. If `FALSE`, open Excel and return `NULL` without
  reading the file back.

## Value

Edited `data.frame`, or `NULL` when `read_line = FALSE`.
