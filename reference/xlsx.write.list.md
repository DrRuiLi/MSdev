# Write list of data frames to Excel

Writes a list of data frames (or other objects) to a single Excel
workbook, each element as a separate sheet. Sheet names are taken from
list names, with invalid characters replaced by underscores.

## Usage

``` r
xlsx.write.list(df.list, file)
```

## Arguments

- df.list:

  List of data frames or objects to write.

- file:

  File path for the output Excel workbook.

## Value

Invisible NULL.
