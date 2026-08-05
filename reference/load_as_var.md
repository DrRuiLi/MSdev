# Load RData File as Variable

Loads an RData file and returns the single variable it contains. The
function validates that exactly one variable exists in the file before
loading.

## Usage

``` r
load_as_var(file_to_load)
```

## Arguments

- file_to_load:

  Character string specifying the path to the RData file to load.

## Value

The data object stored in the RData file. The type depends on what was
saved in the file.
