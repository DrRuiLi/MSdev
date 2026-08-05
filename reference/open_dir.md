# Open Directory in File Explorer

Opens a directory in the system's file explorer (Windows Explorer on
Windows). If a file path is provided, opens the containing directory.

## Usage

``` r
open_dir(x = getwd())
```

## Arguments

- x:

  Character string specifying the path to a directory or file. Defaults
  to the current working directory.

## Value

The normalized path that was opened, returned invisibly.
