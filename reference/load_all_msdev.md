# Load MSdev-related development packages

Load source trees for `MSdev`, `MSCC`, and `MSIP` via
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html).
Package paths under OneDrive are resolved with
`get_dir_expand_from_onedrive()` so the same call works across machines
that share the same OneDrive layout.

## Usage

``` r
load_all_msdev()
```

## Value

Invisibly returns the character vector of package directories loaded.

## Examples

``` r
if (FALSE) { # \dontrun{
load_all_msdev()
} # }
```
