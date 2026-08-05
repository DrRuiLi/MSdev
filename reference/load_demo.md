# Load bundled demo data objects

Load pre-built demo objects from the MSdev demo project directory.
Dispatches to
[`MSdev_load()`](https://drruili.github.io/MSdev/reference/MSdev_IO.md)
for `"MSdev"` and [`readRDS()`](https://rdrr.io/r/base/readRDS.html) for
other demo types. `"xcms"` is an alias for `"XcmsExperiment"`
(preferred). Use `"XCMSnExp"` for the legacy on-disk `XCMSnExp` demo.

## Usage

``` r
load_demo(
  demo = c("MSdev", "XcmsExperiment", "xcms", "XCMSnExp", "SummarizedExperiment",
    "data.se", "Spectra", "sp")
)
```

## Arguments

- demo:

  Character. One of `"MSdev"`, `"XcmsExperiment"`, `"xcms"`,
  `"XCMSnExp"`, `"SummarizedExperiment"`, `"data.se"`, `"Spectra"`, or
  `"sp"`.

## Value

The loaded demo object (type depends on `demo`).

## Examples

``` r
if (FALSE) { # \dontrun{
msdev.demo <- load_demo("MSdev")
xcms.demo <- load_demo("xcms")
xe.demo <- load_demo("XcmsExperiment")
} # }
```
