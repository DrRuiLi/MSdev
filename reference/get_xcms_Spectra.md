# get_xcms_Spectra

**Deprecated.** Use `ProtGenerics::spectra(xcms.xcms)` (and optionally
[`Spectra::filterMsLevel()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
/
[`Spectra::filterPolarity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html))
instead.

Previously: build a
[Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html) object
from the raw files behind an `XCMSnExp`, filtered to the polarity of
`xcms.xcms`. Scan IDs from `get_xcms_scan_Stat()` are assigned as
`spectraNames` and the `scan_id` spectra variable.

## Usage

``` r
get_xcms_Spectra(xcms.xcms)
```

## Arguments

- xcms.xcms:

  An `XCMSnExp` object.

## Value

A `Spectra` object with `scan_id` matching xcms scan metadata.
