# Get XCMS Parameters via Autotuner

Uses the Autotuner algorithm to automatically estimate XCMS peak
detection parameters from an xcms object.

Estimate XCMS peak detection parameters from an xcms object via
Autotuner.

## Usage

``` r
get_xcms_Autotuner(xcms.xcms)
```

## Arguments

- xcms.xcms:

  An xcmsSet object (or similar) containing chromatographic data.

## Value

A list of XCMS parameters as returned by Autotuner's `returnParams`.

## Functions

- `get_xcms_Autotuner()`: get xcms parameters via Autotuner
