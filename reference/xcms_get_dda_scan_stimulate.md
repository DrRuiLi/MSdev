# Stimulate DDA cycle and assign MS2 to feature

Stimulate DDA cycle and assign ms2 to feature, now just support single
file

## Usage

``` r
xcms_get_dda_scan_stimulate(xcms.xcms, dynamic_time = 60)
```

## Arguments

- xcms.xcms:

  xcms object

- dynamic_time:

  minimum time (seconds) before re-acquiring same feature, default 60

## Value

xcms.xcms object with MS2 assignment
