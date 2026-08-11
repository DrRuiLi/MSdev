# Match theoretical isotopes to xcms feature definitions

Match theoretical isotopes to xcms feature definitions

## Usage

``` r
match_isotopes_to_featuredef(
  isotopes_table,
  featuredef,
  mz.ppm = 10,
  rt.tol = 10
)
```

## Arguments

- isotopes_table:

  Data frame from MSCC isotope pattern helpers (needs column `m.z`;
  optional `rt`).

- featuredef:

  xcms feature definition data frame (`mzmed`, optional
  `rtmed`/`feature_id`).

- mz.ppm:

  m/z tolerance in ppm.

- rt.tol:

  RT tolerance in seconds if both sides have RT.

## Value

Matched isotope-feature table.
