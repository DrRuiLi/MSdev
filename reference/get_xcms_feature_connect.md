# TODO: unfinished. Build isotope mass shift grid for multi-tracer

Generates all label-count combinations across tracers with their mass
shifts.

Builds all directed feature pairs from `featureDefinitions` whose
retention times differ by less than `rt.tol`. Used as the RT filter
before isotope mass-shift matching.

## Usage

``` r
get_xcms_feature_connect(xcms.xcms, rt.tol = 5)
```

## Arguments

- xcms.xcms:

  `XCMSnExp` with feature definitions.

- rt.tol:

  Retention time tolerance in seconds (default 5).

- iso_ele:

  Character vector of isotope element strings.

- max_label:

  Named numeric vector of maximum labels per tracer.

## Value

Data.frame with columns for each tracer's label count, total.count, and
mass.shift.

Data frame with integer `from`/`to` row indices, feature ids, m/z, rt,
`mz.diff`, `rt.diff`, and `mz.mean` (mean m/z of the pair).
