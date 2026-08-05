# XCMS isotopologue helpers

Screens isotopologue peaks based on m/z and retention time differences,
assigns isotopologue groups and seeds, and records results in
featureDefinitions. Uses graph-based clustering to identify isotopologue
networks.

TODO: unfinished. Similar to `xcms_get_feature_isotopologues` but
supports multiple isotope labels simultaneously (e.g., `[13]C` and
`[15]N`).

Calculates isotopologue-to-seed ratios and determines traced
isotopologues (label-enriched isotopologues) using one of two methods:

- `untraced_compare` (legacy): compare traced and untraced sample
  sources.

- `natural_based`: compare observed ratio to theoretical natural isotope
  ratio derived from
  [`MSCC::chemform_isotopes_pattern_enviPat()`](https://rdrr.io/pkg/MSCC/man/chemform_isotopes_pattern_enviPat.html).

Results are written to featureDefinitions as `is_labeled` and
`Ratio_to_seed_*` columns (same output contract as the legacy function).

Calculates the fraction of isotopologue intensities relative to their
seed feature intensities for each sample. Returns a matrix of fractions
without natural abundance adjustment.

Identify isotopologues, compute labeling ratios, and isotopologue
fractions.

## Usage

``` r
xcms_get_feature_isotopologues(
  xcms.xcms,
  iso_ele = "[13]C",
  max_label = 10,
  ppm = 10,
  rt.tol = 5,
  net.degree.ratio = 0.3
)

xcms_get_feature_isotopologues_multi_tracer(
  xcms.xcms,
  iso_ele = c("[13]C", "[15]N"),
  max_label = c(`[13]C` = 30, `[15]N` = 10),
  ppm = 5,
  rt.tol = 5,
  net.degree.ratio = 0.3
)

xcms_get_feature_traced_isotopologue(
  xcms.xcms,
  iso_ele = "[13]C",
  method = c("untraced_compare", "natural_based")[1],
  ...
)

get_xcms_iso_fraction(xcms.xcms)
```

## Arguments

- xcms.xcms:

  XCMSnExp object with isotopologue assignments (iso_seed column).

- iso_ele:

  Isotope element string (e.g., `"[13]C"`) used for labeling.

- max_label:

  Named numeric vector of maximum labels per tracer, names must match
  `iso_ele`.

- ppm:

  Mass accuracy tolerance in ppm (default 5).

- rt.tol:

  Retention time tolerance in seconds (default 5).

- net.degree.ratio:

  Ratio threshold for network degree to assign isotopologue seeds
  (default 0.3).

- method:

  Labeling method: `"untraced_compare"` or `"natural_based"`. (Legacy
  aliases `"method1"` / `"method2"` are also accepted.)

- ...:

  Additional arguments passed to internal functions.

## Value

XCMSnExp object with featureDefinitions updated with iso_seed,
iso_count, and iso_connection_group columns.

XCMSnExp object with featureDefinitions updated with iso_seed,
iso_count, iso_connection_group, and per-tracer iso_count\_\* columns.

XCMSnExp object with featureDefinitions updated with is_labeled column
and Ratio_to_seed\_\* columns.

Matrix with rows as features and columns as samples, containing
intensity ratios to seed features.

## Functions

- `xcms_get_feature_isotopologues()`: identify isotopologues

- `xcms_get_feature_isotopologues_multi_tracer()`: identify
  isotopologues with multiple isotope tracers

- `xcms_get_feature_traced_isotopologue()`: calculate
  traced-isotopologue labeling ratios

- `get_xcms_iso_fraction()`: calculate isotopologue fractions
