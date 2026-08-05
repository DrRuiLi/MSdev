# XCMS MS1 annotation helpers

Matches features in an XCMSnExp object to compounds in a CompoundDb
database using m/z and retention time tolerance. Calculates adduct
masses for each compound and finds matches within specified ppm error.
Results are stored as candidate lists in featureDefinitions.

Match xcms features to compound databases by MS1 m/z and retention time.

## Usage

``` r
xcms_get_feature_ms1_candidate(
  xcms.xcms,
  cpdb,
  mz.ppm = 10,
  rt.tol = Inf,
  selected_adduct = MSCC::adduct.table$Adduct,
  ...
)
```

## Arguments

- xcms.xcms:

  XCMSnExp object with feature definitions.

- cpdb:

  CompoundDb object containing compound database.

- mz.ppm:

  Numeric. Mass accuracy tolerance in parts per million (default 10).

- rt.tol:

  Numeric. Retention time tolerance in seconds (default Inf, no RT
  filtering).

- selected_adduct:

  Character vector of adducts to consider (default from
  MSCC::adduct.table\$Adduct).

- ...:

  Additional arguments passed to internal functions.

## Value

XCMSnExp object with featureDefinitions updated with candidate.id,
candidate.formula, candidate.adduct, and candidate.mz columns.

## Functions

- `xcms_get_feature_ms1_candidate()`: match features to compound
  database
