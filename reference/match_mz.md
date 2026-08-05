# Match Ions by M/Z Only

Matches two lists of ions based solely on mass-to-charge ratio (m/z).
Returns the closest match for each ion in the first list, handling
multiple potential matches.

## Usage

``` r
match_mz(mz1, mz2, mz.ppm = 10)
```

## Arguments

- mz1:

  Numeric vector of m/z values for the first ion list to match.

- mz2:

  Numeric vector of m/z values for the second ion list to match against.

- mz.ppm:

  Numeric tolerance for m/z matching in parts per million (ppm). Default
  is 10.

## Value

A numeric vector of the same length as mz1, containing indices of
matched ions in mz2. NA values indicate no match was found within the
specified tolerance.
