# Update legacy MSdev object to current slot layout

Migrate `MSdev` objects saved with older package versions (e.g. missing
`advancedAna` slot, `statData` attribute) to the current S4 definition.

## Usage

``` r
.update_MSdev_object(object)
```

## Arguments

- object:

  MSdev object

## Value

Updated MSdev object with all current slots populated
