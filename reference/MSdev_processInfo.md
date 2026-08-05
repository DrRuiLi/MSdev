# List process history of xcmsData

List all process history entries for a specific xcmsData element in an
MSdev object.

## Usage

``` r
MSdev_processInfo(object, target = "PositiveMS1")
```

## Arguments

- object:

  MSdev object

- target:

  character. Name of the xcmsData element (default "PositiveMS1").

## Value

data.frame of process history, or NULL if not available.
