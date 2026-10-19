# Resolve an MSdev project path

Return an absolute path for a key in `projectInfo$paths`. Paths are
stored relative to `projectDir`. Old objects that have no `paths` entry
are resolved from the existing absolute fields.

## Usage

``` r
MSdev_path(object, key)
```

## Arguments

- object:

  MSdev object

- key:

  One of `rawData`, `msData`, `msdevFile`, `ms1`, `ms2`.

## Value

Absolute path, with forward slashes.
