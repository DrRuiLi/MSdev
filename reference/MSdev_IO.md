# MSdev input and output

save `MSdev` using `qs::qsave()` and qs::qread()\`

Load an MSdev object from a file (tries `qs::qread()`,
[`readRDS()`](https://rdrr.io/r/base/readRDS.html), then
[`load()`](https://rdrr.io/r/base/load.html)). Use
[`.update_MSdev_object`](https://drruili.github.io/MSdev/reference/dot-update_MSdev_object.md)
on objects saved with older MSdev versions.

## Usage

``` r
MSdev_save(object, file = object@projectInfo$MSdevFile)

MSdev_load(file_to_load)
```

## Arguments

- object:

  MSdev

- file_to_load:

  file path

## Value

MSdev

MSdev

## Functions

- `MSdev_save()`: MSdev_save

  this function ....

- `MSdev_load()`: load
