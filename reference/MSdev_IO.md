# MSdev input and output

Save `MSdev` with
[`qs2::qs_save()`](https://rdrr.io/pkg/qs2/man/qs_save.html) / load with
[`qs2::qs_read()`](https://rdrr.io/pkg/qs2/man/qs_read.html).

Load an MSdev object from a file (tries
[`qs2::qs_read()`](https://rdrr.io/pkg/qs2/man/qs_read.html), legacy
`qs::qread()` if installed,
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
