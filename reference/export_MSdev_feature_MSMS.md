# Export MS/MS spectrum for a feature

Export a PNG of the MS/MS spectrum for a given feature (experimental vs
reference mirror when a reference spectrum is present).

## Usage

``` r
export_MSdev_feature_MSMS(MSdev.obj, feature_id, out.dir, cpdb = NULL)
```

## Arguments

- MSdev.obj:

  MSdev object

- feature_id:

  Character string specifying the feature ID

- out.dir:

  Output directory path

- cpdb:

  Optional `CompDb` used to fetch reference MS2 spectra.

## Value

NULL (writes files to disk)
