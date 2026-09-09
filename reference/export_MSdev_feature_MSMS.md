# Export MS/MS spectrum and chromatogram for a feature

Export PNG images of the MS/MS spectrum and chromatogram for a given
feature.

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
