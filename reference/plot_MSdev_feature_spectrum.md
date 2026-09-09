# Plot MS/MS spectrum for a feature

Plot experimental and reference MS/MS spectra for a given feature, with
annotation details. Experimental MS2 is taken from legacy `annotation`
when present, else from `ms2_id` / `MS2_Spectra`, else matched from xcms
MS2 spectra. Reference MS2 is taken from legacy `refSpec` when present,
else from CompDb spectra whose precursor m/z matches the feature
(`compound_id` is preferred when annotated).

## Usage

``` r
plot_MSdev_feature_spectrum(MSdev.obj, feature.id, cpdb = NULL)
```

## Arguments

- MSdev.obj:

  MSdev object

- feature.id:

  Character string specifying the feature ID

- cpdb:

  Optional `CompDb` object. Opened from
  `object@projectInfo$CompoundDB_path` when `NULL`.

## Value

Invisible `TRUE` if a spectrum was plotted, `FALSE` otherwise.
