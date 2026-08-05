# Store MS1 feature purity matrix in XcmsExperiment qdata

Compute a feature-by-sample MS1 purity matrix and store it as assay
`"purity_matrix"` in `MsExperiment::qdata(xcms.xcms)`.

Requires an `XcmsExperiment`. MS1 spectra are taken from
`ProtGenerics::spectra(xcms.xcms)` filtered with
`Spectra::filterMsLevel(1L)`. The matrix is computed by
[`get_xcms_feature_purity_matrix`](https://drruili.github.io/MSdev/reference/get_xcms_feature_purity_matrix.md).
If `qdata` is missing, it is created with
[`xcms::quantify()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
(assay `"raw"`), then assay `"purity_matrix"` is added. Does not write
aggregated purity into `featureDefinitions`.

## Usage

``` r
xcms_get_feature_purity(xcms.xcms, ppm = 10, isolation_half_window = 0.2)
```

## Arguments

- xcms.xcms:

  `XcmsExperiment` with grouped features.

- ppm:

  numeric, ppm tolerance for m/z window.

- isolation_half_window:

  numeric, half isolation window (m/z).

## Value

`XcmsExperiment` with `qdata` containing assay `"purity_matrix"`
(aligned to existing `qdata` rows/columns).

## See also

[`get_xcms_feature_purity_matrix`](https://drruili.github.io/MSdev/reference/get_xcms_feature_purity_matrix.md)
