# Score feature MS2 candidates, parallel only for features with spectra

Features with no experimental MS2 are scored on the main process.
Features with both experimental and reference spectra are sent to
`BPPARAM` as one job each.

## Usage

``` r
.score_ms2_features(
  feature_id,
  candidate_id,
  candidate_adduct,
  sp.exp,
  sp.ref,
  BPPARAM = BiocParallel::SerialParam()
)
```
