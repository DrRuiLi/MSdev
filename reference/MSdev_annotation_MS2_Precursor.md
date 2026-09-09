# Annotate MS2 precursor peaks using a compound database

Annotate `object@advancedAna$MS2_Precursor` (from
[`MSdev_get_peak_table_from_spectra`](https://drruili.github.io/MSdev/reference/MSdev_get_peak_table_from_spectra.md))
with CompDb MS1 candidates and MS2 spectral scores. Same CompDb engine
as
[`MSdev_annotation`](https://drruili.github.io/MSdev/reference/MSdev_workflow.md),
but the targets are spectra-derived precursor peaks, not xcms MS1
features. Isotope-pattern scoring is skipped (no MS1 intensity matrix).
Requires `ms2_id` on the peak table linking to `MS2_Spectra` `sp_id` /
spectraNames.

## Usage

``` r
MSdev_annotation_MS2_Precursor(
  object,
  cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite",
  ppm = 10,
  weight_mz = 0.2,
  weight_ms2 = 0.8,
  ...
)
```

## Arguments

- object:

  MSdev object with `advancedAna$MS2_Precursor`

- cpdb_path:

  path to CompoundDb SQLite database

- ppm:

  m/z tolerance in parts per million for candidate matching

- weight_mz:

  weight for m/z error score (default `0.2`)

- weight_ms2:

  weight for MS2 similarity score (default `0.8`)

- ...:

  additional arguments passed to annotation helpers

## Value

MSdev object with annotated `advancedAna$MS2_Precursor`

## Details

`MS2_Precursor` is not from xcms peak picking. It is built by grouping
MS2 spectra by precursor m/z (`ppm`) and RT gap (`rt_tol`) into rows
such as `MS2P000001`, each with `mzmed`/`rtmed` and `ms2_id` pointing at
those spectra.

For each polarity the function:

1.  matches precursor `mzmed` to CompDb adduct m/z
    (`fdf_get_ms1_candidate`);

2.  scores experimental MS2 vs CompDb reference spectra
    (`fdf_get_ms2_score`, ndotproduct);

3.  sets `score.isopattern` to zeros (no MS1 intensity matrix);

4.  picks the best candidate with default weights `weight_mz = 0.2`,
    `weight_ms2 = 0.8`, `weight_isopattern = 0`.

Results (`compound_id`, `adduct`, `score`, CompDb name/formula/smiles,
...) are written back to `advancedAna$MS2_Precursor`.
`projectInfo$CompoundDB_path` is also set. If `MS2_Precursor` is
missing, a message asks to run
[`MSdev_get_peak_table_from_spectra`](https://drruili.github.io/MSdev/reference/MSdev_get_peak_table_from_spectra.md)
first and the object is returned unchanged.

Difference from
[`MSdev_annotation`](https://drruili.github.io/MSdev/reference/MSdev_workflow.md):

|  |  |  |
|----|----|----|
|  | `MSdev_annotation` | `MSdev_annotation_MS2_Precursor` |
| What | xcms MS1 features (`PositiveMS1` / `NegativeMS1`) | MS2 precursor groups (`advancedAna$MS2_Precursor`) |
| Written to | `featureDefinitions` on the xcms object | `advancedAna$MS2_Precursor` |
| MS2 link | `ms2_id` from [`MSdev_assign_MS2`](https://drruili.github.io/MSdev/reference/MSdev_workflow.md) | `ms2_id` from grouping the MS2 spectra themselves |
| Isotope score | optional (`calc_isopattern_score`) | never |
| Default weights | mz 0.1 / MS2 0.7 / iso 0.2 | mz 0.2 / MS2 0.8 / iso 0 |
| Helpers | `xcms_get_feature_*` | `fdf_get_*` (same logic on a data.frame) |

The two functions do not overwrite each other.

## See also

[`MSdev_annotation`](https://drruili.github.io/MSdev/reference/MSdev_workflow.md),
[`MSdev_get_peak_table_from_spectra`](https://drruili.github.io/MSdev/reference/MSdev_get_peak_table_from_spectra.md)
