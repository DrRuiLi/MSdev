# Build MS2 spectra target list (mz/rt windows)

Extract MS2 spectra from `object@spectra$MS2_Spectra` and summarise them
into a target table with `mz`, `rt`, `rtmin`, `rtmax`. If spectra are
assigned to features (`feature_id`) and feature definitions contain MSIP
annotations (`compound_id`, `iso_form`), targets can be grouped at
compound or isotopologue level.

## Usage

``` r
get_MSdev_spectra_target_list(
  object,
  prefer = c("assigned_feature", "all_ms2"),
  group_by = c("compound_iso", "compound", "feature", "none"),
  rt_expand = 0,
  mz_col = c("isolationWindowTargetMz", "precursorMz"),
  rt_col = c("rtime")
)
```

## Arguments

- object:

  MSdev object.

- prefer:

  character, either `"assigned_feature"` (default; only MS2 with
  non-missing `feature_id`) or `"all_ms2"`.

- group_by:

  character, grouping strategy. One of `"compound_iso"`, `"compound"`,
  `"feature"`, `"none"`. When requested metadata are missing, it falls
  back to the next available level.

- rt_expand:

  numeric, seconds added to both sides of the final RT window.

- mz_col:

  candidate mz columns in spectraData; first match is used. Default
  prefers `isolationWindowTargetMz`, then `precursorMz`.

- rt_col:

  candidate RT columns in spectraData; first match is used. Default is
  `rtime` (Spectra / XcmsExperiment naming).

## Value

data.frame with columns `mz, rt, rtmin, rtmax` plus grouping columns
when available.
