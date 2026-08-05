# Estimate Precursor Purity for Mass Spectra

Estimates the purity of the precursor ion for MS/MS spectra. For MS2
spectra, calculates the proportion of precursor intensity within the
isolation window that belongs to the targeted precursor m/z. For MS1
spectra, estimates purity using an associated MS1 scan.

## Usage

``` r
Spectra_get_purity(sp, msLevel = 2, sp.ms1 = NULL)
```

## Arguments

- sp:

  A `Spectra` object for which to estimate purity.

- msLevel:

  Integer. MS level to process. `1` for MS1 spectra (requires `sp.ms1`),
  `2` for MS2 spectra. Default is `2`.

- sp.ms1:

  Optional. A `Spectra` object containing MS1 spectra used for purity
  calculation when `msLevel = 1`. If `NULL`, MS1 spectra are imported
  from the raw data files.

## Value

The input `Spectra` object with additional spectrum variables: For MS2:
`ms2_purity` (proportion of precursor in isolation window). For MS1:
`ms1.intensity`, `ms2.purity`, and `ms2.ppm`.

## Details

Spectra_get_purity
