# Estimate Noise Level for Spectra

Estimates the noise level for each spectrum using a density-based
method. The noise is estimated as the intensity value at the maximum of
the log10-transformed intensity density distribution. Also calculates
signal-to-noise ratio (SNR).

## Usage

``` r
Spectra_get_noise(sp)
```

## Arguments

- sp:

  A `Spectra` object.

## Value

The input `Spectra` object with two additional spectrum variables:
`noise` (estimated noise level) and `snr` (signal-to-noise ratio,
calculated as base peak intensity divided by noise).

## Details

Spectra_get_noise
