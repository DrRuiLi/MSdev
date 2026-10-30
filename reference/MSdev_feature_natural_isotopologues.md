# Natural isotopologue ladders per polarity

Per-polarity wrapper around
[`xcms_feature_natural_isotopologues`](https://drruili.github.io/MSdev/reference/xcms_feature_natural_isotopologues.md).
Run after
[`MSdev_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/MSdev_feature_EIC_similarity.md).

## Usage

``` r
MSdev_feature_natural_isotopologues(
  object,
  sim_threshold = 0.9,
  ppm = 5,
  max_atoms = c(S = 3, Cl = 3, Br = 3, K = 1),
  intensity_factor_chon = 6,
  intensity_factor_hetero = 1.5,
  min_rel = 0.001
)
```

## Arguments

- object:

  MSdev object.

- sim_threshold:

  numeric(1). Pairs need similarity strictly above this.

- ppm:

  numeric(1). Tolerance of the m/z gap, ppm of the mean m/z.

- max_atoms:

  Named integer. Maximum S, Cl, Br, K atoms per ion.

- intensity_factor_chon:

  numeric(1). Allowed factor over the ceiling for C/H/N/O-only
  compositions.

- intensity_factor_hetero:

  numeric(1). Allowed factor over the ceiling when S/Cl/Br/K is
  involved.

- min_rel:

  numeric(1). Minimum theoretical ceiling for a k to be tried.

## Value

MSdev object; each polarity's MS1 xcms object carries the
`isotopologue_*` columns and `otherData(xcms)$Isotopologue_Global`.
