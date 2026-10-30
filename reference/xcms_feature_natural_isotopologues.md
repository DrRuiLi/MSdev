# Natural isotopologue ladders for xcms features

Finds natural-abundance isotopologues (13C, 2H, 15N, 18O, 34S, 37Cl,
81Br, 41K) among all features of one polarity, without chemical
annotation or feature groups. Only isotopologues (different isotopic
composition, different mass) are detected; isotopomers are not
distinguishable by MS1 mass.

A feature pair is a candidate when its aggregated EIC similarity
(`otherData(xcms)$EIC_Similarity_agg`) is above `sim_threshold`. The
lighter feature is the parent; the m/z gap must match `k` times an
isotope mass difference within `ppm` (standard ppm of the mean m/z). The
abundance ceiling is `prod(choose(n, k) * abundance^k)` with
`n = floor(mz / major mass)` for C/H/N/O and `n = max_atoms` for
S/Cl/Br/K. A row is used only if the ceiling is at least `min_rel`
(detectable k for that m/z) and the observed `into` ratio is at most
ceiling times the factor (`intensity_factor_chon` or
`intensity_factor_hetero`). 41K is used only in positive mode.

Parallel-ladder check: an S/Cl/Br step whose parent is itself an
isotopologue is kept only if the monoisotopic root has the same step
(same isotope and k) to a feature with higher median `into`. Dropped
edges have `pass_parallel = FALSE` in `pairs`.

Each heavy feature keeps one parent, ordered by isotope priority (13C,
then 15N/18O/2H, then 34S/37Cl/81Br, then 41K), \|ppm\| and similarity.
Parents that are not linked to each other make the feature `ambiguous`.
Parents are followed to the monoisotopic root; the cumulative
composition is rechecked against the root (mass, `min_rel`, intensity),
otherwise the direct parent is kept (`chain_fallback`).

Constants come from `MSCC::elem_table`.

## Usage

``` r
xcms_feature_natural_isotopologues(
  xcms.xcms,
  polarity = c("Positive", "Negative"),
  sim_threshold = 0.9,
  ppm = 5,
  max_atoms = c(S = 3, Cl = 3, Br = 3, K = 1),
  intensity_factor_chon = 6,
  intensity_factor_hetero = 1.5,
  min_rel = 0.001
)
```

## Arguments

- xcms.xcms:

  xcms object with `otherData(xcms)$EIC_Similarity_agg` (see
  [`xcms_feature_EIC_similarity`](https://drruili.github.io/MSdev/reference/xcms_feature_EIC_similarity.md)).

- polarity:

  `"Positive"` or `"Negative"`.

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

Updated `xcms.xcms`. `featureDefinitions` gains `isotopologue_status`
(`"monoisotopic"`, `"isotopologue"`, `"ambiguous"`, `"chain_fallback"`
or `NA`), `isotopologue_anchor` (own id for monoisotopic roots),
`isotopologue_formula`, `isotopologue_ppm`, `isotopologue_sim`.
`otherData(xcms)$Isotopologue_Global` holds `params`, `summary`,
`pairs`, `links`, `conflicts`.
