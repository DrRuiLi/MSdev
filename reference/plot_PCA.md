# Plot Principal Component Analysis (PCA)

Performs PCA and creates a scatter plot with optional ellipses.

## Usage

``` r
plot_PCA(
  pca.matrix,
  pca.group,
  force_ellipse = T,
  show_ellipse = T,
  showlabel = F,
  col = NULL
)
```

## Arguments

- pca.matrix:

  A numeric matrix of features (samples as rows).

- pca.group:

  Factor or vector indicating group membership for each sample.

- force_ellipse:

  Logical, whether to force ellipse calculation for groups with few
  samples.

- show_ellipse:

  Logical, whether to show confidence ellipses.

- showlabel:

  Logical, whether to label points with sample names.

- col:

  Optional named vector of colors for groups.

## Value

A ggplot object.
