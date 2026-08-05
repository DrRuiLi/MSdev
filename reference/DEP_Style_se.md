# [`DEP`](https://rdrr.io/pkg/DEP/man/DEP.html) styled [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html) and related analysis

DEP styled SummarizedExperiment and related analysis

DEP list contrast.

Add significance rejections to a SummarizedExperiment based on p-values
and fold changes. Reference to
[`add_rejections`](https://rdrr.io/pkg/DEP/man/add_rejections.html),
which does not support significance without p-adjustment. This function
is a supplementary implementation.

Adjust p-values for multiple testing across all contrasts.

Wrapper of [`test_diff`](https://rdrr.io/pkg/DEP/man/test_diff.html) for
differential analysis between conditions.

Filter significant features based on differential analysis results.

Get differential analysis table for a specified contrast.

Create a volcano plot for differential analysis results.

Create a volcano plot colored by lipid class for lipidomics data.

Plot log2 fold changes grouped by lipid class.

Create a heatmap of expression data using ComplexHeatmap.

Export SummarizedExperiment data to an Excel file with sample info and
compound data.

Create a PCA plot of expression data.

Perform pathway enrichment analysis using Hypergeometric test or Global
test.

Perform gene set enrichment analysis using enrichR.

Adjust sample intensities by weight column.

Perform Kruskal-Wallis test (non-parametric ANOVA) across conditions.

Create a boxplot with jitter for a single feature across conditions.

Impute missing values with the mean of each sample.

Filter features based on missing value ratio within each group.

Filter features based on QC relative standard deviation.

Calculate relative standard deviation (RSD) for QC samples.

Perform preprocessing pipeline: filter missing values, filter QC RSD,
normalize, and impute.

Plot boxplots of expression data before and after normalization.

Get a named color vector for groups in a SummarizedExperiment.

Get significant feature IDs from differential analysis results.

Import data from a MetaboExplorer (ME) result file into a
SummarizedExperiment.

Remove QC and Blank samples from a SummarizedExperiment.

## Usage

``` r
get_MSdev_DEP_se(
  object,
  from = c("metabolite.se", "feature.se"),
  preprocess = T,
  ...
)

DEP_list_contrast(data.se)

DEP_add_rejections(data.se, p.adjust = T, p = 0.05, lfc = 0.5)

DEP_p_adjust(data.se, p.adjust.method = "fdr")

DEP_check_sig(data.se)

DEP_test_diff(se, type, ...)

DEP_filter_significant(
  data.se,
  contrast = DEP_list_contrast(data.se)[1],
  top = Inf
)

DEP_get_diff_table(
  data.se,
  contrast = DEP_list_contrast(data.se)[1],
  keep.all = F
)

DEP_plot_volcano(
  data.se,
  contrast = DEP_list_contrast(data.se)[1],
  show.label = T,
  label.top = 10,
  label.max.char = 15
)

DEP.plot.volcano.lipidomic(
  data.se,
  contrast = DEP_list_contrast(data.se)[1],
  p.adjust = F,
  show.label = T
)

DEP.plot.lfc.lipid.class(
  data.se,
  contrast = DEP_list_contrast(data.se)[1],
  p.adjust = F
)

DEP_plot_heatmap(data.se, feature_id = NULL, ...)

DEP_export_data(data.se, file_path)

DEP_plot_PCA(
  data.se,
  col.group = get_DEP_se_group_color(data.se),
  showlabel = F,
  ...
)

DEP_pathway_enrich(
  data.se,
  contrast,
  method = c("HyperTest", "GlobalTest"),
  filter_Metabolism = F
)

DEP_pathway_enrich_gene(data.se, contrast, database = c("KEGG_2021_Human"))

se_adjuset_by_weight(data.se)

DEP_test_ANOVA(data.se)

DEP_plot_single_bar(data.se, id)

DEP_impute_mean(data.se)

DEP_filter_miss(data.se, group.miss.ratio = 0.3)

DEP_filter_QC_RSD(data.se, QC_RSD = 0.3)

DEP_get_QC_RSD(data.se)

DEP_preprocess(
  data.se,
  group.miss.ratio = 0.3,
  QC_RSD = 0.3,
  keep_before_norm = F
)

DEP_plot_normalization(se, ...)

get_DEP_se_group_color(se)

get_DEP_se_sig_feature(data.diff, contrast = DEP_list_contrast(data.diff)[1])

get_DEP_se_from_ME_result(ME_file)

DEP_remove_QC(data.se, remove_QC = T, remove_Blank = T)
```

## Arguments

- ...:

  Additional SummarizedExperiment objects to include in the plot.

- data.se:

  A SummarizedExperiment object with sample.type column.

- p.adjust:

  Logical indicating whether to use adjusted p-values (default FALSE).

- p:

  Numeric threshold for p-value (or adjusted p-value if p.adjust=TRUE).

- lfc:

  Numeric threshold for absolute log2 fold change.

- p.adjust.method:

  Character string specifying the p-value adjustment method (default
  "fdr").

- se:

  A SummarizedExperiment object with group column.

- type:

  Character string specifying the type of test (default "all").

- contrast:

  Character string specifying the contrast (default first contrast, or
  "all" for all contrasts).

- top:

  Integer specifying the maximum number of significant features to
  return (default Inf).

- keep.all:

  Logical indicating whether to keep all rowData columns (default
  FALSE).

- show.label:

  Logical indicating whether to label significant features (default
  TRUE).

- label.top:

  Integer specifying the number of top features to label (default 10).

- label.max.char:

  Integer specifying maximum character length for labels (default 15).

- feature_id:

  Optional character vector of feature IDs to include (default NULL for
  all).

- file_path:

  Character string specifying the file path for the Excel output.

- col.group:

  Named character vector of colors for groups (default from
  `get_DEP_se_group_color`).

- showlabel:

  Logical indicating whether to show sample labels (default FALSE).

- method:

  Character string specifying enrichment method: "HyperTest" or
  "GlobalTest".

- filter_Metabolism:

  Logical indicating whether to filter to metabolism pathways only
  (default FALSE).

- database:

  Character string specifying the enrichR database (default
  "KEGG_2021_Human").

- id:

  Character string specifying the feature ID to plot.

- group.miss.ratio:

  Numeric threshold for missing value ratio per group (default 0.3).

- QC_RSD:

  Numeric threshold for QC RSD filtering (default 0.3).

- keep_before_norm:

  Logical indicating whether to keep data before normalization as an
  additional assay (default FALSE).

- data.diff:

  A SummarizedExperiment object with significance results.

- ME_file:

  Character string specifying the path to the ME Excel file.

- remove_QC:

  Logical indicating whether to remove QC samples (default TRUE).

- remove_Blank:

  Logical indicating whether to remove Blank samples (default TRUE).

## Value

A character vector of contrast names (e.g., "condA_vs_condB").

A SummarizedExperiment object with added significance columns.

A SummarizedExperiment object with updated p.adjust columns.

A SummarizedExperiment object with differential analysis results.

A SummarizedExperiment object containing only significant features.

A data frame with differential analysis results (log2 fold change,
p-values, significance).

A ggplot2 volcano plot object.

A ggplot2 volcano plot object colored by lipid class.

A ggplot2 plot showing log2 fold changes grouped by lipid class.

A ComplexHeatmap object.

Invisible NULL. The function writes an Excel file as a side effect.

A ggplot2 PCA plot object.

A data frame with pathway enrichment results, or a list of data frames
if contrast="all".

A data frame with enrichment results, or a list of data frames if
contrast="all".

A SummarizedExperiment object with adjusted assay values.

A SummarizedExperiment object with added p.kruskal and p.kruskal.fdr
columns in rowData.

A ggplot2 boxplot object.

A SummarizedExperiment object with imputed values.

A SummarizedExperiment object with filtered features.

A SummarizedExperiment object with filtered features.

A SummarizedExperiment object with added qc_rsd column in rowData.

A SummarizedExperiment object after preprocessing.

A ggplot2 boxplot showing normalization effects.

A named character vector of colors where names are group names.

A character vector of significant feature IDs.

A SummarizedExperiment object with DEP-style formatting.

A SummarizedExperiment object with QC and/or Blank samples removed.

## Functions

- `get_MSdev_DEP_se()`: get DEP style
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  from MSdev

- `DEP_list_contrast()`: list all contrast in SummarizedExperiment

- `DEP_add_rejections()`: Add significant,Reference to
  [`add_rejections`](https://rdrr.io/pkg/DEP/man/add_rejections.html),which
  not support significant with out p adjust, this function as
  supplymentary

- `DEP_p_adjust()`: multiple test

- `DEP_check_sig()`: check if `DEP_add_rejections` performed

- `DEP_test_diff()`: warpper of DEP::test_diff

- `DEP_filter_significant()`: filter significant feature

- `DEP_get_diff_table()`: get differential table

- `DEP_plot_volcano()`: plot volcano

- `DEP.plot.volcano.lipidomic()`: plot volcano with lipid class

- `DEP.plot.lfc.lipid.class()`: plot lfc-class

- `DEP_plot_heatmap()`: plot heatmap

- `DEP_export_data()`: export data, wirte coldata and rowdata to excel

- `DEP_plot_PCA()`: plot PCA

- `DEP_pathway_enrich()`: plot pathway enrich

- `DEP_pathway_enrich_gene()`: plot gene pathway enrich

- `se_adjuset_by_weight()`: adjust by weight

- `DEP_test_ANOVA()`: ANOVA test

- `DEP_plot_single_bar()`: plot bar plot for feature

- `DEP_impute_mean()`: impute with mean

- `DEP_filter_miss()`: filter feature with miss value

- `DEP_filter_QC_RSD()`: filter feature with QC RSD

- `DEP_get_QC_RSD()`: calculate RSD of QC

- `DEP_preprocess()`: filter miss, filter QC rsd, normalization,
  imputation

- `DEP_plot_normalization()`: update `plot_normalization`, add group
  color

- `get_DEP_se_group_color()`: get a vector of group color

- `get_DEP_se_sig_feature()`: get significant feature

- `get_DEP_se_from_ME_result()`: import data from ME result

- `DEP_remove_QC()`: remove QC and Blank
