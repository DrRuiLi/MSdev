# MSdev

LC-MS data analysis shell for project orchestration, xcms-based feature processing, spectra handling, annotation, and downstream analysis.

Chemistry and molecule-graph primitives live in companion package [MSCC](https://github.com/DrRuiLi/MSCC); `MSdev` consumes them for network / atom-transfer workflows.

## Installation

```         
# install.packages("remotes")
remotes::install_github("DrRuiLi/MSdev")
```

Requires R ≥ 4.3 and Bioconductor packages such as `xcms`, `Spectra`, and `ProtGenerics`. Install Bioconductor dependencies first if needed:

```         
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("xcms", "Spectra", "ProtGenerics", "MsExperiment"))
```

## What it does

| Area | Capabilities |
|--------------------|----------------------------------------------------|
| Project shell | `MSdev` S4 object: sample info, xcms data, spectra, annotation, advanced analysis |
| Feature processing | Peak picking / grouping / RT alignment via xcms; polarity-aware MS1 containers |
| Spectra | Extract MS1/MS2, assign MS2 to features (`sp_id`), on-disk storage |
| Feature grouping | Stock xcms compounding (`EicSimilarityParam`) and custom EIC similarity (`MSdev_group_feature_EIC`) |
| Annotation | MS2 scoring and formula / library glue |
| Other workflows | DDA / pseudo-MS2, MRM, statistics, pathway helpers, metabolic flux network (via MSCC) |

## Typical untargeted workflow

```         
library(MSdev)

object <- MSdev(rawDataDir = "path/to/raw")

object <- MSdev_checkSampleInfo(object)
object <- MSdev_msConvert(object)
object <- MSdev_xcmsProcessing(object)
object <- MSdev_extract_Spectra(object)   # also assigns MS2 to features
object <- MSdev_group_feature_EIC(object) # optional custom EIC grouping
object <- MSdev_annotation(object)
```

Prefer accessors over slot digging:

- `get_MSdev_Spectra(object, msLevel = …, polarity = …)`
- `get_MSdev_Chromatogram(object, polarity = …)`

## Documentation

Package vignettes and articles (pkgdown):

- Get started: `vignettes/MSdev.Rmd`
- Untargeted workflow notes: `vignettes/MSdev_untargeted_workflow.Rmd`
- Articles under `vignettes/articles/`:
  - Spectra backends (`MsBackend*`)
  - Fast chromatogram extraction (xcms vs MSdev triad)
  - Feature grouping with `EicSimilarityParam`
  - Rdisop mass decomposition (`decomposeMass`)

Build the site locally with pkgdown when configured:

```         
pkgdown::build_site()
```

## Related packages

| Package | Role |
|-------------------------------------------|-----------------------------|
| [MSCC](https://github.com/DrRuiLi/MSCC) | Chemistry, molecule graphs, adduct / formula helpers |
| [xcms](https://bioconductor.org/packages/xcms/) | Peak detection, correspondence, feature compounding |
| [Spectra](https://bioconductor.org/packages/Spectra/) | Spectrum containers and backends |

## Author

**Rui Li** ([ORCID](https://orcid.org/0000-0002-3199-287X)) — [rli\@sinh.ac.cn](mailto:rli@sinh.ac.cn){.email}

## License

MIT © Rui Li
