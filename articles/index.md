# Articles

### Get started

- [MSdev untargeted
  workflow](https://drruili.github.io/MSdev/articles/MSdev_untargeted_workflow.md):

  End-to-end untargeted LC-MS (lipidomics-style) pipeline: convert, xcms
  processing, spectra extraction, CompoundDb annotation, stats, export,
  and differential figures by sample.type.

- [Get started with
  MSdev](https://drruili.github.io/MSdev/articles/MSdev.md):

  Overview of the MSdev S4 project object, core slots, preferred
  accessors, and a minimal untargeted LC-MS pipeline.

### Technical notes

- [dev_xcms](https://drruili.github.io/MSdev/articles/dev_xcms.md):

- [Rdisop mass decomposition
  (\`decomposeMass\`)](https://drruili.github.io/MSdev/articles/Rdisop_decomposeMass.md):

  Quick reference for how Rdisop turns an exact mass (or isotope
  pattern) into candidate sum formulas, and how that fits a
  formula-calculator workflow.

- [Spectra backends
  (\`MsBackend\*\`)](https://drruili.github.io/MSdev/articles/Spectra_backends.md):

  Quick reference for choosing and using Spectra backends in MSdev /
  Bioconductor LC-MS workflows.

- [Fast chromatogram extraction (xcms vs MSdev
  triad)](https://drruili.github.io/MSdev/articles/xcms_chromatogram_extraction.md):

  How xcms::chromatogram / featureChromatograms / chromPeakChromatograms
  extract EICs, why that path is slow at scale, and how MSdev’s
  get_xcms\_\*\_chromatogram triad speeds up the same work.

- [Feature grouping with
  \`EicSimilarityParam\`](https://drruili.github.io/MSdev/articles/xcms-feature-group-EicSimilarityParam.md):

  xcms / MsFeatures feature compounding (RT → abundance → EIC shape),
  focusing on EIC extraction and the similarity matrix, and how MSdev’s
  stock path compares to custom MSdev_group_feature_EIC().

### Developer notes
