# MSdev — AI wayfinding (`AI_CONTEXT.md`)

Concise index for future sessions. **Do not treat this file as authoritative API detail**; confirm behavior in `man/` (Rd) and roxygen in `R/`.

---

## 1. Core project goal

**MSdev** is an R/Bioconductor-oriented toolkit for **untargeted and targeted LC–MS workflows**: raw-to-feature processing (xcms), spectral handling, annotation, statistics, and **stable-isotope / isotopologue analysis** (PAVE / TRACE-style networks, MSIP targeted isotopologues, CFM-linked atom mapping). The package title/description frame it as **MS data manipulation and advanced algorithms** built around a central **`MSdev` S4 project object**.

---

## 2. Core architecture & file index

| Area | Primary `R/` files | Read first if… |
|------|---------------------|----------------|
| **Project shell & I/O** | `MSdev-class.R` (slots), `MSdev-function.R` (`MSdev()`, `MSdev_save`/`MSdev_load`, export, spectra, annotation orchestration) | Wrong/missing slots, save path, project dirs, `sampleInfo` sync |
| **Reference workflow (glue)** | `MSdev-workflow.R` | Understanding intended **call order** (demo-style pipelines) |
| **xcms / MS1 core** | `MSdev-function.R` (`MSdev_xcmsProcessing`, DDA/MRM dispatch), **`dev_xcms.R`** (large: `xcmsProcessingMS1`, filters, helpers) | Peak picking, grouping, RT correction, polarity split |
| **DDA / MS2 linkage** | `DDA-function.R`, `DDA_Mine_function.R`, `DDA_mine-workflow.R`, `Pseudo-workflow.R` | DDA simulation, MS2–feature assignment, mining workflows |
| **MRM** | `MRM-function.R`, `MRM-WorkFlow.R` | MRM chromatograms, MRM-specific xcms path |
| **Experiment metadata** | `MS_Exp-class.R`, `MS_Exp-function.R` | Chromatography / instrument metadata object |
| **PAVE (dual-label C/N tracing)** | `PAVE.R` (`get_PAVE_from_MSdev`, atom counting), `PAVE2.R`, `PAVE2-fun.R` | Tracer columns, PAVE1 vs PAVE2 logic, `advancedAna$PAVE` |
| **TRACE** | `TRACE.R` | Alternative/heavy **network-based** isotope tracing on xcms features |
| **MSIP (targeted isotopologues)** | `MSIP-function.R`, `MSIP_xcms_processing.targeted` (see exports), `MSIP-foundation.R`, `MSIP-scoring.R`, `MSIPCoreData.R`, `MSIPCoreData_function.R`, `MSIPFragmentMap.R`, `MSIPAtomMap.R`, `MSIP_Isotopomer_SE.R`, `MSIP-class.R` | Compound table → isotopologues, CFM, fragment/atom maps, `SummarizedExperiment` outputs |
| **CFM / structure** | `dev_CFM.R` (and CFM-related MSIP helpers) | CFM prediction, annotation pipelines |
| **Molecule graphs** | `Molecule_igraph-class.R`, `Molecule_igraph-functions.R`, `Molecule_vis.R` | Graph chemistry, visualization |
| **Reaction / flux (network)** | `Reaction_atom_transfer.R`, `Metabolic_flux_network.R` | Atom transfer matrices, network-level analyses |
| **Statistics** | `StatisticFunction.R`, `MSdev-Sta_function.R`, **`dev_DEP.R`** (DEP wrappers), `dev_caret.R` | PCA/ANOVA/diff, DEP proteomics-style stats |
| **Shiny (MSIP)** | `Shiny_MSIP.R`, `Shiny_MSIP_UI.R`, `Shiny_MSIP_SERVER.R`, `Shiny_MSIP_FUN.R`, `MSIP_server.R` | Interactive UI bugs, session logic |
| **On-disk / large data** | `onDiskData.R` | Memory-efficient data paths |
| **Generics / ChemmineR** | `0_class.R` | `atom()`, `vdata`/`edata` generics on graph/SDF-like objects |
| **Utilities (broad)** | `dev_base.R`, `dev_string.R`, `dev_plot.R`, `dev_tidyverse.R`, `dev_openxlsx.R`, `dev_mzR.R`, `dev_Spectra.R`, `dev_xcms.R`, … | **Generic helpers only** — confirm caller in higher-level module before changing behavior |

**Vignettes (`vignettes/`)** — narrative entry points: `MSdev.Rmd`, `MSdev_untargeted_workflow.Rmd`, `MSIP_Workflow.Rmd`, `MSIP_note.Rmd`, `MSIPCoreData.Rmd`, `MSIPAtomMap.Rmd`, `MSIP_Structural_Elucidation_Evaluation.Rmd`.

**Package metadata:** root `DESCRIPTION`, `NAMESPACE` (exports).

---

## 3. Core data flow / state

| Stage | Form |
|-------|------|
| **Input** | Raw vendor files under `rawDataDir` (see `MSdev()`); **`sampleInfo`** data frame (paths, polarity, sample types); optional **`MS_Exp`** in `experimentInfo` |
| **Core container** | **`MSdev`**: `projectInfo`, `sampleInfo`, `experimentInfo`, **`xcmsData`** (named list, typically `PositiveMS1` / `NegativeMS1` `XCMSnExp`), **`spectra`**, **`annotation`**, **`advancedAna`** (MSIP, PAVE, TRACE outputs, etc.) |
| **Typical untargeted chain** | `MSdev_checkSampleInfo` → `MSdev_msConvert` → `MSdev_xcmsProcessing` → spectra extraction/matching → annotation → stats → export/save (see `MSdev-workflow.R`) |
| **Outputs** | **`qs`**-serialized object (`MSdev_save` / `MSdev_load`), tables/plots/Excel from export and stat functions; MSIP may yield **`MSIPIsotopologueData`** / related S4 (`MSIP-class.R`) |

---

## 4. Interaction guidelines (for AI agents)

1. **Read `man/*.Rd` first** for exported functions (this repo builds **~100+** help pages from roxygen). If Rd is stale vs `R/`, treat **source roxygen** as ground truth and suggest `devtools::document()`.
2. **Avoid drive-by edits to `dev_xcms.R` and other large `dev_*` layers** unless the bug is localized; prefer the smallest surface (`MSdev-function.R`, domain file) after reading callers.
3. **Respect S4 slots** on `MSdev`, `MSIPMetaboliteData`, `MSIPIsotopologueData`, `MS_Exp`, etc.; extend via documented patterns rather than ad hoc slot names.
4. **Do not modify `AI_CONTEXT.md` autonomously** unless the user explicitly asks; propose updates in chat if architecture drifts.

---

*Last generated from repository layout (`R/` ≈ 63 files) + `DESCRIPTION` / `NAMESPACE`; not a substitute for reading the modules above.*
