# MSdev — AI wayfinding (`AI_CONTEXT.md`)

Concise index for future sessions. Do not treat this file as authoritative API detail; confirm behavior in `man/` (Rd), roxygen comments, and current `R/` sources.

---

## 1. Core project goal

`MSdev` is now focused on the LC-MS processing shell around the `MSdev` S4 object: project orchestration, xcms-based feature processing, spectra handling, annotation glue, DDA/MRM utilities, network/flux wrappers, statistics, and export/report helpers.

Chemistry and molecule-graph primitives were refactored into `MSCC` and are consumed from there.

---

## 2. Refactor status and package boundaries

### What moved out

- Chemistry/structure modules (for example CFM, RXN mapper, molecule graph classes/functions, related wrappers) were migrated out of this package.
- Local generic definitions for graph accessors (`vdata`, `vdata<-`, `edata`, `edata<-`) were removed from `MSdev`.

### What remains in MSdev

- `MSdev` still defines package-specific S4 classes and workflows for feature processing/analysis.
- `MSdev` provides methods for `Metabolic_flux_network` that rely on imported `MSCC` graph generics.

### Cross-package contract

- `MSCC` defines and exports graph generics/accessors.
- `MSdev` imports them (including replace functions via raw namespace import) in `R/0_mscc_graph_generics_imports.R`.
- Method registration for `Metabolic_flux_network` accessors is done in `.onLoad()` in `R/0_mscc_graph_generics_imports.R` to avoid load-order issues.

---

## 3. Current architecture map (`R/`)

| Area | Primary files | Read first when… |
|------|---------------|------------------|
| Project shell and orchestration | `MSdev-class.R`, `MSdev-function.R`, `MSdev-workflow.R`, `Demo.R` | Understanding core object lifecycle, setup, and top-level API |
| xcms / feature processing | `dev_xcms.R`, `MSdev-function.R`, `onDiskData.R` | Debugging peak picking/grouping/RT alignment or feature extraction |
| DDA / pseudo-MS2 workflow | `DDA-function.R`, `DDA_Mine_function.R`, `DDA_mine-workflow.R`, `Pseudo-workflow.R` | Following DDA simulation/mining flow |
| MRM | `MRM-function.R`, `MRM-WorkFlow.R` | Investigating targeted chromatogram pipeline |
| Network / atom-transfer wrappers | `Metabolic_flux_network.R`, `Reaction_atom_transfer.R`, `0_mscc_graph_generics_imports.R` | Tracing graph accessors, reaction transfer, and MSCC integration points |
| Statistics and downstream analysis | `StatisticFunction.R`, `MSdev-Sta_function.R`, `dev_DEP.R`, `dev_caret.R`, `dev_MetaboSignal.R`, `dev_FELLA.R`, `dev_KEGG.R` | Fixing modeling/DE/pathway/stat output behavior |
| Visualization and plotting | `dev_plot.R` | Plot behavior and formatting |
| Utility layers | `dev_base.R`, `dev_string.R`, `dev_tidyverse.R`, `dev_math.R`, `dev_others.R`, `dev_openxlsx.R`, `dev_RStudio.R` | Shared helpers and local utilities |
| Spectra/instrument helpers | `dev_Spectra.R`, `dev_mzR.R`, `dev_MSInstrument.R`, `MS_Exp-function.R`, `MS_exp-class.R` | Spectra parsing, metadata, and instrument-specific behavior |

---

## 4. Core data flow / state

| Stage | Form |
|-------|------|
| Input | Raw files under `rawDataDir`; `sampleInfo` table; optional `MS_Exp` metadata |
| Core container | `MSdev` S4 with slots `projectInfo`, `processingInfo`, `sampleInfo`, `experimentInfo`, `xcmsData`, `spectra`, `annotation`, `advancedAna` |
| Typical chain | `MSdev_checkSampleInfo` -> `MSdev_msConvert` -> `MSdev_xcmsProcessing` -> spectra extraction/matching -> annotation -> stats/export -> save/load |
| Outputs | Serialized project objects and tabular/plot exports |

---

## 5. Editing guidance for future sessions

1. Read `AI_CONTEXT.md` first, then open the smallest domain file that owns the behavior.
2. For graph/generic issues, check `R/0_mscc_graph_generics_imports.R` and `NAMESPACE` imports from `MSCC` before touching network files.
3. Avoid broad edits in large utility files (especially `dev_xcms.R`) unless caller boundaries are clear.
4. If architecture shifts again (new migrations between `MSdev` and `MSCC`), update this file only after explicit user request.

---

Last refreshed after refactor to `MSCC` ownership for graph/chemistry primitives. This file is a navigation aid, not a substitute for reading source.
