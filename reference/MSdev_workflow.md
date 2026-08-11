# Add new sample files to MSdev object

Add raw data files from a directory to the MSdev object, converting them
as needed.

Perform peak detection and grouping using xcms functions based on
acquisition mode (FS/DDA/MRM).

Read raw MS data files and create xcms experiments for both polarities,
storing them in `object@xcmsData`. Also extracts and stores MS1 and MS2
spectra in `object@spectra`.

Set custom xcms parameters for peak detection and grouping in an MSdev
object. These parameters will be used during xcms processing.

Plot a heatmap using ComplexHeatmap showing whether each feature's peak
is detected (present) or not detected (absent) in each sample.

Open the sample information data frame in Excel for manual editing, then
update the MSdev object.

Read an edited sample-info Excel file into an MSdev object, then refresh
project metadata and (when present) xcms `pData`.

Convert raw data files to mzML format using MSconvertR, updating sample
information.

File-based spectra importer for cases without (or after) xcms init.
Reads raw data files, splits by MS level, optionally evaluates
noise/purity, and stores on-disk data in `object@spectra`.

On the main FS/DDA path, prefer `MSdev_xcmsProcessing` /
`MSdev_get_xcms`, which copy Spectra from the freshly read
`MsExperiment` (filter / `setBackend` only; xcms data unchanged) so
files are not read twice. Use this function when xcms is absent,
skipped, or already locked (e.g. after `MSdev_add_sample`).

Feature matching is not performed here; call `MSdev_assign_MS2` after
xcms feature definitions exist (also done at the end of
`MSdev_xcmsProcessing`).

Match MS2 spectra in `object@spectra$MS2_Spectra` to xcms features using
precursor m/z and retention time tolerances. Stores character `sp_id`
vectors in `featureDefinitions$ms2_id` and writes `feature_id` onto the
MS2 spectra.

Perform feature annotation using a CompoundDb database, including MS1
candidate search, MS2 scoring, and isotope pattern scoring. MS2 spectra
are taken from `object@spectra$MS2_Spectra` and selected via character
`featureDefinitions$ms2_id` (`sp_id` / spectraNames).

Annotate `object@advancedAna$MS2_Precursor` (from
[`MSdev_get_peak_table_from_spectra`](https://drruili.github.io/MSdev/reference/MSdev_get_peak_table_from_spectra.md))
with CompDb MS1 candidates and MS2 spectral scores. Isotope-pattern
scoring is skipped (no MS1 intensity matrix). Requires `ms2_id` on the
peak table linking to `MS2_Spectra` `sp_id` / spectraNames.

Build `advancedAna$feature.se` via `MSdev_get_Se`, then retrieve
compound information, filter based on scores, and optionally build
candidate and metabolite SummarizedExperiment objects.

Extract feature intensity data from xcms processing results, store as
`advancedAna$feature.se`, without compound database lookup or
metabolite/candidate filtering.

## Usage

``` r
MSdev_add_sample(object, raw.data.dir = object@projectInfo$rawDataDir)

MSdev_xcmsProcessing(
  object,
  BPPARAM = BiocParallel::SnowParam(workers = 4, progressbar = TRUE),
  ...
)

MSdev_get_xcms(object)

MSdev_set_param(object, findChromPeaks = NULL, groupChromPeaks = NULL)

plot_MSdev_sample_peaks(object, target = "PositiveMS1", top_n = Inf)

MSdev_checkSampleInfo(object, interactive = TRUE)

MSdev_import_sampleinfo(object, file = NULL, sheet = 1)

MSdev_msConvert(
  object,
  format.to = "mzML",
  BPPARAM = BiocParallel::SnowParam(workers = max(1L, parallel::detectCores() - 1L),
    progressbar = TRUE)
)

MSdev_extract_Spectra(object, rt.tol = 10, eval.noise = F, eval.ms1 = F)

MSdev_assign_MS2(object, rt.tol = 10, ppm = 20)

MSdev_annotation(
  object,
  cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite",
  calc_isopattern_score = F,
  ppm = 10,
  BPPARAM = SerialParam(progressbar = T),
  ...
)

MSdev_annotation_MS2_Precursor(
  object,
  cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite",
  ppm = 10,
  weight_mz = 0.2,
  weight_ms2 = 0.8,
  ...
)

MSdev_get_Stat(
  object,
  keys = c("name", "formula", "kegg_id", "inchikey", "lipidclass"),
  score_thresh = 0.5,
  rt_bin = NA,
  polarity_paired = T,
  candi = F,
  metabolite = T,
  ...
)

MSdev_get_Se(object, polarity_paired = TRUE, ...)
```

## Arguments

- object:

  MSdev object

- raw.data.dir:

  file path

- BPPARAM:

  BiocParallel backend passed to `MSconvertR::msConvert`.

- ...:

  additional arguments passed to
  [`get_xcms_feature_se`](https://drruili.github.io/MSdev/reference/xcms_extension_feature.md)

- findChromPeaks:

  xcms parameter object for peak detection. If `NULL` (default), uses
  the value from `get_MSdev_param(object)`.

- groupChromPeaks:

  xcms parameter object for peak grouping. If `NULL` (default), uses the
  value from `get_MSdev_param(object)`.

- target:

  character. xcmsData element: "PositiveMS1" or "NegativeMS1"
  (shorthand: "pos"/"neg").

- top_n:

  integer. Maximum number of features to plot (default Inf). Features
  are sorted by detection rate.

- interactive:

  Logical; if TRUE (default) use interactive Excel editing. If FALSE,
  write `sample.info.xlsx` under `object@projectInfo$projectDir`, then
  ask the user to edit it and press ENTER to confirm (and read it back).

- file:

  path to the Excel file. Default:
  `file.path(object@projectInfo$projectDir, "sample.info.xlsx")`.

- sheet:

  sheet name or index passed to
  [`openxlsx::read.xlsx`](https://rdrr.io/pkg/openxlsx/man/read.xlsx.html)
  (default `1`).

- format.to:

  target format (default "mzML")

- rt.tol:

  retention time tolerance (seconds)

- eval.noise:

  logical, whether to evaluate noise in MS2 spectra

- eval.ms1:

  logical, whether to evaluate purity using MS1 scans

- ppm:

  m/z tolerance in parts per million for candidate matching

- cpdb_path:

  path to CompoundDb SQLite database

- calc_isopattern_score:

  logical, whether to calculate isotope pattern scores

- weight_mz:

  weight for m/z error score (default `0.2`)

- weight_ms2:

  weight for MS2 similarity score (default `0.8`)

- keys:

  character vector of compound database keys to retrieve (e.g., "name",
  "formula")

- score_thresh:

  minimum annotation score threshold

- rt_bin:

  retention time binning width (seconds)

- polarity_paired:

  logical; if `TRUE`, retain only samples present in both polarities
  (column intersection). If `FALSE`, missing polarity objects are filled
  with an empty `SummarizedExperiment` and positive and negative feature
  sets are row-bound without requiring shared samples.

- candi:

  logical, whether to include all candidates in output

- metabolite:

  logical, whether to filter for metabolite features only

## Value

MSdev

MSdev

MSdev object with xcmsData and spectra populated

MSdev object with updated parameters

ComplexHeatmap object

MSdev a MSdev object

MSdev object with updated `sampleInfo`

MSdev object with converted files

MSdev object with spectra stored in `object@spectra`

MSdev object with updated feature-spectra assignments

MSdev object with annotation results

MSdev object with annotated `advancedAna$MS2_Precursor`

MSdev object with advancedAna populated

MSdev object with `advancedAna$feature.se` populated.

## Details

Per polarity, matching is delegated to `get_Spectra_ms2_feature_id`:

1.  Candidate pairs via
    [`match_mz_rt`](https://drruili.github.io/MSdev/reference/match_mz_rt.md)
    on feature `mzmed`/`rtmed` vs MS2 `precursorMz`/`rtime` (`ppm`,
    `rt.tol`).

2.  For each MS2, keep the feature with smallest RT error to `rtmed`,
    breaking ties by m/z error (not by m/z alone).

This is median-based matching, not
[`xcms::featureSpectra`](https://rdrr.io/pkg/xcms/man/featureSpectra.html)
peak-box matching. Progress is reported as counts/percentages of MS2
assigned and features with at least one MS2.

## Functions

- `MSdev_add_sample()`: add samples

- `MSdev_xcmsProcessing()`: use xcms to Processing data

- `MSdev_get_xcms()`: Initialize xcms data

- `MSdev_set_param()`: set xcms parameters

- `plot_MSdev_sample_peaks()`: plot peak presence heatmap across samples

- `MSdev_checkSampleInfo()`: manually check sampleInfo using excel

- `MSdev_import_sampleinfo()`: import sampleInfo from Excel

- `MSdev_msConvert()`: convert raw files

- `MSdev_extract_Spectra()`: Extract all spectra, split to MS1 and MS2,
  store as onDiskData

- `MSdev_assign_MS2()`: assign MS2 spectra to features

- `MSdev_annotation()`: annotation

- `MSdev_annotation_MS2_Precursor()`: annotate MS2 precursor peaks

- `MSdev_get_Stat()`: extract statistical data

- `MSdev_get_Se()`: extract feature SummarizedExperiment

## See also

`get_Spectra_ms2_feature_id`,
[`match_mz_rt`](https://drruili.github.io/MSdev/reference/match_mz_rt.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a new MSdev object
msdev <- MSdev(rawDataDir = "path/to/raw/data")

# Set custom CentWave parameters
cwp <- xcms::CentWaveParam(
  ppm = 10,
  peakwidth = c(5, 20),
  snthresh = 100,
  prefilter = c(3, 100)
)

# Set custom grouping parameters
gpp <- xcms::PeakDensityParam(
  sampleGroups = "A",
  bw = 5,
  minFraction = 0.6,
  binSize = 0.015
)

# Apply parameters to MSdev object
msdev <- MSdev_set_param(msdev, findChromPeaks = cwp, groupChromPeaks = gpp)
} # }
```
