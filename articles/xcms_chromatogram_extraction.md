# Fast chromatogram extraction (xcms vs MSdev triad)

Extracted-ion chromatograms (EICs) are central to feature inspection,
EIC-based grouping (`EicSimilarityParam` / `MSdev_group_feature_EIC`),
and plotting. xcms provides three related entry points; MSdev mirrors
them with a faster shared engine.

This article summarizes:

1.  How the **xcms** extraction stack works  
2.  How the **MSdev triad** is layered  
3.  **Where** the speedup comes from (and what it cannot change)

Implementation lives in `R/dev_xcms.R`. Project-level caching uses
`MSdev_get_feature_chrom()`.

------------------------------------------------------------------------

## 1. The three extraction problems

| Question | xcms | MSdev |
|----|----|----|
| Extract EICs for arbitrary mz–RT boxes | `chromatogram()` | `get_xcms_chromatogram()` |
| Exact EIC for each chromatographic peak (owning sample) | `chromPeakChromatograms()` | `get_xcms_peaks_chromatogram()` |
| One shared box per feature, applied to selected samples | `featureChromatograms()` | `get_xcms_feature_chromatogram()` |

``` mermaid
flowchart TD
  raw[Raw spectra files]
  raw --> chrom["chromatogram / get_xcms_chromatogram"]
  chrom --> peaks["chromPeakChromatograms / get_xcms_peaks_chromatogram"]
  chrom --> feats["featureChromatograms / get_xcms_feature_chromatogram"]
```

Peaks and features wrappers only **choose boxes**; the heavy work is
always “load spectra → aggregate intensity in each box.”

------------------------------------------------------------------------

## 2. How xcms works

### 2.1 Layering

On current xcms (`XcmsExperiment` path):

``` text
featureChromatograms / chromPeakChromatograms
        │
        ▼
  build mz–RT boxes  (featureArea or each chromPeak’s window)
        │
        ▼
  chromatogram()  →  .mse_chromatogram()
        │
        ▼
  .chromatograms_for_peaks()   ← pure-R kernel
        │
        ▼
  (feature path) attach chromPeaks into each XChromatogram cell
```

`chromatogram()` is the primitive. The other two only decide **which**
boxes to extract and **how** to arrange the result.

### 2.2 What each function returns

**`chromatogram(object, mz, rt)`**

- Input: user-supplied mz / rt matrices (one box per row).  
- Output: regions × samples chromatograms.

**`featureChromatograms(object, features, …)`**

- Needs correspondence (`featureDefinitions`).  
- One **shared** mz–RT box per feature (default: min/max of that
  feature’s peaks via `featureArea`).  
- Same box is applied to **every** selected sample → rows = features,
  columns = samples.  
- Important: the EIC in a sample is **not** exactly that sample’s peak;
  it is the feature-global box.  
- Then attaches feature-linked `chromPeaks` (needed later for
  `removeIntensity(..., "outside_chromPeak")`).

**`chromPeakChromatograms(object, peaks, …)`**

- Needs peak picking only.  
- Each peak uses **its own** `mzmin/mzmax/rtmin/rtmax`, extracted only
  from the **owning sample**.  
- Layout: peaks × 1 column (exact peak shape).

### 2.3 The slow kernel

Per sample (inside `.chromatograms_for_peaks`):

``` text
for each box i:                              # feature-outer
  keep <- scans whose RT is in box i         # ~all scans if full RT
  for each kept scan:
    intensity <- max/sum(peaks in mz window)
  → build Chromatogram S4 object
```

Cost ≈ **`N_boxes × N_scans × cost(mz_filter)`** in interpreted R, plus
S4 construction per cell.

`featureChromatograms` adds a second nested loop over features × samples
to slot `chromPeaks` into each cell.

Sample parallelism (`chunkSize` / `BPPARAM`) only helps **across
files**. It does **not** fix the per-file feature-outer kernel. With
default `chunkSize = 2`, many workers sit idle.

``` mermaid
flowchart TD
  boxes["N mz-rt boxes"]
  chunk["Chunk samples chunkSize default 2"]
  load["Load peaksData for chunk"]
  kernel[".chromatograms_for_peaks"]
  featLoop["for each box"]
  scanLoop["for each scan in RT window: mz filter + aggregate"]
  s4["Chromatogram S4 per cell"]
  attach["feature path: Nf x Ns peak attach"]
  boxes --> chunk --> load --> kernel
  kernel --> featLoop --> scanLoop --> s4 --> attach
```

------------------------------------------------------------------------

## 3. How MSdev works

### 3.1 Triad API

``` r

# arbitrary boxes (engine)
get_xcms_chromatogram(object, mz, rt, aggregationFun = "max", BPPARAM = ...)

# peaks (chromPeakChromatograms analogue)
get_xcms_peaks_chromatogram(
  xcms, peaks.id, all.sample = FALSE,
  rt.range = c("expand", "identity", "all"), expandRt = 15
)

# features (featureChromatograms analogue)
get_xcms_feature_chromatogram(
  xcms, feature.id = NULL, sample = c("maxo", "all", ...),
  rt = c("expand", "identity", "all"),
  attachPeaks = TRUE, BPPARAM = ...
)
```

Old names `get_xcms_peaks_chrom` / `get_xcms_feature_chrom` were
**removed** (no aliases).

### 3.2 Shared engine (`get_xcms_chromatogram`)

``` mermaid
flowchart TD
  boxes["N mz/rt boxes"]
  split["One task per sample file"]
  par["bplapply BPPARAM"]
  once["Load MS1 peaksData + rtime once"]
  mat["Fill intensity matrix scans x boxes"]
  wrap["Build Chromatogram objects once"]
  cbind["Column-bind across files"]
  boxes --> split --> par --> once --> mat --> wrap --> cbind
```

Per file:

``` text
pd[[s]] = peaks for scan s;  rt[s] = retention time

# numeric matrix first (no S4 yet)
I[n_scans, N_boxes] <- max/sum intensity in each mz window
  (skip per-box RT walk when all boxes use full RT)

# then wrap I[, j] into Chromatogram objects
# combine files into MChromatograms / XChromatograms
```

Wrappers only prepare boxes:

- **Peaks:** one box per `chromPeaks` row; default extract in owning
  sample only.  
- **Features:** one shared box per feature (`peakMz*` / `peakRt*` when
  present, else min/max of `peakidx`); optional bulk chromPeak attach.

### 3.3 Project cache

`MSdev_get_feature_chrom()` calls:

``` r

get_xcms_feature_chromatogram(
  xcms.xcms,
  feature.id = fid,
  sample = "all",
  rt = "all",          # full-run EICs when needed
  attachPeaks = TRUE,
  BPPARAM = BPPARAM
)
```

and stores `Positive_Chromatograms` / `Negative_Chromatograms` on disk.
`MSdev_group_feature_EIC()` reuses those chromatograms instead of
re-calling `featureChromatograms` per group.

------------------------------------------------------------------------

## 4. How the extractor speeds up

Same disk I/O (mzML/CDF must be read). Gains are on the **post-load
kernel** and **object assembly**.

| Bottleneck in xcms | MSdev change |
|----|----|
| Feature-outer rewalk of scans | One peaks load per file; matrix fill |
| R `vapply` + `between` per box | Tight numeric loops; full-RT fast path |
| `Chromatogram` S4 inside inner loop | Defer S4 until matrix is done |
| `chunkSize = 2` underuses workers | One task = one full file |
| Nested `Nf × Ns` peak-attach S4 | Bulk assign (or skip if traces only) |
| Full `.mse_chromatogram` machinery | Thin path: Spectra → matrix → chromatograms |

### Expected order-of-magnitude (wall time, ours / xcms)

| Scenario | Typical ratio |
|----|----|
| Few boxes (≤50), narrow RT | ~0.8–1.0× (I/O-dominated) |
| Many features, peak-window RT | ~0.3–0.5× (~2–3× faster) |
| Many features, **full RT** | ~0.05–0.2× on CPU; **~3–10×** wall if I/O is ~half the job |
| Peak-level (`chromPeakChromatograms` scale) | ~0.5–0.8× |

Rough split of xcms full-RT feature extract time:

``` text
~40–70%  disk / Spectra peaksData
~20–50%  .chromatograms_for_peaks R kernel
~5–20%   Chromatogram S4 + peak attach
```

Even a much faster kernel cannot beat disk. Gains grow with more
features, wider/full RT windows, and warmer Spectra caches.

**Target for `MSdev_get_feature_chrom` (full RT, both polarities):**
wall time ≤ about **1/3** of
[`xcms::featureChromatograms`](https://rdrr.io/pkg/xcms/man/featureChromatograms.html),
with max-aggregated traces matching within numerical noise.

### Quick local timing sketch

``` r

register(SerialParam())
fids <- 1:50
pids <- 1:50

system.time(get_xcms_feature_chromatogram(xcms, fids, sample = "all", rt = "expand"))
system.time(xcms::featureChromatograms(xcms, features = fids, expandRt = 15))

system.time(get_xcms_peaks_chromatogram(xcms, pids, rt.range = "expand"))
system.time(xcms::chromPeakChromatograms(
  xcms, peaks = rownames(chromPeaks(xcms))[pids], expandRt = 15
))
```

------------------------------------------------------------------------

## 5. Choosing the right function

| Need | Prefer |
|----|----|
| Custom mz–RT regions | `get_xcms_chromatogram` |
| True peak shape in detecting sample | `get_xcms_peaks_chromatogram` |
| Feature overlays / EIC grouping across samples | `get_xcms_feature_chromatogram` |
| Full project cache for both polarities | `MSdev_get_feature_chrom` |

For EIC grouping with `onlyPeak = TRUE`, feature-level EICs still use a
**shared** box; peak attach +
`removeIntensity(..., "outside_chromPeak")` restricts correlation to
peak windows.

------------------------------------------------------------------------

## 6. Related reading

- xcms:
  [`featureChromatograms`](https://sneumann.github.io/xcms/reference/featureChromatograms.html),
  [`chromPeakChromatograms`](https://sneumann.github.io/xcms/reference/chromPeakChromatograms.html)
- MSdev article: [Feature grouping with
  `EicSimilarityParam`](https://drruili.github.io/MSdev/articles/xcms-feature-group-EicSimilarityParam.md)
- Code: `R/dev_xcms.R` (`get_xcms_chromatogram`,
  `get_xcms_peaks_chromatogram`, `get_xcms_feature_chromatogram`);
  `R/MSdev-function.R` (`MSdev_get_feature_chrom`)
