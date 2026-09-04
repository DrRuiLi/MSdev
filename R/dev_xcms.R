
#' @title  get_features_from_xcms
#' @description extract feature data from xcms::XCMSnExp,
#'  calculate RSD of QC and Sample
#'  ( note this rely on character "QC" and "Sample" in `sampleNames(xcms.xcms)` )
#' @param xcms.xcms XCMSnExp object
#'
#' @return xcms a SummarizedExperiment subject
#' @export
#'

get_features_from_xcms <- function(xcms.xcms,missing = NA){

  xcms.sum <- xcms::quantify(xcms.xcms,missing = missing )
  feature.def <- SummarizedExperiment::rowData(xcms.sum)%>%
    tibble::as_tibble()

  feature.matrix <- SummarizedExperiment::assay(xcms.sum)
  rsd <- function(x){sd(x,na.rm =  T)/mean(x , na.rm = T)}
  feature.matrix.qc <- feature.matrix[,which(grepl("QC",colnames(feature.matrix)))]
  feature.matrix.sample <- feature.matrix[,which(grepl("Sample",colnames(feature.matrix)))]
  if(sum(grepl("QC",colnames(feature.matrix)))>1){

    feature.def$qc_rsd <- apply(feature.matrix.qc, 1, rsd)
  }else(

    feature.def$qc_rsd <- 0
  )

  if(sum(grepl("Sample",colnames(feature.matrix)))>1){

    feature.def$sample_rsd <- apply(feature.matrix.sample, 1, rsd)
  }else(

    feature.def$sample_rsd <- 0
  )
  feature.def$med_intensity <- apply(feature.matrix , 1 ,median,na.rm =T)
  SummarizedExperiment::rowData(xcms.sum) <-feature.def
  return(xcms.sum)
}


#' @description Xcms feature se.
#' @describeIn xcms_extension_feature extract feature data from xcms, convert to SummarizedExperiment
#' @param xcms.xcms xcms object
#' @param missing how missing values should be reported. Allowed values are NA (the default), a numeric or missing = "rowmin_half". The latter replaces any NA with half of the row's minimal (non-missing) value.
#'
#' @returns SummarizedExperiment
#' @export
#'
get_xcms_feature_se <- function(xcms.xcms,...){

  pol <- c("0" = "neg","1" = "pos")

  xcms.xcms <- xcms_get_feature_stat(xcms.xcms)
  sample.info <- Biobase::pData(xcms.xcms)
  rownames(sample.info) <- sample.info$sample.name

  featuredef <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::mutate(xcms_feature_id = feature_id,
                  feature_id = paste0(feature_id ,"_",pol[as.character(polarity)] ))
  rownames(featuredef) <- featuredef$feature_id

  featureval <- xcms::featureValues(xcms.xcms,...)
  colnames(featureval) <- Biobase::pData(xcms.xcms)$sample.name
  rownames(featureval) <- featuredef$feature_id

  feature.se <- SummarizedExperiment::SummarizedExperiment(assays = featureval,
                       rowData = featuredef,
                       colData =sample.info
                      )
  return(feature.se)


}


get_chrom_peaks_shape_score <- function(chrom,
                                        peak.id = chrom@chromPeakData@rownames){
  peak.id = chrom@chromPeakData@rownames
  peak.id <- peak.id[1]
  peaks.data <- xcms::chromPeaks(chrom)[peak.id,,drop = FALSE]
  rtime <- xcms::rtime(chrom)
  int <- xcms::intensity(chrom)

  rtime <- rtime[!is.na(int)]
  int <- int[!is.na(int)]

  int.fit <-peak.gasssian.fit(rtime,
                    peak.apex.intensity = peaks.data[1,"maxo"],
                    peak.apex.rt = peaks.data[1,"rt"],
                    peak.half.width = min(peaks.data[1,"rtmax"]-peaks.data[1,"rt"],
                                          peaks.data[1,"rt"]-peaks.data[1,"rtmin"])/2)
  int[is.na(int)] <- 0


  #sum(abs(int-int.fit)/sum(int.fit))
  #cor(int,int.fit)
  #sqrt(mean((int-int.fit)^2))/mean(int.fit)
  r2 <- 1-sum((int-int.fit)^2)/sum( (int-mean(int))^2 )
  r2.adj <- 1-(1-r2)*(length(int)-1)/(length(int)-2)
  r2.adj
}


get_xchrom_peak_score <- function(xchrom){

  peak.info <- xcms::chromPeaks(xchrom)
  peaks.no <- nrow(peak.info)


  for (i in 1:peaks.no) {


  }


}


sub_xchrom <-  function (x, i, j, drop = TRUE)
{
  if (missing(i) && missing(j))
    return(x)
  if (missing(i))
    i <- seq_len(nrow(x))
  if (missing(j))
    j <- seq_len(ncol(x))
  if (is.logical(i))
    i <- which(i)
  if (is.logical(j))
    j <- which(j)
  if (length(i) > 1 || length(j) > 1)
    drop <- FALSE
  if (length(i) == 1 && length(j) == 1 && drop)
    return(x@.Data[i, j, drop = TRUE][[1]])
  cpeaks_orig <- xcms::chromPeaks(x)
  fts_orig <- xcms::featureDefinitions(x)
  ph <- x@.processHistory
  pd <- x@phenoData
  fd <- x@featureData
  xclass <- class(x)
  x <- as(x@.Data[i = i, j = j, drop = FALSE], xclass)
  pd <- pd[j, ]
  Biobase::pData(pd) <- droplevels(pData(pd))
  x@phenoData <- pd
  fd <- fd[i, ]
  Biobase::pData(fd) <- droplevels(pData(fd))
  x@featureData <- fd
  if (nrow(fts_orig)) {
  cpeaks_sub <- xcms::chromPeaks(x)
    fts <- vector("list", length(i))
    for (el in seq_along(i)) {
      fts_row <- fts_orig[fts_orig$row == i[el], ,
                          drop = FALSE]
      if (nrow(fts_row)) {
        fts_row$row <- el
        fts_row <- xcms:::.subset_features_on_chrom_peaks(fts_row,
                                                          cpeaks_orig, cpeaks_sub)
        fts[[el]] <- fts_row
      }
      else fts[[el]] <- DataFrame()
    }
    x@featureDefinitions <- do.call(rbind, fts)
  }

  x
}

get_xchroms_peaks_count <- function(xchroms){

  peaks.info <- xcms::chromPeaks(xchroms)%>%
    as.data.frame()%>%
    dplyr::group_by(row,column)%>%
    dplyr::mutate(peaks.no = n(),
                  idx = (row -1)*dim(xchroms)[2]+column )%>%
    dplyr::distinct(row,column ,peaks.no,idx)
  peaks.count.matrix <-rep(0,length(xchroms))
  peaks.count.matrix[peaks.info$idx] <- peaks.info$peaks.no
  peaks.count.matrix <- matrix(peaks.count.matrix,
                               nrow = dim(xchroms)[1],
                               byrow = T)%>%
    `colnames<-`(1: dim(xchroms)[2])
  return(peaks.count.matrix)
}

## ---------------------------------------------------------------------------
## Fast chromatogram triad (engine + peaks + features)
## ---------------------------------------------------------------------------

#' Normalize mz/rt inputs to two-column numeric matrices with nrow == n
#' @noRd
.normalize_mz_rt_boxes <- function(mz, rt, n = NULL) {
  if (is.null(dim(mz))) {
    mz <- matrix(as.numeric(mz), ncol = 2L, byrow = TRUE)
  } else {
    mz <- as.matrix(mz)[, 1:2, drop = FALSE]
  }
  storage.mode(mz) <- "numeric"
  if (is.null(dim(rt))) {
    rt <- matrix(as.numeric(rt), ncol = 2L, byrow = TRUE)
  } else {
    rt <- as.matrix(rt)[, 1:2, drop = FALSE]
  }
  storage.mode(rt) <- "numeric"
  if (is.null(n)) {
    n <- max(nrow(mz), nrow(rt))
  }
  if (nrow(mz) == 1L && n > 1L) {
    mz <- mz[rep(1L, n), , drop = FALSE]
  }
  if (nrow(rt) == 1L && n > 1L) {
    rt <- rt[rep(1L, n), , drop = FALSE]
  }
  if (nrow(mz) != n || nrow(rt) != n) {
    stop("'mz' and 'rt' must have the same number of rows (or be a single range)")
  }
  colnames(mz) <- c("mzmin", "mzmax")
  colnames(rt) <- c("rtmin", "rtmax")
  list(mz = mz, rt = rt)
}

#' Assign featureDefinitions compatibly for XCMSnExp / XcmsExperiment
#' @noRd
.xcms_featureDefinitions_replace <- function(object, value) {
  value <- as.data.frame(value, stringsAsFactors = FALSE)
  if (inherits(object, "XcmsExperiment") || inherits(object, "MsExperiment")) {
    xcms::featureDefinitions(object) <- value
  } else {
    xcms::featureDefinitions(object) <- S4Vectors::DataFrame(value)
  }
  object
}

#' Number of sample files in an xcms / MsExperiment object
#' @noRd
.xcms_nfiles <- function(object) {
  if (inherits(object, "MsExperiment") || inherits(object, "XcmsExperiment")) {
    length(object)
  } else {
    length(MSnbase::fileNames(object))
  }
}

#' Subset to one sample/file
#' @noRd
.xcms_filter_file <- function(object, i) {
  if (inherits(object, "MsExperiment") || inherits(object, "XcmsExperiment")) {
    object[as.integer(i)]
  } else {
    MSnbase::filterFile(object, as.integer(i))
  }
}

#' Basename sample labels
#' @noRd
.xcms_sample_names <- function(object) {
  if (inherits(object, "MsExperiment") || inherits(object, "XcmsExperiment")) {
    fn <- tryCatch(xcms::fileNames(object), error = function(e) NULL)
    if (is.null(fn) || !length(fn)) {
      fn <- tryCatch(MsExperiment::sampleData(object)$dataOrigin, error = function(e) NULL)
    }
    if (is.null(fn) || !length(fn)) {
      return(paste0("sample", seq_len(.xcms_nfiles(object))))
    }
    basename(as.character(fn))
  } else {
    basename(MSnbase::fileNames(object))
  }
}

#' Load MS1 peaks + rtime once for a single-file object
#' @noRd
.get_ms1_peaks_rtime_one_file <- function(object_one_file) {
  if (inherits(object_one_file, "XcmsExperiment") ||
      inherits(object_one_file, "MsExperiment")) {
    sp <- ProtGenerics::spectra(object_one_file)
    sp <- Spectra::filterMsLevel(sp, 1L)
    list(
      rt = as.numeric(Spectra::rtime(sp)),
      pd = as.list(Spectra::peaksData(
        sp,
        columns = c("mz", "intensity"),
        f = factor(),
        return.type = "list",
        BPPARAM = BiocParallel::SerialParam()
      ))
    )
  } else {
    fn <- MSnbase::fileNames(object_one_file)
    if (length(fn) != 1L) {
      stop("Expected a single-file xcms object")
    }
    obj_ms1 <- tryCatch(
      MSnbase::filterMsLevel(object_one_file, 1L),
      error = function(e) object_one_file
    )
    rt_obj <- as.numeric(MSnbase::rtime(obj_ms1))
    sp <- Spectra::Spectra(fn, backend = Spectra::MsBackendMzR())
    sp <- Spectra::filterMsLevel(sp, 1L)
    pd <- as.list(Spectra::peaksData(
      sp,
      columns = c("mz", "intensity"),
      f = factor(),
      return.type = "list",
      BPPARAM = BiocParallel::SerialParam()
    ))
    rt_sp <- as.numeric(Spectra::rtime(sp))
    if (length(rt_obj) == length(pd)) {
      list(rt = rt_obj, pd = pd)
    } else {
      list(rt = rt_sp, pd = pd)
    }
  }
}

#' Fill intensity matrix \[n_scans x n_boxes\] (numeric; S4 deferred)
#' @noRd
.eic_intensity_matrix <- function(pd, rt, mzmin, mzmax, rtmin, rtmax,
                                  aggregationFun = "max") {
  ns <- length(pd)
  nb <- length(mzmin)
  I <- matrix(NA_real_, nrow = ns, ncol = nb)
  use_max <- identical(aggregationFun, "max")
  rt_lo <- min(rt, na.rm = TRUE)
  rt_hi <- max(rt, na.rm = TRUE)
  full_rt <- all(rtmin <= rt_lo + 1e-6) && all(rtmax >= rt_hi - 1e-6)

  if (full_rt) {
    for (s in seq_len(ns)) {
      p <- pd[[s]]
      if (is.null(p) || !nrow(p)) {
        next
      }
      mzs <- p[, 1L]
      ints <- p[, 2L]
      for (j in seq_len(nb)) {
        sel <- mzs >= mzmin[j] & mzs <= mzmax[j]
        if (any(sel)) {
          I[s, j] <- if (use_max) max(ints[sel]) else sum(ints[sel])
        }
      }
    }
  } else {
    for (j in seq_len(nb)) {
      keep <- which(rt >= rtmin[j] & rt <= rtmax[j])
      if (!length(keep)) {
        next
      }
      for (s in keep) {
        p <- pd[[s]]
        if (is.null(p) || !nrow(p)) {
          next
        }
        sel <- p[, 1L] >= mzmin[j] & p[, 1L] <= mzmax[j]
        if (any(sel)) {
          I[s, j] <- if (use_max) max(p[sel, 2L]) else sum(p[sel, 2L])
        }
      }
    }
  }
  I
}

#' Build one-file MChromatograms from intensity matrix
#' @noRd
.chromatograms_one_file_from_matrix <- function(I, rt, mz, rt_box,
                                                aggregationFun,
                                                fromFile,
                                                sample_name,
                                                row_names = NULL) {
  nb <- nrow(mz)
  chroms <- vector("list", nb)
  for (j in seq_len(nb)) {
    keep <- which(rt >= rt_box[j, 1L] & rt <= rt_box[j, 2L])
    chroms[[j]] <- MSnbase::Chromatogram(
      rtime = rt[keep],
      intensity = if (length(keep)) I[keep, j] else numeric(),
      mz = mz[j, ],
      filterMz = mz[j, ],
      fromFile = as.integer(fromFile),
      aggregationFun = aggregationFun,
      msLevel = 1L
    )
  }
  mat <- matrix(chroms, nrow = nb, ncol = 1L)
  if (!is.null(row_names)) {
    rownames(mat) <- row_names
  }
  colnames(mat) <- sample_name
  fd <- Biobase::AnnotatedDataFrame(data.frame(
    mzmin = mz[, 1L],
    mzmax = mz[, 2L],
    rtmin = rt_box[, 1L],
    rtmax = rt_box[, 2L],
    row.names = if (is.null(row_names)) seq_len(nb) else row_names,
    stringsAsFactors = FALSE
  ))
  phd <- Biobase::AnnotatedDataFrame(data.frame(
    sampleName = sample_name,
    row.names = sample_name,
    stringsAsFactors = FALSE
  ))
  ## MChromatograms() does matrix(data, ...); pass ncol to avoid flattening
  MSnbase::MChromatograms(
    mat,
    ncol = 1L,
    phenoData = phd,
    featureData = fd
  )
}

#' Combine per-file chromatogram objects by column
#' @noRd
.combine_chromatograms_files <- function(chrom_list) {
  chrom_list <- chrom_list[!vapply(chrom_list, is.null, logical(1))]
  if (!length(chrom_list)) {
    return(NULL)
  }
  if (length(chrom_list) == 1L) {
    return(chrom_list[[1L]])
  }
  ## Avoid base::cbind (drops MChromatograms to a plain matrix)
  mats <- lapply(chrom_list, function(x) x@.Data)
  nr <- vapply(mats, nrow, integer(1))
  if (length(unique(nr)) != 1L) {
    stop("Per-file chromatograms have inconsistent row counts")
  }
  mat <- do.call(cbind, mats)
  fd <- Biobase::featureData(chrom_list[[1L]])
  pds <- lapply(chrom_list, function(x) {
    pd <- Biobase::pData(x)
    ## ensure unique rownames = colnames of each block
    rn <- colnames(x)
    if (!is.null(rn) && length(rn) == nrow(pd)) {
      rownames(pd) <- rn
    }
    pd
  })
  pd <- do.call(rbind, pds)
  phd <- Biobase::AnnotatedDataFrame(pd)
  colnames(mat) <- rownames(pd)
  ## MChromatograms() re-runs matrix(data, ...); must pass ncol
  MSnbase::MChromatograms(
    mat,
    ncol = ncol(mat),
    phenoData = phd,
    featureData = fd
  )
}

#' @describeIn xcms_extension_chromatogram extract chromatograms
#' @description Fast per-file EIC extractor. Loads MS1 peaks once per sample,
#'   fills an intensity matrix for all mz-rt boxes, then wraps
#'   \code{MSnbase::Chromatogram} objects. Drop-in style replacement for
#'   \code{xcms::chromatogram()} for rectangular mz/rt region extraction.
#' @param object XCMSnExp / XcmsExperiment (or single-file subset).
#' @param mz numeric matrix with columns mzmin, mzmax (one row per EIC).
#' @param rt numeric matrix with columns rtmin, rtmax (one row per EIC), or a
#'   single range recycled to all rows of \code{mz}.
#' @param aggregationFun character(1). Intensity aggregation within the m/z
#'   window per scan: \code{"max"} (default) or \code{"sum"}.
#' @param BPPARAM BiocParallel backend; one task per sample file.
#' @param msLevel integer; kept for API compatibility (MS1 extraction).
#' @param ... ignored (compatibility with older callers).
#'
#' @return \code{MChromatograms} with rows = regions, columns = samples.
#' @export
get_xcms_chromatogram <- function(object,
                                  mz,
                                  rt,
                                  aggregationFun = "max",
                                  BPPARAM = SerialParam(),
                                  msLevel = 1L,
                                  ...) {
  aggregationFun <- match.arg(aggregationFun, c("max", "sum"))
  boxes <- .normalize_mz_rt_boxes(mz, rt)
  mz <- boxes$mz
  rt_box <- boxes$rt
  nfiles <- .xcms_nfiles(object)
  if (!nfiles) {
    stop("'object' has no sample files")
  }
  snames <- .xcms_sample_names(object)
  row_names <- rownames(mz)
  if (is.null(row_names)) {
    row_names <- rownames(rt_box)
  }

  chrom_list <- bplapply(seq_len(nfiles), function(i) {
    register(SerialParam())
    one <- .xcms_filter_file(object, i)
    peaks <- .get_ms1_peaks_rtime_one_file(one)
    I <- .eic_intensity_matrix(
      peaks$pd, peaks$rt,
      mz[, 1L], mz[, 2L],
      rt_box[, 1L], rt_box[, 2L],
      aggregationFun = aggregationFun
    )
    .chromatograms_one_file_from_matrix(
      I, peaks$rt, mz, rt_box,
      aggregationFun = aggregationFun,
      fromFile = i,
      sample_name = snames[i],
      row_names = row_names
    )
  }, BPPARAM = BPPARAM)

  .combine_chromatograms_files(chrom_list)
}

#' @describeIn xcms_extension_chromatogram extract chromatogram for a peak
#' @description Extract EICs for chromatographic peaks (xcms
#'   \code{chromPeakChromatograms} analogue). Uses
#'   \code{\link{get_xcms_chromatogram}} as the engine.
#' @param xcms.xcms XCMSnExp / XcmsExperiment with chromPeaks.
#' @param peaks.id character or numeric peak IDs / indices.
#' @param selected_sample Sample selection. \code{NULL} (default) extracts each
#'   peak only in its owning sample (one-column result, chromPeakChromatograms
#'   contract). \code{"all"} extracts every peak box in all samples. Integer
#'   indices or sample name(s) extract in those samples only.
#' @param rt.range one of \code{c("all","identity","expand")}.
#' @param expandRt numeric seconds added on each side when \code{rt.range="expand"}.
#' @param aggregationFun passed to \code{get_xcms_chromatogram}.
#' @param BPPARAM BiocParallel backend.
#'
#' @return MChromatograms / column-bound chromatograms.
#' @export
get_xcms_peaks_chromatogram <- function(xcms.xcms,
                                        peaks.id,
                                        selected_sample = NULL,
                                        rt.range = c("expand", "identity", "all"),
                                        expandRt = 15,
                                        aggregationFun = "max",
                                        BPPARAM = SerialParam()) {
  rt.range <- match.arg(rt.range)
  peaks.data <- xcms::chromPeaks(xcms.xcms)
  if (is.numeric(peaks.id)) {
    peaks.id <- rownames(peaks.data)[peaks.id]
  }
  peaks.data <- peaks.data[peaks.id, , drop = FALSE]
  n <- nrow(peaks.data)
  rt_all <- range(MSnbase::rtime(xcms.xcms), na.rm = TRUE)
  rtr <- switch(
    rt.range,
    "all" = matrix(rt_all, nrow = n, ncol = 2L, byrow = TRUE),
    "expand" = t(apply(peaks.data[, c("rtmin", "rtmax"), drop = FALSE], 1L,
                       expand_range, add = expandRt)),
    "identity" = peaks.data[, c("rtmin", "rtmax"), drop = FALSE]
  )
  mzr <- peaks.data[, c("mzmin", "mzmax"), drop = FALSE]
  rownames(mzr) <- rownames(peaks.data)
  rownames(rtr) <- rownames(peaks.data)

  owning_only <- is.null(selected_sample)
  if (!owning_only) {
    if (identical(selected_sample, "all")) {
      xcms.sub <- xcms.xcms
    } else {
      sn <- .xcms_sample_names(xcms.xcms)
      pdata_names <- tryCatch(
        as.character(Biobase::pData(xcms.xcms)$sample.name),
        error = function(e) NULL
      )
      sample_idx <- .resolve_selected_sample(selected_sample, sn, pdata_names)
      xcms.sub <- .xcms_filter_file(xcms.xcms, sample_idx)
    }
    return(get_xcms_chromatogram(
      xcms.sub,
      mz = mzr,
      rt = rtr,
      aggregationFun = aggregationFun,
      BPPARAM = BPPARAM
    ))
  }

  ## chromPeakChromatograms contract: peaks x 1 (owning sample only)
  sample_ids <- as.integer(peaks.data[, "sample"])
  by_sample <- split(seq_len(n), sample_ids)
  chrom_parts <- vector("list", length(by_sample))
  names(chrom_parts) <- names(by_sample)
  for (si in names(by_sample)) {
    idx <- by_sample[[si]]
    one <- .xcms_filter_file(xcms.xcms, as.integer(si))
    ch <- get_xcms_chromatogram(
      one,
      mz = mzr[idx, , drop = FALSE],
      rt = rtr[idx, , drop = FALSE],
      aggregationFun = aggregationFun,
      BPPARAM = SerialParam()
    )
    ## force single logical column name; keep peak rownames
    colnames(ch) <- "1"
    chrom_parts[[si]] <- ch
  }
  ## stack rows in original peaks.id order
  out_list <- vector("list", n)
  for (si in names(by_sample)) {
    idx <- by_sample[[si]]
    ch <- chrom_parts[[si]]
    for (k in seq_along(idx)) {
      out_list[[idx[k]]] <- ch[k, 1, drop = TRUE]
    }
  }
  mat <- matrix(out_list, nrow = n, ncol = 1L)
  rownames(mat) <- rownames(peaks.data)
  colnames(mat) <- "1"
  fd <- Biobase::AnnotatedDataFrame(data.frame(
    mzmin = mzr[, 1L],
    mzmax = mzr[, 2L],
    rtmin = rtr[, 1L],
    rtmax = rtr[, 2L],
    sample_index = sample_ids,
    row.names = rownames(peaks.data),
    stringsAsFactors = FALSE
  ))
  phd <- Biobase::AnnotatedDataFrame(data.frame(
    sampleName = "1",
    row.names = "1",
    stringsAsFactors = FALSE
  ))
  res <- MSnbase::MChromatograms(
    mat,
    ncol = 1L,
    phenoData = phd,
    featureData = fd
  )
  ## attach chromPeaks when possible
  if (methods::is(res, "MChromatograms")) {
    res <- as(res, "XChromatograms")
    pks <- peaks.data
    pkd <- tryCatch(
      as(xcms::chromPeakData(xcms.xcms)[rownames(peaks.data), , drop = FALSE], "DataFrame"),
      error = function(e) S4Vectors::DataFrame(row.names = rownames(peaks.data))
    )
    tmp <- res@.Data
    for (i in seq_len(n)) {
      slot(tmp[i, 1L][[1L]], "chromPeaks", check = FALSE) <- pks[i, , drop = FALSE]
      slot(tmp[i, 1L][[1L]], "chromPeakData", check = FALSE) <- pkd[i, , drop = FALSE]
    }
    slot(res, ".Data", check = FALSE) <- tmp
  }
  res
}

#' Build feature mz/rt area matrix
#' @noRd
.feature_mz_rt_boxes <- function(xcms.xcms, features.data, rt = c("expand", "identity", "all"),
                                 expandRt = 15, mz.expand = 0) {
  rt <- match.arg(rt)
  n <- nrow(features.data)
  has_peak_rt <- all(c("peakRtMin", "peakRtMax") %in% colnames(features.data))
  has_peak_mz <- all(c("peakMzMin", "peakMzMax") %in% colnames(features.data))

  if (!has_peak_mz || !has_peak_rt) {
    area <- tryCatch(
      as.matrix(xcms::featureArea(
        xcms.xcms,
        features = rownames(features.data)
      )),
      error = function(e) NULL
    )
    if (is.null(area)) {
      ## Fallback for older XCMSnExp without featureArea method
      pks <- xcms::chromPeaks(xcms.xcms)
      area <- do.call(rbind, lapply(features.data$peakidx, function(z) {
        p <- pks[z, , drop = FALSE]
        c(
          mzmin = min(p[, "mzmin"], na.rm = TRUE),
          mzmax = max(p[, "mzmax"], na.rm = TRUE),
          rtmin = min(p[, "rtmin"], na.rm = TRUE),
          rtmax = max(p[, "rtmax"], na.rm = TRUE)
        )
      }))
      rownames(area) <- rownames(features.data)
    }
    mzr <- area[, c("mzmin", "mzmax"), drop = FALSE]
    rtr_id <- area[, c("rtmin", "rtmax"), drop = FALSE]
  } else {
    mzr <- as.matrix(features.data[, c("peakMzMin", "peakMzMax"), drop = FALSE])
    rtr_id <- as.matrix(features.data[, c("peakRtMin", "peakRtMax"), drop = FALSE])
  }
  colnames(mzr) <- c("mzmin", "mzmax")
  colnames(rtr_id) <- c("rtmin", "rtmax")

  rt_all <- range(MSnbase::rtime(xcms.xcms), na.rm = TRUE)
  rtr <- switch(
    rt,
    "all" = matrix(rt_all, nrow = n, ncol = 2L, byrow = TRUE),
    "expand" = t(apply(rtr_id, 1L, expand_range, add = expandRt)),
    "identity" = rtr_id
  )
  colnames(rtr) <- c("rtmin", "rtmax")

  if (mz.expand > 0) {
    mz_range <- mzr[, 2L] - mzr[, 1L]
    mzr[, 1L] <- mzr[, 1L] - mz_range * mz.expand
    mzr[, 2L] <- mzr[, 2L] + mz_range * mz.expand
  }
  rownames(mzr) <- rownames(features.data)
  rownames(rtr) <- rownames(features.data)
  list(mz = mzr, rt = rtr)
}

#' Attach feature chromPeaks into XChromatograms (bulk)
#' @noRd
.attach_feature_chrom_peaks <- function(chrs, xcms.xcms, feature_ids) {
  chrs <- as(chrs, "XChromatograms")
  pks <- xcms::chromPeaks(xcms.xcms)
  pkd <- tryCatch(
    as(xcms::chromPeakData(xcms.xcms), "DataFrame"),
    error = function(e) NULL
  )
  fdef <- xcms::featureDefinitions(xcms.xcms)
  if ("feature_id" %in% colnames(fdef)) {
    f_idx <- match(feature_ids, fdef$feature_id)
    if (anyNA(f_idx)) {
      f_idx <- match(feature_ids, rownames(fdef))
    }
  } else {
    f_idx <- match(feature_ids, rownames(fdef))
  }
  if (anyNA(f_idx)) {
    warning("Some features could not be mapped for chromPeak attachment")
  }
  nf <- nrow(chrs)
  ns <- ncol(chrs)
  pks_empty <- pks[integer(), , drop = FALSE]
  pkd_empty <- if (is.null(pkd)) {
    S4Vectors::DataFrame()
  } else {
    pkd[integer(), , drop = FALSE]
  }
  tmp <- chrs@.Data
  for (i in seq_len(nf)) {
    ii <- f_idx[i]
    if (is.na(ii)) {
      for (j in seq_len(ns)) {
        slot(tmp[i, j][[1L]], "chromPeaks", check = FALSE) <- pks_empty
        if (!is.null(pkd)) {
          slot(tmp[i, j][[1L]], "chromPeakData", check = FALSE) <- pkd_empty
        }
      }
      next
    }
    idx <- fdef$peakidx[[ii]]
    smpl <- as.integer(pks[idx, "sample"])
    for (j in seq_len(ns)) {
      keep <- smpl == j
      if (any(keep)) {
        slot(tmp[i, j][[1L]], "chromPeaks", check = FALSE) <-
          pks[idx[keep], , drop = FALSE]
        if (!is.null(pkd)) {
          slot(tmp[i, j][[1L]], "chromPeakData", check = FALSE) <-
            pkd[idx[keep], , drop = FALSE]
        }
      } else {
        slot(tmp[i, j][[1L]], "chromPeaks", check = FALSE) <- pks_empty
        if (!is.null(pkd)) {
          slot(tmp[i, j][[1L]], "chromPeakData", check = FALSE) <- pkd_empty
        }
      }
    }
  }
  slot(chrs, ".Data", check = FALSE) <- tmp
  chrs
}

#' @describeIn xcms_extension_chromatogram extract chromatograms for features
#' @description Extract EICs for features (xcms \code{featureChromatograms}
#'   analogue) via \code{\link{get_xcms_chromatogram}}. One shared mz-rt box
#'   per feature is applied to each selected sample.
#' @param xcms.xcms XCMSnExp / XcmsExperiment with featureDefinitions.
#' @param feature.id character/numeric feature IDs (default all).
#' @param selected_sample Sample selection. \code{"maxo"} (default) uses the
#'   sample with highest mean feature value; \code{"all"} uses all samples;
#'   integer indices or sample name(s) select those samples.
#' @param rt one of \code{c("all","expand","identity")}.
#' @param expandRt seconds added each side when \code{rt="expand"}.
#' @param mz.expand fraction of mz width to expand on each side.
#' @param aggregationFun passed to \code{get_xcms_chromatogram}.
#' @param attachPeaks logical; attach feature chromPeaks into
#'   \code{XChromatograms} (needed for \code{removeIntensity(..., "outside_chromPeak")}).
#' @param BPPARAM BiocParallel backend.
#'
#' @return XChromatograms. The parent \code{featureDefinitions} rows for the
#'   extracted features are stored in the \code{featureDefinitions} slot
#'   (access with \code{obj@featureDefinitions}), so the object is
#'   self-describing (carries \code{feature_id}, \code{mzmed}, \code{rtmed},
#'   and \code{peakRtMin}/\code{peakRtMax} when present).
#' @export
get_xcms_feature_chromatogram <- function(xcms.xcms,
                                          feature.id = NULL,
                                          selected_sample = "maxo",
                                          rt = c("expand", "identity", "all"),
                                          expandRt = 15,
                                          mz.expand = 0,
                                          aggregationFun = "max",
                                          attachPeaks = TRUE,
                                          BPPARAM = SerialParam(progressbar = TRUE)) {
  rt <- match.arg(rt)
  fdef_all <- xcms::featureDefinitions(xcms.xcms)
  if (is.null(feature.id)) {
    if ("feature_id" %in% colnames(fdef_all)) {
      feature.id <- as.character(fdef_all$feature_id)
    } else {
      feature.id <- rownames(fdef_all)
    }
  }
  if (is.numeric(feature.id)) {
    feature.id <- rownames(fdef_all)[feature.id]
  }
  feature.id <- as.character(feature.id)
  if ("feature_id" %in% colnames(fdef_all)) {
    ii <- match(feature.id, as.character(fdef_all$feature_id))
    if (anyNA(ii)) {
      ii2 <- match(feature.id, rownames(fdef_all))
      ii[is.na(ii)] <- ii2[is.na(ii)]
    }
    if (anyNA(ii)) {
      stop("Unknown feature.id: ", paste(feature.id[is.na(ii)], collapse = ", "))
    }
    features.data <- fdef_all[ii, , drop = FALSE]
  } else {
    features.data <- fdef_all[feature.id, , drop = FALSE]
  }
  rownames(features.data) <- feature.id

  if (identical(selected_sample, "maxo")) {
    features.val <- xcms::featureValues(xcms.xcms, missing = "rowmin_half")
    features.val <- features.val[rownames(features.data), , drop = FALSE]
    orig_idx <- which.max(colMeans(features.val, na.rm = TRUE))
    xcms.sub <- .xcms_filter_file(xcms.xcms, orig_idx)
  } else if (identical(selected_sample, "all") || is.null(selected_sample)) {
    orig_idx <- NULL
    xcms.sub <- xcms.xcms
  } else {
    sn <- .xcms_sample_names(xcms.xcms)
    pdata_names <- tryCatch(
      as.character(Biobase::pData(xcms.xcms)$sample.name),
      error = function(e) NULL
    )
    orig_idx <- .resolve_selected_sample(selected_sample, sn, pdata_names)
    xcms.sub <- .xcms_filter_file(xcms.xcms, orig_idx)
  }

  boxes <- .feature_mz_rt_boxes(
    xcms.xcms, features.data,
    rt = rt, expandRt = expandRt, mz.expand = mz.expand
  )
  x.chrom <- get_xcms_chromatogram(
    xcms.sub,
    mz = boxes$mz,
    rt = boxes$rt,
    aggregationFun = aggregationFun,
    BPPARAM = BPPARAM
  )
  if (isTRUE(attachPeaks)) {
    ## sample indices in subset must map to original sample numbers for peak attach
    if (identical(selected_sample, "all") || is.null(selected_sample)) {
      x.chrom <- .attach_feature_chrom_peaks(x.chrom, xcms.xcms, feature.id)
    } else {
      x.chrom <- as(x.chrom, "XChromatograms")
      ## remap: subset has 1..k columns; chromPeaks sample ids are original
      pks <- xcms::chromPeaks(xcms.xcms)
      pkd <- tryCatch(as(xcms::chromPeakData(xcms.xcms), "DataFrame"), error = function(e) NULL)
      fdef <- xcms::featureDefinitions(xcms.xcms)
      if ("feature_id" %in% colnames(fdef)) {
        f_idx <- match(feature.id, as.character(fdef$feature_id))
      } else {
        f_idx <- match(feature.id, rownames(fdef))
      }
      tmp <- x.chrom@.Data
      for (i in seq_along(feature.id)) {
        ii <- f_idx[i]
        if (is.na(ii)) next
        idx <- fdef$peakidx[[ii]]
        for (jc in seq_len(ncol(x.chrom))) {
          oj <- orig_idx[min(jc, length(orig_idx))]
          keep <- as.integer(pks[idx, "sample"]) == oj
          if (any(keep)) {
            slot(tmp[i, jc][[1L]], "chromPeaks", check = FALSE) <-
              pks[idx[keep], , drop = FALSE]
            if (!is.null(pkd)) {
              slot(tmp[i, jc][[1L]], "chromPeakData", check = FALSE) <-
                pkd[idx[keep], , drop = FALSE]
            }
          }
        }
      }
      slot(x.chrom, ".Data", check = FALSE) <- tmp
    }
  }
  ## Store parent featureDefinitions on the chroms (self-describing)
  if (!inherits(x.chrom, "XChromatograms")) {
    x.chrom <- as(x.chrom, "XChromatograms")
  }
  fdf <- S4Vectors::DataFrame(features.data[rownames(x.chrom), , drop = FALSE])
  x.chrom@featureDefinitions <- fdf
  x.chrom
}


xcms_get_peak_fill <- function(xcms.xcms){

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)
  rt.na <- apply(xcms.peaks,1,function(x){is.na(x["rt"])})
  into.na <- apply(xcms.peaks,1,function(x){is.na(x["into"])})

  ### fill rt with mean
  xcms.peaks[rt.na,"rt"] <- apply(xcms.peaks[rt.na,c("rtmax","rtmin"),drop =F],
                                  1,mean)
  ### fill into with coef maxo peak with
  xcms.peaks.stat <- get_xcms_peaks_stat(xcms.xcms)%>%
    dplyr::mutate(coef = into/maxo/peakWidth)
  into.coef <- median(xcms.peaks.stat$coef,na.rm = T)
  xcms.peaks[into.na,"into"] <- apply(xcms.peaks[into.na,,drop =F],
                                  1,function(x){x["maxo"]*(x["rtmax"]-x["rtmin"])*into.coef})


  ###return
  xcms.ph <- xcms::processHistory(xcms.xcms)
  xcms::chromPeaks(xcms.xcms) <- xcms.peaks
  xcms.xcms@.processHistory <- xcms.ph
  return(xcms.xcms)
}


#' @describeIn xcms_extension_feature_group group features
#' @description Groups features from an XCMSnExp object using multiple criteria: similarity in retention time, abundance (intensity) correlation, and EIC (extracted ion chromatogram) correlation. RT/abundance grouping uses MsFeatures; EIC similarity uses xcms::EicSimilarityParam.
#' @param xcms.xcms XCMSnExp object containing feature definitions.
#' @param diffRt numeric. Maximum allowed retention time difference for grouping by SimilarRtimeParam. If NULL, retention time grouping is skipped. Default is 5.
#' @param intCor numeric. Threshold for abundance similarity (correlation) grouping using AbundanceSimilarityParam. If NULL, intensity correlation grouping is skipped. Default is 0.5.
#' @param eicCor numeric. Threshold for EIC similarity grouping using xcms::EicSimilarityParam. If NULL, EIC correlation grouping is skipped. Default is 0.5.
#'
#' @return XCMSnExp object with featureGroups column added or updated.
#' @export
#'

xcms_get_feature_group <- function(xcms.xcms,
                                   diffRt = 5,
                                   intCor = 0.5,
                                   eicCor = 0.5){

  xcms::featureGroups(xcms.xcms) <- NA
  register(SnowParam(progressbar = T))
  if (!is.null(diffRt)) {
    message_with_time(" group by SimilarRtimeParam")
    xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                   param = MsFeatures::SimilarRtimeParam(diffRt,
                                                             groupFun = groupHclust ))
    message(length(unique(xcms::featureGroups(xcms.xcms)))," feature group")
  }
  if (!is.null(intCor)) {
    message_with_time(" group by AbundanceSimilarityParam")
    xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                    param = MsFeatures::AbundanceSimilarityParam(threshold = intCor,
                                                                     transform = log2 ),
                                    filled = TRUE)
    message(length(unique(xcms::featureGroups(xcms.xcms)))," feature group")
  }
  if (!is.null(eicCor)) {
    register(SerialParam())
    message_with_time(" group by EicSimilarityParam")
    xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                    param = xcms::EicSimilarityParam(threshold = eicCor,
                                                               n = 2))
    message(length(unique(xcms::featureGroups(xcms.xcms)))," feature group")
  }

  return(xcms.xcms)
}


get_chrom_peaks_gaussian_fit <- function(xchrom){

  peaks.info <- xcms::chromPeaks(xchrom)[1,,drop = F]
  peaks.data <- get_chroms_data(xchrom)
  nls(formula = intensity ~ gaussian_functioin(rt ,a,b,c),
      data = peaks.data,control = nls.control(warnOnly = T),
      start = list(a = peaks.info[,"maxo"],
                   b = peaks.info[,"rt"],
                   c = mean(diff(peaks.info[,c("rtmin","rtmax")])))) -> gaussian.fit

}


#' @title get_chroms_data
#' @description extract chomatogram data to a data.frame
#' @param xchrom XChromatograms
#'
#' @return xcms
#' @export
#'

get_chroms_data <- function(xchrom){

  .extract.chrom <- function(i,j){
    this.chrom <- xchrom[i,j]
    data.frame(
      rt = rtime(this.chrom),
      intensity =intensity(this.chrom),
      row =i,col = j
    )
  }
  if (class(xchrom) %in% c("XChromatogram","Chromatogram")) {
    xchrom <- XChromatograms(list(xchrom))
  }
  bp.matrix <- expand.grid(1:nrow(xchrom),1:ncol(xchrom))
  xchrom.data <- BiocParallel::bpmapply(.extract.chrom,
                                        bp.matrix[,1],bp.matrix[,2],
                         BPPARAM = BiocParallel::SerialParam(progressbar = F),SIMPLIFY=F)%>%
    do.call("rbind",.)

  return(xchrom.data)

}


#' @describeIn xcms_extension_chromatogram convert retention time units
#' @description Changes the retention time units of XChromatograms objects. In some situations (e.g., SRM data from Thermo), retention times are recorded in minutes, which can cause errors during peak detection. This function converts between seconds and minutes.
#' @param xchroms XChromatograms or MChromatograms object.
#' @param unit_to Target unit: "s" (seconds) multiplies by 60, "m" (minutes) divides by 60. Default is "s".
#' @param BPPARAM BiocParallel backend for parallel processing. Default is BatchtoolsParam.
#'
#' @return XChromatograms object with converted retention times.
#' @export
#'
XChromatograms_rt_unit <- function(xchroms,unit_to = "s",
                                   BPPARAM = BatchtoolsParam(progressbar = T,log = F,
                                                             registryargs = batchtoolsRegistryargs(packages = c("MSnbase")))){


  unit.mulit <- switch(unit_to,
            "s" = 60,
            "m" = 1/60)
  rtime.max <- max(rtime(xchroms[1,1]))
  f <- function(i,xchroms,unit.mulit){

    i.position <- which(matrix(1:length(xchroms),
                               nrow = dim(xchroms)[1],
                               byrow = T)==i,arr.ind = T)
    x <- i.position[1]
    y <- i.position[2]
    xchroms[x,y]@rtime <- xchroms[x,y]@rtime * unit.mulit
    xchroms[x,y]

  }

  pda <- pData(xchroms)
  fda <- fData(xchroms)

  xchroms.trans <- bplapply(1:length(xchroms),FUN = f,
                              xchroms = xchroms,
                              unit.mulit = unit.mulit,
                              BPPARAM = BPPARAM)
  xchroms.trans <- XChromatograms(xchroms.trans,
                                    nrow = dim(xchroms)[1],byrow = T)

  message( "rtime value ", round(rtime.max,0), " change to ", round(max(rtime(xchroms.trans[1,1])),0))

  pda -> pData(xchroms.trans)
  fda -> fData(xchroms.trans)

  return(xchroms.trans)


}


#' @describeIn xcms_extension_chromatogram subset feature chromatograms by peak RT
#' @description Crop each feature row of an \code{XChromatograms} object to its
#'   \code{peakRtMin}/\code{peakRtMax} window read from the
#'   \code{featureDefinitions} slot. Optional \code{expandRt} (seconds each
#'   side) widens the window; if the resulting width is still below
#'   \code{min_width}, both sides are padded equally until the width reaches
#'   \code{min_width}. Requires \code{peakRtMin} and \code{peakRtMax} in
#'   \code{obj@featureDefinitions} (e.g. after
#'   \code{\link{get_xcms_feature_chromatogram}} on an object that has run
#'   \code{\link{xcms_get_feature_def_stat}}). \code{NA} intensities inside the
#'   window are set to \code{0} so every feature carries a full baseline over
#'   its window; because all features share the sample's master RT grid, this
#'   keeps RTs aligned across features for pairwise comparison.
#' @param xchroms XChromatograms with feature rownames and a
#'   \code{featureDefinitions} slot containing \code{peakRtMin},
#'   \code{peakRtMax}.
#' @param expandRt numeric(1). Seconds added on each side of
#'   \code{[peakRtMin, peakRtMax]}. Default \code{0} (no expansion).
#' @param min_width numeric(1). Minimum RT window width (seconds) after
#'   \code{expandRt}. If the window is shorter, pad both sides equally.
#'   Default \code{0} (no minimum).
#'
#' @return XChromatograms with each feature's EIC cropped to its RT window and
#'   \code{NA} intensities filled with \code{0}. The \code{phenoData},
#'   \code{featureData}, and \code{featureDefinitions} slots are carried over
#'   from the input unchanged.
#' @export
#'
XChromatograms_subset_feature <- function(xchroms,
                                          expandRt = 0,
                                          min_width = 0) {
  if (!inherits(xchroms, "XChromatograms")) {
    stop("`xchroms` must be an XChromatograms object (with a featureDefinitions slot)")
  }
  if (!is.numeric(expandRt) || length(expandRt) != 1L || !is.finite(expandRt)) {
    stop("`expandRt` must be a finite numeric(1)")
  }
  if (!is.numeric(min_width) || length(min_width) != 1L ||
      !is.finite(min_width) || min_width < 0) {
    stop("`min_width` must be a non-negative finite numeric(1)")
  }

  fdef <- xchroms@featureDefinitions
  if (!all(c("peakRtMin", "peakRtMax") %in% colnames(fdef))) {
    stop(
      "featureDefinitions slot must contain 'peakRtMin' and 'peakRtMax' ",
      "(run get_xcms_feature_chromatogram after xcms_get_feature_def_stat)"
    )
  }
  if (nrow(fdef) != nrow(xchroms)) {
    stop("featureDefinitions rows (", nrow(fdef),
         ") do not match chromatogram rows (", nrow(xchroms), ")")
  }

  nr <- nrow(xchroms)
  nc <- ncol(xchroms)
  peak_rt_min <- as.numeric(fdef$peakRtMin)
  peak_rt_max <- as.numeric(fdef$peakRtMax)
  chrom_list <- vector("list", nr * nc)
  k <- 0L
  for (i in seq_len(nr)) {
    rmin <- peak_rt_min[i] - expandRt
    rmax <- peak_rt_max[i] + expandRt
    if (!is.finite(rmin) || !is.finite(rmax)) {
      stop("Non-finite peakRtMin/peakRtMax for row ", i,
           " (", rownames(fdef)[i], ")")
    }
    w <- rmax - rmin
    if (is.finite(w) && w < min_width) {
      pad <- (min_width - w) / 2
      rmin <- rmin - pad
      rmax <- rmax + pad
    }
    for (j in seq_len(nc)) {
      k <- k + 1L
      ch <- MSnbase::filterRt(xchroms[i, j], rt = c(rmin, rmax))
      ## Extraction keeps every scan RT in the window (NA where no signal);
      ## fill NA intensities with 0 so a baseline exists across the whole
      ## window. All features share the master RT grid, so filling per feature
      ## keeps RTs aligned across features for pairwise comparison.
      inten <- MSnbase::intensity(ch)
      inten[!is.finite(inten)] <- 0
      ch@intensity <- inten
      chrom_list[[k]] <- ch
    }
  }

  out <- xcms::XChromatograms(chrom_list, nrow = nr, byrow = TRUE)
  ## Assign metadata via slots: pData<-/fData<- each validate their rownames
  ## against the other slot, so neither can be set first (chicken-and-egg).
  ## Slot assignment skips validObject; phenoData/featureData carried from the
  ## input are already mutually consistent.
  out@phenoData <- xchroms@phenoData
  out@featureData <- xchroms@featureData
  out@featureDefinitions <- fdef
  dimnames(out@.Data) <- list(rownames(xchroms), colnames(xchroms))
  out
}


#' @describeIn xcms_extension_chromatogram fill chromatograms with fewer than two data points
#' @description When using xcms::findChromPeaks, chromatograms with fewer than two data points cause errors. This function identifies such chromatograms and adds a duplicate point (time +1, intensity 0) to ensure at least two points exist.
#' @param xchroms XChromatograms object to be checked and filled.
#'
#' @return XChromatograms object with chromatograms having at least two data points.
#' @export
#'

XChromatograms_fill_2point <- function(xchroms){

  rt.point <- sapply(1:length(xchroms), function(x){
    x.position <- arrayInd(x,.dim = dim(xchroms))
    length(rtime(xcms.xcms[x.position[1],x.position[2]]))
  })

  tofill <- which(rt.point <2)%>%
    arrayInd(.dim = dim(xchroms))

  f <- function(chrom){
    if (length(rtime(chrom)) ==1) {
      chrom@rtime <- c(chrom@rtime , chrom@rtime +1)
      chrom@intensity <- c(chrom@intensity ,0)

    }
    chrom
  }
  xchroms.filled <- bplapply(unlist(xchroms),f ,BPPARAM = BatchtoolsParam(progressbar = T,
                                                             registryargs = batchtoolsRegistryargs(packages = c("MSnbase"))))
  xchroms.filled <- XChromatograms(xchroms.filled,
                                  nrow = dim(xchroms)[1])
  return(xchroms.filled)

}


#' @describeIn xcms_extension_chromatogram plot chromatograms
#' @description Plots XChromatograms data as line plots, with options to normalize intensities to 0-1 range, offset chromatograms for clarity, and customize colors. Returns a ggplot object.
#' @param xchroms XChromatograms object to plot.
#' @param norm logical. If TRUE, normalize intensities to 0-1 range (default TRUE).
#' @param move logical. If TRUE, offset chromatograms by index for better visibility (default TRUE).
#' @param color_by Character indicating grouping for coloring: "column" (by sample) or "row" (by feature). Default is "column".
#' @param color_f Optional character vector of colors for groups. If NULL, uses distinctColorPalette.
#' @param label_df Optional data frame with columns x, y, label for adding text labels via ggrepel.
#'
#' @return ggplot object.
#' @export
#'

plot_XChromatograms <- function(xchroms ,
                                norm = T,
                                move = T,
                                color_by = c("column","row"),
                                color_f = NULL,
                                label_df = NULL){

  color_by = match.arg(color_by)

  if (norm) {
    xchroms <- MSnbase::normalise(xchroms)
    chrom.data <- get_chroms_data(xchroms)%>%
      dplyr::mutate(intensity = intensity*100)
  }else{
    chrom.data <- get_chroms_data(xchroms)
  }


  chrom.data <- chrom.data%>%
    dplyr::mutate(group_idx = case_when(color_by == "column"~col,
                                      color_by == "row"~row,
                                      T~col))
  if (is.null(color_f)) {
    color_f <-paste0("Peaks_",num2str( unique(chrom.data$group_idx)))
  }
  chrom.data <- chrom.data%>%
    dplyr::mutate(peaks.origin = paste0("peak_",num2str(row),"_",num2str(col)),
                  peaks.origin = factor(peaks.origin,level = unique(peaks.origin)))%>%
    dplyr::group_by(peaks.origin)%>%
    dplyr::mutate(peaks.idx =cur_group_id(),
                  color = color_f[group_idx]
    )%>%
    dplyr::ungroup()



    if (move) {
      chrom.data <- chrom.data%>%
        dplyr::mutate(rt = rt +peaks.idx*3,
                      intensity = intensity+peaks.idx*3)
    }


  ggplot(chrom.data)+
    geom_line(aes(x = rt , y = intensity ,
                  group = peaks.idx,
                  col = color),
              linewidth = 0.5,alpha = 0.8)+
    #scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(chrom.data$peaks.idx))))+
    ggsci::scale_fill_npg()+
    labs(x = "Retention time", y = "Intensity",col = "peaks")+
    theme_bw()->p

  if (!is.null(label_df)) {
    p <- p+ggrepel::geom_text_repel(data = label_df,
                                    aes(x= x,y=y,label = label),
                                    hjust = 0)

  }

  return(p)



}


#' @describeIn xcms_extension_feature calculate feature statistics
#' @description Extracts and adds median retention time, signal-to-noise ratio, and maximum intensity for each feature. While xcms::featureDefinitions() provides median mz and rt, this function calculates median values across all peaks within a feature: peakRtMin, peakRtMax, peakWidth, peakMzMin, peakMzMax, peakSN, peakMaxo, and polarity.
#' @param xcms.xcms XCMSnExp object with feature definitions and chromPeaks.
#'
#' @return XCMSnExp object with updated featureDefinitions containing additional statistics.
#' @export
#'
xcms_get_feature_def_stat <- function(xcms.xcms){

  feature.def <- xcms::featureDefinitions(xcms.xcms)
  peaks.data <- xcms::chromPeaks(xcms.xcms)

  .xcmsPeakDataMed <- function(x,peaks.data,key = "rtmax",fun = "median"){
    if (!key%in% colnames(peaks.data)) {
      return(NA)

    }
    x.peaks.data <- peaks.data[x,,drop=F]
    x.peaks.data <- x.peaks.data[!grepl(pattern = "CPM",x = rownames(x.peaks.data)),,drop=F]
    peak.key.value <- x.peaks.data[,key]
    eval(call(fun,peak.key.value,na.rm =T))
  }

  feature.def$peakRtMin <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,key = "rtmin",fun = "min")
  feature.def$peakRtMax <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,key = "rtmax",fun = "max")
  feature.def$peakWidth <- feature.def$peakRtMax-feature.def$peakRtMin
  feature.def$peakMzMin <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"mzmin",fun = "min")
  feature.def$peakMzMax <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"mzmax",fun = "max")
  feature.def$peakSN <-  sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"sn")
  feature.def$peakMaxo <-  sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"maxo")
  feature.def$polarity <- unique(.xcms_polarity(xcms.xcms))
  feature.def.df <- as.data.frame(feature.def)%>%
    dplyr::mutate(feature_id = rownames(.),
                  .before = mzmed)
  xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, feature.def.df)
  return(xcms.xcms)

}



xcms_get_feature_val_stat <- function(xcms.xcms) {

  xcms.pdata <- Biobase::pData(xcms.xcms)
  if (!"sampleNames" %in% colnames(xcms.pdata)) {
    xcms.pdata$sampleNames <- Biobase::sampleNames(xcms.xcms)
  }
  featureval <- xcms::featureValues(xcms.xcms)
  if (!is.null(colnames(featureval)) &&
      !any(xcms.pdata$sampleNames %in% colnames(featureval)) &&
      ncol(featureval) == nrow(xcms.pdata)) {
    colnames(featureval) <- xcms.pdata$sampleNames
  }
  if ("sample.type" %in% colnames(xcms.pdata)) {
    qc.rsd <- featureval[,xcms.pdata%>%dplyr::filter(sample.type =="QC")%>%
                 dplyr::pull(sampleNames),drop=F]%>%
      apply(1,function(x){
        sd(x,na.rm =T)/mean(x,na.rm=T)
      })
    sample.rsd <- featureval[,xcms.pdata%>%
                               dplyr::filter(sample.type =="Sample")%>%
                 dplyr::pull(sampleNames),drop=F]%>%
      apply(1,function(x){
        sd(x,na.rm =T)/mean(x,na.rm=T)
      })
  }else{
    qc.rsd <- NA
    sample.rsd <- NA

  }

  fdf <- as.data.frame(xcms::featureDefinitions(xcms.xcms))
  fdf$qc_rsd <- qc.rsd
  fdf$sample_rsd <- sample.rsd
  xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, fdf)
  return(xcms.xcms)
}


xcms_get_feature_stat <- function(xcms.xcms){
  xcms.xcms <- xcms.xcms %>%
    xcms_get_feature_def_stat()%>%
    xcms_get_feature_val_stat()
  return(xcms.xcms)
}


#' @describeIn xcms_extension_isotopologue identify isotopologues
#' @description Screens isotopologue peaks based on m/z and retention time differences, assigns isotopologue groups and seeds, and records results in featureDefinitions. Uses graph-based clustering to identify isotopologue networks.
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param iso_ele Isotope element string (e.g., `"[13]C"`) for mass difference calculation.
#' @param max_label Maximum number of isotope labels to consider (default 10).
#' @param ppm Mass accuracy tolerance in ppm (default 10).
#' @param rt.tol Retention time tolerance in seconds for grouping (default 5).
#' @param net.degree.ratio Ratio threshold for network degree to assign isotopologue seeds (default 0.3).
#'
#' @return XCMSnExp object with featureDefinitions updated with iso_seed, iso_count, and iso_connection_group columns.
#' @export
#'
xcms_get_feature_isotopologues <- function(xcms.xcms,
                                           iso_ele = "[13]C",
                                           max_label = 10,
                                           ppm = 10,
                                           rt.tol = 5,
                                           net.degree.ratio = 0.3){



  fdf.iso.connect <- get_xcms_feature_iso_connection(xcms.xcms,iso_ele,max_label,
                                                     ppm,rt.tol )


  ### assign isotopologues
  {

    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
      as.data.frame()
    xcms.fdf[,paste0("iso_seed")] <- NA
    xcms.fdf[,paste0("iso_count")] <- NA
    fdf.iso.igraph <- igraph::graph_from_data_frame(fdf.iso.connect)
    #fdf.iso.igraph <- igraph_filter_vertex(fdf.iso.igraph,degree(fdf.iso.igraph)>2)
    node.group <- igraph::components(fdf.iso.igraph)$membership
    xcms.fdf <- as.data.frame(xcms.fdf)
    rownames(xcms.fdf) <-xcms.fdf$feature_id
    message( length(unique(na.omit(node.group)))," iso-group"  )
    message( (length(node.group))," iso-features"  )

    for (i in seq_along(unique(node.group))) {

      #message(i)
      this.nodes <- names(which(node.group==i))
      this.fdf <- xcms.fdf[this.nodes,]
      this.iso <- fdf.iso.connect %>%
        dplyr::filter(from%in%this.nodes | to %in% this.nodes)

      this.igraph <- igraph::graph_from_data_frame(this.iso,vertices =this.fdf[,1:7] )
      #visNetwork::visIgraph(this.igraph)
      this.iso.assign <- get_iso_net_assign(this.igraph,net.degree.ratio = net.degree.ratio)
      xcms.fdf[names(this.iso.assign$iso.seed),
              "iso_seed"] <- this.iso.assign$iso.seed
      xcms.fdf[names(this.iso.assign$iso_count),
               "iso_count"] <- this.iso.assign$iso_count
      xcms.fdf[this.nodes,"iso_connection_group"] <- i
      this.fdf <- xcms.fdf[this.nodes,]


    }





  }


  ### save to featuredef
  {

    xcms.fdf.temp <- xcms::featureDefinitions(xcms.xcms)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    xcms.fdf.temp[,"iso_seed"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_seed"]
    xcms.fdf.temp[,"iso_count"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_count"]
    xcms.fdf.temp[,"iso_connection_group"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_connection_group"]
    xcms.fdf.temp -> xcms::featureDefinitions(xcms.xcms)
    message("Get ",
            sum(!is.na(xcms.fdf.temp[,"iso_count"])),
            " isotopologues")

  }

  return(xcms.xcms)

}

#' @describeIn xcms_extension_isotopologue identify isotopologues with multiple isotope tracers
#' @description TODO: unfinished. Similar to `xcms_get_feature_isotopologues` but supports multiple isotope labels simultaneously (e.g., `[13]C` and `[15]N`).
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param iso_ele Character vector of isotope element strings (e.g., `c("[13]C","[15]N")`).
#' @param max_label Named numeric vector of maximum labels per tracer, names must match `iso_ele`.
#' @param ppm Mass accuracy tolerance in ppm (default 5).
#' @param rt.tol Retention time tolerance in seconds (default 5).
#' @param net.degree.ratio Ratio threshold for network degree to assign isotopologue seeds (default 0.3).
#' @return XCMSnExp object with featureDefinitions updated with iso_seed, iso_count, iso_connection_group, and per-tracer iso_count_* columns.
#' @export
#'
xcms_get_feature_isotopologues_multi_tracer <- function(xcms.xcms,
                                            iso_ele = c("[13]C","[15]N"),
                                            max_label = c("[13]C" = 30,"[15]N" = 10),
                                            ppm = 5,
                                            rt.tol = 5,
                                            net.degree.ratio = 0.3){

  ### find multi-tracer connections
  {
    fdf.iso.connect <- get_xcms_feature_iso_connection_multi_tracer(
      xcms.xcms,
      iso_ele = iso_ele,
      max_label = max_label,
      ppm = ppm,
      rt.tol = rt.tol
    )
  }

  ### assign isotopologues
  {
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms) %>%
      as.data.frame()
    iso_ele_clean <- MSCC::get_ele_uniso(iso_ele)
    for (el in iso_ele_clean) {
      xcms.fdf[, paste0("iso_count_", el)] <- NA
    }
    xcms.fdf[, "iso_seed"] <- NA
    xcms.fdf[, "iso_count"] <- NA

    if (nrow(fdf.iso.connect) == 0) {
      message("No isotopologue connections found")
      fdf.iso.igraph <- igraph::make_empty_graph(n = 0)
      node.group <- integer(0)
    } else {
      fdf.iso.igraph <- igraph::graph_from_data_frame(fdf.iso.connect)
      node.group <- igraph::components(fdf.iso.igraph)$membership
    }

    xcms.fdf <- as.data.frame(xcms.fdf)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    message(length(unique(na.omit(node.group))), " iso-group")
    message(length(node.group), " iso-features")

    for (i in seq_along(unique(node.group))) {
      this.nodes <- names(which(node.group == i))
      this.fdf <- xcms.fdf[this.nodes, ]
      this.iso <- fdf.iso.connect[from %in% this.nodes | to %in% this.nodes]
      this.igraph <- igraph::graph_from_data_frame(this.iso, vertices = this.fdf[, 1:7])
      this.iso.assign <- get_iso_net_assign(this.igraph, net.degree.ratio = net.degree.ratio)

      xcms.fdf[names(this.iso.assign$iso.seed), "iso_seed"] <- this.iso.assign$iso.seed
      xcms.fdf[names(this.iso.assign$iso_count), "iso_count"] <- this.iso.assign$iso_count
      xcms.fdf[this.nodes, "iso_connection_group"] <- i

      ### assign per-tracer counts from the seed
      this.connect <- fdf.iso.connect[from %in% this.nodes & to %in% this.nodes]
      for (el in iso_ele_clean) {
        col_name <- paste0("closest.iso.count_", el)
        if (col_name %in% colnames(this.connect)) {
          seed.fid <- unique(na.omit(xcms.fdf[this.nodes, "iso_seed"]))
          if (length(seed.fid) == 1) {
            seed.idx <- which(this.connect$from == seed.fid | this.connect$to == seed.fid)
            if (length(seed.idx) > 0) {
              for (node in this.nodes) {
                if (node == seed.fid) {
                  xcms.fdf[node, paste0("iso_count_", el)] <- 0
                  next
                }
                edge.idx <- which((this.connect$from == node & this.connect$to == seed.fid) |
                                    (this.connect$from == seed.fid & this.connect$to == node))
                if (length(edge.idx) > 0) {
                  xcms.fdf[node, paste0("iso_count_", el)] <- abs(this.connect[edge.idx[1], col_name])
                }
              }
            }
          }
        }
      }
    }
  }

  ### save to featureDefinitions
  {
    xcms.fdf.temp <- xcms::featureDefinitions(xcms.xcms)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    xcms.fdf.temp[, "iso_seed"] <- xcms.fdf[xcms.fdf.temp$feature_id, "iso_seed"]
    xcms.fdf.temp[, "iso_count"] <- xcms.fdf[xcms.fdf.temp$feature_id, "iso_count"]
    xcms.fdf.temp[, "iso_connection_group"] <- xcms.fdf[xcms.fdf.temp$feature_id, "iso_connection_group"]
    for (el in iso_ele_clean) {
      xcms.fdf.temp[, paste0("iso_count_", el)] <- xcms.fdf[xcms.fdf.temp$feature_id, paste0("iso_count_", el)]
    }
    xcms.fdf.temp -> xcms::featureDefinitions(xcms.xcms)
    message("Get ",
            sum(!is.na(xcms.fdf.temp[, "iso_count"])),
            " isotopologues")
  }

  return(xcms.xcms)

}


#' @title TODO: unfinished. Build isotope mass shift grid for multi-tracer
#' @description Generates all label-count combinations across tracers with their mass shifts.
#' @param iso_ele Character vector of isotope element strings.
#' @param max_label Named numeric vector of maximum labels per tracer.
#' @return Data.frame with columns for each tracer's label count, total.count, and mass.shift.
#'
#' @title Pair xcms features within retention-time tolerance
#' @description Builds all directed feature pairs from `featureDefinitions` whose
#'   retention times differ by less than `rt.tol`. Used as the RT filter before
#'   isotope mass-shift matching.
#' @param xcms.xcms \code{XCMSnExp} with feature definitions.
#' @param rt.tol Retention time tolerance in seconds (default 5).
#' @return Data frame with integer `from`/`to` row indices, feature ids, m/z, rt,
#'   `mz.diff`, `rt.diff`, and `mz.mean` (mean m/z of the pair).
#' @export

get_xcms_feature_connect <- function(xcms.xcms,rt.tol = 5){


  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)

  {
    rt <- xcms.fdf$rtmed
    x <- data.table(rt = rt, start = rt,end = rt)
    y <- x[,.(id = seq_along(rt),rt,start = rt - rt.tol, end = rt + rt.tol)]
    data.table::setkey(y,start,end)
    rtm <- data.table::foverlaps(x, y, type="any", which=TRUE)
    rtm <- rtm[,.(xid,yid = y$id[yid])]
  }


  {
    xcms.net <- rtm[,.(from = xid, to = yid)]
    xcms.net <- xcms.net[from < to ][
      , rt.diff := (xcms.fdf$rtmed[to]-xcms.fdf$rtmed[from]) ][
        abs(rt.diff) < rt.tol,][
          , c("from.mz","to.mz") := .( xcms.fdf$mzmed[from], xcms.fdf$mzmed[to])][
            ,c("mz.diff","mz.mean") := .(to.mz-from.mz, (from.mz+to.mz)/2)]
  }


  return(xcms.net)



}




get_xcms_feature_isotope_grid_multi_tracer <- function(iso_ele, max_label) {

  iso.chemforms <- character()
  iso.counts <- numeric()

  for (i in seq_along(iso_ele)) {
    ele <- iso_ele[i]
    elem.symbol <- stringr::str_extract(ele, "[[:alpha:]]+")
    chemform <- paste0(ele, 1, elem.symbol, -1)
    iso.chemforms[i] <- chemform
    iso.counts[i] <- MSCC::chemform_mz(chemform, 0)
  }

  label.ranges <- lapply(max_label, function(m) -m:m)
  grid <- expand.grid(label.ranges)
  colnames(grid) <- iso_ele

  mass.shift <- as.matrix(grid) %*% iso.counts

  data.frame(
    grid,
    total.count = rowSums(abs(grid)),
    mass.shift = as.vector(mass.shift)
  )
}


get_xcms_feature_iso_connection_multi_tracer <- function(xcms.xcms,
                                                          iso_ele,
                                                          max_label,
                                                          ppm = 10,
                                                          rt.tol = 5) {

  fdf.connect <- get_xcms_feature_connect(xcms.xcms, rt.tol = rt.tol)

  iso.grid <- get_xcms_feature_isotope_grid_multi_tracer(iso_ele, max_label)

  match.res <- match_mz_foverlaps(mz1 = fdf.connect$mz.diff,
                                   mz2 = iso.grid$mass.shift,
                                   ppm.base = fdf.connect$mz.mean,
                                   ppm = ppm)

  iso_ele_clean <- MSCC::get_ele_uniso(iso_ele)
  for (i in seq_along(iso_ele)) {
    fdf.connect[match.res$ion1, paste0("closest.iso.count_", iso_ele_clean[i])] <- iso.grid[match.res$ion2, iso_ele[i]]
  }
  fdf.connect[match.res$ion1, "closest.iso.count"] <- iso.grid$total.count[match.res$ion2]
  fdf.connect[match.res$ion1, "closest.iso.mz"] <- iso.grid$mass.shift[match.res$ion2]
  fdf.connect[match.res$ion1, "mz.error"] <- abs(match.res$mz.ppm)

  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  fdf.iso.connect <- fdf.connect[!is.na(closest.iso.count) & closest.iso.count != 0]
  data.table::setorder(fdf.iso.connect, mz.error)
  fdf.iso.connect <- fdf.iso.connect[!duplicated(fdf.iso.connect, by = c("from", "to"))]
  fdf.iso.connect[, `:=`(from = xcms.fdf$feature_id[from],
                          to = xcms.fdf$feature_id[to])]

  return(fdf.iso.connect)

}


get_xcms_feature_iso_connection <- function(xcms.xcms,
                                            iso_ele,
                                            max_label = 10,
                                            ppm = 10,
                                            rt.tol = 5){


  {

    fdf.connect <- get_xcms_feature_connect(xcms.xcms, rt.tol = rt.tol)

    #isotope <- "[13]C"
    iso.chemform <- paste0(iso_ele,1,
                           str_extract(string = iso_ele,pattern = "[[:alpha:]]+"),-1)
    iso.count <- -max_label:max_label
    iso.mz <- MSCC::chemform_mz(iso.chemform,0)*iso.count

    closest.iso.count <- sapply(fdf.connect$mz.diff, function(x){
      iso.count[which.min(abs(iso.mz-x))]
    } )
    fdf.connect <- fdf.connect%>%
      #rowwise()%>%
      dplyr::mutate(closest.iso.count = closest.iso.count)%>%
      dplyr::mutate(closest.iso.mz = iso.mz[match(closest.iso.count,iso.count)],
                    mz.error = abs(mz.diff-closest.iso.mz),
                    is.iso = mz.error/(from.mz+to.mz)*2 < ppm*1e-6)%>%
      dplyr::ungroup()

    fdf.iso.connect <- fdf.connect%>%
      dplyr::filter(is.iso,closest.iso.count != 0)%>%
      dplyr::group_by(from,closest.iso.count)%>%
      dplyr::slice_min(mz.error)%>%
      dplyr::ungroup()%>%
      dplyr::group_by(to,closest.iso.count)%>%
      dplyr::slice_min(mz.error)%>%
      dplyr::ungroup()%>%
      dplyr::group_by(from,closest.iso.count)%>%
      dplyr::slice_min(mz.error)%>%
      dplyr::ungroup()%>%
      dplyr::mutate(from = from.fid,to = to.fid)

  }

  return(fdf.iso.connect)

}



#' @describeIn xcms_extension_isotopologue calculate traced-isotopologue labeling ratios
#' @description Calculates isotopologue-to-seed ratios and determines traced
#' isotopologues (label-enriched isotopologues) using one of two methods:
#' \itemize{
#'   \item \code{untraced_compare} (legacy): compare traced and untraced sample sources.
#'   \item \code{natural_based}: compare observed ratio to theoretical natural isotope
#'   ratio derived from \code{MSCC::chemform_isotopes_pattern_enviPat()}.
#' }
#' Results are written to featureDefinitions as \code{is_labeled} and
#' \code{Ratio_to_seed_*} columns (same output contract as the legacy function).
#' @param xcms.xcms XCMSnExp object with isotopologue assignments.
#' @param iso_ele Isotope element string (e.g., `"[13]C"`) used for labeling.
#' @param method Labeling method: \code{"untraced_compare"} or
#'   \code{"natural_based"}. (Legacy aliases \code{"method1"} /
#'   \code{"method2"} are also accepted.)
#' @param ... Additional arguments passed to internal functions.
#'
#' @return XCMSnExp object with featureDefinitions updated with is_labeled column and Ratio_to_seed_* columns.
#' @export
#'
xcms_get_feature_traced_isotopologue <- function(xcms.xcms,
                                                 iso_ele = "[13]C",
                                                 method = c("untraced_compare", "natural_based")[1],
                                                 ...){

  method <- as.character(method)[1]
  if (identical(method, "method1")) method <- "untraced_compare"
  if (identical(method, "method2")) method <- "natural_based"
  method <- match.arg(method, c("untraced_compare", "natural_based"))

  ### feature data
  {
    xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
    xcms.se <- xcms.se[,xcms.se$sample.type=="Sample"]
  }

  ### calc iso ratio to seed and aggregate by sample.source
  {
    xcms.ratio.to.seed <- get_xcms_iso_fraction(xcms.xcms)
    xcms.ratio.to.seed <- apply(xcms.ratio.to.seed, 1,
                                function(x){
                                  mean_f(x, f = xcms.se$sample.source, simplify = F, na.rm = T)
                                }) %>%
      do.call(bind_rows, .) %>%
      as.matrix()
  }

  ### determine traced-isotopologue labels
  {
    is.iso <- xcms.se$isotope_tracer %in% iso_ele
    sample.source.iso <- unique(xcms.se$sample.source[is.iso])
    sample.source.uniso <- unique(xcms.se$sample.source[!is.iso])

    if (method == "untraced_compare") {
      if (length(sample.source.uniso) == 0) {
        cli::cli_alert_warning("No untraced sample.source found; cannot determine traced isotopologues with {.code untraced_compare}.")
        is_labeled <- rep(NA, nrow(xcms.ratio.to.seed))
      } else {
      is_labeled <- apply(xcms.ratio.to.seed, 1, function(x){
        any(x[sample.source.iso] > mean(x[sample.source.uniso], na.rm = TRUE), na.rm = TRUE)
      })
      }
    } else {
      # method2: compare observed ratio to theoretical natural ratio.
      xcms.fdf <- xcms::featureDefinitions(xcms.xcms) %>% as.data.frame()

      # Expected natural ratio for one feature row.
      .expected_ratio <- function(formula, iso_count, iso_ele, thresh = 1e-6) {
        if (is.na(formula) || !nzchar(formula) || is.na(iso_count)) return(NA_real_)
        pat <- tryCatch(
          MSCC::chemform_isotopes_pattern_enviPat(formula, thresh = thresh),
          error = function(e) NULL
        )
        if (is.null(pat) || !nrow(pat) || !all(c("isotope_element", "abundance") %in% colnames(pat))) {
          return(NA_real_)
        }
        ie <- as.character(pat$isotope_element)
        ie[is.na(ie)] <- ""
        keep <- ie == "" | grepl(iso_ele, ie, fixed = TRUE)
        pat <- pat[keep, , drop = FALSE]
        if (!nrow(pat)) return(NA_real_)
        iso_chr <- as.character(pat$isotope_element)
        iso_chr[is.na(iso_chr)] <- ""
        iso_num <- ifelse(
          iso_chr == "",
          0L,
          suppressWarnings(as.integer(gsub("[^0-9]", "", gsub(iso_ele, "", iso_chr, fixed = TRUE))))
        )
        iso_num[is.na(iso_num)] <- 0L
        a0 <- sum(pat$abundance[iso_num == 0], na.rm = TRUE)
        ak <- sum(pat$abundance[iso_num == as.integer(iso_count)], na.rm = TRUE)
        if (!is.finite(a0) || a0 <= 0) return(NA_real_)
        ak / a0
      }

      fids <- if ("feature_id" %in% colnames(xcms.fdf)) as.character(xcms.fdf$feature_id) else rownames(xcms.fdf)
      formula_vec <- if ("formula" %in% colnames(xcms.fdf)) as.character(xcms.fdf$formula) else rep(NA_character_, nrow(xcms.fdf))
      iso_count_vec <- if ("iso_count" %in% colnames(xcms.fdf)) suppressWarnings(as.integer(xcms.fdf$iso_count)) else rep(NA_integer_, nrow(xcms.fdf))
      expected <- vapply(seq_len(nrow(xcms.fdf)), function(i) {
        .expected_ratio(formula_vec[i], iso_count_vec[i], iso_ele)
      }, numeric(1))
      names(expected) <- fids

      # If no traced sample.source is available, fallback to all sample.source.
      use.sources <- sample.source.iso
      if (!length(use.sources)) {
        use.sources <- colnames(xcms.ratio.to.seed)
      }
      use.sources <- intersect(use.sources, colnames(xcms.ratio.to.seed))

      is_labeled <- rep(NA, nrow(xcms.ratio.to.seed))
      if (length(use.sources)) {
        is_labeled <- vapply(seq_len(nrow(xcms.ratio.to.seed)), function(i) {
          obs <- xcms.ratio.to.seed[i, use.sources, drop = TRUE]
          fid <- rownames(xcms.ratio.to.seed)[i]
          exp_i <- expected[fid]
          exp_i <- if (length(exp_i)) exp_i[[1]] else NA_real_
          if (!is.finite(exp_i)) return(NA)
          any(obs > exp_i, na.rm = TRUE)
        }, logical(1))
      }
    }
  }

  ### import to xcms
  {
    xcms.fda <- xcms::featureDefinitions(xcms.xcms)
    xcms.fda$is_labeled <- is_labeled
    colnames(xcms.ratio.to.seed) <- paste0("Ratio_to_seed_", colnames(xcms.ratio.to.seed))
    xcms.fda[, colnames(xcms.ratio.to.seed)] <- xcms.ratio.to.seed
    xcms.fda -> xcms::featureDefinitions(xcms.xcms)
    message("Get ",
            sum(xcms.fda$is_labeled, na.rm = TRUE),
            " isotope label")
  }

  return(xcms.xcms)
}

#' @title Xcms Get Feature Isotope Label (Deprecated)
#' @description Deprecated wrapper of \code{xcms_get_feature_traced_isotopologue()}.
#' @param xcms.xcms XCMSnExp object with isotopologue assignments.
#' @param iso_ele Isotope element string (e.g., `"[13]C"`).
#' @param ... Additional arguments.
#' @return XCMSnExp object.
#' @export
xcms_get_feature_isotope_label <- function(xcms.xcms,
                                           iso_ele = "[13]C",
                                           ...){
  .Deprecated("xcms_get_feature_traced_isotopologue",
              package = "MSdev",
              msg = "xcms_get_feature_isotope_label is deprecated. Use xcms_get_feature_traced_isotopologue().")
  xcms_get_feature_traced_isotopologue(xcms.xcms = xcms.xcms,
                                       iso_ele = iso_ele,
                                       method = "untraced_compare",
                                       ...)
}


get_xcms_isotopologues_report <- function(xcms.xcms){



}


#' @describeIn xcms_extension_isotopologue calculate isotopologue fractions
#' @description Calculates the fraction of isotopologue intensities relative to their seed feature intensities for each sample. Returns a matrix of fractions without natural abundance adjustment.
#' @param xcms.xcms XCMSnExp object with isotopologue assignments (iso_seed column).
#'
#' @return Matrix with rows as features and columns as samples, containing intensity ratios to seed features.
#' @export
#'
get_xcms_iso_fraction <- function(xcms.xcms){


  ### feature data
  {
    xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
    #xcms.se <- xcms.se[,xcms.se$sample.type=="Sample"]
     xcms.rda <- rowData(xcms.se)%>%
      as.data.frame()
    xcms.val <- assay(xcms.se)
  }

  ###calc iso ratio to seed
  {
    xcms.ratio.to.seed <- xcms.val
    xcms.ratio.to.seed[,] <-NA
    xcms.fseed <- xcms.rda$iso_seed %>%
      unique()%>%na.omit()
    for (i in seq_unique(xcms.fseed)) {
      this.fid <- xcms.fseed[i]
      this.iso <- xcms.rda%>%
        dplyr::filter(iso_seed %in% this.fid)
      this.matrix <- xcms.val[this.iso$feature_id,,drop = F]
      this.matrix <- t(t(this.matrix)/this.matrix[this.fid,])
      xcms.ratio.to.seed[rownames(this.matrix),] <- this.matrix
    }

    xcms.ratio.to.seed[is.nan(xcms.ratio.to.seed)] <- 0

  }
  return(xcms.ratio.to.seed)

}


#' @describeIn xcms_extension_annotation match features to compound database
#' @description Matches features in an XCMSnExp object to compounds in a CompoundDb database using m/z and retention time tolerance. Calculates adduct masses for each compound and finds matches within specified ppm error. Results are stored as candidate lists in featureDefinitions.
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param cpdb CompoundDb object containing compound database.
#' @param mz.ppm Numeric. Mass accuracy tolerance in parts per million (default 10).
#' @param rt.tol Numeric. Retention time tolerance in seconds (default Inf, no RT filtering).
#' @param selected_adduct Character vector of adducts to consider (default from MSCC::adduct.table$Adduct).
#' @param ... Additional arguments passed to internal functions.
#'
#' @return XCMSnExp object with featureDefinitions updated with candidate.id, candidate.formula, candidate.adduct, and candidate.mz columns.
#' @export
#'

xcms_get_feature_ms1_candidate <- function(xcms.xcms ,
                                           cpdb,
                                           mz.ppm= 10,
                                           rt.tol = Inf,
                                           selected_adduct = MSCC::adduct.table$Adduct,
                                           ...){


  ### calc adduct and filter range
  cpdbt <- compounds(cpdb, columns = CompoundDb::compoundVariables(cpdb,includeId =T))
  #cpdbt <- dplyr::filter(cpdbt,lipidclass=="VAE")
  if ("has_sp"%in% colnames(cpdbt))  cpdbt <- cpdbt[cpdbt$has_sp>0,]
  cpdbt$formula <- MSCC::chemform_formate(cpdbt$formula)

  adducts <-  MSCC::chemform_adduct_check(selected_adduct)%>%
    dplyr::mutate(polarity = case_when(Ion_mode == "negative"~0,T~1))%>%
    dplyr::filter(polarity %in% polarity(xcms.xcms))
  cp.adduct <- MSCC::chemform_adduct(cpdbt$formula,
                                     adducts$adduct.formated,
                                     value = "all" )
  cp.adduct <- cp.adduct%>%
    dplyr::mutate(compound_id=cpdbt$compound_id[id]  )%>%
    dplyr::filter( findInterval(chemform.adduct.mz,
                                mzrange(xcms.xcms))==1)

  ### match database
  xcms.featuredef <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()

  matched.df <- match_mz_foverlaps(mz1 = xcms.featuredef$mzmed,
                            mz2 = cp.adduct$chemform.adduct.mz,
                            ppm = mz.ppm)
  matched.df2 <- cbind( matched.df,cp.adduct[matched.df$ion2,])
  xcms.featuredef$candidate.id <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$compound_id[as.numeric(idx)]
  })
  xcms.featuredef$candidate.formula <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform[as.numeric(idx)]
  })
  xcms.featuredef$candidate.adduct <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$adduct[as.numeric(idx)]
  })
  xcms.featuredef$candidate.mz <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform.adduct.mz[as.numeric(idx)]
  })

  xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, xcms.featuredef)

  return(xcms.xcms)

}


xcms_get_feature_ms2_score <- function(xcms.xcms ,
                                       cpdb,
                                       sp.ms2,
                                       ...){



  ### no ms2
  {
    if (length(sp.ms2)==0) {
      xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
      xcms.fdf$score.ms2 <- lapply(xcms.fdf$candidate.id,function(x){
        rep(0,length(x))
      })

      xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, xcms.fdf)
      return(xcms.xcms)
    }

  }


  ### load spectra database
  {
    Spectra_database <- Spectra::Spectra(cpdb)
    Spectra_database <- Spectra_set_MEM_backend(Spectra_database)
    Spectra_database <- filterPolarity(Spectra_database,
                                       unique(polarity(xcms.xcms)))

    #Spectra::spectraNames(Spectra_database) <- Spectra_database$compound_id
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  }

  ### sp process
  {
    Spectra_database <- Spectra_database%>%
      filterSpectra_below_PrecursorMz()%>%
      normalizeSpectra(norm_to = "max")%>%
      filterSpectraIntensity(ratio = 0.05)%>%
      Spectra::applyProcessing()
    if(length(sp.ms2)!=0){
      sp.ms2 <- sp.ms2%>%
        filterSpectra_below_PrecursorMz()%>%
        normalizeSpectra(norm_to = "max")%>%
        filterSpectraIntensity(ratio = 0.05)%>%
        Spectra_set_MEM_backend()%>%
        Spectra::applyProcessing()
      #if ("from_iso" %in% spectraVariables(sp.ms2)) {
      #  sp.ms2 <- sp.ms2[!sp.ms2$from_iso]
      #}
    }
    if (!all(unlist(xcms.fdf$candidate.id)%in%
             Spectra_database$compound_id)) {
      sp.empty <- makeEmptySpectra(compound_id= setdiff(unlist(xcms.fdf$candidate.id),
                                                        Spectra_database$compound_id))
      Spectra_database <- c(Spectra_database,sp.empty)
    }
  }

  ### sp ms2 split (ms2_id = character sp_id / spectraNames)
  {
    sp.exp <- sapply(1:nrow(xcms.fdf),function(i){

      x <- xcms.fdf$ms2_id[[i]]
      if (!length(x)) return(NULL)
      sp.ms2[match(x,Spectra::spectraNames(sp.ms2))]
      #if (length(x)==0) {
      #  sp <- makeSpectra(xcms.fdf$mzmed[i],
      #                    xcms.fdf$rtmed[i])
      #}else
      #  sp <- list(sp.ms2[x])
      #return(sp)
    })
    names(sp.exp) <- xcms.fdf$feature_id


  }

  ### sp ref split
  {
    sp.split.df <- lapply(1:nrow(xcms.fdf),
                          function( i ){
                            i.candi <- xcms.fdf$candidate.id[[i]]
                            i.adduct <- xcms.fdf$candidate.adduct[[i]]
                            if (!length(i.candi)) return(NULL)
                            i.df <- lapply(i.candi,
                                           function(x){which(Spectra_database$compound_id==x)}
                            )%>%
                              `names<-`(xcms.fdf$candidate.id[[i]])%>%
                              unlist_to_df(name_to = "compound_id",
                                           value_to = "sp_id")
                            i.df$adduct <- i.adduct[match(i.df$compound_id,i.candi)]
                            i.df
                          })%>%
      data.table::rbindlist(use.names = T,idcol = "feature_id")
    sp.ref <- Spectra_database[sp.split.df$sp_id]
    sp.ref$adduct <- sp.split.df$adduct
    sp.ref <- split(sp.ref,xcms.fdf$feature_id[sp.split.df$feature_id])



  }

  ### output all candidate score
  {
    .f <- function(expSpec,refSpec,...){
      if (is.null(expSpec)) {
        if (is.null(refSpec)) {
          return(NULL)
        }else{
          x <- unique(paste0(refSpec$compound_id,"_",refSpec$adduct))
          y <- rep(0,length(x))
          names(y )<-x
          return(y)
        }
      }
      if (is.null(refSpec)) return(NULL)
      scorem <- Spectra::compareSpectra(expSpec,refSpec,,FUN = MsCoreUtils::ndotproduct, m = 2)
      dim(scorem) <- c(length(expSpec),length(refSpec))
      scorem[is.infinite(scorem)|is.na(scorem )] <- 0
      scores <- apply(scorem,2,max,na.rm=T)
      (mean_f(scores,paste0(refSpec$compound_id,"_",refSpec$adduct)))
    }
    xcms.fdf$score.ms2 <- BiocParallel::bplapply(1:length(sp.exp),
                                                       function(i){
        fid <-xcms.fdf$feature_id[i]
      s <- .f(expSpec = sp.exp[[i]], refSpec = sp.ref[[fid]])
      s[paste0(xcms.fdf$candidate.id[[i]],"_",
               xcms.fdf$candidate.adduct[[i]])]%>%
        unname()
    },BPPARAM = BiocParallel::SerialParam(
      progressbar = T))


  }


  xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, xcms.fdf)

  return(xcms.xcms)

}

xcms_get_feature_isopattern_score <- function(xcms.xcms,
                                              ppm = 10,
                                              calc_isopattern_score = T,
                                              BPPARAM = SerialParam(progressbar = T)){

  ### data to calc isopattern
  {
    xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
    xcms.se <- xcms.se[,is.na(xcms.se$isotope_tracer)]
    xcms.se <- xcms.se[,!xcms.se$sample.type%in% "Blank"]
  }

  ### calculate iso-pattern score
  {

    if (calc_isopattern_score) {

      iso.score <- bplapply(seq_len(nrow(xcms.se)),
             FUN = function(i,xcms.se,ppm){
               xcms.fdf <- SummarizedExperiment::rowData(xcms.se)
               xcms.se.temp <- xcms.se[MSdev:::between.range(xcms.fdf$rtmed,
                                                     c(xcms.fdf$rtmed[i]-10,
                                                       xcms.fdf$rtmed[i]+10)),]
               formulas <-mapply( MSCC::chemform_adduct,
                                 chemform = xcms.fdf$candidate.formula[[i]],
                                 adduct =  xcms.fdf$candidate.adduct[[i]],
                                 value = "chemform")
               MSCC::get_isopattern_score(formula = formulas,
                                    mzs = rowData(xcms.se.temp)$mzmed,
                                    int_matrix = assay(xcms.se.temp),
                                    ppm = ppm)
             },xcms.se=xcms.se,ppm=ppm,
             BPPARAM =BPPARAM)

    }else{iso.score <- 0}



  }


  ### return
  {
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
    xcms.fdf$score.isopattern <-iso.score
    xcms.fdf -> xcms::featureDefinitions(xcms.xcms)
    return(xcms.xcms)
  }


}



get_xcms_feature_all_candidate <- function(xcms.xcms){

  xcms.fdf <- get_xcms_feature_definitions(xcms.xcms)
  candi.rda <- xcms.fdf%>%
    dplyr::mutate(candidate.n = lengths(candidate))
  candi.rda.split <- candi.rda[rep(candi.rda$feature_id,
                                   candi.rda$candidate.n),]%>%
    dplyr::group_by(feature_id)%>%
    dplyr::mutate(temp_id = 1:n())%>%
    dplyr::rowwise()%>%
    dplyr::mutate(MSDB_id = candidate[[temp_id]],
                  adduct = candidate.adduct[[temp_id]],
                  mz_ref = candidate.mz[[temp_id]],
                  score = candidate.score[[temp_id]])%>%
    dplyr::ungroup()%>%
    dplyr::select(-c(candidate,candidate.adduct,candidate.mz,candidate.score))

  return(candi.rda.split)

}

xcms_get_feature_annotation <- function(xcms.xcms,
                                        cpdb,
                                        cpdb.keys = c("name","formula","smiles"),
                                        weight_mz = 0.1,
                                        weight_ms2 = 0.7,
                                        weight_isopattern =0.2,
                                        ...){


  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  xcms.fdf$compound_id <- NA
  xcms.fdf$adduct <- NA
  xcms.fdf$score <- NA
  xcms.fdf$mz_ref <- NA



  {
    xcms.candi.dt <-lapply(which(lengths(xcms.fdf$candidate.id)!=0), function(i) {
        data.table(
          feature_id = xcms.fdf$feature_id[i],
          mz = xcms.fdf$mzmed[i],
          candidate.id = xcms.fdf$candidate.id[[i]],
          candidate.adduct = xcms.fdf$candidate.adduct[[i]],
          candidate.formula = xcms.fdf$candidate.formula[[i]],
          candidate.mz = xcms.fdf$candidate.mz[[i]],
          score.ms2 =  xcms.fdf$score.ms2[[i]],
          score.isopattern = xcms.fdf$score.isopattern[[i]]
        )
      }) %>% rbindlist
    xcms.candi.dt <- xcms.candi.dt[
      ,score.isopattern:= ifelse(is.na(score.isopattern),0,score.isopattern)][
        ,score.isopattern:= ifelse(is.nan(score.isopattern),0,score.isopattern)][
          ,score.mz := 1- abs(candidate.mz-mz)/mz *1e6 / 20 ][ ### (1-ppm/10)
            , score.mz := ifelse(score.mz < 0, 0, score.mz)][
              ,score := score.isopattern * weight_isopattern + score.ms2*weight_ms2+score.mz * weight_mz]

    xcms.candi.dt.max <- xcms.candi.dt[, .SD[which.max(score)], by = feature_id]
    data.table::setnames(xcms.candi.dt.max,
                          old =  c("candidate.id","candidate.adduct","candidate.formula","candidate.mz"),
                          new = c("compound_id","adduct","formula","mz_ref") )
    data.table::setkey(xcms.candi.dt.max,feature_id)

    xcms.candi <- xcms.candi.dt.max[xcms.fdf$feature_id]
    xcms.fdf$compound_id <- xcms.candi$compound_id
    xcms.fdf$adduct <- xcms.candi$adduct
    xcms.fdf$formula <- xcms.candi$formula
    xcms.fdf$score <- xcms.candi$score
    xcms.fdf$mz_ref <- xcms.candi$mz_ref
    ## keep per-candidate score lists for candidate.se expansion
    xcms.fdf$candidate.score.ms2 <- xcms.fdf$score.ms2
    xcms.fdf$candidate.score.isopattern <- xcms.fdf$score.isopattern
    ## best-candidate scalars for feature / metabolite export
    xcms.fdf$score.ms2 <- xcms.candi$score.ms2
    xcms.fdf$score.isopattern <- xcms.candi$score.isopattern
    xcms.fdf$score.mz <- xcms.candi$score.mz

  }


  ### Comp info
  {

    dbinfo <- get_CompDb_info(cpdb,
                              xcms.fdf$compound_id,
                              keys = cpdb.keys)
    xcms.fdf[,colnames(dbinfo)] <- dbinfo

  }


  xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, xcms.fdf)

  return(xcms.xcms)



}


#' Match theoretical isotopes to xcms feature definitions
#'
#' @param isotopes_table Data frame from MSCC isotope pattern helpers
#'   (needs column \code{m.z}; optional \code{rt}).
#' @param featuredef xcms feature definition data frame (\code{mzmed},
#'   optional \code{rtmed}/\code{feature_id}).
#' @param mz.ppm m/z tolerance in ppm.
#' @param rt.tol RT tolerance in seconds if both sides have RT.
#'
#' @return Matched isotope-feature table.
#' @keywords internal
match_isotopes_to_featuredef <- function(isotopes_table, featuredef, mz.ppm = 10, rt.tol = 10) {
  if (is.null(isotopes_table) || !nrow(isotopes_table) || is.null(featuredef) || !nrow(featuredef)) {
    return(data.frame())
  }
  iso <- as.data.frame(isotopes_table, stringsAsFactors = FALSE)
  fdf <- as.data.frame(featuredef, stringsAsFactors = FALSE)
  if (!"m.z" %in% colnames(iso) || !"mzmed" %in% colnames(fdf)) return(data.frame())
  if (!"feature_id" %in% colnames(fdf)) fdf$feature_id <- rownames(fdf)
  if (!"rtmed" %in% colnames(fdf)) fdf$rtmed <- NA_real_
  if (!"rt" %in% colnames(iso)) iso$rt <- NA_real_

  out <- lapply(seq_len(nrow(iso)), function(i) {
    mz <- suppressWarnings(as.numeric(iso$m.z[i]))
    rt <- suppressWarnings(as.numeric(iso$rt[i]))
    if (!is.finite(mz)) return(NULL)
    mz_err <- abs(suppressWarnings(as.numeric(fdf$mzmed)) - mz) / mz
    hit <- mz_err <= mz.ppm * 1e-6
    if (is.finite(rt)) {
      rt_err <- abs(suppressWarnings(as.numeric(fdf$rtmed)) - rt)
      hit <- hit & (rt_err <= rt.tol)
    }
    idx <- which(hit)
    if (!length(idx)) return(NULL)
    cbind(
      iso[rep(i, length(idx)), , drop = FALSE],
      fdf[idx, , drop = FALSE],
      mz_error_ppm = mz_err[idx] * 1e6,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(data.frame())
  dplyr::bind_rows(out)
}

#' Append feature intensity matrix to isotope matches
#'
#' @param matched_table Output from \code{match_isotopes_to_featuredef()}.
#' @param featureval Numeric matrix of feature intensities (e.g. from
#'   xcms \code{featureValues()}), with feature IDs as row names.
#'
#' @return Data frame with matched rows and joined intensity columns.
#' @keywords internal
match_isotopes_to_featureval <- function(matched_table, featureval) {
  if (is.null(matched_table) || !nrow(matched_table)) return(data.frame())
  out <- as.data.frame(matched_table, stringsAsFactors = FALSE)
  if (is.null(featureval) || !is.matrix(featureval) || !nrow(featureval)) return(out)
  fv <- as.data.frame(featureval, stringsAsFactors = FALSE)
  fv$feature_id <- rownames(featureval)
  val_cols <- setdiff(colnames(fv), "feature_id")
  colnames(fv)[match(val_cols, colnames(fv))] <- paste0("featureval_", val_cols)
  dplyr::left_join(out, fv, by = "feature_id")
}

#' MS1 candidate matching on a featureDefinitions-style data.frame
#' @param fdf data.frame with at least \code{mzmed} (and preferably
#'   \code{feature_id})
#' @param cpdb CompoundDb object
#' @param polarity integer polarity/polarities to keep for adducts
#' @param mz.ppm m/z tolerance in ppm
#' @param mz_range numeric length-2 m/z range for filtering DB adducts; default
#'   \code{range(fdf$mzmed)}
#' @param selected_adduct adduct names
#' @param ... unused
#' @return \code{fdf} with candidate.* list columns
#' @keywords internal
fdf_get_ms1_candidate <- function(fdf,
                                  cpdb,
                                  polarity,
                                  mz.ppm = 10,
                                  mz_range = NULL,
                                  selected_adduct = MSCC::adduct.table$Adduct,
                                  ...) {
  fdf <- as.data.frame(fdf, stringsAsFactors = FALSE)
  if (!nrow(fdf)) {
    return(fdf)
  }
  if (!"feature_id" %in% names(fdf)) {
    fdf$feature_id <- rownames(fdf)
  }
  polarity <- as.integer(unique(polarity))
  if (is.null(mz_range)) {
    mz_range <- range(as.numeric(fdf$mzmed), na.rm = TRUE)
  }

  cpdbt <- compounds(cpdb, columns = CompoundDb::compoundVariables(cpdb, includeId = TRUE))
  if ("has_sp" %in% colnames(cpdbt)) {
    cpdbt <- cpdbt[cpdbt$has_sp > 0, ]
  }
  cpdbt$formula <- MSCC::chemform_formate(cpdbt$formula)

  polarity_keep <- polarity
  adducts <- MSCC::chemform_adduct_check(selected_adduct) %>%
    dplyr::mutate(polarity = dplyr::case_when(Ion_mode == "negative" ~ 0L, TRUE ~ 1L)) %>%
    dplyr::filter(.data$polarity %in% polarity_keep)

  empty_candi <- function(fdf) {
    fdf$candidate.id <- replicate(nrow(fdf), character(0), simplify = FALSE)
    fdf$candidate.formula <- replicate(nrow(fdf), character(0), simplify = FALSE)
    fdf$candidate.adduct <- replicate(nrow(fdf), character(0), simplify = FALSE)
    fdf$candidate.mz <- replicate(nrow(fdf), numeric(0), simplify = FALSE)
    fdf
  }

  if (!nrow(adducts) || !nrow(cpdbt)) {
    return(empty_candi(fdf))
  }

  cp.adduct <- MSCC::chemform_adduct(cpdbt$formula,
                                     adducts$adduct.formated,
                                     value = "all")
  cp.adduct <- cp.adduct %>%
    dplyr::mutate(compound_id = cpdbt$compound_id[id]) %>%
    dplyr::filter(findInterval(chemform.adduct.mz, mz_range) == 1)

  if (!nrow(cp.adduct)) {
    return(empty_candi(fdf))
  }

  matched.df <- match_mz_foverlaps(
    mz1 = fdf$mzmed,
    mz2 = cp.adduct$chemform.adduct.mz,
    ppm = mz.ppm
  )
  fdf$candidate.id <- lapply(seq_len(nrow(fdf)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$compound_id[as.numeric(idx)]
  })
  fdf$candidate.formula <- lapply(seq_len(nrow(fdf)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform[as.numeric(idx)]
  })
  fdf$candidate.adduct <- lapply(seq_len(nrow(fdf)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$adduct[as.numeric(idx)]
  })
  fdf$candidate.mz <- lapply(seq_len(nrow(fdf)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform.adduct.mz[as.numeric(idx)]
  })
  fdf
}


#' MS2 spectral similarity scores on a featureDefinitions-style data.frame
#' @param fdf data.frame with \code{feature_id}, \code{candidate.*}, \code{ms2_id}
#' @param cpdb CompoundDb object
#' @param sp.ms2 Spectra object (experimental MS2)
#' @param polarity integer polarity for filtering CompDb Spectra
#' @param ... unused
#' @return \code{fdf} with \code{score.ms2} list column
#' @keywords internal
fdf_get_ms2_score <- function(fdf,
                              cpdb,
                              sp.ms2,
                              polarity,
                              ...) {
  fdf <- as.data.frame(fdf, stringsAsFactors = FALSE)
  if (!nrow(fdf)) {
    return(fdf)
  }
  if (!"feature_id" %in% names(fdf)) {
    fdf$feature_id <- rownames(fdf)
  }
  if (!"ms2_id" %in% names(fdf)) {
    fdf$ms2_id <- replicate(nrow(fdf), character(0), simplify = FALSE)
  }
  if (!"candidate.id" %in% names(fdf)) {
    fdf$score.ms2 <- lapply(seq_len(nrow(fdf)), function(i) numeric(0))
    return(fdf)
  }

  if (!length(sp.ms2)) {
    fdf$score.ms2 <- lapply(fdf$candidate.id, function(x) rep(0, length(x)))
    return(fdf)
  }

  polarity <- as.integer(unique(polarity))
  Spectra_database <- Spectra::Spectra(cpdb)
  Spectra_database <- Spectra_set_MEM_backend(Spectra_database)
  Spectra_database <- ProtGenerics::filterPolarity(Spectra_database, polarity)

  Spectra_database <- Spectra_database %>%
    filterSpectra_below_PrecursorMz() %>%
    normalizeSpectra(norm_to = "max") %>%
    filterSpectraIntensity(ratio = 0.05) %>%
    Spectra::applyProcessing()
  sp.ms2 <- sp.ms2 %>%
    filterSpectra_below_PrecursorMz() %>%
    normalizeSpectra(norm_to = "max") %>%
    filterSpectraIntensity(ratio = 0.05) %>%
    Spectra_set_MEM_backend() %>%
    Spectra::applyProcessing()

  missing_ids <- setdiff(unlist(fdf$candidate.id), Spectra_database$compound_id)
  if (length(missing_ids)) {
    sp.empty <- makeEmptySpectra(compound_id = missing_ids)
    Spectra_database <- c(Spectra_database, sp.empty)
  }

  sp_names <- Spectra::spectraNames(sp.ms2)
  sp.exp <- lapply(seq_len(nrow(fdf)), function(i) {
    x <- fdf$ms2_id[[i]]
    if (is.null(x) || !length(x)) {
      return(NULL)
    }
    idx <- match(as.character(x), sp_names)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
      return(NULL)
    }
    sp.ms2[idx]
  })
  names(sp.exp) <- fdf$feature_id

  sp.split.list <- lapply(seq_len(nrow(fdf)), function(i) {
    i.candi <- fdf$candidate.id[[i]]
    i.adduct <- fdf$candidate.adduct[[i]]
    if (!length(i.candi)) {
      return(NULL)
    }
    i.df <- lapply(i.candi, function(x) which(Spectra_database$compound_id == x)) %>%
      `names<-`(fdf$candidate.id[[i]]) %>%
      unlist_to_df(name_to = "compound_id", value_to = "sp_id")
    i.df$adduct <- i.adduct[match(i.df$compound_id, i.candi)]
    i.df
  })
  sp.split.df <- data.table::rbindlist(sp.split.list, use.names = TRUE, idcol = "feature_id")
  if (!nrow(sp.split.df)) {
    fdf$score.ms2 <- lapply(fdf$candidate.id, function(x) rep(0, length(x)))
    return(fdf)
  }

  sp.ref <- Spectra_database[sp.split.df$sp_id]
  sp.ref$adduct <- sp.split.df$adduct
  sp.ref <- split(sp.ref, fdf$feature_id[sp.split.df$feature_id])

  .f <- function(expSpec, refSpec, ...) {
    if (is.null(expSpec)) {
      if (is.null(refSpec)) {
        return(NULL)
      }
      x <- unique(paste0(refSpec$compound_id, "_", refSpec$adduct))
      y <- rep(0, length(x))
      names(y) <- x
      return(y)
    }
    if (is.null(refSpec)) {
      return(NULL)
    }
    scorem <- Spectra::compareSpectra(expSpec, refSpec, , FUN = MsCoreUtils::ndotproduct, m = 2)
    dim(scorem) <- c(length(expSpec), length(refSpec))
    scorem[is.infinite(scorem) | is.na(scorem)] <- 0
    scores <- apply(scorem, 2, max, na.rm = TRUE)
    mean_f(scores, paste0(refSpec$compound_id, "_", refSpec$adduct))
  }

  fdf$score.ms2 <- BiocParallel::bplapply(
    seq_along(sp.exp),
    function(i) {
      fid <- fdf$feature_id[i]
      s <- .f(expSpec = sp.exp[[i]], refSpec = sp.ref[[fid]])
      unname(s[paste0(fdf$candidate.id[[i]], "_", fdf$candidate.adduct[[i]])])
    },
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
  )
  fdf
}


#' Final annotation pick on a featureDefinitions-style data.frame
#' @param fdf data.frame with candidate.* and score.ms2
#' @param cpdb CompoundDb object
#' @param cpdb.keys CompDb columns to attach
#' @param weight_mz,weight_ms2,weight_isopattern score weights
#' @param ... unused
#' @return annotated \code{fdf}
#' @keywords internal
fdf_get_feature_annotation <- function(fdf,
                                       cpdb,
                                       cpdb.keys = c("name", "formula", "smiles"),
                                       weight_mz = 0.1,
                                       weight_ms2 = 0.7,
                                       weight_isopattern = 0.2,
                                       ...) {
  fdf <- as.data.frame(fdf, stringsAsFactors = FALSE)
  if (!nrow(fdf)) {
    return(fdf)
  }
  if (!"feature_id" %in% names(fdf)) {
    fdf$feature_id <- rownames(fdf)
  }
  fdf$compound_id <- NA_character_
  fdf$adduct <- NA_character_
  fdf$score <- NA_real_
  fdf$mz_ref <- NA_real_

  if (!"score.isopattern" %in% names(fdf)) {
    fdf$score.isopattern <- lapply(fdf$candidate.id, function(x) rep(0, length(x)))
  }

  candi_idx <- which(lengths(fdf$candidate.id) != 0)
  if (length(candi_idx)) {
    xcms.candi.dt <- lapply(candi_idx, function(i) {
      iso <- fdf$score.isopattern[[i]]
      if (is.null(iso) || length(iso) != length(fdf$candidate.id[[i]])) {
        iso <- rep(0, length(fdf$candidate.id[[i]]))
      }
      ms2 <- fdf$score.ms2[[i]]
      if (is.null(ms2) || length(ms2) != length(fdf$candidate.id[[i]])) {
        ms2 <- rep(0, length(fdf$candidate.id[[i]]))
      }
      data.table::data.table(
        feature_id = fdf$feature_id[i],
        mz = fdf$mzmed[i],
        candidate.id = fdf$candidate.id[[i]],
        candidate.adduct = fdf$candidate.adduct[[i]],
        candidate.formula = fdf$candidate.formula[[i]],
        candidate.mz = fdf$candidate.mz[[i]],
        score.ms2 = ms2,
        score.isopattern = iso
      )
    }) %>% data.table::rbindlist()

    xcms.candi.dt <- xcms.candi.dt[
      , score.isopattern := ifelse(is.na(score.isopattern), 0, score.isopattern)][
      , score.isopattern := ifelse(is.nan(score.isopattern), 0, score.isopattern)][
      , score.mz := 1 - abs(candidate.mz - mz) / mz * 1e6 / 20][
      , score.mz := ifelse(score.mz < 0, 0, score.mz)][
      , score := score.isopattern * weight_isopattern +
          score.ms2 * weight_ms2 + score.mz * weight_mz]

    xcms.candi.dt.max <- xcms.candi.dt[, .SD[which.max(score)], by = feature_id]
    data.table::setnames(
      xcms.candi.dt.max,
      old = c("candidate.id", "candidate.adduct", "candidate.formula", "candidate.mz"),
      new = c("compound_id", "adduct", "formula", "mz_ref")
    )
    data.table::setkey(xcms.candi.dt.max, feature_id)

    xcms.candi <- xcms.candi.dt.max[fdf$feature_id]
    fdf$compound_id <- xcms.candi$compound_id
    fdf$adduct <- xcms.candi$adduct
    fdf$formula <- xcms.candi$formula
    fdf$score <- xcms.candi$score
    fdf$mz_ref <- xcms.candi$mz_ref
    ## keep per-candidate score lists for candidate.se expansion
    fdf$candidate.score.ms2 <- fdf$score.ms2
    fdf$candidate.score.isopattern <- fdf$score.isopattern
    ## best-candidate scalars for feature / metabolite export
    fdf$score.ms2 <- xcms.candi$score.ms2
    fdf$score.isopattern <- xcms.candi$score.isopattern
    fdf$score.mz <- xcms.candi$score.mz
  }

  dbinfo <- get_CompDb_info(cpdb, fdf$compound_id, keys = cpdb.keys)
  fdf[, colnames(dbinfo)] <- dbinfo
  fdf
}


get_xcms_feature_definitions <- function(xcms.xcms){
  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::select(
      !c(mzmin,mzmax,rtmin,rtmax,npeaks,peakidx)
    )
  return(xcms.fdf)

}


find_xcms_feature <- function(xcms.xcms,mz = 100,ppm = 10){

  fdf <- xcms::featureDefinitions(xcms.xcms)
  mzr <- mz.range.ppm(mz,ppm)
  fdf[between(fdf$mzmed,mzr[1],mzr[2] ), ]%>%
    as_tibble()
}


find_xcms_peaks <- function(xcms.xcms, mz = 100, ppm = 10) {

  peaks <- xcms::chromPeaks(xcms.xcms)
  mzr <- mz.range.ppm(mz, ppm)
  peaks[between(peaks[,"mz"], mzr[1], mzr[2]), ]
}


#' @describeIn xcms_extension_plot plot peaks distribution
#' @description export peaks data by xcms::chromPeaks and plot by ggplot2
#'
#' @param xcms.xcms XCMSnExp object
#' @param plot.title title
#' @param type `"o"`, for geom_point, `"l"`, for geom_segment
#'
#' @return ggplot object
#' @export
#'

plot_xcms_peaks_distribution <- function(xcms.xcms,plot.title = "Peaks distribution",type = "o"){

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::mutate(as.data.frame(xcms::chromPeakData(xcms.xcms)),
                  peak_id = rownames(.),
                  merged = grepl(peak_id,pattern = "CPM"))%>%
    dplyr::filter(!is.na(maxo),
                  rtmax-rtmin <60,
                  !merged)
  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( xcms::processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  if (type == "o") {
    ggplot(xcms.peaks)+
      geom_point(aes(x = rt,y=mz,
                     col = log10(maxo),
                     alpha = log10(maxo)/10,
                     size = (rtmax-rtmin)),
      )+
      scale_size_area(max_size = 8)+
      labs(title = plot.title,
           subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                             "; SN = ",xcms.findpeak.param@snthresh,
                             "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
           col = "Log10\n(Intensity)",
           size = "Peak width",
           x = "Retention time",
           y = "mz")+
      guides(alpha = "none")+
      scale_color_gradientn(breaks = c(0,3,6,9),
                            labels = c(0,3,6,9),
                            limits = c(0,9),
                            values = c(0,2,4,7,9)/9,
                            colors = c("white","white","yellow","red","red"))+
      theme_bw()+
      theme(text = element_text(size = 8))->peaks.dis.plot

  }else if(type == "l"){
    ggplot(xcms.peaks)+
      geom_segment(aes(x = rtmin , xend = rtmax , y = mz, yend = mz,col = log10(maxo)),
                   size = 0.6
      )+
      labs(title = plot.title,
           subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                             "; SN = ",xcms.findpeak.param@snthresh,
                             "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
           col = "Log10(Intensity)",
           x = "Retention time",
           y = "mz")+
      guides(alpha = "none")+
      scale_color_gradientn(breaks = c(0,3,6,9),
                            labels = c(0,3,6,9),
                            limits = c(0,9),
                            values = c(0,2,4,7,9)/9,
                            colors = c("white","white","yellow","red","red"))+
      theme_bw()+
      theme(text = element_text(size = 8))->peaks.dis.plot
    peaks.dis.plot


  }
  return(peaks.dis.plot)


}



#' @describeIn xcms_extension_plot plot features distribution
#' @description Visualizes the distribution of detected features in a 2D space of retention time (x-axis) vs m/z (y-axis). Point size represents peak width, color represents log10 intensity. Includes peak detection parameters in subtitle.
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param plot.title Character title for the plot (default "Features distribution").
#'
#' @return ggplot object.
#' @export
#'

plot_xcms_features_distribution <-
  function(xcms.xcms, plot.title = "Features distribution") {
    xcms.features <- xcms::featureDefinitions(xcms.xcms) %>%
      as.data.frame() %>%
      mutate(mz = mzmed, rt = rtmed)
    xcms.features$maxo <-
      apply(xcms::featureValues(xcms.xcms, value = "maxo"), 1, median, na.rm = T)

    xcms.process.type <-
      xcms::processHistory(xcms.xcms) %>% sapply(xcms::processType)
    xcms.findpeak.param <-
      xcms:: processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]] %>%
      xcms::processParam()
    ### generate scale
    maxo.range <- quantile(log10((xcms.features$maxo)))
    peak.witdh.range <- quantile(xcms.features$peakWidth)

    ggplot(xcms.features) +
      geom_point(aes(
        x = rt,
        y = mz,
        col = log10(maxo),
        alpha = log10(maxo) / 10,
        size =peakWidth
      ),) +
      scale_size_continuous(breaks  = peak.witdh.range[2:4],
                            labels = round(peak.witdh.range[2:4]),
                            range = c(0,10)) +
      xlim(c(0, 800)) +
      labs(
        title = plot.title,
        subtitle = paste0(
          "ppm = ",
          xcms.findpeak.param@ppm,
          "; SN = ",
          xcms.findpeak.param@snthresh,
          "; prefilter = (",
          paste0(xcms.findpeak.param@prefilter, collapse = ","),
          ")"
        ),
        col = "Log10(Intensity)",
        size = "Peak width",
        x = "Retention time",
        y = "mz"
      ) +
      guides(alpha = "none") +
      scale_color_gradientn(breaks = c(0,3,6,9),
                            labels = c(0,3,6,9),
                            limits = c(0,9),
                            values = c(0,2,4,7,9)/9,
                            colors = c("white","white","yellow","red","red"))+
      theme_bw()+
      theme(text = element_text(size = 8)) -> peaks.dis.plot
    peaks.dis.plot
    return(peaks.dis.plot)





  }


xcms_remove_feature_var <- function(xcms.xcms,var){

  xcms.fdf<- xcms::featureDefinitions(xcms.xcms)
  var.selected <- setdiff(colnames(xcms.fdf),var)
  xcms.fdf <- xcms.fdf[,var.selected]
  xcms.fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)
}


#' @description extract Chromatogram from xcms according to feature's mz range and plot
#' @describeIn xcms_extension_plot plot feature chromatogram
#' @param xcms.xcms XCMSnExp object
#' @param feature.id feature id
#' @param sampleNames sample names to include
#'
#' @return ggplot object
#' @export
#'

plot_xcms_feature_chromatogram <- function(xcms.xcms ,feature.id, sampleNames =NULL ){

  ### select samples
  all.sample.names <- Biobase::sampleNames(xcms.xcms)
  xcms.sample.info <- as.data.frame(Biobase::pData(xcms.xcms), stringsAsFactors = FALSE)
  if (!"sampleNames" %in% colnames(xcms.sample.info)) {
    xcms.sample.info$sampleNames <- all.sample.names
  }
  rownames(xcms.sample.info) <- all.sample.names

  if (is.null(sampleNames)) {
    sampleNames <- all.sample.names
  }
  sampleNames <- as.character(sampleNames)
  sampleNames <- intersect(sampleNames, all.sample.names)
  if (!length(sampleNames)) {
    stop("No valid sampleNames found in xcms.xcms.")
  }

  xcms.sample.info <- xcms.sample.info[sampleNames, , drop = FALSE]
  if (nrow(xcms.sample.info) > 5) {
    if ("group" %in% colnames(xcms.sample.info) && any(!is.na(xcms.sample.info$group))) {
      xcms.sample.info.sub <- xcms.sample.info %>%
        dplyr::group_by(group) %>%
        dplyr::slice_sample(n = 1) %>%
        as.data.frame(stringsAsFactors = FALSE)
    } else {
      xcms.sample.info.sub <- xcms.sample.info[seq_len(5), , drop = FALSE]
    }
  } else {
    xcms.sample.info.sub <- xcms.sample.info
  }

  sample.idx <- which(all.sample.names %in% xcms.sample.info.sub$sampleNames)
  if (!length(sample.idx)) {
    stop("No samples selected for chromatogram extraction.")
  }
  xcms.sub <- .xcms_filter_file(xcms.xcms, sample.idx)

  ### mz / rt from feature peaks
  xcms.fdef <- xcms::featureDefinitions(xcms.xcms)
  if (is.numeric(feature.id)) {
    feature.id <- rownames(xcms.fdef)[feature.id]
  }
  feature.id <- as.character(feature.id)[1]
  if (is.na(feature.id) || !nzchar(feature.id) || !(feature.id %in% rownames(xcms.fdef))) {
    stop("feature.id does not exist in featureDefinitions(xcms.xcms).")
  }

  xcms.feature <- xcms.fdef[feature.id, , drop = FALSE]
  peak.idx <- xcms.feature$peakidx[[1]]
  if (is.null(peak.idx) || !length(peak.idx)) {
    stop("Selected feature has no linked peaks (empty peakidx).")
  }
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)[peak.idx, , drop = FALSE]
  mz.range <- as.numeric(c(min(xcms.peaks[, "mzmin"], na.rm = TRUE),
                           max(xcms.peaks[, "mzmax"], na.rm = TRUE)))
  rt.range <- as.numeric(c(min(xcms.peaks[, "rtmin"], na.rm = TRUE),
                           max(xcms.peaks[, "rtmax"], na.rm = TRUE)))
  if (length(mz.range) != 2 || length(rt.range) != 2 ||
      any(!is.finite(mz.range)) || any(!is.finite(rt.range))) {
    stop("Failed to derive finite mz/rt ranges from feature peaks.")
  }

  xcms.chrom <- xcms::chromatogram(xcms.sub,
                                   mz = mz.range,
                                   rt = rt.range,
                                   BPPARAM = BiocParallel::SerialParam())
  if (!methods::is(xcms.chrom, "XChromatograms")) {
    stop("Chromatogram extraction failed for selected feature and samples.")
  }

  xcms.chrom.data <- get_chroms_data(xcms.chrom)%>%
    dplyr::mutate(group = Biobase::sampleNames(xcms.sub)[col])

  rt_all <- if (inherits(xcms.sub, "MsExperiment") || inherits(xcms.sub, "XcmsExperiment")) {
    as.numeric(Spectra::rtime(ProtGenerics::spectra(xcms.sub)))
  } else {
    as.numeric(MSnbase::rtime(xcms.sub))
  }
  ggplot(xcms.chrom.data)+
    geom_line(aes(x = rt,y = intensity , col = group))+
    xlim(c(min(rt_all, na.rm = TRUE), max(rt_all, na.rm = TRUE)))+
    labs(col = "",x = "Retention time", y = "Intensity",
         title = paste0(feature.id),
         subtitle = paste0( "mz: ",round(mz.range[1],4),
                            " ~ ",round(mz.range[2],4),
                            "\nrt: ",round(rt.range[1],2),
                            " ~ ",round(rt.range[2],2) ))+
    theme_bw()+
    theme(text = element_text(size = 8))


}



plot_xcms_peaks_mzerror_density <- function(xcms.xcms,
                                            plot.title = "Peak mz error distribution"){

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    mutate(ppm = (mzmax-mzmin)/mz*1e6,
           mz_diff = mzmax-mzmin,
           peak_width = rtmax-rtmin)

  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( xcms::processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  ggplot(xcms.peaks,aes(x = mz , y = ppm)) +
    stat_density_2d(aes(fill= after_stat(level)),
                    contour = T,
                    geom = "polygon",bins = 100)+
    geom_point(size = 0.1,alpha = 0.1)+
    scale_fill_gradient(low="#00000001",high = "red")+
    scale_x_continuous(expand = c(0.1,0.1))+
    scale_y_continuous(expand = c(0.1,0.1))+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" )
         )+
    guides(fill = "none")+
    theme_bw()+
    theme(text = element_text(size = 8))

}



#' @description plot scans number of MS1 levels in each peak, note that to many peaks will lead to stuck,
#' apply `filterFile` to decrease peaks count
#' @describeIn xcms_extension_plot plot MS1 scan counts for peaks
#' @param xcms.xcms XCMSnExp object should be a `XCMSnExp` object after `findChromPeaks`
#' @param plot.title title
#'
#' @return ggplot object
#' @export
#'

plot_xcms_peaks_ms1_scans <- function(xcms.xcms,plot.title = "Peaks Sans of MS1"){

  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( xcms::processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()
  xcms.scans <- get_xcms_scan_Stat(xcms.xcms)%>%
    dplyr::filter(msLevel== 1)
  peaks_scans <- function(x,xcms.scans){
    sum(x["rtmax"] > xcms.scans$rtime  & x["rtmin"] < xcms.scans$rtime )

  }
  xcms.peaks$ms1_scans_no <- apply(xcms.peaks ,1,peaks_scans , xcms.scans)
  ggplot(xcms.peaks)+
    geom_segment(aes(x = rtmin , xend = rtmax , y = ms1_scans_no, yend = ms1_scans_no,col = log10(maxo)),
                 size = 0.6
    )+
    geom_hline(yintercept = 7)+
    geom_boxplot(aes( x = max(rt)*1.2 , y =ms1_scans_no),width = diff(range(xcms.peaks$rt))*0.1)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
         col = "Log10(Intensity)",
         x = "Retention time",
         y = "Scan count of MS1 in each peak")+
    scale_y_log10(breaks = c(1,2,3,4,5,6,7,8,10,20))+
    guides(alpha = "none")+
    scale_color_gradientn(breaks = c(0,3,6,9),
                          labels = c(0,3,6,9),
                          limits = c(0,log10(median(xcms.peaks$maxo))*1.5),
                          values = c(0,2,4,7,9)/9,
                          na.value = "red",
                          colors = c("white","white","yellow","red","red"))+
    theme_bw()+
    theme(text = element_text(size = 8))->peaks.dis.plot
  peaks.dis.plot

}



#' @description Dot-plot MS1 scan frequency along retention time. Scans are
#' counted in successive \code{rt_window}-wide RT bins (per file); frequency is
#' \code{scan_count / rt_window}.
#' @describeIn xcms_extension_plot plot MS1 scan frequency vs RT
#' @param xcms XCMSnExp / XcmsExperiment object
#' @param rt_window positive numeric; RT window width (same unit as retention
#'   time, typically seconds)
#' @param plot.title title
#'
#' @return ggplot object
#' @export
#'
plot_xcms_ms1_scan_freq <- function(xcms, rt_window = 5,
                                    plot.title = "MS1 Scan Frequency") {

  if (!is.numeric(rt_window) || length(rt_window) != 1L ||
      !is.finite(rt_window) || rt_window <= 0) {
    stop("'rt_window' must be a single positive finite number.", call. = FALSE)
  }

  xcms.scans <- get_xcms_scan_Stat(xcms) %>%
    dplyr::filter(as.integer(msLevel) == 1L)

  if (!nrow(xcms.scans)) {
    stop("No MS1 scans found in 'xcms'.", call. = FALSE)
  }

  ## MS1 only: count MS1 scans per RT window (MS2 excluded)
  freq_df <- xcms.scans %>%
    dplyr::group_by(fileIdx) %>%
    dplyr::mutate(
      rt0 = min(rtime, na.rm = TRUE),
      rt_bin = floor((rtime - rt0) / rt_window)
    ) %>%
    dplyr::group_by(fileIdx, rt_bin, rt0) %>%
    dplyr::summarise(
      scan_count = dplyr::n(),
      rt = dplyr::first(rt0) + (dplyr::first(rt_bin) + 0.5) * rt_window,
      .groups = "drop"
    ) %>%
    dplyr::mutate(freq = scan_count / rt_window)

  rt_span <- diff(range(freq_df$rt, na.rm = TRUE))
  if (!is.finite(rt_span) || rt_span <= 0) {
    rt_span <- rt_window
  }

  ggplot(freq_df) +
    geom_point(aes(x = rt, y = freq, col = freq), size = 0.6) +
    geom_boxplot(aes(x = max(rt) * 1.2, y = freq),
                 width = rt_span * 0.1) +
    labs(title = plot.title,
         subtitle = paste0("rt_window = ", rt_window, " (MS1 only)"),
         col = "MS1 freq",
         x = "Retention time",
         y = "MS1 scan frequency (count / rt_window)") +
    guides(alpha = "none") +
    scale_color_gradientn(
      breaks = c(0, 1, 2, 3, 4, 5),
      labels = c(0, 1, 2, 3, 4, 5),
      limits = c(0, 5),
      colors = ggsci::pal_bs5("indigo")(6),
      na.value = "#666666"
    ) +
    theme_bw() +
    theme(text = element_text(size = 8))
}



#' @description Visualizes the number of MS2 scans that overlap each chromatographic peak based on retention time and m/z ranges. Produces a scatter plot with jitter, violin distribution, and counts of peaks with 0-5 MS2 scans.
#' @param xcms.xcms XCMSnExp object with detected peaks and MS2 scans.
#' @param plot.title Character title for the plot (default "Peaks Sans of MS2").
#'
#' @describeIn xcms_extension_plot plot MS2 scan counts for peaks
#' @return ggplot object.
#' @export
#'

plot_xcms_peaks_ms2_scans <- function(xcms.xcms,plot.title = "Peaks Sans of MS2"){

  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( xcms::processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()
  xcms.scans <- get_xcms_scan_Stat(xcms.xcms)%>%
    dplyr::filter(msLevel== 2)
  peaks_scans <- function(x,xcms.scans){
    sum(x["rtmax"] > xcms.scans$rtime  & x["rtmin"] < xcms.scans$rtime&
          x["mzmax"] > xcms.scans$precursorMz  & x["mzmin"] < xcms.scans$precursorMz)

  }
  xcms.peaks$ms2_scans_no <- apply(xcms.peaks ,1,peaks_scans , xcms.scans)
  ms2_scans_table <- table(xcms.peaks$ms2_scans_no)
  ggplot(xcms.peaks)+
    geom_jitter(aes(x = rt, y = ms2_scans_no, col = log10(maxo)),
                 size = 0.6
    )+
    #geom_hline(yintercept = 7)+
    geom_violin(aes( x = max(rt)*1.2 , y =ms2_scans_no),width = diff(range(xcms.peaks$rt))*0.1)+
    geom_text(aes(x =  max(rt)*1.3,y = 0,label = ms2_scans_table["0"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 1,label = ms2_scans_table["1"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 2,label = ms2_scans_table["2"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 3,label = ms2_scans_table["3"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 4,label = ms2_scans_table["4"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 5,label = ms2_scans_table["5"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.4,y = 5,label = ""),size = 2.67,)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ,"\n",
                           sum(xcms.peaks$ms2_scans_no > 0)," / ",length(xcms.peaks$ms2_scans_no),
                           " ( ",sprintf("%.2f",sum(xcms.peaks$ms2_scans_no > 0)/length(xcms.peaks$ms2_scans_no)*100),"% )"),
         col = "Log10(Intensity)",
         x = "Retention time",
         y = "Scan count of MS2 in each peak")+
    guides(alpha = "none")+
    scale_color_gradientn(breaks = c(0,3,6,9),
                          labels = c(0,3,6,9),
                          limits = c(0,log10(median(xcms.peaks$maxo))*1.5),
                          values = c(0,2,4,7,9)/9,
                          na.value = "red",
                          colors = c("white","white","yellow","red","red"))+
    theme_bw()+
    theme(text = element_text(size = 8))->peaks.dis.plot
  peaks.dis.plot

}



plot_xcms_ms2_distribution <- function(xcms.xcms,plot.title = "MS2 Precursor distribution" ){

 scan.data <- get_xcms_scan_Stat(xcms.xcms)%>%
   dplyr::filter(msLevel==2)

 ms1.rt <- get_xcms_scan_Stat(xcms.xcms)%>%
   dplyr::filter(msLevel==1)%>%
   dplyr::pull(rtime)

 ggplot(scan.data)+
   #geom_vline(xintercept = ms1.rt,linewidth = 0.05,col = "black")+
   geom_point(aes(x = rtime,y= precursorMz,
                  col = log10(precursorIntensity)),
   )+
   labs(title = plot.title,
        col = "Log10\n(Intensity)",
        size = "Peak width",
        x = "Retention time",
        y = "mz")+
   guides(alpha = "none")+
   scale_color_gradientn(breaks = c(0,3,6,9),
                         labels = c(0,3,6,9),
                         limits = c(0,9),
                         values = c(0,2,4,7,9)/9,
                         colors = c("white","white","yellow","red","red"))+
   theme_bw()+
   theme(text = element_text(size = 8))->peaks.dis.plot

 open_plot_win(peaks.dis.plot,width = 25,height = 5)
 peaks.dis.plot
}



plot_xcms_peaks_SN_distribution <- function(xcms.xcms,plot.title = "Peaks SNR(Signal to Noise Ratio)"){


  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( xcms::processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()

  ggplot(xcms.peaks)+
    geom_jitter(aes(x = rt, y = log10(sn), col = log10(maxo)),
                size = 0.6
    )+
    #geom_hline(yintercept = 7)+
    geom_violin(aes( x = max(rt)*1.2 , y =log10(sn)),width = diff(range(xcms.peaks$rt))*0.1)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
         col = "Log10(Intensity)",
         x = "Retention time",
         y = "log10(SNR)")+
    guides(alpha = "none")+
    scale_color_gradientn(breaks = c(0,3,6,9),
                          labels = c(0,3,6,9),
                          limits = c(0,log10(median(xcms.peaks$maxo))*1.5),
                          values = c(0,2,4,7,9)/9,
                          na.value = "red",
                          colors = c("white","white","yellow","red","red"))+
    theme_bw()+
    theme(text = element_text(size = 8))->peaks.dis.plot
  peaks.dis.plot



}

#' @description extract EIC according to peaks' mzrange and rtrange,
#' note that if multiple sample in xcms object, only first sample will be extracted
#' @describeIn xcms_extension_plot plot chromatogram for a peak
#' @param xcms.xcms XCMSnExp object
#' @param peak_id peak id
#' @param rt expansion range for rt
#'
#' @return ggplot object
#' @export
#'

plot_xcms_peaks_Chromatogram <- function(xcms.xcms,peak_id,rt = "expand"){

  peaks.data <- xcms::chromPeaks(xcms.xcms)[peak_id,,drop = F]
  peak_id <- rownames(peaks.data)
  mz.range <- c(peaks.data[,c("mzmin","mzmax")])
  rt.range <- c(peaks.data[,c("rtmin","rtmax")])
  xcms.chrom <- get_xcms_peaks_chromatogram(xcms.xcms,
                                     peaks.id = peak_id,
                                     rt.range = rt)
  chrom.data <- get_chroms_data(xcms.chrom)%>%
    dplyr::mutate(fill = rt > min(rt.range)&rt <max(rt.range),
                  sample = Biobase::sampleNames(xcms.xcms)[col]
                 )%>%
    dplyr::filter(!is.na(intensity))

  ggplot(chrom.data)+
    geom_line(aes(x = rt,y = intensity,col =sample),linetype = 1)+
    geom_area(aes(x = rt,y = intensity, fill = sample),
              stat = "identity",alpha = 0.1)+
    #geom_point(aes(x = rt, y = fit))+
    scale_fill_manual(values = c("FALSE" = "transparent","TRUE" = "grey"))+
    labs(title = paste0(peak_id),
         subtitle = paste0("mz:",paste0(sprintf("%.4f",range(mz.range)),collapse = " - "), ";     ",
                           "rt:",paste0(sprintf("%.2f",range(rt.range)),collapse = " - "),"\n",
                           "mz range = ",sprintf("%.2f",mean(diff(range(mz.range))/mz.range)*1e6)," ppm;     ",
                           "peak width = ", sprintf("%.2f",diff(range(rt.range))),"\n"
                          # "shape score = ",get_chrom_peaks_shape_score(xcms.chrom[1,1])
                           ),
         x = "Retention time")+
    guides(fill = "none")+
    theme_bw()+
    theme(text = element_text(size = 8))



}



chromPeaks_Sta <- function(xcms.xcms){

  xcms.peaks.info <- xcms::chromPeaks(xcms.xcms)
  xcms.peaks.xchroms <- get_xcms_peaks_chromatogram(xcms.xcms,
                                             1:nrow(xcms.peaks.info))




}


#' Build per-spectrum scan table for MsExperiment / XcmsExperiment
#' @noRd
.get_xcms_scan_table <- function(xcms.xcms) {
  if (inherits(xcms.xcms, "MsExperiment") || inherits(xcms.xcms, "XcmsExperiment")) {
    sp <- ProtGenerics::spectra(xcms.xcms)
    if (!length(sp)) {
      return(data.frame(
        fileIdx = integer(),
        msLevel = integer(),
        rtime = numeric(),
        polarity = integer(),
        stringsAsFactors = FALSE
      ))
    }
    sd <- as.data.frame(Spectra::spectraData(sp), stringsAsFactors = FALSE)
    origins <- as.character(sd$dataOrigin)
    if (!length(origins) || all(is.na(origins))) {
      origins <- as.character(sd$dataStorage)
    }
    file_levels <- unique(origins)
    sd$fileIdx <- as.integer(match(origins, file_levels))
    if (!"rtime" %in% names(sd)) {
      sd$rtime <- as.numeric(Spectra::rtime(sp))
    }
    if (!"msLevel" %in% names(sd)) {
      sd$msLevel <- as.integer(Spectra::msLevel(sp))
    }
    if (!"polarity" %in% names(sd)) {
      sd$polarity <- as.integer(Spectra::polarity(sp))
    }
    sd$msLevel <- as.integer(sd$msLevel)
    sd$polarity <- as.integer(sd$polarity)
    sd$rtime <- as.numeric(sd$rtime)
    return(sd)
  }
  stop(
    "OnDiskMSnExp / XCMSnExp scan tables are no longer supported; ",
    "use MsExperiment / XcmsExperiment (Spectra columns: precursorMz, rtime).",
    call. = FALSE
  )
}

#' @title Build xcms centWave roiList from mz/rt targets
#' @description
#' Construct a \code{roiList} accepted by \code{xcms::CentWaveParam(roiList = ...)}.
#' Input a matrix/data.frame with columns \code{mz} and \code{rt} (seconds). For each
#' target, \code{mzmin/mzmax} are calculated using ppm tolerance and the RT window
#' \code{rtmin/rtmax} is mapped to scan indices.
#'
#' Scan indices are computed per file from MS1 spectra, then expanded to the
#' union range across files and \strong{clamped} to
#' \code{[1, min(n_MS1)]} so one shared \code{roiList} stays valid for every
#' file in \code{findChromPeaks()}. Unclamped unions can overflow shorter runs
#' and raise \code{Error in scanrange} inside \code{.centWave_orig}.
#'
#' Note: \code{scmin}/\code{scmax} must be finite integer scan indices
#' (not \code{c(0, Inf)}). \code{centWave} uses them in arithmetic
#' (\code{N <- scmax - scmin + 1}) before clipping to each file's
#' \code{length(scantime)}.
#'
#' @param mzrt matrix/data.frame with columns \code{mz} and \code{rt}.
#' @param xcms.xcms \code{XCMSnExp}, \code{MsExperiment}, or \code{XcmsExperiment}
#'   used to map RT to scan indices.
#' @param ppm numeric, ppm tolerance for mz window.
#' @param rt_tol numeric, RT tolerance in seconds.
#' @param ion_mode optional integer 1 (positive) or 0 (negative). If NULL, inferred
#'   from scan polarity; must be unique.
#'
#' @return list of ROI objects (each ROI is a list with \code{scmin, scmax, mzmin, mzmax, length, intensity}).
#' @export
get_xcms_roi_list <- function(mzrt,
                              xcms.xcms,
                              ppm = 10,
                              rt_tol = 30,
                              ion_mode = NULL) {

  if (is.null(mzrt) || length(mzrt) == 0) return(list())
  mzrt <- as.data.frame(mzrt, stringsAsFactors = FALSE)
  if (!all(c("mz", "rt") %in% colnames(mzrt))) {
    stop("mzrt must have columns: mz, rt")
  }
  mzrt$mz <- suppressWarnings(as.numeric(mzrt$mz))
  mzrt$rt <- suppressWarnings(as.numeric(mzrt$rt))
  mzrt <- mzrt[is.finite(mzrt$mz) & is.finite(mzrt$rt), , drop = FALSE]
  if (!nrow(mzrt)) return(list())

  fdat <- .get_xcms_scan_table(xcms.xcms)
  if (is.null(fdat) || nrow(fdat) == 0) {
    stop("xcms.xcms has empty scan table; cannot derive scan indices for roiList.")
  }
  need_cols <- c("fileIdx", "msLevel", "rtime", "polarity")
  if (!all(need_cols %in% colnames(fdat))) {
    stop("xcms.xcms scan table must contain columns: ", paste(need_cols, collapse = ", "), ".")
  }

  if (is.null(ion_mode)) {
    ion_mode <- unique(as.integer(fdat$polarity[as.integer(fdat$msLevel) == 1L]))
    ion_mode <- ion_mode[!is.na(ion_mode)]
    if (length(ion_mode) != 1) {
      stop("Cannot infer ion_mode (multiple polarities in MS1 scans). Provide ion_mode = 0/1.")
    }
  }
  ion_mode <- as.integer(ion_mode)

  files <- seq_len(.xcms_nfiles(xcms.xcms))
  ms1_rt_by_file <- lapply(files, function(fi) {
    idx <- which(
      as.integer(fdat$fileIdx) == as.integer(fi) &
        as.integer(fdat$msLevel) == 1L &
        as.integer(fdat$polarity) == ion_mode
    )
    if (!length(idx)) return(NULL)
    rt <- as.numeric(fdat$rtime[idx])
    rt[order(rt)]
  })
  n_ms1 <- vapply(ms1_rt_by_file, function(x) {
    if (is.null(x)) 0L else length(x)
  }, integer(1))
  if (!any(n_ms1 > 0L)) {
    stop("No MS1 spectra found for ion_mode=", ion_mode, " while building roiList.")
  }
  min_n_ms1 <- min(n_ms1[n_ms1 > 0L])
  if (length(unique(n_ms1[n_ms1 > 0L])) > 1L) {
    message(sprintf(
      "get_xcms_roi_list: MS1 scan counts differ across files (min=%d, max=%d); clamping ROI scmax to min.",
      min_n_ms1, max(n_ms1)
    ))
  }

  roi_list <- vector("list", length = nrow(mzrt))
  kept <- logical(nrow(mzrt))

  for (i in seq_len(nrow(mzrt))) {
    mz <- mzrt$mz[[i]]
    rt <- mzrt$rt[[i]]
    mzmin <- mz - mz * ppm / 1e6
    mzmax <- mz + mz * ppm / 1e6
    rtmin <- rt - rt_tol
    rtmax <- rt + rt_tol

    scmins <- integer()
    scmaxs <- integer()
    for (fi in files) {
      rts <- ms1_rt_by_file[[fi]]
      if (is.null(rts)) next
      in_rt <- which(rts >= rtmin & rts <= rtmax)
      if (!length(in_rt)) next
      scmins <- c(scmins, min(in_rt))
      scmaxs <- c(scmaxs, max(in_rt))
    }

    if (!length(scmins)) next
    scmin <- max(1L, as.integer(min(scmins)))
    scmax <- min(as.integer(max(scmaxs)), as.integer(min_n_ms1))
    if (!is.finite(scmin) || !is.finite(scmax) || scmax < scmin) next
    roi_list[[i]] <- list(
      scmin = as.integer(scmin),
      scmax = as.integer(scmax),
      mzmin = as.numeric(mzmin),
      mzmax = as.numeric(mzmax),
      length = as.integer(scmax - scmin + 1L),
      intensity = 0
    )
    kept[[i]] <- TRUE
  }

  roi_list[kept]
}


xcms_filter_peaks_NA <- function(xcms.xcms, verbose = TRUE) {

  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) {
      message("xcms_filter_peaks_NA: no chromPeaks to filter")
    }
    return(xcms.xcms)
  }

  mz_col <- NULL
  if ("mz" %in% colnames(pks)) mz_col <- "mz"
  if (is.null(mz_col) && "mzmed" %in% colnames(pks)) mz_col <- "mzmed"
  if (is.null(mz_col)) {
    if (isTRUE(verbose)) {
      message("xcms_filter_peaks_NA: no mz column found in chromPeaks")
    }
    return(xcms.xcms)
  }

  mzv <- as.numeric(pks[, mz_col])
  bad <- is.na(mzv) | is.nan(mzv)
  n_bad <- sum(bad)
  n_total <- nrow(pks)
  ratio <- if (n_total == 0) 0 else n_bad / n_total

  if (isTRUE(verbose)) {
    message(sprintf(
      "xcms_filter_peaks_NA: %d/%d (%.4f) chromPeaks have mz NA/NaN",
      n_bad, n_total, ratio
    ))
  }

  if (n_bad > 0) {
    pks2 <- pks[!bad, , drop = FALSE]
    xcms::chromPeaks(xcms.xcms) <- pks2
    # Keep chromPeakData synchronized with chromPeaks for XcmsExperiment.
    # The chromPeaks<- setter updates only @chromPeaks, which can desync
    # row counts and later trigger out-of-bounds indexing in xcms internals.
    cpd <- tryCatch(xcms::chromPeakData(xcms.xcms), error = function(e) NULL)
    if (!is.null(cpd) && nrow(cpd) > 0 && nrow(cpd) != nrow(pks2)) {
      rn <- rownames(pks2)
      if (!is.null(rn) && length(rn) == nrow(pks2)) {
        idx <- match(rn, rownames(cpd))
        if (!anyNA(idx)) {
          xcms::chromPeakData(xcms.xcms) <- as.data.frame(cpd[idx, , drop = FALSE])
        } else if (nrow(cpd) >= nrow(pks2)) {
          xcms::chromPeakData(xcms.xcms) <- as.data.frame(cpd[seq_len(nrow(pks2)), , drop = FALSE])
        }
      } else if (nrow(cpd) >= nrow(pks2)) {
        xcms::chromPeakData(xcms.xcms) <- as.data.frame(cpd[seq_len(nrow(pks2)), , drop = FALSE])
      }
    }
  }

  return(xcms.xcms)
}

filter_xcms_chromPeaks_mz_width <- function(xcms.xcms, ppm = 20, verbose = TRUE) {
  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) message("filter_xcms_chromPeaks_mz_width: no chromPeaks to filter")
    return(xcms.xcms)
  }
  need <- c("mz", "mzmin", "mzmax")
  if (!all(need %in% colnames(pks))) {
    if (isTRUE(verbose)) {
      message("filter_xcms_chromPeaks_mz_width: required columns not found (mz/mzmin/mzmax)")
    }
    return(xcms.xcms)
  }

  mz <- suppressWarnings(as.numeric(pks[, "mz"]))
  mzmin <- suppressWarnings(as.numeric(pks[, "mzmin"]))
  mzmax <- suppressWarnings(as.numeric(pks[, "mzmax"]))
  mz_width_ppm <- (mzmax - mzmin) / pmax(mz, 1e-12) * 1e6
  bad <- !is.finite(mz_width_ppm) | (mz_width_ppm > ppm)
  n_bad <- sum(bad, na.rm = TRUE)
  n_total <- nrow(pks)

  if (isTRUE(verbose)) {
    message(sprintf(
      "filter_xcms_chromPeaks_mz_width: %d/%d chromPeaks removed (mz width > %.2f ppm); this is aggressive and may remove real peaks",
      n_bad, n_total, ppm
    ))
  }
  if (n_bad > 0) {
    xcms::chromPeaks(xcms.xcms) <- pks[!bad, , drop = FALSE]
  }
  xcms.xcms
}

#' @title Fix overly wide xcms chromPeaks mz window
#' @description
#' For chromPeaks with abnormal m/z width larger than \code{ppm}, this function
#' does not remove peaks. Instead, it recalculates and replaces \code{mzmin} and
#' \code{mzmax} around peak center \code{mz} so final width equals the target
#' ppm window.
#'
#' @param xcms.xcms XCMSnExp object.
#' @param ppm numeric ppm threshold/target width, default \code{20}.
#' @param verbose logical, print summary message.
#'
#' @return XCMSnExp object with fixed \code{chromPeaks} m/z windows.
#' @export
fix_xcms_chromPeaks_mz_width <- function(xcms.xcms, ppm = 20, verbose = TRUE) {
  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) message("fix_xcms_chromPeaks_mz_width: no chromPeaks to fix")
    return(xcms.xcms)
  }
  need <- c("mz", "mzmin", "mzmax")
  if (!all(need %in% colnames(pks))) {
    if (isTRUE(verbose)) {
      message("fix_xcms_chromPeaks_mz_width: required columns not found (mz/mzmin/mzmax)")
    }
    return(xcms.xcms)
  }

  mz <- suppressWarnings(as.numeric(pks[, "mz"]))
  mzmin <- suppressWarnings(as.numeric(pks[, "mzmin"]))
  mzmax <- suppressWarnings(as.numeric(pks[, "mzmax"]))
  bad_mz <- !is.finite(mz) | mz <= 0
  mz[bad_mz] <- (mzmin[bad_mz] + mzmax[bad_mz]) / 2
  bad_mz <- !is.finite(mz) | mz <= 0
  mz_width_ppm <- (mzmax - mzmin) / pmax(mz, 1e-12) * 1e6
  to_fix <- !bad_mz & is.finite(mz_width_ppm) & (mz_width_ppm > ppm)
  n_fix <- sum(to_fix, na.rm = TRUE)
  n_total <- nrow(pks)

  if (n_fix > 0) {
    half <- mz[to_fix] * as.numeric(ppm) / 2e6
    pks[to_fix, "mzmin"] <- mz[to_fix] - half
    pks[to_fix, "mzmax"] <- mz[to_fix] + half
    xcms::chromPeaks(xcms.xcms) <- pks
  }
  if (isTRUE(verbose)) {
    message(sprintf(
      "fix_xcms_chromPeaks_mz_width: %d/%d chromPeaks updated (mz width reset to %.2f ppm)",
      n_fix, n_total, ppm
    ))
  }
  xcms.xcms
}


#' Fill NA `beta_cor` / `beta_snr` via `chromPeakSummary()`.
#' Finite CentWave scores are never overwritten (merged peaks are NA after
#' refine and are the intended fill targets).
#' @noRd
recalc_xcms_chromPeaks_beta <- function(xcms.xcms, BPPARAM, chunkSize = NULL,
                                        verbose = TRUE) {
  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) {
      message("recalc_xcms_chromPeaks_beta: no chromPeaks to score")
    }
    return(xcms.xcms)
  }
  n_na_before <- if ("beta_cor" %in% colnames(pks)) {
    sum(!is.finite(as.numeric(pks[, "beta_cor"])))
  } else {
    nrow(pks)
  }
  n_workers <- tryCatch(
    as.integer(BiocParallel::bpnworkers(BPPARAM)),
    error = function(e) 2L
  )
  if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 2L
  if (is.null(chunkSize)) chunkSize <- n_workers
  beta <- tryCatch(
    as.matrix(xcms::chromPeakSummary(
      xcms.xcms,
      xcms::BetaDistributionParam(),
      BPPARAM = BPPARAM,
      chunkSize = as.integer(chunkSize)
    )),
    error = function(e) {
      if (isTRUE(verbose)) {
        message("recalc_xcms_chromPeaks_beta: chromPeakSummary failed; ",
                conditionMessage(e))
      }
      NULL
    }
  )
  if (is.null(beta) || !nrow(beta)) {
    return(xcms.xcms)
  }
  if (is.null(rownames(beta)) || is.null(rownames(pks)) ||
      (nrow(beta) == nrow(pks) && identical(rownames(beta), rownames(pks)))) {
    if (nrow(beta) != nrow(pks)) {
      if (isTRUE(verbose)) {
        message("recalc_xcms_chromPeaks_beta: row count mismatch; skip write")
      }
      return(xcms.xcms)
    }
    idx <- seq_len(nrow(pks))
  } else {
    idx <- match(rownames(pks), rownames(beta))
  }
  n_filled <- 0L
  for (col in c("beta_cor", "beta_snr")) {
    if (!col %in% colnames(beta)) next
    vals <- if (col %in% colnames(pks)) {
      as.numeric(pks[, col])
    } else {
      rep(NA_real_, nrow(pks))
    }
    fill <- !is.na(idx) & !is.finite(vals)
    if (any(fill)) {
      newv <- as.numeric(beta[idx[fill], col])
      vals[fill] <- newv
      if (identical(col, "beta_cor")) {
        n_filled <- sum(is.finite(newv))
      }
    }
    if (col %in% colnames(pks)) {
      pks[, col] <- vals
    } else {
      pks <- cbind(pks, matrix(vals, ncol = 1L, dimnames = list(rownames(pks), col)))
    }
  }
  xcms::chromPeaks(xcms.xcms) <- pks
  if (isTRUE(verbose)) {
    n_na_after <- sum(!is.finite(as.numeric(pks[, "beta_cor"])))
    message(sprintf(
      "recalc_xcms_chromPeaks_beta: beta_cor NA %d -> %d / %d chromPeaks (%d filled, finite CentWave scores kept)",
      n_na_before, n_na_after, nrow(pks), n_filled
    ))
  }
  xcms.xcms
}

#' Filter chromPeaks by `beta_cor` (drop NA after re-calc).
#' @noRd
filter_xcms_chromPeaks_beta_cor <- function(xcms.xcms, thresh, verbose = TRUE) {
  if (is.null(thresh)) {
    return(xcms.xcms)
  }
  thresh <- as.numeric(thresh)
  if (length(thresh) != 1L || !is.finite(thresh)) {
    stop("beta_cor_thresh must be a single finite number")
  }
  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) {
      message("filter_xcms_chromPeaks_beta_cor: no chromPeaks to filter")
    }
    return(xcms.xcms)
  }
  if (!"beta_cor" %in% colnames(pks)) {
    if (isTRUE(verbose)) {
      message("filter_xcms_chromPeaks_beta_cor: beta_cor column missing; skip")
    }
    return(xcms.xcms)
  }
  beta_cor <- as.numeric(pks[, "beta_cor"])
  n_na <- sum(!is.finite(beta_cor))
  n_pass <- sum(is.finite(beta_cor) & beta_cor >= thresh)
  n_below <- sum(is.finite(beta_cor) & beta_cor < thresh)
  keep <- is.finite(beta_cor) & beta_cor >= thresh
  n_total <- nrow(pks)
  if (isTRUE(verbose)) {
    message(sprintf(
      "filter_xcms_chromPeaks_beta_cor: %d/%d peaks with beta_cor >= %.3f were kept (%d NA removed, %d below thresh removed)",
      n_pass, n_total, thresh, n_na, n_below
    ))
  }
  if (sum(!keep) > 0) {
    xcms.xcms <- xcms::filterChromPeaks(xcms.xcms, keep = keep)
  }
  xcms.xcms
}

#' Enable CentWave `verboseBetaColumns` when the param slot exists.
#' @noRd
.enable_verboseBetaColumns <- function(param) {
  if (methods::.hasSlot(param, "verboseBetaColumns")) {
    param@verboseBetaColumns <- TRUE
  }
  param
}

#' @title xcmsProcessingMS1
#' @description Import `msDataFiles`, filter `ion_mode`, find peaks using `centWaveParam`, correct RT, group peaks using `peaksGroup`, fill peaks by xcms at MS1 Level
#' @param msDataFiles `char` ms file (full) paths
#' @param ion_mode to filter ion_mode, 1: positive, 0: negative, import when scans with both pos and neg
#' @param peaksGroup `vector` to xcms::PeakGroupsParam(sampleGroups), should contain "QC"
#' @param centWaveParam xcms::CentWaveParam()
#' @param beta_cor_thresh optional numeric; if set, fill NA `beta_cor`
#'   after merge with \code{chromPeakSummary()} (finite CentWave scores are
#'   kept), then drop chromPeaks below this value (remaining NA scores are
#'   dropped). Default \code{NULL} skips fill and filtering.
#'
#' @return xcms
#' @export

xcmsProcessingMS1 <- function(xcms.xcms,
                              ion_mode = NA,
                              xcms_param = list(
                                findChromPeaks = xcms::CentWaveParam(),
                                groupChromPeaks = xcms::PeakDensityParam(sampleGroups = "A")
                              ),
                              adjustRT = T,
                              chromPeaks_fix_mz_ppm = NULL,
                              chromPeaks_max_mz_ppm = NULL,
                              beta_cor_thresh = NULL,
                              BPPARAM  = BiocParallel::SnowParam(
                                workers = max(1L, floor(parallel::detectCores() / 3)),
                                progressbar = T
                              ),
                              ...){



  if (is.na(ion_mode)) {
    ion_mode <- unique(.xcms_polarity(xcms.xcms))
    if (length(ion_mode)!=1) {
      stop("MS1 scans contain both positive and negative, please check")
    }
  }

  xcms.xcms <- ProtGenerics::filterPolarity(xcms.xcms, ion_mode)

  xcms_param$findChromPeaks <- .enable_verboseBetaColumns(
    xcms_param$findChromPeaks
  )

  ### Find peaks
  message_with_time(" Find peaks...")
  xcms.xcms<-xcms::findChromPeaks(xcms.xcms,
                            param = xcms_param$findChromPeaks,
                            chunkSize = -1,
                            BPPARAM  = BPPARAM,...)
  message_with_time(" ", nrow(xcms::chromPeaks(xcms.xcms)), " chromPeaks found")

  message_with_time(" Merge neighboring peaks...")
  mpp <- xcms::MergeNeighboringPeaksParam(expandRt = 3,minProp = 0.5,ppm =  xcms_param$findChromPeaks@ppm)
  xcms.xcms <- xcms::refineChromPeaks(xcms.xcms, mpp,
                                      BPPARAM  = BPPARAM)

  message_with_time(" Filter chromPeaks...")
  xcms.xcms <- xcms_filter_peaks_NA(xcms.xcms)
  if (!is.null(chromPeaks_fix_mz_ppm)) {
    xcms.xcms <- fix_xcms_chromPeaks_mz_width(
      xcms.xcms,
      ppm = as.numeric(chromPeaks_fix_mz_ppm)
    )
  }
  if (!is.null(chromPeaks_max_mz_ppm)) {
    xcms.xcms <- filter_xcms_chromPeaks_mz_width(
      xcms.xcms,
      ppm = as.numeric(chromPeaks_max_mz_ppm)
    )
  }
  if (!is.null(beta_cor_thresh)) {
    message_with_time(" Fill NA chromPeak beta_cor...")
    xcms.xcms <- recalc_xcms_chromPeaks_beta(
      xcms.xcms,
      BPPARAM = BPPARAM
    )
    xcms.xcms <- filter_xcms_chromPeaks_beta_cor(
      xcms.xcms,
      thresh = beta_cor_thresh
    )
  }

  ### adujust RT
  if(adjustRT){

    message_with_time(" Group peaks (for RT adjustment)...")
    peaksGroup <- Biobase::pData(xcms.xcms)$sample.type
    peak.density.param <- xcms::PeakDensityParam(sampleGroups = peaksGroup,
                                                 minFraction = 0.4,bw = 30,
                                                 binSize = 0.015)
    xcms.xcms <- xcms::groupChromPeaks(xcms.xcms,param = peak.density.param)

    if (length(Biobase::sampleNames(xcms.xcms))>1) {
      message_with_time(" Adjust RT...")
      if (sum(peaksGroup=="QC") <2 ) {
        rt.adjust.param <- xcms::PeakGroupsParam(minFraction = 0.4,
                                                 #subset = which(peaksGroup == "QC"),
                                                 subsetAdjust = "previous",span = 0.4)
        xcms.xcms <- xcms::adjustRtime(xcms.xcms,param = rt.adjust.param)
      }else{
        ### adjust based on QC
        rt.adjust.param <- xcms::PeakGroupsParam(minFraction = 0.4,
                                                 subset = which(peaksGroup == "QC"),
                                                 subsetAdjust = "average",span = 0.4)
        xcms.xcms <- xcms::adjustRtime(xcms.xcms,param = rt.adjust.param)
      }
    }

  }



  ### group peaks
  message_with_time(" Group peaks...")
  peak.density.param <- xcms_param$groupChromPeaks
  peak.density.param@sampleGroups <- Biobase::pData(xcms.xcms)$sample.type
  xcms.xcms <- xcms::groupChromPeaks(xcms.xcms,param = peak.density.param)
  #xcms.xcms <- xcms_filter_feature_mz_rsd(xcms.xcms,rsd.ppm = 2)
  xcms.xcms <- xcms_get_feature_wmean(xcms.xcms)
  message_with_time(" ",nrow(xcms::featureDefinitions(xcms.xcms))," feature found")


  ### fill peaks
  message_with_time(" Fill peaks...")
  fill_param <- if (inherits(xcms.xcms, "XcmsExperiment") ||
                    inherits(xcms.xcms, "MsExperiment")) {
    xcms::ChromPeakAreaParam()
  } else {
    xcms::FillChromPeaksParam()
  }
  xcms.xcms <- xcms::fillChromPeaks(xcms.xcms,
                                    chunkSize = length(sampleNames(xcms.xcms)),
                                    param = fill_param,
                                    BPPARAM = BPPARAM
                                    )



  return(xcms.xcms)




}

xcmsProcessingMRM <- function(msDataFiles, peaksGroup =NA,
                              centWaveParam ){

  xcms.mrm <- readSRMData(msDataFiles)
  xcms.peaks <- findChromPeaks()


}







matchSpectra_Features <- function(xcmsFeatureDef, spec){

  .matchSP <- function(x,xcmsFeatureDef,
                       mz_ppm = 10,
                       rt_tol = 10){
    mz <- x[["precursorMz"]]%>%as.numeric()
    rt <- x[["rtime"]]%>%as.numeric()
    mzError <- abs((mz - xcmsFeatureDef$mzmed)/mz*1e6)
    rtError <- abs((rt- xcmsFeatureDef$rtmed)/rt)
    feature_id <- rownames(xcmsFeatureDef)[mzError < mz_ppm &rtError < rt_tol]
    #ifelse(length(feature_id)==0, NA,feature_id)


  }
  spec.data <- as.data.frame(Spectra::spectraData(spec))
  feature_id <- apply(spec.data, 1,.matchSP,xcmsFeatureDef)
  spec$feature_id <- feature_id
  spec

}






#' @describeIn xcms_extension_plot plot feature intensity
#' @description plot feature's intensity, ordered by `Biobase::pData(xcms.xcms)$analysis.time.positive` or
#'  `Biobase::pData(xcms.xcms)$analysis.time.negative`
#'
#' @param xcms.xcms XCMSnExp object
#' @param feature_id_to_show feature id to plot
#'
#' @return ggplot object
#' @export
#'

plot_xcms_feature_intensity <- function(xcms.xcms , feature_id_to_show ){

  ion_mode <- unique(as.integer(.xcms_polarity(xcms.xcms)))
  if (length(ion_mode) != 1L || is.na(ion_mode)) {
    stop("Cannot determine unique polarity for xcms.xcms")
  }
  if (ion_mode==1) {
    sample.info <- Biobase::pData(xcms.xcms)%>%
      dplyr::arrange(analysis.time.positive)%>%
      dplyr::mutate(sample.type = factor(sample.type,levels = c("Blank","QC","Sample")),
                    injecton.order = 1:nrow(.))
  }else{

    sample.info <- Biobase::pData(xcms.xcms)%>%
      dplyr::arrange(analysis.time.negative)%>%
      dplyr::mutate(sample.type = factor(sample.type,levels = c("Blank","QC","Sample")),
                    injecton.order = 1:nrow(.))
  }
  features <- xcms::featureValues(xcms.xcms)%>%
    as.data.frame()%>%
    rownames_to_column("feature_id")%>%
    dplyr::filter(feature_id %in% feature_id_to_show )%>%
    dplyr::select(sample.info$sampleNames)%>%as.numeric()

  sample.info$intensity <- features
  sample.info$intensity[is.na(sample.info$intensity )] <- 0
  ggplot(sample.info,aes(x = injecton.order , y = intensity , col = sample.type,na.rm =T))+
    geom_point(size = 0.5)+
    scale_color_manual(values = c("grey","#66CAB7","#EE8E5B"))+
    theme_bw()+
    theme(text = element_text(size = 8))


}




get_xcms_scan_Stat <- function(xcms.xcms){

  if (inherits(xcms.xcms, "MsExperiment") || inherits(xcms.xcms, "XcmsExperiment")) {
    sp <- ProtGenerics::spectra(xcms.xcms)
    sd <- as.data.frame(Spectra::spectraData(sp), stringsAsFactors = FALSE)
    origins <- as.character(sd$dataOrigin)
    if (!length(origins) || all(is.na(origins))) {
      origins <- as.character(sd$dataStorage)
    }
    file_levels <- unique(origins)
    sd$fileIdx <- match(origins, file_levels)
    sd$spIdx <- as.integer(seq_len(nrow(sd)))
    if (!"rtime" %in% names(sd)) {
      sd$rtime <- as.numeric(Spectra::rtime(sp))
    }
    if (!"msLevel" %in% names(sd)) {
      sd$msLevel <- as.integer(Spectra::msLevel(sp))
    }
    if (!"tic" %in% names(sd)) {
      sd$tic <- tryCatch(as.numeric(Spectra::tic(sp)), error = function(e) NA_real_)
    }
    xcms.fdata <- sd %>%
      dplyr::mutate(fileStr = num2str(fileIdx),
                    spStr = num2str(spIdx)) %>%
      dplyr::group_by(fileStr) %>%
      dplyr::mutate(x = 2 - msLevel,
                    ms1_no = cumsum(x)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(ms1_no_str = num2str(ms1_no)) %>%
      dplyr::group_by(fileIdx) %>%
      dplyr::arrange(fileIdx, rtime) %>%
      dplyr::mutate(scan_time = c(diff(rtime), 0),
                    ms1_group = paste0(fileStr, "_", ms1_no_str)) %>%
      dplyr::group_by(ms1_group) %>%
      dplyr::mutate(ms2_count = sum(msLevel == 2),
                    ms1_group_rt = min(rtime),
                    cycle_time = max(rtime) - min(rtime)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fileStr) %>%
      dplyr::mutate(cycle_time = c(diff(ms1_group_rt), 0)) %>%
      dplyr::group_by(ms1_group) %>%
      dplyr::mutate(cycle_time = max(cycle_time)) %>%
      dplyr::mutate(scan_id = paste0("scan_", fileStr, "_", spStr)) %>%
      dplyr::ungroup() %>%
      dplyr::select(scan_id, ms1_no, ms1_group, ms1_group_rt,
                    ms2_count, cycle_time, scan_time, dplyr::everything(),
                    -c(x, fileStr, spStr, ms1_no_str)) %>%
      as.data.frame()
    rownames(xcms.fdata) <- xcms.fdata$scan_id
    return(xcms.fdata)
  }

  stop(
    "OnDiskMSnExp / XCMSnExp scan tables are no longer supported; ",
    "use MsExperiment / XcmsExperiment (Spectra columns: precursorMz, rtime).",
    call. = FALSE
  )
}


xcms_get_scan_Stat <- function(xcms.xcms){
  xcms.fdata <- get_xcms_scan_Stat(xcms.xcms )
  if (inherits(xcms.xcms, "MsExperiment") || inherits(xcms.xcms, "XcmsExperiment")) {
    ## No writable Biobase::fData on MsExperiment/XcmsExperiment; return scan table only.
    attr(xcms.xcms, "scan_stat") <- xcms.fdata
    return(xcms.xcms)
  }
  fData(xcms.xcms) <- xcms.fdata
  return(xcms.xcms)
}




#' @description Plot MS1 total ion chromatograms (TIC) for an XCMSnExp object,
#' colored by sample group from `Biobase::pData(xcms.xcms)$group`.
#' @describeIn xcms_extension_plot plot TIC
#'
#' @param xcms.xcms XCMSnExp object
#' @param col.group named character vector of colors for groups. If `NULL`,
#'   Blank/QC use fixed colors and remaining groups use `ggsci::pal_aaas()`
#'   (or an interpolated palette when there are more than 10 groups)
#' @param title plot title
#'
#' @return ggplot object
#' @export
#'
plot_xcms_TIC <- function(xcms.xcms,col.group = NULL,title = "TIC"){


  xcms.pdata <- Biobase::pData(xcms.xcms)
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  if (!"tic" %in% names(xcms.scan) || all(is.na(xcms.scan$tic))) {
    if (inherits(xcms.xcms, "MsExperiment") || inherits(xcms.xcms, "XcmsExperiment")) {
      xcms.scan$tic <- as.numeric(Spectra::tic(ProtGenerics::spectra(xcms.xcms)))
    } else {
      xcms.scan$tic <- as.numeric(MSnbase::tic(xcms.xcms))
    }
  }
  xcms.scan <- xcms.scan%>%
    dplyr::mutate(group = xcms.pdata$group[fileIdx])%>%
    dplyr::filter(msLevel==1)

  if (is.null(col.group)) {
    groups <- setdiff(unique(xcms.scan$group), c("Blank", "QC"))
    n_groups <- length(groups)
    sample_cols <- if (n_groups <= 10) {
      ggsci::pal_aaas()(max(n_groups, 1))
    } else {
      grDevices::colorRampPalette(ggsci::pal_aaas()(10))(n_groups)
    }
    col.group <- c("grey", "#38C291", sample_cols)
    names(col.group) <- c("Blank", "QC", groups)
  }

  sample_count <- length(unique(xcms.scan$fileIdx))
  line_alpha <- max(1 / max(sample_count, 1), 0.1)

  ggplot(xcms.scan)+
    geom_line(aes(x = rtime , y = tic,
                  col = group,
                  group=fileIdx),alpha = line_alpha)+
    scale_color_manual(values = col.group)+
    #scale_y_log10()+
    labs(title = title,x = "Retention Time", y = "Intensity", col = "")+
    theme_classic()->p
  p



}

plot_xcms_adjustedRT <- function(xcms.xcms){


  xcms.pdata <- Biobase::pData(xcms.xcms)%>%
    dplyr::arrange(ExpTime)%>%
    dplyr::mutate(injection_order = 1:n())
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.scan <- xcms.scan%>%
    dplyr::mutate(adrt = xcms::adjustedRtime(xcms.xcms),
                  group = xcms.pdata$group[fileIdx],
                  injection_order = xcms.pdata$injection_order[fileIdx])%>%
    dplyr::filter(msLevel==1)

 #col.scale <- c("grey","#38C291",ggsci::pal_aaas()(10))
 #names(col.scale) <- c("Blank","QC",
 #                      setdiff(unique(xcms.scan$group),c("Blank","QC")))

  ggplot(xcms.scan)+
    geom_line(aes(x = rtime ,
                  y = adrt-rtime,
                  col = injection_order,group=fileIdx))+
    scale_color_gradient(low = "#FFD700",high = "#EE0000")+
    labs(title = "Retention Time adjust",
         x = "Retention Time", y = "Adjusted Error",
         col = "Injection")+
    theme_classic()->p
  p



}


plot_xcms_scan <- function(xcms.xcms){

  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
 #xcms.scan$precursorMz <- estimatePrecursorIntensity(xcms.xcms,
 #                                        BPPARAM = BatchtoolsParam(progressbar = T,
 #                                                                  log = F))


  ggplot(xcms.scan)+
    geom_point(aes(x = rtime ,
                   y = precursorMz ,
                   col = log10(precursorIntensity)))


}


get_xcms_MS_report <- function(xcms.xcms ,
                               file.path){


  file.path <- "d:/temp/xcms.report.pdf"
  p.tic <- plot_xcms_TIC(xcms.xcms)
  p.rtadj <- plot_xcms_adjustedRT(xcms.xcms )
  p.feature.dis <- plot_xcms_features_distribution(xcms.xcms)


  ### scan
  #plot_xcms_peaks_ms1_scans(xcms.xcms)
  plot_xcms_peaks_ms2_scans(xcms.xcms)



  p1 <- p.tic/p.rtad+(p.feature.dis)

  pdf("d:/temp/aaa.pdf")
  plot(p1)
  dev.off()

  }



get_xcms_peaks_stat <- function(xcms.xcms){

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    rownames_to_column("peak_id")%>%
    dplyr::mutate(peakWidth = rtmax-rtmin,
                  mzWidth = mzmax-mzmin,
                  mzError = mzWidth/mz*1e6)
  return(xcms.peaks)

}



get_xcms_centwave_tune <- function(xcms.xcms,
                                   iteration = 10){

  cwp <- xcms::CentWaveParam(peakwidth = c(5,20),
                       verboseColumns=T,fitgauss = T)
  xcms.tune.df <- data.frame(
    No = 1:iteration,
    ppm = cwp@ppm,
    pwmin = cwp@peakwidth[1],
    pwmax = cwp@peakwidth[2],
    snthresh = cwp@snthresh,
    prefilter = cwp@prefilter[1],
    prefilter.int = cwp@prefilter[2],
    peaks.no = NA,
    mze.range = NA,
    pw.range = NA
  )
  xcms.tune.list <- list()
  for (i in 1:iteration) {

    message("run ",i," ", Sys.time())
    show(cwp)

    ### record param
    xcms.tune.df$ppm[i] <- cwp@ppm
    xcms.tune.df$pwmin[i] <- cwp@peakwidth[1]
    xcms.tune.df$pwmax[i] <- cwp@peakwidth[2]
    xcms.tune.df$snthresh[i] <- cwp@snthresh
    xcms.tune.df$prefilter[i] <- cwp@prefilter[1]
    xcms.tune.df$prefilter.int[i] <- cwp@prefilter[2]


    ### run
    xcms.interation <- findChromPeaks(
      xcms.xcms,
      param = cwp,
      BPPARAM  =SerialParam()
    )
    xcms.peaks <- get_xcms_peaks_stat(xcms.interation)%>%
      dplyr::mutate(scan.no = scmax-scmin)
    xcms.peaks.high.sn <- xcms.peaks%>%
      dplyr::slice_max(sn , n = round(nrow(.)*0.1))

    ### result
    pw <- xcms.peaks.high.sn$peakWidth
    mze <- xcms.peaks.high.sn$mzError
    xcms.tune.df$peaks.no[i] <- nrow(xcms.peaks)
    xcms.tune.df$mze.range[i] <- paste0(quantile(mze,c(0.05,0.95))%>%round,collapse = "-")
    xcms.tune.df$pw.range[i] <- paste0(quantile(pw,c(0.25,0.75))%>%round,collapse = "-")

    ### update param
    cwp@ppm <- quantile(mze,0.5)
    cwp@peakwidth <- quantile(pw,c(0.25,0.75))
    #cwp@snthresh <- quantile(xcms.peaks$sn,0.05)
    #cwp@prefilter[1] <- quantile(xcms.peaks$scan.no,0.05)

    ### record
    xcms.tune.list[[i]] <- xcms.interation


  }





}



#' @title get_xcms_Spectra
#' @description
#' **Deprecated.** Use \code{ProtGenerics::spectra(xcms.xcms)} (and optionally
#' \code{Spectra::filterMsLevel()} / \code{Spectra::filterPolarity()}) instead.
#'
#' Previously: build a [Spectra::Spectra] object from the raw files behind an
#' `XCMSnExp`, filtered to the polarity of `xcms.xcms`. Scan IDs from
#' `get_xcms_scan_Stat()` are assigned as `spectraNames` and the
#' `scan_id` spectra variable.
#'
#' @param xcms.xcms An `XCMSnExp` object.
#'
#' @return A `Spectra` object with `scan_id` matching xcms scan metadata.
#' @export
get_xcms_Spectra <- function(xcms.xcms){
  .Deprecated("ProtGenerics::spectra")

  xcms.files <- filepaths(xcms.xcms)
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  message_with_time("loading Spectra from files...")
  xcms.sp <- Spectra::Spectra(xcms.files,
                         backend = Spectra::MsBackendMemory(),
                         BPPARAM = SerialParam(progressbar = T))%>%
    filterPolarity(unique(polarity(xcms.xcms)))
  Spectra::spectraNames(xcms.sp) <- xcms.sp$scan_id <- xcms.scan$scan_id
  return(xcms.sp)

}

## Biobase-compatible accessors for MsExperiment / XcmsExperiment
## (XCMSnExp already provides these via MSnbase inheritance)

#' Polarity vector for XCMSnExp / MsExperiment / XcmsExperiment
#' @noRd
.xcms_polarity <- function(xcms.xcms) {
  if (inherits(xcms.xcms, "MsExperiment") || inherits(xcms.xcms, "XcmsExperiment")) {
    sp <- ProtGenerics::spectra(xcms.xcms)
    if (!length(sp)) return(integer())
    return(as.integer(Spectra::polarity(sp)))
  }
  as.integer(ProtGenerics::polarity(xcms.xcms))
}

#' @importFrom Biobase fData
#' @importFrom methods setMethod
setMethod("fData", "MsExperiment", function(object) {
  .get_xcms_scan_table(object)
})

#' @importFrom ProtGenerics polarity
setMethod("polarity", "MsExperiment", function(object) {
  .xcms_polarity(object)
})

#' Filter MsExperiment / XcmsExperiment spectra by polarity
#'
#' Uses \code{MsExperiment::filterSpectra} so sample–spectra links stay
#' consistent. Covers \code{XcmsExperiment} via inheritance. Does not sync
#' \code{chromPeaks} (use before peak detection).
#'
#' @noRd
#' @importFrom ProtGenerics filterPolarity
#' @importFrom methods setMethod
setMethod(
  "filterPolarity", "MsExperiment",
  function(object, polarity = integer()) {
    if (!length(ProtGenerics::spectra(object))) {
      return(object)
    }
    MsExperiment::filterSpectra(object, Spectra::filterPolarity, polarity)
  }
)

#' @importFrom Biobase pData
#' @importFrom methods setMethod setReplaceMethod
setMethod("pData", "MsExperiment", function(object) {
  as.data.frame(MsExperiment::sampleData(object), stringsAsFactors = FALSE)
})

#' @importFrom Biobase pData<-
setReplaceMethod("pData", "MsExperiment", function(object, value) {
  if (is.null(value)) {
    stop("'value' must be a data.frame or DataFrame")
  }
  value <- as.data.frame(value, stringsAsFactors = FALSE)
  MsExperiment::sampleData(object) <- S4Vectors::DataFrame(value, check.names = FALSE)
  object
})

#' @importFrom Biobase sampleNames
setMethod("sampleNames", "MsExperiment", function(object) {
  sd <- MsExperiment::sampleData(object)
  ## Prefer explicit sample id columns (rownames often = basename(files))
  for (col in c("sampleNames", "sample.name", "sample_name")) {
    if (col %in% colnames(sd) && !all(is.na(sd[[col]]))) {
      return(as.character(sd[[col]]))
    }
  }
  rn <- rownames(sd)
  if (!is.null(rn) && length(rn) && !all(rn == as.character(seq_len(nrow(sd))))) {
    return(as.character(rn))
  }
  fn <- tryCatch(xcms::fileNames(object), error = function(e) NULL)
  if (!is.null(fn) && length(fn)) {
    return(basename(as.character(fn)))
  }
  as.character(seq_len(length(object)))
})

#' @importFrom Biobase sampleNames<-
setReplaceMethod("sampleNames", "MsExperiment", function(object, value) {
  value <- as.character(value)
  sd <- MsExperiment::sampleData(object)
  if (length(value) != nrow(sd)) {
    stop("'value' length must equal number of samples (", nrow(sd), ")")
  }
  ## Keep rownames as file basenames when possible; store ids in columns
  sd$sampleNames <- value
  if ("sample.name" %in% colnames(sd)) {
    sd$sample.name <- value
  }
  MsExperiment::sampleData(object) <- sd
  object
})

#' @importFrom xcms filepaths
#' @export
setMethod(f = "filepaths",
                       signature = "XCMSnExp",
                       definition = function(object) {
                         paste0(BiocGenerics::dirname(object), "/", Biobase::sampleNames(object))
                       })

#' @importFrom xcms filepaths
#' @export
setMethod(f = "filepaths",
                       signature = "MsExperiment",
                       definition = function(object) {
                         fn <- tryCatch(xcms::fileNames(object), error = function(e) NULL)
                         if (is.null(fn) || !length(fn)) {
                           fn <- tryCatch(
                             as.character(MsExperiment::sampleData(object)$dataOrigin),
                             error = function(e) NULL
                           )
                         }
                         if (is.null(fn) || !length(fn)) {
                           stop("Cannot determine file paths for MsExperiment object")
                         }
                         as.character(fn)
                       })

#' @importFrom xcms mzrange
#' @export
setMethod(f = mzrange,
                       signature = "XCMSnExp",
                       definition = function(object) {
                         xcms.fdata <- fData(object)
                         return(c(min(xcms.fdata$scanWindowLowerLimit, na.rm = TRUE),
                                  max(xcms.fdata$scanWindowUpperLimit, na.rm = TRUE)))
                       })

#' @importFrom xcms mzrange
#' @export
setMethod(f = mzrange,
                       signature = "MsExperiment",
                       definition = function(object) {
                         sp <- ProtGenerics::spectra(object)
                         sd <- as.data.frame(Spectra::spectraData(sp))
                         lo <- if ("scanWindowLowerLimit" %in% names(sd)) {
                           sd$scanWindowLowerLimit
                         } else {
                           NA_real_
                         }
                         hi <- if ("scanWindowUpperLimit" %in% names(sd)) {
                           sd$scanWindowUpperLimit
                         } else {
                           NA_real_
                         }
                         if (all(is.na(lo)) || all(is.na(hi))) {
                           mz <- tryCatch(Spectra::mz(sp), error = function(e) NULL)
                           if (!is.null(mz) && length(mz)) {
                             mz_num <- unlist(mz, use.names = FALSE)
                             return(c(min(mz_num, na.rm = TRUE), max(mz_num, na.rm = TRUE)))
                           }
                           return(c(NA_real_, NA_real_))
                         }
                         c(min(lo, na.rm = TRUE), max(hi, na.rm = TRUE))
                       })



get_xcms_precursor_intensity <- function(xcms.xcms,...){

  estimatePrecursorIntensity(xcms.xcms,
                             method = "previous",...)

}




plotly_feature_span <- function(xcms.fdf){

  plot_ly(xcms.fdf)%>%
    add_segments(x = ~rtmin,xend = ~rtmax,
                 y = ~mzmed , yend = ~mzmed)

}



simulate_dda <- function(xcms.fdf,
                         ms1.time = 0.6,
                         ms2.time = 0.6,
                         topn = 10){

  t <- 0
  scan.df <- data.frame(spIdx = 1,
                        msLevel= 1,
                        rtime= 0,
                        precursorMz= NA,
                        ion_id = NA
  )
  while(t < max(xcms.fdf$rtmax)){

    ion.to.ms2 <- which(t < xcms.fdf$rtmax&t>xcms.fdf$rtmin)
    ion.to.ms2 <- na.omit(ion.to.ms2[order(xcms.fdf$peakMaxo[ion.to.ms2],decreasing = T)[1:topn]])
    ms2.scan <- data.frame(spIdx = rep(NA,length(ion.to.ms2)),
                           msLevel= rep(2,length(ion.to.ms2)),
                           rtime= t + seq_along(ion.to.ms2)*ms2.time,
                           precursorMz = xcms.fdf$mzmed[ion.to.ms2],
                           ion_id = ion.to.ms2
    )
    scan.df <- rbind(scan.df,ms2.scan)
    t <- t + ms2.time*length(ion.to.ms2)

    ms1.scan <- data.frame(spIdx = NA,
                           msLevel= 1,
                           rtime = t + ms1.time,
                           precursorMz= NA,
                           ion_id = NA
    )
    scan.df <- rbind(scan.df,ms1.scan)
    t <-  t + ms1.time
  }


  return(scan.df)

}


simulate_prm <- function(xcms.fdf,
                         total.time = 28*60,
                         ms2.time = 0.265){


  t <- 0
  scan.df <- data.frame(spIdx = NULL,
                        msLevel= NULL,
                        rtime= NULL,
                        precursorMz= NULL,
                        ion_id = NULL
  )
  while(t < max(xcms.fdf$rtmax)){

    ion.to.ms2 <- which(t < xcms.fdf$rtmax&t>xcms.fdf$rtmin)
    ion.to.ms2 <- na.omit(ion.to.ms2)
    ms2.scan <- data.frame(spIdx = rep(NA,length(ion.to.ms2)),
                           msLevel= rep(2,length(ion.to.ms2)),
                           rtime= t + seq_along(ion.to.ms2)*ms2.time,
                           precursorMz = xcms.fdf$mzmed[ion.to.ms2],
                           ion_id = ion.to.ms2
    )
    scan.df <- rbind(scan.df,ms2.scan)
    t <- t + ifelse(ms2.time*length(ion.to.ms2)>0,ms2.time*length(ion.to.ms2),0.1)


  }



  return(scan.df)
}


#' @title Compute MS1 purity matrix for xcms features
#' @description
#' Calculate a feature-by-sample MS1 purity matrix by extracting, for each feature
#' in \code{xcms.xcms}, the closest MS1 scan (by retention time) from \code{xcms.ms1.sp}
#' in each sample file (matched by \code{Spectra::dataOrigin}). Purity is calculated
#' within an isolation window around the feature m/z.
#'
#' If \code{xcms.ms1.sp} is not provided, MS1 spectra are taken via
#' \code{ProtGenerics::spectra(xcms.xcms)} and filtered to MS level 1.
#'
#' Optimized path: closest-scan lookup via \code{findInterval}, a single
#' \code{peaksData()} extraction, and purity computed on peak matrices
#' (no per-feature \code{Spectra} subsetting).
#'
#' @param xcms.xcms \code{XCMSnExp} with grouped features (must have \code{featureDefinitions}).
#' @param xcms.ms1.sp Optional MS1 \code{Spectra} covering the same files as \code{xcms.xcms}.
#'   If missing or \code{NULL}, taken from \code{ProtGenerics::spectra()}.
#' @param ppm numeric, ppm tolerance for m/z window.
#' @param isolation_half_window numeric, half isolation window (m/z).
#'
#' @return numeric matrix with rows = \code{feature_id}, columns = sample files (by \code{dataOrigin}).
#' @export
get_xcms_feature_purity_matrix <- function(xcms.xcms,
                                           xcms.ms1.sp = NULL,
                                           ppm = 5,
                                           isolation_half_window = 0.2){

  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  n_features <- nrow(xcms.fdf)
  rtmed <- xcms.fdf$rtmed
  mzmed <- xcms.fdf$mzmed
  ppm_tol <- ppm * 1e-6

  if (missing(xcms.ms1.sp) || is.null(xcms.ms1.sp)) {
    message_with_time("xcms.ms1.sp not provided; extracting MS1 via spectra()...")
    xcms.ms1.sp <- ProtGenerics::spectra(xcms.xcms) %>%
      Spectra::filterMsLevel(1L)
  }

  ### closest MS1 scan index per feature x file
  {
    sp.rt <- Spectra::rtime(xcms.ms1.sp)
    sp.origin <- xcms.ms1.sp$dataOrigin
    sp.idx.split <- split(seq_along(xcms.ms1.sp), sp.origin)
    sp.rt.split <- split(sp.rt, sp.origin)
    origins <- names(sp.idx.split)
    n_files <- length(origins)

    f.sp.idx <- matrix(NA_integer_, nrow = n_features, ncol = n_files)
    for (j in seq_len(n_files)) {
      idx_j <- sp.idx.split[[j]]
      rt_j <- sp.rt.split[[j]]
      ord <- order(rt_j)
      rt_sorted <- rt_j[ord]
      n_j <- length(rt_sorted)
      if (n_j == 1L) {
        f.sp.idx[, j] <- idx_j[1L]
      } else {
        pos <- findInterval(rtmed, rt_sorted)
        pos <- pmax(1L, pmin(pos, n_j))
        pos2 <- pmin(pos + 1L, n_j)
        best <- ifelse(abs(rt_sorted[pos] - rtmed) <= abs(rt_sorted[pos2] - rtmed),
                       pos, pos2)
        f.sp.idx[, j] <- idx_j[ord[best]]
      }
    }
  }

  ### purity from peak matrices (single peaksData call)
  {
    message_with_time("extracting peaks data...")
    all_pd <- as.list(Spectra::peaksData(xcms.ms1.sp))

    message_with_time("calculating MS1 purity per feature x sample...")
    ms1_purity_matrix <- matrix(0, nrow = n_features, ncol = n_files)
    colnames(ms1_purity_matrix) <- basename(origins)
    rownames(ms1_purity_matrix) <- xcms.fdf$feature_id

    pb <- get_progress_bar(n_features)
    for (i in seq_len(n_features)) {
      ion_mz <- mzmed[i]
      for (j in seq_len(n_files)) {
        x <- all_pd[[f.sp.idx[i, j]]]
        mz.diff <- abs(x[, 1L] - ion_mz)
        num <- sum(x[mz.diff / ion_mz < ppm_tol, 2L])
        den <- sum(x[mz.diff < isolation_half_window, 2L])
        p <- num / den
        if (is.finite(p)) ms1_purity_matrix[i, j] <- p
      }
      pb$tick()
    }
  }

  return(ms1_purity_matrix)

}



#' @title Store MS1 feature purity matrix in XcmsExperiment qdata
#' @description
#' Compute a feature-by-sample MS1 purity matrix and store it as assay
#' \code{"purity_matrix"} in \code{MsExperiment::qdata(xcms.xcms)}.
#'
#' Requires an \code{XcmsExperiment}. MS1 spectra are taken from
#' \code{ProtGenerics::spectra(xcms.xcms)} filtered with
#' \code{Spectra::filterMsLevel(1L)}. The matrix is computed by
#' \code{\link{get_xcms_feature_purity_matrix}}. If \code{qdata} is missing,
#' it is created with \code{xcms::quantify()} (assay \code{"raw"}), then
#' assay \code{"purity_matrix"} is added. Does not write aggregated purity
#' into \code{featureDefinitions}.
#'
#' @param xcms.xcms \code{XcmsExperiment} with grouped features.
#' @param ppm numeric, ppm tolerance for m/z window.
#' @param isolation_half_window numeric, half isolation window (m/z).
#'
#' @return \code{XcmsExperiment} with \code{qdata} containing assay
#'   \code{"purity_matrix"} (aligned to existing \code{qdata} rows/columns).
#'
#' @seealso \code{\link{get_xcms_feature_purity_matrix}}
#' @export
xcms_get_feature_purity <- function(xcms.xcms,
                                    ppm = 10,
                                    isolation_half_window = 0.2) {

  if (!inherits(xcms.xcms, "XcmsExperiment")) {
    stop("'xcms.xcms' must be an XcmsExperiment (qdata is not available on XCMSnExp).")
  }

  ### MS1 spectra from xcms object
  message_with_time("extracting MS1 Spectra from xcms...")
  xcms.ms1.sp <- ProtGenerics::spectra(xcms.xcms) %>%
    Spectra::filterMsLevel(1L)

  ### calc ms1_purity_matrix
  message_with_time("computing feature purity matrix...")
  ms1_purity_matrix <- get_xcms_feature_purity_matrix(
    xcms.xcms,
    xcms.ms1.sp = xcms.ms1.sp,
    ppm = ppm,
    isolation_half_window = isolation_half_window
  )

  ### ensure qdata (SummarizedExperiment)
  se <- MsExperiment::qdata(xcms.xcms)
  if (is.null(se)) {
    message_with_time("qdata missing; creating via quantify()...")
    se <- xcms::quantify(xcms.xcms)
  }
  if (!inherits(se, "SummarizedExperiment")) {
    stop("'qdata(xcms.xcms)' must be a SummarizedExperiment (or NULL).")
  }

  ### align purity matrix to qdata dimensions
  {
    se_rows <- rownames(se)
    if (is.null(se_rows) || !length(se_rows) || all(se_rows == as.character(seq_len(nrow(se))))) {
      rd <- SummarizedExperiment::rowData(se)
      if ("feature_id" %in% colnames(rd)) {
        se_rows <- as.character(rd$feature_id)
      } else {
        se_rows <- as.character(seq_len(nrow(se)))
      }
    }
    pm_rows <- rownames(ms1_purity_matrix)
    row_idx <- match(se_rows, pm_rows)
    if (anyNA(row_idx)) {
      if (nrow(ms1_purity_matrix) == nrow(se)) {
        row_idx <- seq_len(nrow(se))
      } else {
        stop("Cannot align purity matrix rows to qdata (feature ids / dimensions mismatch).")
      }
    }

    se_cols <- colnames(se)
    if (is.null(se_cols) || !length(se_cols)) {
      se_cols <- basename(as.character(xcms::fileNames(xcms.xcms)))
    }
    pm_cols <- colnames(ms1_purity_matrix)
    col_idx <- match(se_cols, pm_cols)
    if (anyNA(col_idx)) {
      ## try basename(fileNames) vs sampleNames bridge
      fn_base <- basename(as.character(xcms::fileNames(xcms.xcms)))
      sn <- tryCatch(as.character(Biobase::sampleNames(xcms.xcms)), error = function(e) fn_base)
      alt <- se_cols
      map_fn <- match(se_cols, sn)
      alt[!is.na(map_fn)] <- fn_base[map_fn[!is.na(map_fn)]]
      col_idx <- match(alt, pm_cols)
    }
    if (anyNA(col_idx)) {
      if (ncol(ms1_purity_matrix) == ncol(se)) {
        col_idx <- seq_len(ncol(se))
      } else {
        stop("Cannot align purity matrix columns to qdata (sample names / dimensions mismatch).")
      }
    }

    aligned <- ms1_purity_matrix[row_idx, col_idx, drop = FALSE]
    dimnames(aligned) <- list(rownames(se), colnames(se))
  }

  SummarizedExperiment::assay(se, "purity_matrix") <- aligned
  MsExperiment::qdata(xcms.xcms) <- se
  xcms.xcms
}



cbind_Chromatograms <- function(...){

  chrom.list <- list(...)
  if (length(chrom.list) == 1) {
    return(chrom.list[[1]])
  }
  chrom.featureDefinitions <- sapply(chrom.list,
                                     function(x)x@featureDefinitions)[[1]]
  chrom.phenoData <- sapply(chrom.list,
                                     function(x)x@phenoData)%>%
    do.call(Biobase::combine,.)
  chrom.featureData <- sapply(chrom.list,
                                     function(x)x@featureData)%>%
    do.call(Biobase::combine,.)
  chrom.Data<- sapply(chrom.list,
                              function(x)x@.Data)
  chrom.processHistory<- sapply(chrom.list,
                              function(x)x@.processHistory)[1]
  xchrom <-xcms::XChromatograms(chrom.Data,
                     ncol = ncol(chrom.Data),
                     phenoData = chrom.phenoData)
  xchrom@featureDefinitions <- chrom.featureDefinitions
  xchrom@featureData <- chrom.featureData
  xchrom@.processHistory <- chrom.processHistory
  return(xchrom)
}



get_xcms_quantify_MSIP <- function(xcms.xcms){

  xcms::quantify(xcms.xcms,missing = 1,method="max",value = "into")

}




xcms_from_ms2_spectra <- function(sp.ms2 ,
                                  sample.info,
                                  ppm = 10,
                                  peak_width = 60){


  ### assign
  {
    sample.info <- sample.info%>%
      dplyr::mutate(sample = as.numeric(factor(msData.files)))

    sp.peaks.df <- data.frame(
      mz = precursorMz(sp.ms2),
      rt = rtime(sp.ms2),
      sample.files = dataOrigin(sp.ms2),
      into = sp.ms2$totIonCurrent
    )%>%
      dplyr::mutate(sample = match_path(sample.files,
                                        sample.info$msData.files),
                    sample = sample.info$sample[sample])

    sp.peaks.matrix <- sp.peaks.df%>%
      dplyr::mutate(mzmin = mz,
                    mzmax = mz,
                    rtmin = rt-peak_width/2,
                    rtmax = rt+peak_width/2)%>%
      dplyr::select(any_of(c("mz","mzmin","mzmax",
                             "rt","rtmin","rtmax",
                             "into","intb","maxo",
                             "sn","sample")))%>%
      as.matrix()
    sp.peaks.data <- sp.peaks.df%>%
      dplyr::mutate(ms_level = 1,
                    ms_level = as.integer(ms_level),
                    is_filled = F)%>%
      S4Vectors::DataFrame()

    ion_df <- do_groupChromPeaks_density(sp.peaks.df,
                                         bw = peak_width,
                                         sampleGroups = sample.info$sample.source,
                                         binSize = 0.001,
                                         ppm = ppm)

    ion_table <-ion_df %>%
      dplyr::mutate(feature_id = paste0("FTS",num2str(1:n())),
                    .before = mzmed)

  }


  ### simulate xcms class
  {

    MsFeatureData <- new("MsFeatureData",
                         chromPeaks = sp.peaks.matrix,
                         chromPeakData = sp.peaks.data,
                         featureDefinitions =  S4Vectors::DataFrame(ion_table))

    XCMSnExp <- new("XCMSnExp")
    XCMSnExp@msFeatureData <- MsFeatureData
  }


  return(XCMSnExp)

}


xcms_get_feature_adduct_connection <- function(xcms.xcms,rt.tol = 5,ppm = 10){


  pol <- unique(polarity(xcms.xcms))
  adduct.diff <- MSCC::get_adduct_mass_diff(pol)
  adduct.diff <- adduct.diff[order(adduct.diff$mass_diff),]
  #xcms.xcms <- xcms_get_feature_group(xcms.xcms,diffRt = 10,intCor = NULL,eicCor = NULL)
  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)


  ### Construct connection
  {




    # Generate all connection
    # filter rt diff
    # calc mz diff
    # pre-filter mz.diff
    # calc mz.mean
    {

      xcms.net <- expand.grid(
        from = 1:nrow(xcms.fdf),
        to = 1:nrow(xcms.fdf)
      )
      xcms.net <- data.table::as.data.table(xcms.net)
      xcms.net <- xcms.net[from != to ][
        , rt.diff := abs(xcms.fdf$rtmed[to]-xcms.fdf$rtmed[from]) ][
        rt.diff < rt.tol,][
          , c("mz.diff") := .( xcms.fdf$mzmed[to] - xcms.fdf$mzmed[from])][
            mz.diff > min(adduct.diff$mass_diff)&mz.diff < max(adduct.diff$mass_diff) ][
              ,mz.mean := xcms.fdf$mzmed[to] + xcms.fdf$mzmed[from]   ]
    }


    # match mz.diff to adduct.diff
    match.df <- match_mz_foverlaps(mz1 = xcms.net$mz.diff,mz2 = adduct.diff$mass_diff,
                                   ppm.base = xcms.net$mz.mean,ppm = ppm)


    # add adduct.diff data
    xcms.net.matched <- xcms.net[match.df$ion1,][
      ,c("adduct.diff.idx","mz.ppm"):= .(match.df$ion2,match.df$mz.ppm)
    ][mz.ppm  < ppm,][ ### connect within ppm
      ,c("adduct.mass.diff","from.adduct","to.adduct"):= .(adduct.diff$mass_diff[adduct.diff.idx],
                                              adduct.diff$adduct.from[adduct.diff.idx],
                                              adduct.diff$adduct.to[adduct.diff.idx])
    ]


  }


  ###
  {

    xcms.net.matched$label  <- paste0(xcms.net.matched$from.adduct," to ",xcms.net.matched$to.adduct)
    xcms.fdf.ig <- igraph::graph_from_data_frame(xcms.net.matched)
    visIgraph(igraph_filter_distance(xcms.fdf.ig,from = "11959",1))%>%
      visEdges(smooth = T)


    node.group <- igraph::components(xcms.fdf.ig)$membership

    ig <- igraph_filter_vertex(xcms.fdf.ig , which(node.group==1))

    eda <- edata(ig)
    vda <- vdata(ig)
    plot_density(eda$rt.diff)

  }

}


plotly_xcms_feature_group <- function(xcms.xcms){


  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()
  #ggplot(xcms.fdf)+
  #  geom_point(aes(x = rtmed , y = mzmed, col = feature_group))+
  #  theme(legend.position = "none")

  plotly::plot_ly(xcms.fdf)%>%
    add_markers(x = ~rtmed, y = ~mzmed, color = ~feature_group)

}


#' @title Update feature mz/rt using peak-intensity weighted means
#' @description Recomputes `mzmed` and `rtmed` in `xcms::featureDefinitions(xcms.xcms)`
#' using peak-level `mz`/`rt` weighted by peak `maxo` (maximum intensity) from
#' `xcms::chromPeaks(xcms.xcms)`.
#'
#' @param xcms.xcms An `xcms::XCMSnExp` object with feature definitions and chromPeaks.
#'
#' @return An updated `xcms::XCMSnExp` object where `featureDefinitions(object)$mzmed`
#' and `featureDefinitions(object)$rtmed` are replaced by the intensity-weighted means
#' across each feature's constituent peaks.
#'
#' @export
xcms_get_feature_wmean <- function(xcms.xcms){

  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  xcms.pks <- xcms::chromPeaks(xcms.xcms)

  wrt <- sapply(xcms.fdf$peakidx,function(x){

    x.rt <- xcms.pks[x,"rt"]
    x.int <- xcms.pks[ x, "maxo"]
    weighted.mean(x.rt,x.int)
  })

  wmz <- sapply(xcms.fdf$peakidx,function(x){

    x.mz <- xcms.pks[x,"mz"]
    x.int <- xcms.pks[ x, "maxo"]
    weighted.mean(x.mz,x.int)
  })

  xcms.fdf$mzmed <- wmz
  xcms.fdf$rtmed <- wrt
  xcms.fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)
}


xcms_filter_feature_mz_rsd <- function(xcms.xcms, rsd.ppm = 2){


  fdf <- xcms::featureDefinitions(xcms.xcms)
  ch <- xcms::chromPeaks(xcms.xcms)
  mz.sd <- sapply(fdf$peakidx,function(x){
    sd(ch[x,'mz'])/mean(ch[x,'mz']) * 1e6
  })
  #plot_density(mz.sd)
  fdf <- fdf[mz.sd < rsd.ppm,]
  fdf$feature_id <- paste0("FT",num2str(1:nrow(fdf)))
  rownames(fdf) <- fdf$feature_id
  fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)

}
xcms_filter_feature_rt_rsd <- function(xcms.xcms, rt.shift = 5 ){


  fdf <- xcms::featureDefinitions(xcms.xcms)
  ch <- xcms::chromPeaks(xcms.xcms)
  rt.sd <- sapply(fdf$peakidx,function(x){
    sd(ch[x,'rt'])/mean(ch[x,'rt'])
  })
  #plot_density(mz.sd)
  fdf <- fdf[rt.sd < rt.shift,]
  fdf$feature_id <- paste0("FT",num2str(1:nrow(fdf)))
  rownames(fdf) <- fdf$feature_id
  fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)

}

#' Extract summed chromatogram trace for ggplot XIC
#' @noRd
.extract_xcms_xic_chromatogram <- function(xcms.filt, mzr, rtr) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required for XIC plots.", call. = FALSE)
  }
  chr <- xcms::chromatogram(xcms.filt, mz = mzr, rt = rtr)
  chr.sp <- chr[[1L]]
  rt_vals <- as.numeric(MSnbase::rtime(chr.sp))
  int_vals <- as.numeric(MSnbase::intensity(chr.sp))
  int_vals[!is.finite(int_vals)] <- 0
  data.frame(rt = rt_vals, intensity = int_vals, stringsAsFactors = FALSE)
}

#' Extract m/z–RT points for ggplot XIC
#' @noRd
.extract_xcms_xic_points <- function(xcms.filt) {
  sps <- spectra(xcms.filt)
  if (!length(sps)) {
    return(data.frame(
      rt = numeric(),
      mz = numeric(),
      intensity = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- lapply(sps, function(sp) {
    mz_vals <- mz(sp)
    if (!length(mz_vals)) {
      return(NULL)
    }
    data.frame(
      rt = rep(rtime(sp), length(mz_vals)),
      mz = mz_vals,
      intensity = intensity(sp),
      stringsAsFactors = FALSE
    )
  })
  pts <- do.call(rbind, out)
  if (is.null(pts) || !nrow(pts)) {
    return(data.frame(
      rt = numeric(),
      mz = numeric(),
      intensity = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  pts[is.finite(pts$rt) & is.finite(pts$mz) & is.finite(pts$intensity), , drop = FALSE]
}

#' @describeIn xcms_extension_plot ggplot2 XIC plot matching xcms \code{plot(type = \"XIC\")}
#'
#' Upper panel: extracted ion chromatogram (intensity vs retention time).
#' Lower panel: m/z vs retention time with points coloured by intensity.
#'
#' @param xcms.filt \code{XCMSnExp} after \code{filterRt()} and \code{filterMz()}.
#' @param mzr Optional m/z range used for extraction (for axis limits).
#' @param rtr Optional RT range used for extraction (for axis limits).
#' @param title Optional plot title.
#' @param subtitle Optional subtitle.
#' @param base_size Base font size.
#' @param return.data If \code{TRUE}, return a list with \code{plot}, \code{p_chr},
#'   \code{p_mz}, \code{chrom}, and \code{points}.
#'
#' @returns A patchwork object with two panels, or a list when
#'   \code{return.data = TRUE}.
#' @export
plot_xcms_xic <- function(
    xcms.filt,
    mzr = NULL,
    rtr = NULL,
    title = NULL,
    subtitle = NULL,
    base_size = 6,
    return.data = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_xcms_xic().", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for plot_xcms_xic().", call. = FALSE)
  }
  if (is.null(mzr) || is.null(rtr)) {
    stop("mzr and rtr must be supplied.", call. = FALSE)
  }
  mzr <- as.numeric(mzr)
  rtr <- as.numeric(rtr)
  if (length(mzr) != 2L || length(rtr) != 2L || any(!is.finite(c(mzr, rtr)))) {
    stop("mzr and rtr must be finite length-2 numeric vectors.", call. = FALSE)
  }

  chrom.df <- .extract_xcms_xic_chromatogram(xcms.filt, mzr = mzr, rtr = rtr)
  pts.df <- .extract_xcms_xic_points(xcms.filt)

  p.chr <- ggplot2::ggplot(chrom.df, ggplot2::aes(x = .data$rt, y = .data$intensity)) +
    ggplot2::geom_line(linewidth = 0.35, colour = "black") +
    ggplot2::scale_x_continuous(limits = rtr, expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::labs(x = NULL, y = "Intensity") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(2, 4, 0, 6, "pt"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 0, l = 0))
    )

  if (nrow(pts.df)) {
    p.mz <- ggplot2::ggplot(
      pts.df,
      ggplot2::aes(x = .data$rt, y = .data$mz, colour = .data$intensity)
    ) +
      ggplot2::geom_point(size = 0.45, alpha = 0.85) +
      ggplot2::scale_colour_gradientn(
        colours = grDevices::topo.colors(64),
        na.value = NA
      ) +
      ggplot2::scale_x_continuous(limits = rtr, expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
      ggplot2::scale_y_continuous(limits = mzr, expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
      ggplot2::labs(x = "Retention time", y = expression(italic(m/z)), colour = "Intensity") +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 4, 2, 6, "pt"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 0, l = 0)),
        legend.position = "right",
        legend.key.height = grid::unit(0.35, "cm"),
        legend.key.width = grid::unit(0.15, "cm"),
        legend.title = ggplot2::element_text(size = base_size * 0.85),
        legend.text = ggplot2::element_text(size = base_size * 0.75)
      )
  } else {
    p.mz <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = mean(rtr), y = mean(mzr), label = "No MS points", size = base_size / 2) +
      ggplot2::scale_x_continuous(limits = rtr) +
      ggplot2::scale_y_continuous(limits = mzr) +
      ggplot2::labs(x = "Retention time", y = expression(italic(m/z))) +
      ggplot2::theme_bw(base_size = base_size)
  }

  p.all <- p.chr / p.mz +
    patchwork::plot_layout(heights = c(1, 1.1)) +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.15, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = base_size * 0.95, hjust = 0.5, lineheight = 0.95)
      )
    )

  if (return.data) {
    return(list(
      plot = p.all,
      p_chr = p.chr,
      p_mz = p.mz,
      chrom = chrom.df,
      points = pts.df,
      mzr = mzr,
      rtr = rtr
    ))
  }
  p.all
}
