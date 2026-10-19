#' @title Resolve an MSdev project path
#' @description Return an absolute path for a key in \code{projectInfo$paths}.
#'   Paths are stored relative to \code{projectDir}. Old objects that have no
#'   \code{paths} entry are resolved from the existing absolute fields.
#' @param object MSdev object
#' @param key One of \code{rawData}, \code{msData}, \code{msdevFile},
#'   \code{ms1}, \code{ms2}.
#' @return Absolute path, with forward slashes.
#' @export
MSdev_path <- function(object, key) {
  stopifnot(inherits(object, "MSdev"))
  key <- as.character(key)[[1L]]
  known <- c("rawData", "msData", "msdevFile", "ms1", "ms2")
  if (!key %in% known) {
    stop(
      "Unknown MSdev path key: ", key,
      ". Expected one of: ", paste(known, collapse = ", "),
      call. = FALSE
    )
  }
  rel <- NULL
  paths <- object@projectInfo$paths
  if (!is.null(paths) && !is.null(paths[[key]]) && length(paths[[key]]) &&
      nzchar(as.character(paths[[key]])[[1L]])) {
    rel <- as.character(paths[[key]])[[1L]]
  } else {
    rel <- .msdev_default_relative(object, key)
  }
  .msdev_join_root(object@projectInfo$projectDir, rel)
}


#' @title Point an MSdev object at a new project directory
#' @description Rewrite absolute paths that sit under the saved
#'   \code{projectDir} so they sit under \code{projectDir}. Sample files,
#'   \code{onDiskData} locations, Spectra \code{dataOrigin} / \code{dataStorage},
#'   and xcms experiment file paths are updated. Files are not moved; the
#'   project tree is assumed to already be at the new location.
#'
#'   \code{rawDataDir} overrides the raw-data directory when it is not inside
#'   the project directory. \code{compoundDb} sets
#'   \code{projectInfo$CompoundDB_path} without prefix replacement.
#'   After rewriting, \code{projectInfo$paths} holds paths relative to
#'   \code{projectDir}.
#' @param object MSdev object
#' @param projectDir New project directory.
#' @param rawDataDir Optional new raw-data directory. Used when raw data does
#'   not live under \code{projectDir}.
#' @param compoundDb Optional CompoundDb SQLite file for this machine.
#' @return The same MSdev object with updated paths.
#' @export
MSdev_relocate <- function(object,
                           projectDir,
                           rawDataDir = NULL,
                           compoundDb = NULL) {
  stopifnot(inherits(object, "MSdev"))
  if (missing(projectDir) || !length(projectDir) ||
      !nzchar(as.character(projectDir)[[1L]])) {
    stop("'projectDir' is required.", call. = FALSE)
  }

  old_project <- .msdev_path_prefix(object@projectInfo$projectDir)
  new_project <- .msdev_path_prefix(projectDir)
  if (!nzchar(new_project)) {
    stop("'projectDir' is empty.", call. = FALSE)
  }
  if (!dir.exists(new_project)) {
    dir.create(new_project, recursive = TRUE, showWarnings = FALSE)
  }
  new_project <- .msdev_path_prefix(new_project)

  maps <- list()
  if (!is.null(rawDataDir) && length(rawDataDir) &&
      nzchar(as.character(rawDataDir)[[1L]])) {
    old_raw <- .msdev_path_prefix(object@projectInfo$rawDataDir)
    new_raw <- .msdev_path_prefix(rawDataDir)
    if (nzchar(old_raw) && nzchar(new_raw)) {
      if (!dir.exists(new_raw)) {
        dir.create(new_raw, recursive = TRUE, showWarnings = FALSE)
      }
      maps[[length(maps) + 1L]] <- list(old = old_raw, new = .msdev_path_prefix(new_raw))
    }
  }
  if (nzchar(old_project)) {
    maps[[length(maps) + 1L]] <- list(old = old_project, new = new_project)
  }

  if (length(maps) && !.msdev_maps_noop(maps)) {
    for (key in c("projectDir", "rawDataDir", "msDataDir", "MSdevFile")) {
      val <- object@projectInfo[[key]]
      if (is.character(val) && length(val) == 1L && !is.na(val) && nzchar(val)) {
        object@projectInfo[[key]] <- .msdev_rewrite_paths(val, maps)
      }
    }
    object@projectInfo$projectDir <- new_project

    if (is.data.frame(object@sampleInfo) && ncol(object@sampleInfo)) {
      for (col in intersect(c("raw.files", "msData.files"), names(object@sampleInfo))) {
        object@sampleInfo[[col]] <- .msdev_rewrite_column(object@sampleInfo[[col]], maps)
      }
    }

    if (is.list(object@spectra)) {
      object@spectra <- .msdev_rewrite_ondisk_tree(object@spectra, maps)
    }
    if (is.list(object@xcmsData)) {
      object@xcmsData <- .msdev_rewrite_ondisk_tree(object@xcmsData, maps)
      for (i in seq_along(object@xcmsData)) {
        object@xcmsData[[i]] <- .msdev_rewrite_xcms(object@xcmsData[[i]], maps)
      }
    }
  } else {
    object@projectInfo$projectDir <- new_project
  }

  if (!is.null(compoundDb) && length(compoundDb) &&
      nzchar(as.character(compoundDb)[[1L]])) {
    object@projectInfo$CompoundDB_path <- .msdev_path_prefix(compoundDb)
  }

  object <- .msdev_fill_paths(object)
  message_with_time("Relocated projectDir: ", old_project, " -> ", object@projectInfo$projectDir)
  object
}


.msdev_slash <- function(x) {
  if (is.null(x) || !length(x)) {
    return(x)
  }
  gsub("\\\\", "/", as.character(x), fixed = FALSE)
}

.msdev_path_prefix <- function(path) {
  if (is.null(path) || !length(path) || is.na(path[[1L]]) ||
      !nzchar(as.character(path[[1L]]))) {
    return("")
  }
  path <- tryCatch(
    normalizePath(path[[1L]], winslash = "/", mustWork = FALSE),
    error = function(e) as.character(path[[1L]])
  )
  path <- gsub("\\\\", "/", path)
  sub("/+$", "", path)
}

.msdev_same_path <- function(a, b) {
  a <- .msdev_path_prefix(a)
  b <- .msdev_path_prefix(b)
  if (!nzchar(a) || !nzchar(b)) {
    return(FALSE)
  }
  if (.Platform$OS.type == "windows") {
    return(tolower(a) == tolower(b))
  }
  identical(a, b)
}

.msdev_is_absolute <- function(path) {
  path <- .msdev_slash(path)
  grepl("^(/|[A-Za-z]:/)", path)
}

.msdev_under_root <- function(path, root) {
  path <- .msdev_slash(path)
  root <- .msdev_path_prefix(root)
  if (!nzchar(root) || !length(path)) {
    return(logical(length(path)))
  }
  path_cmp <- path
  root_cmp <- root
  if (.Platform$OS.type == "windows") {
    path_cmp <- tolower(path)
    root_cmp <- tolower(root)
  }
  path_cmp == root_cmp | startsWith(path_cmp, paste0(root_cmp, "/"))
}

.msdev_relative <- function(path, root) {
  if (is.null(path) || !length(path) || is.na(path[[1L]]) ||
      !nzchar(as.character(path[[1L]]))) {
    return("")
  }
  path_n <- .msdev_path_prefix(path)
  root_n <- .msdev_path_prefix(root)
  if (!nzchar(path_n) || !nzchar(root_n) || !.msdev_under_root(path_n, root_n)) {
    return(path_n)
  }
  if (.msdev_same_path(path_n, root_n)) {
    return(".")
  }
  substring(path_n, nchar(root_n) + 2L)
}

.msdev_join_root <- function(root, rel) {
  root <- .msdev_path_prefix(root)
  if (is.null(rel) || !length(rel) || is.na(rel[[1L]]) || !nzchar(rel[[1L]]) ||
      identical(rel[[1L]], ".")) {
    return(root)
  }
  rel <- .msdev_slash(rel[[1L]])
  if (.msdev_is_absolute(rel)) {
    return(.msdev_path_prefix(rel))
  }
  .msdev_slash(file.path(root, rel))
}

.msdev_default_relative <- function(object, key) {
  projectDir <- object@projectInfo$projectDir
  abs <- switch(
    key,
    rawData = object@projectInfo$rawDataDir,
    msData = object@projectInfo$msDataDir,
    msdevFile = object@projectInfo$MSdevFile,
    ms1 = {
      p <- object@spectra$MS1_Spectra
      if (inherits(p, "onDiskData") && length(p@path) && nzchar(p@path)) {
        p@path
      } else {
        file.path(projectDir, "MS1_Spectra.rds")
      }
    },
    ms2 = {
      p <- object@spectra$MS2_Spectra
      if (inherits(p, "onDiskData") && length(p@path) && nzchar(p@path)) {
        p@path
      } else {
        file.path(projectDir, "MS2_Spectra.rds")
      }
    },
    NULL
  )
  if (is.null(abs) || !length(abs) || is.na(abs[[1L]]) || !nzchar(as.character(abs[[1L]]))) {
    stop("MSdev object has no path for '", key, "'.", call. = FALSE)
  }
  .msdev_relative(abs[[1L]], projectDir)
}

.msdev_fill_paths <- function(object) {
  projectDir <- .msdev_path_prefix(object@projectInfo$projectDir)
  object@projectInfo$projectDir <- projectDir
  ms1 <- NULL
  ms2 <- NULL
  if (is.list(object@spectra)) {
    ms1 <- object@spectra$MS1_Spectra
    ms2 <- object@spectra$MS2_Spectra
  }
  abs_map <- c(
    rawData = .msdev_chr1(object@projectInfo$rawDataDir),
    msData = .msdev_chr1(object@projectInfo$msDataDir),
    msdevFile = .msdev_chr1(object@projectInfo$MSdevFile),
    ms1 = if (inherits(ms1, "onDiskData") && length(ms1@path) && nzchar(ms1@path)) {
      ms1@path
    } else {
      file.path(projectDir, "MS1_Spectra.rds")
    },
    ms2 = if (inherits(ms2, "onDiskData") && length(ms2@path) && nzchar(ms2@path)) {
      ms2@path
    } else {
      file.path(projectDir, "MS2_Spectra.rds")
    }
  )
  object@projectInfo$paths <- vapply(
    abs_map,
    function(p) .msdev_relative(p, projectDir),
    character(1),
    USE.NAMES = TRUE
  )
  object
}

.msdev_chr1 <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) {
    return("")
  }
  as.character(x[[1L]])
}

.msdev_rewrite_one <- function(path, maps) {
  path_s <- .msdev_slash(path)
  path_n <- .msdev_path_prefix(path_s)
  if (!nzchar(path_n)) {
    return(path_s)
  }
  for (m in maps) {
    if (.msdev_under_root(path_n, m$old)) {
      return(.msdev_join_root(m$new, .msdev_relative(path_n, m$old)))
    }
  }
  path_n
}

.msdev_rewrite_paths <- function(paths, maps) {
  if (is.null(paths)) {
    return(paths)
  }
  if (is.factor(paths)) {
    paths <- as.character(paths)
  }
  if (!is.character(paths)) {
    return(paths)
  }
  if (!length(paths)) {
    return(paths)
  }
  nm <- names(paths)
  out <- vapply(paths, function(p) {
    if (is.na(p) || !nzchar(p)) {
      return(if (is.na(p)) NA_character_ else p)
    }
    .msdev_rewrite_one(p, maps)
  }, character(1), USE.NAMES = FALSE)
  names(out) <- nm
  out
}

.msdev_rewrite_column <- function(x, maps) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (!is.character(x)) {
    return(x)
  }
  .msdev_rewrite_paths(x, maps)
}

.msdev_maps_noop <- function(maps) {
  if (!length(maps)) {
    return(TRUE)
  }
  all(vapply(maps, function(m) .msdev_same_path(m$old, m$new), logical(1)))
}

.msdev_rewrite_spectra_paths <- function(sp, maps) {
  if (!inherits(sp, "Spectra") || !length(sp)) {
    return(sp)
  }
  ## dataStorage<- / dataOrigin<- are methods for MsBackend*, not Spectra.
  ## Spectra::$<- writes those spectra variables on either backend.
  ds <- tryCatch(Spectra::dataStorage(sp), error = function(e) NULL)
  if (is.character(ds) && length(ds) == length(sp)) {
    new_ds <- .msdev_rewrite_paths(ds, maps)
    if (!identical(as.character(ds), as.character(new_ds))) {
      sp$dataStorage <- new_ds
    }
  }
  orig <- tryCatch(Spectra::dataOrigin(sp), error = function(e) NULL)
  if (is.character(orig) && length(orig) == length(sp)) {
    new_orig <- .msdev_rewrite_paths(orig, maps)
    if (!identical(as.character(orig), as.character(new_orig))) {
      sp$dataOrigin <- new_orig
    }
  }
  sp
}

.msdev_spectra_paths_equal <- function(a, b) {
  if (!inherits(a, "Spectra") || !inherits(b, "Spectra") || length(a) != length(b)) {
    return(FALSE)
  }
  ds_ok <- identical(
    as.character(Spectra::dataStorage(a)),
    as.character(Spectra::dataStorage(b))
  )
  orig_ok <- identical(
    as.character(Spectra::dataOrigin(a)),
    as.character(Spectra::dataOrigin(b))
  )
  ds_ok && orig_ok
}

.msdev_stored_changed <- function(a, b) {
  if (inherits(a, "Spectra") && inherits(b, "Spectra")) {
    return(!.msdev_spectra_paths_equal(a, b))
  }
  if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.data.frame(b) &&
      length(a) == length(b)) {
    return(any(vapply(
      seq_along(a),
      function(i) .msdev_stored_changed(a[[i]], b[[i]]),
      logical(1)
    )))
  }
  FALSE
}

.msdev_rewrite_stored <- function(data, maps) {
  if (inherits(data, "Spectra")) {
    return(.msdev_rewrite_spectra_paths(data, maps))
  }
  if (is.list(data) && !is.data.frame(data)) {
    data[] <- lapply(data, .msdev_rewrite_stored, maps = maps)
    return(data)
  }
  data
}

.msdev_rewrite_ondisk <- function(object, maps) {
  if (!inherits(object, "onDiskData")) {
    return(object)
  }
  new_path <- unname(.msdev_rewrite_paths(object@path, maps))
  object@path <- new_path
  if (!length(new_path) || is.na(new_path) || !nzchar(new_path) || !file.exists(new_path)) {
    if (length(new_path) && !is.na(new_path) && nzchar(new_path) && !file.exists(new_path)) {
      message_with_time("onDiskData file is not at the new path: ", new_path)
    }
    return(object)
  }
  data <- onDiskData_retrieve(object)
  data2 <- .msdev_rewrite_stored(data, maps)
  if (.msdev_stored_changed(data, data2)) {
    object <- onDiskData_update(object, data2)
  }
  object
}

.msdev_rewrite_ondisk_tree <- function(x, maps) {
  if (inherits(x, "onDiskData")) {
    return(.msdev_rewrite_ondisk(x, maps))
  }
  if (inherits(x, "MsExperiment") || inherits(x, "XCMSnExp") || inherits(x, "OnDiskMSnExp")) {
    return(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    x[] <- lapply(x, .msdev_rewrite_ondisk_tree, maps = maps)
    return(x)
  }
  x
}

.msdev_rewrite_sample_table <- function(object, maps) {
  sd <- tryCatch(MsExperiment::sampleData(object), error = function(e) NULL)
  if (is.null(sd) || !nrow(sd)) {
    return(object)
  }
  cols <- intersect(
    c("raw.files", "msData.files", "dataOrigin", "dataStorage"),
    colnames(sd)
  )
  if (!length(cols)) {
    return(object)
  }
  for (col in cols) {
    rewritten <- .msdev_rewrite_column(as.character(sd[[col]]), maps)
    if (length(rewritten) == nrow(sd)) {
      sd[[col]] <- rewritten
    }
  }
  MsExperiment::`sampleData<-`(object, value = sd)
}

.msdev_rewrite_mse_spectra <- function(object, maps) {
  sp <- tryCatch(ProtGenerics::spectra(object), error = function(e) NULL)
  if (is.null(sp) || !inherits(sp, "Spectra") || !length(sp)) {
    return(object)
  }
  sp <- .msdev_rewrite_spectra_paths(sp, maps)
  updated <- tryCatch(
    ProtGenerics::`spectra<-`(object, value = sp),
    error = function(e) e
  )
  if (inherits(updated, "error")) {
    if (methods::.hasSlot(object, "spectra")) {
      object@spectra <- sp
      return(object)
    }
    warning("Could not update spectra paths: ", conditionMessage(updated), call. = FALSE)
    return(object)
  }
  updated
}

.msdev_rewrite_experiment_files <- function(object, maps) {
  if (!methods::.hasSlot(object, "experimentFiles")) {
    return(object)
  }
  ef <- object@experimentFiles
  if (is.null(ef) || !length(ef)) {
    return(object)
  }
  nms <- names(ef)
  if (is.null(nms)) {
    nms <- as.character(seq_along(ef))
  }
  tryCatch({
    for (nm in nms) {
      val <- ef[[nm]]
      if (is.character(val)) {
        ef[[nm]] <- .msdev_rewrite_paths(val, maps)
      }
    }
    object@experimentFiles <- ef
    object
  }, error = function(e) object)
}

.msdev_rewrite_msnexp_files <- function(object, maps) {
  if (methods::.hasSlot(object, "processingData")) {
    pd <- object@processingData
    if (methods::.hasSlot(pd, "files") && length(pd@files)) {
      pd@files <- .msdev_rewrite_paths(pd@files, maps)
      object@processingData <- pd
    }
  }
  pdata <- tryCatch(Biobase::pData(object), error = function(e) NULL)
  if (is.data.frame(pdata) && ncol(pdata)) {
    cols <- intersect(
      c("raw.files", "msData.files", "dataOrigin", "dataStorage"),
      names(pdata)
    )
    if (length(cols)) {
      for (col in cols) {
        pdata[[col]] <- .msdev_rewrite_column(pdata[[col]], maps)
      }
      object <- tryCatch(
        Biobase::`pData<-`(object, value = pdata),
        error = function(e) object
      )
    }
  }
  object
}

.msdev_rewrite_xcms <- function(object, maps) {
  if (is.null(object) || identical(object, NA)) {
    return(object)
  }
  if (inherits(object, "MsExperiment")) {
    object <- .msdev_rewrite_mse_spectra(object, maps)
    object <- .msdev_rewrite_sample_table(object, maps)
    object <- .msdev_rewrite_experiment_files(object, maps)
    return(object)
  }
  if (inherits(object, "XCMSnExp") || inherits(object, "OnDiskMSnExp")) {
    return(.msdev_rewrite_msnexp_files(object, maps))
  }
  object
}

.msdev_load_relocate <- function(object, file_to_load) {
  saved <- object@projectInfo$projectDir
  saved_ok <- is.character(saved) && length(saved) == 1L && !is.na(saved) &&
    nzchar(saved) && dir.exists(saved)
  if (isTRUE(saved_ok)) {
    return(object)
  }
  loaded_dir <- tryCatch(
    .msdev_path_prefix(dirname(file_to_load)),
    error = function(e) ""
  )
  if (!nzchar(loaded_dir) || .msdev_same_path(saved, loaded_dir)) {
    return(object)
  }
  MSdev_relocate(object, projectDir = loaded_dir)
}

.msdev_resolve_cpdb_path <- function(object = NULL, cpdb_path = NULL) {
  candidates <- character()
  add <- function(x) {
    if (is.null(x) || !length(x) || is.na(x[[1L]])) {
      return()
    }
    x <- as.character(x)[[1L]]
    if (nzchar(x)) {
      candidates <<- c(candidates, x)
    }
  }
  add(cpdb_path)
  if (!is.null(object) && inherits(object, "MSdev")) {
    add(object@projectInfo$CompoundDB_path)
  }
  add(getOption("MSdev.compoundDb", default = NULL))
  env <- Sys.getenv("MSDEV_COMPOUNDDB", unset = "")
  add(env)
  candidates <- unique(candidates)
  for (p in candidates) {
    if (file.exists(p)) {
      return(.msdev_path_prefix(p))
    }
  }
  NULL
}

.msdev_require_cpdb_path <- function(object = NULL, cpdb_path = NULL) {
  resolved <- .msdev_resolve_cpdb_path(object, cpdb_path)
  if (is.null(resolved)) {
    stop(
      "CompoundDb path not found. Pass `cpdb_path`, set options(MSdev.compoundDb = ...), or set the MSDEV_COMPOUNDDB environment variable.",
      call. = FALSE
    )
  }
  resolved
}
