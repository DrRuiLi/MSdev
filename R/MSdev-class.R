setClass("MSdev",
         slots = list(
           projectInfo = "list",
           processingInfo = "list",
           sampleInfo = "data.frame",
           experimentInfo ="MS_Exp",
           xcmsData = "list",
           spectra = "list",
           annotation = "list",
           advancedAna = "list"
           ))




#' @importFrom magrittr `%>%` `%<>%`
#' @import dplyr ggplot2 patchwork
#' @importFrom BiocParallel SnowParam SerialParam bpmapply bplapply snowWorkers register
#' @importFrom tibble tibble as.tibble rownames_to_column remove_rownames
#' @importFrom S4Vectors isEmpty
#' @importFrom Biobase pData fData
#' @import ProtGenerics
#' @importFrom CompoundDb compounds
#' @importFrom stringr str_extract str_split str_extract_all
#' @importFrom plotly layout add_bars
#' @importFrom SummarizedExperiment rowData colData assay `rowData<-` `colData<-` `assay<-`
#' @importFrom plotly add_markers plot_ly add_pie add_segments event_data add_text add_lines
#' @importFrom visNetwork visNetwork visLayout visOptions visInteraction visEvents visEdges visGetEdges visNetworkProxy visUpdateEdges
#' @importFrom igraph as_adjacency_matrix distances edge.attributes degree V E `V<-` `E<-`
#' @importFrom ComplexHeatmap draw Legend rowAnnotation
#' @importFrom grid gpar grid.rect grid.circle
#' @importFrom data.table rbindlist data.table
#' @importFrom readr write_csv

#' @title Create an MSdev object
#' @description Create a new MSdev object for mass spectrometry data analysis.
#' This function initializes the project structure, detects raw data format,
#' and creates sample information from the raw data directory.
#'
#' @param rawDataDir path to directory containing raw mass spectrometry data files.
#'   Required. Supported extensions: \code{.wiff}, \code{.raw}, \code{.mzXML},
#'   \code{.mzML}, \code{.lcd}; matching is case-insensitive.
#' @param projectDir project directory for storing processed data, defaults to parent of rawDataDir.
#'   This is the single relocatable project root. Other project paths are also
#'   stored relative to it in \code{projectInfo$paths}.
#' @param experimentInfo MS_Exp object containing experiment metadata
#'
#' @return MSdev object
#' @export
MSdev <- function(rawDataDir,
                  projectDir = dirname(rawDataDir),
                  experimentInfo = MS_Exp()){
  .Object <- new("MSdev")
  {

    ### check param
    {
      if (missing(rawDataDir) || !length(rawDataDir) ||
          !nzchar(as.character(rawDataDir)[[1L]])) {
        stop("'rawDataDir' is required.", call. = FALSE)
      }
      if (!dir.exists(rawDataDir)) {
        message("rawDataDir do not exist")
      }
      rawDataDir <- .msdev_path_prefix(rawDataDir)
      if (!dir.exists(projectDir)) {
        dir.create(projectDir, recursive = TRUE, showWarnings = FALSE)
      }
      projectDir <- .msdev_path_prefix(projectDir)
      msDataDir <- .msdev_path_prefix(file.path(rawDataDir, "msData"))
      msdevFile <- .msdev_path_prefix(
        file.path(projectDir, paste0("MSdev", date_suffix(), ".Rdata"))
      )
    }

    .Object@projectInfo = list(
      creatTime = Sys.time(),
      projectDir = projectDir,
      rawDataDir = rawDataDir,
      rawDataFormat = NULL,
      rawDataFileCount = NULL,
      msDataDir = msDataDir,
      experimentTime = NULL,
      msAcquisition = "Unknow",
      msLevel= NA,
      polarity = -1,
      sampleCount = NA,
      MSdevFile = msdevFile
    )
    .Object <- .msdev_fill_paths(.Object)
    .Object@experimentInfo = experimentInfo

    ### check raw data
    {
      rowDataFile <- dir(.Object@projectInfo$rawDataDir, recursive = TRUE)
      .Object@projectInfo$rawDataFormat <- dplyr::case_when(
        any(grepl("\\.wiff$", rowDataFile, ignore.case = TRUE)) ~ ".wiff",
        any(grepl("\\.raw$", rowDataFile, ignore.case = TRUE)) ~ ".raw",
        any(grepl("\\.mzXML$", rowDataFile, ignore.case = TRUE)) ~ ".mzXML",
        any(grepl("\\.mzML$", rowDataFile, ignore.case = TRUE)) ~ ".mzML",
        any(grepl("\\.lcd$", rowDataFile, ignore.case = TRUE)) ~ ".lcd",
        TRUE ~ NA_character_
      )
      if (is.na(.Object@projectInfo$rawDataFormat)) {
        stop(
          "No supported raw data files found under ",
          .Object@projectInfo$rawDataDir,
          " (expected .wiff, .raw, .mzXML, .mzML, or .lcd).",
          call. = FALSE
        )
      }
      rowDataFile <- dir(
        path = .Object@projectInfo$rawDataDir,
        full.names = TRUE,
        recursive = TRUE
      )
      rowDataFile <- rowDataFile[
        grepl(
          paste0("\\", .Object@projectInfo$rawDataFormat, "$"),
          rowDataFile,
          ignore.case = TRUE
        )
      ]
      .Object@projectInfo$rawDataFileCount <- length(rowDataFile)
      .Object@projectInfo$experimentTime <- max(file.info(rowDataFile)$mtime)

    }
    .Object@sampleInfo <- get_MS_sampleinfo(.Object@projectInfo$rawDataDir,
                                            rawDataFormat = .Object@projectInfo$rawDataFormat,
                                            verbose = F)
    .Object <- MSdev_get_MSinfo(.Object)
    .Object <- .updateProjectInfoFromSampleInfo(.Object )
    .Object@processingInfo$readInRawData$done <- T
    MSdev_save(.Object)
    .Object

  }

}





setMethod("show" ,
  "MSdev",
  definition = function(object) {
    project_info <- object@projectInfo
    cat(paste0("MSdev project"))

    if (isTRUE(object@processingInfo$readInRawData$done)) {

      cat(paste0("\n",
                          'Total ',
                          project_info$sampleCount," Samples, ",
                          project_info$rawDataFileCount," *",project_info$rawDataFormat," files"))

      print(object@projectInfo$rawDatafiles  )

    }




  }
)


#' @title List process history of xcmsData
#' @description List all process history entries for a specific xcmsData element in an MSdev object.
#' @param object MSdev object
#' @param target character. Name of the xcmsData element (default "PositiveMS1").
#' @return data.frame of process history, or NULL if not available.
#' @export
#'
MSdev_processInfo <- function(object, target = "PositiveMS1") {

  xcms.xcms <- object@xcmsData[[target]]
  if (is.null(xcms.xcms) || identical(xcms.xcms, NA)  ) {
    message("xcmsData$", target, " is NULL or NA")
    return(NULL)
  }

  ph <- xcms::processHistory(xcms.xcms)
  if (length(ph) == 0) {
    message("No process history for xcmsData$", target)
    return(NULL)
  }

 php <- lapply(seq_along(ph), function(i) {
    p <- ph[[i]]
    #print(xcms::processType(p))
    p.param <- xcms::processParam(p)
    if (xcms::processType(p) == "Retention time correction") {
      p.param@peakGroupsMatrix <- matrix()
    }
    #print(p.param)
    return(p.param)
  })
 names(php) <- lapply(seq_along(ph), function(i) {
   xcms::processType( ph[[i]])
 })
 return(php)

}

