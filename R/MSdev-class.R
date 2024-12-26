setClass("MSdev",
         slots = list(
           projectInfo = "list",
           processingInfo = "list",
           sampleInfo = "data.frame",
           experimentInfo ="MS_Exp",
           xcmsData = "list",
           spectra = "list",
           annotation = "list",
           statData = "list"
           ))




#' @importFrom magrittr `%>%` `%<>%`
#' @import  tidyverse dplyr ggplot2 patchwork
#' @importFrom xcms filepaths mzrange plotSpec noise filterFile dirname chromPeaks
#' featureDefinitions polarity `featureDefinitions<-` featureValues featureChromatograms
#' @importFrom MSnbase readMSData sampleNames fileNames
#' @importFrom BiocParallel SnowParam SerialParam bpmapply bplapply snowWorkers
#' @importFrom tibble tibble as.tibble rownames_to_column remove_rownames
#' @importFrom S4Vectors isEmpty
#' @importFrom Biobase pData fData
#' @importFrom Spectra msLevel Spectra MsBackendMemory `spectraNames<-`
#' peaksData spectraNames filterPolarity rtime dataOrigin
#' @importFrom ProtGenerics filterMsLevel filterPolarity collisionEnergy
#' @importFrom CompoundDb compounds
#' @importFrom ChemmineR atomcountMA MF `cid<-` cid validSDF canonicalize
#' atomblock bonds bondblock rings
#' @importFrom MSCC chemform_adduct_check
#' @importFrom stringr str_extract
#' @importFrom plotly layout
#' @importFrom SummarizedExperiment rowData colData assay `rowData<-` `colData<-` `assay<-`
#' @importFrom shiny shinyApp fluidPage fluidRow
#' checkboxInput selectInput navbarPage
#' column h1 h3 h4 br wellPanel
#' @importFrom plotly add_markers plot_ly add_pie add_segments
#' event_data add_text add_lines
#' @importFrom visNetwork visNetwork visLayout
#' @importFrom igraph as_adjacency_matrix distances edge.attributes
#' degree V E
#' @importFrom ComplexHeatmap draw Legend rowAnnotation
#' @importFrom grid gpar grid.rect grid.circle
MSdev <- function(rawDataDir =
                    "C:/Users/91879/OneDrive/Code/R/Projecct/2022.1.8_MS.demo/Demo/raw.data",
                  projectDir = dirname(rawDataDir),
                  experimentInfo = MS_Exp()){
  .Object <- new("MSdev")
  {

    ### check param
    {
      if (!dir.exists(rawDataDir)) {
        message("rawDataDir do not exist")

      }
      rawDataDir <- normalizePath(rawDataDir)
      projectDir <- normalizePath(projectDir)
      dir.create(projectDir,recursive = T,showWarnings = F)
    }

    .Object@projectInfo = list(
      creatTime = Sys.time(),
      projectDir =projectDir,
      rawDataDir =rawDataDir,
      rawDataFormat = NULL,
      rawDataFileCount = NULL,
      msDataDir =paste0(rawDataDir , "/msData"),
      experimentTime = NULL,
      msAcquisition = "Unknow",
      msLevel= NA,
      polarity = -1,
      sampleCount = NA,
      MSdevFile = paste0(projectDir , "/MSdev",date_suffix(),".Rdata")
    )
    .Object@experimentInfo = experimentInfo

    ### check raw data
    {
      rowDataFile <- dir(.Object@projectInfo$rawDataDir)
      .Object@projectInfo$rawDataFormat =
        dplyr::case_when(any(grepl(pattern = ".wiff$", x = rowDataFile))~".wiff",
                         any(grepl(pattern = ".raw$", x = rowDataFile))~".raw",
                         any(grepl(pattern = ".lcd$", x = rowDataFile))~".lcd")
      rowDataFile <- dir(path = .Object@projectInfo$rawDataDir,
                         pattern =paste0(.Object@projectInfo$rawDataFormat,"$") ,
                         full.names = T)
      .Object@projectInfo$rawDataFileCount <- length(rowDataFile)
      .Object@projectInfo$experimentTime <- max(file.info(rowDataFile)$mtime)

    }
    .Object@sampleInfo <- get_MS_sampleinfo(.Object@projectInfo$rawDataDir,
                                            rawDataFormat = .Object@projectInfo$rawDataFormat,
                                            verbose = F)
    .Object <- .updateProjectInfoFromSampleInfo(.Object )
    .Object@processingInfo$readInRawData$done <- T
    MSdev_save(.Object)
    .Object

  }

}





setMethod(
  "show" ,
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

