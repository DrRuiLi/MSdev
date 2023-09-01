setClass("MSU",
                slots = list(
                    projectInfo = "list",
                    sampleInfo = "data.frame",
                    experimentInfo = "MS_Exp",
                    xcmsData = "list",
                    spectra = "list",
                    stat = "list"
                  )
                )



MSU <- function(rawDataDir =
                    "C:/Users/91879/OneDrive/Documents/Code/R/Projecct/2022.1.8_MS.demo/Demo/raw.data",
                  projectDir = dirname(rawDataDir),
                  experimentInfo = MS_Exp()){
  .Object <- new("MSU")


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
      msDataDir =paste0(projectDir , "/msData"),
      experimentTime = NULL,
      msAcquisition = "DDA",
      sampleCount = NA,
      MSUFile = paste0(projectDir , "/MSU",date_suffix(),".Rdata")
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
                                            .Object@projectInfo$rawDataFormat)
    .Object@projectInfo$sampleCount <- length(unique(.Object@sampleInfo$sample.name))
    saveMSU(.Object)
    .Object



}


setMethod(
  "show" ,
  "MSU",
  definition = function(object) {
    project_info <- object@projectInfo
    cat(paste0("MSU project"))

    cat(paste0("\n",
                 'Total ',
                 project_info$sampleCount," Samples, ",
                 project_info$rawDataFileCount," *",
                 project_info$rawDataFormat," files"))

      print(object@projectInfo$rawDatafiles  )





  }
)









