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



MSdev <- function(rawDataDir =
                    "C:/Users/91879/OneDrive/Documents/Code/R/Projecct/2022.1.8_MS.demo/Demo/raw.data",
                  projectDir =
                    "C:/Users/91879/OneDrive/Documents/Code/R/Projecct/2022.1.8_MS.demo/Demo3",
                  experimentInfo = MS_Exp()){
  new("MSdev",rawDataDir,projectDir,experimentInfo)

}



setMethod("initialize" , "MSdev",
          function(.Object,
                   rawDataDir,
                   projectDir,
                   experimentInfo
                   ){

            ### check param
            {
              if (!dir.exists(rawDataDir)) {
                message("rawDataDir do not exist")

              }
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
              msAcquisition = "unknown",
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
            .Object <- readInRawData(.Object)
            saveMSdev(.Object)
            .Object

          })



setMethod(
  "show" ,
  "MSdev",
  definition = function(object) {
    project_info <- object@projectInfo
    show_info <-paste0("MSdev project")

    if (isTRUE(object@processingInfo$readInRawData$done)) {

      show_info <- paste0(show_info,"\n",
                          'Total ',
                          project_info$sampleCount," Samples, ",
                          project_info$rawDataFileCount," *",project_info$rawDataFormat," files")
    }
    cat(show_info)


  }
)

