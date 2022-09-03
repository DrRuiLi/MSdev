
saveMSdev <- function(object){
  MSdev <- object
  save(MSdev, file =  object@projectInfo$MSdevFile)
  invisible(MSdev)
}

#' @title readInRawData
#' @description read in ms raw data from `object@projectInfo$msDataDir`
#' and generate a table `sampleInfo`
#' note this function read in file according to their file names
#' @param object a `MSdev` object
#' @details discriminate sample type and ion mode according char in file names
#'
#' grep "pos" and "neg" for ion mode
#'
#' grep "blank", "blk" and "QC" for sample type, other samples will regard as "Sample"
#' @return a `MSdev` object
#' @export
#'
#' @examples
readInRawData <- function(object){

  projectInfo <- object@projectInfo
  msData.dir <- projectInfo$msDataDir
  raw.files <- dir(path = projectInfo$rawDataDir,
                     pattern = paste0(projectInfo$rawDataFormat,"$"),
                     full.names = T)
  sample.info <- data.frame(raw.files = raw.files)%>%
    dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
    dplyr::mutate(polarity = case_when(grepl(pattern = "pos", x = raw.files,ignore.case = T)~"raw.file.positive",
                                       grepl(pattern = "neg", x = raw.files,ignore.case = T)~"raw.file.negative",
                                       T~"error"),
                  sample.abbreviation= gsub(pattern = paste0("pos|neg|",projectInfo$rawDataFormat,"$"),
                                            x = basename(raw.files) ,
                                            ignore.case = T,
                                            replacement = ""),
                  sample.abbreviation = tolower(sample.abbreviation))%>%
    pivot_wider(names_from = polarity , values_from = raw.files)%>%
    dplyr::mutate(analysis.time.positive = file.info(raw.file.positive)$mtime,
                  analysis.time.negative = file.info(raw.file.negative)$mtime)%>%
    dplyr::mutate(sample.type = case_when(grepl(pattern = "QC",x = sample.abbreviation,ignore.case = T)~ "QC",
                                          grepl(pattern = "blank|blk",x = sample.abbreviation,ignore.case = T)~ "Blank",
                                          T~"Sample"),
                  .before = sample.abbreviation)%>%
    dplyr::mutate(sample.name = paste0(sample.type,str_pad(1:nrow(.), ceiling(log10(nrow(.))),pad = "0")),
                  .before = sample.type)%>%
    dplyr::mutate(msData.file.positive = case_when(is.na(raw.file.positive)~raw.file.positive,
                                                 T~paste0(msData.dir,"pos/",sample.name,".mzML")),
                  msData.file.negative = case_when(is.na(raw.file.negative)~raw.file.positive,
                                                 T~paste0(msData.dir,"neg/",sample.name,".mzML")))%>%
    dplyr::select(sample.name,sample.type,sample.abbreviation,
                  raw.file.positive,raw.file.negative,
                  analysis.time.positive,analysis.time.negative,
                  msData.file.positive,msData.file.negative)

  ### save
  {
    object@sampleInfo <- sample.info
    object@processingInfo$readInRawData$done <- T
    object@projectInfo$sampleCount <- sum(sample.info$sample.type=="Sample")


    }
  object


}




