
saveMSdev <- function(object){
  MSdev <- object
  save(MSdev, file =  object@projectInfo$MSdevFile)
  invisible(MSdev)
}

#' @title readInRawData
#' @description read in ms raw data from `object@projectInfo$msDataDir`
#' and generate a table `sampleInfo`
#' note this function read in file according to their file names, those contaion both "pos" and "neg"
#' will be regarded as two files
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

  ### generate sampleInfo
  {
  .select_char <- function(char_vector){
    if (all(is.na(char_vector))) {
      return(NA)
    }
    max(char_vector,na.rm = T)
  }
  sample.info <- data.frame(raw.files = raw.files)%>%
    dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
    dplyr::mutate(raw.file.positive = case_when(grepl(pattern = "pos", x = raw.files,ignore.case = T)~raw.files),
                  raw.file.negative = case_when(grepl(pattern = "neg", x = raw.files,ignore.case = T)~raw.files))%>%
    dplyr::mutate(sample.abbreviation= gsub(pattern = paste0("pos|neg|",projectInfo$rawDataFormat,"$"),
                                            x = basename(raw.files) ,
                                            ignore.case = T,
                                            replacement = ""),
                  sample.abbreviation = tolower(sample.abbreviation))%>%
    dplyr::group_by(sample.abbreviation)%>%
    dplyr::mutate(raw.file.positive = .select_char (raw.file.positive),
                  raw.file.negative = .select_char (raw.file.negative))%>%
    dplyr::ungroup()%>%
    dplyr::distinct(sample.abbreviation,.keep_all = T)%>%
    dplyr::mutate(analysis.time.positive =as.character( file.info(raw.file.positive)$mtime),
                  analysis.time.negative = as.character(file.info(raw.file.negative)$mtime))%>%
    dplyr::mutate(sample.type = case_when(grepl(pattern = "QC",x = sample.abbreviation,ignore.case = T)~ "QC",
                                          grepl(pattern = "blank|blk",x = sample.abbreviation,ignore.case = T)~ "Blank",
                                          T~"Sample"),
                  .before = sample.abbreviation)%>%
    dplyr::mutate(sample.name = paste0(sample.type,str_pad(1:nrow(.), ceiling(log10(nrow(.))),pad = "0")),
                  .before = sample.type)%>%
    dplyr::mutate(msData.file.positive = case_when(is.na(raw.file.positive)~raw.file.positive,
                                                   T~paste0(msData.dir,"/pos/",sample.name,".mzML")),
                  msData.file.negative = case_when(is.na(raw.file.negative)~raw.file.positive,
                                                   T~paste0(msData.dir,"/neg/",sample.name,".mzML")))%>%
    dplyr::mutate(group = sample.type,weight = 1 ,
                  xcmsProcessing = "Both")%>%
    dplyr::select(sample.name,sample.type,group , weight,
                  sample.abbreviation,
                  raw.file.positive,raw.file.negative,
                  analysis.time.positive,analysis.time.negative,
                  msData.file.positive,msData.file.negative,
                  xcmsProcessing)

  }
  ### save
  {
    object@sampleInfo <- sample.info
    object <- .updateProjectInfoFromSampleInfo(object )
    object@processingInfo$readInRawData$done <- T


    }
  object


}

#' @title checkSampleInfo
#' @description manually check sampleInfo using excel
#' @param object a `MSdev` object
#'
#' @return a `MSdev` object
#' @export
#'
#' @examples
checkSampleInfo <- function(object){

  sampleInfo <- object@sampleInfo
  sampleInfo <- edit_df_in_excel(sampleInfo)
  ### save
  {
    object@sampleInfo <- sampleInfo
    object <- .updateProjectInfoFromSampleInfo(object )

  }

  object
}

.updateProjectInfoFromSampleInfo <- function(object){

  object@sampleInfo -> sampleInfo
  object@projectInfo$sampleCount <- sum(sampleInfo$sample.type=="Sample")
  object@projectInfo$rawDataFileCount <- sum(!is.na(sampleInfo$raw.file.positive),
                                             !is.na(sampleInfo$raw.file.negative))
  sampleInfo <- object@sampleInfo%>%
    dplyr::mutate(pos = ifelse(is.na(raw.file.positive),NA,"positive"),
                  neg = ifelse(is.na(raw.file.negative),NA,"negative"))%>%
    pivot_longer(c("pos","neg"),names_to = "ion_mode",values_to = "p")
  object@projectInfo$rawDatafiles <- table(sampleInfo$p,sampleInfo$sample.type)

  object
}


msConvert <- function(object){

  ### convert
  if (object@projectInfo$rawDataFormat == ".raw") {
    msconvert_raw2mzML(raw.files  = object@sampleInfo$raw.file.positive,
                       mzML.files = object@sampleInfo$msData.file.positive,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

    msconvert_raw2mzML(raw.files  = object@sampleInfo$raw.file.negative,
                       mzML.files = object@sampleInfo$msData.file.negative,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

  }else if(object@projectInfo$rawDataFormat == ".wiff"){
      stop("not supprot wiff file convert")

  }

  ### return
  {
    object@processingInfo$rawDataConvert <- list(
      done = T,
      time = Sys.time(),
      rawFormat =object@projectInfo$rawDataFormat,
      msDataFormat =".mzML"

    )
    saveMSdev(object )
    return(object)

  }
}

xcmsProcessing_fullscan_DDA <- function(object){

  object@xcmsData$positiveMS1  <- xcmsProcessingMS1(msDataFiles = object@sampleInfo$msData.file.positive,
                                                       ion_mode = 1,
                                                       peaksGroup =object@sampleInfo$sample.type,
                                                       centWaveParam =xcms::CentWaveParam(ppm = 5,
                                                                                          peakwidth = c(5,50))
  )
  object@xcmsData$negativeMS1  <- xcmsProcessingMS1(msDataFiles = object@sampleInfo$msData.file.negative,
                                                       ion_mode = 0,
                                                       peaksGroup =object@sampleInfo$sample.type,
                                                       centWaveParam =xcms::CentWaveParam(ppm = 5,
                                                                                          peakwidth = c(5,50))
  )

  object





}




extractSpectra_fullscan_DDA <- function(object){

  sampleInfo <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing %in% c("both","MS2"))
  spectra.pos <- Spectra::Spectra(sampleInfo$msData.file.positive)%>%
    filterMsLevel(2)%>%
    setBackend(MsBackendDataFrame())
  spectra.neg <- Spectra::Spectra(sampleInfo$msData.file.negative)%>%
    filterMsLevel(2)%>%
    setBackend(MsBackendDataFrame())

  object@spectra$positiveMS2 <- spectra.pos
  object@spectra$negativeMS2 <- spectra.neg
  return(object)

}


featureSpectra_fullscan_DDA <- function(object){

  .matchSP <- function(x,spectras,
                       mz_ppm = 10,
                       rt_tol = 10){
    mz <- x$mzmed
    rt <- x$rtmed
    mzError <- abs((mz - spectras$precursorMz)/mz*1e6)
    rtError <- abs((rt- spectras$rtime)/rt)
    matchedspectras <- which(mzError < mz_ppm &rtError < rt_tol)
    if (length(matchedspectras)==0) {
      return(NA)
    }else{
      return(matchedspectras)
    }


  }

  object@spectra$positiveFeatureMS2Map <- apply(featureDefinitions(object@xcmsData$PositiveMS1), 1, .matchSP , object@spectra$positiveMS2 )
  object@spectra$negativeFeatureMS2Map <- apply(featureDefinitions(object@xcmsData$negativeMS1), 1, .matchSP , object@spectra$negativeMS2 )
  object
}

annotateMSdev <- function(object){

  .annotateMSdev <- function(xcmsFeature , ion_mode){


  }




}
