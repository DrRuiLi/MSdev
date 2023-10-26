
#' check SampleInfo in excel
#' @description manually check sampleInfo using excel
#' @param object a `MSdev` object
#'
#' @return a `MSdev` object
#' @export
#'
#' @examples
MSdev_checkSampleInfo <- function(object){

  sampleInfo <- object@sampleInfo
  sampleInfo <- edit_df_in_excel(sampleInfo)
  ### save
  {
    object@sampleInfo <- sampleInfo
    object <- .updateProjectInfoFromSampleInfo(object )

  }

  object
}




#' @title msConvert_MSdev
#'
#' @param object
#'
#' @return
#' @export
#' @importFrom BiocParallel  bplapply
#' @examples
#'

MSdev_msConvert<- function(object){


  ### filter files
  {
    sample.info <- object@sampleInfo%>%
      dplyr::mutate(raw.exist = file.exists(raw.files),
                    ms.exist = file.exists(msData.files))%>%
      dplyr::filter(raw.exist,!ms.exist)

  }

  ### convert

  if (nrow(sample.info)) {
    MSconvertR::msConvert2mzML(raw.files  = object@sampleInfo$raw.files,
                               mzML.files = object@sampleInfo$msData.files,
                               BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))


  }

  object <- get_MSdev_MSinfo(object)

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



MSdev_extract_Spectra <- function(object){

  sampleInfo <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing %in% c("Both","MS2"))

  if (nrow(sampleInfo)==0) {
    sp <- Spectra::Spectra()
  } else {
    sp <- Spectra::Spectra(na.omit(sampleInfo$msData.files),backend = Spectra::MsBackendDataFrame())%>%
      filterMsLevel(2)
    sp$sp_id <- paste0("MS2_SP",num2str(1:length(sp)))
    Spectra::spectraNames(sp) <- sp$sp_id
  }

  object@spectra$MS2_Spectra <- sp
  return(object)
}




MSdev_match_Spectra_to_feature <- function(object){


  object@spectra$MS2_Spectra$ms2_matched_feature <-NA
  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    sp.ms2 <- object@spectra$MS2_Spectra%>%
      ProtGenerics::filterPolarity(i)
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
      as.data.frame()
    sp.ms2.data <- Spectra::spectraData(sp.ms2)%>%
      as.data.frame()%>%
      dplyr::mutate(precursorMZ = precursorMz,
                    retentionTime = rtime)%>%
      get_xcms_ms2_feature_id(featuredef = xcms.fdf)

    ### update MS2_Spectra
    sp.ms2.total <-object@spectra$MS2_Spectra %>%
      Spectra::spectraData()%>%
      as.data.frame()%>%
      dplyr::mutate(ms2_matched_feature = case_when(
        polarity==i ~ sp.ms2.data[sp_id,]$ms2_matched_feature,
        T~ms2_matched_feature
      ))
    Spectra::spectraData(object@spectra$MS2_Spectra ) <- DataFrame(sp.ms2.total)

    ### update xcms featuredef
    xcms.fdf$ms2_id <- sapply(xcms.fdf$feature_id,
                              function(i){
                                sp_id <- sp.ms2.data$sp_id[which(sp.ms2.data$ms2_matched_feature==i)]
                                return(sp_id)
                              }
    )
    xcms::featureDefinitions(xcms.xcms) <- xcms.fdf%>%
      DataFrame()
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]


  }

  return(object)


}

MSdev_annotation <- function(object,db.path = "d:/MSdb/msdb.HMDB.Rdata"){


  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.xcms <- get_xcms_feature_ms1_candidate(xcms.xcms,
                                   "d:/MSdb/msdb.HMDB.Rdata")
    xcms.xcms <- get_xcms_feature_ms2_score(xcms.xcms ,
                                            db.path = "d:/MSdb/msdb.HMDB.Rdata",
                                            object@spectra$MS2_Spectra)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]

  }
  object@projectInfo$MSdbPath <- db.path
  return(object)


}
