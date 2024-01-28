#' @title getmsExpTime
#' @description get MS experiment time by mzR::runInfo()$startTimeStamp
#' @param msdataFiles
#'
#' @return a data.frame with two column: msDataFiles, ExpTime
#' @export
#'
#' @examples
getmsExpTime<- function(msDataFiles){
  msdata <- data.frame(files = msDataFiles,
                       ExpTime = NA)
  for (i in 1:length(msDataFiles)) {
    #cat(i,"\n")
    if (is.na(msDataFiles[i])) {
      msdata$ExpTime[i] <- NA
      next

    }
    msmzR <- mzR::openMSfile(msDataFiles[i])
    msdata$ExpTime[i] <- mzR::runInfo(msmzR)$startTimeStamp
    mzR::close(msmzR)
  }
  return(msdata)

}







#' get_MSinfo_mzR
#' get MS experiment time by mzR
#'
#' @param msDataFiles
#'
#' @return
#' @export
#'
#' @examples
get_MSinfo_mzR<- function(msDataFiles){

  msdata <- data.frame(files = msDataFiles,
                       ExpTime = NA,
                       msLevels = NA,
                       polarity = NA,
                       manufacturer=NA,
                       model =NA
                       )
  for (i in 1:length(msDataFiles)) {
    #cat(i,"\n")
    if (is.na(msDataFiles[i])|!file.exists(msDataFiles[i])) {
      next

    }
    msmzR <- mzR::openMSfile(msDataFiles[i])
    msrunInfo <- mzR::runInfo(msmzR)
    msinstru <- mzR::instrumentInfo(msmzR)
    msheader <- mzR::header(msmzR)
    msdata$ExpTime[i] <- msrunInfo$startTimeStamp
    msdata$msLevels[i] <- case_when(is_empty(msrunInfo$msLevels)~NA,
                                 T~max(msrunInfo$msLevels))
    msdata$polarity[i] <- case_when(is_empty(unique(msheader$polarity))~NA,
                                    T~paste0(unique(msheader$polarity),collapse = ";"))
    msdata$manufacturer[i] <- msinstru$manufacturer
    msdata$model[i] <- msinstru$model

    mzR::close(msmzR)
  }
  return(msdata)

}
