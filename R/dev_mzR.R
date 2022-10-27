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
