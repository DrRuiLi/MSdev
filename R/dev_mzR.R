#' @title getmsExpTime
#' @description get MS experiment time by mzR::runInfo()$startTimeStamp
#' @param msdataFiles file path
#'
#' @return a data.frame with two column: msDataFiles, ExpTime
#' @export
#'

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
#' @title Get Msinfo Mzr
#' @description MSinfo mzR.
#' @param msDataFiles file path
#'
#' @return mzR info
#' @export
#'

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
    tr <- try(msmzR <- mzR::openMSfile(msDataFiles[i]),silent  = T)
    if (is.character(tr)&&grepl("Error",tr)) next
    msrunInfo <- mzR::runInfo(msmzR)
    msinstru <- mzR::instrumentInfo(msmzR)
    msheader <- mzR::header(msmzR)
    msdata$ExpTime[i] <- msrunInfo$startTimeStamp
    msdata$msLevels[i] <- dplyr::case_when(
      rlang::is_empty(msrunInfo$msLevels) ~ NA_character_,
      TRUE ~ {
        lv <- sort(unique(as.integer(msrunInfo$msLevels)))
        lv <- lv[!is.na(lv)]
        if (!length(lv)) {
          NA_character_
        } else if (all(c(1L, 2L) %in% lv)) {
          "1;2"
        } else {
          as.character(lv[[length(lv)]])
        }
      }
    )
    msdata$polarity[i] <- case_when(rlang::is_empty(unique(msheader$polarity))~NA,
                                    T~paste0(unique(msheader$polarity),collapse = ";"))
    msdata$manufacturer[i] <- msinstru$manufacturer
    msdata$model[i] <- msinstru$model

    mzR::close(msmzR)

  }
  return(msdata)

}
