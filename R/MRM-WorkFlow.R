MRM_WorkFlow <- function(rawDataDir){

  rawDataDir <- "d:/2023.08.22.Pseudo/Result/"

  msdev.mrm <- MSdev(rawDataDir = rawDataDir )
  ### project info
  {
    msdev.mrm@projectInfo[["msAcquisition"]] <- "SRM"
    msdev.mrm@experimentInfo <-MSdev::MS_Experiment[9]
  }
  msdev.mrm <- checkSampleInfo(msdev.mrm)
  msdev.mrm <- msConvert_MSdev(msdev.mrm)
  msdev.mrm <- xcmsProcessingMSdev(msdev.mrm)





}
