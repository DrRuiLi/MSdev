QE_DDA_Mine_workflow <- function(){

  ### run after QC FS, once
  data.dir <- "d:/DDA.mine.test/pos"
  msdev.qe <- MSdev(rawDataDir = data.dir)
  msdev.qe <- MSdev_load("d:/DDA.mine.test/MSdev_2025_06_18.Rdata")
  msdev.qe <- MSdev_msConvert(msdev.qe)
  msdev.qe <- MSdev_xcmsProcessing(msdev.qe)

  msdev.qe@statData <- list()
  msdev.qe <- MSdev_get_Inclusion_Queue(msdev.qe)

  ### run after every time DDA acquired
  msdev.qe <- MSdev_get_Inclusion_List(msdev.qe)
  msdev.qe <- MSdev_add_sample(msdev.qe,raw.data.dir = "d:/DDA.mine.test/pos")
  msdev.qe <- MSdev_get_MS2acquisitionStat(msdev.qe)

  table(msdev.qe@statData$DDA_mine_queue_Positive$acquired)





}
