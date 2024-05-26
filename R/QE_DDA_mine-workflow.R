QE_DDA_Mine_workflow <- function(){


  ### run after QC FS, once
  data.dir <- "d:/DDA_mine_test/Pos/rawData/"
  msdev.qe <- MSdev(rawDataDir = data.dir)
  msdev.qe <- msConvert_MSdev(msdev.qe)
  msdev.qe <- xcmsProcessingMSdev(msdev.qe)
  msdev.qe <- get_MSdev_Inclusion_Queue(msdev.qe)
  ### run after every time DDA acquired
  msdev.qe <- get_MSdev_Inclusion_List(msdev.qe)
  msdev.qe <- get_MSdev_newSamples(msdev.qe)
  msdev.qe <- get_MSdev_MS2acquisitionStat(msdev.qe)

  table(msdev.qe@statData$DDA.mine.queue.Positive$acquired.time)






}
