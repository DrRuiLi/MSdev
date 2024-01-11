pseudo_workflow <- function(){

  #msdev.pseudo <- MSdev("../../Projecct/2023.07.12.Pseudo.dev/msdata/")
  msdev.pseudo <- load_as_var("../../Projecct/2023.07.12.Pseudo.dev/MSdev_2023_07_18.Rdata")
  msdev.pseudo <- msConvert_MSdev(msdev.pseudo)
  msdev.pseudo <- checkSampleInfo(msdev.pseudo)
  msdev.pseudo <- xcmsProcessingMSdev(msdev.pseudo)

  msdev.pseudo <- extractSpectra_fullscan_DDA(msdev.pseudo)
  msdev.pseudo <- featureSpectra_fullscan_DDA(msdev.pseudo)
  msdev.pseudo <- featureCandidate(msdev.pseudo,
                                   mz.ppm = 20,
                                   spectraDatabase = "d:/MSdb.2023.05.30/MSdb_Spectra_database_2023_06_06.rda")

  msdev.pseudo <- annotateMSdev(msdev.pseudo)
  msdev.pseudo <- getStaDataMSdev(msdev.pseudo)









}




