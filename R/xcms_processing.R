#' @title xcms_procecing
#' @description read in files in sample.info$mzML.files. note the sample were reordered in xcms object by analysis time, may not identical to sample.name, should be correspondended by mzML.file
#' @param sample.info
#' @import xcms
#' @import MsFeatures
#' @return
#' @export
#'
#' @examples
xcms_processing <- function(sample.info, ion_mode ){

  sample.info.xcms <- sample.info %>%
    dplyr::mutate(mzML.file = eval(parse(text = paste0("mzML.file.",ion_mode))),
                  raw.file = eval(parse(text = paste0("raw.file.",ion_mode))),
                  analysis.time = eval(parse(text = paste0("analysis.time.",ion_mode))),
                  ion_mode = ion_mode)%>%
    dplyr::select("sample.name" ,           "sample.type"        , "ion_mode",   "sample.abbreviation" ,
                  "mzML.file"   ,           "raw.file"           ,    "analysis.time")%>%
    dplyr::filter(!is.na(mzML.file))%>%
    dplyr::mutate(injection.order = rank(analysis.time),.after = sample.name)%>%
    dplyr::arrange(injection.order)


  message(Sys.time()," xcms data importing... ",ion_mode)
  xcms.xcms <-  readMSData(sample.info.xcms$mzML.file, mode = "onDisk")
  pData(xcms.xcms) <- cbind(pData(xcms.xcms),sample.info.xcms)

  message(Sys.time()," xcms peak finding...",ion_mode)
  centwave.param <- CentWaveParam(peakwidth = c(5,50),
                                  prefilter = c(3,100),
                                  snthresh = 10,
                                  ppm = 20)
  xcms.xcms<-findChromPeaks(xcms.xcms,
                            param = centwave.param)

  message(Sys.time()," xcms rime adjusting... ",ion_mode)
  peak.density.param <- PeakDensityParam(sampleGroups = sample.info.xcms$sample.type,
                                         minFraction = 0.8,bw = 30,
                                         binSize = 0.015)
  xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)
  peak.group.param <- PeakGroupsParam(minFraction = 0.85,
                                      subset = which(sample.info.xcms$sample.type == "QC"),
                                      subsetAdjust = "average",span = 0.4)

  xcms.xcms <- adjustRtime(xcms.xcms,param = peak.group.param)

  message(Sys.time()," xcms peak grouping... ",ion_mode)
  peak.density.param <- PeakDensityParam(sampleGroups = sample.info.xcms$sample.type,
                                         minFraction = 0.5,bw = 30,
                                         binSize = 0.015)

  xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)
  xcms.xcms <- fillChromPeaks(xcms.xcms,param = FillChromPeaksParam())


  message(Sys.time()," xcms features grouping... ",ion_mode)

  return(xcms.xcms)



  }
