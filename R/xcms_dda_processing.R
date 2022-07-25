xcms_dda_processing <- function(ms.ana ,  polarity ){


  message(Sys.time()," xcms processing ",polarity)

  ### data to xcms
  {
    sample.info.xcms <- ms.ana$sample.info %>%
      dplyr::mutate(mzML.file = eval(parse(text = paste0("mzML.file.",polarity))),
                    raw.file = eval(parse(text = paste0("raw.file.",polarity))),
                    analysis.time = eval(parse(text = paste0("analysis.time.",polarity))),
                    polarity = polarity)%>%
      dplyr::select("sample.name" ,           "sample.type"        , "polarity",   "sample.abbreviation" ,
                    "mzML.file"   ,           "raw.file"           ,    "analysis.time")%>%
      dplyr::filter(!is.na(mzML.file))%>%
      dplyr::mutate(injection.order = rank(analysis.time),.after = sample.name)%>%
      dplyr::arrange(injection.order)



  }

  xcms.processing <- ms.ana[["processing.info"]][["xcms.processing"]][[polarity]]

  ###read ms data
  if (!isTRUE(xcms.processing$readmsdata$done ) ){
    xcms.xcms <-  readMSData(sample.info.xcms$mzML.file, mode = "onDisk")
    xcms.processing$readmsdata$done <- T
    xcms.processing$readmsdata$time <- Sys.time()
  }else{
    xcms.xcms <- ms.ana[[ paste0("xcms.",polarity)]]

  }


  ###find peaks
  if (!isTRUE(xcms.processing$findpeak$done ) ){

    centwave.param <- CentWaveParam(peakwidth = c(5,30),
                                    prefilter = c(3,100),
                                    ppm = 20)
    xcms.xcms<-findChromPeaks(xcms.xcms,
                              param = centwave.param)

    xcms.processing$findpeak$done <- T
    xcms.processing$findpeak$param <- centwave.param
    xcms.processing$findpeak$time <- Sys.time()

  }else{
    xcms.xcms <- ms.ana[[ paste0("xcms.",polarity)]]

  }

  ###adjust retention time

  if (!isTRUE(xcms.processing$adjustRT$done ) ){

    peak.density.param <- PeakDensityParam(sampleGroups = sample.info.xcms$sample.type,
                                           minFraction = 0.8,bw = 30,
                                           binSize = 0.015)

    xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)

    peak.group.param <- PeakGroupsParam(minFraction = 0.85,
                                        subset = which(sample.info.xcms$sample.type == "QC"),
                                        subsetAdjust = "average",span = 0.4)

    xcms.xcms <- adjustRtime(xcms.xcms,param = peak.group.param)


    xcms.processing$adjustRT$done <- T
    xcms.processing$adjustRT$param <- peak.group.param
    xcms.processing$adjustRT$time <- Sys.time()

  }else{
    xcms.xcms <- ms.ana[[ paste0("xcms.",polarity)]]

  }



  ### group features
  if (!isTRUE(xcms.processing$correspondence$done ) ){

    peak.density.param <- PeakDensityParam(sampleGroups = sample.info.xcms$sample.type,
                                           minFraction = 0.5,bw = 30,
                                           binSize = 0.015)

    xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)
    xcms.xcms <- fillChromPeaks(xcms.xcms,param = FillChromPeaksParam())


    xcms.processing$correspondence$done <- T
    xcms.processing$correspondence$param <- peak.density.param
    xcms.processing$correspondence$time <- Sys.time()

  }else{
    xcms.xcms <- ms.ana[[ paste0("xcms.",polarity)]]

  }


  ms.ana[["processing.info"]][["xcms.processing"]][[polarity]] <- xcms.processing
  ms.ana[[ paste0("xcms.",polarity)]]  <- xcms.xcms
  save(ms.ana,file = ms.ana$processing.info$project.info$ms.ana.file)

  return(ms.ana)


}

