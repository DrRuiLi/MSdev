# Wed Oct 18 21:18:36 2023 ------------------------------
### create demo

msdev.demo <- MSdev()
msdev.demo <- checkSampleInfo(msdev.demo)
msdev.demo <- msConvert_MSdev(msdev.demo)
msdev.demo <- xcmsProcessingMSdev(msdev.demo)

xcms.xcms <- msdev.demo@xcmsData$PositiveMS1
save(xcms.xcms,file = paste0(msdev.demo@projectInfo$projectDir,"/XCMSnExp_2023_10_19.rda"))

# Thu Oct 19 10:36:46 2023 ------------------------------
xcms.xcms <- filterFile(xcms.xcms , 2)

xcms.peaks.centwave.iso <- findChromPeaks(xcms.xcms,
                             param = CentWavePredIsoParam(),
                             BPPARAM = SerialParam())%>%
  groupChromPeaks(param = PeakDensityParam(1))

xcms.peaks.centwave <- findChromPeaks(xcms.xcms,
                                      param = CentWaveParam(),
                                      BPPARAM = SerialParam())%>%
  groupChromPeaks(param = PeakDensityParam(1))


a <- chromPeaks(xcms.peaks.centwave)%>%
  as.data.frame()%>%
  dplyr::mutate(cluster_ion(mz,rt))

b <- chromPeaks(xcms.peaks.centwave.iso)%>%
  as.data.frame()%>%
  dplyr::mutate(cluster_ion(mz,rt))


c <- featureDefinitions_PeakSta(xcms.peaks.centwave)



xcms.xcms %>%
  filterMz( mz.range.ppm(828.55040,25)) %>%
  filterRt(c(150,200))%>%
  plot(type = "XIC")




