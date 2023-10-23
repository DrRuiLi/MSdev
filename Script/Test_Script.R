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




msdev.demo <- MSdev()
msdev.demo <- checkSampleInfo(msdev.demo)
msdev.demo <- msConvert_MSdev(msdev.demo)
msdev.demo <- xcmsProcessingMSdev(msdev.demo)
msdev.demo <- extract_Spectra_MSdev(msdev.demo)
msdev.demo <- match_Spectra_to_feature_MSdev(msdev.demo)


xcms.xcms <- load_demo("XCMS")
get_xcms_MS_report()



a <- readxl::read_excel("aaa.xlsx")

plot.data <- a%>%
  dplyr::mutate(date = as.Date(date))%>%
  dplyr::group_by(name)%>%
  dplyr::mutate(group = case_when(date <max(date)~"HUA",
                                  T~"Gout"))%>%
  dplyr::select(!id)%>%
  pivot_wider(names_from = "group",
              values_from = "date")

ggplot(plot.data)+
  geom_segment(aes(x = HUA,xend = Gout,
                   y= name,yend= name))+
  geom_point(aes(x = HUA,
                   y= name,
                 col = "#91D1C2"),show.legend = F)+
  geom_point(aes(x = Gout,
                 y= name,
                 col = "#E64B35"),show.legend = F)+
  labs( x = "Date", y = "", title = "Sample Collection before and after Gout")+
  theme_bw()+
  theme(axis.text.y = element_blank())->p

open_ggplot_win(p,6,8)


QE.POS <- load_as_var("d:/2023.09.02.DDA.mine/Pos/MSdev_2023_09_02.back.Rdata")
QE.POS.list <-QE.POS@statData[["DDA.mine.list.Positive"]][["DDA.mine.list001"]]











