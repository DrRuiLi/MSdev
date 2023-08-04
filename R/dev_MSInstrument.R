export_QE_InclusionList_From_xcmsFeature <- function(xcms.xcms){

  feature.def <- featureDefinitions_PeakSta(xcms.xcms)
  ion_mode <- unique(fData(xcms.xcms)$polarity)
  if (length(ion_mode)>1) {
    stop("multiple polarity")

  }else{

    ion_mode <- ifelse(ion_mode==1,"Positive","Negative")
  }
  inclusion.list <- QEinclusionListTemplate[rep(1,nrow(feature.def)),]
  inclusion.list$`Mass [m/z]` <- feature.def$mzmed
  inclusion.list$Polarity <- ion_mode
  inclusion.list$`Start [min]` <- feature.def$peakRtMin/60
  inclusion.list$`End [min]` <- feature.def$peakRtMax/60
  inclusion.list$Comment <- "Export from xcms"
  inclusion.list$maxo <-feature.def$peakMaxo%>%log10
  inclusion.list
}

export_QE_ExclusionList_From_xcmsPeaks <- function(xcms.xcms,peak.count.thresh = 10){


  raw.file <- fileNames(xcms.xcms)
  if (length(raw.file)>1) {
    raw.file <- paste0(dirname(raw.file[1]),"/multiple_file")
  }
  ion_mode <- unique(fData(xcms.xcms)$polarity)
  if (length(ion_mode)>1) {
    stop("multiple polarity")

  }else{

    ion_mode <- ifelse(ion_mode==1,"Positive","Negative")
  }
  ms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    arrange(mz)%>%
    mutate(mzdiff = c(diff(mz),0),
           mz.group = groupMz(mz)$mz.group)

  ms.peaks.count <- ms.peaks %>%
    group_by(mz.group)%>%
    mutate(mz.group.count=length(mz.group) )%>%
    ungroup()%>%
    filter(mz.group.count > peak.count.thresh)
  ggplot(ms.peaks.count)+
    geom_segment(aes(x = rtmin,xend =  rtmax , y = mz,yend = mz,col = mz.group),size = 0.1)+
    labs(title = paste0("Total ",length(unique(ms.peaks.count$mz.group))," ions"),
         x = "rt")+
    guides(col = "none")+
    theme_bw()+
    theme(text = element_text(size = 8))->gp
  gp
  peaks.to.exclusion <-ms.peaks.count%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate(rt.min = min(rtmin),
                  rt.max = max(rtmax),
                  mz = median(mz))%>%
    dplyr::ungroup()%>%
    dplyr::select(mz.group,mz,rt.min,rt.max)%>%
    dplyr::distinct(mz.group,.keep_all = T)
  exclusion.list <- QEExclusionListTemplate[rep(1,nrow(peaks.to.exclusion)),]%>%
    dplyr::mutate(`Mass [m/z]` = peaks.to.exclusion$mz,
                  `Start [min]` = peaks.to.exclusion$rt.min/60,
                  `End [min]` = peaks.to.exclusion$rt.max/60,
                  Polarity = ion_mode,
                  Comment =  "Export from xcms")
  export::graph2png(gp,file =paste0(dirname(raw.file),"/ExclusionList_From_",basename(raw.file),".png") ,
                    width =5,height = 4)
  write_csv(exclusion.list,file =paste0(dirname(raw.file),"/ExclusionList_From_",basename(raw.file),".csv"))
  return(gp)

}
