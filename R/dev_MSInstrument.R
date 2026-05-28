export_QE_InclusionList_From_xcmsFeature <- function(xcms.xcms){

  feature.def <- featureDefinitions_PeakSta(xcms.xcms)
  ion_mode <- unique(Biobase::fData(xcms.xcms)$polarity)
  if (length(ion_mode)>1) {
    stop("multiple polarity")

  }else{

    ion_mode <- ifelse(ion_mode==1,"Positive","Negative")
  }
  inclusion.list <- QEinclusionListTemplate[rep(1,nrow(feature.def)),]
  inclusion.list$`Mass [m/z]` <- feature.def$mzmed
  inclusion.list$Polarity <- ion_mode
  inclusion.list$`Start [min]` <- (feature.def$peakRtMin-10)%>%
    ifelse(. < 0 ,0,.)/60
  inclusion.list$`End [min]` <- (feature.def$peakRtMax+10) %>%
    ifelse(. >max(rtime(xcms.xcms)),max(rtime(xcms.xcms)),.)/60
  inclusion.list$Comment <- "Export from MSdev"
  return(inclusion.list)
}

export_QE_ExclusionList_From_xcmsPeaks <- function(xcms.xcms,peak.count.thresh = 10){


  raw.file <- MSnbase::fileNames(xcms.xcms)
  if (length(raw.file)>1) {
    raw.file <- paste0(dirname(raw.file[1]),"/multiple_file")
  }
  ion_mode <- unique(fData(xcms.xcms)$polarity)
  if (length(ion_mode)>1) {
    stop("multiple polarity")

  }else{

    ion_mode <- ifelse(ion_mode==1,"Positive","Negative")
  }
  ms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
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


get_MRM_list <- function(feature_def_sta){


  median.peak.width <- 5
  time.max <-max(xcms.features.transition$peakRtMax)

  feature_def_sta.calc <- feature_def_sta%>%
    dplyr::mutate(peakRtMin = peakRtMin-10,
                  peakRtMax = peakRtMax+10,
                  peakRtMin = case_when(peakRtMin <0~0,
                                        T~peakRtMin),
                  peakRtMax = case_when(peakRtMax >time.max~time.max,
                                        T~peakRtMax))
  #median.peak.width <- median(srm.list$peakRtMax- srm.list$peakRtMin)/2

  time.seq <- seq(0,time.max+median.peak.width,median.peak.width)
  scan.time <- 1
  dwt.list <- list()
  for (i in 1:(length(time.seq)-1)) {
    time.window <- time.seq[c(i,i+1)]

    dwt.list[[i]] <- feature_def_sta.calc%>%
      dplyr::mutate(in.window = peakRtMin < time.window[2] & peakRtMax > time.window[1] )%>%
      dplyr::group_by(in.window)%>%
      dplyr::mutate(ion_count = n( ),
                    part =1/log(peakMaxo),
                    dwt = scan.time *part/sum(part) *1000,
                    dwt = case_when(in.window~dwt,T~NA))


  }

  dwt.matrix <- sapply(dwt.list,function(x){
    x$dwt
  },USE.NAMES = T)
  feature_def_sta.calc$dwt.mean <- apply(dwt.matrix , 1 , mean,na.rm =T)

  time.window.count <- sapply(dwt.list,function(x){
    min(x$ion_count)
  })


  srm.list <- feature_def_sta.calc%>%
    dplyr::mutate(dwt = case_when(dwt.mean > 50~50,
                                  T~dwt.mean),
                  peakRtMin = peakRtMin/60,
                  peakRtMax= peakRtMax/60
                  )%>%
    dplyr::select(Compound = feature.id,
                  `Start Time (min)` =peakRtMin,
                  `End Time (min)`=peakRtMax,
                  `Precursor (m/z)`=precursorMz,
                  `Product (m/z)` = productMz,
                  `Collision Energy (V)` = collisionEnergy,
                  `Dwell Time (ms)` = dwt,
                  everything(),
                  -peakidx
                  )


  return(srm.list)



}




QE_list_2feature_def <- function(table_to_trans,keep = T ){

  # featureDefinitions_PeakSta(MSdev@xcmsData$PositiveMS1)->table_to_trans
  # table_to_trans <- dda.acq.list
  var.map <-c(mzmed = "Mass [m/z]",
              peakRtMin = "Start [min]" ,
              peakRtMax = "End [min]",
              collisionEnergy = "(N)CE",
              feature_id = "Comment",
              polarity = "Polarity")
  if (any(colnames(table_to_trans) %in% var.map) ) {
    type <- "QE_list"
    rt_multi <- 60

    var.map.matched <- names(var.map)[var.map%in% colnames(table_to_trans)]
    names(var.map.matched) <- var.map[var.map%in% colnames(table_to_trans)]
    table_transed <- table_to_trans
    table_transed <- dplyr_rename_to(table_transed,names(var.map.matched),var.map.matched)



  }else if (any(colnames(table_to_trans) %in% names(var.map))) {

    type <- "feature_df"
    rt_multi <- 1/60
    var.map.matched <- var.map[names(var.map)%in% colnames(table_to_trans)]
    table_transed <- table_to_trans
  }else{
    stop("var not match")
  }

  if ("rtmin" %in% colnames(table_transed))
    table_transed$peakRtMin <- table_transed$peakRtMin *rt_multi
  if ("rtmax" %in% colnames(table_transed))
    table_transed$peakRtMax <- table_transed$peakRtMax *rt_multi
  if ("polarity" %in% colnames(table_transed)){
    if (type =="feature_df") {
      table_transed$polarity <- ifelse(table_transed$polarity == 0,"Negative","Positive")
    }
    if (type =="QE_list") {
      table_transed$polarity <- ifelse(table_transed$polarity == "Positive",1,0)
    }
  }


  table_transed <- dplyr_rename_to(table_transed,old = names(var.map.matched),new = var.map.matched)
  if (!keep)
    table_transed <- table_transed[,var.map.matched]
  return(table_transed)



}


