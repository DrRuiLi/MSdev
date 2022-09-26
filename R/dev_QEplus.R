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
  inclusion.list
}


