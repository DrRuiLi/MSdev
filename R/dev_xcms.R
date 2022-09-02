#' @title  get_features_from_xcms
#' @description extract feature data from xcms::XCMSnExp,
#'  calculate RSD of QC and Sample
#'  ( note this rely on character "QC" and "Sample" in `sampleNames(xcms.xcms)` )
#' @param xcms.xcms
#'
#' @return
#' @export
#'
#' @examples
get_features_from_xcms <- function(xcms.xcms){

  xcms.sum <- quantify(xcms.xcms,missing = "rowmin_half")
  feature.def <- SummarizedExperiment::rowData(xcms.sum)%>%
    tibble::as_tibble()

  feature.matrix <- SummarizedExperiment::assay(xcms.sum)
  rsd <- function(x){sd(x,na.rm =  T)/mean(x , na.rm = T)}
  feature.matrix.qc <- feature.matrix[,which(grepl("QC",colnames(feature.matrix)))]
  feature.matrix.sample <- feature.matrix[,which(grepl("Sample",colnames(feature.matrix)))]
  feature.def$qc_rsd <- apply(feature.matrix.qc, 1, rsd)
  feature.def$sample_rsd <- apply(feature.matrix.sample, 1, rsd)
  feature.def$med_intensity <- apply(feature.matrix , 1 ,median)
  SummarizedExperiment::rowData(xcms.sum) <-feature.def
  return(xcms.sum)
}










#' @title plot_xcms_peaks_distribution
#' @description export peaks data by xcms::chromPeaks and plot by ggplot2
#'
#' @param xcms.xcms
#' @param plot.title
#' @param type `"o"`, for geom_point, `"l"`, for geom_segment
#'
#' @return
#' @export
#'
#' @examples
plot_xcms_peaks_distribution <- function(xcms.xcms,plot.title = "Peaks distribution",type = "o"){

  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::filter(!is.na(maxo))
  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  if (type == "o") {
    ggplot(xcms.peaks)+
      geom_point(aes(x = rt,y=mz,
                     col = log10(maxo),
                     alpha = log10(maxo)/10,
                     size = (rtmax-rtmin)),
      )+
      scale_size_area(max_size = 8)+
      labs(title = plot.title,
           subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                             "; SN = ",xcms.findpeak.param@snthresh,
                             "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
           col = "Log10\n(Intensity)",
           size = "Peak width",
           x = "Retention time",
           y = "mz")+
      guides(alpha = "none")+
      scale_color_gradientn(breaks = c(0,3,6,9),
                            labels = c(0,3,6,9),
                            limits = c(0,9),
                            values = c(0,2,4,7,9)/9,
                            colors = c("white","white","yellow","red","red"))+
      theme_bw()+
      theme(text = element_text(size = 8))->peaks.dis.plot

  }else if(type == "l"){
    ggplot(xcms.peaks)+
      geom_segment(aes(x = rtmin , xend = rtmax , y = mz, yend = mz,col = log10(maxo)),
                   size = 0.6
      )+
      labs(title = plot.title,
           subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                             "; SN = ",xcms.findpeak.param@snthresh,
                             "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
           col = "Log10(Intensity)",
           x = "Retention time",
           y = "mz")+
      guides(alpha = "none")+
      scale_color_gradientn(breaks = c(0,3,6,9),
                            labels = c(0,3,6,9),
                            limits = c(0,9),
                            values = c(0,2,4,7,9)/9,
                            colors = c("white","white","yellow","red","red"))+
      theme_bw()+
      theme(text = element_text(size = 8))->peaks.dis.plot
    peaks.dis.plot


  }
  return(peaks.dis.plot)


}



#' @title plot_xcms_peaks_distribution
#' @description plot_xcms_peaks_distribution
#' @param xcms.xcms
#' @param plot.title
#'
#' @return
#' @export
#'
#' @examples
plot_xcms_features_distribution <-
  function(xcms.xcms, plot.title = "Features distribution") {
    xcms.features <- featureDefinitions(xcms.xcms) %>%
      as.data.frame() %>%
      mutate(mz = mzmed, rt = rtmed)
    xcms.features$maxo <-
      apply(featureValues(xcms.xcms, value = "maxo"), 1, mean, na.rm = T)

    xcms.process.type <-
      processHistory(xcms.xcms) %>% sapply(processType)
    xcms.findpeak.param <-
      processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]] %>%
      processParam()
    ggplot(xcms.features) +
      geom_point(aes(
        x = rt,
        y = mz,
        col = log10(maxo),
        alpha = log10(maxo) / 10,
        size = (rtmax - rtmin)
      ),) +
      scale_size(range = c(2,6)) +
      xlim(c(0, 800)) +
      labs(
        title = plot.title,
        subtitle = paste0(
          "ppm = ",
          xcms.findpeak.param@ppm,
          "; SN = ",
          xcms.findpeak.param@snthresh,
          "; prefilter = (",
          paste0(xcms.findpeak.param@prefilter, collapse = ","),
          ")"
        ),
        col = "Log10(Intensity)",
        size = "Peak width",
        x = "Retention time",
        y = "mz"
      ) +
      guides(alpha = "none") +
      scale_color_gradientn(breaks = c(0,3,6,9),
                            labels = c(0,3,6,9),
                            limits = c(0,9),
                            values = c(0,2,4,7,9)/9,
                            colors = c("white","white","yellow","red","red"))+
      theme(text = element_text(size = 8)) -> peaks.dis.plot
    return(peaks.dis.plot)





  }







plot_xcms_peaks_mzerror_density <- function(xcms.xcms,
                                            plot.title = "Peak mz error distribution"){

  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    mutate(ppm = (mzmax-mzmin)/mz*1e6,
           mz_diff = mzmax-mzmin,
           peak_width = rtmax-rtmin)

  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  ggplot(xcms.peaks,aes(x = mz , y = ppm)) +
    stat_density_2d(aes(fill= ..level..),
                    contour = T,
                    geom = "polygon",bins = 100)+
    geom_point(size = 0.1,alpha = 0.1)+
    scale_fill_gradient(low="#00000001",high = "red")+
    scale_x_continuous(expand = c(0.1,0.1))+
    scale_y_continuous(expand = c(0.1,0.1))+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" )
         )+
    guides(fill = "none")+
    theme_bw()+
    theme(text = element_text(size = 8))

}



#' @title plot_xcms_peaks_ms1_scans
#' @description plot scans number of MS1 levels in each peak, note that to many peaks will lead to stuck,
#' apply `filterFile` to decrease peaks count
#' @param xcms.xcms should be a `XCMSnExp` object after `findChromPeaks`
#' @param plot.title
#'
#' @return
#' @export
#'
#' @examples
plot_xcms_peaks_ms1_scans <- function(xcms.xcms,plot.title = "Peaks Sans of MS1"){

  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()
  xcms.scans <- fData(xcms.xcms)%>%
    dplyr::filter(msLevel== 1)
  peaks_scans <- function(x,xcms.scans){
    sum(x["rtmax"] > xcms.scans$retentionTime  & x["rtmin"] < xcms.scans$retentionTime )

  }
  xcms.peaks$ms1_scans_no <- apply(xcms.peaks ,1,peaks_scans , xcms.scans)
  ggplot(xcms.peaks)+
    geom_segment(aes(x = rtmin , xend = rtmax , y = ms1_scans_no, yend = ms1_scans_no,col = log10(maxo)),
                 size = 0.6
    )+
    geom_hline(yintercept = 7)+
    geom_boxplot(aes( x = max(rt)*1.2 , y =ms1_scans_no),width = diff(range(xcms.peaks$rt))*0.1)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
         col = "Log10(Intensity)",
         x = "Retention time",
         y = "Scan count of MS1 in each peak")+
    scale_y_log10(breaks = c(1,2,3,4,5,6,7,8,10,20))+
    guides(alpha = "none")+
    scale_color_gradientn(breaks = c(0,3,6,9),
                          labels = c(0,3,6,9),
                          limits = c(0,log10(median(xcms.peaks$maxo))*1.5),
                          values = c(0,2,4,7,9)/9,
                          na.value = "red",
                          colors = c("white","white","yellow","red","red"))+
    theme_bw()+
    theme(text = element_text(size = 8))->peaks.dis.plot
  peaks.dis.plot

}



plot_xcms_peaks_ms2_scans <- function(xcms.xcms,plot.title = "Peaks Sans of MS2"){

  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()
  xcms.scans <- fData(xcms.xcms)%>%
    dplyr::filter(msLevel== 2)
  peaks_scans <- function(x,xcms.scans){
    sum(x["rtmax"] > xcms.scans$retentionTime  & x["rtmin"] < xcms.scans$retentionTime&
          x["mzmax"] > xcms.scans$precursorMZ  & x["mzmin"] < xcms.scans$precursorMZ)

  }
  xcms.peaks$ms2_scans_no <- apply(xcms.peaks ,1,peaks_scans , xcms.scans)
  ggplot(xcms.peaks)+
    geom_jitter(aes(x = rt, y = ms2_scans_no, col = log10(maxo)),
                 size = 0.6
    )+
    #geom_hline(yintercept = 7)+
    geom_violin(aes( x = max(rt)*1.2 , y =ms2_scans_no),width = diff(range(xcms.peaks$rt))*0.1)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
         col = "Log10(Intensity)",
         x = "Retention time",
         y = "Scan count of MS2 in each peak")+
    guides(alpha = "none")+
    scale_color_gradientn(breaks = c(0,3,6,9),
                          labels = c(0,3,6,9),
                          limits = c(0,log10(median(xcms.peaks$maxo))*1.5),
                          values = c(0,2,4,7,9)/9,
                          na.value = "red",
                          colors = c("white","white","yellow","red","red"))+
    theme_bw()+
    theme(text = element_text(size = 8))->peaks.dis.plot
  peaks.dis.plot

}



plot_xcms_peaks_SN_distribution <- function(xcms.xcms,plot.title = "Peaks SNR(Signal to Noise Ratio)"){


  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()

  ggplot(xcms.peaks)+
    geom_jitter(aes(x = rt, y = log10(sn), col = log10(maxo)),
                size = 0.6
    )+
    #geom_hline(yintercept = 7)+
    geom_violin(aes( x = max(rt)*1.2 , y =log10(sn)),width = diff(range(xcms.peaks$rt))*0.1)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
         col = "Log10(Intensity)",
         x = "Retention time",
         y = "log10(SNR)")+
    guides(alpha = "none")+
    scale_color_gradientn(breaks = c(0,3,6,9),
                          labels = c(0,3,6,9),
                          limits = c(0,log10(median(xcms.peaks$maxo))*1.5),
                          values = c(0,2,4,7,9)/9,
                          na.value = "red",
                          colors = c("white","white","yellow","red","red"))+
    theme_bw()+
    theme(text = element_text(size = 8))->peaks.dis.plot
  peaks.dis.plot



}

plot_xcms_peaks_Chromatogram <- function(xcms.xcms,peak_id){

  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()
  xcms.peaks <- xcms.peaks[peak_id,]
  mz.range <-c(xcms.peaks$mzmin,xcms.peaks$mzmax)
  rt.range <-c(xcms.peaks$rtmin,xcms.peaks$rtmax)

  xcms.chrom <- chromatogram(xcms.xcms ,
                             mz = mz.range,
                             rt = rt.range+diff(rt.range)*c(-1.5,1.5)
                             )
  xcms.chrom <- xcms.chrom[1,1]
  chrom.data <- data.frame(rt = rtime(xcms.chrom),
                           intensity = intensity(xcms.chrom))%>%
    #dplyr::filter(!is.na(intensity))%>%
    dplyr::mutate(fill = rt > min(rt.range)&rt <max(rt.range))

  ggplot(chrom.data)+
    geom_line(aes(x = rt,y = intensity),linetype = 1)+
    geom_area(aes(x = rt,y = intensity, fill = fill))+
    scale_fill_manual(values = c("FALSE" = "transparent","TRUE" = "grey"))+
    labs(title = paste0(xcms.chrom@chromPeakData@rownames),
         subtitle = paste0("mz:",paste0(sprintf("%.5f",xcms.chrom@mz),collapse = " - "), ";     ",
                           "rt:",paste0(sprintf("%.2f",rt.range),collapse = " - "),"\n",
                           "mz error = ",sprintf("%.2f",mean(diff(xcms.chrom@mz)/xcms.chrom@mz)*1e6)," ppm;     ",
                           "peak width = ", sprintf("%.2f",diff(rt.range))
                           ),
         x = "Retention time")+
    guides(fill = "none")+
    theme_bw()+
    theme(text = element_text(size = 8))



}

























