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
#' @param xcms.xcms
#' @param plot.title
#'
#' @return
#' @export
#'
#' @examples
plot_xcms_peaks_distribution <- function(xcms.xcms,plot.title = "Peaks distribution"){

  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::filter(!is.na(intb))
  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  ggplot(xcms.peaks)+
    geom_point(aes(x = rt,y=mz,
                   col = log10(intb),
                   alpha = log10(intb)/10,
                   size = (rtmax-rtmin)),
    )+
    scale_size_area(max_size = 8)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ),
         col = "Log10(Intensity)",
         size = "Peak width",
         x = "Retention time",
         y = "mz")+
    guides(alpha = "none")+
    scale_color_gradientn(breaks = c(0,3,6,10),
                          colors = c("white","orange","red","red"))+
    theme(text = element_text(size = 8))->peaks.dis.plot
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
    xcms.features$intb <-
      apply(featureValues(xcms.xcms, value = "intb"), 1, mean, na.rm = T)

    xcms.process.type <-
      processHistory(xcms.xcms) %>% sapply(processType)
    xcms.findpeak.param <-
      processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]] %>%
      processParam()
    ggplot(xcms.features) +
      geom_point(aes(
        x = rt,
        y = mz,
        col = log10(intb),
        alpha = log10(intb) / 10,
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
      scale_color_gradientn(breaks = c(0, 3, 6, 10),
                            colors = c("white", "orange", "red", "red")) +
      theme(text = element_text(size = 8)) -> peaks.dis.plot
    return(peaks.dis.plot)





  }
