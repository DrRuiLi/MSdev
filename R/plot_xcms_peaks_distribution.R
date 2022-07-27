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
    as.data.frame()
  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  ggplot(xcms.peaks)+
    geom_point(aes(x = rt,y=mz,
                   col = log10(intb),
                   alpha = log10(intb)/10,
                   size = (rtmax-rtmin)),
    )+
    scale_size_area(max_size = 5)+
    xlim(c(0,800))+
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
