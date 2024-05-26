#' xcms_get_dda_scan_stimulate
#' stimulate DDA cycle and assign ms2 to feature, now just support single file
#' @param xcms.scan xcms.scan
#' @param feature_def feature_def
#' @param dynamic_time dynamic_time
#'
#' @return  xcms.xcms
#' @export
#'

xcms_get_dda_scan_stimulate <- function(xcms.xcms ,
                              dynamic_time = 60){

  xcms.xcms <- xcms_get_scan_Stat(xcms.xcms)
  xcms.scan <- fData(xcms.xcms)%>%
    dplyr::mutate(cycle_ion_count = NA,
                  cycle_ion_ms2_count = NA,
                  ms2_hit = F,
                  ms2_feature_id = NA,
                  temp_var = scan_id)%>%
    remove_rownames()%>%
    column_to_rownames("temp_var")
  feature_def <- featureDefinitions(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::mutate(temp_var = feature_id)%>%
    remove_rownames()%>%
    column_to_rownames("temp_var")



  xcms.scan.ms1 <- xcms.scan %>%
    dplyr::filter(msLevel ==1)
  xcms.scan.ms2 <-xcms.scan %>%
    dplyr::filter(msLevel ==2)
  feature_def.temp <- feature_def%>%
    dplyr::mutate(ms2_acquire_count = 0,
                  ms2_acquired = F)
  for (i in 1:nrow(xcms.scan.ms1)) {

    cycle.rt <- c(xcms.scan.ms1$ms1_group_rt[i],
                  xcms.scan.ms1$ms1_group_rt[i]+xcms.scan.ms1$cycle_time[i])

    cycle.feature <- feature_def.temp%>%
      dplyr::filter(peakRtMin < cycle.rt[1]&peakRtMax > cycle.rt[1])


    cycle.ms2 <-xcms.scan.ms2 %>%
      dplyr::filter(ms1_group %in% xcms.scan.ms1$ms1_group[i])

    if (nrow(cycle.feature)!=0&nrow(cycle.ms2)!=0) {

      mz.error <- abs(matrixSub(cycle.ms2$precursorMZ,cycle.feature$mzmed))/cycle.ms2$precursorMZ
      mz.min.error <- apply(mz.error, 1, which.min)
      mz.hit <- apply(mz.error, 1, min) < 1e-5

      cycle.ms2 <-cycle.ms2  %>%
        dplyr::mutate(ms2_hit = mz.hit,
                      ms2_feature_id = cycle.feature$feature_id[mz.min.error],
                      ms2_feature_id = case_when(ms2_hit~ms2_feature_id,
                                             !ms2_hit~NA))

      #message(i,";",sum(feature_def.temp$ms2_acquired))

    }

    cycle.feature <- cycle.feature %>%
      dplyr::mutate(ms2_acquired = feature_id %in% cycle.ms2$ms2_feature_id)

    feature_def.temp[rownames(cycle.feature),]  <- feature_def.temp[rownames(cycle.feature),]%>%
      dplyr::mutate(ms2_acquired = case_when(ms2_acquired~T,
                                             !ms2_acquired ~cycle.feature$ms2_acquired,
                                             T~ms2_acquired)) %>%
      dplyr::mutate(peakRtMin= case_when(ms2_acquired ~  cycle.rt[1] + dynamic_time,
                                         T~peakRtMin),
                    ms2_acquire_count = case_when(
                      !ms2_acquired ~ ms2_acquire_count + 1,
                      ms2_acquired ~ ms2_acquire_count + 1,
                      T~ms2_acquire_count))

    xcms.scan[rownames(cycle.ms2),] <- xcms.scan[rownames(cycle.ms2),] %>%
      dplyr::mutate(ms2_feature_id = cycle.ms2$ms2_feature_id,
                    ms2_hit = cycle.ms2$ms2_hit)

    xcms.scan[xcms.scan.ms1$scan_id[i],]$cycle_ion_count <- nrow(cycle.feature)
    xcms.scan[xcms.scan.ms1$scan_id[i],]$cycle_ion_ms2_count <- sum(cycle.feature$ms2_acquired)



  }
  feature_def$ms2_acquire_count <-feature_def.temp$ms2_acquire_count
  feature_def$ms2_acquired <- feature_def.temp$ms2_acquired

  featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(feature_def)
  rownames(xcms.scan) <- rownames(fData(xcms.xcms))
  fData(xcms.xcms) <- xcms.scan
  return(xcms.xcms)


}




#' xcms_get_dda_ms2_assignment
#'
#'  match xcms scan of msLevel 2 to feature Definitions
#'  based on precursorMZ, retentionTime in `xcms.scan`;
#'  peakRtMin, peakRtMax, feature_id in `featuredef`
#'
#' @param xcms.xcms xcms
#' @param featuredef featuredef
#'
#' @return xcms
#' @export
#'

xcms_get_dda_ms2_assignment <- function(xcms.xcms){

  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  featuredef <- featureDefinitions(xcms.xcms)%>%as.data.frame()
  match.df <- match_mz_rt(featuredef$mzmed,featuredef$rtmed,
                          xcms.scan$precursorMZ,
                          xcms.scan$retentionTime)
  match.df <- match.df%>%
    dplyr::mutate(featuredef[ion1,],
                  scan_id = xcms.scan$scan_id[ion2],
                  scan_rt = xcms.scan$retentionTime[ion2])%>%
    dplyr::filter(scan_rt < peakRtMax&scan_rt>peakRtMin )%>%
    dplyr::group_by(scan_id)%>%
    dplyr::slice_min(mz.error)



  xcms.scan$ms2_feature_id <- sapply(
    1:nrow(xcms.scan),
    function(i){
      fdf <- match.df[match.df$ion2==i,]
      fid <- fdf$feature_id[which.min(fdf$mz.error)]
      if (length(fid)==0) {
        return(NA)
      }
      return(fid)
    }
  )
  xcms.scan$ms2_hit <- !is.na(xcms.scan$ms2_feature_id)

  featuredef$ms2_scan_id <- sapply(featuredef$feature_id,
                                   function(fid){
                                     xcms.scan$scan_id[which(xcms.scan$ms2_feature_id==fid)]
                                   })

  fData(xcms.xcms) <- xcms.scan
  featureDefinitions(xcms.xcms) <- featuredef%>%
    dplyr::mutate(ms2_acquire_count = sapply(ms2_scan_id,length) ,
                  ms2_acquired= ms2_acquire_count>0)%>%
    S4Vectors::DataFrame()


  return(xcms.xcms)



}


plot_xcms_dda_acquisition <- function(xcms.xcms) {

  if (!"ms2_feature_id" %in%colnames(xcms.xcms@featureData@data)) {
    xcms.xcms <- xcms_get_dda_ms2_assignment(xcms.xcms)
  }
  xcms.scan.ms2 <- fData(xcms.xcms)%>%
    dplyr::filter(msLevel == 2)

  feature_def <- featureDefinitions(xcms.xcms)%>%
    as.data.frame()

  ggplot()+
    geom_segment(data = feature_def,
                 aes(x = peakRtMin , xend = peakRtMax,
                     y= mzmed , yend = mzmed , col = ms2_acquired))+
    geom_point(data = xcms.scan.ms2, aes(x = retentionTime , y = precursorMZ,
                                         fill = ms2_hit),pch = 21)+
    scale_x_continuous(breaks = seq(0,max(xcms.scan.ms2$retentionTime),60))+
    #xlim(c(0,150))+
    #ylim(c(0,1200))+
    labs(x = "Acquisition windows",
         y = "Mz",
         col = "MS2 acquired",
         fill = "MS2 Hit")+
    theme_bw()


}



plot_xcms_dda_cycle_stimulate <- function( xcms.xcms,
                            topn =20){


  xcms.scan <- fData(xcms.xcms)

  plot.data <- xcms.scan%>%
    dplyr::filter(msLevel==1)%>%
    dplyr::mutate(cycle_ion_not_acqed_count = cycle_ion_count-cycle_ion_ms2_count)%>%
    pivot_longer(c(cycle_ion_ms2_count, cycle_ion_not_acqed_count),
                 names_to = "type",
                 values_to = "count")%>%
    dplyr::mutate(type = factor(type,level = c("cycle_ion_not_acqed_count",
                                               "cycle_ion_ms2_count")))

  ggplot(plot.data)+
    geom_bar(aes(x = retentionTime ,
                 y = count,
                 fill = type),
             alpha = 0.5,
             width = 2,
             stat = "identity",
             position = "stack",
             show.legend = F)+
    #geom_text(aes(x = max(plot.data$retentionTime)*0.8,
    #              y = max(plot.data$count)*0.9,
    #              label = "Features"),
    #          hjust=0,
    #          col = "#F7857D",
    #          check_overlap = T)+
    #geom_text(aes(x = max(plot.data$retentionTime)*0.8,
    #              y = max(plot.data$count)*0.8,
    #              label = "MS2 Acquired"),
    #          hjust=0,
    #          col = "#3EC6C9",
    #          check_overlap = T)+
    labs(x = "Retention Time",
         y = "Count in single DDA cycle")+
    ggbreak::scale_y_break(breaks = c(20,20),
                           scales = 0.5,
                           space = 0.1)+
    theme_bw()->p

  p


}


plot_xcms_dda_feature_stat <- function(xcms.xcms){

  feature_def <- featureDefinitions(xcms.xcms)%>%
    as.data.frame()
  ggplot(feature_def,
         aes( x = 1 , fill = ms2_acquired))+
    geom_bar(show.legend = F)+
    geom_text(aes(label =  after_stat(count)),
              stat = "count")+
    coord_polar(theta = "y")+
    theme_void()->p1

  ggplot(feature_def,
         aes( y = peakWidth,
                          x = ms2_acquired,
                          col = ms2_acquired))+
    geom_boxplot(alpha = 1,outlier.size = 0,
                 show.legend = F)+
    geom_jitter(size = 0.1,alpha = 0.3,
                show.legend = F)+
    labs(x = "MS2 acquired" , y = "Peak width")+
    ylim(c(0,50))+
    theme_classic()+
    theme(axis.text.x  = element_blank(),
          axis.title.x = element_blank())->p2

  ggplot(feature_def,
         aes( y = log10(peakMaxo),
              x = ms2_acquired,
              col = ms2_acquired))+
    geom_boxplot(alpha = 1,outlier.size = 0,
                 show.legend = F)+
    geom_jitter(size = 0.1,alpha = 0.3,
                show.legend = F)+
    labs(x = "MS2 acquired" , y = "Log10 intensity")+
    theme_classic()+
    theme(axis.text.x  = element_blank(),
          axis.title.x = element_blank())->p3


  ggplot(feature_def)+
    stat_density2d(aes(x = log10(peakMaxo),
                       y = peakWidth,
                       col = ms2_acquired),
                   n = 100)+
    labs(col = "MS2 acquired" ,
         y = "Peak width",
         x = "Log10 intensity")+
    #xlim(c(4,8))+
    ylim(c(0,50))+
    theme_bw()+
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(fill = "transparent"))->p4


  (p1+(p2+p3))/p4+
    plot_layout(widths = c(2, 1,1),
                heights = c(1, 2))

}

plot_xcms_dda_cycle_stat <- function(xcms.xcms){


  xcms.scan <- fData(xcms.xcms)
  xcms.scan.ms1 <- dplyr::filter(xcms.scan,msLevel == 1)
  ggplot(xcms.scan.ms1)+
    geom_histogram(aes(y = ms2_count),
                   binwidth = 1,
                   fill = "#1E80DE",col = "white")+
    labs(y = "MS2 count",x = "MS1 cycle count")+
    #scale_y_reverse()+
    scale_x_reverse()+
    ggsci::scale_color_gsea()+
    theme_bw()+theme(
      plot.margin =unit(c(0,0,0,0),"inch")
      )->p1

  ggplot(xcms.scan.ms1)+
    geom_jitter(aes(x = retentionTime ,
                   y = ms2_count,
                   col = log10(totIonCurrent)),
                width = 0,
                height = 0.1,
                alpha = 0.8,
                size = 1)+
    labs(y = "",col = "Log10 TIC")+
    ggsci::scale_color_gsea()+
    theme_bw()+
    theme(axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          plot.margin =unit(c(0,0,0,0),"inch"),
      axis.line.y.left  = element_blank())->p2


  p1+p2+plot_layout(widths = c(1,4))

}


