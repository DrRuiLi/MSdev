#' @title  get_features_from_xcms
#' @description extract feature data from xcms::XCMSnExp,
#'  calculate RSD of QC and Sample
#'  ( note this rely on character "QC" and "Sample" in `sampleNames(xcms.xcms)` )
#' @param xcms.xcms
#'
#' @return a SummarizedExperiment subject
#' @export
#'
#' @examples
get_features_from_xcms <- function(xcms.xcms,missing = NA){

  xcms.sum <- quantify(xcms.xcms,missing = missing )
  feature.def <- SummarizedExperiment::rowData(xcms.sum)%>%
    tibble::as_tibble()

  feature.matrix <- SummarizedExperiment::assay(xcms.sum)
  rsd <- function(x){sd(x,na.rm =  T)/mean(x , na.rm = T)}
  feature.matrix.qc <- feature.matrix[,which(grepl("QC",colnames(feature.matrix)))]
  feature.matrix.sample <- feature.matrix[,which(grepl("Sample",colnames(feature.matrix)))]
  if(sum(grepl("QC",colnames(feature.matrix)))>1){

    feature.def$qc_rsd <- apply(feature.matrix.qc, 1, rsd)
  }else(

    feature.def$qc_rsd <- 0
  )

  if(sum(grepl("Sample",colnames(feature.matrix)))>1){

    feature.def$sample_rsd <- apply(feature.matrix.sample, 1, rsd)
  }else(

    feature.def$sample_rsd <- 0
  )
  feature.def$med_intensity <- apply(feature.matrix , 1 ,median,na.rm =T)
  SummarizedExperiment::rowData(xcms.sum) <-feature.def
  return(xcms.sum)
}


get_chrom_peaks_shape_score <- function(chrom,
                                        peak.id = chrom@chromPeakData@rownames){
  peak.id = chrom@chromPeakData@rownames
  peak.id <- peak.id[1]
  peaks.data <- chromPeaks(chrom)[peak.id,,drop = F]
  rtime <- rtime(chrom)
  int <- intensity(chrom)

  rtime <- rtime[!is.na(int)]
  int <- int[!is.na(int)]

  int.fit <-peak.gasssian.fit(rtime,
                    peak.apex.intensity = peaks.data[1,"maxo"],
                    peak.apex.rt = peaks.data[1,"rt"],
                    peak.half.width = min(peaks.data[1,"rtmax"]-peaks.data[1,"rt"],
                                          peaks.data[1,"rt"]-peaks.data[1,"rtmin"])/2)
  int[is.na(int)] <- 0


  #sum(abs(int-int.fit)/sum(int.fit))
  #cor(int,int.fit)
  #sqrt(mean((int-int.fit)^2))/mean(int.fit)
  r2 <- 1-sum((int-int.fit)^2)/sum( (int-mean(int))^2 )
  r2.adj <- 1-(1-r2)*(length(int)-1)/(length(int)-2)
  r2.adj
}


get_xchrom_peak_score <- function(xchrom){

  peak.info <- chromPeaks(xchrom)
  peaks.no <- nrow(peak.info)


  for (i in 1:peaks.no) {


  }


}

get_xchroms_peaks_count <- function(xchroms){

  peaks.info <- chromPeaks(xchroms)%>%
    as.data.frame()%>%
    dplyr::group_by(row,column)%>%
    dplyr::mutate(peaks.no = n(),
                  idx = (row -1)*dim(xchroms)[2]+column )%>%
    dplyr::distinct(row,column ,peaks.no,idx)
  peaks.count.matrix <-rep(0,length(xchroms))
  peaks.count.matrix[peaks.info$idx] <- peaks.info$peaks.no
  peaks.count.matrix <- matrix(peaks.count.matrix,
                               nrow = dim(xchroms)[1],
                               byrow = T)%>%
    `colnames<-`(1: dim(xchroms)[2])
  return(peaks.count.matrix)
}

#' @title get_xcms_peaks_chrom
#' @description
#' extract chromatograph from XCMSnExp,
#' if `all.sample` = F, only the samples, in which given peaks.id are detected will be return,
#' else extract from all samples
#'
#' @param xcms.xcms
#' @param peaks.id
#' @param all.sample
#' @param rt one of c("all","identity","expand")
#'
#' @return XChromatograms
#' @import xcms
#' @export
#'
#' @examples
get_xcms_peaks_chrom <- function(xcms.xcms,
                                 peaks.id ,
                                 all.sample =F,
                                 rt.range = "expand"){

  peaks.data <- xcms::chromPeaks(xcms.xcms)
  if(is.numeric(peaks.id)) { peaks.id <-rownames(peaks.data)[peaks.id]}
  peaks.data <- peaks.data[peaks.id,,drop=F]
  peaks.data[,c("rtmin","rtmax")] <- switch (rt.range,
                                             "all" = c(min(rtime(xcms.xcms)),max(rtime(xcms.xcms))),
                                             "expand" = apply(peaks.data[,c("rtmin","rtmax"),drop =F],1,expand_range,add = 15)%>%t,
                                             "identity" = peaks.data[,c("rtmin","rtmax"),drop =F]
  )

  if (all.sample){
    x.chrom <- extract_chrom(xcms.xcms,
                             mzr = peaks.data[,c("mzmin","mzmax")],
                             rtr = peaks.data[,c("rtmin","rtmax")],
                             sample = sampleNames(xcms.xcms))
  }else{
    peaks.id.split <- split(rownames(peaks.data),
                              f = peaks.data[,"sample"])
    x.chrom.split <- bplapply(peaks.id.split,
                              FUN = function(x,xcms.x,peaks.data){

      x <- peaks.data[x,,drop =F]
      extract_chrom(xcms.x,
                    mzr = x[,c("mzmin","mzmax")],
                    rtr = x[,c("rtmin","rtmax")],
                    sample = x[,"sample"])
    },xcms.x=xcms.xcms ,peaks.data=peaks.data,
    BPPARAM = SerialParam( progressbar = T))

    x.chrom <- XChromatograms(unlist(x.chrom.split))

    x.chrom
  }



  return(x.chrom)
}



#' extract_chrom
#'
#' @param xcms.xcms
#' @param rtr
#' @param mzr
#' @param sample
#'
#' @return
#' @export
#'
#' @examples
extract_chrom <- function(xcms.xcms,
                          rtr,
                          mzr,
                          sample){
  xcms.xcms <- MSnbase::filterFile(xcms.xcms,sample)
  xcms::chromatogram(xcms.xcms,
                     rt = rtr,
                     mz = mzr,
                     BPPARAM = BiocParallel::SerialParam()
                     )
}




get_xcms_feature_chrom <- function(xcms.xcms,
                                    feature.id,
                                    sample = "maxo",
                                    rt = "expand"){

  features.data <- featureDefinitions(xcms.xcms)
  features.val <- featureValues(xcms.xcms,missing = "rowmin_half")
  if(is.numeric(feature.id)) { feature.id <-rownames(features.data)[feature.id]}
  features.data <- features.data[feature.id,,drop=F]
  features.val <- features.val[feature.id,,drop =F]
  if (sample == "maxo")  {
    xcms.sub <- MSnbase::filterFile(xcms.xcms,
                                    which.max(features.val))
  }else if (sample == "all")  {
    xcms.sub <- xcms.xcms
  }else {
    xcms.sub <- MSnbase::filterFile(xcms.xcms, sample )
  }

    bp <-BiocParallel::SerialParam(progressbar = F)


  rtr <- switch (rt,
                 "all" = c(min(rtime(xcms.sub)),max(rtime(xcms.sub))),
                 "expand" = t(apply(features.data[,c("peakRtMin","peakRtMax"),drop =F],1,expand_range,add = 15)),
                 "identity" = features.data[,c("peakRtMin","peakRtMax"),drop =F]

  )
  x.chrom <-  xcms::chromatogram(xcms.sub,
                                 mz = features.data[,c("peakMzMin","peakMzMax")]%>%as.matrix(),
                                 rt = rtr,
                                 aggregationFun = "max",
                                 BPPARAM = bp)
  return(x.chrom)
}


get_chrom_peaks_gaussian_fit <- function(xchrom){

  peaks.info <- chromPeaks(xchrom)[1,,drop = F]
  peaks.data <- get_chroms_data(xchrom)
  nls(formula = intensity ~ gaussian_functioin(rt ,a,b,c),
      data = peaks.data,control = nls.control(warnOnly = T),
      start = list(a = peaks.info[,"maxo"],
                   b = peaks.info[,"rt"],
                   c = mean(diff(peaks.info[,c("rtmin","rtmax")])))) -> gaussian.fit

}


#' @title get_chroms_data
#' @description extract chomatogram data to a data.frame
#' @param xchrom
#'
#' @return
#' @export
#'
#' @examples
get_chroms_data <- function(xchrom){

  .extract.chrom <- function(i,j){
    this.chrom <- xchrom[i,j]
    data.frame(
      rt = rtime(this.chrom),
      intensity =intensity(this.chrom),
      row =i,col = j
    )
  }
  if (class(xchrom) %in% c("XChromatogram","Chromatogram")) {
    xchrom <- XChromatograms(list(xchrom))
  }
  bp.matrix <- expand.grid(1:nrow(xchrom),1:ncol(xchrom))
  xchrom.data <- BiocParallel::bpmapply(.extract.chrom,
                                        bp.matrix[,1],bp.matrix[,2],
                         BPPARAM = BiocParallel::SerialParam(progressbar = F),SIMPLIFY=F)%>%
    do.call("rbind",.)

  return(xchrom.data)

}


#' XChromatograms_rt_unit
#'
#'  change rtime units, in some situation (such as SRM data from Thermo), rtime are recorded with unit "m",
#'  this will lead to error when findChrompeaks
#'
#' @param xchroms `XChromatograms` or `MChromatograms` object
#' @param unit_to "s" or "m", "s": rtime*60; "m": rtime/60
#'
#' @return
#' @export
#'
#' @examples
#'
XChromatograms_rt_unit <- function(xchroms,unit_to = "s"){


  unit.mulit <- switch(unit_to,
            "s" = 60,
            "m" = 1/60)
  rtime.max <- max(rtime(xchroms[1,1]))
  f <- function(i,xchroms,unit.mulit){

    i.position <- which(matrix(1:length(xchroms),
                               nrow = dim(xchroms)[1],
                               byrow = T)==i,arr.ind = T)
    x <- i.position[1]
    y <- i.position[2]
    xchroms[x,y]@rtime <- xchroms[x,y]@rtime * unit.mulit
    xchroms[x,y]

  }

  xchroms.trans <- bplapply(1:length(xchroms),FUN = f,
                              xchroms = xchroms,
                              unit.mulit = unit.mulit,
                              BPPARAM = BatchtoolsParam(progressbar = T,
                                                        registryargs = batchtoolsRegistryargs(packages = c("MSnbase"))))
  xchroms.trans <- XChromatograms(xchroms.trans,
                                    nrow = dim(xchroms)[1],byrow = T)

  message( "rtime value ", round(rtime.max,0), " change to ", round(max(rtime(xchroms.trans[1,1])),0))

  return(xchroms.trans)


}



#' XChromatograms_fill_2point
#'
#' when xcms::findChromPeaks(), if any Chromatogram contain less than 2 point, this will lead to error
#'
#' @param xchroms
#'
#' @return
#' @export
#'
#' @examples
XChromatograms_fill_2point <- function(xchroms){

  rt.point <- sapply(1:length(xchroms), function(x){
    x.position <- arrayInd(x,.dim = dim(xchroms))
    length(rtime(xcms.xcms[x.position[1],x.position[2]]))
  })

  tofill <- which(rt.point <2)%>%
    arrayInd(.dim = dim(xchroms))

  f <- function(chrom){
    if (length(rtime(chrom)) ==1) {
      chrom@rtime <- c(chrom@rtime , chrom@rtime +1)
      chrom@intensity <- c(chrom@intensity ,0)

    }
    chrom
  }
  xchroms.filled <- bplapply(unlist(xchroms),f ,BPPARAM = BatchtoolsParam(progressbar = T,
                                                             registryargs = batchtoolsRegistryargs(packages = c("MSnbase"))))
  xchroms.filled <- XChromatograms(xchroms.filled,
                                  nrow = dim(xchroms)[1])
  return(xchroms.filled)

}


#' plot_XChromatograms
#'
#' @param xchrom
#' @param norm
#' @param move
#'
#' @return
#' @export
#'
#' @examples
plot_XChromatograms <- function(xchroms , norm = T,move = T,column.group = 1:ncol(xchroms)){


  if (norm) {
    xchroms <- normalise(xchroms)
    chrom.data <- get_chroms_data(xchroms)%>%
      dplyr::mutate(peaks.origin = paste0("peak_",num2str(row),"_",column.group[col]),
                    peaks.origin = factor(peaks.origin,level = unique(peaks.origin)))%>%
      dplyr::group_by(peaks.origin)%>%
      dplyr::mutate(peaks.idx =cur_group_id(),
                    intensity = intensity*100
      )%>%
      dplyr::ungroup()
  }else{
    chrom.data <- get_chroms_data(xchroms)%>%
      dplyr::mutate(peaks.origin = paste0("peak_",num2str(row),"_",column.group[col]),
                    peaks.origin = factor(peaks.origin,level = unique(peaks.origin)))%>%
      dplyr::group_by(peaks.origin)%>%
      dplyr::mutate(peaks.idx =cur_group_id(),
                    #intensity = case_when(is.na(intensity)~ 0 ,
                    #                      T~intensity)
      )%>%
      dplyr::ungroup()

  }


    if (move) {
      chrom.data <- chrom.data%>%
        dplyr::mutate(rt = rt +peaks.idx*3,
                      intensity = intensity+peaks.idx*3)
    }







  ggplot(chrom.data)+
    geom_line(aes(x = rt , y = intensity , col = peaks.origin),
              linewidth = 0.5,alpha = 0.8)+
    scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(chrom.data$peaks.idx))))+
    labs(x = "Retention time", y = "Intensity")+
    theme_bw()



}


#' @title featureDefinitions_PeakSta
#' @description extract features' median rt, sn and maxo,
#' `xcms::featureDefinitions()` return a `DataFrame`, in which rtmin, rtmax, rtmed was median of `xcms::chromPeaks()$rt`,
#' but not the median range of peaks. peakRtMin, peakRtMax, peakSN, peakMaxo are median of all peaks in a feature
#'
#' @param xcms.xcms
#'
#' @return
#' @export
#'
featureDefinitions_PeakSta<- function(xcms.xcms){

  feature.def <- featureDefinitions(xcms.xcms)
  feature.val <- featureValues(xcms.xcms)
  peaks.data <- chromPeaks(xcms.xcms)

  .xcmsPeakDataMed <- function(x,peaks.data,key = "rtmax",fun = "median"){
    x.peaks.data <- peaks.data[x,,drop=F]
    x.peaks.data <- x.peaks.data[!grepl(pattern = "CPM",x = rownames(x.peaks.data)),,drop=F]
    peak.key.value <- x.peaks.data[,key]
    eval(call(fun,peak.key.value,na.rm =T))
  }


  feature.def$peakRtMin <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,key = "rtmin",fun = "min")
  feature.def$peakRtMax <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,key = "rtmax",fun = "max")
  feature.def$peakWidth <- feature.def$peakRtMax-feature.def$peakRtMin
  feature.def$peakMzMin <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"mzmin",fun = "min")
  feature.def$peakMzMax <- sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"mzmax",fun = "max")
  feature.def$peakSN <-  sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"sn")
  feature.def$peakMaxo <-  sapply(feature.def$peakidx,.xcmsPeakDataMed,peaks.data,"maxo")
  feature.def$polarity <- polarity(xcms.xcms)%>%unique()
  feature.def.df <- as.data.frame(feature.def)%>%
    dplyr::mutate(feature_id = rownames(.),
                  .before = mzmed)
  feature.def.df
}


get_xcms_feature_def_stat <- function(xcms.xcms){

  xcms.def <- featureDefinitions_PeakSta(xcms.xcms)%>%
    DataFrame()
  featureDefinitions(xcms.xcms)<-xcms.def
  return(xcms.xcms)

}



get_xcms_feature_val_stat <- function(xcms.xcms) {

  xcms.pdata <- Biobase::pData(xcms.xcms)
  featureval <- featureValues(xcms.xcms)
  if ("sample.type" %in% colnames(xcms.pdata)) {
    qc.rsd <- featureval[,xcms.pdata%>%dplyr::filter(sample.type =="QC")%>%
                 dplyr::pull(sampleNames)]%>%
      apply(1,function(x){
        sd(x,na.rm =T)/mean(x,na.rm=T)
      })
    sample.rsd <- featureval[,xcms.pdata%>%
                               dplyr::filter(sample.type =="Sample")%>%
                 dplyr::pull(sampleNames)]%>%
      apply(1,function(x){
        sd(x,na.rm =T)/mean(x,na.rm=T)
      })
  }else{
    qc.rsd <- NA
    sample.rsd <- NA

  }

  fdf <- featureDefinitions(xcms.xcms)
  fdf$qc_rsd <- qc.rsd
  fdf$sample_rsd <- sample.rsd
  fdf -> featureDefinitions(xcms.xcms)
  return(xcms.xcms)
}


get_xcms_feature_stat <- function(xcms.xcms){
  xcms.xcms <- xcms.xcms %>%
    get_xcms_feature_def_stat()%>%
    get_xcms_feature_val_stat()
  return(xcms.xcms)
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
    dplyr::mutate(as.data.frame(chromPeakData(xcms.xcms)))%>%
    dplyr::filter(!is.na(maxo),
                  rtmax-rtmin <60,
                  !merged)
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
      apply(featureValues(xcms.xcms, value = "maxo"), 1, median, na.rm = T)

    xcms.process.type <-
      processHistory(xcms.xcms) %>% sapply(processType)
    xcms.findpeak.param <-
      processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]] %>%
      processParam()
    ### generate scale
    maxo.range <- quantile(log10((xcms.features$maxo)))
    peak.witdh.range <- quantile(xcms.features$peakWidth)

    ggplot(xcms.features) +
      geom_point(aes(
        x = rt,
        y = mz,
        col = log10(maxo),
        alpha = log10(maxo) / 10,
        size =peakWidth
      ),) +
      scale_size_continuous(breaks  = peak.witdh.range[2:4],
                            labels = round(peak.witdh.range[2:4]),
                            range = c(0,10)) +
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
      theme_bw()+
      theme(text = element_text(size = 8)) -> peaks.dis.plot
    peaks.dis.plot
    return(peaks.dis.plot)





  }


#' @title plot_xcms_feature_chromatogram
#' @description extract Chromatogram from xcms according to feature's mz range and plot
#' @param xcms.xcms
#' @param feature.id
#' @param sampleNames
#'
#' @return
#' @export
#'
#' @examples
plot_xcms_feature_chromatogram <- function(xcms.xcms ,feature.id, sampleNames =NULL ){

  ### select samples
  xcms.sample.info <- Biobase::pData(xcms.xcms)
  if (is.null(sampleNames)) {
    sampleNames <- xcms.sample.info$sampleNames
  }
  xcms.sample.info <- xcms.sample.info[sampleNames,,drop=F]
  if (length(sampleNames) > 5) {
    if (!is.null(xcms.sample.info$group)) {
      xcms.sample.info.sub <- xcms.sample.info%>%
        dplyr::group_by(group)%>%
        dplyr::slice_sample(n=1)
    }
  }else{
    xcms.sample.info.sub <- xcms.sample.info
  }
  xcms.sub <- filterFile(xcms.xcms,which(Biobase::sampleNames(xcms.xcms)%in% xcms.sample.info.sub$sampleNames))
  ### mz
  xcms.feature <- featureDefinitions(xcms.xcms)[feature.id,]
  feature.id<-rownames(xcms.feature)
  xcms.peaks <- chromPeaks(xcms.xcms)[xcms.feature$peakidx[[1]],,drop = F]
  mz.range <- c(min(xcms.peaks[,"mzmin"]),
                max(xcms.peaks[,"mzmax"]))
  rt.range <- c(min(xcms.peaks[,"rtmin"]),
                max(xcms.peaks[,"rtmax"]))
  xcms.chrom <- chromatogram(xcms.sub , mz = mz.range,
                            #rt =c(xcms.feature$rtmin,xcms.feature$rtmax),
                            #rt = rt.range,
                             adjustedRtime  =F)

  xcms.chrom.data <- get_intensity_rtime_df_from_XChromatogram(xcms.chrom)%>%
    dplyr::mutate(group = sample.name)

  ggplot(xcms.chrom.data)+
    geom_line(aes(x = rt,y = intensity , col = group))+
    xlim(c(min(rtime(xcms.sub)),max(rtime(xcms.sub))))+
    labs(col = "",x = "Retention time", y = "Intensity",
         title = paste0(feature.id),
         subtitle = paste0( "mz: ",round(mz.range[1],6),
                            " ~ ",round(mz.range[2],6),
                            "\nrt: ",round(rt.range[1],2),
                            " ~ ",round(rt.range[2],2) ))+
    theme_bw()+
    theme(text = element_text(size = 8))


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
  xcms.scans <- get_xcms_scan_Stat(xcms.xcms)%>%
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
  ms2_scans_table <- table(xcms.peaks$ms2_scans_no)
  ggplot(xcms.peaks)+
    geom_jitter(aes(x = rt, y = ms2_scans_no, col = log10(maxo)),
                 size = 0.6
    )+
    #geom_hline(yintercept = 7)+
    geom_violin(aes( x = max(rt)*1.2 , y =ms2_scans_no),width = diff(range(xcms.peaks$rt))*0.1)+
    geom_text(aes(x =  max(rt)*1.3,y = 0,label = ms2_scans_table["0"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 1,label = ms2_scans_table["1"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 2,label = ms2_scans_table["2"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 3,label = ms2_scans_table["3"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 4,label = ms2_scans_table["4"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.3,y = 5,label = ms2_scans_table["5"]),size = 2.67,hjust = 0)+
    geom_text(aes(x =  max(rt)*1.4,y = 5,label = ""),size = 2.67,)+
    labs(title = plot.title,
         subtitle = paste0("ppm = ",xcms.findpeak.param@ppm,
                           "; SN = ",xcms.findpeak.param@snthresh,
                           "; prefilter = (",paste0(xcms.findpeak.param@prefilter,collapse = ","),")" ,"\n",
                           sum(xcms.peaks$ms2_scans_no > 0)," / ",length(xcms.peaks$ms2_scans_no),
                           " ( ",sprintf("%.2f",sum(xcms.peaks$ms2_scans_no > 0)/length(xcms.peaks$ms2_scans_no)*100),"% )"),
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



plot_xcms_ms2_distribution <- function(xcms.xcms,plot.title = "MS2 Precursor distribution" ){

 scan.data <- fData(xcms.xcms)%>%
   dplyr::filter(msLevel==2)

 ms1.rt <- fData(xcms.xcms)%>%
   dplyr::filter(msLevel==1)%>%
   dplyr::pull(retentionTime)

 ggplot(scan.data)+
   geom_vline(xintercept = ms1.rt,linewidth = 0.05,col = "black")+
   geom_point(aes(x = retentionTime,y= precursorMZ,
                  col = log10(precursorIntensity)),
   )+
   labs(title = plot.title,
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
 peaks.dis.plot
 open_ggplot_win(peaks.dis.plot,width = 25,height = 5)

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

#' @title plot_xcms_peaks_Chromatogram
#' @description extract EIC according to peaks' mzrange and rtrange,
#' note that if multiple sample in xcms object, only first sample will be extracted
#'
#' @param xcms.xcms
#' @param peak_ids
#' @param rt_expand foldchange to expand rt range
#'
#' @return
#' @export
#'
#' @examples
plot_xcms_peaks_Chromatogram <- function(xcms.xcms,peak_id,rt = "expand"){

  peaks.data <- chromPeaks(xcms.xcms)[peak_id,,drop = F]
  peak_id <- rownames(peaks.data)
  mz.range <- c(peaks.data[,c("mzmin","mzmax")])
  rt.range <- c(peaks.data[,c("rtmin","rtmax")])
  xcms.chrom <- get_xcms_peaks_chrom(xcms.xcms,peaks.id = peak_id,rt = rt)
  chrom.data <- get_chroms_data(xcms.chrom)%>%
    dplyr::mutate(fill = rt > min(rt.range)&rt <max(rt.range),
                  fit = gaussian_functioin(rt,
                                          peak.apex.intensity = peaks.data[1,"maxo"],
                                          peak.apex.rt = peaks.data[1,"rt"],
                                          peak.half.width = min(peaks.data[1,"rtmax"]-peaks.data[1,"rt"],
                                                                peaks.data[1,"rt"]-peaks.data[1,"rtmin"])/2
                  ))

  ggplot(chrom.data)+
    geom_line(aes(x = rt,y = intensity),linetype = 1)+
    geom_area(aes(x = rt,y = intensity, fill = fill))+
    geom_point(aes(x = rt, y = fit))+
    scale_fill_manual(values = c("FALSE" = "transparent","TRUE" = "grey"))+
    labs(title = paste0(peak_id),
         subtitle = paste0("mz:",paste0(sprintf("%.5f",mz.range),collapse = " - "), ";     ",
                           "rt:",paste0(sprintf("%.2f",rt.range),collapse = " - "),"\n",
                           "mz error = ",sprintf("%.2f",mean(diff(mz.range)/mz.range)*1e6)," ppm;     ",
                           "peak width = ", sprintf("%.2f",diff(rt.range)),"\n",
                           "shape score = ",get_chrom_peaks_shape_score(xcms.chrom[1,1])
                           ),
         x = "Retention time")+
    guides(fill = "none")+
    theme_bw()+
    theme(text = element_text(size = 8))



}



chromPeaks_Sta <- function(xcms.xcms){

  xcms.peaks.info <- xcms::chromPeaks(xcms.xcms)
  xcms.peaks.xchroms <- get_xcms_peaks_chrom(xcms.xcms,
                                             1:nrow(xcms.peaks.info))




}


#' @title xcmsProcessingMS1
#' @description Import `msDataFiles`, filter `ion_mode`, find peaks using `centWaveParam`, correct RT, group peaks using `peaksGroup`, fill peaks by xcms at MS1 Level
#' @param msDataFiles `char` ms file (full) paths
#' @param ion_mode to filter ion_mode, 1: positive, 0: negative, import when scans with both pos and neg
#' @param peaksGroup `vector` to PeakGroupsParam(sampleGroups), should contain "QC"
#' @param centWaveParam xcms::CentWaveParam()
#'
#' @return
#' @export
#' @import xcms
#' @examples
xcmsProcessingMS1 <- function(msDataFiles,ion_mode = NA,peaksGroup =NA,
                              centWaveParam = xcms::CentWaveParam(ppm = 20,
                                                                  peakwidth = c(5,20),
                                                                  snthresh = 10,
                                                                  fitgauss = T,
                                                                  verboseColumns = T,
                                                                  prefilter = c(3,100))){
  if (!length(msDataFiles)) {
    return(NA)
  }
  xcms.xcms <-  MSnbase::readMSData(msDataFiles, mode = "onDisk")
  if (is.na(ion_mode)) {
    ion_mode <- polarity(xcms.xcms )%>%unique()
    if (length(ion_mode)!=1) {
      stop("MS1 scans contain both positive and negative, please check")
    }
  }
  if (all(is.na(peaksGroup))) {
    peaksGroup<- rep("A",length(msDataFiles))

  }
  xcms.xcms <- ProtGenerics::filterPolarity(xcms.xcms , ion_mode)
  message(Sys.time()," Find peaks...")
  xcms.xcms<-xcms::findChromPeaks(xcms.xcms,
                            param = centWaveParam,
                            BPPARAM  = BiocParallel::SerialParam(progressbar = T))
  #mpp <- xcms::MergeNeighboringPeaksParam(expandRt = 2.5,minProp = 0.5)
  #xcms.xcms <- xcms::refineChromPeaks(xcms.xcms, mpp,
  #                                    BPPARAM  = BiocParallel::SerialParam(progressbar = T))
  message(Sys.time()," Adjust RT...")
  peak.density.param <- xcms::PeakDensityParam(sampleGroups = peaksGroup,
                                         minFraction = 0.4,bw = 30,
                                         binSize = 0.015)
  xcms.xcms <- xcms::groupChromPeaks(xcms.xcms,param = peak.density.param)



  if (length(oligoClasses::sampleNames(xcms.xcms))>1) {
    if (sum(peaksGroup=="QC") <2 ) {
      #rt.adjust.param <- ObiwarpParam()
     # xcms.xcms <- adjustRtime(xcms.xcms,param = rt.adjust.param)
    }else{
      ### adjust based on QC
      rt.adjust.param <- PeakGroupsParam(minFraction = 0.4,
                                          subset = which(peaksGroup == "QC"),
                                          subsetAdjust = "average",span = 0.4)
      xcms.xcms <- adjustRtime(xcms.xcms,param = rt.adjust.param)
    }
  }

  message(Sys.time()," Group peaks...")

  peak.density.param <- PeakDensityParam(sampleGroups =peaksGroup,
                                         minFraction = 0.4,bw = 30,
                                         binSize = 0.015)
  xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)
  xcms.xcms <- fillChromPeaks(xcms.xcms,param = FillChromPeaksParam())
  xcms.xcms <- get_xcms_feature_stat(xcms.xcms )
  return(xcms.xcms)




}

xcmsProcessingMRM <- function(msDataFiles, peaksGroup =NA,
                              centWaveParam ){

  xcms.mrm <- readSRMData(msDataFiles)
  xcms.peaks <- findChromPeaks()


}







matchSpectra_Features <- function(xcmsFeatureDef, spec){

  .matchSP <- function(x,xcmsFeatureDef,
                       mz_ppm = 10,
                       rt_tol = 10){
    mz <- x[["precursorMz"]]%>%as.numeric()
    rt <- x[["rtime"]]%>%as.numeric()
    mzError <- abs((mz - xcmsFeatureDef$mzmed)/mz*1e6)
    rtError <- abs((rt- xcmsFeatureDef$rtmed)/rt)
    feature_id <- rownames(xcmsFeatureDef)[mzError < mz_ppm &rtError < rt_tol]
    #ifelse(length(feature_id)==0, NA,feature_id)


  }
  spec.data <- as.data.frame(Spectra::spectraData(spec))
  feature_id <- apply(spec.data, 1,.matchSP,xcmsFeatureDef)
  spec$feature_id <- feature_id
  spec

}






#' @title plot_xcms_feature_intensity
#' @description plot feature's intensity, ordered by `Biobase::pData(xcms.xcms)$analysis.time.positive` or
#'  `Biobase::pData(xcms.xcms)$analysis.time.negative`
#'
#' @param xcms.xcms
#' @param feature_id_to_show
#'
#' @return
#' @export
#'
#' @examples
plot_xcms_feature_intensity <- function(xcms.xcms , feature_id_to_show ){

  ion_mode <- unique(fData(xcms.xcms)$polarity)
  if (ion_mode==1) {
    sample.info <- Biobase::pData(xcms.xcms)%>%
      dplyr::arrange(analysis.time.positive)%>%
      dplyr::mutate(sample.type = factor(sample.type,levels = c("Blank","QC","Sample")),
                    injecton.order = 1:nrow(.))
  }else{

    sample.info <- Biobase::pData(xcms.xcms)%>%
      dplyr::arrange(analysis.time.negative)%>%
      dplyr::mutate(sample.type = factor(sample.type,levels = c("Blank","QC","Sample")),
                    injecton.order = 1:nrow(.))
  }
  features <- featureValues(xcms.xcms)%>%
    as.data.frame()%>%
    rownames_to_column("feature_id")%>%
    dplyr::filter(feature_id %in% feature_id_to_show )%>%
    dplyr::select(sample.info$sampleNames)%>%as.numeric()

  sample.info$intensity <- features
  sample.info$intensity[is.na(sample.info$intensity )] <- 0
  ggplot(sample.info,aes(x = injecton.order , y = intensity , col = sample.type,na.rm =T))+
    geom_point(size = 0.5)+
    scale_color_manual(values = c("grey","#66CAB7","#EE8E5B"))+
    theme_bw()+
    theme(text = element_text(size = 8))


}




get_xcms_scan_Stat <- function(xcms.xcms){

  xcms.fdata <-Biobase::fData(xcms.xcms)%>%
    dplyr::mutate(fileStr = num2str(fileIdx),
                  spStr = num2str(spIdx)
                  )%>%
    dplyr::group_by(fileStr)%>%
    dplyr::mutate(x = 2-msLevel,
                  ms1_no = cumsum(x))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(ms1_no_str = num2str(ms1_no))%>%
    dplyr::group_by(fileIdx)%>%
    dplyr::arrange(fileIdx,retentionTime  )%>%
    dplyr::mutate(scan_time = c(diff(retentionTime),0),
                  ms1_group = paste0(fileStr,"_",ms1_no_str))%>%
    dplyr::group_by(ms1_group)%>%
    dplyr::mutate(ms2_count = sum(msLevel==2),
                  ms1_group_rt = min(retentionTime),
                  cycle_time = max(retentionTime)-min(retentionTime))%>%
    dplyr::ungroup()%>%
    dplyr::group_by(fileStr)%>%
    dplyr::mutate(cycle_time = c(diff(ms1_group_rt),0))%>%
    dplyr::group_by(ms1_group)%>%
    dplyr::mutate(cycle_time = max(cycle_time))%>%
    dplyr::mutate(scan_id = paste0(fileStr,"_",spStr))%>%
    dplyr::ungroup()%>%
    dplyr::select(scan_id,ms1_no,ms1_group,ms1_group_rt,
                  ms2_count,cycle_time,scan_time,everything(),
                  -c(x,fileStr,spStr,ms1_no_str))

  xcms.fdata
}


#' get_xcms_scan_feature_id
#'
#'  match xcms scan of msLevel 2 to feature Definitions
#'  based on precursorMZ, retentionTime in `xcms.scan`;
#'  peakRtMin, peakRtMax, feature_id in `featuredef`
#'
#' @param xcms.scan
#' @param featuredef
#'
#' @return
#' @export
#'
#' @examples
get_xcms_scan_feature_id <- function(xcms.scan,
                                     featuredef ){

  assign_ms2 <- function(pmz,prt,xcms.featuredef){


    mz.pass <- abs(xcms.featuredef$mzmed - pmz)/pmz <1e-5
    rt.pass <- (xcms.featuredef$peakRtMax > prt)&(xcms.featuredef$peakRtMin < prt)

    x <- xcms.featuredef$feature_id[mz.pass&rt.pass]
    if (length(x)==0) {
      return(NA)

    }
    return(paste0(x,collapse = ";"))

  }

  if (sum(xcms.scan$msLevel >1) > 10000) {
    bpp<-BiocParallel::SnowParam(progressbar = T)
  }else{
    bpp<-BiocParallel::SerialParam(progressbar = T)

  }

  xcms.scan <- xcms.scan%>%
    dplyr::filter(msLevel == 2)%>%
    dplyr::mutate(ms2_matched_feature =
                    BiocParallel::bpmapply(assign_ms2 ,
                             pmz = precursorMZ,
                             prt = retentionTime,
                             MoreArgs = list(xcms.featuredef = featuredef),
                             BPPARAM = bpp))

  return(xcms.scan)



}


plot_xcms_TIC <- function(xcms.xcms){


  xcms.pdata <- Biobase::pData(xcms.xcms)
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.scan <- xcms.scan%>%
    dplyr::mutate(tic = MSnbase::tic(xcms.xcms),
                  group = xcms.pdata$group[fileIdx])%>%
    dplyr::filter(msLevel==1)

  col.scale <- c("grey","#38C291",ggsci::pal_aaas()(10))
  names(col.scale) <- c("Blank","QC",
                        setdiff(unique(xcms.scan$group),c("Blank","QC")))

  ggplot(xcms.scan)+
    geom_line(aes(x = retentionTime , y = tic,col = group,group=fileIdx))+
    scale_color_manual(values = col.scale)+
    labs(title = "TIC",x = "Retention Time", y = "Intensity", col = "")+
    theme_classic()->p
  p



}

plot_xcms_adjustedRT <- function(xcms.xcms){


  xcms.pdata <- Biobase::pData(xcms.xcms)%>%
    dplyr::arrange(ExpTime)%>%
    dplyr::mutate(injection_order = 1:n())
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.scan <- xcms.scan%>%
    dplyr::mutate(adrt = adjustedRtime(xcms.xcms),
                  group = xcms.pdata$group[fileIdx],
                  injection_order = xcms.pdata$injection_order[fileIdx])%>%
    dplyr::filter(msLevel==1)

 #col.scale <- c("grey","#38C291",ggsci::pal_aaas()(10))
 #names(col.scale) <- c("Blank","QC",
 #                      setdiff(unique(xcms.scan$group),c("Blank","QC")))

  ggplot(xcms.scan)+
    geom_line(aes(x = retentionTime ,
                  y = adrt-retentionTime,
                  col = injection_order,group=fileIdx))+
    scale_color_gradient(low = "#FFD700",high = "#EE0000")+
    labs(title = "Retention Time adjust",
         x = "Retention Time", y = "Adjusted Error",
         col = "Injection")+
    theme_classic()->p
  p



}


plot_xcms_scan <- function(xcms.xcms){

  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.scan$precursorMZ <- estimatePrecursorIntensity(xcms.xcms,
                                          BPPARAM = BatchtoolsParam(progressbar = T,
                                                                    log = F))


  ggplot(xcms.scan)+
    geom_point(aes(x = retentionTime ,
                   y = precursorMZ ,
                   col = log10(precursorIntensity)))


}


get_xcms_MS_report <- function(xcms.xcms){


  p.tic <- plot_xcms_TIC(xcms.xcms)
  p.rtadj <- plot_xcms_adjustedRT(xcms.xcms )
  p.feature.dis <-plot_xcms_features_distribution(xcms.xcms)


  ### scan
  #plot_xcms_peaks_ms1_scans(xcms.xcms)
  plot_xcms_peaks_ms2_scans(xcms.xcms)



  p1 <- p.tic/p.rtad+(p.feature.dis)

  pdf("d:/temp/aaa.pdf")
  plot(p1)
  dev.off()

  }



get_ms2_feature <- function(precursorMz,
                            rtime ,
                            feature_def
                            ){

  assign_ms2_list <- function(pmz,prt ,feature_def){

    feature_def %>%
      dplyr::mutate(mz.ppm = abs(mzmed-pmz)/pmz,
                    rt.mean = (peakRtMin+peakRtMax)/2,
                    rt.error = abs(prt - rt.mean))%>%
      dplyr::filter( mz.ppm< 1e-5,
                    prt < peakRtMax,
                    prt > peakRtMin)%>%
      dplyr::slice_min(rt.error)%>%
      dplyr::slice_min(mz.ppm)%>%
      dplyr::pull( feature_id)->x
    if (length(x)==0) {
      return(NA)

    }
    return(x)

  }

  ion.df <- data.frame(mz = precursorMz,
                       rt = rtime)%>%
    dplyr::rowwise()%>%
    dplyr::mutate(feature_id = assign_ms2_list(mz,rt,feature_def))%>%
    dplyr::ungroup()


  return(ion.df$feature_id)


}







