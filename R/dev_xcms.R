#' @title  get_features_from_xcms
#' @description extract feature data from xcms::XCMSnExp,
#'  calculate RSD of QC and Sample
#'  ( note this rely on character "QC" and "Sample" in `sampleNames(xcms.xcms)` )
#' @param xcms.xcms XCMSnExp object
#'
#' @return xcms a SummarizedExperiment subject
#' @export
#'

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


get_xcms_feature_se <- function(xcms.xcms,...){

  pol <- c("0" = "neg","1" = "pos")

  xcms.xcms <- xcms_get_feature_stat(xcms.xcms)
  sample.info <- Biobase::pData(xcms.xcms)
  rownames(sample.info) <- sample.info$sample.name

  featuredef <- featureDefinitions(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::mutate(xcms_feature_id = feature_id,
                  feature_id = paste0(feature_id ,"_",pol[as.character(polarity)] ))
  rownames(featuredef) <- featuredef$feature_id

  featureval <- featureValues(xcms.xcms,...)
  colnames(featureval) <- Biobase::pData(xcms.xcms)$sample.name
  rownames(featureval) <- featuredef$feature_id

  feature.se <- SummarizedExperiment(assays = featureval,
                       rowData = featuredef,
                       colData =sample.info
                      )
  return(feature.se)


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
#' @param xcms.xcms XCMSnExp object
#' @param peaks.id char of num
#' @param all.sample should all samples included
#' @param rt one of c("all","identity","expand")
#'
#' @return xcms XChromatograms
#' @import xcms
#' @export
#'

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
    BPPARAM = SerialParam( progressbar = F))

    x.chrom <- XChromatograms(unlist(x.chrom.split))

    x.chrom
  }



  return(x.chrom)
}



xcms_get_peak_fill <- function(xcms.xcms){

  xcms.peaks <- chromPeaks(xcms.xcms)
  rt.na <- apply(xcms.peaks,1,function(x){is.na(x["rt"])})
  into.na <- apply(xcms.peaks,1,function(x){is.na(x["into"])})

  ### fill rt with mean
  xcms.peaks[rt.na,"rt"] <- apply(xcms.peaks[rt.na,c("rtmax","rtmin"),drop =F],
                                  1,mean)
  ### fill into with coef maxo peak with
  xcms.peaks.stat <- get_xcms_peaks_stat(xcms.xcms)%>%
    dplyr::mutate(coef = into/maxo/peakWidth)
  into.coef <- median(xcms.peaks.stat$coef,na.rm = T)
  xcms.peaks[into.na,"into"] <- apply(xcms.peaks[into.na,,drop =F],
                                  1,function(x){x["maxo"]*(x["rtmax"]-x["rtmin"])*into.coef})


  ###return
  chromPeaks(xcms.xcms) <- xcms.peaks
  return(xcms.xcms)
}


#' xcms feature group
#'
#' @param xcms.xcms XCMSnExp object
#'
#' @return xcms
#' @export
#' @import MsFeatures
#'

xcms_get_feature_group <- function(xcms.xcms,
                                   diffRt = 5,
                                   intCor = 0.5,
                                   eicCor = 0.5){

  featureGroups(xcms.xcms) <- NA
  register(SnowParam(progressbar = T))
  if (!is.null(diffRt)) {
    message(Sys.time()," group by SimilarRtimeParam")
    xcms.xcms <- groupFeatures(xcms.xcms,
                                   param = SimilarRtimeParam(diffRt,
                                                             groupFun = groupHclust ))
    message(length(unique(featureGroups(xcms.xcms)))," feature group")
  }
  if (!is.null(intCor)) {
    message(Sys.time()," group by AbundanceSimilarityParam")
    xcms.xcms <- groupFeatures(xcms.xcms,
                                    param = AbundanceSimilarityParam(threshold = intCor,
                                                                     transform = log2 ),
                                    filled = TRUE)
    message(length(unique(featureGroups(xcms.xcms)))," feature group")
  }
  if (!is.null(eicCor)) {
    register(SerialParam())
    message(Sys.time()," group by EicSimilarityParam")
    xcms.xcms <- groupFeatures(xcms.xcms,
                                    param = EicSimilarityParam(threshold = eicCor,
                                                               n=2))
    message(length(unique(featureGroups(xcms.xcms)))," feature group")
  }

  return(xcms.xcms)
}




#' extract_chrom
#'
#' @param xcms.xcms XCMSnExp object
#' @param rtr rt range
#' @param mzr MZ range
#' @param sample
#'
#' @return xcms
#' @export
#'

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
  if ("maxo"%in% sample  )  {
    xcms.sub <- MSnbase::filterFile(xcms.xcms,
                                    which.max(features.val))
  }else if ("all"%in% sample  )  {
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
#' @param xchrom XChromatograms
#'
#' @return xcms
#' @export
#'

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
#' @return xcms
#' @export
#'

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
#' @param xchroms XChromatograms
#'
#' @return xcms
#' @export
#'

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
#' @param xchrom XChromatograms
#' @param norm norm to 0-1
#' @param move move step
#'
#' @return xcms
#' @export
#'

plot_XChromatograms <- function(xchroms ,
                                norm = T,
                                move = T,
                                color_by = c("column","row"),
                                color_f = NULL,
                                label_df = NULL){

  color_by = match.arg(color_by)

  if (norm) {
    xchroms <- normalise(xchroms)
    chrom.data <- get_chroms_data(xchroms)%>%
      dplyr::mutate(intensity = intensity*100)
  }else{
    chrom.data <- get_chroms_data(xchroms)
  }


  chrom.data <- chrom.data%>%
    dplyr::mutate(group_idx = case_when(color_by == "column"~col,
                                      color_by == "row"~row,
                                      T~col))
  if (is.null(color_f)) {
    color_f <-paste0("Peaks_",num2str( unique(chrom.data$group_idx)))
  }
  chrom.data <- chrom.data%>%
    dplyr::mutate(peaks.origin = paste0("peak_",num2str(row),"_",num2str(col)),
                  peaks.origin = factor(peaks.origin,level = unique(peaks.origin)))%>%
    dplyr::group_by(peaks.origin)%>%
    dplyr::mutate(peaks.idx =cur_group_id(),
                  color = color_f[group_idx]
    )%>%
    dplyr::ungroup()



    if (move) {
      chrom.data <- chrom.data%>%
        dplyr::mutate(rt = rt +peaks.idx*3,
                      intensity = intensity+peaks.idx*3)
    }


  ggplot(chrom.data)+
    geom_line(aes(x = rt , y = intensity ,
                  group = peaks.idx,
                  col = color),
              linewidth = 0.5,alpha = 0.8)+
    scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(chrom.data$peaks.idx))))+
    labs(x = "Retention time", y = "Intensity",col = "peaks")+
    theme_bw()->p

  if (!is.null(label_df)) {
    p <- p+ggrepel::geom_text_repel(data = label_df,
                                    aes(x= x,y=y,label = label),
                                    hjust = 0)

  }

  return(p)



}


#'  stat featureDefinitions based on chrompeaks
#' @description extract features' median rt, sn and maxo,
#' `xcms::featureDefinitions()` return a `DataFrame`, in which rtmin, rtmax, rtmed was median of `xcms::chromPeaks()$rt`,
#' but not the median range of peaks. peakRtMin, peakRtMax, peakSN, peakMaxo are median of all peaks in a feature
#'
#' @param xcms.xcms XCMSnExp object
#'
#' @return xcms
#' @export
#'
xcms_get_feature_def_stat <- function(xcms.xcms){

  feature.def <- featureDefinitions(xcms.xcms)
  peaks.data <- chromPeaks(xcms.xcms)

  .xcmsPeakDataMed <- function(x,peaks.data,key = "rtmax",fun = "median"){
    if (!key%in% colnames(peaks.data)) {
      return(NA)

    }
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
  featureDefinitions(xcms.xcms)<-feature.def.df%>%S4Vectors::DataFrame()
  return(xcms.xcms)

}



xcms_get_feature_val_stat <- function(xcms.xcms) {

  xcms.pdata <- Biobase::pData(xcms.xcms)
  featureval <- featureValues(xcms.xcms)
  if ("sample.type" %in% colnames(xcms.pdata)) {
    qc.rsd <- featureval[,xcms.pdata%>%dplyr::filter(sample.type =="QC")%>%
                 dplyr::pull(sampleNames),drop=F]%>%
      apply(1,function(x){
        sd(x,na.rm =T)/mean(x,na.rm=T)
      })
    sample.rsd <- featureval[,xcms.pdata%>%
                               dplyr::filter(sample.type =="Sample")%>%
                 dplyr::pull(sampleNames),drop=F]%>%
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


xcms_get_feature_stat <- function(xcms.xcms){
  xcms.xcms <- xcms.xcms %>%
    xcms_get_feature_def_stat()%>%
    xcms_get_feature_val_stat()
  return(xcms.xcms)
}


xcms_get_feature_isotopologues <- function(xcms.xcms,
                                           isotope = "[13]C",
                                           max_label = 10,
                                           ppm = 5,
                                           net.degree.ratio= 0.5){

  ### calc mz diff within group
  {

    xcms.fdf <- featureDefinitions(xcms.xcms)%>%
      as.data.frame()%>%
      dplyr::arrange(feature_group)%>%
      dplyr::mutate(temp.id = 1:n())%>%
      dplyr::group_by(feature_group)%>%
      dplyr::mutate(group.feature.count = n(),
                    temp2.id = 1:n())

    from.id = rep(xcms.fdf$temp.id,times = xcms.fdf$group.feature.count)
    to.id = xcms.fdf%>%
      dplyr::slice(rep(temp2.id,unique(group.feature.count)))

    fdf.connect <- data.frame(from = from.id,
                              to = to.id$temp.id,
                              feature.group = to.id$feature_group)%>%
      dplyr::filter(from != to)%>%
      dplyr::mutate(from.fid = xcms.fdf$feature_id[from],
                    from.rt = xcms.fdf$rtmed[from],
                    from.mz = xcms.fdf$mzmed[from],
                    to.fid = xcms.fdf$feature_id[to],
                    to.rt = xcms.fdf$rtmed[to],
                    to.mz = xcms.fdf$mzmed[to],
                    mz.diff = to.mz - from.mz)

    #isotope <- "[13]C"
    iso.chemform <- paste0(isotope,1,str_extract(string = isotope,pattern = "[[:alpha:]]+"),-1)
    iso.count <- -max_label:max_label
    iso.mz <- chemform_mz(iso.chemform,0)*iso.count

    fdf.connect <- fdf.connect%>%
      rowwise()%>%
      dplyr::mutate(closest.iso.count = iso.count[which.min(abs(iso.mz-mz.diff))])%>%
      dplyr::mutate(closest.iso.mz = iso.mz[match(closest.iso.count,iso.count)],
                    mz.error = abs(mz.diff-closest.iso.mz),
                    is.iso = mz.error/mean(from.mz,to.mz) < ppm*1e-6)%>%
      dplyr::ungroup()

    fdf.iso.connect <- fdf.connect%>%
      dplyr::filter(is.iso,closest.iso.count != 0)%>%
      dplyr::group_by(from,closest.iso.count)%>%
      dplyr::slice_min(mz.error)%>%
      dplyr::ungroup()%>%
      dplyr::group_by(to,closest.iso.count)%>%
      dplyr::slice_min(mz.error)%>%
      dplyr::ungroup()%>%
      dplyr::group_by(from,closest.iso.count)%>%
      dplyr::slice_min(mz.error)%>%
      dplyr::ungroup()%>%
      dplyr::mutate(from = from.fid,to = to.fid)

  }



  ### assign isotope
  {

    iso.colname <- paste0(str_extract(string = isotope,pattern = "[[:alpha:]]+"),
                          str_extract(string = isotope,pattern = "[[:digit:]]+"))
    xcms.fdf[,paste0(iso.colname,"_seed")] <- NA
    xcms.fdf[,paste0(iso.colname,"_count")] <- NA
    fdf.iso.igraph <- igraph::graph_from_data_frame(fdf.iso.connect)
    node.group <- igraph::components(fdf.iso.igraph)$membership
    xcms.fdf <- as.data.frame(xcms.fdf)
    rownames(xcms.fdf) <-xcms.fdf$feature_id
    for (i in seq_along(unique(node.group))) {


      this.nodes <- names(which(node.group==i))
      this.iso <- fdf.iso.connect %>%
        dplyr::filter(from%in%this.nodes | to %in% this.nodes)

      this.igraph <- igraph::graph_from_data_frame(this.iso)
      #visNetwork::visIgraph(this.igraph)
      to.delete <- degree(this.igraph)<(length(this.nodes)-1)*2*net.degree.ratio
      #message(i," ",sum(to.delete)," of ",length(this.nodes)," nodes remove")
      this.igraph.sub <- delete.vertices(this.igraph,to.delete )
      #visNetwork::visIgraph(this.igraph.sub)
      this.dis <- igraph::distances(this.igraph.sub,mode = "out",
                                    weights = edge.attributes(this.igraph.sub)$closest.iso.count )
      this.dis[this.dis<0] <- 0
      dis.sum <- apply(this.dis,1,sum)
      seed.fid <- names(which.max(dis.sum))
      dis.to.seed <- this.dis[names(which.max(dis.sum)),]
      xcms.fdf[names(dis.to.seed),
               paste0(iso.colname,"_seed")] <- seed.fid
      xcms.fdf[names(dis.to.seed),
               paste0(iso.colname,"_count")] <- unname(dis.to.seed)
      #message(sum(is.na(xcms.fdf$feature_id)))

    }





  }


  ### save to featuredef
  {

    xcms.fdf.temp <- featureDefinitions(xcms.xcms)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    xcms.fdf.temp[,paste0(iso.colname,"_seed")] <- xcms.fdf[rownames(xcms.fdf.temp),paste0(iso.colname,"_seed")]
    xcms.fdf.temp[,paste0(iso.colname,"_count")] <- xcms.fdf[rownames(xcms.fdf.temp),paste0(iso.colname,"_count")]
    xcms.fdf.temp -> featureDefinitions(xcms.xcms)
    message("Get ",
            sum(!is.na(xcms.fdf.temp[,paste0(iso.colname,"_count")])),
            " isotopologues")

  }

  return(xcms.xcms)

}

xcms_get_feature_isotope_label <- function(xcms.xcms,
                                           isotope = "[13]C",
                                           ...){


  ### feature data
  {
    xcms.se <- get_xcms_feature_se(xcms.xcms,
                                   missing = 1)
    rownames(xcms.se) <- rowData(xcms.se)$xcms_feature_id
    iso.colname <- paste0(str_extract(string = isotope,pattern = "[[:alpha:]]+"),
                          str_extract(string = isotope,pattern = "[[:digit:]]+"))
    xcms.rda <- rowData(xcms.se)%>%
      as.data.frame()
    xcms.rda$iso_seed <- xcms.rda[,paste0(iso.colname,"_seed")]
    xcms.rda$iso_count <- xcms.rda[,paste0(iso.colname,"_count")]
    xcms.val <- assay(xcms.se)
  }

  ###calc iso ratio to seed
  {
    xcms.ratio.to.seed <- xcms.val
    xcms.ratio.to.seed[,] <-NA
    xcms.fseed <- xcms.rda$iso_seed %>%
      unique()%>%na.omit()
    for (i in seq_unique(xcms.fseed)) {
      this.fid <- xcms.fseed[i]
      this.iso <- xcms.rda%>%
        dplyr::filter(iso_seed %in% this.fid)
      this.matrix <- xcms.val[this.iso$xcms_feature_id,]
      this.matrix <- t(t(this.matrix)/this.matrix[this.fid,])
      xcms.ratio.to.seed[rownames(this.matrix),] <- this.matrix
    }

    xcms.ratio.to.seed[is.nan(xcms.ratio.to.seed)] <- 0

  }

  ### compare labeled and unlabeled
  {
    is.iso <- xcms.se$isotope_label%in% isotope
    iso.stat <- apply(xcms.ratio.to.seed,1, function(x){
      labeled.mean <- mean(x[is.iso])
      unlabeled.mean <- mean(x[!is.iso])
      p.t.test <- t.test_dev(x[is.iso],x[!is.iso])
      data.frame(labeled.mean,
                 unlabeled.mean,
                 p.t.test
      )
    })%>%
      data.table::rbindlist(idcol = "feature_id")
    iso.stat <- iso.stat%>%
      dplyr::mutate(is_labeled = (labeled.mean > unlabeled.mean&
                                    p.t.test < 0.05))


  }


  ### import to xcms
  {
    xcms.fda <- featureDefinitions(xcms.xcms)
    xcms.fda$is_labeled <- iso.stat$is_labeled
    xcms.fda -> featureDefinitions(xcms.xcms)
    message("Get ",
            sum(xcms.fda$is_labeled ),
            " isotope label")
  }

  return(xcms.xcms)
}


#' match feature to database based on mz and rt
#' MSDB_id in db with mz error < mz.ppm will b
#'
#' @param xcms.xcms XCMSnExp object
#' @param cpdb compoundDb
#' @param mz.ppm num
#' @param rt.tol num
#'
#' @return xcms
#' @export
#'

xcms_get_feature_ms1_candidate <- function(xcms.xcms ,
                                           cpdb,
                                           mz.ppm= 10,
                                           rt.tol = Inf,
                                           selected_adduct = MSCC::adduct.table$Adduct,
                                           ...){


  ### calc adduct and filter range
  cpdbt <- compounds(cpdb, columns = CompoundDb::compoundVariables(cpdb,includeId =T))
  cpdbt <- cpdbt[cpdbt$has_sp>0,]
  cpdbt$formula <- MSCC::chemform_formate(cpdbt$formula)

  adducts <- chemform_adduct_check(selected_adduct)%>%
    dplyr::mutate(polarity = case_when(Ion_mode == "negative"~0,T~1))%>%
    dplyr::filter(polarity %in% polarity(xcms.xcms))
  cp.adduct <- MSCC::chemform_adduct(cpdbt$formula,
                                     adducts$adduct.formated )
  cp.adduct <- cp.adduct%>%
    dplyr::mutate(compound_id=cpdbt$compound_id[id]  )%>%
    dplyr::filter( findInterval(chemform.adduct.mz,
                                mzrange(xcms.xcms))==1)

  ### match database
  xcms.featuredef <- featureDefinitions(xcms.xcms)%>%
    as.data.frame()

  matched.df <- match_mz_rt(mz1 = xcms.featuredef$mzmed,
                            mz2 = cp.adduct$chemform.adduct.mz,
                            mz.ppm = mz.ppm)
  xcms.featuredef$candidate <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$compound_id[as.numeric(idx)]
  })
  xcms.featuredef$candidate.adduct <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$adduct[as.numeric(idx)]
  })
  xcms.featuredef$candidate.mz <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform.adduct.mz[as.numeric(idx)]
  })

  featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(xcms.featuredef)

  return(xcms.xcms)

}


xcms_get_feature_ms2_score <- function(xcms.xcms ,
                                       cpdb,
                                       sp.ms2,
                                       ...){


  ### load spectra database
  Spectra_database <- Spectra(cpdb)
  Spectra_database <- get_Spectra_MEM_backend(Spectra_database)
  Spectra_database <- filterPolarity(Spectra_database,
                                     unique(polarity(xcms.xcms)))
  #spectraNames(Spectra_database) <- Spectra_database$compound_id
  xcms.fdf <- featureDefinitions(xcms.xcms)

  ### sp process
  Spectra_database <- Spectra_database%>%
    filterSpectra_below_PrecursorMz()%>%
    normalizeSpectra(norm_to = "max")%>%
    filterSpectraIntensity(ratio = 0.05)%>%
      applyProcessing()

  if(length(sp.ms2)!=0){
    sp.ms2 <- sp.ms2%>%
      filterSpectra_below_PrecursorMz()%>%
      normalizeSpectra(norm_to = "max")%>%
      filterSpectraIntensity(ratio = 0.05)%>%
      get_Spectra_MEM_backend()%>%
      applyProcessing()

  }



  sp.exp <- sapply(1:nrow(xcms.fdf),function(i){

    x <- xcms.fdf$ms2_id[[i]]
    if (length(x)==0) {
      sp <- makeSpectra(xcms.fdf$mzmed[i],
                        xcms.fdf$rtmed[i])
    }else
    sp <- list(sp.ms2[x])
    return(sp)
  })
  if (!all(unlist(xcms.fdf$candidate)%in% Spectra_database$compound_id)) {
    sp.empty <- makeEmptySpectra(compound_id= setdiff(unlist(xcms.fdf$candidate),
                                                      Spectra_database$compound_id))
    Spectra_database <- c(Spectra_database,sp.empty)
  }
  sp.ref.list <- split(Spectra_database,
                       Spectra_database$compound_id)
  sp.ref <- bplapply(1:nrow(xcms.fdf), function(i){

    cp_id <- xcms.fdf$candidate[[i]]
    sp.temp <- sp.ref.list[cp_id]
    if (length(sp.temp)==0) return(NULL)
    for (j in 1:length(cp_id)) {
      sp.temp[[j]]$adduct <- xcms.fdf$candidate.adduct[[i]][j]
      sp.temp[[j]]$precursorMz <- xcms.fdf$candidate.mz[[i]][j]
    }
    sp.temp <- do.call(what = "c",args = unname(sp.temp))
    return(sp.temp)
  },BPPARAM = SerialParam(progressbar = T))



  ### output all candidate score
  {
    .f <- function(expSpec,refSpec,...){
      if (is.null(refSpec)) {
        return(NULL)
      }
      scorem <- compareSpectra(expSpec,refSpec,...)
      dim(scorem) <- c(length(expSpec),length(refSpec))
      scorem[is.infinite(scorem)|is.na(scorem )] <- 0
      scores <- apply(scorem,2,max,na.rm=T)
      unname(mean_f(scores,paste0(refSpec$compound_id,"_",refSpec$adduct)))
    }
    xcms.fdf$candidate.score <- BiocParallel::bplapply(1:length(sp.exp),
                                                       function(i){
      .f(expSpec = sp.exp[[i]], refSpec = sp.ref[[i]])
    },BPPARAM = BiocParallel::SerialParam(
      progressbar = T))


  }


  featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(xcms.fdf)

  return(xcms.xcms)

}


get_xcms_feature_all_candidate <- function(xcms.xcms){

  xcms.fdf <- get_xcms_feature_definitions(xcms.xcms)
  candi.rda <- xcms.fdf%>%
    dplyr::mutate(candidate.n = lengths(candidate))
  candi.rda.split <- candi.rda[rep(candi.rda$feature_id,
                                   candi.rda$candidate.n),]%>%
    dplyr::group_by(feature_id)%>%
    dplyr::mutate(temp_id = 1:n())%>%
    dplyr::rowwise()%>%
    dplyr::mutate(MSDB_id = candidate[[temp_id]],
                  adduct = candidate.adduct[[temp_id]],
                  mz_ref = candidate.mz[[temp_id]],
                  score = candidate.score[[temp_id]])%>%
    dplyr::ungroup()%>%
    dplyr::select(-c(candidate,candidate.adduct,candidate.mz,candidate.score))

  return(candi.rda.split)

}

xcms_get_feature_annotation <- function(xcms.xcms,
                                        ...){


  xcms.fdf <- featureDefinitions(xcms.xcms)
  xcms.fdf$compound_id <- NA
  xcms.fdf$adduct <- NA
  xcms.fdf$score <- NA
  xcms.fdf$mz_ref <- NA
  xcms.fdf$rt_ref <- NA
  for (i in 1:nrow(xcms.fdf)) {

    this.mz <- xcms.fdf$mzmed[i]
    candi.compound_id <- xcms.fdf$candidate[[i]]
    candi.mz <- xcms.fdf$candidate.mz[[i]]
    candi.adduct <- xcms.fdf$candidate.adduct[[i]]
    candi.score <- xcms.fdf$candidate.score[[i]]

    if (length(candi.compound_id)==0) next

    ### score
    score.mz <- abs(candi.mz-this.mz)
    score.mz <- 1-score.mz/max(score.mz)
    score.ms2 <- candi.score
    score <- score.mz * 0.2 + score.ms2*0.8
    selected <- which.max(score)

    ### info
    xcms.fdf$compound_id[i] <- candi.compound_id[selected]
    xcms.fdf$adduct[i] <- candi.adduct[selected]
    xcms.fdf$score[i] <- candi.score[selected]
    xcms.fdf$mz_ref[i] <- candi.mz[selected]
    #xcms.fdf$rt_ref[i] <- candi.msdbid[selected]

  }



  featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(xcms.fdf)

  return(xcms.xcms)



}


get_xcms_feature_definitions <- function(xcms.xcms){
  xcms.fdf <- featureDefinitions(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::select(
      !c(mzmin,mzmax,rtmin,rtmax,npeaks,peakidx)
    )
  return(xcms.fdf)

}
#' @title plot_xcms_peaks_distribution
#' @description export peaks data by xcms::chromPeaks and plot by ggplot2
#'
#' @param xcms.xcms XCMSnExp object
#' @param plot.title title
#' @param type `"o"`, for geom_point, `"l"`, for geom_segment
#'
#' @return xcms
#' @export
#'

plot_xcms_peaks_distribution <- function(xcms.xcms,plot.title = "Peaks distribution",type = "o"){

  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::mutate(as.data.frame(chromPeakData(xcms.xcms)),
                  peak_id = rownames(.),
                  merged = grepl(peak_id,pattern = "CPM"))%>%
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
#' @param xcms.xcms XCMSnExp object
#' @param plot.title title
#'
#' @return xcms
#' @export
#'

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
#' @param xcms.xcms XCMSnExp object
#' @param feature.id feature id
#' @param sampleNames
#'
#' @return xcms
#' @export
#'

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
  xcms.chrom <- extract_chrom(xcms.sub ,
                              mzr = mz.range,
                              rtr = rt.range)

  xcms.chrom.data <- get_chroms_data(xcms.chrom)%>%
    dplyr::mutate(group = sampleNames[col])

  ggplot(xcms.chrom.data)+
    geom_line(aes(x = rt,y = intensity , col = group))+
    xlim(c(min(rtime(xcms.sub)),max(rtime(xcms.sub))))+
    labs(col = "",x = "Retention time", y = "Intensity",
         title = paste0(feature.id),
         subtitle = paste0( "mz: ",round(mz.range[1],4),
                            " ~ ",round(mz.range[2],4),
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
    stat_density_2d(aes(fill= after_stat(level)),
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
#' @param xcms.xcms XCMSnExp object should be a `XCMSnExp` object after `findChromPeaks`
#' @param plot.title title
#'
#' @return xcms
#' @export
#'

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



#' Title
#'
#' @param xcms.xcms XCMSnExp object
#' @param plot.title title
#'
#' @return xcms
#' @export
#' @import xcms
#'

plot_xcms_peaks_ms2_scans <- function(xcms.xcms,plot.title = "Peaks Sans of MS2"){

  xcms.process.type <- processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    processParam()
  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()
  xcms.scans <- Biobase::fData(xcms.xcms)%>%
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

 open_ggplot_win(peaks.dis.plot,width = 25,height = 5)
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

#' @title plot_xcms_peaks_Chromatogram
#' @description extract EIC according to peaks' mzrange and rtrange,
#' note that if multiple sample in xcms object, only first sample will be extracted
#'
#' @param xcms.xcms XCMSnExp object
#' @param peak_ids peaks id
#' @param rt_expand foldchange to expand rt range
#'
#' @return xcms
#' @export
#'

plot_xcms_peaks_Chromatogram <- function(xcms.xcms,peak_id,rt = "expand"){

  peaks.data <- chromPeaks(xcms.xcms)[peak_id,,drop = F]
  peak_id <- rownames(peaks.data)
  mz.range <- c(peaks.data[,c("mzmin","mzmax")])
  rt.range <- c(peaks.data[,c("rtmin","rtmax")])
  xcms.chrom <- get_xcms_peaks_chrom(xcms.xcms,
                                     peaks.id = peak_id,
                                     rt = rt)
  chrom.data <- get_chroms_data(xcms.chrom)%>%
    dplyr::mutate(fill = rt > min(rt.range)&rt <max(rt.range),
                  sample = sampleNames(xcms.xcms)[row]
                 )%>%
    dplyr::filter(!is.na(intensity))

  ggplot(chrom.data)+
    geom_line(aes(x = rt,y = intensity,col =sample),linetype = 1)+
    geom_area(aes(x = rt,y = intensity, fill = sample),
              stat = "identity",alpha = 0.1)+
    #geom_point(aes(x = rt, y = fit))+
    scale_fill_manual(values = c("FALSE" = "transparent","TRUE" = "grey"))+
    labs(title = paste0(peak_id),
         subtitle = paste0("mz:",paste0(sprintf("%.4f",range(mz.range)),collapse = " - "), ";     ",
                           "rt:",paste0(sprintf("%.2f",range(rt.range)),collapse = " - "),"\n",
                           "mz range = ",sprintf("%.2f",mean(diff(range(mz.range))/mz.range)*1e6)," ppm;     ",
                           "peak width = ", sprintf("%.2f",diff(range(rt.range))),"\n"
                          # "shape score = ",get_chrom_peaks_shape_score(xcms.chrom[1,1])
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
#' @return xcms
#' @export
#' @import xcms

xcmsProcessingMS1 <- function(xcms.xcms,
                              ion_mode = NA,
                              xcms_param = NULL,
                              ...){



  if (is.na(ion_mode)) {
    ion_mode <- polarity(xcms.xcms )%>%unique()
    if (length(ion_mode)!=1) {
      stop("MS1 scans contain both positive and negative, please check")
    }
  }

  xcms.xcms <- ProtGenerics::filterPolarity(xcms.xcms , ion_mode)


  ### Find peaks
  message(Sys.time()," Find peaks...")
  xcms.xcms<-xcms::findChromPeaks(xcms.xcms,
                            param = xcms_param$findChromPeaks,
                            BPPARAM  = BiocParallel::SnowParam(progressbar = T))
  xcms.xcms <- xcms_get_peak_fill(xcms.xcms)
  #mpp <- xcms::MergeNeighboringPeaksParam(expandRt = 2.5,minProp = 0.5)
  #xcms.xcms <- xcms::refineChromPeaks(xcms.xcms, mpp,
  #                                    BPPARAM  = BiocParallel::SerialParam(progressbar = T))

  ### adujust RT
  message(Sys.time()," Adjust RT...")
  peaksGroup <- Biobase::pData(xcms.xcms)$sample.type
  peak.density.param <- xcms::PeakDensityParam(sampleGroups = peaksGroup,
                                         minFraction = 0.4,bw = 30,
                                         binSize = 0.015)
  xcms.xcms <- xcms::groupChromPeaks(xcms.xcms,param = peak.density.param)



  if (length(oligoClasses::sampleNames(xcms.xcms))>1) {
    if (sum(peaksGroup=="QC") <2 ) {
      rt.adjust.param <- PeakGroupsParam(minFraction = 0.4,
                                         #subset = which(peaksGroup == "QC"),
                                         subsetAdjust = "previous",span = 0.4)
      xcms.xcms <- adjustRtime(xcms.xcms,param = rt.adjust.param)
    }else{
      ### adjust based on QC
      rt.adjust.param <- PeakGroupsParam(minFraction = 0.4,
                                          subset = which(peaksGroup == "QC"),
                                          subsetAdjust = "average",span = 0.4)
      xcms.xcms <- adjustRtime(xcms.xcms,param = rt.adjust.param)
    }
  }


  ### group peaks
  message(Sys.time()," Group peaks...")
  peak.density.param <- xcms_param$groupChromPeaks
  peak.density.param@sampleGroups <- Biobase::pData(xcms.xcms)$sample.type
  xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)
  message(Sys.time()," ",nrow(featureDefinitions(xcms.xcms))," feature found")
  xcms.xcms <- fillChromPeaks(xcms.xcms,param = FillChromPeaksParam())



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
#' @param xcms.xcms XCMSnExp object
#' @param feature_id_to_show
#'
#' @return xcms
#' @export
#'

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
    dplyr::mutate(scan_id = paste0("scan_",fileStr,"_",spStr))%>%
    dplyr::ungroup()%>%
    dplyr::select(scan_id,ms1_no,ms1_group,ms1_group_rt,
                  ms2_count,cycle_time,scan_time,everything(),
                  -c(x,fileStr,spStr,ms1_no_str))%>%
    as.data.frame()
  suppressMessages(row.names(xcms.fdata) <- rownames(fData(xcms.xcms)))

  xcms.fdata
}


xcms_get_scan_Stat <- function(xcms.xcms){
  xcms.fdata <- get_xcms_scan_Stat(xcms.xcms )
  fData(xcms.xcms) <- xcms.fdata
  return(xcms.xcms)
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
    geom_line(aes(x = retentionTime , y = tic,
                  col = group,
                  group=fileIdx))+
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
 #xcms.scan$precursorMZ <- estimatePrecursorIntensity(xcms.xcms,
 #                                        BPPARAM = BatchtoolsParam(progressbar = T,
 #                                                                  log = F))


  ggplot(xcms.scan)+
    geom_point(aes(x = retentionTime ,
                   y = precursorMZ ,
                   col = log10(precursorIntensity)))


}


get_xcms_MS_report <- function(xcms.xcms ,
                               file.path){


  file.path <- "d:/temp/xcms.report.pdf"
  p.tic <- plot_xcms_TIC(xcms.xcms)
  p.rtadj <- plot_xcms_adjustedRT(xcms.xcms )
  p.feature.dis <- plot_xcms_features_distribution(xcms.xcms)


  ### scan
  #plot_xcms_peaks_ms1_scans(xcms.xcms)
  plot_xcms_peaks_ms2_scans(xcms.xcms)



  p1 <- p.tic/p.rtad+(p.feature.dis)

  pdf("d:/temp/aaa.pdf")
  plot(p1)
  dev.off()

  }



get_xcms_peaks_stat <- function(xcms.xcms){

  xcms.peaks <- chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    rownames_to_column("peak_id")%>%
    dplyr::mutate(peakWidth = rtmax-rtmin,
                  mzWidth = mzmax-mzmin,
                  mzError = mzWidth/mz*1e6)
  return(xcms.peaks)

}



get_xcms_centwave_tune <- function(xcms.xcms,
                                   iteration = 10){

  cwp <- CentWaveParam(peakwidth = c(5,20),
                       verboseColumns=T,fitgauss = T)
  xcms.tune.df <- data.frame(
    No = 1:iteration,
    ppm = cwp@ppm,
    pwmin = cwp@peakwidth[1],
    pwmax = cwp@peakwidth[2],
    snthresh = cwp@snthresh,
    prefilter = cwp@prefilter[1],
    prefilter.int = cwp@prefilter[2],
    peaks.no = NA,
    mze.range = NA,
    pw.range = NA
  )
  xcms.tune.list <- list()
  for (i in 1:iteration) {

    message("run ",i," ", Sys.time())
    show(cwp)

    ### record param
    xcms.tune.df$ppm[i] <- cwp@ppm
    xcms.tune.df$pwmin[i] <- cwp@peakwidth[1]
    xcms.tune.df$pwmax[i] <- cwp@peakwidth[2]
    xcms.tune.df$snthresh[i] <- cwp@snthresh
    xcms.tune.df$prefilter[i] <- cwp@prefilter[1]
    xcms.tune.df$prefilter.int[i] <- cwp@prefilter[2]


    ### run
    xcms.interation <- findChromPeaks(
      xcms.xcms,
      param = cwp,
      BPPARAM  =SerialParam()
    )
    xcms.peaks <- get_xcms_peaks_stat(xcms.interation)%>%
      dplyr::mutate(scan.no = scmax-scmin)
    xcms.peaks.high.sn <- xcms.peaks%>%
      dplyr::slice_max(sn , n = round(nrow(.)*0.1))

    ### result
    pw <- xcms.peaks.high.sn$peakWidth
    mze <- xcms.peaks.high.sn$mzError
    xcms.tune.df$peaks.no[i] <- nrow(xcms.peaks)
    xcms.tune.df$mze.range[i] <- paste0(quantile(mze,c(0.05,0.95))%>%round,collapse = "-")
    xcms.tune.df$pw.range[i] <- paste0(quantile(pw,c(0.25,0.75))%>%round,collapse = "-")

    ### update param
    cwp@ppm <- quantile(mze,0.5)
    cwp@peakwidth <- quantile(pw,c(0.25,0.75))
    #cwp@snthresh <- quantile(xcms.peaks$sn,0.05)
    #cwp@prefilter[1] <- quantile(xcms.peaks$scan.no,0.05)

    ### record
    xcms.tune.list[[i]] <- xcms.interation


  }





}



get_xcms_Spectra <- function(xcms.xcms){

  xcms.files <- paste0(dirname(xcms.xcms),"/",sampleNames(xcms.xcms))
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.sp <- Spectra(xcms.files,
                         backend = MsBackendDataFrame(),
                         BPPARAM = SerialParam(progressbar = T))%>%
    filterPolarity(unique(polarity(xcms.xcms)))
  spectraNames(xcms.sp) <- xcms.sp$scan_id <- xcms.scan$scan_id
  return(xcms.sp)

}


setMethod(f = "filepaths",
          signature = "XCMSnExp",
          definition = function(object){
  paste0(dirname(object),"/",sampleNames(object))
})
setMethod(f = "mzrange",
          signature = "XCMSnExp",
          definition = function(object){
            xcms.fdata <- fData(object)
            return(c(min(xcms.fdata$scanWindowLowerLimit,na.rm = T),
                     max(xcms.fdata$scanWindowUpperLimit,na.rm = T)))
          })



get_xcms_precursor_intensity <- function(xcms.xcms,...){

  estimatePrecursorIntensity(xcms.xcms,
                             method = "previous",...)

}




plotly_feature_span <- function(xcms.fdf){

  plot_ly(xcms.fdf)%>%
    add_segments(x = ~rtmin,xend = ~rtmax,
                 y = ~mzmed , yend = ~mzmed)

}



simulate_dda <- function(xcms.fdf,
                         ms1.time = 0.6,
                         ms2.time = 0.6,
                         topn = 10){

  t <- 0
  scan.df <- data.frame(spIdx = 1,
                        msLevel= 1,
                        retentionTime= 0,
                        precursorMZ= NA,
                        ion_id = NA
  )
  while(t < max(xcms.fdf$rtmax)){

    ion.to.ms2 <- which(t < xcms.fdf$rtmax&t>xcms.fdf$rtmin)
    ion.to.ms2 <- na.omit(ion.to.ms2[order(xcms.fdf$peakMaxo[ion.to.ms2],decreasing = T)[1:topn]])
    ms2.scan <- data.frame(spIdx = rep(NA,length(ion.to.ms2)),
                           msLevel= rep(2,length(ion.to.ms2)),
                           retentionTime= t + seq_along(ion.to.ms2)*ms2.time,
                           precursorMZ = xcms.fdf$mzmed[ion.to.ms2],
                           ion_id = ion.to.ms2
    )
    scan.df <- rbind(scan.df,ms2.scan)
    t <- t + ms2.time*length(ion.to.ms2)

    ms1.scan <- data.frame(spIdx = NA,
                           msLevel= 1,
                           retentionTime = t + ms1.time,
                           precursorMZ= NA,
                           ion_id = NA
    )
    scan.df <- rbind(scan.df,ms1.scan)
    t <-  t + ms1.time
  }


  return(scan.df)

}
