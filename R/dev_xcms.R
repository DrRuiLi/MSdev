
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

  xcms.sum <- xcms::quantify(xcms.xcms,missing = missing )
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


#' @description Xcms feature se.
#' @describeIn xcms_extenstion extract feature data from xcms, convert to SummarizedExperiment
#' @title get_xcms_feature_se
#' @param xcms.xcms xcms object
#' @param missing how missing values should be reported. Allowed values are NA (the default), a numeric or missing = "rowmin_half". The latter replaces any NA with half of the row's minimal (non-missing) value.
#'
#' @returns SummarizedExperiment
#' @export
#'
get_xcms_feature_se <- function(xcms.xcms,...){

  pol <- c("0" = "neg","1" = "pos")

  xcms.xcms <- xcms_get_feature_stat(xcms.xcms)
  sample.info <- Biobase::pData(xcms.xcms)
  rownames(sample.info) <- sample.info$sample.name

  featuredef <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::mutate(xcms_feature_id = feature_id,
                  feature_id = paste0(feature_id ,"_",pol[as.character(polarity)] ))
  rownames(featuredef) <- featuredef$feature_id

  featureval <- xcms::featureValues(xcms.xcms,...)
  colnames(featureval) <- Biobase::pData(xcms.xcms)$sample.name
  rownames(featureval) <- featuredef$feature_id

  feature.se <- SummarizedExperiment::SummarizedExperiment(assays = featureval,
                       rowData = featuredef,
                       colData =sample.info
                      )
  return(feature.se)


}


get_chrom_peaks_shape_score <- function(chrom,
                                        peak.id = chrom@chromPeakData@rownames){
  peak.id = chrom@chromPeakData@rownames
  peak.id <- peak.id[1]
  peaks.data <- xcms::chromPeaks(chrom)[peak.id,,drop = FALSE]
  rtime <- xcms::rtime(chrom)
  int <- xcms::intensity(chrom)

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

  peak.info <- xcms::chromPeaks(xchrom)
  peaks.no <- nrow(peak.info)


  for (i in 1:peaks.no) {


  }


}


sub_xchrom <-  function (x, i, j, drop = TRUE)
{
  if (missing(i) && missing(j))
    return(x)
  if (missing(i))
    i <- seq_len(nrow(x))
  if (missing(j))
    j <- seq_len(ncol(x))
  if (is.logical(i))
    i <- which(i)
  if (is.logical(j))
    j <- which(j)
  if (length(i) > 1 || length(j) > 1)
    drop <- FALSE
  if (length(i) == 1 && length(j) == 1 && drop)
    return(x@.Data[i, j, drop = TRUE][[1]])
  cpeaks_orig <- xcms::chromPeaks(x)
  fts_orig <- xcms::featureDefinitions(x)
  ph <- x@.processHistory
  pd <- x@phenoData
  fd <- x@featureData
  xclass <- class(x)
  x <- as(x@.Data[i = i, j = j, drop = FALSE], xclass)
  pd <- pd[j, ]
  Biobase::pData(pd) <- droplevels(pData(pd))
  x@phenoData <- pd
  fd <- fd[i, ]
  Biobase::pData(fd) <- droplevels(pData(fd))
  x@featureData <- fd
  if (nrow(fts_orig)) {
  cpeaks_sub <- xcms::chromPeaks(x)
    fts <- vector("list", length(i))
    for (el in seq_along(i)) {
      fts_row <- fts_orig[fts_orig$row == i[el], ,
                          drop = FALSE]
      if (nrow(fts_row)) {
        fts_row$row <- el
        fts_row <- xcms:::.subset_features_on_chrom_peaks(fts_row,
                                                          cpeaks_orig, cpeaks_sub)
        fts[[el]] <- fts_row
      }
      else fts[[el]] <- DataFrame()
    }
    x@featureDefinitions <- do.call(rbind, fts)
  }

  x
}

get_xchroms_peaks_count <- function(xchroms){

  peaks.info <- xcms::chromPeaks(xchroms)%>%
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

#' @describeIn xcms_extenstion extract chromatogram for a peak
#' @title get_xcms_peaks_chrom
#' @description
#' extract chromatograph from XCMSnExp,
#' if `all.sample` = F, only the samples, in which given peaks.id are detected will be return,
#' else extract from all samples
#'
#' @param xcms.xcms XCMSnExp object
#' @param peaks.id char or num
#' @param all.sample should all samples be included
#' @param rt one of c("all","identity","expand")
#'
#' @return XChromatograms object
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
    x.chrom <- get_xcms_chromatogram(xcms.xcms,
                             mz = peaks.data[,c("mzmin","mzmax")],
                             rt = peaks.data[,c("rtmin","rtmax")])
  }else{
    peaks.id.split <- split(rownames(peaks.data),
                              f = peaks.data[,"sample"])
    x.chrom.split <- bplapply(peaks.id.split,
                              FUN = function(x,xcms.x,peaks.data){

      x <- peaks.data[x,,drop =F]
      xcms.sub <- MSnbase::filterFile(xcms.x, as.integer(unique(x[,"sample"])))
      get_xcms_chromatogram(xcms.sub,
                    mz = x[,c("mzmin","mzmax")],
                    rt = x[,c("rtmin","rtmax")])
    },xcms.x=xcms.xcms ,peaks.data=peaks.data,
    BPPARAM = SerialParam( progressbar = F))

    x.chrom <- do.call(cbind_Chromatograms, x.chrom.split)

    x.chrom
  }



  return(x.chrom)
}



xcms_get_peak_fill <- function(xcms.xcms){

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)
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
  xcms.ph <- xcms::processHistory(xcms.xcms)
  xcms::chromPeaks(xcms.xcms) <- xcms.peaks
  xcms.xcms@.processHistory <- xcms.ph
  return(xcms.xcms)
}


#' @describeIn xcms_extenstion group features
#' @title Group xcms Features
#' @description Groups features from an XCMSnExp object using multiple criteria: similarity in retention time, abundance (intensity) correlation, and EIC (extracted ion chromatogram) correlation. The grouping is performed sequentially using the MsFeatures package functions.
#' @param xcms.xcms XCMSnExp object containing feature definitions.
#' @param diffRt numeric. Maximum allowed retention time difference for grouping by SimilarRtimeParam. If NULL, retention time grouping is skipped. Default is 5.
#' @param intCor numeric. Threshold for abundance similarity (correlation) grouping using AbundanceSimilarityParam. If NULL, intensity correlation grouping is skipped. Default is 0.5.
#' @param eicCor numeric. Threshold for EIC similarity grouping using EicSimilarityParam. If NULL, EIC correlation grouping is skipped. Default is 0.5.
#'
#' @return XCMSnExp object with featureGroups column added or updated.
#' @export
#'

xcms_get_feature_group <- function(xcms.xcms,
                                   diffRt = 5,
                                   intCor = 0.5,
                                   eicCor = 0.5){

  xcms::featureGroups(xcms.xcms) <- NA
  register(SnowParam(progressbar = T))
  if (!is.null(diffRt)) {
    message_with_time(" group by SimilarRtimeParam")
    xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                   param = MsFeatures::SimilarRtimeParam(diffRt,
                                                             groupFun = groupHclust ))
    message(length(unique(xcms::featureGroups(xcms.xcms)))," feature group")
  }
  if (!is.null(intCor)) {
    message_with_time(" group by AbundanceSimilarityParam")
    xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                    param = MsFeatures::AbundanceSimilarityParam(threshold = intCor,
                                                                     transform = log2 ),
                                    filled = TRUE)
    message(length(unique(featureGroups(xcms.xcms)))," feature group")
  }
  if (!is.null(eicCor)) {
    register(SerialParam())
    message_with_time(" group by EicSimilarityParam")
    xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                    param = MsFeatures::EicSimilarityParam(threshold = eicCor,
                                                               n=2))
    message(length(unique(xcms::featureGroups(xcms.xcms)))," feature group")
  }

  return(xcms.xcms)
}




#' @describeIn xcms_extenstion extract chromatograms
#' @title Get Xcms Chromatogram
#' @description Extracts chromatograms from each file in the XCMSnExp object using xcms::chromatogram and combines them into a single XChromatograms object. This function iterates over each file, extracts chromatograms with provided parameters, and returns a combined chromatograms object.
#' @param object XCMSnExp object from which to extract chromatograms.
#' @param BPPARAM BiocParallel backend for parallel processing. Default is SerialParam().
#' @param ... Additional arguments passed to xcms::chromatogram, such as rt (retention time range) and mz (m/z range).
#'
#' @return XChromatograms object containing extracted chromatograms for all files.
#' @export
#'
get_xcms_chromatogram <- function(object,
                                  BPPARAM= SerialParam(),
                                  ...){


  xcms.split <- lapply(seq_along(MSnbase::fileNames(object)),
                       function(x) MSnbase::filterFile(object,x))

  xcms.chrom <- bplapply(xcms.split,
                         function(x,...){
    #message(nrow(rt))
    register(SerialParam())
    y <-NA
    try(y <- xcms::chromatogram(x,BPPARAM = BPPARAM,...)
    )
    return(y)
  },..., BPPARAM = BPPARAM)


  do.call(cbind_Chromatograms,xcms.chrom)

}


get_xcms_feature_chrom <- function(xcms.xcms,
                                    feature.id = xcms::featureDefinitions(xcms.xcms)$feature_id,
                                    sample = "maxo",
                                    rt = "expand",
                                    mz.expand = 0,
                                    BPPARAM = SerialParam()){

  features.data <- xcms::featureDefinitions(xcms.xcms)
  features.val <- xcms::featureValues(xcms.xcms,missing = "rowmin_half")
  if(is.numeric(feature.id)) { feature.id <-rownames(features.data)[feature.id]}
  features.data <- features.data[feature.id,,drop=F]
  features.val <- features.val[feature.id,,drop =F]
  if ("maxo"%in% sample  )  {
    xcms.sub <- MSnbase::filterFile(xcms.xcms,
                                    which.max(colMeans(features.val)))
  }else if ("all"%in% sample  )  {
    xcms.sub <- xcms.xcms
  }else {
    xcms.sub <- MSnbase::filterFile(xcms.xcms, sample )
  }



  rtr <- switch (rt,
                 "all" = c(min(rtime(xcms.sub)),max(rtime(xcms.sub))),
                 "expand" = t(apply(features.data[,c("peakRtMin","peakRtMax"),drop =F],1,expand_range,add = 15)),
                 "identity" = features.data[,c("peakRtMin","peakRtMax"),drop =F]

  )
  mzr <- features.data[,c("peakMzMin","peakMzMax")]%>%as.matrix()
  if (mz.expand > 0) {
    mz_range <- mzr[,2] - mzr[,1]
    mzr[,1] <- mzr[,1] - mz_range * mz.expand
    mzr[,2] <- mzr[,2] + mz_range * mz.expand
  }

  x.chrom <-  get_xcms_chromatogram(xcms.sub,
                                  mz = mzr,
                                  rt = rtr,
                                  BPPARAM = BPPARAM)
  return(x.chrom)
}


get_chrom_peaks_gaussian_fit <- function(xchrom){

  peaks.info <- xcms::chromPeaks(xchrom)[1,,drop = F]
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


#' @describeIn xcms_extenstion convert retention time units
#' @title XChromatograms Rt Unit
#' @description Changes the retention time units of XChromatograms objects. In some situations (e.g., SRM data from Thermo), retention times are recorded in minutes, which can cause errors during peak detection. This function converts between seconds and minutes.
#' @param xchroms XChromatograms or MChromatograms object.
#' @param unit_to Target unit: "s" (seconds) multiplies by 60, "m" (minutes) divides by 60. Default is "s".
#' @param BPPARAM BiocParallel backend for parallel processing. Default is BatchtoolsParam.
#'
#' @return XChromatograms object with converted retention times.
#' @export
#'
XChromatograms_rt_unit <- function(xchroms,unit_to = "s",
                                   BPPARAM = BatchtoolsParam(progressbar = T,log = F,
                                                             registryargs = batchtoolsRegistryargs(packages = c("MSnbase")))){


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

  pda <- pData(xchroms)
  fda <- fData(xchroms)

  xchroms.trans <- bplapply(1:length(xchroms),FUN = f,
                              xchroms = xchroms,
                              unit.mulit = unit.mulit,
                              BPPARAM = BPPARAM)
  xchroms.trans <- XChromatograms(xchroms.trans,
                                    nrow = dim(xchroms)[1],byrow = T)

  message( "rtime value ", round(rtime.max,0), " change to ", round(max(rtime(xchroms.trans[1,1])),0))

  pda -> pData(xchroms.trans)
  fda -> fData(xchroms.trans)

  return(xchroms.trans)


}



#' @describeIn xcms_extenstion fill chromatograms with fewer than two data points
#' @title XChromatograms Fill 2point
#' @description When using xcms::findChromPeaks, chromatograms with fewer than two data points cause errors. This function identifies such chromatograms and adds a duplicate point (time +1, intensity 0) to ensure at least two points exist.
#' @param xchroms XChromatograms object to be checked and filled.
#'
#' @return XChromatograms object with chromatograms having at least two data points.
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


#' @describeIn xcms_extenstion plot chromatograms
#' @title Plot XChromatograms
#' @description Plots XChromatograms data as line plots, with options to normalize intensities to 0-1 range, offset chromatograms for clarity, and customize colors. Returns a ggplot object.
#' @param xchroms XChromatograms object to plot.
#' @param norm logical. If TRUE, normalize intensities to 0-1 range (default TRUE).
#' @param move logical. If TRUE, offset chromatograms by index for better visibility (default TRUE).
#' @param color_by Character indicating grouping for coloring: "column" (by sample) or "row" (by feature). Default is "column".
#' @param color_f Optional character vector of colors for groups. If NULL, uses distinctColorPalette.
#' @param label_df Optional data frame with columns x, y, label for adding text labels via ggrepel.
#'
#' @return ggplot object.
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
    xchroms <- MSnbase::normalise(xchroms)
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
    #scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(chrom.data$peaks.idx))))+
    ggsci::scale_fill_npg()+
    labs(x = "Retention time", y = "Intensity",col = "peaks")+
    theme_bw()->p

  if (!is.null(label_df)) {
    p <- p+ggrepel::geom_text_repel(data = label_df,
                                    aes(x= x,y=y,label = label),
                                    hjust = 0)

  }

  return(p)



}


#' @describeIn xcms_extenstion calculate feature statistics
#' @title Xcms Get Feature Def Stat
#' @description Extracts and adds median retention time, signal-to-noise ratio, and maximum intensity for each feature. While xcms::featureDefinitions() provides median mz and rt, this function calculates median values across all peaks within a feature: peakRtMin, peakRtMax, peakWidth, peakMzMin, peakMzMax, peakSN, peakMaxo, and polarity.
#' @param xcms.xcms XCMSnExp object with feature definitions and chromPeaks.
#'
#' @return XCMSnExp object with updated featureDefinitions containing additional statistics.
#' @export
#'
xcms_get_feature_def_stat <- function(xcms.xcms){

  feature.def <- xcms::featureDefinitions(xcms.xcms)
  peaks.data <- xcms::chromPeaks(xcms.xcms)

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
  xcms::featureDefinitions(xcms.xcms) <- feature.def.df %>% S4Vectors::DataFrame()
  return(xcms.xcms)

}



xcms_get_feature_val_stat <- function(xcms.xcms) {

  xcms.pdata <- Biobase::pData(xcms.xcms)
  featureval <- xcms::featureValues(xcms.xcms)
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

  fdf <- xcms::featureDefinitions(xcms.xcms)
  fdf$qc_rsd <- qc.rsd
  fdf$sample_rsd <- sample.rsd
  fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)
}


xcms_get_feature_stat <- function(xcms.xcms){
  xcms.xcms <- xcms.xcms %>%
    xcms_get_feature_def_stat()%>%
    xcms_get_feature_val_stat()
  return(xcms.xcms)
}


#' @describeIn xcms_extenstion identify isotopologues
#' @title Xcms Get Feature Isotopologues
#' @description Screens isotopologue peaks based on m/z and retention time differences, assigns isotopologue groups and seeds, and records results in featureDefinitions. Uses graph-based clustering to identify isotopologue networks.
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param iso_ele Isotope element string (e.g., `"[13]C"`) for mass difference calculation.
#' @param max_label Maximum number of isotope labels to consider (default 10).
#' @param ppm Mass accuracy tolerance in ppm (default 10).
#' @param rt.tol Retention time tolerance in seconds for grouping (default 5).
#' @param net.degree.ratio Ratio threshold for network degree to assign isotopologue seeds (default 0.3).
#'
#' @return XCMSnExp object with featureDefinitions updated with iso_seed, iso_count, and iso_connection_group columns.
#' @export
#'
xcms_get_feature_isotopologues <- function(xcms.xcms,
                                           iso_ele = "[13]C",
                                           max_label = 10,
                                           ppm = 10,
                                           rt.tol = 5,
                                           net.degree.ratio = 0.3){



  fdf.iso.connect <- get_xcms_feature_iso_connection(xcms.xcms,iso_ele,max_label,
                                                     ppm,rt.tol )


  ### assign isotopologues
  {

    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
      as.data.frame()
    xcms.fdf[,paste0("iso_seed")] <- NA
    xcms.fdf[,paste0("iso_count")] <- NA
    fdf.iso.igraph <- igraph::graph_from_data_frame(fdf.iso.connect)
    #fdf.iso.igraph <- igraph_filter_vertex(fdf.iso.igraph,degree(fdf.iso.igraph)>2)
    node.group <- igraph::components(fdf.iso.igraph)$membership
    xcms.fdf <- as.data.frame(xcms.fdf)
    rownames(xcms.fdf) <-xcms.fdf$feature_id
    message( length(unique(na.omit(node.group)))," iso-group"  )
    message( (length(node.group))," iso-features"  )

    for (i in seq_along(unique(node.group))) {

      #message(i)
      this.nodes <- names(which(node.group==i))
      this.fdf <- xcms.fdf[this.nodes,]
      this.iso <- fdf.iso.connect %>%
        dplyr::filter(from%in%this.nodes | to %in% this.nodes)

      this.igraph <- igraph::graph_from_data_frame(this.iso,vertices =this.fdf[,1:7] )
      #visNetwork::visIgraph(this.igraph)
      this.iso.assign <- get_iso_net_assign(this.igraph,net.degree.ratio = net.degree.ratio)
      xcms.fdf[names(this.iso.assign$iso.seed),
              "iso_seed"] <- this.iso.assign$iso.seed
      xcms.fdf[names(this.iso.assign$iso_count),
               "iso_count"] <- this.iso.assign$iso_count
      xcms.fdf[this.nodes,"iso_connection_group"] <- i
      this.fdf <- xcms.fdf[this.nodes,]


    }





  }


  ### save to featuredef
  {

    xcms.fdf.temp <- xcms::featureDefinitions(xcms.xcms)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    xcms.fdf.temp[,"iso_seed"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_seed"]
    xcms.fdf.temp[,"iso_count"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_count"]
    xcms.fdf.temp[,"iso_connection_group"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_connection_group"]
    xcms.fdf.temp -> xcms::featureDefinitions(xcms.xcms)
    message("Get ",
            sum(!is.na(xcms.fdf.temp[,"iso_count"])),
            " isotopologues")

  }

  return(xcms.xcms)

}

#' @describeIn xcms_extenstion identify isotopologues with multiple isotope tracers
#' @title Xcms Get Feature Isotopologues Multi Tracer
#' @description TODO: unfinished. Similar to `xcms_get_feature_isotopologues` but supports multiple isotope labels simultaneously (e.g., `[13]C` and `[15]N`).
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param iso_ele Character vector of isotope element strings (e.g., `c("[13]C","[15]N")`).
#' @param max_label Named numeric vector of maximum labels per tracer, names must match `iso_ele`.
#' @param ppm Mass accuracy tolerance in ppm (default 5).
#' @param rt.tol Retention time tolerance in seconds (default 5).
#' @param net.degree.ratio Ratio threshold for network degree to assign isotopologue seeds (default 0.3).
#' @return XCMSnExp object with featureDefinitions updated with iso_seed, iso_count, iso_connection_group, and per-tracer iso_count_* columns.
#' @export
#'
xcms_get_feature_isotopologues_multi_tracer <- function(xcms.xcms,
                                            iso_ele = c("[13]C","[15]N"),
                                            max_label = c("[13]C" = 30,"[15]N" = 10),
                                            ppm = 5,
                                            rt.tol = 5,
                                            net.degree.ratio = 0.3){

  ### find multi-tracer connections
  {
    fdf.iso.connect <- get_xcms_feature_iso_connection_multi_tracer(
      xcms.xcms,
      iso_ele = iso_ele,
      max_label = max_label,
      ppm = ppm,
      rt.tol = rt.tol
    )
  }

  ### assign isotopologues
  {
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms) %>%
      as.data.frame()
    iso_ele_clean <- get_ele_uniso(iso_ele)
    for (el in iso_ele_clean) {
      xcms.fdf[, paste0("iso_count_", el)] <- NA
    }
    xcms.fdf[, "iso_seed"] <- NA
    xcms.fdf[, "iso_count"] <- NA

    if (nrow(fdf.iso.connect) == 0) {
      message("No isotopologue connections found")
      fdf.iso.igraph <- igraph::make_empty_graph(n = 0)
      node.group <- integer(0)
    } else {
      fdf.iso.igraph <- igraph::graph_from_data_frame(fdf.iso.connect)
      node.group <- igraph::components(fdf.iso.igraph)$membership
    }

    xcms.fdf <- as.data.frame(xcms.fdf)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    message(length(unique(na.omit(node.group))), " iso-group")
    message(length(node.group), " iso-features")

    for (i in seq_along(unique(node.group))) {
      this.nodes <- names(which(node.group == i))
      this.fdf <- xcms.fdf[this.nodes, ]
      this.iso <- fdf.iso.connect[from %in% this.nodes | to %in% this.nodes]
      this.igraph <- igraph::graph_from_data_frame(this.iso, vertices = this.fdf[, 1:7])
      this.iso.assign <- get_iso_net_assign(this.igraph, net.degree.ratio = net.degree.ratio)

      xcms.fdf[names(this.iso.assign$iso.seed), "iso_seed"] <- this.iso.assign$iso.seed
      xcms.fdf[names(this.iso.assign$iso_count), "iso_count"] <- this.iso.assign$iso_count
      xcms.fdf[this.nodes, "iso_connection_group"] <- i

      ### assign per-tracer counts from the seed
      this.connect <- fdf.iso.connect[from %in% this.nodes & to %in% this.nodes]
      for (el in iso_ele_clean) {
        col_name <- paste0("closest.iso.count_", el)
        if (col_name %in% colnames(this.connect)) {
          seed.fid <- unique(na.omit(xcms.fdf[this.nodes, "iso_seed"]))
          if (length(seed.fid) == 1) {
            seed.idx <- which(this.connect$from == seed.fid | this.connect$to == seed.fid)
            if (length(seed.idx) > 0) {
              for (node in this.nodes) {
                if (node == seed.fid) {
                  xcms.fdf[node, paste0("iso_count_", el)] <- 0
                  next
                }
                edge.idx <- which((this.connect$from == node & this.connect$to == seed.fid) |
                                    (this.connect$from == seed.fid & this.connect$to == node))
                if (length(edge.idx) > 0) {
                  xcms.fdf[node, paste0("iso_count_", el)] <- abs(this.connect[edge.idx[1], col_name])
                }
              }
            }
          }
        }
      }
    }
  }

  ### save to featureDefinitions
  {
    xcms.fdf.temp <- xcms::featureDefinitions(xcms.xcms)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    xcms.fdf.temp[, "iso_seed"] <- xcms.fdf[xcms.fdf.temp$feature_id, "iso_seed"]
    xcms.fdf.temp[, "iso_count"] <- xcms.fdf[xcms.fdf.temp$feature_id, "iso_count"]
    xcms.fdf.temp[, "iso_connection_group"] <- xcms.fdf[xcms.fdf.temp$feature_id, "iso_connection_group"]
    for (el in iso_ele_clean) {
      xcms.fdf.temp[, paste0("iso_count_", el)] <- xcms.fdf[xcms.fdf.temp$feature_id, paste0("iso_count_", el)]
    }
    xcms.fdf.temp -> xcms::featureDefinitions(xcms.xcms)
    message("Get ",
            sum(!is.na(xcms.fdf.temp[, "iso_count"])),
            " isotopologues")
  }

  return(xcms.xcms)

}


#' @title TODO: unfinished. Build isotope mass shift grid for multi-tracer
#' @description Generates all label-count combinations across tracers with their mass shifts.
#' @param iso_ele Character vector of isotope element strings.
#' @param max_label Named numeric vector of maximum labels per tracer.
#' @return Data.frame with columns for each tracer's label count, total.count, and mass.shift.
#'
#' @title Pair xcms features within retention-time tolerance
#' @description Builds all directed feature pairs from `featureDefinitions` whose
#'   retention times differ by less than `rt.tol`. Used as the RT filter before
#'   isotope mass-shift matching.
#' @param xcms.xcms \code{XCMSnExp} with feature definitions.
#' @param rt.tol Retention time tolerance in seconds (default 5).
#' @return Data frame with integer `from`/`to` row indices, feature ids, m/z, rt,
#'   `mz.diff`, `rt.diff`, and `mz.mean` (mean m/z of the pair).
#' @export

get_xcms_feature_connect <- function(xcms.xcms,rt.tol = 5){


  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)

  {
    rt <- xcms.fdf$rtmed
    x <- data.table(rt = rt, start = rt,end = rt)
    y <- x[,.(id = seq_along(rt),rt,start = rt - rt.tol, end = rt + rt.tol)]
    data.table::setkey(y,start,end)
    rtm <- data.table::foverlaps(x, y, type="any", which=TRUE)
    rtm <- rtm[,.(xid,yid = y$id[yid])]
  }


  {
    xcms.net <- rtm[,.(from = xid, to = yid)]
    xcms.net <- xcms.net[from < to ][
      , rt.diff := (xcms.fdf$rtmed[to]-xcms.fdf$rtmed[from]) ][
        abs(rt.diff) < rt.tol,][
          , c("from.mz","to.mz") := .( xcms.fdf$mzmed[from], xcms.fdf$mzmed[to])][
            ,c("mz.diff","mz.mean") := .(to.mz-from.mz, (from.mz+to.mz)/2)]
  }


  return(xcms.net)



}




get_xcms_feature_isotope_grid_multi_tracer <- function(iso_ele, max_label) {

  iso.chemforms <- character()
  iso.counts <- numeric()

  for (i in seq_along(iso_ele)) {
    ele <- iso_ele[i]
    elem.symbol <- stringr::str_extract(ele, "[[:alpha:]]+")
    chemform <- paste0(ele, 1, elem.symbol, -1)
    iso.chemforms[i] <- chemform
    iso.counts[i] <- MSCC::chemform_mz(chemform, 0)
  }

  label.ranges <- lapply(max_label, function(m) -m:m)
  grid <- expand.grid(label.ranges)
  colnames(grid) <- iso_ele

  mass.shift <- as.matrix(grid) %*% iso.counts

  data.frame(
    grid,
    total.count = rowSums(abs(grid)),
    mass.shift = as.vector(mass.shift)
  )
}


get_xcms_feature_iso_connection_multi_tracer <- function(xcms.xcms,
                                                          iso_ele,
                                                          max_label,
                                                          ppm = 10,
                                                          rt.tol = 5) {

  fdf.connect <- get_xcms_feature_connect(xcms.xcms, rt.tol = rt.tol)

  iso.grid <- get_xcms_feature_isotope_grid_multi_tracer(iso_ele, max_label)

  match.res <- match_mz_foverlaps(mz1 = fdf.connect$mz.diff,
                                   mz2 = iso.grid$mass.shift,
                                   ppm.base = fdf.connect$mz.mean,
                                   ppm = ppm)

  iso_ele_clean <- get_ele_uniso(iso_ele)
  for (i in seq_along(iso_ele)) {
    fdf.connect[match.res$ion1, paste0("closest.iso.count_", iso_ele_clean[i])] <- iso.grid[match.res$ion2, iso_ele[i]]
  }
  fdf.connect[match.res$ion1, "closest.iso.count"] <- iso.grid$total.count[match.res$ion2]
  fdf.connect[match.res$ion1, "closest.iso.mz"] <- iso.grid$mass.shift[match.res$ion2]
  fdf.connect[match.res$ion1, "mz.error"] <- abs(match.res$mz.ppm)

  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  fdf.iso.connect <- fdf.connect[!is.na(closest.iso.count) & closest.iso.count != 0]
  data.table::setorder(fdf.iso.connect, mz.error)
  fdf.iso.connect <- fdf.iso.connect[!duplicated(fdf.iso.connect, by = c("from", "to"))]
  fdf.iso.connect[, `:=`(from = xcms.fdf$feature_id[from],
                          to = xcms.fdf$feature_id[to])]

  return(fdf.iso.connect)

}


get_xcms_feature_iso_connection <- function(xcms.xcms,
                                            iso_ele,
                                            max_label = 10,
                                            ppm = 10,
                                            rt.tol = 5){


  {

    fdf.connect <- get_xcms_feature_connect(xcms.xcms, rt.tol = rt.tol)

    #isotope <- "[13]C"
    iso.chemform <- paste0(iso_ele,1,
                           str_extract(string = iso_ele,pattern = "[[:alpha:]]+"),-1)
    iso.count <- -max_label:max_label
    iso.mz <- MSCC::chemform_mz(iso.chemform,0)*iso.count

    closest.iso.count <- sapply(fdf.connect$mz.diff, function(x){
      iso.count[which.min(abs(iso.mz-x))]
    } )
    fdf.connect <- fdf.connect%>%
      #rowwise()%>%
      dplyr::mutate(closest.iso.count = closest.iso.count)%>%
      dplyr::mutate(closest.iso.mz = iso.mz[match(closest.iso.count,iso.count)],
                    mz.error = abs(mz.diff-closest.iso.mz),
                    is.iso = mz.error/(from.mz+to.mz)*2 < ppm*1e-6)%>%
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

  return(fdf.iso.connect)

}



#' @describeIn xcms_extenstion calculate traced-isotopologue labeling ratios
#' @title Xcms Get Feature Traced Isotopologue
#' @description Calculates isotopologue-to-seed ratios and determines traced
#' isotopologues (label-enriched isotopologues) using one of two methods:
#' \itemize{
#'   \item \code{untraced_compare} (legacy): compare traced and untraced sample sources.
#'   \item \code{natural_based}: compare observed ratio to theoretical natural isotope
#'   ratio derived from \code{MSCC::chemform_isotopes_pattern_enviPat()}.
#' }
#' Results are written to featureDefinitions as \code{is_labeled} and
#' \code{Ratio_to_seed_*} columns (same output contract as the legacy function).
#' @param xcms.xcms XCMSnExp object with isotopologue assignments.
#' @param iso_ele Isotope element string (e.g., `"[13]C"`) used for labeling.
#' @param method Labeling method: \code{"untraced_compare"} or
#'   \code{"natural_based"}. (Legacy aliases \code{"method1"} /
#'   \code{"method2"} are also accepted.)
#' @param ... Additional arguments passed to internal functions.
#'
#' @return XCMSnExp object with featureDefinitions updated with is_labeled column and Ratio_to_seed_* columns.
#' @export
#'
xcms_get_feature_traced_isotopologue <- function(xcms.xcms,
                                                 iso_ele = "[13]C",
                                                 method = c("untraced_compare", "natural_based")[1],
                                                 ...){

  method <- as.character(method)[1]
  if (identical(method, "method1")) method <- "untraced_compare"
  if (identical(method, "method2")) method <- "natural_based"
  method <- match.arg(method, c("untraced_compare", "natural_based"))

  ### feature data
  {
    xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
    xcms.se <- xcms.se[,xcms.se$sample.type=="Sample"]
  }

  ### calc iso ratio to seed and aggregate by sample.source
  {
    xcms.ratio.to.seed <- get_xcms_iso_fraction(xcms.xcms)
    xcms.ratio.to.seed <- apply(xcms.ratio.to.seed, 1,
                                function(x){
                                  mean_f(x, f = xcms.se$sample.source, simplify = F, na.rm = T)
                                }) %>%
      do.call(bind_rows, .) %>%
      as.matrix()
  }

  ### determine traced-isotopologue labels
  {
    is.iso <- xcms.se$isotope_tracer %in% iso_ele
    sample.source.iso <- unique(xcms.se$sample.source[is.iso])
    sample.source.uniso <- unique(xcms.se$sample.source[!is.iso])

    if (method == "untraced_compare") {
      if (length(sample.source.uniso) == 0) {
        cli::cli_alert_warning("No untraced sample.source found; cannot determine traced isotopologues with {.code untraced_compare}.")
        is_labeled <- rep(NA, nrow(xcms.ratio.to.seed))
      } else {
      is_labeled <- apply(xcms.ratio.to.seed, 1, function(x){
        any(x[sample.source.iso] > mean(x[sample.source.uniso], na.rm = TRUE), na.rm = TRUE)
      })
      }
    } else {
      # method2: compare observed ratio to theoretical natural ratio.
      xcms.fdf <- xcms::featureDefinitions(xcms.xcms) %>% as.data.frame()

      # Expected natural ratio for one feature row.
      .expected_ratio <- function(formula, iso_count, iso_ele, thresh = 1e-6) {
        if (is.na(formula) || !nzchar(formula) || is.na(iso_count)) return(NA_real_)
        pat <- tryCatch(
          MSCC::chemform_isotopes_pattern_enviPat(formula, thresh = thresh),
          error = function(e) NULL
        )
        if (is.null(pat) || !nrow(pat) || !all(c("isotope_element", "abundance") %in% colnames(pat))) {
          return(NA_real_)
        }
        ie <- as.character(pat$isotope_element)
        ie[is.na(ie)] <- ""
        keep <- ie == "" | grepl(iso_ele, ie, fixed = TRUE)
        pat <- pat[keep, , drop = FALSE]
        if (!nrow(pat)) return(NA_real_)
        iso_chr <- as.character(pat$isotope_element)
        iso_chr[is.na(iso_chr)] <- ""
        iso_num <- ifelse(
          iso_chr == "",
          0L,
          suppressWarnings(as.integer(gsub("[^0-9]", "", gsub(iso_ele, "", iso_chr, fixed = TRUE))))
        )
        iso_num[is.na(iso_num)] <- 0L
        a0 <- sum(pat$abundance[iso_num == 0], na.rm = TRUE)
        ak <- sum(pat$abundance[iso_num == as.integer(iso_count)], na.rm = TRUE)
        if (!is.finite(a0) || a0 <= 0) return(NA_real_)
        ak / a0
      }

      fids <- if ("feature_id" %in% colnames(xcms.fdf)) as.character(xcms.fdf$feature_id) else rownames(xcms.fdf)
      formula_vec <- if ("formula" %in% colnames(xcms.fdf)) as.character(xcms.fdf$formula) else rep(NA_character_, nrow(xcms.fdf))
      iso_count_vec <- if ("iso_count" %in% colnames(xcms.fdf)) suppressWarnings(as.integer(xcms.fdf$iso_count)) else rep(NA_integer_, nrow(xcms.fdf))
      expected <- vapply(seq_len(nrow(xcms.fdf)), function(i) {
        .expected_ratio(formula_vec[i], iso_count_vec[i], iso_ele)
      }, numeric(1))
      names(expected) <- fids

      # If no traced sample.source is available, fallback to all sample.source.
      use.sources <- sample.source.iso
      if (!length(use.sources)) {
        use.sources <- colnames(xcms.ratio.to.seed)
      }
      use.sources <- intersect(use.sources, colnames(xcms.ratio.to.seed))

      is_labeled <- rep(NA, nrow(xcms.ratio.to.seed))
      if (length(use.sources)) {
        is_labeled <- vapply(seq_len(nrow(xcms.ratio.to.seed)), function(i) {
          obs <- xcms.ratio.to.seed[i, use.sources, drop = TRUE]
          fid <- rownames(xcms.ratio.to.seed)[i]
          exp_i <- expected[fid]
          exp_i <- if (length(exp_i)) exp_i[[1]] else NA_real_
          if (!is.finite(exp_i)) return(NA)
          any(obs > exp_i, na.rm = TRUE)
        }, logical(1))
      }
    }
  }

  ### import to xcms
  {
    xcms.fda <- xcms::featureDefinitions(xcms.xcms)
    xcms.fda$is_labeled <- is_labeled
    colnames(xcms.ratio.to.seed) <- paste0("Ratio_to_seed_", colnames(xcms.ratio.to.seed))
    xcms.fda[, colnames(xcms.ratio.to.seed)] <- xcms.ratio.to.seed
    xcms.fda -> xcms::featureDefinitions(xcms.xcms)
    message("Get ",
            sum(xcms.fda$is_labeled, na.rm = TRUE),
            " isotope label")
  }

  return(xcms.xcms)
}

#' @title Xcms Get Feature Isotope Label (Deprecated)
#' @description Deprecated wrapper of \code{xcms_get_feature_traced_isotopologue()}.
#' @param xcms.xcms XCMSnExp object with isotopologue assignments.
#' @param iso_ele Isotope element string (e.g., `"[13]C"`).
#' @param ... Additional arguments.
#' @return XCMSnExp object.
#' @export
xcms_get_feature_isotope_label <- function(xcms.xcms,
                                           iso_ele = "[13]C",
                                           ...){
  .Deprecated("xcms_get_feature_traced_isotopologue",
              package = "MSdev",
              msg = "xcms_get_feature_isotope_label is deprecated. Use xcms_get_feature_traced_isotopologue().")
  xcms_get_feature_traced_isotopologue(xcms.xcms = xcms.xcms,
                                       iso_ele = iso_ele,
                                       method = "untraced_compare",
                                       ...)
}


get_xcms_isotopologues_report <- function(xcms.xcms){



}


#' @describeIn xcms_extenstion calculate isotopologue fractions
#' @title Get Xcms Iso Fraction
#' @description Calculates the fraction of isotopologue intensities relative to their seed feature intensities for each sample. Returns a matrix of fractions without natural abundance adjustment.
#' @param xcms.xcms XCMSnExp object with isotopologue assignments (iso_seed column).
#'
#' @return Matrix with rows as features and columns as samples, containing intensity ratios to seed features.
#' @export
#'
get_xcms_iso_fraction <- function(xcms.xcms){


  ### feature data
  {
    xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
    #xcms.se <- xcms.se[,xcms.se$sample.type=="Sample"]
     xcms.rda <- rowData(xcms.se)%>%
      as.data.frame()
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
      this.matrix <- xcms.val[this.iso$feature_id,,drop = F]
      this.matrix <- t(t(this.matrix)/this.matrix[this.fid,])
      xcms.ratio.to.seed[rownames(this.matrix),] <- this.matrix
    }

    xcms.ratio.to.seed[is.nan(xcms.ratio.to.seed)] <- 0

  }
  return(xcms.ratio.to.seed)

}


#' @describeIn xcms_extenstion match features to compound database
#' @title Xcms Get Feature Ms1 Candidate
#' @description Matches features in an XCMSnExp object to compounds in a CompoundDb database using m/z and retention time tolerance. Calculates adduct masses for each compound and finds matches within specified ppm error. Results are stored as candidate lists in featureDefinitions.
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param cpdb CompoundDb object containing compound database.
#' @param mz.ppm Numeric. Mass accuracy tolerance in parts per million (default 10).
#' @param rt.tol Numeric. Retention time tolerance in seconds (default Inf, no RT filtering).
#' @param selected_adduct Character vector of adducts to consider (default from MSCC::adduct.table$Adduct).
#' @param ... Additional arguments passed to internal functions.
#'
#' @return XCMSnExp object with featureDefinitions updated with candidate.id, candidate.formula, candidate.adduct, and candidate.mz columns.
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
  #cpdbt <- dplyr::filter(cpdbt,lipidclass=="VAE")
  if ("has_sp"%in% colnames(cpdbt))  cpdbt <- cpdbt[cpdbt$has_sp>0,]
  cpdbt$formula <- MSCC::chemform_formate(cpdbt$formula)

  adducts <-  MSCC::chemform_adduct_check(selected_adduct)%>%
    dplyr::mutate(polarity = case_when(Ion_mode == "negative"~0,T~1))%>%
    dplyr::filter(polarity %in% polarity(xcms.xcms))
  cp.adduct <- MSCC::chemform_adduct(cpdbt$formula,
                                     adducts$adduct.formated,
                                     value = "all" )
  cp.adduct <- cp.adduct%>%
    dplyr::mutate(compound_id=cpdbt$compound_id[id]  )%>%
    dplyr::filter( findInterval(chemform.adduct.mz,
                                mzrange(xcms.xcms))==1)

  ### match database
  xcms.featuredef <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()

  matched.df <- match_mz_foverlaps(mz1 = xcms.featuredef$mzmed,
                            mz2 = cp.adduct$chemform.adduct.mz,
                            ppm = mz.ppm)
  matched.df2 <- cbind( matched.df,cp.adduct[matched.df$ion2,])
  xcms.featuredef$candidate.id <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$compound_id[as.numeric(idx)]
  })
  xcms.featuredef$candidate.formula <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform[as.numeric(idx)]
  })
  xcms.featuredef$candidate.adduct <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$adduct[as.numeric(idx)]
  })
  xcms.featuredef$candidate.mz <- sapply(1:nrow(xcms.featuredef),function(i){
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform.adduct.mz[as.numeric(idx)]
  })

  xcms::featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(xcms.featuredef)

  return(xcms.xcms)

}


xcms_get_feature_ms2_score <- function(xcms.xcms ,
                                       cpdb,
                                       sp.ms2,
                                       ...){



  ### no ms2
  {
    if (length(sp.ms2)==0) {
      xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
      xcms.fdf$score.ms2 <- lapply(xcms.fdf$candidate.id,function(x){
        rep(0,length(x))
      })

      xcms::featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(xcms.fdf)
      return(xcms.xcms)
    }

  }


  ### load spectra database
  {
    Spectra_database <- Spectra(cpdb)
    Spectra_database <- Spectra_set_MEM_backend(Spectra_database)
    Spectra_database <- filterPolarity(Spectra_database,
                                       unique(polarity(xcms.xcms)))

    #spectraNames(Spectra_database) <- Spectra_database$compound_id
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  }

  ### sp process
  {
    Spectra_database <- Spectra_database%>%
      filterSpectra_below_PrecursorMz()%>%
      normalizeSpectra(norm_to = "max")%>%
      filterSpectraIntensity(ratio = 0.05)%>%
      Spectra::applyProcessing()
    if(length(sp.ms2)!=0){
      sp.ms2 <- sp.ms2%>%
        filterSpectra_below_PrecursorMz()%>%
        normalizeSpectra(norm_to = "max")%>%
        filterSpectraIntensity(ratio = 0.05)%>%
        Spectra_set_MEM_backend()%>%
        Spectra::applyProcessing()
      #if ("from_iso" %in% spectraVariables(sp.ms2)) {
      #  sp.ms2 <- sp.ms2[!sp.ms2$from_iso]
      #}
    }
    if (!all(unlist(xcms.fdf$candidate.id)%in%
             Spectra_database$compound_id)) {
      sp.empty <- makeEmptySpectra(compound_id= setdiff(unlist(xcms.fdf$candidate.id),
                                                        Spectra_database$compound_id))
      Spectra_database <- c(Spectra_database,sp.empty)
    }
  }

  ### sp ms2 split
  {
    sp.exp <- sapply(1:nrow(xcms.fdf),function(i){

      x <- xcms.fdf$ms2_id[[i]]
      if (!length(x)) return(NULL)
      sp.ms2[match(x,spectraNames(sp.ms2))]
      #if (length(x)==0) {
      #  sp <- makeSpectra(xcms.fdf$mzmed[i],
      #                    xcms.fdf$rtmed[i])
      #}else
      #  sp <- list(sp.ms2[x])
      #return(sp)
    })
    names(sp.exp) <- xcms.fdf$feature_id


  }

  ### sp ref split
  {
    sp.split.df <- lapply(1:nrow(xcms.fdf),
                          function( i ){
                            i.candi <- xcms.fdf$candidate.id[[i]]
                            i.adduct <- xcms.fdf$candidate.adduct[[i]]
                            if (!length(i.candi)) return(NULL)
                            i.df <- lapply(i.candi,
                                           function(x){which(Spectra_database$compound_id==x)}
                            )%>%
                              `names<-`(xcms.fdf$candidate.id[[i]])%>%
                              unlist_to_df(name_to = "compound_id",
                                           value_to = "sp_id")
                            i.df$adduct <- i.adduct[match(i.df$compound_id,i.candi)]
                            i.df
                          })%>%
      data.table::rbindlist(use.names = T,idcol = "feature_id")
    sp.ref <- Spectra_database[sp.split.df$sp_id]
    sp.ref$adduct <- sp.split.df$adduct
    sp.ref <- split(sp.ref,xcms.fdf$feature_id[sp.split.df$feature_id])



  }

  ### output all candidate score
  {
    .f <- function(expSpec,refSpec,...){
      if (is.null(expSpec)) {
        if (is.null(refSpec)) {
          return(NULL)
        }else{
          x <- unique(paste0(refSpec$compound_id,"_",refSpec$adduct))
          y <- rep(0,length(x))
          names(y )<-x
          return(y)
        }
      }
      if (is.null(refSpec)) return(NULL)
      scorem <- Spectra::compareSpectra(expSpec,refSpec,,FUN = MsCoreUtils::ndotproduct, m = 2)
      dim(scorem) <- c(length(expSpec),length(refSpec))
      scorem[is.infinite(scorem)|is.na(scorem )] <- 0
      scores <- apply(scorem,2,max,na.rm=T)
      (mean_f(scores,paste0(refSpec$compound_id,"_",refSpec$adduct)))
    }
    xcms.fdf$score.ms2 <- BiocParallel::bplapply(1:length(sp.exp),
                                                       function(i){
        fid <-xcms.fdf$feature_id[i]
      s <- .f(expSpec = sp.exp[[i]], refSpec = sp.ref[[fid]])
      s[paste0(xcms.fdf$candidate.id[[i]],"_",
               xcms.fdf$candidate.adduct[[i]])]%>%
        unname()
    },BPPARAM = BiocParallel::SerialParam(
      progressbar = T))


  }


  xcms::featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(xcms.fdf)

  return(xcms.xcms)

}

xcms_get_feature_isopattern_score <- function(xcms.xcms,
                                              ppm = 10,
                                              calc_isopattern_score = T,
                                              BPPARAM = SerialParam(progressbar = T)){

  ### data to calc isopattern
  {
    xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
    xcms.se <- xcms.se[,is.na(xcms.se$isotope_tracer)]
    xcms.se <- xcms.se[,!xcms.se$sample.type%in% "Blank"]
  }

  ### calculate iso-pattern score
  {

    if (calc_isopattern_score) {

      iso.score <- bplapply(seq_len(nrow(xcms.se)),
             FUN = function(i,xcms.se,ppm){
               xcms.fdf <- SummarizedExperiment::rowData(xcms.se)
               xcms.se.temp <- xcms.se[MSdev:::between.range(xcms.fdf$rtmed,
                                                     c(xcms.fdf$rtmed[i]-10,
                                                       xcms.fdf$rtmed[i]+10)),]
               formulas <-mapply( MSCC::chemform_adduct,
                                 chemform = xcms.fdf$candidate.formula[[i]],
                                 adduct =  xcms.fdf$candidate.adduct[[i]],
                                 value = "chemform")
               get_isopattern_score(formula = formulas,
                                    mzs = rowData(xcms.se.temp)$mzmed,
                                    int_matrix = assay(xcms.se.temp),
                                    ppm = ppm)
             },xcms.se=xcms.se,ppm=ppm,
             BPPARAM =BPPARAM)

    }else{iso.score <- 0}



  }


  ### return
  {
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
    xcms.fdf$score.isopattern <-iso.score
    xcms.fdf -> xcms::featureDefinitions(xcms.xcms)
    return(xcms.xcms)
  }


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
                                        cpdb,
                                        cpdb.keys = c("name","formula","smiles"),
                                        weight_mz = 0.1,
                                        weight_ms2 = 0.7,
                                        weight_isopattern =0.2,
                                        ...){


  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  xcms.fdf$compound_id <- NA
  xcms.fdf$adduct <- NA
  xcms.fdf$score <- NA
  xcms.fdf$mz_ref <- NA



  {
    xcms.candi.dt <-lapply(which(lengths(xcms.fdf$candidate.id)!=0), function(i) {
        data.table(
          feature_id = xcms.fdf$feature_id[i],
          mz = xcms.fdf$mzmed[i],
          candidate.id = xcms.fdf$candidate.id[[i]],
          candidate.adduct = xcms.fdf$candidate.adduct[[i]],
          candidate.formula = xcms.fdf$candidate.formula[[i]],
          candidate.mz = xcms.fdf$candidate.mz[[i]],
          score.ms2 =  xcms.fdf$score.ms2[[i]],
          score.isopattern = xcms.fdf$score.isopattern[[i]]
        )
      }) %>% rbindlist
    xcms.candi.dt <- xcms.candi.dt[
      ,score.isopattern:= ifelse(is.na(score.isopattern),0,score.isopattern)][
        ,score.isopattern:= ifelse(is.nan(score.isopattern),0,score.isopattern)][
          ,score.mz := 1- abs(candidate.mz-mz)/mz *1e6 / 20 ][ ### (1-ppm/10)
            , score.mz := ifelse(score.mz < 0, 0, score.mz)][
              ,score := score.isopattern * weight_isopattern + score.ms2*weight_ms2+score.mz * weight_mz]

    xcms.candi.dt.max <- xcms.candi.dt[, .SD[which.max(score)], by = feature_id]
    data.table::setnames(xcms.candi.dt.max,
                          old =  c("candidate.id","candidate.adduct","candidate.formula","candidate.mz"),
                          new = c("compound_id","adduct","formula","mz_ref") )
    data.table::setkey(xcms.candi.dt.max,feature_id)

    xcms.candi <- xcms.candi.dt.max[xcms.fdf$feature_id]
    xcms.fdf$compound_id <- xcms.candi$compound_id
    xcms.fdf$adduct <- xcms.candi$adduct
    xcms.fdf$formula <- xcms.candi$formula
    xcms.fdf$score <- xcms.candi$score
    xcms.fdf$mz_ref <- xcms.candi$mz_ref
    xcms.fdf$score.ms2 <- xcms.candi$score.ms2
    xcms.fdf$score.isopattern <- xcms.candi$score.isopattern
    xcms.fdf$score.mz <- xcms.candi$score.mz

  }


  ### Comp info
  {

    dbinfo <- get_CompDb_info(cpdb,
                              xcms.fdf$compound_id,
                              keys = cpdb.keys)
    xcms.fdf[,colnames(dbinfo)] <- dbinfo

  }


  xcms::featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(xcms.fdf)

  return(xcms.xcms)



}


get_xcms_feature_definitions <- function(xcms.xcms){
  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::select(
      !c(mzmin,mzmax,rtmin,rtmax,npeaks,peakidx)
    )
  return(xcms.fdf)

}


find_xcms_feature <- function(xcms.xcms,mz = 100,ppm = 10){

  fdf <- xcms::featureDefinitions(xcms.xcms)
  mzr <- mz.range.ppm(mz,ppm)
  fdf[between(fdf$mzmed,mzr[1],mzr[2] ), ]%>%
    as_tibble()
}


find_xcms_peaks <- function(xcms.xcms, mz = 100, ppm = 10) {

  peaks <- xcms::chromPeaks(xcms.xcms)
  mzr <- mz.range.ppm(mz, ppm)
  peaks[between(peaks[,"mz"], mzr[1], mzr[2]), ]
}


#' @describeIn xcms_extenstion plot peaks distribution
#' @title plot_xcms_peaks_distribution
#' @description export peaks data by xcms::chromPeaks and plot by ggplot2
#'
#' @param xcms.xcms XCMSnExp object
#' @param plot.title title
#' @param type `"o"`, for geom_point, `"l"`, for geom_segment
#'
#' @return ggplot object
#' @export
#'

plot_xcms_peaks_distribution <- function(xcms.xcms,plot.title = "Peaks distribution",type = "o"){

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    dplyr::mutate(as.data.frame(chromPeakData(xcms.xcms)),
                  peak_id = rownames(.),
                  merged = grepl(peak_id,pattern = "CPM"))%>%
    dplyr::filter(!is.na(maxo),
                  rtmax-rtmin <60,
                  !merged)
  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
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



#' @describeIn xcms_extenstion plot features distribution
#' @title Plot Xcms Features Distribution
#' @description Visualizes the distribution of detected features in a 2D space of retention time (x-axis) vs m/z (y-axis). Point size represents peak width, color represents log10 intensity. Includes peak detection parameters in subtitle.
#' @param xcms.xcms XCMSnExp object with feature definitions.
#' @param plot.title Character title for the plot (default "Features distribution").
#'
#' @return ggplot object.
#' @export
#'

plot_xcms_features_distribution <-
  function(xcms.xcms, plot.title = "Features distribution") {
    xcms.features <- xcms::featureDefinitions(xcms.xcms) %>%
      as.data.frame() %>%
      mutate(mz = mzmed, rt = rtmed)
    xcms.features$maxo <-
      apply(xcms::featureValues(xcms.xcms, value = "maxo"), 1, median, na.rm = T)

    xcms.process.type <-
      xcms::processHistory(xcms.xcms) %>% sapply(processType)
    xcms.findpeak.param <-
      xcms:: processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]] %>%
      xcms::processParam()
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


xcms_remove_feature_var <- function(xcms.xcms,var){

  xcms.fdf<- xcms::featureDefinitions(xcms.xcms)
  var.selected <- setdiff(colnames(xcms.fdf),var)
  xcms.fdf <- xcms.fdf[,var.selected]
  xcms.fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)
}


#' @title plot_xcms_feature_chromatogram
#' @description extract Chromatogram from xcms according to feature's mz range and plot
#' @describeIn xcms_extenstion plot feature chromatogram
#' @param xcms.xcms XCMSnExp object
#' @param feature.id feature id
#' @param sampleNames sample names to include
#'
#' @return ggplot object
#' @export
#'

plot_xcms_feature_chromatogram <- function(xcms.xcms ,feature.id, sampleNames =NULL ){

  ### select samples
  all.sample.names <- Biobase::sampleNames(xcms.xcms)
  xcms.sample.info <- as.data.frame(Biobase::pData(xcms.xcms), stringsAsFactors = FALSE)
  if (!"sampleNames" %in% colnames(xcms.sample.info)) {
    xcms.sample.info$sampleNames <- all.sample.names
  }
  rownames(xcms.sample.info) <- all.sample.names

  if (is.null(sampleNames)) {
    sampleNames <- all.sample.names
  }
  sampleNames <- as.character(sampleNames)
  sampleNames <- intersect(sampleNames, all.sample.names)
  if (!length(sampleNames)) {
    stop("No valid sampleNames found in xcms.xcms.")
  }

  xcms.sample.info <- xcms.sample.info[sampleNames, , drop = FALSE]
  if (nrow(xcms.sample.info) > 5) {
    if ("group" %in% colnames(xcms.sample.info) && any(!is.na(xcms.sample.info$group))) {
      xcms.sample.info.sub <- xcms.sample.info %>%
        dplyr::group_by(group) %>%
        dplyr::slice_sample(n = 1) %>%
        as.data.frame(stringsAsFactors = FALSE)
    } else {
      xcms.sample.info.sub <- xcms.sample.info[seq_len(5), , drop = FALSE]
    }
  } else {
    xcms.sample.info.sub <- xcms.sample.info
  }

  sample.idx <- which(all.sample.names %in% xcms.sample.info.sub$sampleNames)
  if (!length(sample.idx)) {
    stop("No samples selected for chromatogram extraction.")
  }
  xcms.sub <- MSnbase::filterFile(xcms.xcms, sample.idx)

  ### mz / rt from feature peaks
  xcms.fdef <- xcms::featureDefinitions(xcms.xcms)
  if (is.numeric(feature.id)) {
    feature.id <- rownames(xcms.fdef)[feature.id]
  }
  feature.id <- as.character(feature.id)[1]
  if (is.na(feature.id) || !nzchar(feature.id) || !(feature.id %in% rownames(xcms.fdef))) {
    stop("feature.id does not exist in featureDefinitions(xcms.xcms).")
  }

  xcms.feature <- xcms.fdef[feature.id, , drop = FALSE]
  peak.idx <- xcms.feature$peakidx[[1]]
  if (is.null(peak.idx) || !length(peak.idx)) {
    stop("Selected feature has no linked peaks (empty peakidx).")
  }
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)[peak.idx, , drop = FALSE]
  mz.range <- as.numeric(c(min(xcms.peaks[, "mzmin"], na.rm = TRUE),
                           max(xcms.peaks[, "mzmax"], na.rm = TRUE)))
  rt.range <- as.numeric(c(min(xcms.peaks[, "rtmin"], na.rm = TRUE),
                           max(xcms.peaks[, "rtmax"], na.rm = TRUE)))
  if (length(mz.range) != 2 || length(rt.range) != 2 ||
      any(!is.finite(mz.range)) || any(!is.finite(rt.range))) {
    stop("Failed to derive finite mz/rt ranges from feature peaks.")
  }

  xcms.chrom <- xcms::chromatogram(xcms.sub,
                                   mz = mz.range,
                                   rt = rt.range,
                                   BPPARAM = BiocParallel::SerialParam())
  if (!methods::is(xcms.chrom, "XChromatograms")) {
    stop("Chromatogram extraction failed for selected feature and samples.")
  }

  xcms.chrom.data <- get_chroms_data(xcms.chrom)%>%
    dplyr::mutate(group = Biobase::sampleNames(xcms.sub)[col])

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

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    mutate(ppm = (mzmax-mzmin)/mz*1e6,
           mz_diff = mzmax-mzmin,
           peak_width = rtmax-rtmin)

  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
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
#' @describeIn xcms_extenstion plot MS1 scan counts for peaks
#' @param xcms.xcms XCMSnExp object should be a `XCMSnExp` object after `findChromPeaks`
#' @param plot.title title
#'
#' @return ggplot object
#' @export
#'

plot_xcms_peaks_ms1_scans <- function(xcms.xcms,plot.title = "Peaks Sans of MS1"){

  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
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



#' Plot number of MS2 scans overlapping each chromatographic peak
#'
#' @title Plot Xcms Peaks Ms2 Scans
#' @description Visualizes the number of MS2 scans that overlap each chromatographic peak based on retention time and m/z ranges. Produces a scatter plot with jitter, violin distribution, and counts of peaks with 0-5 MS2 scans.
#' @param xcms.xcms XCMSnExp object with detected peaks and MS2 scans.
#' @param plot.title Character title for the plot (default "Peaks Sans of MS2").
#'
#' @describeIn xcms_extenstion plot MS2 scan counts for peaks
#' @return ggplot object.
#' @export
#'

plot_xcms_peaks_ms2_scans <- function(xcms.xcms,plot.title = "Peaks Sans of MS2"){

  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
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
   #geom_vline(xintercept = ms1.rt,linewidth = 0.05,col = "black")+
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

 open_plot_win(peaks.dis.plot,width = 25,height = 5)
 peaks.dis.plot
}



plot_xcms_peaks_SN_distribution <- function(xcms.xcms,plot.title = "Peaks SNR(Signal to Noise Ratio)"){


  xcms.process.type <- xcms::processHistory(xcms.xcms) %>% sapply( xcms::processType )
  xcms.findpeak.param <- xcms::processHistory(xcms.xcms)[[which(xcms.process.type == "Peak detection")]]%>%
    xcms::processParam()
  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
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
#' @describeIn xcms_extenstion plot chromatogram for a peak
#' @param xcms.xcms XCMSnExp object
#' @param peak_id peak id
#' @param rt expansion range for rt
#'
#' @return ggplot object
#' @export
#'

plot_xcms_peaks_Chromatogram <- function(xcms.xcms,peak_id,rt = "expand"){

  peaks.data <- xcms::chromPeaks(xcms.xcms)[peak_id,,drop = F]
  peak_id <- rownames(peaks.data)
  mz.range <- c(peaks.data[,c("mzmin","mzmax")])
  rt.range <- c(peaks.data[,c("rtmin","rtmax")])
  xcms.chrom <- get_xcms_peaks_chrom(xcms.xcms,
                                     peaks.id = peak_id,
                                     rt.range = rt)
  chrom.data <- get_chroms_data(xcms.chrom)%>%
    dplyr::mutate(fill = rt > min(rt.range)&rt <max(rt.range),
                  sample = Biobase::sampleNames(xcms.xcms)[col]
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


#' @title Build xcms centWave roiList from mz/rt targets
#' @description
#' Construct a \code{roiList} accepted by \code{xcms::CentWaveParam(roiList = ...)}.
#' Input a matrix/data.frame with columns \code{mz} and \code{rt} (seconds). For each
#' target, \code{mzmin/mzmax} are calculated using ppm tolerance and the RT window
#' \code{rtmin/rtmax} is mapped to scan indices. Because scan numbers differ across
#' samples/files, the returned ROI uses the **union scan range** (min start to max end)
#' across all files for that RT window.
#'
#' @param mzrt matrix/data.frame with columns \code{mz} and \code{rt}.
#' @param xcms.xcms \code{XCMSnExp} object used to map RT to scan indices.
#' @param ppm numeric, ppm tolerance for mz window.
#' @param rt_tol numeric, RT tolerance in seconds.
#' @param ion_mode optional integer 1 (positive) or 0 (negative). If NULL, inferred
#'   from \code{Biobase::fData(xcms.xcms)$polarity}; must be unique.
#'
#' @return list of ROI objects (each ROI is a list with \code{scmin, scmax, mzmin, mzmax, length, intensity}).
#' @export
get_xcms_roi_list <- function(mzrt,
                              xcms.xcms,
                              ppm = 10,
                              rt_tol = 30,
                              ion_mode = NULL) {

  if (is.null(mzrt) || length(mzrt) == 0) return(list())
  mzrt <- as.data.frame(mzrt, stringsAsFactors = FALSE)
  if (!all(c("mz", "rt") %in% colnames(mzrt))) {
    stop("mzrt must have columns: mz, rt")
  }
  mzrt$mz <- suppressWarnings(as.numeric(mzrt$mz))
  mzrt$rt <- suppressWarnings(as.numeric(mzrt$rt))
  mzrt <- mzrt[is.finite(mzrt$mz) & is.finite(mzrt$rt), , drop = FALSE]
  if (!nrow(mzrt)) return(list())

  fdat <- Biobase::fData(xcms.xcms)
  if (is.null(fdat) || nrow(fdat) == 0) {
    stop("xcms.xcms has empty fData; cannot derive scan indices for roiList.")
  }
  need_cols <- c("fileIdx", "msLevel", "retentionTime", "polarity")
  if (!all(need_cols %in% colnames(fdat))) {
    stop("xcms.xcms fData must contain columns: ", paste(need_cols, collapse = ", "), ".")
  }

  if (is.null(ion_mode)) {
    ion_mode <- unique(fdat$polarity[fdat$msLevel == 1])
    ion_mode <- ion_mode[!is.na(ion_mode)]
    if (length(ion_mode) != 1) {
      stop("Cannot infer ion_mode (multiple polarities in MS1 scans). Provide ion_mode = 0/1.")
    }
  }

  files <- seq_along(MSnbase::fileNames(xcms.xcms))
  ms1_rt_by_file <- lapply(files, function(fi) {
    idx <- which(fdat$fileIdx == fi & fdat$msLevel == 1 & fdat$polarity == ion_mode)
    if (!length(idx)) return(NULL)
    rt <- fdat$retentionTime[idx]
    rt[order(rt)]
  })

  roi_list <- vector("list", length = nrow(mzrt))
  kept <- logical(nrow(mzrt))

  for (i in seq_len(nrow(mzrt))) {
    mz <- mzrt$mz[[i]]
    rt <- mzrt$rt[[i]]
    mzmin <- mz - mz * ppm / 1e6
    mzmax <- mz + mz * ppm / 1e6
    rtmin <- rt - rt_tol
    rtmax <- rt + rt_tol

    scmins <- integer()
    scmaxs <- integer()
    for (fi in files) {
      rts <- ms1_rt_by_file[[fi]]
      if (is.null(rts)) next
      in_rt <- which(rts >= rtmin & rts <= rtmax)
      if (!length(in_rt)) next
      scmins <- c(scmins, min(in_rt))
      scmaxs <- c(scmaxs, max(in_rt))
    }

    if (!length(scmins)) next
    scmin <- min(scmins)
    scmax <- max(scmaxs)
    roi_list[[i]] <- list(
      scmin = as.integer(scmin),
      scmax = as.integer(scmax),
      mzmin = as.numeric(mzmin),
      mzmax = as.numeric(mzmax),
      length = as.integer(scmax - scmin + 1L),
      intensity = 0
    )
    kept[[i]] <- TRUE
  }

  roi_list[kept]
}


xcms_filter_peaks_NA <- function(xcms.xcms, verbose = TRUE) {

  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) {
      message("xcms_filter_peaks_NA: no chromPeaks to filter")
    }
    return(xcms.xcms)
  }

  mz_col <- NULL
  if ("mz" %in% colnames(pks)) mz_col <- "mz"
  if (is.null(mz_col) && "mzmed" %in% colnames(pks)) mz_col <- "mzmed"
  if (is.null(mz_col)) {
    if (isTRUE(verbose)) {
      message("xcms_filter_peaks_NA: no mz column found in chromPeaks")
    }
    return(xcms.xcms)
  }

  mzv <- as.numeric(pks[, mz_col])
  bad <- is.na(mzv) | is.nan(mzv)
  n_bad <- sum(bad)
  n_total <- nrow(pks)
  ratio <- if (n_total == 0) 0 else n_bad / n_total

  if (isTRUE(verbose)) {
    message(sprintf(
      "xcms_filter_peaks_NA: %d/%d (%.4f) chromPeaks have mz NA/NaN",
      n_bad, n_total, ratio
    ))
  }

  if (n_bad > 0) {
    pks2 <- pks[!bad, , drop = FALSE]
    xcms::chromPeaks(xcms.xcms) <- pks2
  }

  return(xcms.xcms)
}

filter_xcms_chromPeaks_mz_width <- function(xcms.xcms, ppm = 20, verbose = TRUE) {
  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) message("filter_xcms_chromPeaks_mz_width: no chromPeaks to filter")
    return(xcms.xcms)
  }
  need <- c("mz", "mzmin", "mzmax")
  if (!all(need %in% colnames(pks))) {
    if (isTRUE(verbose)) {
      message("filter_xcms_chromPeaks_mz_width: required columns not found (mz/mzmin/mzmax)")
    }
    return(xcms.xcms)
  }

  mz <- suppressWarnings(as.numeric(pks[, "mz"]))
  mzmin <- suppressWarnings(as.numeric(pks[, "mzmin"]))
  mzmax <- suppressWarnings(as.numeric(pks[, "mzmax"]))
  mz_width_ppm <- (mzmax - mzmin) / pmax(mz, 1e-12) * 1e6
  bad <- !is.finite(mz_width_ppm) | (mz_width_ppm > ppm)
  n_bad <- sum(bad, na.rm = TRUE)
  n_total <- nrow(pks)

  if (isTRUE(verbose)) {
    message(sprintf(
      "filter_xcms_chromPeaks_mz_width: %d/%d chromPeaks removed (mz width > %.2f ppm); this is aggressive and may remove real peaks",
      n_bad, n_total, ppm
    ))
  }
  if (n_bad > 0) {
    xcms::chromPeaks(xcms.xcms) <- pks[!bad, , drop = FALSE]
  }
  xcms.xcms
}

#' @title Fix overly wide xcms chromPeaks mz window
#' @description
#' For chromPeaks with abnormal m/z width larger than \code{ppm}, this function
#' does not remove peaks. Instead, it recalculates and replaces \code{mzmin} and
#' \code{mzmax} around peak center \code{mz} so final width equals the target
#' ppm window.
#'
#' @param xcms.xcms XCMSnExp object.
#' @param ppm numeric ppm threshold/target width, default \code{20}.
#' @param verbose logical, print summary message.
#'
#' @return XCMSnExp object with fixed \code{chromPeaks} m/z windows.
#' @export
fix_xcms_chromPeaks_mz_width <- function(xcms.xcms, ppm = 20, verbose = TRUE) {
  pks <- xcms::chromPeaks(xcms.xcms)
  if (is.null(pks) || length(pks) == 0 || nrow(pks) == 0) {
    if (isTRUE(verbose)) message("fix_xcms_chromPeaks_mz_width: no chromPeaks to fix")
    return(xcms.xcms)
  }
  need <- c("mz", "mzmin", "mzmax")
  if (!all(need %in% colnames(pks))) {
    if (isTRUE(verbose)) {
      message("fix_xcms_chromPeaks_mz_width: required columns not found (mz/mzmin/mzmax)")
    }
    return(xcms.xcms)
  }

  mz <- suppressWarnings(as.numeric(pks[, "mz"]))
  mzmin <- suppressWarnings(as.numeric(pks[, "mzmin"]))
  mzmax <- suppressWarnings(as.numeric(pks[, "mzmax"]))
  bad_mz <- !is.finite(mz) | mz <= 0
  mz[bad_mz] <- (mzmin[bad_mz] + mzmax[bad_mz]) / 2
  bad_mz <- !is.finite(mz) | mz <= 0
  mz_width_ppm <- (mzmax - mzmin) / pmax(mz, 1e-12) * 1e6
  to_fix <- !bad_mz & is.finite(mz_width_ppm) & (mz_width_ppm > ppm)
  n_fix <- sum(to_fix, na.rm = TRUE)
  n_total <- nrow(pks)

  if (n_fix > 0) {
    half <- mz[to_fix] * as.numeric(ppm) / 2e6
    pks[to_fix, "mzmin"] <- mz[to_fix] - half
    pks[to_fix, "mzmax"] <- mz[to_fix] + half
    xcms::chromPeaks(xcms.xcms) <- pks
  }
  if (isTRUE(verbose)) {
    message(sprintf(
      "fix_xcms_chromPeaks_mz_width: %d/%d chromPeaks updated (mz width reset to %.2f ppm)",
      n_fix, n_total, ppm
    ))
  }
  xcms.xcms
}


#' @title xcmsProcessingMS1
#' @description Import `msDataFiles`, filter `ion_mode`, find peaks using `centWaveParam`, correct RT, group peaks using `peaksGroup`, fill peaks by xcms at MS1 Level
#' @param msDataFiles `char` ms file (full) paths
#' @param ion_mode to filter ion_mode, 1: positive, 0: negative, import when scans with both pos and neg
#' @param peaksGroup `vector` to xcms::PeakGroupsParam(sampleGroups), should contain "QC"
#' @param centWaveParam xcms::CentWaveParam()
#'
#' @return xcms
#' @export

xcmsProcessingMS1 <- function(xcms.xcms,
                              ion_mode = NA,
                              xcms_param = list(
                                findChromPeaks = xcms::CentWaveParam(),
                                groupChromPeaks = xcms::PeakDensityParam(sampleGroups = "A")
                              ),
                              adjustRT = T,
                              chromPeaks_fix_mz_ppm = NULL,
                              chromPeaks_max_mz_ppm = NULL,
                              BPPARAM  = BiocParallel::SnowParam(workers = 4,progressbar = T),
                              ...){



  if (is.na(ion_mode)) {
    ion_mode <- polarity(xcms.xcms )%>%unique()
    if (length(ion_mode)!=1) {
      stop("MS1 scans contain both positive and negative, please check")
    }
  }

  xcms.xcms <- ProtGenerics::filterPolarity(xcms.xcms , ion_mode)


  ### Find peaks
  message_with_time(" Find peaks...")
  xcms.xcms<-xcms::findChromPeaks(xcms.xcms,
                            param = xcms_param$findChromPeaks,
                            BPPARAM  = BPPARAM,...)

  xcms.xcms <- xcms_filter_peaks_NA(xcms.xcms)
  if (!is.null(chromPeaks_fix_mz_ppm)) {
    xcms.xcms <- fix_xcms_chromPeaks_mz_width(
      xcms.xcms,
      ppm = as.numeric(chromPeaks_fix_mz_ppm)
    )
  }
  if (!is.null(chromPeaks_max_mz_ppm)) {
    xcms.xcms <- filter_xcms_chromPeaks_mz_width(
      xcms.xcms,
      ppm = as.numeric(chromPeaks_max_mz_ppm)
    )
  }
  #xcms.xcms <- xcms_get_peak_fill(xcms.xcms)
  #mpp <- xcms::MergeNeighboringPeaksParam(expandRt = 2.5,minProp = 0.5)
  #xcms.xcms <- xcms::refineChromPeaks(xcms.xcms, mpp,
  #                                    BPPARAM  = BiocParallel::SerialParam(progressbar = T))

  ### adujust RT
  if(adjustRT){

    message_with_time(" Adjust RT...")
    peaksGroup <- Biobase::pData(xcms.xcms)$sample.type
    peak.density.param <- xcms::PeakDensityParam(sampleGroups = peaksGroup,
                                                 minFraction = 0.4,bw = 30,
                                                 binSize = 0.015)
    xcms.xcms <- xcms::groupChromPeaks(xcms.xcms,param = peak.density.param)



    if (length(Biobase::sampleNames(xcms.xcms))>1) {
      if (sum(peaksGroup=="QC") <2 ) {
        rt.adjust.param <- xcms::PeakGroupsParam(minFraction = 0.4,
                                                 #subset = which(peaksGroup == "QC"),
                                                 subsetAdjust = "previous",span = 0.4)
        xcms.xcms <- xcms::adjustRtime(xcms.xcms,param = rt.adjust.param)
      }else{
        ### adjust based on QC
        rt.adjust.param <- xcms::PeakGroupsParam(minFraction = 0.4,
                                                 subset = which(peaksGroup == "QC"),
                                                 subsetAdjust = "average",span = 0.4)
        xcms.xcms <- xcms::adjustRtime(xcms.xcms,param = rt.adjust.param)
      }
    }

  }



  ### group peaks
  message_with_time(" Group peaks...")
  peak.density.param <- xcms_param$groupChromPeaks
  peak.density.param@sampleGroups <- Biobase::pData(xcms.xcms)$sample.type
  xcms.xcms <- xcms::groupChromPeaks(xcms.xcms,param = peak.density.param)
  #xcms.xcms <- xcms_filter_feature_mz_rsd(xcms.xcms,rsd.ppm = 2)
  xcms.xcms <- xcms_get_feature_wmean(xcms.xcms)
  message_with_time(" ",nrow(xcms::featureDefinitions(xcms.xcms))," feature found")
  xcms.xcms <- xcms::fillChromPeaks(xcms.xcms,param = xcms::FillChromPeaksParam())



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






#' @describeIn xcms_extenstion plot feature intensity
#' @title plot_xcms_feature_intensity
#' @description plot feature's intensity, ordered by `Biobase::pData(xcms.xcms)$analysis.time.positive` or
#'  `Biobase::pData(xcms.xcms)$analysis.time.negative`
#'
#' @param xcms.xcms XCMSnExp object
#' @param feature_id_to_show feature id to plot
#'
#' @return ggplot object
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
  features <- xcms::featureValues(xcms.xcms)%>%
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




plot_xcms_TIC <- function(xcms.xcms,col.group = NULL,title = "TIC"){


  xcms.pdata <- Biobase::pData(xcms.xcms)
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.scan <- xcms.scan%>%
    dplyr::mutate(tic = MSnbase::tic(xcms.xcms),
                  group = xcms.pdata$group[fileIdx])%>%
    dplyr::filter(msLevel==1)

  if (is.null(col.group)) {
    col.group <- c("grey","#38C291",ggsci::pal_aaas()(10))
    names(col.group) <- c("Blank","QC",
                          setdiff(unique(xcms.scan$group),c("Blank","QC")))
  }

  ggplot(xcms.scan)+
    geom_line(aes(x = retentionTime , y = tic,
                  col = group,
                  group=fileIdx),alpha = 0.1)+
    scale_color_manual(values = col.group)+
    #scale_y_log10()+
    labs(title = title,x = "Retention Time", y = "Intensity", col = "")+
    theme_classic()->p
  p



}

plot_xcms_adjustedRT <- function(xcms.xcms){


  xcms.pdata <- Biobase::pData(xcms.xcms)%>%
    dplyr::arrange(ExpTime)%>%
    dplyr::mutate(injection_order = 1:n())
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.scan <- xcms.scan%>%
    dplyr::mutate(adrt = xcms::adjustedRtime(xcms.xcms),
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

  xcms.peaks <- xcms::chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    rownames_to_column("peak_id")%>%
    dplyr::mutate(peakWidth = rtmax-rtmin,
                  mzWidth = mzmax-mzmin,
                  mzError = mzWidth/mz*1e6)
  return(xcms.peaks)

}



get_xcms_centwave_tune <- function(xcms.xcms,
                                   iteration = 10){

  cwp <- xcms::CentWaveParam(peakwidth = c(5,20),
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

  xcms.files <- paste0(dirname(xcms.xcms),"/",Biobase::sampleNames(xcms.xcms))
  xcms.scan <- get_xcms_scan_Stat(xcms.xcms)
  xcms.sp <- Spectra(xcms.files,
                         backend = MsBackendMemory(),
                         BPPARAM = SerialParam(progressbar = T))%>%
    filterPolarity(unique(polarity(xcms.xcms)))
  spectraNames(xcms.sp) <- xcms.sp$scan_id <- xcms.scan$scan_id
  return(xcms.sp)

}

#' @importFrom ProtGenerics polarity


#' @importFrom xcms filepaths
#' @export
setMethod(f = "filepaths",
                       signature = "XCMSnExp",
                       definition = function(object) {
                         paste0(dirname(object), "/", Biobase::sampleNames(object))
                       })
#' @importFrom xcms mzrange
#' @export
setMethod(f = mzrange,
                       signature = "XCMSnExp",
                       definition = function(object) {
                         xcms.fdata <- fData(object)
                         return(c(min(xcms.fdata$scanWindowLowerLimit, na.rm = TRUE),
                                  max(xcms.fdata$scanWindowUpperLimit, na.rm = TRUE)))
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


simulate_prm <- function(xcms.fdf,
                         total.time = 28*60,
                         ms2.time = 0.265){


  t <- 0
  scan.df <- data.frame(spIdx = NULL,
                        msLevel= NULL,
                        retentionTime= NULL,
                        precursorMZ= NULL,
                        ion_id = NULL
  )
  while(t < max(xcms.fdf$rtmax)){

    ion.to.ms2 <- which(t < xcms.fdf$rtmax&t>xcms.fdf$rtmin)
    ion.to.ms2 <- na.omit(ion.to.ms2)
    ms2.scan <- data.frame(spIdx = rep(NA,length(ion.to.ms2)),
                           msLevel= rep(2,length(ion.to.ms2)),
                           retentionTime= t + seq_along(ion.to.ms2)*ms2.time,
                           precursorMZ = xcms.fdf$mzmed[ion.to.ms2],
                           ion_id = ion.to.ms2
    )
    scan.df <- rbind(scan.df,ms2.scan)
    t <- t + ifelse(ms2.time*length(ion.to.ms2)>0,ms2.time*length(ion.to.ms2),0.1)


  }



  return(scan.df)
}


#' @title Compute MS1 purity matrix for xcms features
#' @description
#' Calculate a feature-by-sample MS1 purity matrix by extracting, for each feature
#' in \code{xcms.xcms}, the closest MS1 scan (by retention time) from \code{xcms.ms1.sp}
#' in each sample file (matched by \code{Spectra::dataOrigin}). Purity is calculated
#' by \code{get_spectra_ion_purity()} within an isolation window around the feature m/z.
#'
#' Note: this function requires MS1 \code{Spectra}. It does **not** import spectra
#' from \code{xcms.xcms} automatically. Use \code{object@spectra$MS1_Spectra} from
#' \code{MSdev_extract_Spectra()} (recommended) or build MS1 Spectra yourself.
#'
#' @param xcms.xcms \code{XCMSnExp} with grouped features (must have \code{featureDefinitions}).
#' @param xcms.ms1.sp MS1 \code{Spectra} object covering the same files as \code{xcms.xcms}.
#' @param ppm numeric, ppm tolerance for m/z window.
#' @param isolation_half_window numeric, half isolation window (m/z).
#'
#' @return numeric matrix with rows = \code{feature_id}, columns = sample files (by \code{dataOrigin}).
#' @export
get_xcms_feature_purity_matrix <- function(xcms.xcms,
                                           xcms.ms1.sp,
                                           ppm = 10,
                                           isolation_half_window = 0.2){

  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)

  if (missing(xcms.ms1.sp) || is.null(xcms.ms1.sp)) {
    stop("xcms.ms1.sp (MS1 Spectra) must be provided.")
  }

  ### calc ms1_purity by ms1.sp
  {
    sp.rt <- rtime(xcms.ms1.sp)
    sp.origin <- xcms.ms1.sp$dataOrigin
    sp.idx.split <- split(seq_along(xcms.ms1.sp),sp.origin)
    f.sp.idx <- lapply(seq_len(nrow(xcms.fdf)),function(x){
      rt.diff <- abs(sp.rt-xcms.fdf$rtmed[x])
      rt.diff <- split(rt.diff,f = sp.origin)
      rt.idx <- unname(sapply(rt.diff,which.min))
      sapply(seq_along(rt.idx),function(i){
        sp.idx.split[[i]][rt.idx[i] ]
      })
    })
    f.sp <- lapply(f.sp.idx,function(x){xcms.ms1.sp[x]})
    message_with_time("calculating MS1 purity...")
    ms1_purity <- bplapply(seq_len(nrow(xcms.fdf)),function(x){
      get_spectra_ion_purity(f.sp[[x]],xcms.fdf$mzmed[x],ppm,isolation_half_window)
    },BPPARAM = SerialParam(progressbar = T))
    ms1_purity_matrix <- do.call(rbind,ms1_purity)
    rownames(ms1_purity_matrix) <-xcms.fdf$feature_id
  }

  return(ms1_purity_matrix)

}



xcms_get_feature_purity <- function(xcms.xcms,
                                    xcms.ms1.sp,
                                    ms1_purity_matrix = NULL,
                                    split_source = F,
                                    selected.sample =which( pData(xcms.xcms)$sample.type !="Blank"),
                                    FUN = max,
                                    ppm = 10 ,
                                    isolation_half_window = 0.2
                                    ){

  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)

  ### calc ms1_purity_matrix
  {
    if (is.null(ms1_purity_matrix)) {
      ms1_purity_matrix <- get_xcms_feature_purity_matrix(xcms.xcms,
                                                          xcms.ms1.sp = xcms.ms1.sp,
                                                          ppm = ppm,
                                                          isolation_half_window = isolation_half_window)
    }
  }

  ### aggregate purity
  {

    xcms.pdata <- pData(xcms.xcms)[selected.sample,]
    xcms.pm <- ms1_purity_matrix[,xcms.pdata$sampleNames,drop = F]
    if (split_source) {
      xcms.pm <- apply(xcms.pm,1,mean_f,xcms.pdata$sample.source,simplify =F)%>%
        do.call(rbind,.)
    }

    ms1_purity <- apply(xcms.pm,1,FUN,na.rm =T)

  }

  xcms.fdf$ms1_purity <- ms1_purity
  xcms.fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)
}



cbind_Chromatograms <- function(...){

  chrom.list <- list(...)
  if (length(chrom.list) == 1) {
    return(chrom.list[[1]])
  }
  chrom.featureDefinitions <- sapply(chrom.list,
                                     function(x)x@featureDefinitions)[[1]]
  chrom.phenoData <- sapply(chrom.list,
                                     function(x)x@phenoData)%>%
    do.call(Biobase::combine,.)
  chrom.featureData <- sapply(chrom.list,
                                     function(x)x@featureData)%>%
    do.call(Biobase::combine,.)
  chrom.Data<- sapply(chrom.list,
                              function(x)x@.Data)
  chrom.processHistory<- sapply(chrom.list,
                              function(x)x@.processHistory)[1]
  xchrom <-xcms::XChromatograms(chrom.Data,
                     ncol = ncol(chrom.Data),
                     phenoData = chrom.phenoData)
  xchrom@featureDefinitions <- chrom.featureDefinitions
  xchrom@featureData <- chrom.featureData
  xchrom@.processHistory <- chrom.processHistory
  return(xchrom)
}



get_xcms_quantify_MSIP <- function(xcms.xcms){

  xcms::quantify(xcms.xcms,missing = 1,method="max",value = "into")

}




xcms_from_ms2_spectra <- function(sp.ms2 ,
                                  sample.info,
                                  ppm = 10,
                                  peak_width = 60){


  ### assign
  {
    sample.info <- sample.info%>%
      dplyr::mutate(sample = as.numeric(factor(msData.files)))

    sp.peaks.df <- data.frame(
      mz = precursorMz(sp.ms2),
      rt = rtime(sp.ms2),
      sample.files = dataOrigin(sp.ms2),
      into = sp.ms2$totIonCurrent
    )%>%
      dplyr::mutate(sample = match_path(sample.files,
                                        sample.info$msData.files),
                    sample = sample.info$sample[sample])

    sp.peaks.matrix <- sp.peaks.df%>%
      dplyr::mutate(mzmin = mz,
                    mzmax = mz,
                    rtmin = rt-peak_width/2,
                    rtmax = rt+peak_width/2)%>%
      dplyr::select(any_of(c("mz","mzmin","mzmax",
                             "rt","rtmin","rtmax",
                             "into","intb","maxo",
                             "sn","sample")))%>%
      as.matrix()
    sp.peaks.data <- sp.peaks.df%>%
      dplyr::mutate(ms_level = 1,
                    ms_level = as.integer(ms_level),
                    is_filled = F)%>%
      S4Vectors::DataFrame()

    ion_df <- do_groupChromPeaks_density(sp.peaks.df,
                                         bw = peak_width,
                                         sampleGroups = sample.info$sample.source,
                                         binSize = 0.001,
                                         ppm = ppm)

    ion_table <-ion_df %>%
      dplyr::mutate(feature_id = paste0("FTS",num2str(1:n())),
                    .before = mzmed)

  }


  ### simulate xcms class
  {

    MsFeatureData <- new("MsFeatureData",
                         chromPeaks = sp.peaks.matrix,
                         chromPeakData = sp.peaks.data,
                         featureDefinitions =  S4Vectors::DataFrame(ion_table))

    XCMSnExp <- new("XCMSnExp")
    XCMSnExp@msFeatureData <- MsFeatureData
  }


  return(XCMSnExp)

}


xcms_get_feature_adduct_connection <- function(xcms.xcms,rt.tol = 5,ppm = 10){


  pol <- unique(polarity(xcms.xcms))
  adduct.diff <- get_adduct_mass_diff(pol)
  adduct.diff <- adduct.diff[order(adduct.diff$mass_diff),]
  #xcms.xcms <- xcms_get_feature_group(xcms.xcms,diffRt = 10,intCor = NULL,eicCor = NULL)
  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)


  ### Construct connection
  {




    # Generate all connection
    # filter rt diff
    # calc mz diff
    # pre-filter mz.diff
    # calc mz.mean
    {

      xcms.net <- expand.grid(
        from = 1:nrow(xcms.fdf),
        to = 1:nrow(xcms.fdf)
      )
      xcms.net <- data.table::as.data.table(xcms.net)
      xcms.net <- xcms.net[from != to ][
        , rt.diff := abs(xcms.fdf$rtmed[to]-xcms.fdf$rtmed[from]) ][
        rt.diff < rt.tol,][
          , c("mz.diff") := .( xcms.fdf$mzmed[to] - xcms.fdf$mzmed[from])][
            mz.diff > min(adduct.diff$mass_diff)&mz.diff < max(adduct.diff$mass_diff) ][
              ,mz.mean := xcms.fdf$mzmed[to] + xcms.fdf$mzmed[from]   ]
    }


    # match mz.diff to adduct.diff
    match.df <- match_mz_foverlaps(mz1 = xcms.net$mz.diff,mz2 = adduct.diff$mass_diff,
                                   ppm.base = xcms.net$mz.mean,ppm = ppm)


    # add adduct.diff data
    xcms.net.matched <- xcms.net[match.df$ion1,][
      ,c("adduct.diff.idx","mz.ppm"):= .(match.df$ion2,match.df$mz.ppm)
    ][mz.ppm  < ppm,][ ### connect within ppm
      ,c("adduct.mass.diff","from.adduct","to.adduct"):= .(adduct.diff$mass_diff[adduct.diff.idx],
                                              adduct.diff$adduct.from[adduct.diff.idx],
                                              adduct.diff$adduct.to[adduct.diff.idx])
    ]


  }


  ###
  {

    xcms.net.matched$label  <- paste0(xcms.net.matched$from.adduct," to ",xcms.net.matched$to.adduct)
    xcms.fdf.ig <- igraph::graph_from_data_frame(xcms.net.matched)
    visIgraph(igraph_filter_distance(xcms.fdf.ig,from = "11959",1))%>%
      visEdges(smooth = T)


    node.group <- igraph::components(xcms.fdf.ig)$membership

    ig <- igraph_filter_vertex(xcms.fdf.ig , which(node.group==1))

    eda <- edata(ig)
    vda <- vdata(ig)
    plot_density(eda$rt.diff)

  }

}


plotly_xcms_feature_group <- function(xcms.xcms){


  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
    as.data.frame()
  #ggplot(xcms.fdf)+
  #  geom_point(aes(x = rtmed , y = mzmed, col = feature_group))+
  #  theme(legend.position = "none")

  plotly::plot_ly(xcms.fdf)%>%
    add_markers(x = ~rtmed, y = ~mzmed, color = ~feature_group)

}


#' @title Update feature mz/rt using peak-intensity weighted means
#' @description Recomputes `mzmed` and `rtmed` in `xcms::featureDefinitions(xcms.xcms)`
#' using peak-level `mz`/`rt` weighted by peak `maxo` (maximum intensity) from
#' `xcms::chromPeaks(xcms.xcms)`.
#'
#' @param xcms.xcms An `xcms::XCMSnExp` object with feature definitions and chromPeaks.
#'
#' @return An updated `xcms::XCMSnExp` object where `featureDefinitions(object)$mzmed`
#' and `featureDefinitions(object)$rtmed` are replaced by the intensity-weighted means
#' across each feature's constituent peaks.
#'
#' @export
xcms_get_feature_wmean <- function(xcms.xcms){

  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  xcms.pks <- xcms::chromPeaks(xcms.xcms)

  wrt <- sapply(xcms.fdf$peakidx,function(x){

    x.rt <- xcms.pks[x,"rt"]
    x.int <- xcms.pks[ x, "maxo"]
    weighted.mean(x.rt,x.int)
  })

  wmz <- sapply(xcms.fdf$peakidx,function(x){

    x.mz <- xcms.pks[x,"mz"]
    x.int <- xcms.pks[ x, "maxo"]
    weighted.mean(x.mz,x.int)
  })

  xcms.fdf$mzmed <- wmz
  xcms.fdf$rtmed <- wrt
  xcms.fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)
}


xcms_filter_feature_mz_rsd <- function(xcms.xcms, rsd.ppm = 2){


  fdf <- xcms::featureDefinitions(xcms.xcms)
  ch <- xcms::chromPeaks(xcms.xcms)
  mz.sd <- sapply(fdf$peakidx,function(x){
    sd(ch[x,'mz'])/mean(ch[x,'mz']) * 1e6
  })
  #plot_density(mz.sd)
  fdf <- fdf[mz.sd < rsd.ppm,]
  fdf$feature_id <- paste0("FT",num2str(1:nrow(fdf)))
  rownames(fdf) <- fdf$feature_id
  fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)

}
xcms_filter_feature_rt_rsd <- function(xcms.xcms, rt.shift = 5 ){


  fdf <- xcms::featureDefinitions(xcms.xcms)
  ch <- xcms::chromPeaks(xcms.xcms)
  rt.sd <- sapply(fdf$peakidx,function(x){
    sd(ch[x,'rt'])/mean(ch[x,'rt'])
  })
  #plot_density(mz.sd)
  fdf <- fdf[rt.sd < rt.shift,]
  fdf$feature_id <- paste0("FT",num2str(1:nrow(fdf)))
  rownames(fdf) <- fdf$feature_id
  fdf -> xcms::featureDefinitions(xcms.xcms)
  return(xcms.xcms)

}
