annotateSpectra <- function(expSpec,refSpec){

  to.return <- list(
    mz = median(ProtGenerics::precursorMz(expSpec)),
    rt = median(ProtGenerics::rtime(expSpec)),
    ref.mz = NA,
    ref.rt = NA,
    score = 0,
    compound = NA,
    adduct = NA,
    inchikey = NA,
    kegg.id = NA,
    origin = NA,
    "expSpec" = expSpec
  )
  if (length(refSpec)==0) {
    return(to.return)
  }

  mz.error <- abs(matrixSub(expSpec$precursorMz , refSpec$precursorMz))/to.return$mz*1e6
  mz.score <- 1 - mz.error/20 ###ppm 20

  rt.error <- abs(matrixSub(expSpec$rtime , refSpec$rtime))
  rt.score <- 2- rt.error/20
  rt.score[is.na(rt.score)] <- 0

  sp.score <- Spectra::compareSpectra(expSpec,refSpec)
  sp.score[is.nan(sp.score)] <- 0

  score <- mz.score*0.2 + rt.score*0.2 + sp.score*0.6


    refSpecMatched <- refSpec[ceiling(which.max(score)/length(expSpec))]
    score <- sp.score[which.max(score)]
    to.return$ref.mz <- refSpecMatched$precursorMz
    to.return$ref.rt <- refSpecMatched$rtime
    to.return$score <- score
    to.return$compound <- refSpecMatched$name
    to.return$adduct <- refSpecMatched$adduct
    to.return$inchikey <- refSpecMatched$inchikey
    to.return$kegg.id <- refSpecMatched$kegg.id
    to.return$origin <- refSpecMatched$database
    to.return$refSpec <- refSpecMatched
  return(to.return)




}

annotateSpectraMSdb <- function(expSpec,refSpec){

  to.return <- list(
    mz = median(ProtGenerics::precursorMz(expSpec)),
    rt = median(ProtGenerics::rtime(expSpec)),
    ref.mz = NA,
    ref.rt = NA,
    score = 0,
    MSDB_id =NA,
    "expSpec" = expSpec,
    "refSpec" = refSpec
  )
  if (length(refSpec)==0) {
    return(to.return)
  }

  mz.error <- abs(matrixSub(expSpec$precursorMz , refSpec$precursorMz))/to.return$mz*1e6
  mz.score <- 1 - mz.error/20 ###ppm 20

  rt.error <- abs(matrixSub(expSpec$rtime , refSpec$rtime))
  rt.score <- 2- rt.error/20
  rt.score[is.na(rt.score)] <- 0

  expSpec.sub <-Spectra::filterMzRange(expSpec,c(0,expSpec$precursorMz-0.5))
  refSpec.sub <-Spectra::filterMzRange(refSpec,c(0,expSpec$precursorMz-0.5))
  sp.score <- Spectra::compareSpectra(expSpec.sub,refSpec.sub)
  sp.score[is.nan(sp.score)] <- 0

  score <- mz.score*0.2 + rt.score*0.2 + sp.score*0.6


  refSpecMatched <- refSpec[ceiling(which.max(score)/length(expSpec))]
  score <- sp.score[which.max(score)]
  to.return$ref.mz <- refSpecMatched$precursorMz
  to.return$ref.rt <- refSpecMatched$rtime
  to.return$score <- score
  to.return$MSDB_id <-refSpecMatched$MSDB_id
  to.return$refSpec <- refSpecMatched
  return(to.return)




}



makeSpectra <- function(precursorMz = 0,
                        rtime = 0 ,...){

  Spectra::Spectra(S4Vectors::DataFrame(precursorMz = precursorMz,
                             rtime = rtime,
                             ...),backend = Spectra::MsBackendMemory())

}

makeEmptySpectra <- function(...){
  Spectra::Spectra(S4Vectors::DataFrame(precursorMz = 0,
                                        rtime = 0,
                                        ...),backend = Spectra::MsBackendMemory())

}


filterSpectraIntensity <- function(sp,ratio = 0.05){

  Spectra::filterIntensity( sp, intensity = function(z) z/max(z) > ratio )

}


filterSpectra_below_PrecursorMz <- function(sp){

  if (!length(sp)) return(sp)
  nf <-  function(z,  precursorMz ,...) {
    precursorMz <- ifelse(is.na(precursorMz),Inf,precursorMz)
    idx <- z[, "mz"] < precursorMz-0.5
    z[idx,,drop = F]
  }
  sp <- Spectra::addProcessing(sp,FUN = nf, spectraVariables = "precursorMz")%>%
    Spectra::applyProcessing(BPPARAM = SerialParam())

  return(sp)
}



normalizeSpectra <- function(sp,norm_to = "tic"){

  if (norm_to == "tic") {
    nf <-  function(z, ...) {
      #z[,"maxIntensity"] <- max( z[, "intensity"] )
      z[, "intensity"] <- z[, "intensity"] /
        sum(z[, "intensity"], na.rm = TRUE) * 100
      z
    }
  }
  if (norm_to == "max") {
    nf <-  function(z, ...) {
      #z[,"maxIntensity"] <- max( z[, "intensity"] )
      z[, "intensity"] <- z[, "intensity"] /
        max(z[, "intensity"], na.rm = TRUE) * 100
      z
    }
  }

  sp <- Spectra::addProcessing(sp,nf)

  return(sp)

}

normalizeSpectra_by_precursorIntensity <- function(sp){

  if (!length(sp)) return(sp)
  if(all(sp$precursorIntensity==0)){
    sp$precursorIntensity[sp$precursorIntensity==0] <- sp$totIonCurrent[sp$precursorIntensity==0]
    warning("No precursor intensity, replace with TIC")
  }
  nf <-  function(z,  precursorIntensity ,...) {
    #z[,"maxIntensity"] <- max( z[, "intensity"] )
    z[, "intensity"] <- z[, "intensity"] / precursorIntensity * 100
    z
  }
  sp <- Spectra::addProcessing(sp,FUN = nf, spectraVariables = "precursorIntensity")%>%
    Spectra::applyProcessing(BPPARAM = SerialParam())

  return(sp)

}

#' get_Spectra_data
#'
#' @title Get Spectra Data
#' @description Extracts data from a Spectra object into a data.frame with columns for spectrum ID, mz, intensity, and requested spectrum variables.
#' @param sp A `Spectra` object containing mass spectrometry data.
#' @param var Character vector of spectrum variables to include in the output.
#'   Default is `c("precursorMz", "collisionEnergy")`. Should be valid names
#'   from `Spectra::spectraVariables()`.
#'
#' @return A data.frame with columns: `sp.id` (spectrum identifier), `mz` (mass-to-charge ratio),
#'   `intensity` (peak intensity), and any additional variables specified in `var`.
#' @export
#'

get_Spectra_data <- function(sp,var = c("precursorMz","collisionEnergy")){


  if (!length(sp)) {
    sp.data <- data.frame(matrix(nrow = 0,ncol = 3+length(var)))%>%
      `colnames<-`(c("sp.id","mz","intensity",var))
      return(sp.data)
  }
  sp.data <- Spectra::spectraData(sp)%>%
    as.data.frame()%>%
    dplyr::mutate(id = 1:n())

  sp.peaks.data <- peaksData(sp)%>%
    unname%>%
    #lapply(as.matrix,ncol = 2)%>%
    lapply(data.frame)%>%
    data.table::rbindlist(use.names = F,
                          idcol = "sp.id")

  spec.df <- sp.peaks.data %>%
    dplyr::mutate(
      sp.data[match(sp.id,sp.data$id)   ,var,drop = F]
    )
  if (!is.null(spectraNames(sp)))
    spec.df$sp.id <- spectraNames(sp)[spec.df$sp.id]

  return(spec.df)
}

#' combineSpectra_groupby_ce
#'
#' @title Combine Spectra by Collision Energy
#' @description Combines spectra that share the same collision energy using a TIC-weighted
#'   peak combination method. This is useful for merging technical replicates or
#'   averaged spectra at each collision energy level.
#' @param sp A `Spectra` object to combine.
#' @param minProp Numeric between 0 and 1. Minimum proportion of spectra that must
#'   contain a peak for it to be retained in the combined spectrum. Default is `0.5`.
#' @param ppm Numeric. Parts per million tolerance for grouping peaks. Default is `10`.
#' @param plot Logical. If `TRUE`, generates plots of the combination results. Default is `FALSE`.
#' @param ... Additional arguments passed to [Spectra::combineSpectra()].
#'
#' @return A `Spectra` object with combined spectra for each unique collision energy.
#' @export
#'
combineSpectra_groupby_ce <- function(sp,
                                      minProp = 0.5,
                                      ppm = 10,
                                      plot = F,
                                      ...){

  sp.ce <- Spectra::combineSpectra(sp,
                                   peaks = "intersect",
                                   intensityFun = mean,
                                   FUN = combinePeaksData_tic_weighted,
                 f = sp$collisionEnergy,
                 weighted =T,
                 minProp=minProp,
                 ppm = ppm,
                 ...)
  if(plot){
    for (ce  in unique(collisionEnergy(sp))) {
      this.sp <- sp[collisionEnergy(sp)==ce]
      this.sp.list <- split(this.sp,f = this.sp$sp_id)
      this.p.list <- lapply(this.sp.list,plotSpec)
      p <- ggplot_sum_patchwork(this.p.list)+
        plot_layout(ncol = 2)+
        plot_annotation(title = paste0("collisionEnergy = ",ce))
      open_plot_win(p,10,ceiling(length(this.p.list)/2)*3)
    }


  }

  return(sp.ce)

}


combinePeaksData_tic_weighted <-
  function (x, intensityFun = base::mean, mzFun = base::mean,
          weighted = FALSE, tolerance = 0, ppm = 0, timeDomain = FALSE,
          peaks = c("union", "intersect"), main = integer(), minProp = 0.5,
          ...)
{
  peaks <- match.arg(peaks)
  lenx <- length(x)
  if (lenx == 1)
    return(x[[1]])
  mzs <- lapply(x, "[", , y = 1)
  mzs_lens <- lengths(mzs)
  mzs <- unlist(mzs, use.names = FALSE)
  mz_order <- order(mzs)
  mzs <- mzs[mz_order]
  if (timeDomain)
    mz_groups <- MsCoreUtils::group(sqrt(mzs), tolerance = tolerance,
                       ppm = ppm)
  else mz_groups <- MsCoreUtils::group(mzs, tolerance = tolerance, ppm = ppm)
  ints <- unlist(lapply(x, "[", , y = 2), use.names = FALSE)[mz_order]
  tics <- rep(sapply(lapply(x, "[", , y = 2),sum),mzs_lens)[mz_order]
  if (length(main)) {
    if (main < 1 || main > lenx)
      stop("'main' has to be larger than 1 and smaller than ",
           lenx)
    is_in_main <- rep.int(seq_along(mzs_lens), mzs_lens)[mz_order] ==
      main
    keep <- mz_groups %in% mz_groups[is_in_main]
    mz_groups <- mz_groups[keep]
    mzs <- mzs[keep]
    ints <- ints[keep]
  }
  mzs <- split(mzs, mz_groups)
  ints <- split(ints, mz_groups)
  tics <- split(tics,mz_groups)
  if (peaks == "intersect") {
    keep <- lengths(mzs) >= (lenx * minProp)
    if (any(keep)) {
      mzs <- mzs[keep]
      ints <- ints[keep]
      tics <- tics[keep]
    }
    else return(cbind(mz = numeric(), intensity = numeric()))
  }
  if (weighted) {
    wm <- stats::weighted.mean
    mzs <- mapply(mzs, ints, FUN = function(mz, w) wm(mz,
                                                      w + 1, na.rm = TRUE, USE.NAMES = FALSE))
  }
  else mzs <- vapply1d(mzs, FUN = mzFun)
  if (is.unsorted(mzs))
    stop("m/z values of combined spectrum are not ordered increasingly")
  ints <- mapply(ints, tics, FUN = function(int, w) stats::weighted.mean(int,
                                                    w , na.rm = TRUE, USE.NAMES = FALSE))
 # vapply1d(ints, FUN = intensityFun)
  cbind(mz = mzs, intensity = ints)
}


Spectra_fill_3CE <- function(sp){

  ce.to.fill <- setdiff(c(10,20,40),
                        collisionEnergy(sp))
  sp.to.fill <- sp[1]
  remove_peaks <- function(z){
    z[1,] <- c(0,0)
    z[1,,drop = F]
  }
  sp.to.fill <- addProcessing(sp,remove_peaks)%>%
    Spectra::applyProcessing(BPPARAM = SerialParam())
  sp.to.fill <- sp.to.fill[rep(1,length(ce.to.fill))]
  collisionEnergy(sp.to.fill) <- ce.to.fill
  sp.filled <- c(sp,sp.to.fill)
  return(sp.filled)
}

get_Spectra_ms2_feature_id <- function(sp,
                                    featuredef,
                                    ppm = 5,
                                    rt.tol = 5){



  sp.data <- Spectra::spectraData(sp)%>%
    as.data.frame()%>%
    dplyr::mutate(precursorMZ = precursorMz,
                  retentionTime = rtime)%>%
    dplyr::filter(msLevel == 2)

  match.df <- match_mz_rt(featuredef$mzmed,
                          featuredef$rtmed,
                          sp.data$precursorMZ,
                          sp.data$retentionTime,
                          mz.ppm = ppm,
                          rt.tol = rt.tol)
  match.df <- match.df%>%
    dplyr::mutate(featuredef[ion1,],
                  sp_id = sp.data$sp_id[ion2],
                  sp_rt = sp.data$retentionTime[ion2])%>%
    dplyr::group_by(sp_id)%>%
    dplyr::slice_min(mz.error)%>%
    dplyr::ungroup()



  sp.data$feature_id<- sapply(
    1:nrow(sp.data),
    function(i){
      fid <- match.df$feature_id[match.df$ion2==i]
      if (length(fid)==0) {
        return(NA)
      }
      return(fid)
    }
  )

  return(sp.data)



}

get_Spectra_transition <- function(sp){

  sp <- filterSpectra_below_PrecursorMz(sp)
  sp.data <- get_Spectra_data(sp)%>%
    dplyr::mutate(groupMz(mz,return.type = "data.frame"))%>%
    dplyr::group_by(mz.group )%>%
    dplyr::mutate(mz.group.max = max(intensity),
                  mz.group.n = n())%>%
    dplyr::slice_max(intensity , n = 1)%>%
    dplyr::ungroup()%>%
    dplyr::slice_max(intensity , n = 2)%>%
    dplyr::select( precursorMz , productMz =  mz,
                   collisionEnergy , productRatio = intensity)

  return(sp.data)

}


plot_Spectra_Mirror<- function(sp1,sp2,show.label = "rtime"){

  sp.data <-rbind(get_Spectra_data(sp1)%>%
                    dplyr::mutate(sp.id=1),
                  get_Spectra_data(sp2)%>%
                    dplyr::mutate(sp.id=2))%>%
    dplyr::mutate(groupMz(mz))%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = case_when(sp.id==1~intensity,
                                   T~-intensity))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate(matched = n()==2)

  label.df <- dplyr::filter(sp.data , matched)
  ymax.abs <- max(abs(sp.data$intensity))
  ggplot(sp.data)+
    geom_hline(yintercept = 0 , size = 0.2,col = "grey")+
    geom_segment(aes(x = x,y =y,xend = xend,
                     yend = yend,col = matched),
                 size = 0.2,
                 show.legend = F)+
    geom_point(aes(x = x, y = yend , alpha = matched,col = matched),
               show.legend = F,size = 0.5)+
    geom_text(aes(x = quantile(range(sp.data$x),0.8),
                  y =ymax.abs*0.8,
                  label = paste0(show.label," = \n",sp1[[show.label]])),
              size = 2,
              check_overlap = T)+
    geom_text(aes(x = quantile(range(sp.data$x),0.8),
                  y = -ymax.abs*0.8,
                  label = paste0(show.label," = \n",sp2[[show.label]])),
              size = 2,
              check_overlap = T)+
    ggrepel::geom_text_repel(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                             #segment.size = 0.1,
                             data = label.df)+
    scale_color_manual(values = c(`FALSE` = "grey",`TRUE` = "#80B1D3"))+
    scale_alpha_manual(values = c(`FALSE` = 0,`TRUE` =1))+
    ylim(c(-ymax.abs,ymax.abs))+
    labs(x = "Mz",y = "Intensity")+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))->p
  p


}





#' plot_Spectra
#'
#' @title Plot Mass Spectra
#' @description Creates a ggplot2 bar plot visualization of mass spectra data.
#'   Peaks above 10% of maximum intensity are highlighted in blue.
#' @param sp A `Spectra` object to plot.
#' @param label.top Integer. Number of highest intensity peaks to label with m/z values.
#'   Default is `10`.
#'
#' @return A `ggplot` object displaying the mass spectrum.
#' @export
plot_Spectra<- function(sp,label.top = 10){

  sp.data <-get_Spectra_data(sp)%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = intensity,
                  int.rank = rank(-intensity) ,
                  matched = intensity > max(intensity,na.rm = T)/10 )

  label.df <- dplyr::filter(sp.data , int.rank  < label.top)
  ymax.abs <- max(abs(sp.data$intensity))
  ggplot(sp.data)+
    geom_hline(yintercept = 0 , linewidth = 0.2,col = "grey")+
    geom_segment(aes(x = x,y =y,xend = xend,
                     yend = yend,col = matched),
                 linewidth = 0.2,
                 show.legend = F)+
    geom_point(aes(x = x, y = yend , alpha = matched,col = matched),
               show.legend = F,size = 0.5)+
    ggrepel::geom_text_repel(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                            # segment.size = 0.1,
                             data = label.df)+
    scale_color_manual(values = c(`FALSE` = "grey",`TRUE` = "#80B1D3"))+
    scale_alpha_manual(values = c(`FALSE` = 0,`TRUE` =1))+
    scale_y_continuous(expand = expansion(0,0),
                       limits = c(0,ymax.abs*1.1),
                       labels = scales::scientific)+
    labs(x = "Mz",y = "Intensity")+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))->p
  p


}


plot_Spectra_quality <- function(sp){

  plot_Spectra(sp)+
    geom_hline(yintercept = sp$noise,col = "grey")+
    labs(subtitle  = paste0("Pre.int: ",format(Spectra::precursorIntensity(sp),digit = 3,sci = T) ,"\n",
                           "TIC: ", format(sp$totIonCurrent,digit = 3,sci = T) ,"\n",
                           "Base.int: ", format(sp$basePeakIntensity,digit = 3,sci = T),"\n",
                           "Noise: ", format(sp$noise,digit = 3,sci = T),"\n",
                           "SNR: ",format(sp$snr,digit = 3,sci = T)))

}

plot_Spectra_CE<-function(sp){


  sp.data <- get_Spectra_data(sp)%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = intensity,
                  int.rank = rank(-intensity) ,
                  matched = intensity > max(intensity,na.rm = T)/10 )%>%
    dplyr::group_by(collisionEnergy)%>%
    dplyr::mutate(step = cur_group_id()-1)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(xstep = step * max(x,na.rm = T)/100*5,
                  ystep = step * max(yend,na.rm = T)/100*5,
                  x = x+xstep,
                  xend = xend + xstep,
                  y = y+ystep,
                  yend = yend + ystep,
                  collisionEnergy= factor(collisionEnergy))%>%
    dplyr::mutate(groupMz(mz,return.type = "data.frame"))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate( highlight = length(unique(collisionEnergy)) <
                                        length(levels(collisionEnergy))*0.5,
                   hx = min(x),
                   hxend = max(x),
                   hy = min(y),
                   hyend =max(y))

  col.list <- randomcoloR::distinctColorPalette(length(unique(sp.data$step)))
  label.df <- dplyr::filter(sp.data , int.rank  < 5)
  ymax.abs <- max(abs(sp.data$yend))
  ggplot(sp.data)+
    geom_segment(aes(x = hx,y =hy,xend = hxend,
                     yend = hyend) ,col = "grey",linewidth =1)+
    geom_hline(aes( yintercept= ystep , col = collisionEnergy))+
    geom_segment(aes(x = x,y =y,xend = xend,
                     yend = yend,col = collisionEnergy),
                 linewidth = 0.5,alpha = 0.7,
                 show.legend = F)+
    geom_point(aes(x = x, y = yend ,col = collisionEnergy),
               show.legend = F,size = 0.5)+
    geom_text(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                             #segment.size = 0.1,
                             data = label.df)+
    scale_color_manual(values = col.list)+
    scale_y_continuous(expand = expansion(0,0),lim = c(0,ymax.abs*1.05))+
    labs(x = "Mz",y = "Intensity")+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))->p
  p




}


plot_Spectra_product_CE_curve <- function(sp){

  sp.data <- get_Spectra_data(sp)%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = intensity,
                  int.rank = rank(-intensity) ,
                  matched = intensity > max(intensity,na.rm = T)/10 )%>%
    dplyr::group_by(collisionEnergy)%>%
    dplyr::mutate(step = cur_group_id()-1)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(xstep = step * max(x,na.rm = T)/100*5,
                  ystep = step * max(yend,na.rm = T)/100*5,
                  x = x+xstep,
                  xend = xend + xstep,
                  y = y+ystep,
                  yend = yend + ystep,
                  cex = as.numeric(collisionEnergy),
                  collisionEnergy= factor(collisionEnergy))%>%
    dplyr::mutate(groupMz(mz,return.type = "data.frame"))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate( highlight = length(unique(collisionEnergy)) >
                     length(levels(collisionEnergy))*0.5,
                   hx = min(x),
                   hxend = max(x),
                   hy = min(y),
                   hyend =max(y),
                   mz.center = sprintf("%.4f",mz.center)
                   )%>%
    dplyr::ungroup()%>%
    dplyr::filter(highlight)

  ggplot(sp.data,aes(x = cex , y = intensity ,col = mz.center))+
    geom_point()+
    geom_smooth(formula = y~x,method = "loess")+
    scale_color_random()+
    #scale_y_log10()+
    labs(x = "Collision Energy",y = "Intensity normalized to Precursor",
         col = "Product mz")+
    guides(col = guide_legend(ncol = 3))+
    theme_bw()+
    theme(legend.position = "right")->p
  p


}


plot_Spectra_RT<-function(sp){

  sp.data <- get_Spectra_data(sp,var =  c("precursorMz","rtime"))%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = intensity,
                  int.rank = rank(-intensity) ,
                  matched = intensity > max(intensity,na.rm = T)/10 )%>%
    dplyr::group_by(rtime)%>%
    dplyr::mutate(step = cur_group_id()-1)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(xstep = step * max(x,na.rm = T)/100*5,
                  ystep = step * max(yend,na.rm = T)/100*5,
                  x = x+xstep,
                  xend = xend + xstep,
                  y = y+ystep,
                  yend = yend + ystep,
                  rtime= factor(rtime))%>%
    dplyr::mutate(groupMz(mz,return.type = "data.frame"))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate( highlight = length(unique(rtime)) <
                     length(levels(rtime))*0.5,
                   hx = min(x),
                   hxend = max(x),
                   hy = min(y),
                   hyend =max(y))

  col.list <- randomcoloR::distinctColorPalette(length(unique(sp.data$step)))
  label.df <- dplyr::filter(sp.data , int.rank  < 5)
  ymax.abs <- max(abs(sp.data$yend))
  ggplot(sp.data)+
    geom_segment(aes(x = hx,y =hy,xend = hxend,
                     yend = hyend) ,col = "grey",linewidth =1)+
    geom_hline(aes( yintercept= ystep , col = rtime))+
    geom_segment(aes(x = x,y =y,xend = xend,
                     yend = yend,col = rtime),
                 linewidth = 0.5,alpha = 0.7,
                 show.legend = F)+
    geom_point(aes(x = x, y = yend ,col = rtime),
               show.legend = F,size = 0.5)+
    ggrepel::geom_text_repel(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                             #segment.size = 0.1,
                             data = label.df)+
    scale_color_manual(values = col.list)+
    scale_y_continuous(expand = expansion(0,0),lim = c(0,ymax.abs*1.05))+
    labs(x = "Mz",y = "Intensity",col = "Retention Time")+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))->p
  p




}



plot_Spectra_similarity <- function(sp,col_title = NULL){

  if(length(sp) > 100){
    sp <- sp[sample(seq_along(sp),100)]
    warning("Too many Spectra")
  }

  sp.data <- spectraData(sp)%>%as.data.frame()%>%
    dplyr::arrange(collisionEnergy,rtime)
  sp <- sp[rownames(sp.data)]
  sp.simlilarity <- Spectra::compareSpectra(sp)
  colnames(sp.simlilarity) <- sp.data$rtime%>%round(0)
  Heatmap(sp.simlilarity,
          name = "Spectra\nSimilarity",
          col = circlize::colorRamp2(
            breaks = c(0,0.25,0.5,0.75,1),
            colors = c("#395586","#0094B2",
                       "#FFE8E9","#FF737B",
                       "#E32237")
                          ),
          top_annotation = HeatmapAnnotation(
            "Precursor\nIntensity" = anno_lines(sp.data$precursorIntensity,
                                   #pt_gp = gpar(col=sp.data$collisionEnergy),
                                   add_points = T),
            "Collision\nEnergy" = anno_barplot(sp.data$collisionEnergy,
                              gp = gpar(fill =sp.data$collisionEnergy ),
                              bar_width = 1),
            annotation_height = unit(1.5,"inch"),
            show_annotation_name = T),
          column_title = col_title,
          cluster_columns = F,
          cluster_rows = F,
          show_row_names = F,
          show_column_names = T)

}


plot_Spectra_Precursor_Int <- function(sp,
                                       precursor.ratio = 0.01,
                                       label.top = 10){

  sp.data <-get_Spectra_data(sp)%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = intensity,
                  int.rank = rank(-intensity) ,
                  matched = intensity > max(intensity,na.rm = T)/10 )

  label.df <- dplyr::filter(sp.data , int.rank  < label.top)
  ymax.abs <- max(abs(sp.data$intensity))
  ggplot(sp.data)+
    geom_hline(yintercept = 0 , linewidth = 0.2,col = "grey")+
    geom_segment(aes(x = x,y =y,xend = xend,
                     yend = yend,col = matched),
                 linewidth = 0.2,
                 show.legend = F)+
    geom_point(aes(x = x, y = yend , alpha = matched,col = matched),
               show.legend = F,size = 0.5)+
    geom_hline(
      yintercept =Spectra::precursorIntensity(sp)*precursor.ratio
      #yintercept = median(sp.data$intensity)
               )+
    geom_text(aes(x = quantile(range(x),0.9),
                  y = quantile(range(yend),0.9),
                  label = format(Spectra::precursorIntensity(sp),
                                 digit = 3,scientific=T)),
              check_overlap = T)+
    ggrepel::geom_text_repel(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                             #segment.size = 0.1,
                             data = label.df)+
    scale_color_manual(values = c(`FALSE` = "grey",`TRUE` = "#80B1D3"))+
    scale_alpha_manual(values = c(`FALSE` = 0,`TRUE` =1))+
    scale_y_log10(labels = scales::scientific,
                  expand = expansion(0,0))+
    labs(x = "Mz",y = "Intensity")+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))->p
  p


}

#setMethod("plotSpec",
#          signature = "Spectra",
#          definition = function(object){
#            plot_Spectra(object)
#          })

setMethod("plot",
          signature = "Spectra",
          definition = function(x){
            plot_Spectra(x)
          })

export_Spectra <- function(sp,
                           format = "msp",
                           file ){
  bk = switch(format,
              "msp" = MsBackendMsp::MsBackendMsp(),
              "mona" = MsBackendMsp::MsBackendMsp(),
              "mgf" = MsBackendMgf::MsBackendMgf())
  Spectra::export(sp,
                  file = file,
                  mapping = spectraVariableMapping(bk),
                  backend = bk)


}

export_Spectra_peak_list_for_cfm <- function(sp,
                                   file ){
  readr::write_lines(NULL,file = file)
  for (ce in collisionEnergy(sp)) {
    sp.this <- sp[collisionEnergy(sp)==ce]
    energy.type <- switch(as.character(ce),
           "10" = "energy0",
           "20" = "energy1",
           "40" = "energy2")
    readr::write_lines(energy.type,file = file,append = T)
    data.to.write <- get_Spectra_data(sp.this)%>%
      dplyr::arrange(mz)%>%
      dplyr::mutate(to.write = paste0(mz," ",intensity))%>%
      dplyr::select(to.write)%>%
      as.matrix()
    readr::write_lines(data.to.write,file = file,append = T)

  }


  return(invisible(file))


}

load_Spectra <- function(file) {

  format <- get_file_formate(file)

  bk = switch(format,
              "msp" = MsBackendMsp::MsBackendMsp(),
              "mona" = MsBackendMsp::MsBackendMsp(),
              "mgf" = MsBackendMgf::MsBackendMgf())

  Spectra(file,source = bk ,
          BPPARAM = SerialParam())
}

plot_Spectra_Injection <- function(sp){

  sp.data <- Spectra::spectraData(sp)%>%as.data.frame()
  ggplot(sp.data)+
    geom_point(aes(x = log10(Spectra::precursorIntensity),
                   y = log10(totIonCurrent),
                   col = injectionTime),
               shape = 19,
               size = 2,
               stroke = 0,
               #col = "transparent",
               alpha = 0.3)+
    ggsci::scale_fill_gsea()+
    ggsci::scale_color_gsea()+
    theme_bw()


}




#' get_CFM_data_Spectra
#'
#' @title Convert CFM Data to Spectra Object
#' @description Converts CFM (Competitive Fragmentation Modeling) peak assignment data
#'   to a `Spectra` object with spectra at three collision energy levels (10, 20, 40 eV).
#' @param cfmd A CFM data object containing peak assignment information.
#'   Typically obtained from `read_CFM_xxx` functions.
#'
#' @return A `Spectra` object containing MS/MS spectra at three collision energy levels.
#' @export
#'
get_CFM_data_Spectra <- function(cfmd ){

  cfm.df <- cfmd@peak_assignment
  sp.list<- list()
  for (i in 0:2) {
    this.sp.data <- cfm.df%>%
      dplyr::filter(energy == paste0("energy",i))%>%
      dplyr::distinct(mz,intensity)
    this.sp.df <- S4Vectors::DataFrame(collisionEnergy = switch(as.character(i),
                                                     "0" = 10,
                                                     "1" = 20,
                                                     "2" = 40))
    this.sp.df$mz <- list(this.sp.data$mz)
    this.sp.df$intensity  <-list( this.sp.data$intensity)
    this.sp <- Spectra(this.sp.df)
    sp.list[[ paste0("energy",i)]] <- this.sp
  }
  sp <- do.call(c,unname(sp.list))
  return(sp)


}

get_Spectra_adduct_expand <- function(sp,
                                      selected_adduct = MSCC::adduct.table$Adduct){


  #sp <- Spectra_database
  adduct.table <- MSCC::adduct.table%>%
    dplyr::filter(Adduct%in% selected_adduct)
  adduct.count <- dplyr::count(adduct.table,Ion_mode)
  adduct.pos <- dplyr::pull(filter(adduct.table,Ion_mode=="positive"),Adduct)
  adduct.neg <- dplyr::pull(filter(adduct.table,Ion_mode=="negative"),Adduct)
  sp.data <- spectraData(sp)%>%
    as.data.frame()

  chem_unique <- sp.data%>%
    dplyr::distinct(formula,polarity)%>%
    dplyr::mutate(temp_idx = 1:n(),
                  adduct.count = case_when(polarity == 0~adduct.count$n[1],
                                           polarity == 1~adduct.count$n[2]) )%>%
    dplyr::slice(rep(temp_idx,adduct.count))%>%
    dplyr::group_by(temp_idx)%>%
    dplyr::mutate(adduct = 1:n(),
                  adduct = case_when(polarity == 0~adduct.neg[adduct],
                                     polarity == 1~adduct.pos[adduct]))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(MSCC::chemform_adduct(formula,adduct),
                  match.id = paste0(formula,adduct))

  sp.data <- sp.data%>%
    dplyr::mutate(sp_temp_idx = 1:n(),
                  adduct.count = case_when(polarity == 0~adduct.count$n[1],
                                           polarity == 1~adduct.count$n[2]) )%>%
    dplyr::slice(rep(sp_temp_idx,adduct.count))%>%
    dplyr::group_by(sp_temp_idx)%>%
    dplyr::mutate(adduct = 1:n(),
                  adduct = case_when(polarity == 0~adduct.neg[adduct],
                                     polarity == 1~adduct.pos[adduct]))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(match.id = paste0(formula,adduct),
                  chem_unique[match(match.id , chem_unique$match.id),  ]  )

  sp.expand <- sp[sp.data$sp_temp_idx]
  sp.expand$adduct <- sp.data$adduct
  sp.expand$precursorMz <- sp.data$chemform.adduct.mz

  return(sp.expand)


}

Spectra_set_MEM_backend <- function(sp){

  if(class(sp@backend)=="MsBackendCompDb"){
    sp$centroided <-T
    sp$smoothed <-T
  }

  sp <- Spectra::setBackend(sp,
                  backend = Spectra::MsBackendMemory())
  sp
}




#' plotly_Spectra
#'
#' @title Create Interactive Plotly Mass Spectrum
#' @description Creates an interactive plotly visualization of a single mass spectrum.
#'   Hovering over peaks shows m/z and intensity values. The top N peaks are highlighted.
#' @param sp A `Spectra` object. If multiple spectra are provided, only the first is used.
#' @param label.top Integer. Number of highest intensity peaks to highlight.
#'   Default is `10`.
#' @param show.info Logical. If `TRUE`, displays precursor information (m/z, intensity,
#'   collision energy) in the plot title. Default is `FALSE`.
#'
#' @return A `plotly` object displaying an interactive mass spectrum.
#' @export
#'
plotly_Spectra <- function(sp,label.top = 10,show.info = F){

  if (length(sp)>1) {
      warning("more than 1 spectra input, select the first")
    sp <- sp[1]
  }
  sp.data <-get_Spectra_data(sp)%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = intensity,
                  int.rank = rank(-intensity) ,
                  highlight = intensity > quantile(intensity,(n()-label.top)/n()),
                  size = ifelse(highlight,1,0))

  if(show.info){
    label.to.show <- paste0(
      "precursorMz = ", round(precursorMz(sp),digits = 4),"\n",
      "precursorIntensity = ", format(Spectra::precursorIntensity(sp),digits = 3,sci=T),"\n",
      "collisionEnergy = ", (collisionEnergy(sp)),"\n"
    )
  }else{
    label.to.show <- ""
  }


  plot_ly(sp.data)%>%
    add_segments(x = ~x, xend = ~ xend,
                 y = ~y,yend = ~ yend, color = I("grey"),
                 showlegend = F)%>%
    add_markers(x = ~xend ,y = ~yend,
                size = I(5),
                showlegend = F,
                hovertemplate = "mz:%{x}\nintensity:%{y}<extra></extra>")%>%
    layout(xaxis =list(title =  "mz"),
           title = list(text = label.to.show,
                        x = 0.15,y = 0.9,
                        xanchor = "left"),
           yaxis =list(title =  "intensity",linewidth = 1))




}


#' plotly_Spectra_mirror
#'
#' @title Create Interactive Mirror Plot of Two Spectra
#' @description Creates an interactive plotly mirror visualization comparing two mass spectra.
#'   The first spectrum is displayed on top, the second on bottom (negative intensity).
#'   Matched peaks are highlighted with markers.
#' @param sp1 A `Spectra` object for the top spectrum.
#' @param sp2 A `Spectra` object for the bottom spectrum.
#'
#' @return A `plotly` object displaying a mirror plot of the two spectra.
#' @export
#'
plotly_Spectra_mirror <- function(sp1,sp2){


  sp.data <-rbind(get_Spectra_data(sp1)%>%
                    dplyr::mutate(sp.id=1),
                  get_Spectra_data(sp2)%>%
                    dplyr::mutate(sp.id=2))%>%
    dplyr::mutate(groupMz(mz,return.type = "data.frame"))%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = case_when(sp.id==1~intensity,
                                   T~-intensity))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate(matched = n()==2)


  plot_ly(sp.data)%>%
    add_segments(x = ~x, xend = ~ xend,
                 y = ~y,yend = ~ yend, color = I("grey"),
                 showlegend = F)%>%
    dplyr::filter(matched)%>%
    add_markers(x = ~xend ,y = ~yend,
                size = I(5),
                showlegend = F,
                hovertemplate = "mz:%{x}\nintensity:%{y}<extra></extra>")%>%
    layout(xaxis =list(title =  "mz"),
           yaxis =list(title =  "intensity",linewidth = 1))


}


#' plotly_Spectra_iso_mirror
#'
#' @title Create Interactive Isotopic Mirror Plot
#' @description Creates an interactive plotly mirror visualization comparing a spectrum
#'   with its isotopically labeled counterpart. Peaks are matched based on the expected
#'   isotopic mass difference and highlighted accordingly.
#' @param sp A `Spectra` object for the natural abundance spectrum (top).
#' @param sp.iso A `Spectra` object for the isotopically labeled spectrum (bottom).
#' @param ppm Numeric. Parts per million tolerance for matching peaks. Default is `10`.
#' @param iso_mass_diff Numeric. Expected mass difference between isotopes.
#'   Default is `1.003355` (carbon-13 difference).
#' @param iso_count Integer. Number of isotope labels to account for. Default is `1`.
#'
#' @return A `plotly` object displaying an isotopic mirror plot with matched peaks highlighted.
#' @export
#'
plotly_Spectra_iso_mirror <- function(sp,sp.iso ,
                                      ppm = 10,
                                      iso_mass_diff = 1.003355,
                                      iso_count = 1){


  sp.data <- get_Spectra_data(sp)%>%
    dplyr::mutate(sp.id=1,hl = F)
  sp.iso.data <- get_Spectra_data(sp.iso)%>%
    dplyr::mutate(sp.id=2,hl = F)
  matched.id <- match_mz(sp.data$mz+iso_mass_diff*iso_count,
                         sp.iso.data$mz,mz.ppm = ppm)
  sp.data$hl[!is.na(matched.id)] <-T
  sp.iso.data$hl[(matched.id)] <-T


  sp.data <-rbind(sp.data,sp.iso.data)%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = case_when(sp.id==1~intensity,
                                   T~-intensity))


  plot_ly(sp.data)%>%
    add_segments(x = ~x, xend = ~ xend,
                 y = ~y,yend = ~ yend, color = I("grey"),
                 showlegend = F)%>%
    dplyr::filter(hl)%>%
    add_markers(x = ~xend ,y = ~yend,
                size = I(5),
                showlegend = F,
                hovertemplate = "mz:%{x}\nintensity:%{y}<extra></extra>")%>%
    add_text(x = 1.1*max(precursorMz(sp)),y = 0.8*max(abs(sp.data$yend))*c(1,-1),
             text = paste0("M",c(0,iso_count)),
             showlegend = F,
             size= I(20))%>%
    add_text(x = 1.1*max(precursorMz(sp)),
             y = 0.5*max(abs(sp.data$yend))*c(1,-1),
             text = paste0("Collision Energy\n",c(unique(collisionEnergy(sp)),
                                                  unique(collisionEnergy(sp.iso)))),
             showlegend = F,
             size= I(10))%>%
    layout(xaxis =list(title =  "mz",range = c(0,1.3*max(sp.data$x))),
           yaxis =list(title =  "intensity",linewidth = 1))


}


#' combineSpectra_ce_max_precursor
#'
#' @title Select Spectra with Maximum Precursor Intensity per Collision Energy
#' @description For each unique collision energy in the Spectra object, selects the spectrum
#'   with the highest precursor intensity. Useful for selecting the most representative
#'   spectrum when multiple spectra exist at each collision energy.
#' @param sp A `Spectra` object containing multiple spectra.
#'
#' @return A `Spectra` object with one spectrum per collision energy (the one with
#'   maximum precursor intensity).
#' @export
combineSpectra_ce_max_precursor  <- function(sp){

  sp.data <- spectraData(sp)%>%
    as.data.frame()%>%
    rownames_to_column("sp.name")%>%
    dplyr::group_by(collisionEnergy)%>%
    dplyr::slice_max(precursorIntensity)

  sp[sp.data$sp.name]
}





#' Spectra_get_noise
#'
#' @title Estimate Noise Level for Spectra
#' @description Estimates the noise level for each spectrum using a density-based method.
#'   The noise is estimated as the intensity value at the maximum of the log10-transformed
#'   intensity density distribution. Also calculates signal-to-noise ratio (SNR).
#' @param sp A `Spectra` object.
#'
#' @return The input `Spectra` object with two additional spectrum variables:
#'   `noise` (estimated noise level) and `snr` (signal-to-noise ratio, calculated as
#'   base peak intensity divided by noise).
#' @export
#'
Spectra_get_noise <- function(sp){

  .f <- function(x){

    if (length(x[,2])<2) return(0)
    den <- density(log10(x[,2]))

    ### density max
    {
      id <- which.max(den$y)

    }

    ### density mutation
    {
      #den.y.diff <- diff(den$y)
      #den.y.diff2 <- diff(den.y.diff)
      #den.y.diff[den.y.diff<0] <- -1
      #den.y.diff[den.y.diff>0] <- 1
      #id <- which(diff(den.y.diff)==2)[1]
      #par(mfrow = c(2,1))
      #plot(den$y)
      #abline(v = id)
      #plot(den.y.diff)
      #abline(v = id)
    }
    return(10^den$x[id])
  }

  sp.ms2.noise <- bplapply(peaksData(sp),.f ,
                           BPPARAM = SerialParam(progressbar = T))
  sp$noise <- unlist(sp.ms2.noise)
  sp$snr <- sp$basePeakIntensity/sp$noise
  return(sp)
}


#' Spectra_filter_noise
#'
#' @title Filter Peaks Below Noise Level
#' @description Filters out mass spectral peaks that are below the noise level.
#'   Requires that a `noise` spectrum variable has been previously calculated
#'   (e.g., using [Spectra_get_noise]). Peaks with intensity less than or equal
#'   to the noise level are removed.
#' @param sp A `Spectra` object with a `noise` spectrum variable.
#'
#' @return A filtered `Spectra` object with peaks below the noise level removed.
#' @export
Spectra_filter_noise <- function(sp){

  if (!"noise" %in% spectraVariables(sp)) {
    warning("No noise varibales, skip")
    return(sp)
  }
  nf <-  function(z,  noise ,...) {

    idx <- z[, 2] > noise
    z[idx,,drop = F]
  }
  sp <- Spectra::addProcessing(sp,FUN = nf,
                               spectraVariables = "noise")%>%
    Spectra::applyProcessing(BPPARAM = SerialParam())

  return(sp)
}


#setMethod(noise,"Spectra",
#          definition = function(object )object$noise)


setMethod(sampleNames,"Spectra",
          definition = function(object){
            basename(object$dataOrigin)
          })

#' Spectra_get_purity
#'
#' @title Estimate Precursor Purity for Mass Spectra
#' @description Estimates the purity of the precursor ion for MS/MS spectra.
#'   For MS2 spectra, calculates the proportion of precursor intensity within the
#'   isolation window that belongs to the targeted precursor m/z.
#'   For MS1 spectra, estimates purity using an associated MS1 scan.
#' @param sp A `Spectra` object for which to estimate purity.
#' @param msLevel Integer. MS level to process. `1` for MS1 spectra (requires `sp.ms1`),
#'   `2` for MS2 spectra. Default is `2`.
#' @param sp.ms1 Optional. A `Spectra` object containing MS1 spectra used for
#'   purity calculation when `msLevel = 1`. If `NULL`, MS1 spectra are imported
#'   from the raw data files.
#'
#' @return The input `Spectra` object with additional spectrum variables:
#'   For MS2: `ms2_purity` (proportion of precursor in isolation window).
#'   For MS1: `ms1.intensity`, `ms2.purity`, and `ms2.ppm`.
#' @export
#'
Spectra_get_purity <- function(sp,msLevel = 2,sp.ms1= NULL){

  if (msLevel == 1) {
    sp <- .Spectra_get_purity_ms1(sp,sp.ms1)
  }
  if (msLevel == 2) {
    sp <- .Spectra_get_purity_ms2(sp)
  }

  return(sp)

}


.Spectra_get_purity_ms1 <- function(sp,sp.ms1= NULL){


  sp.raw <- sp
  sp$temp.id = 1:length(sp)

  ###extract sp
  {
    if (is.null(sp.ms1)) {
      message("Import MS1 Spectra")
      ms.files <- unique(dataOrigin(sp))
      sp.ms1 <- Spectra(ms.files,msLevel=1)%>%
        filterMsLevel(1)
    }
    sp.ms2.list <- split(sp,sp$dataOrigin)[unique(sp$dataOrigin)]
    sp.ms1.list <- split(sp.ms1,sp.ms1$dataOrigin)[unique(sp$dataOrigin)]
  }

  ### func
  .f <- function(ms1,ms2){

    ms2$ms1.intensity<-NA
    ms2$ms2.purity<-NA
    ms2$ms2.ppm <- NA
    for (i in 1:length(ms2)) {
      x.sp <- ms2[i]
      precursor.ms1 <- max(which(rtime(ms1)<rtime(x.sp)))
      x <- peaksData(ms1[precursor.ms1])[[1]]
      x <- x[x[,1]>x.sp$isolationWindowLowerMz&
               x[,1]< x.sp$isolationWindowUpperMz,,drop =F]
      #id <- which(round(x[,1],8)==round(x.sp$isolationWindowTargetMz ,8))
      id <- which.min(abs(x[,1]-x.sp$isolationWindowTargetMz ))
      ms2$ms1.intensity[i] <- x[id,2]
      ms2$ms2.purity[i] <- x[id,2]/sum(x[,2])
      ms2$ms2.ppm[i] <- min(abs(x[,1]-x.sp$isolationWindowTargetMz ))/
        x.sp$isolationWindowTargetMz*1e6

    }
    return(ms2)

  }


  ### process
  {
    sp.x <- bpmapply(.f,sp.ms1.list,
                     sp.ms2.list,
                     BPPARAM = SerialParam(progressbar = T))
    sp.x <- do.call(c,unname(sp.x))


    sp.raw$ms1.intensity <- sp.x$ms1.intensity[order(sp.x$temp.id)]
    sp.raw$ms2.purity <- sp.x$ms2.purity[order(sp.x$temp.id)]
    sp.raw$ms2.ppm <- sp.x$ms2.ppm[order(sp.x$temp.id)]
  }

  return(sp.raw)
}

.Spectra_get_purity_ms2 <- function(sp,ppm = 10){

  pmz <- precursorMz(sp)
  iwlm <- isolationWindowLowerMz(sp)
  iwum <- isolationWindowUpperMz(sp)

  spd <- peaksData(sp)

  sp$ms2_purity <- sapply(seq_along(spd),
                          function(x){
                            pd <- spd[[x]]
                            iniw <-  between(pd[,1],iwlm[x],iwum[x])
                            isp <- between.range(pd[,1],
                                                 mz.range.ppm(pmz[x],ppm))
                            purity <- sum(pd[isp,2])/sum(pd[iniw,2])
                            ifelse(is.nan(purity),NA,purity)
                          })

  return(sp)
}



get_spectra_ion_purity <- function(sp,ion_mz,ppm = 10,isolation_half_window = 0.2){

  sp.peaks.data <- peaksData(sp)
  names(sp.peaks.data) <- basename(dataOrigin(sp))
  ion.purity <- sapply(sp.peaks.data,function(x){
    mz.diff <-abs(x[,1]-ion_mz)
    idx.in.window <- which( mz.diff < isolation_half_window)
    idx.target <- which(mz.diff/ion_mz < ppm*1e-6)
    unname(sum(x[idx.target,2])/sum(x[idx.in.window,2]))
  })
  ion.purity[is.infinite(ion.purity)] <- 0
  ion.purity[is.nan(ion.purity)] <- 0
  return(ion.purity)

}


get_spectra_by_name <- function(sp,sp.name){

  idx <- match(sp.name,spectraNames(sp))

  sp[idx]

}
