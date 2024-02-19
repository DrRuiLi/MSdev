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
                             ...),backend = MsBackendMemory())

}

makeEmptySpectra <- function(...){
  Spectra::Spectra(S4Vectors::DataFrame(precursorMz = 0,
                                        rtime = 0,
                                        ...),backend = MsBackendMemory())

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
    applyProcessing(BPPARAM = SerialParam())

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
    applyProcessing(BPPARAM = SerialParam())

  return(sp)

}

#' get_Spectra_data
#'
#' @param sp `Spectra` object
#' @param var any variables of `Spectra::spectraVariables()`
#'
#' @return a dataframe with mz and intensity and vars
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

  mz <- mz(sp)%>%
    as.data.frame()%>%
    dplyr::select(everything(),mz = value)
  int <- intensity(sp)%>%
    as.data.frame()%>%
    dplyr::select(everything(),intensity = value)

  spec.df <- mz %>%
    dplyr::mutate(int,
                  sp.id = group,
                  sp.data[match(sp.id,sp.data$id)   ,var])%>%
    dplyr::select(sp.id ,!all_of(c("group","group_name")))
  return(spec.df)
}

#' combineSpectra_groupby_ce
#'
#' @param sp Spectra
#' @param minProp NUM
#'
#' @return Spectra
#' @export
#'

combineSpectra_groupby_ce <- function(sp,
                                      minProp = 0.5,
                                      ppm = 5,
                                      ...){

  sp.ce <- Spectra::combineSpectra(sp,
                                   peaks = "intersect",
                                   intensityFun = median,
                 f = sp$collisionEnergy,
                 weighted =T,
                 minProp=minProp,
                 ppm = ppm,
                 ...)


}

get_Spectra_ms2_feature_id <- function(sp,
                                    featuredef ){



  sp.data <- Spectra::spectraData(sp)%>%
    as.data.frame()%>%
    dplyr::mutate(precursorMZ = precursorMz,
                  retentionTime = rtime)%>%
    dplyr::filter(msLevel == 2)

  match.df <- match_mz_rt(featuredef$mzmed,featuredef$rtmed,
                          sp.data$precursorMZ,
                          sp.data$retentionTime)
  match.df <- match.df%>%
    dplyr::mutate(featuredef[ion1,],
                  sp_id = sp.data$sp_id[ion2],
                  sp_rt = sp.data$retentionTime[ion2])%>%
    dplyr::filter(sp_rt < peakRtMax&sp_rt>peakRtMin )%>%
    dplyr::group_by(sp_id)%>%
    dplyr::slice_min(mz.error)



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
    dplyr::mutate(groupMz(mz))%>%
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
                             segment.size = 0.1,
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
    geom_hline(yintercept = 0 , size = 0.2,col = "grey")+
    geom_segment(aes(x = x,y =y,xend = xend,
                     yend = yend,col = matched),
                 size = 0.2,
                 show.legend = F)+
    geom_point(aes(x = x, y = yend , alpha = matched,col = matched),
               show.legend = F,size = 0.5)+
    ggrepel::geom_text_repel(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                             segment.size = 0.1,
                             data = label.df)+
    scale_color_manual(values = c(`FALSE` = "grey",`TRUE` = "#80B1D3"))+
    scale_alpha_manual(values = c(`FALSE` = 0,`TRUE` =1))+
    scale_y_continuous(expand = expansion(0,0),
                       #limits = c(0,ymax.abs*1.1),
                       labels = scales::scientific)+
    labs(x = "Mz",y = "Intensity")+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))->p
  p


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
    dplyr::mutate(groupMz(mz))%>%
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
    ggrepel::geom_text_repel(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                             segment.size = 0.1,
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
    dplyr::mutate(groupMz(mz))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate( highlight = length(unique(collisionEnergy)) >
                     length(levels(collisionEnergy))*0.5,
                   hx = min(x),
                   hxend = max(x),
                   hy = min(y),
                   hyend =max(y))%>%
    dplyr::ungroup()%>%
    dplyr::filter(highlight)

  ggplot(sp.data,aes(x = cex , y = intensity ,col = mz.group))+
    geom_point()+
    geom_smooth(formula = y~x,method = "loess")+
    scale_color_manual(values = randomcoloR::distinctColorPalette(30),
                       labels = sprintf("%.4f",unique(sp.data$mz)))+
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
    dplyr::mutate(groupMz(mz))%>%
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
                             segment.size = 0.1,
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
  sp.simlilarity <- compareSpectra(sp)
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
      yintercept =precursorIntensity(sp)*precursor.ratio
      #yintercept = median(sp.data$intensity)
               )+
    geom_text(aes(x = quantile(range(x),0.9),
                  y = quantile(range(yend),0.9),
                  label = format(precursorIntensity(sp),
                                 digit = 3,scientific=T)),
              check_overlap = T)+
    ggrepel::geom_text_repel(aes(x = x, y = yend ,
                                 label = format(mz,digit = 4,nsmall = 4)),
                             size =2,
                             col = "#00000088",
                             segment.size = 0.1,
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

setMethod("plotSpec",
          signature = "Spectra",
          definition = function(object){
            plot_Spectra(object)
          })

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
  write_lines(NULL,file = file)
  for (ce in collisionEnergy(sp)) {
    sp.this <- sp[collisionEnergy(sp)==ce]
    energy.type <- switch(as.character(ce),
           "10" = "energy0",
           "20" = "energy1",
           "40" = "energy2")
    write_lines(energy.type,file = file,append = T)
    data.to.write <- get_Spectra_data(sp.this)%>%
      dplyr::arrange(mz)%>%
      dplyr::mutate(to.write = paste0(mz," ",intensity))%>%
      dplyr::select(to.write)%>%
      as.matrix()
    write_lines(data.to.write,file = file,append = T)

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

  sp.data <- spectraData(sp)%>%as.data.frame()
  ggplot(sp.data)+
    geom_point(aes(x = log10(precursorIntensity),
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



get_Spectra_CFM <- function(sp){



}


#' get_Spectra_from_CFM
#'
#' @param cfm.data from read_CFM_xxx
#'
#' @return Spcetra
#' @export
#'
get_Spectra_from_CFM <- function(cfm.data ){

  cfm.df <- cfm.data$peak_assignment
  sp.list<- list()
  for (i in 0:2) {
    this.sp.data <- cfm.df%>%
      dplyr::filter(energy == paste0("energy",i))%>%
      dplyr::distinct(mz,intensity)
    this.sp.df <- DataFrame(collisionEnergy = switch(as.character(i),
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

get_Spectra_MEM_backend <- function(sp){

  if(class(sp@backend)=="MsBackendCompDb"){
    sp$centroided <-T
    sp$smoothed <-T
  }

  sp <- setBackend(sp,
                  backend = MsBackendMemory())
  sp
}
