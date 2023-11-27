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



makeSpectra <- function(precursorMz ,
                        rtime ,...){

  Spectra::Spectra(S4Vectors::DataFrame(precursorMz = precursorMz,
                             rtime = rtime,
                             ...))

}


filterSpectraIntensity <- function(sp,ratio){

  Spectra::filterIntensity( sp, intensity = function(z) z/max(z) > ratio )

}


filterSpectra_below_PrecursorMz <- function(sp){

  if (!length(sp)) return(sp)
  nf <-  function(z,  precursorMz ,...) {

    idx <- z[, "mz"] < precursorMz-0.5
    z[idx,,drop = F]
  }
  sp <- Spectra::addProcessing(sp,FUN = nf, spectraVariables = "precursorMz")%>%
    applyProcessing(BPPARAM = SerialParam())

  return(sp)
}



normalizeSpectra <- function(sp){

  nf <-  function(z, ...) {
    #z[,"maxIntensity"] <- max( z[, "intensity"] )
    z[, "intensity"] <- z[, "intensity"] /
      max(z[, "intensity"], na.rm = TRUE) * 100
    z
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
#' @examples
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

combineSpectra_groupby_ce <- function(sp,minProp = 0.5){

  sp.ce <- Spectra::combineSpectra(sp,peaks = "intersect",
                                   intensityFun = median,
                 f = sp$collisionEnergy,
                 minProp=minProp,
                 ppm = 5)


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
