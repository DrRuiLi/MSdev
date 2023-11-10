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

  Spectra::filterIntensity( sp, intensity = function(z) z/max(z) > r )

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
  sp$precursorIntensity[sp$precursorIntensity==0] <- sp$totIonCurrent[sp$precursorIntensity==0]
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

  combineSpectra(sp,peaks = "intersect",
                 f = sp$collisionEnergy,
                 minProp=minProp,
                 ppm = 25)


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


plot_Spectra_Mirror<- function(sp1,sp2){

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
    scale_y_continuous(expand = expansion(0,0))+
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



