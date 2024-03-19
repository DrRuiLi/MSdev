#' get_MSdev_isotopologues
#'
#' extract iso-labeled compound info and spectra,
#' filter compound without 3 CE spectra of M0
#'
#' @param object MSdev
#'
#' @return a list of isotopologues
#' @export
#'
get_MSdev_isotopologues <- function(object){

  iso.list <- list()
  sp.ms2 <- object@spectra$MS2_Spectra
  cpdb <- CompoundDb::CompDb(object@projectInfo$CompoundDB_path)
  ### Requirement in xcms features
  ### C13_seed
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.fdf <- get_xcms_feature_definitions(xcms.xcms)
    db.info <- get_CompDb_info(xcms.fdf$compound_id,
                             cpdb = cpdb,
                             keys = c("name","kegg_id",
                                      "formula", "inchikey","smiles"))
    xcms.fdf <- cbind(xcms.fdf,db.info)
    iso <- grep(pattern = "_seed",colnames(xcms.fdf),value = T)
    xcms.fdf[,"iso_seed"] <- xcms.fdf[,iso]
    xcms.fdf[,"iso_count"] <- xcms.fdf[,sub(x = iso,pattern = "_seed",replacement = "_count")]
    xcms.fdf <- xcms.fdf%>%
      dplyr::group_by(iso_seed)%>%
      dplyr::filter(!is.na(iso_seed),
                    n() > 1,
                    rep(any(is_labeled) ,n()))%>%
      dplyr::ungroup()

    seed.id <- unique(xcms.fdf$iso_seed)
    for (j in seq_along(seed.id)) {

      this.fdf <- xcms.fdf%>%
        dplyr::filter(iso_seed == seed.id[j])
      this.seed.df <- this.fdf%>%
        dplyr::filter(feature_id==iso_seed)
      this.seed.fid <- this.seed.df$feature_id
      this.sp <-  sp.ms2[this.fdf$ms2_id%>%unlist()]
      this.list <- list()
      ### compound_info
      {
        this.list$compound_info <- list(name = this.seed.df$name,
                                        compound_id = this.seed.df$compound_id,
                                        formula = this.seed.df$formula,
                                        smiles = this.seed.df$smiles,
                                        mz_ref = this.seed.df$mz_ref,
                                        adduct = this.seed.df$adduct,
                                        polarity = this.seed.df$polarity)
      }
      ### Spectra split
      if (length(this.sp)>0) {
        this.sp.list <- split(this.sp,this.sp$feature_id)
        names(this.sp.list) <- paste0("M",
                                      this.fdf$iso_count[ match(names(this.sp.list),this.fdf$feature_id)])

        this.list <- append(this.list,this.sp.list)
      }

      ### filter
      {
        if (!"M0" %in% names(this.list)) next
        if (!all(c(10,20,40) %in% collisionEnergy(this.list$M0))) next
        if (is.na(this.list$compound_info$smiles)) next
        if (length(this.list)<=2) next

      }

      iso.list[[paste0(seed.id[j],"_",pol)]] <- this.list
    }
  }

  return(iso.list)

}


#' get_isotopologues_CFM_annotation
#'
#' @param iso.list list of cfm data
#'
#' @return  list of cfm data
#' @export
#'
get_isotopologues_CFM_annotation<- function(iso.list,
                                            BPPARAM = SnowParam(
                                              workers = 12,
                                              progressbar = T)){

  ff <- function(x){

    ### combine sp of seed to cfm-annotate
    {
      ### process M0 Spectra
      {
        seed.sp.c <- x$M0%>%
          Spectra_filter_noise()%>%
          normalizeSpectra("tic")%>%
          combineSpectra_groupby_ce(ppm = 10,
                                    minProp = 0.2,
                                    plot = F)
      }
      CFM_result <- NA
      try.return <- try(
        CFM_result <- CFM_annotate(smiles_or_inchi = x$compound_info$smiles,
                                     spectrum_file = seed.sp.c,
                                     ppm_mass_tol = 5,
                                     abs_mass_tol = 0.002,
                                     param_adduct = switch(x$compound_info$polarity,
                                                           "0"="[M-H]-",
                                                           "1"="[M+H]+") )
      )
      if (class(CFM_result) !="CFM_data") return(x)
    }
    x$CFM_annotation <-CFM_result
    return(x)
  }




  iso.cfm <- bplapply(iso.list,ff,
                      BPPARAM = BPPARAM
  )
  annotated <- sapply(iso.cfm,function(x) "CFM_annotation" %in% names(x))
  iso.cfm <- iso.cfm[annotated]
  return(iso.cfm)


}



#' MSdev_find_isotope_label
#'
#' @param object MSdev
#' @param isotope "[13]C"
#' @param max_label 10
#' @param ppm 20
#' @param net.degree.ratio 0.5
#'
#' @return MSdev
#' @export
#'
MSdev_find_isotope_label <- function(object,
                                     isotope = "[13]C",
                                     ...){

  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.xcms <- xcms_get_feature_isotopologues(xcms.xcms,
                                                isotope = isotope,
                                                ...)
    xcms.xcms <- xcms_get_feature_isotope_label(xcms.xcms,
                                                isotope = isotope,
                                                ...)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
  }

  return(object)

}


get_MSdev_iso_acq_list <- function(object){



  acq.list <- list()
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.fdf <- get_xcms_feature_definitions(xcms.xcms)


    ### Comp info
    cpdb <- CompDb(msdev.fs@projectInfo$CompoundDB_path)
    dbinfo <- get_CompDb_info(cpdb,
                              xcms.fdf$compound_id,
                              keys = c("name","formula","smiles"))
    xcms.fdf <- cbind(xcms.fdf,dbinfo[,c("name","formula","smiles")])


    ### calc intensity of iso and un-iso labeled sample
    xcms.pdata <- pData(xcms.xcms)%>%
      dplyr::filter(!sample.type%in% c("Blank"))
    xcms.fv <- featureValues(xcms.xcms,value = "maxo",missing = 1)[,xcms.pdata$sampleNames]
    uniso.mean <- xcms.fv[,xcms.pdata$sampleNames[is.na(xcms.pdata$isotope_label)]]%>%
      apply(1,mean)
    iso.mean <- xcms.fv[,xcms.pdata$sampleNames[!is.na(xcms.pdata$isotope_label)]]%>%
      apply(1,mean)

    xcms.fdf$mean.iso <- log10(iso.mean)
    xcms.fdf$mean.uniso <- log10(uniso.mean)


    ### iso stat
    xcms.fdf.stat <- xcms.fdf%>%
      dplyr::mutate(compound_id = case_when(
        feature_id == C13_seed~ compound_id,
        T~ NA
      ),name = case_when(
        feature_id == C13_seed~ name,
        T ~ NA
      ))%>%
      dplyr::filter(!is.na(C13_seed))%>%
      dplyr::filter(peakMaxo > 1e5)%>%
      dplyr::group_by(C13_seed)%>%
      dplyr::mutate(total.isotopologues = n(),
                    iso.maxo = max(log10(peakMaxo)))%>%
      dplyr::arrange(-iso.maxo,
                     C13_seed,
                     C13_count)%>%
      dplyr::filter(any( is_labeled ),
                    any( !is.na( compound_id ) ) )%>%
      dplyr::ungroup( )


    #edit_df_in_excel(xcms.fdf.stat)
    acq.list[[pol]] <- xcms.fdf.stat


  }
  acq.list



}


#' get_isotopologues_Spectra_process
#'
#' @param iso.list iso.list
#'
#' @return iso.list
#' @export
#' @import  Spectra
get_isotopologues_Spectra_process <- function(iso.list){

  for (i in seq_along(iso.list)) {

    this.iso <- iso.list[[i]]
    this.count <- stringr::str_extract(names(this.iso),"[:digit:]+")%>%
      as.numeric()%>%na.omit()
    for (j in this.count) {
      this.sp <- this.iso[[paste0("M",j)]]

      ### filter sp


      ### combine sp
      this.sp <- normalizeSpectra(this.sp,"max")
      this.sp <- combineSpectra_groupby_ce(this.sp,
                                           minProp = 0.49,
                                           ppm = 10)
      this.sp -> this.iso[[paste0("M",j)]]
    }
    this.iso -> iso.list[[i]]
  }

  return(iso.list)



}



get_isotopologues_label_fraction <- function(iso.list){

  .f <- function(x){

    cfm.anno <- x$CFM_annotation$peak_assignment%>%
      dplyr::mutate(groupMz(mz))%>%
      dplyr::arrange(mz)
    iso.count <- stringr::str_extract(names(x),"[:digit:]+")%>%
      as.numeric()%>%na.omit()
    for (i in iso.count) {

      this.sp <-x[[paste0("M",i)]]
      da <- CFM_annotate_Spectra(this.sp,
                                 CFM_annotation =  x$CFM_annotation,
                                 iso.count = i)

    }





  }



}


