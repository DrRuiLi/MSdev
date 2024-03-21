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
                                        score = this.seed.df$score,
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
#' @import magrittr tidyverse
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
        seed.sp.c <- Spectra_filter_noise(x$M0)
        seed.sp.c <- normalizeSpectra(sp = seed.sp.c,norm_to = "tic")
        seed.sp.c<- combineSpectra_groupby_ce(seed.sp.c,ppm = 10,
                                    minProp = 0.3,
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

    CFM_result <- cfm_data_get_fragment_group(CFM_result)
    CFM_result <- CFM_data_get_igraph(CFM_result)
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

  .f <- function(x.iso.cfm){

    cfmd <- x.iso.cfm$CFM_annotation
    iso.count <- stringr::str_extract(names(x.iso.cfm),"[:digit:]+")%>%
      as.numeric()%>%na.omit()
    atom_ele <-  vdata(cfmd@fragment_igraph$Fragment001)$atom
    names(atom_ele) <- vdata(cfmd@fragment_igraph$Fragment001)$name
    c_ele <- atom_ele[atom_ele=="C"]
    for (i in iso.count) {

      this.sp <-x.iso.cfm[[paste0("M",i)]]
      this.sp <- Spectra_filter_noise(this.sp)
      sp.frag <- CFM_annotate_isotopologues(this.sp,
                                 cfmd  = cfmd,
                                 iso.count = i)%>%
        dplyr::filter(!is.na(iso))
      idx <- split(1:nrow(sp.frag),sp.frag$fragment_group)
      frag.group <- lapply(idx,
             function(x){
              x.df <- sp.frag[x,]
              x.int <- x.df%>%
                dplyr::select(-mz,-collisionEnergy)%>%
                tidyr::pivot_wider(names_from ="iso",
                                   id_cols = "sp.id",
                                   values_from = "intensity")%>%
                tibble::column_to_rownames("sp.id")%>%
                dplyr::select(dplyr::starts_with("M"))%>%
                as.matrix()
              x.int[is.na(x.int)] <- 0
              if (!"M0" %in%colnames(x.int))
                x.int <- cbind(matrix(0,nrow(x.int),1,dimnames = list(NULL,"M0")),x.int)
              x.weight <- rowSums(x.int)
              x.int <- t(apply(x.int,1,function(z) z/sum(z)))
              x.int <- x.int[,order(colnames(x.int)),drop =F]
              x.int <- apply(x.int,2,weighted.mean,m = x.weight)
              return(x.int)
             })

      frag.c <- matrix(ncol = length(c_ele),
                       nrow = length(frag.group),
                       dimnames = list(names(frag.group),
                                       names(c_ele)))
      for (j in seq_along(frag.group)) {

        this.frag.group <- names(frag.group)[j]
        this.frags <- cfmd@fragment_define[cfmd@fragment_define$fragment_group==this.frag.group,]
        this.frag.ratio <-frag.group[[j]]
        this.iso.expectation <- sum(str_extract_num(names(this.frag.ratio))*this.frag.ratio)
        this.frag.atom <- get_cfm_data_fg_atom_map(cfmd,this.frag.group)
        this.frag.c <- this.frag.atom[names(c_ele)]
        this.frag.c <- this.frag.c*this.iso.expectation/sum(this.frag.c)
        this.frag.c <- this.frag.c[this.frag.c!=0]
        frag.c[this.frag.group,names(this.frag.c)] <- this.frag.c
        #vis_sdf_igraph(cfmd@fragment_igraph[[1]],show.label = F,highlight = names(this.frag.atom))
      }

      heatmap_atom_iso_prob(t(frag.c))

    }






  }



}


