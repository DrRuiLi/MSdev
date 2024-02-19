#' get_MSdev_isotopologues
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
      this.ms2 <- this.fdf$ms2_id%>%unlist()
      this.ms2 <-  sp.ms2[this.ms2]
      iso.list[[paste0(seed.id[j],"_",pol)]] <- list(iso.df = this.fdf,
                                                     sp = this.ms2)
    }
  }

  return(iso.list)

}


get_isotopologues_CFM_annotation<- function(iso.list){

  ff <- function(x){

    this.iso.df <- x$iso.df
    this.seed.fid <- this.iso.df%>%
      dplyr::filter(feature_id==iso_seed)%>%
      dplyr::pull(feature_id)
    this.smiles <- this.iso.df%>%
      dplyr::filter(feature_id==iso_seed)%>%
      dplyr::pull(smile)
    this.pol <- unique(this.iso.df$polarity)
    this.sp <- x$sp%>%
      split(.,x$sp$ms2_matched_feature)

    if (is.na(this.smiles)) return(x)



    ### combine sp of seed and cfm-annotate
    {
      seed.sp <- this.sp[[this.seed.fid]]
      if (!all(c(10,20,40) %in% collisionEnergy(seed.sp)))  return(x)
      seed.sp.c <- combineSpectra_groupby_ce(seed.sp)
      try.return <- try(
        x$CFM_result <- CFM_annotate(smiles_or_inchi = this.smiles,
                                     spectrum_file = seed.sp.c,
                                     ppm_mass_tol = 5,
                                     abs_mass_tol =0.002,
                                     param_adduct = switch(this.pol,
                                                           "0"="[M-H]-",
                                                           "1"="[M+H]+") )
      )
    }


    ### cfm result
    {
      for (i in nrow(this.iso.df)) {
        cfm.annotated <- x$CFM_result$peak_assignment%>%
          dplyr::filter(energy %in% paste0("energy",i),
                        !is.na(fragment_id))
        cfm.frag.def <- x$CFM_result$fragment_define
        if (!nrow(cfm.annotated)) next

      }




    }
    return(x)


  }

  iso.cfm <- bplapply(iso.list,ff,
                        BPPARAM = SnowParam(progressbar = T))

  return(iso.cfm)


}

