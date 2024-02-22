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

    this.iso.df <- x$iso.df
    this.seed.df <- this.iso.df%>%
      dplyr::filter(feature_id==iso_seed)
    this.seed.fid <- this.seed.df$feature_id
    this.smiles <- this.seed.df$smiles
    this.pol <- unique(this.iso.df$polarity)
    this.sp <- x$sp%>%
      split(.,x$sp$ms2_matched_feature)

    ### construct iso-cfm data
    {
      this.iso.cfm <- list()
      this.iso.cfm$compound_info <- list(name = this.seed.df$name,
                                         compound_id = this.seed.df$compound_id,
                                         formula = this.seed.df$formula,
                                         smiles = this.seed.df$smiles,
                                         mz_ref = this.seed.df$mz_ref,
                                         adduct = this.seed.df$adduct,
                                         polarity = this.seed.df$polarity)
      }

    if (is.na(this.smiles)) return(this.iso.cfm)



    ### combine sp of seed to cfm-annotate
    {
      seed.sp <- this.sp[[this.seed.fid]]
      if (length(seed.sp)==0) return(this.iso.cfm)
      if (!all(c(10,20,40) %in% collisionEnergy(seed.sp)))  return(this.iso.cfm)
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
      if (is.null(x$CFM_result )) return(this.iso.cfm)
    }


    ### cfm result
    {
      for (i in 1:nrow(this.iso.df)) {

        mlabel <- paste0("M",this.iso.df$iso_count[i])
        this.iso.cfm[[mlabel]]$Spectra <- this.sp[[this.iso.df$feature_id[i]]]

      }
      this.iso.cfm$M0$CFM_annotation <- x$CFM_result



    }
    return(this.iso.cfm)


  }

  iso.cfm <- bplapply(iso.list,ff,
                      BPPARAM = BPPARAM
  )

  return(iso.cfm)


}

