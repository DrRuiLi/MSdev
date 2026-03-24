#' @title Find isotope labels in MSdev object
#' @description Find isotope labels in MSdev object by processing isotopologues
#' and isotope labels for both positive and negative ion modes.
#'
#' @param object MSdev object
#' @param iso_ele isotope element, default `"[13]C"`
#' @param ppm ppm tolerance, default 10
#' @param ... additional arguments passed to xcms_get_feature_isotopologues
#'
#' @return MSdev object with isotope labels added
#' @export
#'
MSdev_find_isotope_label <- function(object,
                                     iso_ele = "[13]C",
                                     ppm = 10,
                                     ...){

  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.xcms <- xcms_get_feature_isotopologues(xcms.xcms,
                                                iso_ele = iso_ele,
                                                ppm = ppm,
                                                ...)
    xcms.xcms <- xcms_get_feature_isotope_label(xcms.xcms,
                                                iso_ele = iso_ele,
                                                ...)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
  }

  return(object)

}


MSIP_get_isotopologues_table <- function(object,
                                         int_thresh = 5,
                                         ppm=10,
                                         extract_chrom = F){



  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]


    ### MS1 purity
    {
      #sample.idx <- !is.na(pData(xcms.xcms)$isotope_tracer)
      sample.idx <- seq_along(filepaths(xcms.xcms))
      #sample.idx <-  pData(xcms.xcms)$group=="FSLowCFullGlu"
      if ( !"ms1_purity"%in%   colnames(xcms.xcms@msFeatureData$featureDefinitions)) {
        xcms.purity.matrix <- get_xcms_feature_purity_matrix(xcms.xcms ,
                                                             ppm = ppm,
                                                             isolation_half_window = 0.2)
        xcms.xcms <- xcms_get_feature_purity(xcms.xcms,
                                             ms1_purity_matrix =xcms.purity.matrix,
                                             split_source  = T,
                                             selected.sample = sample.idx)
        xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
        xcms.purity.matrix <- xcms.purity.matrix[,pData(xcms.xcms)$sample.type!="Blank"]
        xcms.purity.matrix -> object@statData$MSIP$isotopologues_matrix$ms1_purity[[pol]]
      }

      }

    ### calc ratio of iso and un-iso labeled sample
    {
      xcms.pdata <- pData(xcms.xcms)%>%
        dplyr::filter(!sample.type%in% c("Blank"))
      xcms.fdf <- featureDefinitions(xcms.xcms)%>%as.data.frame()
      xcms.fv <- get_xcms_quantify_MSIP(xcms.xcms)
      xcms.fv <- assay(xcms.fv)[,xcms.pdata$sampleNames,drop = F]
      f.traced <- ifelse(is.na(xcms.pdata$isotope_tracer),
                         "int_mean_nontracer","int_mean_tracer")
      int.mean <- apply(xcms.fv,1,mean_f,f = f.traced,simplify = F)%>%do.call(rbind,.)
      int.mean <- cbind(int.mean,matrix(0,nrow = nrow(int.mean),
                                        ncol =length(setdiff(c("int_mean_nontracer","int_mean_tracer"),colnames(int.mean))),
                               dimnames = list(rownames(int.mean),
                                               setdiff(c("int_mean_nontracer","int_mean_tracer"),colnames(int.mean)))))
      xcms.fdf[,colnames(int.mean)] <- int.mean

      ratio.matrix <- get_xcms_iso_fraction(xcms.xcms)
      ratio.matrix -> object@statData$MSIP$isotopologues_matrix$ratio_to_seed[[pol]]
    }


    ### iso stat and filter, marker features to acq
    {

      xcms.fdf.iso <- xcms.fdf%>%
        dplyr::mutate(is_seed = feature_id%in% iso_seed)
      xcms.fdf.iso[!xcms.fdf.iso$is_seed,
                   c("compound_id","name","adduct","score",
                     "mz_ref","rt_ref","formula","smiles")] <- NA
      #xcms.fdf.iso$ms1_purity[xcms.fdf.iso$is_seed] <- 1
      na.unique <- function(x){
        if (all(is.na(x))) return(NA)
        return(unique(na.omit(x)))
      }
      xcms.fdf.selected <- xcms.fdf.iso%>%
        ### filter not assigned iso
        #dplyr::filter(!is.na(iso_seed))%>%
        dplyr::group_by(iso_seed)%>%
        dplyr::mutate(
          is_isotopologues = !is.na(iso_seed)&any(!is.na(compound_id)),
          formula_max_iso_count = get_formula_ele_count(formula))%>%
        dplyr::mutate(
          formula_max_iso_count = na.unique(formula_max_iso_count),
          iso_seed = case_when(iso_count > formula_max_iso_count~ NA,
                               T~iso_seed),
          iso_count = case_when(iso_count > formula_max_iso_count~ NA,
                               T~iso_count))%>%
        dplyr::ungroup()%>%
        dplyr::group_by(iso_seed)%>%
        dplyr::mutate(
          compound_id = na.unique(compound_id),
          name = na.unique(name),
          adduct = na.unique(adduct),
          formula = na.unique(formula),
          smiles = na.unique(smiles),
          selected_to_acq = T)%>%
        dplyr::mutate(
          selected_to_acq = case_when(
            all(is.na(compound_id))~F,
            int_mean_tracer < int_thresh&int_mean_nontracer<int_thresh ~F,
            ms1_purity < 0.8~ F,
            iso_count == formula_max_iso_count~F,
            !is_labeled&!is_seed~ F,
            T~selected_to_acq
          ),
          seed.acq = is_seed&selected_to_acq,
          selected_to_acq = case_when(
            any(seed.acq)~ selected_to_acq,
            T~F
          ),
          selected_to_acq = case_when(
            sum(selected_to_acq)>1~ selected_to_acq,
            T~F
          ))%>%
        dplyr::ungroup()%>%
        dplyr::arrange(#-int_mean_nontracer,
          iso_seed,
          iso_count)%>%
        dplyr::select(any_of(c("feature_id","mzmed","rtmed","rtmin","rtmax","peakMaxo","polarity",
                             "score","iso_seed","iso_count","ms2_id",
                             "is_seed","is_isotopologues","is_labeled","compound_id","adduct",
                             "name","formula","smiles","int_mean_tracer", "int_mean_nontracer",
                             "ms1_purity","selected_to_acq")),
                      grep(pattern = "Ratio_to_seed",colnames(xcms.fdf),value = T))
    }



    object@statData$MSIP$isotopologues_table[[pol]] <- xcms.fdf.selected


    ### extract chrom
    {
      if (extract_chrom) {
        fid <- xcms.fdf.selected%>%
          dplyr::filter(is_isotopologues)%>%
          dplyr::pull(feature_id)
        message_with_time("Extract ",length(fid)," Chromatograms in ",pol)
        xcms.chrom <- featureChromatograms(xcms.xcms,
                                           features = fid,
                                           expandRt = Inf,
                                           filled = T,
                                           BPPARAM = SnowParam(
                                             workers  = min(snowWorkers(), length(fileNames(xcms.xcms))),
                                             progressbar = T)
        )
        xcms.chrom <- onDiskData(xcms.chrom,
                                 path = paste0(object@projectInfo$projectDir,"/",pol,"_Chromatograms.rds"))
        object@xcmsData[[paste0(pol,"_Chromatograms")]] <- xcms.chrom

      }

    }

  }

  return(object)

}

MSIP_assign_MS2 <- function(object,rt.tol = 10){


  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    sp.ms2 <- object@spectra$MS2_Spectra%>%
      ProtGenerics::filterPolarity(i)
    xcms.fdf <- object@statData$MSIP$isotopologues_table[[pol]]
    sp.ms2.data <- get_Spectra_ms2_feature_id(sp.ms2,
                                              xcms.fdf,
                                              ppm = 10,
                                              rt.tol = rt.tol)
    sum(!is.na(sp.ms2.data$feature_id))

   ### update xcms featuredef
    xcms.fdf$ms2_id <- sapply(xcms.fdf$feature_id,
                              function(i){
                                sp_id <- sp.ms2.data$sp_id[which(sp.ms2.data$feature_id==i)]
                                return(sp_id)
                              }
    )

    xcms.fdf -> object@statData$MSIP$isotopologues_table[[pol]]


  }

  return(object)

}


#' @title Get isotopologues data from MSdev
#' @description Extract iso-labeled compound info and spectra, filter compound
#' based on isotope labeling data.
#'
#' @param object MSdev object
#' @param iso_ele isotope element, default retrieved from object
#'
#' @return a list of isotopologues data
#' @export
#'
MSIP_get_isotopologues_data <- function(object,
                                        iso_ele = get_MSdev_iso_ele(object)){


  sp.ms2 <- onDiskData_retrieve(object@spectra$MS2_Spectra)
  sp.ms2$sample.source <- object@sampleInfo$sample.source[
    match_path(dataOrigin(sp.ms2),object@sampleInfo$msData.files)]
  #cpdb <- CompoundDb::CompDb(object@projectInfo$CompoundDB_path)

  ### Requirement in xcms features
  ### iso_seed
  for (i.pol in 0:1) {

    pol <- ifelse(i.pol==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
    samples <- levels(groupStringFactor(xcms.se$sample.source))
    xcms.fdf <- featureDefinitions(xcms.xcms)
    iso.table <- object@statData$MSIP$isotopologues_table[[pol]]
    iso.table$ms2_id <-xcms.fdf$ms2_id[match(iso.table$feature_id,xcms.fdf$feature_id  )]
    iso.table <- iso.table%>%
      dplyr::mutate(ms2_count = lengths(ms2_id),
                    ele_count = get_formula_ele_count(formula,
                                          ele = get_ele_uniso(iso_ele)))%>%
      dplyr::filter(!is.na(iso_seed),
                    iso_count <= ele_count
                    )%>%
      dplyr::group_by(iso_seed)%>%
      dplyr::filter(#rep(any(is_labeled)),
        any(ms2_count)>0,
                    n() > 1)%>%
      dplyr::ungroup()
    seed.id <- unique(iso.table$iso_seed)
    match(xcms.fdf$iso_seed[which((lengths(xcms.fdf$ms2_id)>0)&!is.na(xcms.fdf$iso_seed))],
          seed.id)
    for (i_seed_id in seq_along(seed.id)) {
      this.fdf <- iso.table%>%
        dplyr::filter(iso_seed == seed.id[i_seed_id])
      this.seed.df <- this.fdf%>%
        dplyr::filter(feature_id==iso_seed)
      this.seed.fid <- this.seed.df$feature_id
      this.list <- list()
      ### compound_info
      {
        this.list$compound_info <- list(name = this.seed.df$name,
                                        compound_id = this.seed.df$compound_id,
                                        mz = this.seed.df$mzmed,
                                        rt = this.seed.df$rtmed,
                                        formula = this.seed.df$formula,
                                        smiles = this.seed.df$smiles,
                                        score = this.seed.df$score,
                                        adduct = this.seed.df$adduct,
                                        polarity = this.seed.df$polarity,
                                        feature_id = setNames(
                                          this.fdf$feature_id,
                                          format_isotopologue(this.fdf$iso_count,"M")),
                                        merged = F)

        atom.count <- get_formula_ele_count(this.list$compound_info$formula,get_ele_uniso(iso_ele))
      }



      ### Spectra split
      {
        this.sp <-  lapply(this.fdf$ms2_id,function(x){
          get_spectra_by_name(sp.ms2,x)
        })
        #message(i_seed_id," ",lengths(this.fdf$ms2_id))
        #this.sp <- lapply(this.sp,combineSpectra_groupby_ce)
        names(this.sp) <- paste0("M", this.fdf$iso_count[ match(names(this.sp),
                                                                this.fdf$feature_id)])
        #this.sp <- this.sp[lengths(this.sp)!=0]
        #if (!length(this.sp)) next
        this.sp <- lapply(this.sp , function(x) split(x,x$sample.source))
        sp_count <- lapply(this.sp,function(y){
          y.l <- lengths(y)
          if (length(y.l)) return(y.l)
          return(0)
        })%>%do.call(rbind,.)
        sp_count <- get_matrix_value_fill_with_NA(sp_count,colnames_vec = colnames(sp_count),rownames_vec = names(this.sp))
        this.list$Spectra <- this.sp
        this.list$compound_info$ms2_count <- sp_count
      }

      ### matrix
      {

        ### intensity_matrix
        {

        }
        ### ratio_matrix
        {
          ratio_matrix <- object@statData$MSIP$
            isotopologues_matrix$ratio_to_seed[[pol]][this.fdf$feature_id,]
          xcms.se.temp <-xcms.se[,colnames(ratio_matrix)]
          ratio_matrix <- apply(ratio_matrix,1,mean_f,
                                f = xcms.se.temp$sample.source,simplify = F)%>%
            do.call(rbind,.)
          rownames(ratio_matrix) <- format_isotopologue(this.fdf$iso_count,"M")
          ratio_matrix <- get_matrix_value_fill_with_NA(ratio_matrix,
                                                        colnames_vec = samples,
                                        rownames_vec = format_isotopologue(0:atom.count,"M"))
        }

        ###purity_matrix
        {
          purity_matrix <- object@statData$MSIP$
            isotopologues_matrix$ms1_purity[[pol]][this.fdf$feature_id,]
          xcms.se.temp <-xcms.se[,colnames(purity_matrix)]
          purity_matrix <- apply(purity_matrix,1,mean_f,
                                 f = xcms.se.temp$sample.source,simplify = F)%>%
            do.call(rbind,.)
          rownames(purity_matrix) <- paste0("M",this.fdf$iso_count)
          purity_matrix <- get_matrix_value_fill_with_NA(purity_matrix,
                                                         colnames_vec = samples,
                                                        rownames_vec = format_isotopologue(0:atom.count,"M"))
        }


        ### natural_matrix
        {

          natural_matrix <- get_iso_natural_ratio(formula = this.list$compound_info$formula,
                                                  iso_ele = iso_ele,ratio_matrix = ratio_matrix)
          natural_matrix <- get_matrix_value_fill_with_NA(natural_matrix,
                                                          colnames_vec = samples,
                                                          rownames_vec = format_isotopologue(0:atom.count,"M"))
        }

        this.list$compound_info$ratio_matrix <- ratio_matrix
        this.list$compound_info$purity_matrix <- purity_matrix
        this.list$compound_info$natural_matrix <- natural_matrix

        }


      ### filter
      {
        #if ( length(this.list$Spectra$M0)==0 ) next
        if (is.na(this.list$compound_info$formula)) next
        if (length(this.list$Spectra)<1) next



      }


      ### MSIPMetaboliteData
      {

        object@statData$MSIP$isotopologues_data[[paste0(seed.id[i_seed_id],"_",pol)]] <-
          MSIPMetaboliteData(CompoundInfo = this.list$compound_info,
                           Spectra = this.list$Spectra     )

      }

    }


  }

  #object@statData$MSIP$isotopologues_data <- iso.list
  return(object)

}



MSIP_get_isotopologues_data_fid <- function(object,fid_seed,polarity,
                                            iso_ele = get_MSdev_iso_ele(object)){





  ### from MSIP_get_isotopologues_data
  {
    sp.ms2 <- onDiskData_retrieve(object@spectra$MS2_Spectra)
    sp.ms2$sample.source <- object@sampleInfo$sample.source[
      match_path(dataOrigin(sp.ms2),object@sampleInfo$msData.files)]
    #cpdb <- CompoundDb::CompDb(object@projectInfo$CompoundDB_path)

    ### Requirement in xcms features
    ### iso_seed

    {
      pol <- ifelse(polarity ==0,"Negative","Positive")
      xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
      xcms.se <- get_xcms_quantify_MSIP(xcms.xcms)
      samples <- levels(groupStringFactor(xcms.se$sample.source))
      xcms.fdf <- featureDefinitions(xcms.xcms)
      iso.table <- object@statData$MSIP$isotopologues_table[[pol]]
      iso.table$ms2_id <-xcms.fdf$ms2_id[match(iso.table$feature_id,xcms.fdf$feature_id  )]
      iso.table <- iso.table%>%
        dplyr::mutate(ms2_count = lengths(ms2_id),
                      ele_count = get_formula_ele_count(formula,
                                                        ele = get_ele_uniso(iso_ele)))%>%
        dplyr::filter(!is.na(iso_seed),
                      iso_count <= ele_count
        )

    }



      {

        this.fdf <- iso.table%>%
          dplyr::filter(iso_seed == fid_seed)
        this.seed.df <- this.fdf%>%
          dplyr::filter(feature_id==iso_seed)
        this.seed.fid <- this.seed.df$feature_id
        this.list <- list()
      }
        ### compound_info
        {
          this.list$compound_info <- list(name = this.seed.df$name,
                                          compound_id = this.seed.df$compound_id,
                                          mz = this.seed.df$mzmed,
                                          rt = this.seed.df$rtmed,
                                          formula = this.seed.df$formula,
                                          smiles = this.seed.df$smiles,
                                          score = this.seed.df$score,
                                          adduct = this.seed.df$adduct,
                                          polarity = this.seed.df$polarity,
                                          merged = F)

          atom.count <- get_formula_ele_count(this.list$compound_info$formula,get_ele_uniso(iso_ele))
        }



        ### Spectra split
        {
          this.sp <-  lapply(this.fdf$ms2_id,function(x){
            get_spectra_by_name(sp.ms2,x)
          })
          #message(i_seed_id," ",lengths(this.fdf$ms2_id))
          names(this.sp) <- paste0("M", this.fdf$iso_count[ match(names(this.sp),
                                                                  this.fdf$feature_id)])

          this.sp <- lapply(this.sp , function(x) split(x,x$sample.source))
          sp_count <- lapply(this.sp,function(y){
            y.l <- lengths(y)
            if (length(y.l)) return(y.l)
            return(0)
          })%>%do.call(rbind,.)

          to.add <- setdiff(unique(xcms.se$sample.source),colnames(sp_count))
          sp_count <- get_matrix_value_fill_with_NA(sp_count,
                                                    colnames_vec = unique(xcms.se$sample.source),
                                                    rownames_vec = names(this.sp))
          sp_count[is.na(sp_count)] <- 0

          this.list$Spectra <- this.sp
          this.list$compound_info$ms2_count <-sp_count
        }

        ### matrix
        {

          ### intensity_matrix
          {

          }
          ### ratio_matrix
          {
            ratio_matrix <- object@statData$MSIP$
              isotopologues_matrix$ratio_to_seed[[pol]][this.fdf$feature_id,]
            xcms.se.temp <-xcms.se[,colnames(ratio_matrix)]
            ratio_matrix <- apply(ratio_matrix,1,mean_f,
                                  f = xcms.se.temp$sample.source,simplify = F)%>%
              do.call(rbind,.)
            rownames(ratio_matrix) <- format_isotopologue(this.fdf$iso_count,"M")
            ratio_matrix <- get_matrix_value_fill_with_NA(ratio_matrix,
                                                          colnames_vec = samples,
                                                          rownames_vec = format_isotopologue(0:atom.count,"M"))
          }

          ###purity_matrix
          {
            purity_matrix <- object@statData$MSIP$
              isotopologues_matrix$ms1_purity[[pol]][this.fdf$feature_id,]
            xcms.se.temp <-xcms.se[,colnames(purity_matrix)]
            purity_matrix <- apply(purity_matrix,1,mean_f,
                                   f = xcms.se.temp$sample.source,simplify = F)%>%
              do.call(rbind,.)
            rownames(purity_matrix) <- paste0("M",this.fdf$iso_count)
            purity_matrix <- get_matrix_value_fill_with_NA(purity_matrix,
                                                           colnames_vec = samples,
                                                           rownames_vec = format_isotopologue(0:atom.count,"M"))
          }


          ### natural_matrix
          {

            natural_matrix <- get_iso_natural_ratio(formula = this.list$compound_info$formula,
                                                    iso_ele = iso_ele,ratio_matrix = ratio_matrix)
            natural_matrix <- get_matrix_value_fill_with_NA(natural_matrix,
                                                            colnames_vec = samples,
                                                            rownames_vec = format_isotopologue(0:atom.count,"M"))
          }

          this.list$compound_info$ratio_matrix <- ratio_matrix
          this.list$compound_info$purity_matrix <- purity_matrix
          this.list$compound_info$natural_matrix <- natural_matrix

        }


        ### filter
        {
          #if ( length(this.list$Spectra$M0)==0 ) next
          if (is.na(this.list$compound_info$formula)) next
          if (length(this.list$Spectra)<1) next



        }


        ### MSIPMetaboliteData
        {

          object@statData$MSIP$isotopologues_data[[paste0(fid_seed,"_",pol)]] <-
            MSIPMetaboliteData(CompoundInfo = this.list$compound_info,
                               Spectra = this.list$Spectra     )

        }

      }

  return(object)



}


get_MSdev_iso_ele <- function(object){
  iso_ele <- object@sampleInfo$isotope_tracer
  iso_ele <- unique(na.omit(iso_ele))
  if(length(iso_ele)>1){
    message("Multiple iso_ele, please check")
    iso_ele <- iso_ele[1]
  }
  if(length(iso_ele)==0){
    message("No iso_ele, using [13]C please check")
    iso_ele <- "[13]C"
  }
  return(iso_ele)
}

#' @title Annotate isotopologues with CFM prediction
#' @description Annotate isotopologues using CFM (Collision-induced dissociation
#' Fragmentation simulator) prediction data.
#'
#' @param object MSdev object
#' @param ppm ppm tolerance for matching, default 20
#' @param check_temp if TRUE, check for existing CFM data in temp directory
#' @param BPPARAM BiocParallel backend parameter
#'
#' @return MSdev object with CFM annotation added to isotopologues
#' @export
#'
MSIP_get_isotopologues_CFM_annotation <- function(object,
                                                  ppm = 20,
                                                  check_temp = T,
                                                  BPPARAM = SnowParam(progressbar = T)){


  ff <- function(x){
    #x <- msdev.M1@statData$MSIP$isotopologues_data[["FT06121_Negative"]]
    ### combine sp of seed to cfm-annotate
    {
      ### process M0 Spectra
      if(F){
        ### M0 just annotated for smiles assign
        seed.sp.c <- do.call(c,unname(x$Spectra$M0))
        seed.sp.c <- Spectra_filter_noise(seed.sp.c)
        seed.sp.c <- normalizeSpectra(sp = seed.sp.c,norm_to = "tic")
        seed.sp.c <- combineSpectra_groupby_ce(seed.sp.c,
                                               ppm = ppm,
                                               minProp = 0.3,
                                               plot = F)
        seed.sp.c <- Spectra_fill_3CE(seed.sp.c)
      }

      if (!is.null(x@CompoundInfo[['CFM_annotation']])) {
        message_with_time("cfmd exist, skip")
        return(x)
      }

      CFM_result <- NA
      try.return <- try({

        CFM_result <- get_CFM_data_from_smiles(
          smiles = x@CompoundInfo$smiles,
          compound_id = x@CompoundInfo$compound_id,
          ppm =  ppm,
          adduct = switch(as.character(x@CompoundInfo$polarity),
                                     "0"="[M-H]-",
                                     "1"="[M+H]+"),
          check_temp = check_temp,
          temp_dir = paste0(object@projectInfo$CompoundDB_path,"_cfmd"))

      })
      if (class(CFM_result) !="CFM_data") return(x)
    }

    x@CompoundInfo$CFM_annotation <-CFM_result
    return(x)
  }


   iso.cfm <- bplapply(object@statData$MSIP$isotopologues_data,
                      ff,
                      BPPARAM = BPPARAM )
  annotated <- sapply(iso.cfm,function(x) "CFM_annotation" %in% names(x@CompoundInfo))
  object@statData$MSIP$isotopologues_data <- iso.cfm[annotated]
  object
}



MSIPMetaboliteData_get_CFM_annotation <- function(MSIPMetaboliteData,
                                                  ppm = 20,
                                                  check_temp = T){

  if (!is.null(MSIPMetaboliteData@CompoundInfo[['CFM_annotation']])) {
    message_with_time("cfmd exist, skip")
    return(x)
  }


  CFM_result <- NA
  try.return <- try({

    CFM_result <- get_CFM_data_from_smiles(
      smiles = MSIPMetaboliteData@CompoundInfo$smiles,
      compound_id = MSIPMetaboliteData@CompoundInfo$compound_id,
      ppm =  ppm,
      adduct = switch(as.character(MSIPMetaboliteData@CompoundInfo$polarity),
                      "0"="[M-H]-",
                      "1"="[M+H]+"),
      check_temp = check_temp,
      temp_dir = paste0(object@projectInfo$CompoundDB_path,"_cfmd"))

  })
  if (class(CFM_result) !="CFM_data") return(MSIPMetaboliteData)

  MSIPMetaboliteData@CompoundInfo$CFM_annotation <-CFM_result
  return(x)
}



get_MSdev_iso_acq_list <- function(object,keep = F){

  acq.list <- object@statData$MSIP$isotopologues_table
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    acq.list.pol <- acq.list[[pol]]%>%
      dplyr::filter(selected_to_acq)%>%
      dplyr::mutate(rtmed = case_when(
        is_seed  ~ rtmed,
        T ~ NA
      ))%>%
      dplyr::group_by(iso_seed)%>%
      dplyr::filter(n()>1&any(is_seed))%>%
      dplyr::mutate(rtmed = na.omit(rtmed),
                    rtmin = rtmed - 10,
                    rtmax = rtmed+10
                    )

    acq.list[[pol]] <- QE_list_2feature_def(acq.list.pol,keep = keep)

  }

  acq.df <- do.call(bind_rows,acq.list)
  acq.df.group <- acq.df%>%
    dplyr::filter(is_seed)%>%
    dplyr::group_by(compound_id)%>%
    dplyr::mutate(rt.diff =diff(range(rtmed)),
                  count = n(),
                  is_paired = case_when(
                    n() >1&rt.diff < 10~T,
                    T~F
                  ))



  return(acq.list)
}


MSIP_export_isotopologues_table <- function(object){

  isotopologues_list <- object@statData$MSIP$isotopologues_table

}


MSIP_export_isotopologues_acquisition_list <- function(object,
                                                       keep = T,
                                                       file_path = "acq.list.xlsx"){

  acq.list <- get_MSdev_iso_acq_list(object,keep = keep)

  message_with_time("Export isotopologues acquisition list to ",normalizePath(file_path))
  xlsx.write.list(acq.list,file = file_path)



}


get_iso_net_assign <- function(iso.ig,net.degree.ratio = 0.6){

  iso.fdf <- vdata(iso.ig)
  dis.con <-  as_adjacency_matrix(iso.ig,sparse	 =F ,type = "upper")
  dis.mw <-  distances(iso.ig,mode = "out",
                       weights = edata(iso.ig)$closest.iso.count)
  if (nrow(iso.fdf)<=2) {
    cl <- rep(1,nrow(iso.fdf))
  }else
    cl <- fpc::pamk(dis.con,krange = 1:(nrow(iso.fdf)-1))$pamobject$clustering
  iso.seed <- rep(NA,nrow(iso.fdf))
  names(iso.seed)<-iso.fdf$name
  iso_count<- iso.seed
  for (i.cl in unique(cl)) {
    this.igraph.sub <- igraph_filter_vertex(iso.ig,cl==i.cl)

    igraph::distances(this.igraph.sub,mode = "out",
                      weights = edge.attributes(this.igraph.sub)$closest.iso.count )
    to.delete <- degree(this.igraph.sub)<(length(V(this.igraph.sub))-1)*2*net.degree.ratio
    this.igraph.sub <- igraph::delete.vertices(this.igraph.sub,to.delete )
    #visNetwork::visIgraph(this.igraph.sub)
    this.dis <- igraph::distances(this.igraph.sub,mode = "out",
                                  weights = edge.attributes(this.igraph.sub)$closest.iso.count )
    this.dis[this.dis<0] <- 0
    dis.sum <- apply(this.dis,1,sum)
    seed.fid <- names(which.max(dis.sum))
    dis.to.seed <- this.dis[names(which.max(dis.sum)),]
    iso.seed[names(dis.to.seed)] <- seed.fid
    iso_count[names(dis.to.seed)] <- dis.to.seed
  }

  x <- list(iso.seed=iso.seed,iso_count=iso_count)
  return(x)

}



#' @title Process isotopologues spectra
#' @description Process isotopologues spectra by normalizing and combining
#' spectra for each collision energy.
#'
#' @param iso.list list of isotopologues data
#'
#' @return processed iso.list with normalized and combined spectra
#' @export
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

iso_prob <- function(iso_count=1,C_count = 5){

  M0 <- dbinom(0,C_count,0.0107)
  Mn <- dbinom(iso_count,C_count,0.0107)
  Mn/M0

}


get_isotopomer_atom_matrix <- function(atom.count = 10,
                                       labeled.count = NULL ) {


  m <- atom.count

  if (!is.null(labeled.count)) {

    combos <- combn(atom.count, labeled.count)
    total <- ncol(combos)

    i <- rep(1:total, each = labeled.count )
    j <- as.vector(combos)
    mat <- matrix(0L, nrow = total, ncol = atom.count)
    for (i in seq_len(total)) {
      mat[i, combos[, i]] <- 1L
    }


  }else{

    n <- 2^m
    mat <- matrix(1L, nrow = n, ncol = m)
    for (i in 1:m) {
      repeat_size <- 2^(m - i)
      x <- rep(rep(c(0L, 1L), each = repeat_size), times = 2^(i - 1))
      mat[, i] <- x
    }
  }

  #size_of(mat)
  rownames(mat) <- apply(mat, 1, paste0, collapse = "")

  return(mat)

}

MSIP_solve_isotopologues <- function(object,
                                     process.info = MSIP_solve_computation_evaluate(object,F),
                                     int_thresh = 10^3.8,
                                     certainty_thresh = 0.6,
                                     weight_fun = .intensity_weight,
                            ppm = 20,
                            timeout = 60,
                            BPPARAM = SerialParam(progressbar = T)){



  iso.data <- object@statData$MSIP$isotopologues_data
  .f <- function(i){

    this.fid <- process.info$feature_id[i]
    this.iso.count <- process.info$iso_count[i]
    this.sample <- process.info$samples[i]
    this.sp.count <- process.info$ms2.count[i]


    start_time <- Sys.time()
    message_with_time(this.fid,";",this.sample,";",str_isotope2_num(this.iso.count),";",this.sp.count," Spectra")

    this.natural.ratio <- process.info$natural.ratio[i]
    this.msip.mtbd <- iso.data[[this.fid]]
    this.cfmd <- this.msip.mtbd@CompoundInfo$CFM_annotation
    this.sp.iso <-this.msip.mtbd@Spectra[[str_isotope2_num(this.iso.count)]][[this.sample]]
    #if(!length(this.sp.iso)) return(NULL)
    msip.core <- get_MSIPCoreData(sp.iso = this.sp.iso,
                                  cfmd = this.cfmd,
                                  iso_count_max = this.iso.count,
                                  ppm = ppm)
    msip.core <- MSIPCore_solve(msip.core,
                                int_thresh = int_thresh,
                                weight_fun=weight_fun,
                                certainty_thresh = certainty_thresh)
    if (process.info$traced[i]&!is.na(this.natural.ratio)){
      msip.core <- MSIPCore_correct_natural(msip.core, natural.ratio = this.natural.ratio)
    }
    time_consume <- Sys.time()-start_time
    MSIPIsotopomerMap <- msip.core@solve$MSIPIsotopomerMap
    if (is.null(MSIPIsotopomerMap )){
      n_isotopomers <- 0
    }else{
      n_isotopomers <- length(MSIPIsotopomerMap@isotopomer.defination)
    }

    message_with_time("total of ",n_isotopomers," isotopomers, time consume: ",
                      format(time_consume,digits = 4))

    return(msip.core)


  }
  msip.result.list <- bplapply(1:nrow(process.info),
                               FUN = function(i){
                                 R.utils::withTimeout(.f(i),
                                                      timeout = timeout,
                                                      onTimeout = "silent")
                               },
                 BPPARAM =BPPARAM)


  for (i in 1:nrow(process.info)) {
    this.fid <- process.info$feature_id[i]
    this.iso.count <- process.info$iso_count[i]%>%str_isotope2_num()
    this.sample<- process.info$samples[i]
    if (i < length(msip.result.list)) {
      iso.data[[this.fid]]@MSIPIsotopologueDatas[[this.iso.count]][[this.sample]] <-
        msip.result.list[[i]]
    }
  }

  iso.data -> object@statData$MSIP$isotopologues_data
  return(object)


}


MSIP_drop_isotopologues_tempdata <- function(object){

  isotopologues_data <- object@statData$MSIP$isotopologues_data
  for (i_fid in names(isotopologues_data)) {

    this.iso.data <-isotopologues_data[[i_fid]]
    for (i_iso_count in names(this.iso.data$MSIP_result)) {
      this.sample.data <-this.iso.data$MSIP_result[[i_iso_count]]
      for (i_sample in names(this.sample.data)) {
        this.msip.core.data <- this.sample.data[[i_sample]]
        this.msip.core.data <- MSIPCore_drop(this.msip.core.data)
        this.msip.core.data -> this.sample.data[[i_sample]]
      }
      this.sample.data -> this.iso.data$MSIP_result[[i_iso_count]]

    }
    this.iso.data -> isotopologues_data[[i_fid]]

  }
  isotopologues_data -> object@statData$MSIP$isotopologues_data
  return(object)

}
.get_MSIP_tracer <- function(object){

  x <- object@sampleInfo%>%
    dplyr::filter(!sample.type %in% c("Blank"))%>%
    dplyr::select(sample.source,isotope_tracer)%>%
    dplyr::distinct()%>%
    dplyr::pull(isotope_tracer,name = sample.source)

  return(x)

}




get_iso_natural_ratio <- function(formula, iso_ele, ratio_matrix ){

  thresh = min(ratio_matrix,na.rm = T)
  if (thresh < 1e-5)  thresh<- 1e-3
  iso_pattern <- MSCC::chemform_isotopes_pattern_enviPat(
    formula,
    thresh = thresh)%>%
    dplyr::filter(grepl(iso_ele,isotope_element,fixed = T)|isotope_element=="")%>%
    dplyr::mutate(iso_count = sub(pattern = iso_ele,
                                  x = isotope_element,
                                  replacement = "",fixed = T),
                  iso_count = case_when(iso_count==""~"0",
                                        T~iso_count),
                  iso_count = paste0("M",iso_count))%>%
    dplyr::filter(iso_count%in% rownames(ratio_matrix))%>%
    dplyr::add_row(iso_count = setdiff(rownames(ratio_matrix),.$iso_count),
                   abundance = 0)
  iso_pattern$abundance/100/ratio_matrix[iso_pattern$iso_count,,drop=F]

}



MSIP_build_compoundDB <- function(cpdb = CompDb("C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb"),
                                  compound_id,
                                  db_path){

  #db_path <- tempfile()
  cpt <- CompoundDb::compounds(cpdb,
                               compoundVariables(cpdb,
                                                 includeId = T))
  cpsp <- CompoundDb::Spectra(cpdb)
  cpt <- cpt[cpt$compound_id %in% compound_id,,drop = F]
  cpsp <- cpsp[cpsp$compound_id %in% compound_id]

  cpdb_new <- CompoundDb::emptyCompDb(db_path)
  cpdb_new <- insertCompound(cpdb_new,cpt)
  cpdb_new <- insertSpectra(cpdb_new,cpsp,
                            columns = setdiff(spectraVariables(cpsp),
                                              "synonym"))

  cpdb_new <- CompDb(db_path)
  return(cpdb_new)
}

MSIP_update_compoundDB_from_interest_list <- function(){

  roi <- readxl::read_excel("C:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/MSIP.interest.list.xlsx")
  MSIP_build_compoundDB(compound_id = roi$compound_id,
                        db_path = "c:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/MSDB.ROI.sqlite")

}



get_MSIP_solve_computation_evaluate <- function(object,
                                            include.merged = F,
                                            show_message = F){

  iso.data <- object@statData$MSIP$isotopologues_data
  iso_ele <- get_MSdev_iso_ele(object)
  target_ele <- get_ele_uniso(iso_ele)
  all.sample <- .get_MSIP_tracer(object)%>%names()
  #traced.sample <- names(na.omit(all.sample))


  comp.eval.list <- list()
  for (i in seq_along(iso.data)) {

    this.msip.mtbd <- iso.data[[i]]

    ### consider pos and neg
    cfmd <-this.msip.mtbd@CompoundInfo[[grep("CFM",names(this.msip.mtbd@CompoundInfo),value = T)[1]]]

    if (is.null(cfmd)) next
    cfmd.ig <- get_cfm_data_sdf_igraph(cfmd)
    this.atom <- get_sdf_igraph_atom(cfmd.ig,ele = target_ele)
    this.ele.count <-length(this.atom)
    iso_count <- names(this.msip.mtbd@Spectra)%>%
      format_isotopologue("num")


    natural.ratio.matrix <- this.msip.mtbd@CompoundInfo$natural_matrix
    ms2_count.matrix <- this.msip.mtbd@CompoundInfo$ms2_count
    comp.eval <- expand.grid(
      feature_id = names(iso.data)[i],
      name = this.msip.mtbd@CompoundInfo$name,
      compound_id = this.msip.mtbd@CompoundInfo$compound_id,
      iso_count = iso_count,
      target_ele_count =this.ele.count,
      samples = all.sample,
      stringsAsFactors =F,
      merged = this.msip.mtbd@CompoundInfo$merged
    )%>%
      dplyr::mutate(isotopomer = choose(target_ele_count ,iso_count ))%>%
      dplyr::rowwise()%>%
      dplyr::mutate(natural.ratio =
                      get_matrix_value_fill_with_NA(natural.ratio.matrix,
                                                    str_isotope2_num(iso_count),
                                                    paste0(samples)),
                    ms2.count =
                      get_matrix_value_fill_with_NA(ms2_count.matrix,
                                                    str_isotope2_num(iso_count),
                                                    paste0(samples)),
                    solved = F,
                    FSIS.count = NA,
                    sp.consistency.icc= NA,
                    sp.consistency.cos= NA)%>%
      dplyr::ungroup()%>%
      dplyr::filter(!is.na(ms2.count))
  for (j in 1:nrow(comp.eval)) {

    msip.core <- this.msip.mtbd@MSIPIsotopologueDatas[[str_isotope2_num(comp.eval$iso_count[j])]][[
      comp.eval$samples[j]]]
    comp.eval$solved[j] <-!is.null(msip.core)


    if (comp.eval$solved[j]) {
      comp.eval$sp.consistency.icc[j] <- mean(msip.core@FG_map@FG.data$icc,na.rm = T)
      comp.eval$sp.consistency.cos[j] <- mean(msip.core@FG_map@FG.data$cos,na.rm = T)

      if(!isEmpty(msip.core@solve$MSIPIsotopomerMap)){
        #length(msip.core@solve$MSIPIsotopomerMap@)
        comp.eval$FSIS.count[j] <-length(msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set)
      }
    }
  }
  if (show_message) {

    message(names(iso.data)[i],", Total ",this.ele.count," ", target_ele)
    mes <- paste0("M",comp.eval$iso_count,", isoforms: ",
                  comp.eval$isotopomer)%>%
      paste0(collapse = "\n")
    message(mes)
    message("")
  }




    comp.eval.list[[i]] <-comp.eval

  }

  sample.tracer <- .get_MSIP_tracer(object)
  comp.eval <- do.call(rbind,comp.eval.list)%>%
    #dplyr::filter(ms2.count>0)%>%
    dplyr::mutate(traced = case_when(
      is.na(sample.tracer[samples] )~F,
      T~T))
  if(include.merged){
    comp.eval <- comp.eval
  }else{
    comp.eval <- comp.eval%>%
      dplyr::filter(!merged)
  }

  return(invisible( comp.eval ))

}




MSIP_from_Spectra <- function(object,
                              ppm = 5,
                              peak_width = 30){


  sp.ms2 <- onDiskData_retrieve(object@spectra$MS2_Spectra)
  polarity.index <- c("0" = "Negative",
                      "1" = "Positive")
  for (i in 0:1) {
    polarity.index[i]
    sp.pol <- filterPolarity(sp.ms2,i)
    sample.info.pol <- object@sampleInfo%>%
      dplyr::filter(polarity == i)

    xcms.pol <- xcms_from_ms2_spectra(sp.pol,
                                      sample.info = sample.info.pol,
                                      ppm = ppm,peak_width =peak_width )
    xcms.pol <- xcms_get_feature_isotopologues(xcms.pol,
                                               iso_ele = get_MSdev_iso_ele(object),
                                               ppm = ppm,
                                               rt.tol = peak_width)




  }



}



#' @title Merge positive and negative isotopologues data
#' @description Merge MSIP data from positive and negative ion modes for the same compound.
#' @describeIn MSIP merge msip data of positive and negative
#' @param object MSIP/MSdev object
#' @param rt.tol permitted rt error to merge pos and neg data, default 5 seconds
#'
#' @return MSIP/MSdev object with merged isotopologues
#' @export
#'
MSIP_merge_isotopologues <- function(object,
                       rt.tol = 5){

  isotopologues.data.list <- object@statData$MSIP$isotopologues_data

  isotopologues.to.merge <- get_MSIP_compound_info(object,vars = c("all"))%>%
    dplyr::mutate(feature_id = names(isotopologues.data.list),
                  rt = as.numeric(rt),
                  rt.group = MsCoreUtils::group(rt,
                                                tolerance = rt.tol),
                  cp.group = paste0(compound_id ,"_", rt.group)
                  )%>%
    dplyr::filter(polarity %in%0:1)%>%
    dplyr::group_by(cp.group,rt.group)%>%
    dplyr::filter( any(polarity==0) & any(polarity==1))%>%
    dplyr::mutate(rt.diff = abs(rt - mean(rt)))%>%
    dplyr::group_by(polarity,.add = T)%>%
    dplyr::slice_min(rt.diff)

  to.merge <- unique(isotopologues.to.merge$cp.group)
  for (i in seq_along(to.merge)) {
    this.cp.group <-to.merge[i]
    this.isotopologues <- isotopologues.to.merge %>%
      dplyr::filter(cp.group==this.cp.group)

    iso.data1 <- isotopologues.data.list[[this.isotopologues$feature_id[1]]]
    iso.data2 <- isotopologues.data.list[[this.isotopologues$feature_id[2]]]

    iso.data <- merge_isotopologues_data(iso.data1,iso.data2)
    isotopologues.data.list[[paste0(unique(this.cp.group),"_merged")]]<-
      iso.data

  }


  isotopologues.data.list -> object@statData$MSIP$isotopologues_data

  return(object)

}




merge_isotopologues_data <- function(iso.data1,iso.data2){

  iso.data <- iso.data1

  ### compound.info
  {

    iso.data$compound_info$mz <- NA
    iso.data$compound_info$rt <- mean(iso.data1$compound_info$rt,iso.data2$compound_info$rt )
    iso.data$compound_info$score <- max(iso.data1$compound_info$score ,
                                        iso.data2$compound_info$score )
    iso.data$compound_info$adduct <-paste0(iso.data1$compound_info$adduct,";",
                                           iso.data2$compound_info$adduct)
    iso.data$compound_info$polarity <- paste0(iso.data1$compound_info$polarity,";",
                                              iso.data2$compound_info$polarity)
    iso.data$compound_info$merged <- T
    iso.data$compound_info$ms2_count <- add_matrix(iso.data1$compound_info$ms2_count,
                                                   iso.data2$compound_info$ms2_count)
    iso.data$compound_info$ratio_matrix <- mean_matrix(iso.data1$compound_info$ratio_matrix,
                                                   iso.data2$compound_info$ratio_matrix)
    iso.data$compound_info$purity_matrix <- mean_matrix(iso.data1$compound_info$purity_matrix,
                                                   iso.data2$compound_info$purity_matrix)
    iso.data$compound_info$natural_matrix <- mean_matrix(iso.data1$compound_info$natural_matrix,
                                                       iso.data2$compound_info$natural_matrix)


  }


  ### Spectra
  {
    ms2_count <- iso.data$compound_info$ms2_count
    for (i in rownames(ms2_count)) {
      for (j in colnames(ms2_count)) {
        #message(i,j)
        iso.data$Spectra[[i]][[j]] <- c(
          iso.data1$Spectra[[i]][[j]],
          iso.data2$Spectra[[i]][[j]]
        )
      }
    }

  }


  ### cfmd
  {
    iso.data[[paste0("CFM_annotation",iso.data1$compound_info$polarity)]] <-
      iso.data1$CFM_annotation
    iso.data[[paste0("CFM_annotation",iso.data2$compound_info$polarity)]] <-
      iso.data2$CFM_annotation
    iso.data$CFM_annotation <- NULL
  }


  ### msipcore
  {
    iso.data$MSIP_result <- list()
    for (i in rownames(ms2_count)) {
      for (j in colnames(ms2_count)) {

        MSIPCoreData1 <- iso.data1$MSIP_result[[i]][[j]]
        MSIPCoreData2 <- iso.data2$MSIP_result[[i]][[j]]
        MSIPCoreData <- MSIPCore_merge(MSIPCoreData1,
                                             MSIPCoreData2,
                                             suffix1 = iso.data1$compound_info$polarity,
                                             suffix2 = iso.data2$compound_info$polarity)
        if (!is.null(MSIPCoreData)) {
          MSIPCoreData <- MSIPCore_solve(MSIPCoreData)
        }
        iso.data$MSIP_result[[i]][[j]] <- MSIPCoreData

        }
      }


  }


  ###
  return(iso.data)

}



get_MSIP_M1_Statistic <- function(object){


  non.tracer.sample <- .get_MSIP_tracer(object)%>%
    is.na()%>%
    which()%>%
    names()

  process.info <- MSIP_solve_computation_evaluate(object,
                                                  include.merged = T)
  process.info.m1 <- process.info%>%
    dplyr::filter(iso_count == 1,
                  samples %in%non.tracer.sample,
                  solved)
  process.info.m1$r2 <- NA
  process.info.m1$rmse <- NA
  process.info.m1$is.count <- NA
  iso.list <- object@statData$MSIP$isotopologues_data
  for (i in 1:nrow(process.info.m1)) {

    msip.core <- iso.list[[process.info.m1$feature_id[i]]]$
      MSIP_result[[str_isotope2_num(process.info.m1$iso_count[i])]][[process.info.m1$samples[i]]]
    if (is.null(msip.core)|isEmpty(msip.core)) {
      next
    }
    msip.ifmap <- msip.core@solve$MSIPIsotopomerMap
    x <- lengths(msip.ifmap@solve$isotopomer.set)
    y <- msip.ifmap@solve$isotopomer.set.prob
    if (length(x)) {
      x <- x/sum(x)
      process.info.m1$is.count[i] <- length(x)
      process.info.m1$r2[i]<-  summary(lm(x~y))$r.squared
      process.info.m1$rmse[i] <- sqrt(mean((x - y) ^ 2))
    }
  }


  return(process.info.m1)



}


get_MSIP_fragment_Statistic <- function(object){

  process.info <- MSIP_solve_computation_evaluate(object,
                                                  include.merged = T)
  process.info.fg <- process.info%>%
    dplyr::filter(iso_count>0,solved)
  iso.list <- object@statData$MSIP$isotopologues_data
  process.info.fg$fragment.count <- NA
  process.info.fg$is.count<- NA
  for (i in 1:nrow(process.info.fg)) {

    msip.core <- iso.list[[process.info.fg$feature_id[i]]]$
      MSIP_result[[str_isotope2_num(process.info.fg$iso_count[i])]][[process.info.fg$samples[i]]]
    if (is.null(msip.core)|isEmpty(msip.core)) {
      next
    }
    msip.map <- msip.core@FG_map
    msip.ifmap<- msip.core@solve$MSIPIsotopomerMap
    process.info.fg$fragment.count[i]  <- sum(msip.map@fragment.include)
    process.info.fg$is.count[i] <- length( lengths(msip.ifmap@solve$isotopomer.set))

  }

  return(process.info.fg)

}


get_MSIP_result_Statistic <- function(object,
                                      include.merged = F,
                                      show_message = F){

  iso.data <- object@statData$MSIP$isotopologues_data
  iso_ele <- get_MSdev_iso_ele(object)
  target_ele <- get_ele_uniso(iso_ele)
  all.sample <- .get_MSIP_tracer(object)%>%names()
  #traced.sample <- names(na.omit(all.sample))


  comp.eval.list <- list()
  for (i in seq_along(iso.data)) {

    cfmd <- iso.data[[i]][[grep("CFM",names(iso.data[[i]]),value = T)[1]]]
    if (is.null(cfmd)) next
    cfmd.ig <- get_cfm_data_sdf_igraph(cfmd)
    this.atom <- get_sdf_igraph_atom(cfmd.ig,ele = target_ele)
    this.ele.count <-length(this.atom)
    iso_count <- names(iso.data[[i]]$Spectra)%>%
      str_isotope2_num()


    natural.ratio.matrix <- this.msip.mtbd@CompoundInfo$natural_matrix
    ms2_count.matrix <- this.msip.mtbd@CompoundInfo$ms2_count
    comp.eval <- expand.grid(
      feature_id = names(iso.data)[i],
      name = this.msip.mtbd@CompoundInfo$name,
      compound_id = this.msip.mtbd@CompoundInfo$compound_id,
      iso_count = iso_count,
      target_ele_count =this.ele.count,
      samples = all.sample,
      stringsAsFactors =F,
      merged = this.msip.mtbd@CompoundInfo$merged
    )%>%
      dplyr::mutate(isotopomer = choose(target_ele_count ,iso_count ))%>%
      dplyr::rowwise()%>%
      dplyr::mutate(natural.ratio =
                      get_matrix_value_fill_with_NA(natural.ratio.matrix,
                                                    str_isotope2_num(iso_count),
                                                    paste0(samples)),
                    ms2.count =
                      get_matrix_value_fill_with_NA(ms2_count.matrix,
                                                    str_isotope2_num(iso_count),
                                                    paste0(samples)),
                    solved = F,
                    FSIS.count = NA)%>%
      dplyr::ungroup()%>%
      dplyr::filter(!is.na(ms2.count))
    for (j in 1:nrow(comp.eval)) {
      msip.core <- iso.data[[i]][["MSIP_result"]][[str_isotope2_num(comp.eval$iso_count[j])]][[
        comp.eval$samples[j]]]
      comp.eval$solved[j] <-!is.null(msip.core)
      comp.eval$FSIS.exist <- NA
      comp.eval$isotopomer.exist <- NA
      if (comp.eval$solved[j]) {
        if(!isEmpty(msip.core@solve$MSIPIsotopomerMap)){
          #length(msip.core@solve$MSIPIsotopomerMap@)
          comp.eval$FSIS.count[j] <-length(msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set)
          comp.eval$FSIS.exist[j]  <- sum(msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob>0.00001)
          comp.eval$isotopomer.exist[j] <- sum(msip.core@solve$MSIPIsotopomerMap@isotopomer.probability>0.00001)
        }
      }
    }
    if (show_message) {

      message(names(iso.data)[i],", Total ",this.ele.count," ", target_ele)
      mes <- paste0("M",comp.eval$iso_count,", isoforms: ",
                    comp.eval$isotopomer)%>%
        paste0(collapse = "\n")
      message(mes)
      message("")
    }




    comp.eval.list[[i]] <-comp.eval

  }

  sample.tracer <- .get_MSIP_tracer(object)
  comp.eval <- do.call(rbind,comp.eval.list)%>%
    dplyr::filter(ms2.count>0)%>%
    dplyr::mutate(traced = case_when(
      is.na(sample.tracer[samples] )~F,
      T~T))
  if(include.merged){
    comp.eval <- comp.eval
  }else{
    comp.eval <- comp.eval%>%
      dplyr::filter(!merged)
  }

  return(invisible( comp.eval ))

}



get_MSIP_weight_fun <- function(object){

  msmodel <- object@projectInfo$msModel
  if ("Orbitrap Astral" %in% msmodel) {
    return(.intensity_weight_astral)
  }
  return(.intensity_weight)

}




get_MSIP_Molecule_igraph <- function(object,fraction_thresh = 0.001){

  isotopologues_datas <- object@statData$MSIP$isotopologues_data
  dm <- list(names(isotopologues_datas),unique(object@sampleInfo$sample.source))
  Molecule_igraph_matrix <- matrix(list(),
                                   nrow = lengths(dm)[1],
                                   ncol =lengths(dm)[2],dimnames = dm
  )
  target_ele <- get_MSdev_iso_ele(object)
  for (i in seq_along(isotopologues_datas)) { ### Metabolites
    isotopologues_data <- isotopologues_datas[[i]]
    ratio_matrix <- isotopologues_data@CompoundInfo$ratio_matrix
    natural_matrix <- isotopologues_data@CompoundInfo$natural_matrix
    ratio_adj_matrix <- ratio_matrix*(1-natural_matrix)
    ratio_adj_matrix["M0",] <- 1
    ratio_adj_matrix[ratio_adj_matrix<0] <- 0
    ratio_adj_matrix[is.na(ratio_adj_matrix)] <- 0
    for (j in colnames(Molecule_igraph_matrix)  ) {### samples

      mol.ig <- get_Molecule_igraph_from_smiles(isotopologues_data@CompoundInfo$smiles)
      for (k in rownames(ratio_matrix) ) { ### isotopologues

        message_with_time("i ",i," j ",j," k ",k)
        if (k=="M0") next
        msip.core <- isotopologues_data@MSIPIsotopologueDatas[[k]][[j]]
        if (isEmpty(msip.core)) {

          ### if MSIP not data exist, simulate isotopologue
          mol.ig <- Molecule_igraph_add_isotopologue(
            Molecule_igraph = mol.ig,
            isotopologue = k,
            abundance = ratio_adj_matrix[k,j],
            target_ele =target_ele,
            all_isotopomers = F

          )

        }else{
          isotopomers.ele <- msip.core@solve$MSIPIsotopomerMap@isotopomer.defination
          istp <- rownames(isotopomers.ele)
          isotopomers.prob <- msip.core@solve$MSIPIsotopomerMap@isotopomer.probability[istp]
          isotopomers.FSIS <- msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set%>%
            lapply(function(x){data.frame(isotopomer = x)})%>%
            data.table::rbindlist(idcol = "FSIS")%>%
            dplyr::mutate(FSIS = paste0(k,"_",FSIS))%>%
            dplyr::pull(FSIS,name = isotopomer)
          isotopomers.FSIS <- isotopomers.FSIS[istp]
          isotopomers.abundance <- (isotopomers.prob * ratio_adj_matrix[k,j])[istp]
          ### fliter
          {

            idx.filt <- which(isotopomers.prob > fraction_thresh)
            isotopomers.ele <- isotopomers.ele[idx.filt,,drop = F]
            isotopomers.FSIS <- isotopomers.FSIS[idx.filt]
            isotopomers.abundance <- isotopomers.abundance[idx.filt]

          }
          isotopomers.ele.list <- apply(isotopomers.ele,1,function(x){
            make_vector(target_ele,
                        colnames(isotopomers.ele)[x==1])},simplify = F)

          for (ii in seq_len(nrow(isotopomers.ele))) {
            #message(ii)
            mol.ig <- Molecule_igraph_add_isotopomer(
              mol.ig,
              isotopomer = rownames(isotopomers.ele)[ii],
              iso_vec = isotopomers.ele.list[[ii]],
              FSIS = isotopomers.FSIS[ii],
              abundance = isotopomers.abundance[ii])
          }


        }


      }

      mol.ig <- Molecule_igraph_filter_isotopomers(mol.ig,fraction_thresh = fraction_thresh)

      Molecule_igraph_matrix[i,j] <- list(mol.ig)
    }

  }


  return(Molecule_igraph_matrix)


}

get_MSIPIsotopologueData_Molecule_igraphs  <- function(object,fraction_thresh = 0.001){


  { ### Metabolites
    isotopologues_data <- object
    ratio_matrix <- isotopologues_data@CompoundInfo$ratio_matrix
    natural_matrix <- isotopologues_data@CompoundInfo$natural_matrix
    ratio_adj_matrix <- ratio_matrix*(1-natural_matrix)
    ratio_adj_matrix["M0",] <- 1
    ratio_adj_matrix[ratio_adj_matrix<0] <- 0
    ratio_adj_matrix[is.na(ratio_adj_matrix)] <- 0
    for (j in colnames(Molecule_igraph_matrix)  ) {### samples

      mol.ig <- get_Molecule_igraph_from_smiles(isotopologues_data@CompoundInfo$smiles)
      for (k in rownames(ratio_matrix) ) { ### isotopologues

        message_with_time(" j ",j," k ",k)
        if (k=="M0") next
        msip.core <- isotopologues_data@MSIPIsotopologueDatas[[k]][[j]]
        if (isEmpty(msip.core)) {

          ### if MSIP not data exist, simulate isotopologue
          mol.ig <- Molecule_igraph_add_isotopologue(
            Molecule_igraph = mol.ig,
            isotopologue = k,
            abundance = ratio_adj_matrix[k,j],
            target_ele =target_ele,
            all_isotopomers = F

          )

        }else{
          isotopomers.ele <- msip.core@solve$MSIPIsotopomerMap@isotopomer.defination
          istp <- rownames(isotopomers.ele)
          isotopomers.prob <- msip.core@solve$MSIPIsotopomerMap@isotopomer.probability[istp]
          isotopomers.FSIS <- msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set%>%
            lapply(function(x){data.frame(isotopomer = x)})%>%
            data.table::rbindlist(idcol = "FSIS")%>%
            dplyr::mutate(FSIS = paste0(k,"_",FSIS))%>%
            dplyr::pull(FSIS,name = isotopomer)
          isotopomers.FSIS <- isotopomers.FSIS[istp]
          isotopomers.abundance <- (isotopomers.prob * ratio_adj_matrix[k,j])[istp]
          ### fliter
          {

            idx.filt <- which(isotopomers.prob > fraction_thresh)
            isotopomers.ele <- isotopomers.ele[idx.filt,,drop = F]
            isotopomers.FSIS <- isotopomers.FSIS[idx.filt]
            isotopomers.abundance <- isotopomers.abundance[idx.filt]

          }
          isotopomers.ele.list <- apply(isotopomers.ele,1,function(x){
            make_vector(target_ele,
                        colnames(isotopomers.ele)[x==1])},simplify = F)

          for (ii in seq_len(nrow(isotopomers.ele))) {
            #message(ii)
            mol.ig <- Molecule_igraph_add_isotopomer(
              mol.ig,
              isotopomer = rownames(isotopomers.ele)[ii],
              iso_vec = isotopomers.ele.list[[ii]],
              FSIS = isotopomers.FSIS[ii],
              abundance = isotopomers.abundance[ii])
          }


        }


      }

      mol.ig <- Molecule_igraph_filter_isotopomers(mol.ig,fraction_thresh = fraction_thresh)


    }

  }
}



MSIP_clear_previous_data <- function(object){

  object@statData$MSIP$isotopologues_data <- list()
  return(object)
}



Report_MSIP_raw_data <- function(object,
                                 file = paste0(object@projectInfo$projectDir,"/MSIP.raw.data.pdf"),
                        show_chrom = F){



  ### get all chroms
  {
   if (show_chrom) {
     chroms.list <- object@xcmsData[c("Negative_Chromatograms",
                                      "Positive_Chromatograms")]%>%
       lapply(function(x){

         x <- onDiskData_retrieve(x)
         # pol <- ifelse(polarity(x)==0,"Negative","Positive")
         # rownames(x) <- paste0( featureDefinitions(x)$feature_id,"_", pol)
         return(x)
       })

   }

   # chroms <- rbind(chroms.list[[1]],
   #                 chroms.list[[2]])


  }

  comp.eval <- get_MSIP_solve_computation_evaluate(object )

  comp.eval <- comp.eval%>%
    dplyr::filter(solved,iso_count > 0)



  ### out put
  {

    file.dir <- paste0(object@projectInfo$projectDir,"/MSIP/raw.data")
    dir.create(file.dir,showWarnings = F,recursive = T)
    open_dir(file.dir)
    for (i in 1:nrow(comp.eval)) {

      ### retrieve data
      {

        this.fid <- comp.eval$feature_id[i]
        this.sample <- comp.eval$samples[i]
        this.isotopologue <- comp.eval$iso_count[i]

        this.data <- object@statData$MSIP$isotopologues_data[[this.fid]]


        this.cfmd <- this.data@CompoundInfo$CFM_annotation
        this.msip.core <- this.data@MSIPIsotopologueDatas[[format_isotopologue(this.isotopologue,"M")]][[this.sample]]
        this.sp <- this.msip.core@Spectra_data
        this.mz <- this.data@CompoundInfo$mz + 1.003355484 * this.isotopologue


        this.cp.name <- this.data@CompoundInfo$name
        this.title <- paste0(this.cp.name,", ",
                             format_isotopologue(this.isotopologue,"+"),", ",
                             this.sample)

        this.subtitle <- paste0(this.data@CompoundInfo$adduct,",",
                                "mz = ",str_digit(this.mz,5))


        if (isEmpty(this.msip.core)) next

      }

      {
        file.path <- paste0(file.dir,"/",this.fid,".pdf")
      }

      ### Chrom + MS2 Acq
      {
        if(show_chrom){
          this.chroms <- chroms.list[[this.data@CompoundInfo$polarity+1]]
          this.chrom <- this.chroms[match(sub(x= this.fid,pattern = "_Negative|_Positive",""),
                                          featureDefinitions(this.chroms)$feature_id),]



          p.chrom <- plot_XChromatograms(this.chrom,
                                         color_f = pData(this.chrom)$sample.source,
                                         norm = F,move = F)+
            geom_vline(xintercept = this.data@CompoundInfo$rt,
                       linewidth = 5,
                       col = "grey",alpha = 0.5)





        }else{
          p.chrom <- ggplot()
        }


        this.sp.temp <- Spectra_filter_cfm_annotated(this.sp)%>%
          Spectra_get_totIonCurrent()
        ms2.acq <- Spectra::spectraData(this.sp.temp   )%>%
          as.data.frame()

        p.chrom+
          geom_line(data = ms2.acq,
                    aes(x = rtime ,
                        y = totIonCurrent ,
                        group = collisionEnergy,
                    ),col = "grey")+
          geom_point(data = ms2.acq,
                     aes(x = rtime ,
                         y = totIonCurrent ,
                         fill = collisionEnergy),
                     pch = 21)+
          xlim(this.data@CompoundInfo$rt+c(-60,60))+
          scale_fill_gradient(low = "orange", high = "red")+
          labs(title = this.title,subtitle = this.subtitle)+
          scale_y_log10()+
          theme(legend.position = "nonde")-> p.chrom.acq


      }


      ### Spectra
      {

        sp.temp <- filterSpectraIntensity(this.sp,ratio = 0.01)%>%
          Spectra::filterIntensity(intensity = 1e3)
        p.sp.ce <- plot_Spectra_CE(sp.temp)+
          theme(legend.position = "none")

      }

      ### Fragment
      {

        p.sp.fg <- plot_MSIPCore_spectra_consistency_hm(this.msip.core)



      }

      ### Isotopomer result
      {
        p.isotopomer <- plot_MSIPCore_result(this.msip.core)+
          labs(title = this.title,subtitle = this.subtitle)


      }


      ### Molecule igraph
      {
        this.mig <- get_Molecule_igraph_from_smiles(this.data@CompoundInfo$smiles)
        this.mig <- this.cfmd@fragment_sdf@SDF[[1]]%>%get_Molecule_igraph_from_sdf()
        p.mig <- plot_Molecule_igraph(this.mig,show_id = T,size = 2)

      }


      ### export
      {


        p1 <- p.chrom.acq/p.sp.ce
        #p1 <- p.sp.ce

        export_graph2pdf(p1,file_path = file.path,append = T,width = 5,height = 6)

        export_graph2pdf(p.sp.fg,file_path = file.path,append = T,width = ncol(p.sp.fg@matrix)*0.2+2,height = 6)

        p3 <- p.mig+p.isotopomer+plot_layout(widths = c(2,1))

        export_graph2pdf(p3,file_path = file.path,append = T,width = 15,height = 6)
      }
    }
  }


}



Report_MSIP_isotopomers <- function(object,
                                    file = paste0(object@projectInfo$projectDir,"/MSIP.isotopomer.pdf")){


  MSIP.mol.igs <- object@statData$MSIP$isotopomer_Molecule_igraph
  suppressWarnings(file.remove(file))
  open_dir(dirname(file))

  for (i in 1:nrow(MSIP.mol.igs)) {

    ### info
    {

      fid <- rownames(MSIP.mol.igs)[i]
      message_with_time(fid)
      this.msip.mtblt <-  object@statData$MSIP$isotopologues_data[[fid]]
      cp <- this.msip.mtblt@CompoundInfo$name

      this.info  <-
        paste0(fid,": ",cp,"\n",
               "Formula: ",object@statData$MSIP$isotopologues_data[[fid]]@CompoundInfo$formula,"\n",
               "RT: ",str_digit(object@statData$MSIP$isotopologues_data[[fid]]@CompoundInfo$rt,0),"\n",
               "mz: ",str_digit(object@statData$MSIP$isotopologues_data[[fid]]@CompoundInfo$mz,4),"\n"
               )
      p.info <- ggplot()+
        geom_text(aes(x=0,y=1,label = this.info),vjust = 1,hjust = 0)+
        xlim(c(0,1))+
        ylim(c(0,1))+
        theme_void()

    }

    ### composition of every sample
    {
      mol.igs <- MSIP.mol.igs[i,]
      p.cirs <- list()
      for (i.mig in seq_along(mol.igs)) {
        mig <- mol.igs[[i.mig]]
        p.cirs[[i.mig]] <- plot_Molecular_igraph_isotopomer_circle(mig)+
          labs(title = names(mol.igs)[i.mig])+
          theme(plot.title = element_text(hjust = 0.5))
      }
      p.per.cirs <- ggplot_sum_patchwork(p.cirs)+
        patchwork::plot_layout(nrow = 1)+
        patchwork::plot_annotation(tag_levels  =NULL)
    }

    ### compare of samples
    {
      p.bar <- plot_Molecular_igraphs_isotopomer_bar(mol.igs)
      if (is.null(p.bar)) next
      #p.bar/p.per.cirs
    }


    ### Molecular structure
    {
      p.mig <- plot_Molecule_igraph(mol.igs[[1]],show_id = T,size = 2)
    }



    ### export
    {

      export_graph2pdf(p.info+ p.mig +p.bar,file , width = 20,height = 5,append = T)
      export_graph2pdf(p.per.cirs ,file , width = 20,height = 6,append = T)


    }



  }

  return(invisible())

}


MSIP_get_Molecule_igraph <- function(object,fraction_thresh = 0.001){

  MSIP.mol.igs <- get_MSIP_Molecule_igraph(object,fraction_thresh = fraction_thresh)
  object@statData$MSIP$isotopomer_Molecule_igraph <-MSIP.mol.igs
  return(object)

}
