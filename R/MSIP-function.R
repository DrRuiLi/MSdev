#' @title Import compound table into MSdev
#' @description Read a compound table from file (.xlsx or .csv) and store it in obj@advancedAna$MSIP$compound_table.
#' The table must contain columns: compound_id, name, formula, smiles, rt.
#'
#' @param obj MSdev object
#' @param table.path file path to the compound table (.xlsx or .csv)
#'
#' @return MSdev object with compound table stored in advancedAna
#' @export
#'
MSIP_import_compound_table <- function(obj, table.path) {

  ext <- tolower(tools::file_ext(table.path))

  if (ext == "xlsx") {
    compound_table <- openxlsx::read.xlsx(table.path)
  } else if (ext == "csv") {
    compound_table <- read.csv(table.path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported format: .", ext, ". Accepted formats: .xlsx, .csv")
  }

  required_cols <- c("compound_id", "name", "formula", "smiles", "rt")
  missing_cols <- setdiff(required_cols, colnames(compound_table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  obj@advancedAna$MSIP$compound_table <- compound_table
  return(obj)
}


#' @title Find isotopologues from compound table
#' @description Find isotopologue features for compounds listed in the compound table.
#' Calculates expected m/z for each isotope label combination using formula and isotope mass differences,
#' then matches to xcms features by m/z and rt. Stores result in 'obj@advancedAna$MSIP$isotopologues_table'.
#'
#' @details
#' \lifecycle{deprecated}
#' This function is deprecated. Use \code{MSIP_xcms_processing.targeted()} instead,
#' which handles isotopologue annotation and M+0 injection automatically.
#'
#' @param obj MSdev object with compound table in advancedAna$MSIP$compound_table
#' @param iso_ele isotope element, default "\[13\]C"
#' @param ppm ppm tolerance for m/z matching, default 10
#' @param rt.tol rt tolerance in seconds, default 10
#'
#' @return MSdev object with isotopologues_table stored in advancedAna
#' @export
#'
MSIP_find_isotopologue_from_compound_table <- function(obj,
                                                        iso_ele = "[13]C",
                                                        ppm = 10,
                                                        rt.tol = 10) {
  .Deprecated("MSIP_xcms_processing.targeted",
              package = "MSdev",
              msg = "MSIP_find_isotopologue_from_compound_table is deprecated. \\
Use MSIP_xcms_processing.targeted() which handles isotopologue annotation automatically.")

  compound_table <- obj@advancedAna$MSIP$compound_table
  if (is.null(compound_table)) {
    stop("compound_table not found. Run MSIP_import_compound_table first.")
  }

  for (i in 0:1) {
    pol <- ifelse(i == 0, "Negative", "Positive")
    xcms.xcms <- obj@xcmsData[[paste0(pol, "MS1")]]
    if (is.null(xcms.xcms) || identical(xcms.xcms, NA)) next

    xcms.fdf <- featureDefinitions(xcms.xcms) %>% as.data.frame()
    polarity_val <- i
    target_ele <- get_ele_uniso(iso_ele)

    result_list <- list()
    for (j in seq_len(nrow(compound_table))) {
      cp <- compound_table[j, ]
      cp_formula <- cp$formula

      adduct <- ifelse(i == 1, "[M+H]+", "[M-H]-")
      cp_adduct <- MSCC::chemform_adduct(cp_formula, adduct, value = "all")
      if (is.null(cp_adduct) || nrow(cp_adduct) == 0) next
      seed_mz <- cp_adduct$chemform.adduct.mz[1]

      cp_max_iso <- get_formula_ele_count(cp_formula, target_ele)
      if (is.na(cp_max_iso) || cp_max_iso == 0) next

      ele_counts <- setNames(list(cp_max_iso), target_ele)
      iso_grid <- do.call(MSCC::get_isotope_mass_diff, ele_counts)
      iso_grid$iso_form <- iso_grid$chemform_diff
      iso_grid$iso_mz <- seed_mz + iso_grid$mass_diff

      if (!is.null(cp$rt) && !is.na(cp$rt)) {
        match_res <- match_mz_rt(mz1 = iso_grid$iso_mz,
                                  rt1 = cp$rt,
                                  mz2 = xcms.fdf$mzmed,
                                  rt2 = xcms.fdf$rtmed,
                                  mz.ppm = ppm,
                                  rt.tol = 60)
      } else {
        match_res <- match_mz_rt(mz1 = iso_grid$iso_mz,
                                  mz2 = xcms.fdf$mzmed,
                                  mz.ppm = ppm)
      }

      if (nrow(match_res) == 0) next

      matched_features <- match_res %>%
        dplyr::mutate(iso_grid[ion1,],
                      xcms.fdf[ion2, ],
                      rt.cl = cluster_rt(rtmed,rt.tol = rt.tol))%>%
        dplyr::group_by(rt.cl)%>%
        dplyr::mutate(rt.cl.error = mean(rt.error))%>%
        dplyr::ungroup()%>%
        dplyr::filter(rt.cl.error == min(rt.cl.error))

      seed_fid <- matched_features$feature_id[which.min(matched_features$mass_diff)]
      iso_df <- matched_features %>%
        dplyr::mutate(
          iso_seed = seed_fid,
          is_seed = feature_id == iso_seed,
          is_isotopologues = TRUE,
          compound_id = cp$compound_id,
          name = cp$name,
          formula = cp_formula,
          smiles = cp$smiles,
          rt_ref = ifelse(is.null(cp$rt), NA, cp$rt)
        )

      result_list[[length(result_list) + 1]] <- iso_df
    }

    if (length(result_list) > 0) {
      xcms.fdf.selected <- do.call(rbind, result_list)
      xcms.fdf.selected <- xcms.fdf.selected %>%
        dplyr::arrange(compound_id, iso_form)

      obj@advancedAna$MSIP$isotopologues_table[[pol]] <- xcms.fdf.selected
    }
  }

  return(obj)
}


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
        xcms.purity.matrix -> object@advancedAna$MSIP$isotopologues_matrix$ms1_purity[[pol]]
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
      ratio.matrix -> object@advancedAna$MSIP$isotopologues_matrix$ratio_to_seed[[pol]]
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



    object@advancedAna$MSIP$isotopologues_table[[pol]] <- xcms.fdf.selected


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
                                             workers  = min(snowWorkers(), length(MSnbase::fileNames(xcms.xcms))),
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
    xcms.fdf <- object@advancedAna$MSIP$isotopologues_table[[pol]]
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

    xcms.fdf -> object@advancedAna$MSIP$isotopologues_table[[pol]]


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
#' @return a list of isotopologues data (legacy)
#' @keywords internal
#'
.MSIP_get_isotopologues_data_legacy <- function(object,
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
    iso.table <- object@advancedAna$MSIP$isotopologues_table[[pol]]
    iso.table$ms2_id <-xcms.fdf$ms2_id[match(iso.table$feature_id,xcms.fdf$feature_id  )]
    iso.table <- iso.table%>%
      dplyr::mutate(ms2_count = lengths(ms2_id),
                    ele_count = get_formula_ele_count(formula,
                                          ele = get_ele_uniso(iso_ele)))%>%
      dplyr::filter(!is.na(iso_seed),
                    #iso_count <= ele_count
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
          ratio_matrix <- object@advancedAna$MSIP$
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
          purity_matrix <- object@advancedAna$MSIP$
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

        object@advancedAna$MSIP$isotopologues_data[[paste0(seed.id[i_seed_id],"_",pol)]] <-
          MSIPMetaboliteData(CompoundInfo = this.list$compound_info,
                           Spectra = this.list$Spectra     )

      }

    }


  }

  #object@advancedAna$MSIP$isotopologues_data <- iso.list
  return(object)

}


#' @title Build MSIP isotopologue-level SummarizedExperiment
#' @description
#' Construct \code{\link{MSIPIsotopologueData}} objects (one per compound) and return
#' them as a named list. This is the updated
#' isotopologue-level container: a \code{SummarizedExperiment} with assays
#' \code{intensity}, \code{ratio} and \code{purity} across \code{sample.source}.
#'
#' The function uses existing xcms/MSIP outputs:
#' - xcms feature definitions (\code{featureDefinitions(object@xcmsData[[...]] )})
#'   for feature ↔ compound ↔ isotopologue mapping
#' - optional matrices in \code{object@advancedAna$MSIP$isotopologues_matrix}:
#'   \code{ratio_to_seed} and \code{ms1_purity}
#' - intensities quantified from \code{object@xcmsData} via \code{get_xcms_quantify_MSIP()}.
#'
#' @param object MSdev object.
#' @param compound_id optional character vector. If NULL, build for all compounds in \code{compound_table}.
#' @param iso_ele isotope element, default from \code{get_MSdev_iso_ele(object)}.
#' @param assay_fun function used to aggregate replicates within the same \code{sample.source};
#'   default \code{mean}.
#'
#' @return named list of \code{MSIPIsotopologueData}, one element per \code{compound_id}.
#' @export
get_MSIPIsotopologueData <- function(object,
                                     compound_id = NULL,
                                     iso_ele = get_MSdev_iso_ele(object),
                                     assay_fun = mean) {
  msip <- object@advancedAna[["MSIP"]]
  if (is.null(msip)) {
    stop("object@advancedAna[['MSIP']] is missing.")
  }

  .extract_compound_table <- function(x, max_iter = 12L) {
    need <- c("compound_id", "name", "formula", "smiles", "rt")
    cur <- x
    for (i in seq_len(max_iter)) {
      if (is.data.frame(cur) && all(need %in% colnames(cur))) return(cur)
      if (!is.list(cur)) break
      if ("compound_table" %in% names(cur)) {
        cur <- cur[["compound_table"]]
        next
      }
      if ("MSIP" %in% names(cur)) {
        nxt <- cur[["MSIP"]]
        if (identical(nxt, cur)) break
        cur <- nxt
        next
      }
      break
    }
    NULL
  }

  compound_table <- .extract_compound_table(msip)
  if (is.null(compound_table)) {
    stop("compound_table not found or invalid in object@advancedAna$MSIP. ",
         "Expected columns: compound_id, name, formula, smiles, rt.")
  }
  if (is.null(compound_id)) compound_id <- as.character(compound_table$compound_id)
  compound_id <- unique(as.character(compound_id))

  # helper: aggregate per-sample columns to sample.source
  .agg_to_source <- function(mat, sources, fun = mean) {
    if (is.null(mat) || !nrow(mat) || !ncol(mat)) return(mat)
    sources <- as.character(sources)
    u <- unique(sources)
    out <- sapply(u, function(ss) {
      idx <- which(sources %in% ss)
      if (!length(idx)) return(rep(NA_real_, nrow(mat)))
      apply(mat[, idx, drop = FALSE], 1, fun, na.rm = TRUE)
    })
    out <- as.matrix(out)
    rownames(out) <- rownames(mat)
    colnames(out) <- u
    out
  }

  # helper: pull matrices for one polarity (0/1)
  .get_pol_mats <- function(ion_mode) {
    pol <- ifelse(ion_mode == 0, "Negative", "Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]
    if (is.null(xcms.xcms) || identical(xcms.xcms, NA)) return(NULL)

    se <- tryCatch(
      get_xcms_quantify_MSIP(xcms.xcms),
      error = function(e) NULL
    )
    if (is.null(se)) return(NULL)
    int <- SummarizedExperiment::assay(se)
    colnames(int) <- Biobase::pData(xcms.xcms)$sample.name
    sources <- Biobase::pData(xcms.xcms)$sample.source
    groups <- Biobase::pData(xcms.xcms)$group

    int.src <- .agg_to_source(int, sources, fun = assay_fun)

    # ratio: compute isotopologue fraction (ratio to seed) per feature
    ratio.src <- NULL
    ratio.matrix <- tryCatch(get_xcms_iso_fraction(xcms.xcms), error = function(e) NULL)
    if (!is.null(ratio.matrix)) {
      pdata <- Biobase::pData(xcms.xcms)
      ratio.col <- colnames(ratio.matrix)
      src_ratio <- pdata$sample.source[match(ratio.col, pdata$sample.name)]
      if (all(is.na(src_ratio)) && "sampleNames" %in% colnames(pdata)) {
        src_ratio <- pdata$sample.source[match(ratio.col, pdata$sampleNames)]
      }
      nm_ratio <- pdata$sample.name[match(ratio.col, pdata$sample.name)]
      nm_ratio[is.na(nm_ratio)] <- ratio.col[is.na(nm_ratio)]
      colnames(ratio.matrix) <- nm_ratio
      ratio.src <- .agg_to_source(ratio.matrix, src_ratio, fun = assay_fun)
    }

    # purity: compute MS1 purity matrix (feature-by-sample) from MS1 Spectra, then aggregate
    purity.src <- NULL

    # If MS1 purity is not already available, try to compute from MS1 spectra.
    if (is.null(purity.src)) {
      ms1.sp <- NULL
      if (is.null(object@spectra$MS1_Spectra)) {
        message("[get_MSIPIsotopologueData] ", pol,
                ": object@spectra$MS1_Spectra is NULL; purity will be NA.")
      } else {
        ms1.sp <- tryCatch(onDiskData_retrieve(object@spectra$MS1_Spectra), error = function(e) {
          message("[get_MSIPIsotopologueData] ", pol,
                  ": failed to retrieve MS1 Spectra: ", conditionMessage(e))
          NULL
        })
      }
      if (!is.null(ms1.sp) && inherits(ms1.sp, "Spectra")) {
        ms1.sp <- tryCatch(ProtGenerics::filterPolarity(ms1.sp, ion_mode), error = function(e) ms1.sp)
        if (!length(ms1.sp)) {
          message("[get_MSIPIsotopologueData] ", pol,
                  ": MS1 Spectra has 0 spectra after polarity filter; purity will be NA.")
          ms1.sp <- NULL
        }
      }
      purity.matrix <- tryCatch(
        get_xcms_feature_purity_matrix(
          xcms.xcms,
          xcms.ms1.sp = ms1.sp,
          ppm = 10,
          isolation_half_window = 0.2
        ),
        error = function(e) {
          message("[get_MSIPIsotopologueData] ", pol,
                  ": get_xcms_feature_purity_matrix() failed: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(purity.matrix)) {
        purity.matrix <- purity.matrix[rownames(int), , drop = FALSE]
        pdata <- Biobase::pData(xcms.xcms)
        pur.col <- colnames(purity.matrix)
        src_pur <- pdata$sample.source[match_path(pur.col, pdata$sampleNames)]
        nm_pur <- pdata$sample.name[match_path(pur.col, pdata$sampleNames)]
        nm_pur[is.na(nm_pur)] <- pur.col[is.na(nm_pur)]
        colnames(purity.matrix) <- nm_pur
        purity.src <- .agg_to_source(purity.matrix, src_pur, fun = assay_fun)
      }
    }

    # colData by sample.source
    cda <- data.frame(sample.source = colnames(int.src), stringsAsFactors = FALSE)
    cda$group <- vapply(cda$sample.source, function(ss) {
      g <- unique(as.character(groups[sources %in% ss]))
      g <- g[!is.na(g) & nzchar(g)]
      if (!length(g)) return(NA_character_)
      paste(g, collapse = ";")
    }, character(1))
    rownames(cda) <- cda$sample.source

    list(intensity = int.src, ratio = ratio.src, purity = purity.src, colData = cda)
  }

  # helper: derive iso_count robustly from featureDefinitions
  .derive_iso_count <- function(df, iso_ele = "[13]C") {
    n <- nrow(df)
    iso_count <- rep(NA_integer_, n)
    if ("iso_count" %in% colnames(df)) {
      iso_count <- suppressWarnings(as.integer(df$iso_count))
    }
    if (all(is.na(iso_count)) && "isotopologue_form" %in% colnames(df)) {
      x <- as.character(df$isotopologue_form)
      num <- suppressWarnings(as.integer(gsub(".*?([0-9]+)$", "\\1", x)))
      num[is.na(num) & !is.na(x) & nzchar(x)] <- 0L
      iso_count <- num
    }
    if (all(is.na(iso_count)) && "isotope_element" %in% colnames(df)) {
      x <- as.character(df$isotope_element)
      x0 <- ifelse(is.na(x), "", x)
      x0 <- gsub(iso_ele, "", x0, fixed = TRUE)
      x0 <- gsub("[^0-9]", "", x0)
      num <- suppressWarnings(as.integer(x0))
      num[is.na(num) & !is.na(df$isotope_element)] <- 0L
      iso_count <- num
    }
    if ("feature_id" %in% colnames(df) && "iso_seed" %in% colnames(df)) {
      is_seed <- !is.na(df$iso_seed) & as.character(df$feature_id) == as.character(df$iso_seed)
      iso_count[is_seed & is.na(iso_count)] <- 0L
    }
    iso_count
  }

  # helper: build compound/isotopologue map from xcms featureDefinitions
  .get_iso_map_from_fdf <- function(ion_mode) {
    pol <- ifelse(ion_mode == 0, "Negative", "Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]
    if (is.null(xcms.xcms) || identical(xcms.xcms, NA)) return(NULL)
    fdf <- xcms::featureDefinitions(xcms.xcms) %>% as.data.frame()
    if (!nrow(fdf) || !("feature_id" %in% colnames(fdf))) return(NULL)

    # diagnostic: report when featureDefinitions lacks isotopologue-matching variables
    match_vars <- c("compound_id", "name", "iso_seed",
                    "isotopologue_form", "isotope_element")
    if (!any(match_vars %in% colnames(fdf))) {
      message("[get_MSIPIsotopologueData] ", pol,
              ": featureDefinitions lacks isotopologue-matching variables ",
              "(compound_id/name/iso_seed/isotopologue_form/isotope_element). ",
              "No isotopologue mapping can be derived for this polarity.")
    }

    if (!("compound_id" %in% colnames(fdf))) fdf$compound_id <- NA_character_
    if (!("name" %in% colnames(fdf))) fdf$name <- NA_character_
    if (!("iso_seed" %in% colnames(fdf))) fdf$iso_seed <- NA_character_

    fdf$iso_count <- .derive_iso_count(fdf, iso_ele = iso_ele)

    # propagate compound/name within iso_seed group when only seed has annotation
    if (any(!is.na(fdf$iso_seed))) {
      fdf <- fdf %>%
        dplyr::group_by(iso_seed) %>%
        dplyr::mutate(
          compound_id = {
            x <- unique(stats::na.omit(as.character(compound_id)))
            if (length(x)) x[[1]] else NA_character_
          },
          name = {
            x <- unique(stats::na.omit(as.character(name)))
            if (length(x)) x[[1]] else NA_character_
          }
        ) %>%
        dplyr::ungroup()
    }

    fdf <- fdf %>%
      dplyr::mutate(
        compound_id = as.character(compound_id),
        name = as.character(name),
        feature_id = as.character(feature_id),
        iso_count = suppressWarnings(as.integer(iso_count)),
        polarity = pol
      ) %>%
      dplyr::filter(!is.na(compound_id), !is.na(feature_id), !is.na(iso_count)) %>%
      dplyr::select(feature_id, compound_id, name, iso_count, polarity) %>%
      dplyr::distinct()
    if (!nrow(fdf)) return(NULL)
    fdf
  }

  .build_from_legacy_nested <- function(object, iso_ele = "[13]C") {
    # Legacy payload observed in old saved objects:
    # (A) msip$isotopologue_data[[1]]$compound_table + $isotopomer_data
    # (B) msip$compound_table + msip$isotopomer_data
    legacy_candidates <- list(
      msip[["isotopologue_data"]],
      msip
    )
    cpt <- NULL
    iso_mat_list <- NULL
    for (cand in legacy_candidates) {
      if (is.null(cand)) next
      if (is.list(cand) && all(c("compound_table", "isotopomer_data") %in% names(cand))) {
        cpt <- cand[["compound_table"]]
        iso_mat_list <- cand[["isotopomer_data"]]
      } else if (is.list(cand) && length(cand) && is.list(cand[[1]]) &&
                 all(c("compound_table", "isotopomer_data") %in% names(cand[[1]]))) {
        cpt <- cand[[1]][["compound_table"]]
        iso_mat_list <- cand[[1]][["isotopomer_data"]]
      }
      if (is.data.frame(cpt) && is.list(iso_mat_list)) break
    }
    if (!is.data.frame(cpt) || !is.list(iso_mat_list)) return(NULL)
    if (!all(c("compound_id", "name") %in% colnames(cpt))) return(NULL)

    # sample.source -> group from sampleInfo
    sample.info <- as.data.frame(object@sampleInfo)
    src_group <- NULL
    if (all(c("sample.source", "group") %in% colnames(sample.info))) {
      src_group <- tapply(as.character(sample.info$group),
                          as.character(sample.info$sample.source),
                          function(z) paste(unique(stats::na.omit(z)), collapse = ";"))
    }

    out <- list()
    cids <- intersect(as.character(cpt$compound_id), names(iso_mat_list))
    for (cid in cids) {
      m <- iso_mat_list[[cid]]
      if (is.null(dim(m)) || is.null(rownames(m)) || is.null(colnames(m))) next
      iso_form <- rownames(m)
      sample_src <- colnames(m)
      iso_count <- suppressWarnings(as.integer(gsub(".*?([0-9]+)$", "\\1", iso_form)))
      iso_count[is.na(iso_count)] <- seq_along(iso_form) - 1L
      iso_id <- paste0(cid, "_", iso_form)
      cp_name <- cpt$name[match(cid, cpt$compound_id)]
      cp_name <- ifelse(length(cp_name) && !is.na(cp_name), as.character(cp_name[[1]]), NA_character_)

      rda <- S4Vectors::DataFrame(
        isotopologue_form = iso_form,
        isotopologue_id = iso_id,
        compound_id = rep(cid, length(iso_form)),
        compound_name = rep(cp_name, length(iso_form)),
        label.isotopologue = ifelse(is.na(cp_name), NA_character_, paste0(cp_name, "_M+", iso_count)),
        row.names = iso_id
      )
      cda <- data.frame(
        sample.source = sample_src,
        group = if (is.null(src_group)) rep(NA_character_, length(sample_src)) else as.character(src_group[sample_src]),
        row.names = sample_src,
        stringsAsFactors = FALSE
      )

      .msipcore_intensity <- function(x) {
        # Return a single numeric intensity proxy (TIC sum) from legacy MSIPCoreData-like objects.
        if (is.null(x)) return(NA_real_)
        # Many saved objects store MSIPCoreData cells as data.frame already.
        if (methods::is(x, "MSIPCoreData")) {
          sp <- x@Spectra
          if (inherits(sp, "Spectra")) {
            ints <- Spectra::intensity(sp)
            return(sum(vapply(ints, function(v) sum(v, na.rm = TRUE), numeric(1)), na.rm = TRUE))
          }
          if (is.data.frame(sp) && all(c("sp_id", "intensity") %in% colnames(sp))) {
            return(sum(as.numeric(sp$intensity), na.rm = TRUE))
          }
          return(NA_real_)
        }
        if (inherits(x, "Spectra")) {
          ints <- Spectra::intensity(x)
          return(sum(vapply(ints, function(v) sum(v, na.rm = TRUE), numeric(1)), na.rm = TRUE))
        }
        if (is.data.frame(x)) {
          if ("totIonCurrent" %in% colnames(x)) return(sum(as.numeric(x$totIonCurrent), na.rm = TRUE))
          if ("intensity" %in% colnames(x)) return(sum(as.numeric(x$intensity), na.rm = TRUE))
          return(NA_real_)
        }
        if (is.list(x) && length(x) == 1) {
          return(.msipcore_intensity(x[[1]]))
        }
        NA_real_
      }

      intensity.mat <- matrix(NA_real_, nrow = length(iso_form), ncol = length(sample_src),
                              dimnames = list(iso_id, sample_src))
      purity.mat <- intensity.mat
      for (i_row in seq_along(iso_form)) {
        for (i_col in seq_along(sample_src)) {
          cell <- m[i_row, i_col][[1]]
          intensity.mat[i_row, i_col] <- .msipcore_intensity(cell)
        }
      }

      ratio.mat <- intensity.mat
      if (nrow(intensity.mat) >= 1) {
        denom <- intensity.mat[1, ]
        ratio.mat <- sweep(intensity.mat, 2, denom, `/`)
      }
      out[[cid]] <- MSIPIsotopologueData(
        assays = list(
          intensity.positive = intensity.mat,
          intensity.negative = purity.mat,
          ratio.positive = ratio.mat,
          ratio.negative = purity.mat,
          purity.positive = purity.mat,
          purity.negative = purity.mat
        ),
        rowData = rda,
        colData = S4Vectors::DataFrame(cda)
      )
    }
    if (!length(out)) return(NULL)
    out
  }

  # collect polarity mapping + matrices
  iso_all <- list()
  pol_mats <- list()
  for (ion_mode in 0:1) {
    pol <- ifelse(ion_mode == 0, "Negative", "Positive")
    iso_all[[pol]] <- .get_iso_map_from_fdf(ion_mode)
    pol_mats[[pol]] <- .get_pol_mats(ion_mode)
  }
  iso_all_df <- do.call(rbind, iso_all)
  if (is.null(iso_all_df) || !nrow(iso_all_df)) {
    legacy.out <- .build_from_legacy_nested(object, iso_ele = iso_ele)
    if (!is.null(legacy.out)) {
      return(legacy.out)
    }
    stop("No isotopologue mapping found in xcms featureDefinitions. ",
         "Check that features are annotated with compound_id and isotopologue info.")
  }

  # ensure required columns exist in derived map
  need <- c("compound_id", "name", "iso_count", "feature_id", "polarity")
  miss <- setdiff(need, colnames(iso_all_df))
  if (length(miss)) {
    stop("isotopologue mapping missing required columns: ", paste(miss, collapse = ", "))
  }
  iso_all_df$compound_id <- as.character(iso_all_df$compound_id)
  iso_all_df$name <- as.character(iso_all_df$name)
  iso_all_df$iso_count <- suppressWarnings(as.integer(iso_all_df$iso_count))

  # build per compound
  out <- list()
  for (cid in compound_id) {
    this.df <- iso_all_df[iso_all_df$compound_id %in% cid, , drop = FALSE]
    if (!nrow(this.df)) next

    cp.name <- unique(na.omit(this.df$name))
    cp.name <- if (length(cp.name)) cp.name[[1]] else NA_character_

    max_iso <- max(this.df$iso_count, na.rm = TRUE)
    if (!is.finite(max_iso) || is.na(max_iso)) next
    iso_counts <- 0:max_iso
    iso_form <- paste0(iso_ele, iso_counts)
    iso_id <- paste0(cid, "_", iso_form)
    label_iso <- if (!is.na(cp.name)) paste0(cp.name, "_M+", iso_counts) else rep(NA_character_, length(iso_counts))

    # union sample.sources across available polarity matrices
    cda <- NULL
    all_sources <- character()
    for (pol in c("Negative", "Positive")) {
      pm <- pol_mats[[pol]]
      if (is.null(pm)) next
      all_sources <- union(all_sources, colnames(pm$intensity))
      if (is.null(cda)) cda <- pm$colData
    }
    all_sources <- unique(all_sources)
    if (!length(all_sources)) next
    if (is.null(cda)) {
      cda <- data.frame(sample.source = all_sources, group = NA_character_, row.names = all_sources)
    } else {
      cda <- cda[all_sources, , drop = FALSE]
    }

    intensity.pos <- matrix(NA_real_, nrow = length(iso_counts), ncol = length(all_sources),
                            dimnames = list(iso_id, all_sources))
    intensity.neg <- intensity.pos
    ratio.pos <- intensity.pos
    ratio.neg <- intensity.pos
    purity.pos <- intensity.pos
    purity.neg <- intensity.pos

    for (k in seq_along(iso_counts)) {
      ic <- iso_counts[[k]]

      # Negative
      fids.neg <- this.df$feature_id[this.df$iso_count %in% ic & this.df$polarity %in% "Negative"]
      fids.neg <- unique(na.omit(as.character(fids.neg)))
      if (length(fids.neg) && !is.null(pol_mats$Negative)) {
        pm <- pol_mats$Negative
        fi <- intersect(fids.neg, rownames(pm$intensity))
        if (length(fi)) {
          intensity.neg[k, ] <- apply(pm$intensity[fi, , drop = FALSE], 2, mean, na.rm = TRUE)
        }
        if (!is.null(pm$ratio)) {
          fi <- intersect(fids.neg, rownames(pm$ratio))
          if (length(fi)) ratio.neg[k, ] <- apply(pm$ratio[fi, , drop = FALSE], 2, mean, na.rm = TRUE)
        }
        if (!is.null(pm$purity)) {
          fi <- intersect(fids.neg, rownames(pm$purity))
          if (length(fi)) purity.neg[k, ] <- apply(pm$purity[fi, , drop = FALSE], 2, mean, na.rm = TRUE)
        }
      }

      # Positive
      fids.pos <- this.df$feature_id[this.df$iso_count %in% ic & this.df$polarity %in% "Positive"]
      fids.pos <- unique(na.omit(as.character(fids.pos)))
      if (length(fids.pos) && !is.null(pol_mats$Positive)) {
        pm <- pol_mats$Positive
        fi <- intersect(fids.pos, rownames(pm$intensity))
        if (length(fi)) {
          intensity.pos[k, ] <- apply(pm$intensity[fi, , drop = FALSE], 2, mean, na.rm = TRUE)
        }
        if (!is.null(pm$ratio)) {
          fi <- intersect(fids.pos, rownames(pm$ratio))
          if (length(fi)) ratio.pos[k, ] <- apply(pm$ratio[fi, , drop = FALSE], 2, mean, na.rm = TRUE)
        }
        if (!is.null(pm$purity)) {
          fi <- intersect(fids.pos, rownames(pm$purity))
          if (length(fi)) purity.pos[k, ] <- apply(pm$purity[fi, , drop = FALSE], 2, mean, na.rm = TRUE)
        }
      }
    }

    rda <- S4Vectors::DataFrame(
      isotopologue_form = iso_form,
      isotopologue_id = iso_id,
      compound_id = rep(cid, length(iso_counts)),
      compound_name = rep(cp.name, length(iso_counts)),
      label.isotopologue = label_iso,
      row.names = iso_id
    )

    assays <- list(
      intensity.positive = intensity.pos,
      intensity.negative = intensity.neg,
      ratio.positive = ratio.pos,
      ratio.negative = ratio.neg,
      purity.positive = purity.pos,
      purity.negative = purity.neg
    )
    out[[cid]] <- MSIPIsotopologueData(
      assays = assays,
      rowData = rda,
      colData = S4Vectors::DataFrame(cda)
    )
  }

  return(out)
}


#' @title Populate isotopologue data into MSdev
#' @description
#' Build isotopologue-level \code{MSIPIsotopologueData} objects with
#' \code{\link{get_MSIPIsotopologueData}} and store them in
#' \code{object@advancedAna$MSIP$isotopologue_data}. This function returns the
#' updated \code{MSdev} object.
#'
#' @param object MSdev object.
#' @param ... additional arguments passed to \code{get_MSIPIsotopologueData()}.
#'
#' @return MSdev object with \code{advancedAna$MSIP$isotopologue_data} updated.
#' @export
MSIP_get_isotopologues_data <- function(object, ...) {
  iso.list <- get_MSIPIsotopologueData(object, ...)
  object@advancedAna[["MSIP"]][["isotopologue_data"]] <- iso.list
  object@advancedAna[["MSIP"]][["isotopologues_table"]] <- NULL
  object
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
      iso.table <- object@advancedAna$MSIP$isotopologues_table[[pol]]
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
            ratio_matrix <- object@advancedAna$MSIP$
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
            purity_matrix <- object@advancedAna$MSIP$
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

          object@advancedAna$MSIP$isotopologues_data[[paste0(fid_seed,"_",pol)]] <-
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
    #x <- msdev.M1@advancedAna$MSIP$isotopologues_data[["FT06121_Negative"]]
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


   iso.cfm <- bplapply(object@advancedAna$MSIP$isotopologues_data,
                      ff,
                      BPPARAM = BPPARAM )
  annotated <- sapply(iso.cfm,function(x) "CFM_annotation" %in% names(x@CompoundInfo))
  object@advancedAna$MSIP$isotopologues_data <- iso.cfm[annotated]
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

  acq.list <- object@advancedAna$MSIP$isotopologues_table
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

  isotopologues_list <- object@advancedAna$MSIP$isotopologues_table

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



  iso.data <- object@advancedAna$MSIP$isotopologues_data
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
                                  msipAtomMap = this.cfmd,
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
    MSIPIsotopomerMap <- msip.core@Solve$MSIPIsotopomerMap
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

  iso.data -> object@advancedAna$MSIP$isotopologues_data
  return(object)


}


MSIP_drop_isotopologues_tempdata <- function(object){

  isotopologues_data <- object@advancedAna$MSIP$isotopologues_data
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
  isotopologues_data -> object@advancedAna$MSIP$isotopologues_data
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

  iso.data <- object@advancedAna$MSIP$isotopologues_data
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
    cfmd.ig <- get_MSIPAtomMap_sdf_igraph(cfmd)
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
      comp.eval$sp.consistency.icc[j] <- mean(msip.core@MSIPFragmentMap@FG.data$icc,na.rm = T)
      comp.eval$sp.consistency.cos[j] <- mean(msip.core@MSIPFragmentMap@FG.data$cos,na.rm = T)

      if(!isEmpty(msip.core@Solve$MSIPIsotopomerMap)){
        #length(msip.core@Solve$MSIPIsotopomerMap@)
        comp.eval$FSIS.count[j] <-length(msip.core@Solve$MSIPIsotopomerMap@solve$isotopomer.set)
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

  isotopologues.data.list <- object@advancedAna$MSIP$isotopologues_data

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


  isotopologues.data.list -> object@advancedAna$MSIP$isotopologues_data

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
  iso.list <- object@advancedAna$MSIP$isotopologues_data
  for (i in 1:nrow(process.info.m1)) {

    msip.core <- iso.list[[process.info.m1$feature_id[i]]]$
      MSIP_result[[str_isotope2_num(process.info.m1$iso_count[i])]][[process.info.m1$samples[i]]]
    if (is.null(msip.core)|isEmpty(msip.core)) {
      next
    }
    msip.ifmap <- msip.core@Solve$MSIPIsotopomerMap
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
  iso.list <- object@advancedAna$MSIP$isotopologues_data
  process.info.fg$fragment.count <- NA
  process.info.fg$is.count<- NA
  for (i in 1:nrow(process.info.fg)) {

    msip.core <- iso.list[[process.info.fg$feature_id[i]]]$
      MSIP_result[[str_isotope2_num(process.info.fg$iso_count[i])]][[process.info.fg$samples[i]]]
    if (is.null(msip.core)|isEmpty(msip.core)) {
      next
    }
    msip.map <- msip.core@MSIPFragmentMap
    msip.ifmap<- msip.core@Solve$MSIPIsotopomerMap
    process.info.fg$fragment.count[i]  <- sum(msip.map@fragment.include)
    process.info.fg$is.count[i] <- length( lengths(msip.ifmap@solve$isotopomer.set))

  }

  return(process.info.fg)

}


get_MSIP_result_Statistic <- function(object,
                                      include.merged = F,
                                      show_message = F){

  iso.data <- object@advancedAna$MSIP$isotopologues_data
  iso_ele <- get_MSdev_iso_ele(object)
  target_ele <- get_ele_uniso(iso_ele)
  all.sample <- .get_MSIP_tracer(object)%>%names()
  #traced.sample <- names(na.omit(all.sample))


  comp.eval.list <- list()
  for (i in seq_along(iso.data)) {

    cfmd <- iso.data[[i]][[grep("CFM",names(iso.data[[i]]),value = T)[1]]]
    if (is.null(cfmd)) next
    cfmd.ig <- get_MSIPAtomMap_sdf_igraph(cfmd)
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
        if(!isEmpty(msip.core@Solve$MSIPIsotopomerMap)){
          #length(msip.core@Solve$MSIPIsotopomerMap@)
          comp.eval$FSIS.count[j] <-length(msip.core@Solve$MSIPIsotopomerMap@solve$isotopomer.set)
          comp.eval$FSIS.exist[j]  <- sum(msip.core@Solve$MSIPIsotopomerMap@solve$isotopomer.set.prob>0.00001)
          comp.eval$isotopomer.exist[j] <- sum(msip.core@Solve$MSIPIsotopomerMap@isotopomer.probability>0.00001)
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

  isotopologues_datas <- object@advancedAna$MSIP$isotopologues_data
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
          isotopomers.ele <- msip.core@Solve$MSIPIsotopomerMap@isotopomer.defination
          istp <- rownames(isotopomers.ele)
          isotopomers.prob <- msip.core@Solve$MSIPIsotopomerMap@isotopomer.probability[istp]
          isotopomers.FSIS <- msip.core@Solve$MSIPIsotopomerMap@solve$isotopomer.set%>%
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
          isotopomers.ele <- msip.core@Solve$MSIPIsotopomerMap@isotopomer.defination
          istp <- rownames(isotopomers.ele)
          isotopomers.prob <- msip.core@Solve$MSIPIsotopomerMap@isotopomer.probability[istp]
          isotopomers.FSIS <- msip.core@Solve$MSIPIsotopomerMap@solve$isotopomer.set%>%
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

  object@advancedAna$MSIP$isotopologues_data <- list()
  return(object)
}


#' @title Targeted xcms processing for MSIP
#' @description
#' Perform xcms MS1 processing similar to \code{\link{MSdev_xcmsProcessing}} but
#' with targeted peak detection using \code{xcms::CentWaveParam(roiList = ...)}.
#' The \code{roiList} is constructed from \code{object@advancedAna$MSIP$compound_table}
#' by simulating all possible isotopologues for each compound (default: \code{[13]C})
#' and creating RT/mz ROIs around expected signals.
#'
#' After xcms processing, the function annotates \code{featureDefinitions} with
#' isotopologue metadata (\code{compound_id}, \code{name}, \code{iso_count},
#' \code{iso_form}, \code{iso_seed}) by matching detected features to the
#' theoretical isotopologue grid.
#'
#' If a compound's M+0 (unlabeled) feature is missing from the detected features,
#' a synthetic feature is injected with zero intensity. This ensures downstream
#' functions like \code{\link{get_MSIPIsotopologueData}} can construct a complete
#' isotopologue series for every compound.
#'
#' @param object MSdev object.
#' @param iso_ele character, isotope element, default \code{"[13]C"}.
#' @param mz_ppm numeric, ppm tolerance used to construct ROI mz ranges and
#'   to match detected features to theoretical isotopologues.
#' @param rt_tol numeric, RT tolerance (seconds) used to construct ROI RT ranges
#'   and to match detected features to theoretical isotopologues.
#' @param max_iso integer, optional cap for maximum isotopologue count per compound.
#' @param adjustRT logical, whether to perform retention time adjustment.
#' @param BPPARAM BiocParallel backend passed to \code{xcms::findChromPeaks()}.
#' @param ... passed to \code{xcms::findChromPeaks()}.
#'
#' @return MSdev object with processed \code{object@xcmsData$PositiveMS1} and/or
#' \code{object@xcmsData$NegativeMS1}. The \code{featureDefinitions} of each
#' polarity object are annotated with isotopologue columns
#' (\code{compound_id}, \code{name}, \code{iso_count}, \code{iso_form},
#' \code{iso_seed}). Missing M+0 features are injected as zero-intensity
#' synthetic peaks.
#' @export
MSIP_xcms_processing.targeted <- function(object,
                                         iso_ele = "[13]C",
                                         mz_ppm = 10,
                                         rt_tol = 30,
                                         max_iso = NULL,
                                         adjustRT = F,
                                         BPPARAM = BiocParallel::SnowParam(
                                           workers = 4,
                                           progressbar = TRUE
                                         ),
                                         ...) {

  if (is.null(object@advancedAna$MSIP$compound_table)) {
    stop("compound_table not found in object@advancedAna$MSIP$compound_table. ",
         "Run MSIP_import_compound_table() first.")
  }

  compound_table <- object@advancedAna$MSIP$compound_table
  required_cols <- c("compound_id", "name", "formula", "smiles", "rt")
  missing_cols <- setdiff(required_cols, colnames(compound_table))
  if (length(missing_cols) > 0) {
    stop("compound_table missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Ensure xcmsData exists
  object <- MSdev_get_xcms(object)

  # Filter samples: keep MS1/Both and available msData.files
  sampleInfo <- dplyr::filter(object@sampleInfo, xcmsProcessing %in% c("MS1", "Both")) %>%
    dplyr::filter(!is.na(msData.files))

  if (nrow(sampleInfo) == 0) {
    message_with_time("No MS1/Both samples with msData.files; skip targeted xcms processing.")
    return(object)
  }

  # ---------------------------------------------------------------------------
  # Internal helpers
  # ---------------------------------------------------------------------------
  .as_centwave_with_roi <- function(param_obj, roiList) {
    cwp <- xcms::CentWaveParam()
    if (inherits(param_obj, "CentWaveParam")) {
      # copy common slots if present; avoid failing on version differences
      for (slot_nm in intersect(slotNames(cwp), slotNames(param_obj))) {
        if (slot_nm %in% c("roiList")) next
        try({
          methods::slot(cwp, slot_nm) <- methods::slot(param_obj, slot_nm)
        }, silent = TRUE)
      }
    }
    cwp@roiList <- roiList
    cwp
  }

  .compound_table_to_iso_grid <- function(compound_table,
                                         ion_mode,
                                         iso_ele = "[13]C",
                                         max_iso = NULL) {
    target_ele <- get_ele_uniso(iso_ele)
    adduct <- ifelse(ion_mode == 1, "[M+H]+", "[M-H]-")

    out <- list()
    k <- 0L

    for (j in seq_len(nrow(compound_table))) {
      formula <- compound_table$formula[[j]]
      rt_ref <- suppressWarnings(as.numeric(compound_table$rt[[j]]))
      if (is.na(rt_ref) || !is.finite(rt_ref)) next
      if (is.na(formula) || !nzchar(formula)) next

      cp_adduct <- MSCC::chemform_adduct(formula, adduct, value = "all")
      if (is.null(cp_adduct) || nrow(cp_adduct) == 0) next
      seed_mz <- cp_adduct$chemform.adduct.mz[1]
      if (!is.finite(seed_mz)) next

      max_ele <- get_formula_ele_count(formula, target_ele)
      if (is.na(max_ele) || max_ele <= 0) next
      if (!is.null(max_iso)) max_ele <- min(max_ele, as.integer(max_iso))
      if (max_ele <= 0) next

      ele_counts <- setNames(list(max_ele), target_ele)
      iso_grid <- do.call(MSCC::get_isotope_mass_diff, ele_counts)
      if (is.null(iso_grid) || nrow(iso_grid) == 0) next

      n_iso <- nrow(iso_grid)
      k <- k + 1L
      out[[k]] <- data.frame(
        compound_id = rep(compound_table$compound_id[[j]], n_iso),
        name = rep(compound_table$name[[j]], n_iso),
        iso_form = iso_grid$chemform_diff,
        iso_count = seq(0L, length.out = n_iso),
        mz = seed_mz + iso_grid$mass_diff,
        rt = rep(rt_ref, n_iso),
        stringsAsFactors = FALSE
      )
    }

    if (!length(out)) {
      return(data.frame(compound_id = character(0), name = character(0),
                        iso_form = character(0), iso_count = integer(0),
                        mz = numeric(0), rt = numeric(0),
                        stringsAsFactors = FALSE))
    }
    do.call(rbind, out)
  }

  # Extract mzrt matrix (mz, rt) from iso_grid for roiList construction
  .iso_grid_to_mzrt <- function(iso_grid) {
    if (!nrow(iso_grid)) {
      return(matrix(numeric(), ncol = 2, dimnames = list(NULL, c("mz", "rt"))))
    }
    mzrt <- unique(iso_grid[, c("mz", "rt"), drop = FALSE])
    as.matrix(mzrt)
  }

  # Annotate featureDefinitions with isotopologue info and inject missing M+0
  .annotate_fdf_with_iso <- function(xcms.xcms, iso_grid, ion_mode,
                                     mz_ppm = 10, rt_tol = 30) {
    pol <- ifelse(ion_mode == 0, "Negative", "Positive")
    fdf <- as.data.frame(xcms::featureDefinitions(xcms.xcms))

    if (!nrow(fdf) || !all(c("mzmed", "rtmed") %in% colnames(fdf))) {
      message("[MSIP_xcms_processing.targeted] ", pol,
              ": no featureDefinitions to annotate.")
      return(xcms.xcms)
    }

    # Add annotation columns (initialised to NA)
    if (!("compound_id" %in% colnames(fdf))) fdf$compound_id <- NA_character_
    if (!("name"          %in% colnames(fdf))) fdf$name          <- NA_character_
    if (!("iso_count"     %in% colnames(fdf))) fdf$iso_count     <- NA_integer_
    if (!("iso_form"      %in% colnames(fdf))) fdf$iso_form      <- NA_character_
    if (!("iso_seed"      %in% colnames(fdf))) fdf$iso_seed      <- NA_character_

    # ---- match each row in iso_grid to the closest feature ----
    for (j in seq_len(nrow(iso_grid))) {
      ig <- iso_grid[j, ]
      mz_err <- abs(fdf$mzmed - ig$mz) / ig$mz
      rt_err <- abs(fdf$rtmed - ig$rt)

      hit <- which(mz_err < mz_ppm * 1e-6 & rt_err < rt_tol)
      if (!length(hit)) next

      # pick the closest by m/z error
      best <- hit[which.min(mz_err[hit])]
      fdf$compound_id[best] <- ig$compound_id
      fdf$name[best]          <- ig$name
      fdf$iso_count[best]     <- as.integer(ig$iso_count)
      fdf$iso_form[best]      <- ig$iso_form
    }

    # ---- set iso_seed: for each compound, the feature with iso_count == 0 ----
    compounds <- unique(na.omit(fdf$compound_id))
    for (cid in compounds) {
      sel <- which(fdf$compound_id == cid & !is.na(fdf$iso_count))
      if (!length(sel)) next
      seed_rows <- sel[fdf$iso_count[sel] == 0L]
      seed_fid <- if (length(seed_rows)) rownames(fdf)[seed_rows[[1]]] else NA_character_
      fdf$iso_seed[sel] <- seed_fid
    }

    # ---- inject missing M+0 features ----
    nsamples <- length(MSnbase::fileNames(xcms.xcms))
    pks <- xcms::chromPeaks(xcms.xcms)  # numeric matrix; "sample" col = file index
    n_injected <- 0L

    new_pk_list <- list()
    new_fdf_list <- list()

    for (cid in compounds) {
      # check: does this compound have iso_count == 0?
      sel <- which(fdf$compound_id == cid & !is.na(fdf$iso_count))
      has_m0 <- any(fdf$iso_count[sel] == 0L)
      if (has_m0) next

      # get theoretical M+0 from iso_grid
      m0_rows <- iso_grid[iso_grid$compound_id == cid & iso_grid$iso_count == 0, ]
      if (!nrow(m0_rows)) next
      m0_mz <- m0_rows$mz[[1]]
      m0_rt <- m0_rows$rt[[1]]
      m0_name <- m0_rows$name[[1]]
      m0_form <- m0_rows$iso_form[[1]]

      mz_ppm_half <- m0_mz * mz_ppm / 1e6

      # create one synthetic chromPeak per file (intensity = 0)
      n_existing_pks <- nrow(pks) + length(new_pk_list)
      peak_row_ids <- integer(nsamples)
      for (fi in seq_len(nsamples)) {
        pid <- sprintf("CPM0_%s_%d", cid, fi)
        peak_row_ids[fi] <- n_existing_pks + fi
        new_pk_list[[pid]] <- c(
          mz = m0_mz, mzmin = m0_mz - mz_ppm_half, mzmax = m0_mz + mz_ppm_half,
          rt = m0_rt, rtmin = m0_rt - rt_tol, rtmax = m0_rt + rt_tol,
          into = 0, maxo = 0, sn = NA_real_, sample = fi
        )
      }

      fid_new <- sprintf("CPM0_%s", cid)
      n_injected <- n_injected + 1L

      new_fdf_list[[fid_new]] <- data.frame(
        feature_id = fid_new,
        mzmed = m0_mz, mzmin = m0_mz - mz_ppm_half, mzmax = m0_mz + mz_ppm_half,
        rtmed = m0_rt, rtmin = m0_rt - rt_tol, rtmax = m0_rt + rt_tol,
        npeaks = nsamples,
        peakidx = I(list(peak_row_ids)),
        compound_id = cid, name = m0_name,
        iso_count = 0L, iso_form = m0_form, iso_seed = fid_new,
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }

    if (n_injected > 0L) {
      message("[MSIP_xcms_processing.targeted] ", pol,
              ": injected ", n_injected, " missing M+0 feature(s) with zero intensity.")
    }

    # ---- append synthetic chromPeaks ----
    if (length(new_pk_list)) {
      pk_new_mat <- do.call(rbind, new_pk_list)
      # align columns to match existing chromPeaks matrix exactly
      pk_template <- matrix(NA_real_, nrow = nrow(pk_new_mat), ncol = ncol(pks),
                            dimnames = list(rownames(pk_new_mat), colnames(pks)))
      shared <- intersect(colnames(pk_new_mat), colnames(pks))
      pk_template[, shared] <- pk_new_mat[, shared, drop = FALSE]
      pk_combined <- rbind(pks, pk_template)
      xcms::chromPeaks(xcms.xcms) <- pk_combined

      # update chromPeakData to match new row count
      cpd <- xcms::chromPeakData(xcms.xcms)
      n_new <- nrow(pk_combined) - nrow(cpd)
      if (n_new > 0) {
        cpd_extra <- data.frame(
          msLevel = rep(1L, n_new),
          is_filled = rep(TRUE, n_new),
          stringsAsFactors = FALSE
        )
        # align columns with existing chromPeakData
        for (nm in setdiff(colnames(cpd), colnames(cpd_extra))) {
          cpd_extra[[nm]] <- NA
        }
        cpd_extra <- cpd_extra[, colnames(cpd), drop = FALSE]
        xcms::chromPeakData(xcms.xcms) <- S4Vectors::DataFrame(
          rbind(as.data.frame(cpd), cpd_extra)
        )
      }
    }

    # ---- append new featureDefinitions rows ----
    if (length(new_fdf_list)) {
      fdf_new <- do.call(rbind, new_fdf_list)
      # ensure peakidx exists in original fdf
      if (!"peakidx" %in% colnames(fdf)) {
        fdf$peakidx <- I(rep(list(integer(0)), nrow(fdf)))
      }
      # align columns
      common <- intersect(colnames(fdf), colnames(fdf_new))
      fdf_combined <- rbind(fdf[, common, drop = FALSE], fdf_new[, common, drop = FALSE])
    } else {
      fdf_combined <- fdf
    }

    xcms::featureDefinitions(xcms.xcms) <- S4Vectors::DataFrame(fdf_combined)
    xcms.xcms
  }

  # ---------------------------------------------------------------------------
  # Run targeted MS1 processing for each polarity
  # ---------------------------------------------------------------------------
  polarity.index <- c("0" = "Negative", "1" = "Positive")
  xcms.param <- get_MSdev_param(object)

  for (i in c(0, 1)) {
    sample.info.polarity <- sampleInfo %>%
      dplyr::filter(grepl(i, polarity))
    polarity.tag <- paste0(polarity.index[as.character(i)], "MS1")
    if (!nrow(sample.info.polarity)) next

    xcms.xcms <- filterFile(
      object@xcmsData[[polarity.tag]],
      which(Biobase::pData(object@xcmsData[[polarity.tag]])$sample.name %in%
              sample.info.polarity$sample.name)
    )

    message_with_time("Build roiList from compound_table (", polarity.tag, ") ...")
    iso_grid <- .compound_table_to_iso_grid(
      compound_table = compound_table,
      ion_mode = i,
      iso_ele = iso_ele,
      max_iso = max_iso
    )
    mzrt <- .iso_grid_to_mzrt(iso_grid)
    roiList <- get_xcms_roi_list(
      mzrt = mzrt,
      xcms.xcms = xcms.xcms,
      ppm = mz_ppm,
      rt_tol = rt_tol,
      ion_mode = i
    )
    message_with_time("roiList size: ", length(roiList))

    cwp <- .as_centwave_with_roi(xcms.param$findChromPeaks, roiList = roiList)
    xcms_param_targeted <- list(
      findChromPeaks = cwp,
      groupChromPeaks = xcms.param$groupChromPeaks
    )

    xcms.xcms <- xcmsProcessingMS1(
      xcms.xcms = xcms.xcms,
      ion_mode = i,
      xcms_param = xcms_param_targeted,
      adjustRT = adjustRT,
      BPPARAM = BPPARAM,
      ...
    )

    xcms.xcms <- xcms_get_feature_stat(xcms.xcms)
    message_with_time("Annotate featureDefinitions with isotopologue info (", polarity.tag, ") ...")
    xcms.xcms <- .annotate_fdf_with_iso(xcms.xcms, iso_grid, ion_mode = i,
                                         mz_ppm = mz_ppm, rt_tol = rt_tol)
    object@xcmsData[[polarity.tag]] <- xcms.xcms
  }

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

        this.data <- object@advancedAna$MSIP$isotopologues_data[[this.fid]]


        this.cfmd <- this.data@CompoundInfo$CFM_annotation
        this.msip.core <- this.data@MSIPIsotopologueDatas[[format_isotopologue(this.isotopologue,"M")]][[this.sample]]
        this.sp <- this.msip.core@Spectra
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


  MSIP.mol.igs <- object@advancedAna$MSIP$isotopomer_Molecule_igraph
  suppressWarnings(file.remove(file))
  open_dir(dirname(file))

  for (i in 1:nrow(MSIP.mol.igs)) {

    ### info
    {

      fid <- rownames(MSIP.mol.igs)[i]
      message_with_time(fid)
      this.msip.mtblt <-  object@advancedAna$MSIP$isotopologues_data[[fid]]
      cp <- this.msip.mtblt@CompoundInfo$name

      this.info  <-
        paste0(fid,": ",cp,"\n",
               "Formula: ",object@advancedAna$MSIP$isotopologues_data[[fid]]@CompoundInfo$formula,"\n",
               "RT: ",str_digit(object@advancedAna$MSIP$isotopologues_data[[fid]]@CompoundInfo$rt,0),"\n",
               "mz: ",str_digit(object@advancedAna$MSIP$isotopologues_data[[fid]]@CompoundInfo$mz,4),"\n"
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
  object@advancedAna$MSIP$isotopomer_Molecule_igraph <-MSIP.mol.igs
  return(object)

}


#' @title Get isotopomer data from MSdev object
#' @description Extract and solve isotopomer data for compounds in the compound table.
#' For each compound, calculates m/z of isotopologues, matches MS2 spectra,
#' constructs CFM data and MSIPCoreData, and solves for isotopomer distributions.
#'
#' @param object MSdev object
#' @param mode Mode: "untargeted" (default) or "targeted"
#' @param iso_ele Isotope element, default `"[13]C"`
#' @param ppm PPM tolerance for m/z matching, default 10
#' @param rt.tol RT tolerance in seconds, default 10
#' @param sp_top Optional integer. If provided, keep only the top-N spectra (by
#'   TIC) within each split group before constructing MSIPCoreData. See
#'   `Spectra_filter_TIC()`.
#' @param int_thresh Intensity threshold for MSIPCore_solve, default 10^3.8
#' @param certainty_thresh Certainty threshold for MSIPCore_solve, default 0.6
#' @param weight_fun Weight function for MSIPCore_solve, default .intensity_weight
#' @param iso_count_max Maximum isotopologue count to calculate, default NULL (auto)
#' @param check_temp Check temp directory for cached CFM data, default TRUE
#' @param temp_dir Temporary directory for CFM data cache
#' @param BPPARAM BiocParallel parameter for parallel processing
#'
#' @return MSdev object with isotopomer_data stored in obj@advancedAna$MSIP$isotopomer_data
#' @export
#'
MSIP_get_isotopomer_data <- function(object,
                                      mode = c("untargeted", "targeted")[1],
                                      iso_ele = "[13]C",
                                      ppm = 5,
                                      rt.tol = 30,
                                      sp_top = NULL,
                                      int_thresh = 10^3.8,
                                      certainty_thresh = 0.6,
                                      weight_fun = .intensity_weight,
                                      iso_count_max = NULL,
                                      check_temp = TRUE,
                                      temp_dir = NULL,
                                      BPPARAM = SerialParam(progressbar = TRUE)) {

  # Set default temp_dir if not provided
  if (is.null(temp_dir)) {
    if (is.null(object@projectInfo$CompoundDB_path)) {
      temp_dir <- get_dir_expand_from_onedrive("/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb_cfmd")
    } else {
      temp_dir <- paste0(object@projectInfo$CompoundDB_path, "_cfmd")
    }
  }
  # Save temp_dir into msdev object for reuse
  object@projectInfo$MSIP_cfmd_tempdir <- temp_dir

  # ============================================================================
  # Step 1: Check if compound_table exists
  # ============================================================================
  if (is.null(object@advancedAna$MSIP$compound_table)) {
    stop("compound_table not found in obj@advancedAna$MSIP$compound_table. ",
         "Run MSIP_import_compound_table first.")
  }

  compound_table <- object@advancedAna$MSIP$compound_table

  # ============================================================================
  # Step 2: Check compound_table format
  # ============================================================================
  required_cols <- c("compound_id", "name", "formula", "smiles", "rt")
  missing_cols <- setdiff(required_cols, colnames(compound_table))
  if (length(missing_cols) > 0) {
    stop("compound_table missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ============================================================================
  # Get isotope element info
  # ============================================================================
  target_ele <- get_ele_uniso(iso_ele)
  if (is.null(iso_count_max)) {
    iso_count_max <- 30  # default maximum isotopologue count
  }

  # ============================================================================
  # Get sample sources
  # ============================================================================
  sample_sources <- unique(object@sampleInfo$sample.source)
  if (length(sample_sources) == 0) {
    sample_sources <- "unknown"
  }

  # ============================================================================
  # Initialize result list
  # ============================================================================
  isotopomer_data_list <- list()

  # ============================================================================
  # Load MS2 spectra once outside the loop
  # ============================================================================
  sp.ms2 <- onDiskData_retrieve(object@spectra$MS2_Spectra)
  if (is.null(sp.ms2) || length(sp.ms2) == 0) {
    stop("No MS2 spectra found in object@spectra$MS2_Spectra")
  }

  # ============================================================================
  # Call targeted function if mode is targeted
  # ============================================================================
  if (mode == "targeted") {
    isotopomer_data_list <- get_MSIP_isotopomer_data.targeted(
      sp.ms2 = sp.ms2,
      compound_table = compound_table,
      sample_sources = sample_sources,
      target_ele = target_ele,
      iso_count_max = iso_count_max,
      temp_dir = temp_dir,
      check_temp = check_temp,
      ppm = ppm,
      rt.tol = rt.tol,
      sp_top = sp_top,
      int_thresh = int_thresh,
      certainty_thresh = certainty_thresh,
      weight_fun = weight_fun,
      iso_ele = iso_ele
    )
  }

  # ============================================================================
  # Package result as list of MSIPCoreData matrices
  # ============================================================================
  # ============================================================================
  # Step 5: Store in object@advancedAna$MSIP$isotopomer_data
  # ============================================================================
  object@advancedAna$MSIP$isotopomer_data <- isotopomer_data_list

  message_with_time("MSIP_get_isotopomer_data completed. ",
                   "Processed ", length(isotopomer_data_list), " compounds.")

  return(object)
}


#' @title Targeted MSIP isotopomer data extraction
#' @description Internal function for targeted MSIP isotopomer analysis.
#' Processes each compound in the compound table: calculates m/z of isotopologues,
#' matches MS2 spectra, constructs CFM data and MSIPCoreData, and solves.
#'
#' @param sp.ms2 MS2 spectra object (onDiskData retrieved)
#' @param compound_table Compound table from MSdev object
#' @param sample_sources Unique sample sources
#' @param target_ele Target element (result of get_ele_uniso)
#' @param iso_count_max Maximum isotopologue count
#' @param temp_dir Temp directory for CFM cache
#' @param check_temp Whether to check/use temp cache
#' @param ppm PPM tolerance
#' @param rt.tol RT tolerance
#' @param int_thresh Intensity threshold
#' @param certainty_thresh Certainty threshold
#' @param weight_fun Weight function
#' @param iso_ele Isotope element string
#'
#' @return List of isotopomer data results
#' @keywords internal
#'
get_MSIP_isotopomer_data.targeted <- function(sp.ms2,
                                               compound_table,
                                               sample_sources,
                                               target_ele,
                                               iso_count_max,
                                               temp_dir,
                                               check_temp,
                                               ppm,
                                               rt.tol,
                                               sp_top,
                                               int_thresh,
                                               certainty_thresh,
                                               weight_fun,
                                               iso_ele) {

  # Initialize result list
  isotopomer_data_list <- list()

  .infer_iso_count_max_from_iso_form <- function(iso_form, iso_ele = "[13]C") {
    if (is.null(iso_form) || is.na(iso_form) || !nzchar(iso_form)) return(0)
    if (!grepl(iso_ele, iso_form, fixed = TRUE)) return(0)
    x <- gsub(iso_ele, "", iso_form, fixed = TRUE)
    x <- gsub("[^0-9]", "", x)
    if (!nzchar(x)) return(1)
    as.integer(x)
  }

  # ============================================================================
  # Loop through each compound
  # ============================================================================
  for (i in seq_len(nrow(compound_table))) {
    cp <- compound_table[i, ]
    compound_id <- cp$compound_id
    name <- cp$name
    formula <- cp$formula
    smiles <- cp$smiles
    rt_ref <- cp$rt

    message_with_time("Processing compound: ", compound_id, " (", name, ")")

    # --------------------------------------------------------------------------
    # Step 3a: Calculate mz of isotopologues and adducts (M+H and M-H)
    # --------------------------------------------------------------------------
    # Get base formula mz for both positive and negative adducts
    adducts <- c("[M+H]+", "[M-H]-")
    polarity_vals <- c(1, 0)  # 1 for positive, 0 for negative

    # Calculate base mz for each adduct
    base_mz_list <- list()
    for (idx in seq_along(adducts)) {
      adduct <- adducts[idx]
      cp_adduct <- MSCC::chemform_adduct(formula, adduct, value = "all")
      if (is.null(cp_adduct) || nrow(cp_adduct) == 0) {
        base_mz_list[[adduct]] <- NA
      } else {
        base_mz_list[[adduct]] <- cp_adduct$chemform.adduct.mz[1]
      }
    }

    # Calculate maximum isotopologue count for this compound
    max_iso <- get_formula_ele_count(formula, target_ele)
    if (is.na(max_iso) || max_iso == 0) {
      message_with_time("  Skipping ", compound_id, ": no ", target_ele, " atoms found")
      next
    }

    # Build isotopologue mz grid for each adduct
    isotopologue_mz_list <- list()
    for (idx in seq_along(adducts)) {
      adduct <- adducts[idx]
      seed_mz <- base_mz_list[[adduct]]
      if (is.na(seed_mz)) next

      ele_counts <- setNames(list(max_iso), target_ele)
      iso_grid <- do.call(MSCC::get_isotope_mass_diff, ele_counts)
      iso_grid$iso_mz <- seed_mz + iso_grid$mass_diff
      iso_grid$iso_form <- iso_grid$chemform_diff
      iso_grid$adduct <- adduct
      isotopologue_mz_list[[adduct]] <- iso_grid
    }

    # Combine all isotopologue m/z values
    all_isotopologue_mz <- do.call(rbind, isotopologue_mz_list)
    if (is.null(all_isotopologue_mz) || nrow(all_isotopologue_mz) == 0) {
      message_with_time("  Skipping ", compound_id, ": could not calculate isotopologue m/z")
      next
    }

    # --------------------------------------------------------------------------
    # Step 3b: Check if Spectra exist for the isotopologue by matching mz and rt
    # --------------------------------------------------------------------------
    # Match spectra by mz and rt for each polarity
    matched_spectra <- list()
    for (idx in seq_along(adducts)) {
      adduct <- adducts[idx]
      polarity_val <- polarity_vals[idx]
      iso_grid <- isotopologue_mz_list[[adduct]]
      if (is.null(iso_grid) || nrow(iso_grid) == 0) next

      # Filter spectra by polarity
      sp.pol <- filterPolarity(sp.ms2, polarity_val)
      if (length(sp.pol) == 0) next
      sp.pol.idx <- which(sp.ms2$polarity == polarity_val)

      # Match mz first, then match rt
      match_res <- match_mz_rt(
        mz1 = iso_grid$iso_mz,
        mz2 = sp.pol$isolationWindowTargetMz,
        mz.ppm = ppm
      )
      if (nrow(match_res) == 0) next

      # If rt reference is available, filter by rt
      if (!is.na(rt_ref) && !is.null(rt_ref)) {
        rt2 <- rtime(sp.pol)[match_res$ion2]
        in_rt <- abs(rt2 - rt_ref) <= rt.tol

        if (!any(in_rt)) {
          rt_target <- rt_ref + c(-rt.tol, rt.tol)
          rt_selected <- range(rt2, na.rm = TRUE)
          message_with_time(
            nrow(match_res), " sp targeted to ", compound_id,
            " with rt range: ", format(rt_target[1], digits = 4), " - ", format(rt_target[2], digits = 4),
            ", selected 0 sp with rt range: ", format(rt_selected[1], digits = 4), " - ", format(rt_selected[2], digits = 4)
          )
          next
        }

        if (sum(in_rt) < nrow(match_res)) {
          rt_target <- rt_ref + c(-rt.tol, rt.tol)
          rt_selected <- range(rt2[in_rt], na.rm = TRUE)
          message_with_time(
            nrow(match_res), " sp targeted to ", compound_id,
            " with rt range: ", format(rt_target[1], digits = 4), " - ", format(rt_target[2], digits = 4),
            ", selected ", sum(in_rt), " sp with rt range: ", format(rt_selected[1], digits = 4), " - ", format(rt_selected[2], digits = 4)
          )
        }

        match_res <- match_res[in_rt, , drop = FALSE]
      }

      if (nrow(match_res) == 0) next

      # Add iso_form and iso_mz to match_res from iso_grid
      match_res$iso_form <- iso_grid$iso_form[match_res$ion1]
      match_res$iso_mz <- iso_grid$iso_mz[match_res$ion1]
      match_res$sp_idx <- sp.pol.idx[match_res$ion2]

      # Store matched spectra info
      match_res$adduct <- adduct
      match_res$polarity <- polarity_val
      matched_spectra[[adduct]] <- match_res
    }

    if (length(matched_spectra) == 0) {
      message_with_time("  Skipping ", compound_id, ": no matching spectra found")
      next
    }

    # --------------------------------------------------------------------------
    # Step 3c: Separate spectra by isotopologue and sample.source
    # --------------------------------------------------------------------------
    # Combine matched spectra from all adducts
    all_matched <- do.call(rbind, matched_spectra)

    # Get unique isotopologue forms
    iso_forms <- unique(all_matched$iso_form)

    # Split spectra by isotopologue form and polarity and sample.source
    spectra_by_isotope <- list()
    for (iso_form in iso_forms) {
      # Get spectra indices (in sp.ms2) for this isotopologue across polarities
      sp_indices <- unique(all_matched$sp_idx[all_matched$iso_form == iso_form])
      if (length(sp_indices) == 0) next

      sp.iso <- sp.ms2[sp_indices]

      # Get polarity of these spectra
      sp.polarity <- ifelse(length(sp.iso$polarity) > 0, sp.iso$polarity[1], NA)

      # Split by sample.source
      if (!is.null(sp.iso$sample.source)) {
        sp.split <- split(sp.iso, sp.iso$sample.source)
      } else {
        sp.split <- list("unknown" = sp.iso)
      }

      spectra_by_isotope[[iso_form]] <- list(
        sp = sp.split,
        mz = unique(all_matched$iso_mz[all_matched$iso_form == iso_form]),
        polarity = sp.polarity
      )
    }

    if (length(spectra_by_isotope) == 0) {
      message_with_time("  Skipping ", compound_id, ": no valid isotopologue spectra")
      next
    }

    # --------------------------------------------------------------------------
    # Step 3d: Construct MSIPAtomMap using get_MSIPAtomMap_from_smiles
    # --------------------------------------------------------------------------
    msipAtomMap_list <- list()
    for (idx in seq_along(adducts)) {
      adduct <- adducts[idx]
      polarity_val <- polarity_vals[idx]

      msipAtomMap <- tryCatch({
        get_MSIPAtomMap_from_smiles(
          smiles = smiles,
          compound_id = compound_id,
          ppm = ppm,
          adduct = adduct,
          check_temp = check_temp,
          iso_ele = iso_ele,
          temp_dir = temp_dir
        )
      }, error = function(e) {
        message_with_time("  MSIPAtomMap error for ", compound_id, " ", adduct, ": ", e$message)
        return(NULL)
      })

      if (!is.null(msipAtomMap)) {
        msipAtomMap_list[[adduct]] <- list(
          msipAtomMap = msipAtomMap,
          polarity = polarity_val
        )
      }
    }

    if (length(msipAtomMap_list) == 0) {
      message_with_time("  Skipping ", compound_id, ": could not generate MSIPAtomMap data")
      next
    }

    # --------------------------------------------------------------------------
    # Step 3e & 3f: Construct MSIPCoreData and solve for each isotopologue/sample
    # --------------------------------------------------------------------------
    # Build result matrix: rows = isotopologues ([13]C0, [13]C1, [13]C2...), cols = sample.source
    result_matrix <- matrix(
      list(),
      nrow = length(iso_forms),
      ncol = length(sample_sources),
      dimnames = list(iso_forms, sample_sources)
    )

    # Process each isotopologue and sample combination
    for (iso_form in names(spectra_by_isotope)) {
      # Determine iso_count_max from iso_form & iso_ele (e.g. "[13]C1" -> 1)
      this_iso_count_max <- .infer_iso_count_max_from_iso_form(iso_form, iso_ele = iso_ele)
      if (is.na(this_iso_count_max)) this_iso_count_max <- 0
      this_iso_count_max <- min(this_iso_count_max, max_iso, iso_count_max)

      iso_data <- spectra_by_isotope[[iso_form]]
      sp.by.pol <- iso_data$sp

      # Loop through sample.source
      for (sample_src in names(sp.by.pol)) {
        sp.iso <- sp.by.pol[[sample_src]]
        if (length(sp.iso) == 0) next

        # Get polarity of these spectra
        sp.polarity <- ifelse(length(sp.iso$polarity) > 0, sp.iso$polarity[1], NA)

        # Loop through polarity
        msip.core.pos <- NULL
        msip.core.neg <- NULL

        for (pol_idx in seq_along(polarity_vals)) {
          polarity_val <- polarity_vals[pol_idx]
          adduct <- adducts[pol_idx]
          msipAtomMap_info <- msipAtomMap_list[[adduct]]

          # Get spectra for this polarity
          sp.pol <- sp.iso[sp.iso$polarity == polarity_val]
          if (length(sp.pol) == 0) next
          if (!is.null(sp_top)) {
            sp.pol <- Spectra_filter_TIC(
              sp.pol,
              topN = sp_top,
              split_var = c("sample.source", "collisionEnergy", "polarity")
            )
            if (!length(sp.pol)) next
          }

          tryCatch({
            # Construct MSIPCoreData
            msip.core <- get_MSIPCoreData(
              sp.iso = sp.pol,
              msipAtomMap = msipAtomMap_info$msipAtomMap,
              iso_count_max = this_iso_count_max,
              iso_ele = iso_ele,
              ppm = ppm
            )
            if (polarity_val == 1) {
              msip.core.pos <- msip.core
            } else {
              msip.core.neg <- msip.core
            }

          }, error = function(e) {
            message_with_time("  Error solving ", compound_id, " ", iso_form, " ",
                            adduct, " ", sample_src, ": ", e$message)
          })
        }

        # Merge pos/neg first, then solve
        msip.core.merged <- NULL
        if (!is.null(msip.core.pos) && !is.null(msip.core.neg)) {
          msip.core.merged <- MSIPCore_merge(msip.core.pos, msip.core.neg,
                                             suffix1 = "Positive", suffix2 = "Negative")
        } else if (!is.null(msip.core.pos)) {
          msip.core.merged <- msip.core.pos
        } else if (!is.null(msip.core.neg)) {
          msip.core.merged <- msip.core.neg
        }

        if (!is.null(msip.core.merged) && !isEmpty(msip.core.merged)) {
          msip.core.merged <- MSIPCore_solve(
            msip.core.merged,
            int_thresh = int_thresh,
            certainty_thresh = certainty_thresh,
            weight_fun = weight_fun
          )
          msip.core.merged@Solve$iso_count_max <- this_iso_count_max
          msip.core.merged@Solve$iso_form <- iso_form
          msip.core.merged@Solve$iso_ele <- iso_ele
          msip.core.merged@Solve$ppm <- ppm
          msip.core.merged@Solve$rt_tol <- rt.tol
          message_with_time("  Solved ", compound_id, " ", iso_form, " ", sample_src)
        }

        result_matrix[iso_form, sample_src] <- list(msip.core.merged)
      }
    }

    # --------------------------------------------------------------------------
    # Step 3g: Store compound result
    # --------------------------------------------------------------------------
    # Store in result list (named by compound_id)
    isotopomer_data_list[[compound_id]] <- result_matrix

  }  # End loop through compounds

  return(isotopomer_data_list)
}
