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
                                                ppm = 10,
                                                ...)
    xcms.xcms <- xcms_get_feature_isotope_label(xcms.xcms,
                                                isotope = isotope,
                                                ...)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
  }

  return(object)

}


MSIP_get_isotopologues_table <- function(object,extract_chrom = F){



  object@statData$MSIP <- list()
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]


    ### MS1 purity
    {
      sample.idx <- !is.na(pData(xcms.xcms)$isotope_label)
      #sample.idx <-  pData(xcms.xcms)$group=="FSLowCFullGlu"
      if ( !"ms1_purity"%in%   colnames(xcms.xcms@msFeatureData$featureDefinitions)) {
        xcms.xcms <- xcms_get_feature_purity(xcms.xcms,
                                             grouped = T,
                                             selected.sample = sample.idx)
        xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
      }
      xcms.fdf <- featureDefinitions(xcms.xcms)%>%
        as.data.frame()
      }
    ### Comp info
    {
      cpdb <- CompDb(object@projectInfo$CompoundDB_path)
      dbinfo <- get_CompDb_info(cpdb,
                                xcms.fdf$compound_id,
                                keys = c("name","formula","smiles"))
      xcms.fdf <- cbind(xcms.fdf,dbinfo[,c("name","formula","smiles")])
    }
    ### calc intensity of iso and un-iso labeled sample
    {
      xcms.pdata <- pData(xcms.xcms)%>%
        dplyr::filter(!sample.type%in% c("Blank"))
      xcms.fv <- featureValues(xcms.xcms,value = "maxo",missing = 1)[,xcms.pdata$sampleNames]
      uniso.mean <- xcms.fv[,xcms.pdata$sampleNames[is.na(xcms.pdata$isotope_label)]]%>%
        apply(1,mean)
      iso.mean <- xcms.fv[,xcms.pdata$sampleNames[!is.na(xcms.pdata$isotope_label)]]%>%
        apply(1,mean)

      xcms.fdf$mean.iso <- log10(iso.mean)
      xcms.fdf$mean.uniso <- log10(uniso.mean)


    }
    ### iso stat and filter
    {

      xcms.fdf.iso <- xcms.fdf%>%
        dplyr::mutate(is_seed = feature_id%in% C13_seed )
      xcms.fdf.iso[!xcms.fdf.iso$is_seed,
                   c("compound_id","name","adduct","score",
                     "mz_ref","rt_ref","formula","smiles")] <- NA
      #xcms.fdf.iso$ms1_purity[xcms.fdf.iso$is_seed] <- 1
      na.unique <- function(x){
        if (all(is.na(x))) return(NA)
        return(unique(na.omit(x)))
      }
      xcms.fdf.iso <- xcms.fdf.iso%>%
        ### filter not assigned iso
        #dplyr::filter(!is.na(C13_seed))%>%
        dplyr::group_by(C13_seed)%>%
        dplyr::mutate(is_isotopologues = !is.na(C13_seed)&any(!is.na(compound_id))  )%>%
        ### filter not annotated
        #dplyr::filter(any( !is.na( compound_id ) ),
        #              any( is_labeled ))%>%
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
            mean.iso < 4&mean.uniso<4 ~F,
            ms1_purity < 0.8~ F,
            !is_labeled&!is_seed~ F,
            T~selected_to_acq
          ))%>%
        dplyr::ungroup()%>%
        dplyr::arrange(#-mean.uniso,
          C13_seed,
          C13_count)%>%
        dplyr::select("feature_id","mzmed","rtmed","rtmin","rtmax","peakMaxo","polarity",
                      "score","C13_seed","C13_count",
                      "is_seed","is_isotopologues","is_labeled","compound_id","adduct",
                      "name","formula","smiles","mean.iso", "mean.uniso",
                      grep(pattern = "Ratio_to_seed",colnames(xcms.fdf),value = T),
                      "ms1_purity","selected_to_acq")
    }



    object@statData$MSIP$isotopologues_table[[pol]] <- xcms.fdf.iso


    ### extract chrom
    {
      if (extract_chrom) {
        fid <- xcms.fdf.iso%>%
          dplyr::filter(is_isotopologues)%>%
          dplyr::pull(feature_id)
        xcms.chrom <- featureChromatograms(xcms.xcms,
                                           features = fid,
                                           expandRt = Inf,
                                           filled = T,
                                           BPPARAM = SnowParam(
                                             workers  = min(snowWorkers(), length(fileNames(xcms.xcms))),
                                             progressbar = T)
        )
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


#' get_MSdev_isotopologues
#'
#' extract iso-labeled compound info and spectra,
#' filter compound
#'
#' @param object MSdev
#'
#' @return a list of isotopologues
#' @export
#'
MSIP_get_isotopologues_data <- function(object,iso_ele = "[13]C"){

  iso.list <- list()
  sp.ms2 <- object@spectra$MS2_Spectra
  sp.ms2$sample.source <- object@sampleInfo$sample.source[
    match(sampleNames(sp.ms2),basename(object@sampleInfo$msData.files))]
  #cpdb <- CompoundDb::CompDb(object@projectInfo$CompoundDB_path)

  ### Requirement in xcms features
  ### C13_seed
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.fdf <- object@statData$MSIP$isotopologues_table[[pol]]
    xcms.fdf[,"iso_seed"] <- xcms.fdf[,paste0(trans_iso_ele(iso_ele),"_seed")]
    xcms.fdf[,"iso_count"] <- xcms.fdf[,paste0(trans_iso_ele(iso_ele),"_count")]
    xcms.fdf <- xcms.fdf%>%
      dplyr::group_by(iso_seed)%>%
      dplyr::filter(!is.na(iso_seed),
                    n() > 1,
                    rep(any(is_labeled) ,n()))%>%
      dplyr::ungroup()
    seed.id <- unique(xcms.fdf$iso_seed)
    for (i_seed_id in seq_along(seed.id)) {
      #message(i_seed_id)
      this.fdf <- xcms.fdf%>%
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
                                        polarity = this.seed.df$polarity)
      }
      ### Spectra split
      {
        this.sp <-  lapply(this.fdf$ms2_id,function(x) sp.ms2[x])
        #this.sp <- lapply(this.sp,combineSpectra_groupby_ce)
        names(this.sp) <- paste0("M", this.fdf$iso_count[ match(names(this.sp),
                                                                this.fdf$feature_id)])
        this.sp <- lapply(this.sp , function(x) split(x,x$sample.source))
        this.list$Spectra <- this.sp

      }

      ### Ratio matrix
      {
        ratio_matrix <- this.fdf %>%
          dplyr::select(contains("Ratio_to_seed"))%>%
          as.matrix()
        rownames(ratio_matrix) <- paste0("M",this.fdf$iso_count)
        this.list$compound_info$ratio_matrix <- ratio_matrix
      }


      ### filter
      {
        if ( length(this.list$Spectra$M0)==0 ) next
       # if (!all(c(10,20,40) %in% collisionEnergy(this.list$Spectra$M0))) next
        if (is.na(this.list$compound_info$formula)) next
        uniso.ele <- get_ele_uniso(iso_ele)
        atom.count <- MSCC:::chemform_parse(this.list$compound_info$formula)
        if (!uniso.ele%in% colnames(atom.count)) {
          max.atom <- 0
        }else
          max.atom <- atom.count[,uniso.ele]
        idx <- str_extract_num(names(this.list$Spectra)) <= max.atom
        this.list$Spectra <- this.list$Spectra[idx]

        if (length(this.list$Spectra)<=1) next

        object@statData$MSIP$isotopologues_data[[paste0(seed.id[i_seed_id],"_",pol)]]$compound_info <- this.list$compound_info
        object@statData$MSIP$isotopologues_data[[paste0(seed.id[i_seed_id],"_",pol)]]$Spectra <- this.list$Spectra
      }

      #iso.list[[paste0(seed.id[i_seed_id],"_",pol)]] <- this.list
    }


  }

  #object@statData$MSIP$isotopologues_data <- iso.list
  return(object)

}


#' get_isotopologues_CFM_annotation
#'
#' @param object msdev
#'
#' @return  list of cfm data
#' @import magrittr tidyverse
#' @export
#'
MSIP_get_isotopologues_CFM_annotation <- function(object,
                                                  ppm = 20,
                                                  BPPARAM = SnowParam(progressbar = T)){


  ff <- function(x){
      x <- msdev.combine@statData$MSIP$isotopologues_data[["FT7476_Negative"]]
    ### combine sp of seed to cfm-annotate
    {
      ### process M0 Spectra
      {
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
      CFM_result <- NA
      try.return <- try(
        CFM_result <- CFM_annotate_by_predict(smiles_or_inchi = x$compound_info$smiles,
                                   spectrum_file = seed.sp.c,
                                   ppm_mass_tol = ppm,
                                   abs_mass_tol = 0.005,
                                   param_adduct = switch(as.character(x$compound_info$polarity),
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




  iso.cfm <- bplapply(object@statData$MSIP$isotopologues_data,
                      ff,
                      BPPARAM = BPPARAM )
  annotated <- sapply(iso.cfm,function(x) "CFM_annotation" %in% names(x))
  object@statData$MSIP$isotopologues_data <- iso.cfm[annotated]
  object
}

get_MSdev_iso_acq_list <- function(object){

  acq.list <- object@statData$iso.acq.list
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    acq.list.pol <- acq.list[[pol]]%>%
      dplyr::filter(selected_to_acq)%>%
      dplyr::mutate(rtmed = case_when(
        feature_id == C13_seed ~ rtmed,
        T ~ NA
      ))%>%
      dplyr::group_by(C13_seed)%>%
      dplyr::filter(n()>1)%>%
      dplyr::mutate(rtmed = na.omit(rtmed),
                    rtmin = rtmed - 10,
                    rtmax = rtmed+10
                    )

    acq.list.pol <- acq.list[[pol]] <- QE_list_2feature_def(acq.list.pol)

  }
  return(acq.list)
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
  iso.count<- iso.seed
  for (i.cl in unique(cl)) {
    this.igraph.sub <- igraph_filter_vertex(iso.ig,cl==i.cl)

    igraph::distances(this.igraph.sub,mode = "out",
                      weights = edge.attributes(this.igraph.sub)$closest.iso.count )
    to.delete <- degree(this.igraph.sub)<(length(V(this.igraph.sub))-1)*2*net.degree.ratio
    this.igraph.sub <-delete.vertices(this.igraph.sub,to.delete )
    #visNetwork::visIgraph(this.igraph.sub)
    this.dis <- igraph::distances(this.igraph.sub,mode = "out",
                                  weights = edge.attributes(this.igraph.sub)$closest.iso.count )
    this.dis[this.dis<0] <- 0
    dis.sum <- apply(this.dis,1,sum)
    seed.fid <- names(which.max(dis.sum))
    dis.to.seed <- this.dis[names(which.max(dis.sum)),]
    iso.seed[names(dis.to.seed)] <- seed.fid
    iso.count[names(dis.to.seed)] <- dis.to.seed
  }

  x <- list(iso.seed=iso.seed,iso.count=iso.count)
  return(x)

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



MSIP_get_isotopologues_label_fraction <- function(object,ppm = 20){

  .f <- function(x.iso.cfm){
    #x.iso.cfm <- msdev.combine@statData$MSIP$isotopologues_data$FT0288_Negative
    cfmd <- x.iso.cfm$CFM_annotation
    all.iso.count <- stringr::str_extract(names(x.iso.cfm$Spectra),"[:digit:]+")%>%
      as.numeric()%>%na.omit()
    msip_result <-list()
    natural.ratio.matrix <- get_iso_natural_ratio(
      formula = x.iso.cfm$compound_info$formula,
      iso_ele = "[13]C",
      ratio_matrix = x.iso.cfm$compound_info$ratio_matrix)
    msip_result <- list()
    for (i.iso in all.iso.count) {

     # message(i.iso)
      if (i.iso==0) {
        next
      }
      this.sp <-x.iso.cfm$Spectra[[paste0("M",i.iso)]]
      lengths(this.sp)
      for (i.sample in names(this.sp)) {
        x <- .get_isotopologues_label_fraction(
          sp.iso = this.sp[[i.sample]],
          cfmd = cfmd,
          ppm = ppm,
          iso.count = i.iso,
          natural.ratio = natural.ratio.matrix[paste0("M",i.iso),
                                               paste0("Ratio_to_seed_",i.sample)]
        )
        msip_result[[paste0("M",i.iso)]][[i.sample]] <- x
      }



    }
    x.iso.cfm$MSIP_result <-msip_result

    return(x.iso.cfm)

  }

  object@statData$MSIP$MSIP_result <- bplapply(object@statData$MSIP$isotopologues_data,
           .f,
           BPPARAM = SerialParam(progressbar = T))
  #x.iso.cfm <- object@statData$MSIP$isotopologues_data[[1]]


}


.get_isotopologues_label_fraction <- function(sp.iso,
                                              cfmd,
                                              ppm = 10,
                                              iso.count,
                                              natural.ratio){

  sp.iso <- Spectra_filter_noise(sp.iso)
  #sp.iso <- combineSpectra_groupby_ce(sp.back,
  #                                    minProp = 0.5,plot = T,
  #                                    ppm = 20)
  sp.frag.data <- CFM_annotate_isotopologues(sp.iso,
                                             cfmd  = cfmd,
                                             ppm = ppm,
                                             iso.count = iso.count)
  sp.frag.data <- CFM_spectra_data_int_weight(sp.frag.data,iso.count)
  if (sum(sp.frag.data$sp.id=="combined_sp")==0) return(NA)
  fg.map <- get_frag_group_map(sp.frag.data,cfmd,iso.count = iso.count)
  heatmap.fg.map(fg.map)
  if.map <- get_iso_form_map(fg.map)
  sp.frag.data <- CFM_spectra_data_remove_natural(sp.frag.data,natural.ratio,if.map)
  if (sum(sp.frag.data$sp.id=="combined_sp")==0) return(NA)
  fg.map <- get_frag_group_map(sp.frag.data,cfmd,iso.count = iso.count)
  fg.map <- merge_frag_group_map(fg.map)
  if.map <- get_iso_form_set_map(if.map ,fg.map)
  if.map <- get_iso_form_prob_GLPK(if.map)
  heatmap.ifs.map(if.map)
  c.prob <- get_iso_from_C_prob(if.map, cfmd,iso.count)
  #sum(iso.form.map$iso.form.prob)




  ### vis
  # {
  #   this.dir <- paste0("d:/temp/nad_M",i.iso)
  #   dir.create(this.dir,showWarnings = F)
  #   hm <-  heatmap.fg.map(fg.map)
  #   export::graph2png(hm,
  #                     file = paste0(this.dir,"/Frag.map.png"),
  #                     width = 10,
  #                     height = nrow(fg.map$frag.c.matrix)*0.8)
  #   hm <- heatmap.ifs.map(iso.form.map)
  #   export::graph2png(draw(hm),
  #                     file = paste0(this.dir,"/Iso.form.map.png"),
  #                     width = 10,
  #                     height = nrow(fg.map$frag.c.matrix)*1.5)
  #   p <- vis_sdf_ig_prob(cfmd@fragment_igraph[[1]],c.prob,show.label = T)
  #   saveWidget(p,file = paste0(this.dir,"/Atom.prob.html"))
  # }

  ### SAVE
  {
    return(list(
      sp.data = sp.frag.data,
      fg.map = fg.map,
      if.map = if.map,
      c.prob = c.prob
    ))


  }


}




get_iso_natural_ratio <- function(formula, iso_ele, ratio_matrix  ){

  iso_pattern <- chemform_isotopes_pattern_enviPat(
    formula,
    thresh = min(ratio_matrix))%>%
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


