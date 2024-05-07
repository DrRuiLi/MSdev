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
get_MSdev_isotopologues <- function(object,iso_ele = "[13]C"){

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
    iso <- grep(pattern = "_seed$",colnames(xcms.fdf),value = T)
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
                                        mz = this.seed.df$mzmed,
                                        rt = this.seed.df$rtmed,
                                        formula = this.seed.df$formula,
                                        smiles = this.seed.df$smiles,
                                        score = this.seed.df$score,
                                        adduct = this.seed.df$adduct,
                                        polarity = this.seed.df$polarity)
      }
      ### Spectra split
      if (length(this.sp)>0) {
        this.sp.list <- split(this.sp,this.sp$feature_id)
        names(this.sp.list) <- paste0("M",
                                      this.fdf$iso_count[ match(names(this.sp.list),this.fdf$feature_id)])

        this.list$Spectra <- this.sp.list
      }

      ### filter
      {
        if (!"M0" %in% names(this.list$Spectra)) next
        if (!all(c(10,20,40) %in% collisionEnergy(this.list$Spectra$M0))) next
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
                                            ppm = 20,
                                            BPPARAM = SnowParam(
                                              workers = 12,
                                              progressbar = T)){

  ff <- function(x){

    ### combine sp of seed to cfm-annotate
    {
      ### process M0 Spectra
      {
        ### M0 just annotated for smiles assign
        seed.sp.c <- Spectra_filter_noise(x$Spectra$M0)
        #seed.sp.c <- normalizeSpectra(sp = seed.sp.c,norm_to = "tic")
        seed.sp.c <- combineSpectra_groupby_ce(seed.sp.c,ppm = ppm,
                                    minProp = 0.3,
                                    plot = F)
        seed.sp.c <- Spectra_fill_3CE(seed.sp.c)
      }
      CFM_result <- NA
      try.return <- try(
        CFM_result <- CFM_annotate(smiles_or_inchi = x$compound_info$smiles,
                                     spectrum_file = seed.sp.c,
                                     ppm_mass_tol = ppm,
                                     abs_mass_tol = 0.005,
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
                                                ppm = 5,
                                                ...)
    xcms.xcms <- xcms_get_feature_isotope_label(xcms.xcms,
                                                isotope = isotope,
                                                ...)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
  }

  return(object)

}


MSdev_get_iso_acq_list <- function(object){



  acq.list <- list()
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
      xcms.fdf.iso$ms1_purity[xcms.fdf.iso$is_seed] <- 1
      xcms.fdf.iso <- xcms.fdf.iso%>%
        ### filter not assigned iso
        dplyr::filter(!is.na(C13_seed))%>%
        dplyr::group_by(C13_seed)%>%
        ### filter not annotated
        dplyr::filter(any( !is.na( compound_id ) ),
                      any( is_labeled ))%>%
        dplyr::mutate(
          compound_id = na.omit(compound_id),
          name = na.omit(name),
          adduct = na.omit(adduct),
          formula = na.omit(formula),
          smiles = na.omit(smiles),
          selected_to_acq = T)%>%
        dplyr::mutate(
          selected_to_acq = case_when(
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
                      "is_labeled","compound_id","adduct",
                      "name","formula","smiles","mean.iso",
                      "mean.uniso", "ms1_purity","selected_to_acq")
    }

   # ### assign and split
   # {
   #   scan.df <- simulate_prm(xcms.fdf.stat)
   #   xcms.fdf.assign <- xcms.fdf.stat %>%
   #     dplyr::mutate(scan.count = table(scan.df$ion_id)[as.character(1:n())],
   #                   idx = sample(1:(ceiling(min.scan/min(scan.count))),n(),replace =T))
   #   xcms.fdf.assign <- split(xcms.fdf.assign,xcms.fdf.assign$idx)
#
   # }

    ### trans to QE
    #pol.list <- lapply(xcms.fdf.assign,QE_list_2feature_def,keep=T)

    #edit_df_in_excel(xcms.fdf.stat)
    acq.list[[pol]] <- xcms.fdf.iso


  }
  #unlist(acq.list,recursive = F)

  object@statData$iso.acq.list <- acq.list
  return(object)

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



get_isotopologues_label_fraction <- function(iso.cfm){

  .f <- function(x.iso.cfm){
    #x.iso.cfm <- iso.cfm$FT7478_Positive
    cfmd <- x.iso.cfm$CFM_annotation
    iso.count <- stringr::str_extract(names(x.iso.cfm$Spectra),"[:digit:]+")%>%
      as.numeric()%>%na.omit()
    atom_ele <-  vdata(cfmd@fragment_igraph[[1]])$atom
    names(atom_ele) <- vdata(cfmd@fragment_igraph[[1]])$name
    c_ele <- atom_ele[atom_ele=="C"]
    for (i.iso in iso.count) {

     # message(i.iso)
      if (i.iso==0) {
        next
      }
      this.sp <-x.iso.cfm$Spectra[[paste0("M",i.iso)]]
      this.sp <- Spectra_filter_noise(this.sp)
      sp.frag.data <- CFM_annotate_isotopologues(this.sp,
                                 cfmd  = cfmd,ppm = 10,
                                 iso.count = i.iso)%>%
        dplyr::filter(!is.na(iso))
      if (!nrow(sp.frag.data)) next
      fg.map <- get_frag_group_map(sp.frag.data,cfmd,c_ele,i.iso)
      fg.map <- merge_frag_group_map(fg.map)
      iso.form.map <- get_iso_form_map(fg.map ,atom_prob = F)
      iso.form.map <- get_iso_form_prob_GLPK(iso.form.map)
      c.prob <- get_iso_from_C_prob(iso.form.map, x.iso.cfm$CFM_annotation,i.iso)
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
        x.iso.cfm$msip_data[[paste0("M",i.iso)]] <- list(
          fg.map = fg.map,
          iso.form.map = iso.form.map,
          c.prob = c.prob
        )


      }

    }

    return(x.iso.cfm)






  }

  bplapply(iso.cfm,.f,BPPARAM = SerialParam(progressbar = T))



}






