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
    x.iso.cfm <- iso.cfm$FT7478_Positive
    cfmd <- x.iso.cfm$CFM_annotation
    iso.count <- stringr::str_extract(names(x.iso.cfm),"[:digit:]+")%>%
      as.numeric()%>%na.omit()
    atom_ele <-  vdata(cfmd@fragment_igraph$Fragment001)$atom
    names(atom_ele) <- vdata(cfmd@fragment_igraph$Fragment001)$name
    c_ele <- atom_ele[atom_ele=="C"]
    for (i.iso in iso.count) {

      this.sp <-x.iso.cfm[[paste0("M",i.iso)]]
      this.sp <- Spectra_filter_noise(this.sp)
      sp.frag.data <- CFM_annotate_isotopologues(this.sp,
                                 cfmd  = cfmd,ppm = 10,
                                 iso.count = i.iso)%>%
        dplyr::filter(!is.na(iso))
      ### frag group to label fraction
      {
        fg.idx <- split(1:nrow(sp.frag.data),sp.frag.data$fragment_group)
        frag.iso.matrix <- matrix(
          nrow = length(fg.idx),ncol = i.iso+1,
          dimnames = list(names(fg.idx),paste0("M",0:i.iso)))
        frag.int.sum <- c()
        for (i.fg in seq_along(fg.idx)) {
          x.df <- sp.frag.data[fg.idx[[i.fg]],]
          x.int <- x.df%>%
            dplyr::select(-mz,-collisionEnergy)%>%
            tidyr::pivot_wider(names_from ="iso",
                               id_cols = "sp.id",
                               values_from = "intensity",
                               values_fn = sum)%>%
            tibble::column_to_rownames("sp.id")%>%
            dplyr::select(dplyr::starts_with("M"))%>%
            as.matrix()
          to.add <- setdiff(paste0("M",0:i.iso),colnames(x.int))
          x.int <- cbind(matrix(0,nrow(x.int),length(to.add),
                                dimnames = list(NULL,to.add)),x.int)
          x.int <- x.int[,paste0("M",0:i.iso),drop =F]
          x.int[is.na(x.int)] <- 0
          x.weight <- rowSums(x.int)
          x.int <- t(apply(x.int,1,function(z) z/sum(z)))
          x.int.weighted <- apply(x.int,2,weighted.mean,m = x.weight)
          frag.iso.matrix[i.fg,] <- x.int.weighted
          frag.int.sum[i.fg] <- sum(x.weight)
        }
        }


      ### frag group to C atom prob
      {
        frag.c.matrix <- matrix(ncol = length(c_ele),
                         nrow = length(fg.idx),
                         dimnames = list(names(fg.idx),
                                         names(c_ele)))
        for (i.fg in seq_along(fg.idx)) {

          this.frag.group <- names(fg.idx)[i.fg]
          this.frags <- cfmd@fragment_define[cfmd@fragment_define$fragment_group==this.frag.group,]
          this.frag.ratio <-frag.iso.matrix[i.fg,]
          this.frag.atom <- get_cfm_data_fg_atom_map(cfmd,this.frag.group)
          this.frag.c <- this.frag.atom[names(c_ele)]
          #this.iso.expectation <- sum(str_extract_num(names(this.frag.ratio))*this.frag.ratio)
          #this.frag.c <- this.frag.c*this.iso.expectation/sum(this.frag.c)
          #this.frag.c <- this.frag.c[this.frag.c!=0]
          frag.c.matrix[this.frag.group,names(this.frag.c)] <- this.frag.c
        }

      }

      ### vis
      {
       # h1 <- heatmap_atom_iso_prob(frag.c.matrix)
       # h2 <- heatmap_atom_iso_prob(frag.iso.matrix)
       # h1+h2
       hm <-  heatmap.frag.group.maps(frag.c.matrix ,frag.iso.matrix)
        open_plot_win(hm,10,5)
      }

      ### merge duplicate and complementary
      {
        z <- frag.c.matrix
        z[z>0] <- 1
        ### duplicated
        frag.c.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.c.matrix))
        frag.iso.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.iso.matrix))
        frag.int.sum1 <- c()
        z.split <- split(1:nrow(z),apply(z,1,paste0,collapse = ";"))
        for (i.z in seq_along(z.split)) {
          idx <- z.split[[i.z]]
          this.frag.c <- apply(frag.c.matrix[idx,,drop =F],2,
                               mean,weight = frag.int.sum[idx])
          this.frag.iso <- apply(frag.iso.matrix[idx,,drop =F],2,
                               mean,weight = frag.int.sum[idx])
          frag.c.matrix1 <- rbind(frag.c.matrix1,this.frag.c)
          frag.iso.matrix1 <- rbind(frag.iso.matrix1,this.frag.iso)
          frag.int.sum1 <- c(frag.int.sum1, max(frag.int.sum[idx]) )
          rn <- rownames(frag.iso.matrix)[idx][which.max(frag.int.sum[idx])]
          rownames(frag.c.matrix1)[i.z] <- rownames(frag.iso.matrix1)[i.z]  <- rn
        }
        frag.c.matrix <- frag.c.matrix1
        frag.iso.matrix <- frag.iso.matrix1
        frag.int.sum <- frag.int.sum1
        ### complementary
        z <- frag.c.matrix
        z[z>0] <- 1
        z.comple <- apply(z,1,function(x){
          z1 <- t(t(z)+x)
          apply(z1,1,function(xx){all(xx==1)})
          })
        z.comple <- which(z.comple,arr.ind = T)
        z.split <- sapply(1:nrow(z),function(x){
          x.c <- z.comple[z.comple[,1] == x,2]
          sort( c(x,x.c))
        })
        z.split <- z.split[!duplicated(z.split)]
        frag.c.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.c.matrix))
        frag.iso.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.iso.matrix))
        frag.int.sum1 <- c()
        for (i.z in seq_along(z.split)) {
          idx <- z.split[[i.z]]
          this.frag.c <- apply(frag.c.matrix[idx,,drop =F],2,
                               mean,weight = frag.int.sum[idx])
          this.frag.iso <- apply(frag.iso.matrix[idx,,drop =F],2,
                                 mean,weight = frag.int.sum[idx])
          frag.c.matrix1 <- rbind(frag.c.matrix1,this.frag.c)
          frag.iso.matrix1 <- rbind(frag.iso.matrix1,this.frag.iso)
          frag.int.sum1 <- c(frag.int.sum1, max(frag.int.sum[idx]) )
          rn <- rownames(frag.iso.matrix)[idx][which.max(frag.int.sum[idx])]
          rownames(frag.c.matrix1)[i.z] <- rownames(frag.iso.matrix1)[i.z]  <- rn
        }
        frag.c.matrix <- frag.c.matrix1
        frag.iso.matrix <- frag.iso.matrix1
        frag.int.sum <- frag.int.sum1


      }


      ### probability distribution
      {
        save(frag.c.matrix,frag.iso.matrix,frag.int.sum,file = "temp.rda")

      }

    }






  }



}



