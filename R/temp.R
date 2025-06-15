Spectra_annotate_cfmd <- function(
    sp,
    cfmd,
    iso_ele = "[13]C",
    iso_count_max = 3 ,
    ppm = 10
){


  ### pre-process
  {


    if (!"fragment_group"%in% colnames(cfmd@peak_assignment)) {
      cfmd <- cfm_data_get_fragment_group(cfmd)
    }

    cfm.peaks.data <- cfmd@peak_assignment%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::mutate(collisionEnergy = case_when(energy == "energy0"~10,
                                                energy == "energy1"~20,
                                                energy == "energy2"~40,
      ))

    if (!nrow(cfm.peaks.data))
      cfm.peaks.data <-  cfmd@peak_assignment%>%
      dplyr::mutate(mz = 0)
  }


  ### iso mz
  {

    diff.formula <- paste0(iso_ele,get_ele_uniso(iso_ele),"-1")
    iso.mz.diff <- (0:iso_count_max)*MSCC::chemform_mz(diff.formula)
    mz.labeled.m <- matrixSub(cfm.peaks.data$mz,-iso.mz.diff)%>%
      `colnames<-`(paste0("M",0:iso_count_max))%>%
      as.data.frame()%>%
      dplyr::mutate(fragment_group=cfm.peaks.data$fragment_group)%>%
      dplyr::distinct(M0,.keep_all = T)%>%
      tidyr::pivot_longer(paste0("M",0:iso_count_max),
                          names_to = "iso_count",
                          values_to = "mz")%>%
      dplyr::mutate(iso_count = str_extract_num(iso_count))


  }


  ### annotate
  {



    idx <- lapply(mz(sp),function(x){
      match_mz(x,mz.labeled.m$mz, mz.ppm = ppm) })

    sp$iso_count <- lapply(idx,function(x)mz.labeled.m$iso_count[x])#%>%as("IntegerList")
    sp$fragment_group <- lapply(idx,function(x)mz.labeled.m$fragment_group[x])#%>%as("CharacterList")
    sp$fragment_group_mz <- lapply(idx,function(x)mz.labeled.m$mz[x])#%>%as.array("NumericList")

  }




  return(sp)
}


Spectra_calculate_fragment_iso_ratio <- function(sp){

  sp$fragment_group_int_sum <- lapply(1:length(sp),function(i){

    fg <- sp$fragment_group[[i]]
    isoc <- sp$iso_count[[i]]
    mzi <- peaksData(sp)[[i]]
    fg.i <- rep(NA,length(fg))

    for (i.fg in unique(na.omit (fg))) {
      idx <- which(fg == i.fg)

      fg.i[idx] <- sum(mzi[idx,"intensity"])
    }

    return(fg.i)

  })


  sp$fragment_group_ratio <- lapply(1:length(sp),function(i){

    fg <- sp$fragment_group[[i]]
    isoc <- sp$iso_count[[i]]
    mzi <- peaksData(sp)[[i]]
    fg.r <- rep(NA,length(fg))

    for (i.fg in unique(na.omit (fg))) {
      idx <- which(fg == i.fg)

      fg.r[idx] <- mzi[idx,"intensity"]/sp$fragment_group_int_sum[[i]][idx]
    }

    return(fg.r)

  })


  return(sp)


  }



get_Spectra_fg_ratio_se <- function(sp,iso_count_max = 3){



  ### ratio matrix
  {
    fragment_group_ratio_matrix <- lapply(1:length(sp),function(i){

      fg <- sp$fragment_group[[i]]
      isoc <- sp$iso_count[[i]]
      fg.r <- sp$fragment_group_ratio[[i]]
      fg.f <-  paste0(fg,"_M",isoc)

      names(fg.r) <- fg.f
      fg.r <- na.omit(fg.r)
      fg.r <-sapply(split(fg.r,names(fg.r )),sum)

      return(fg.r)


    })

    fragment_group_ratio_matrix <- do.call(bind_rows,fragment_group_ratio_matrix)%>%as.matrix()

    rownames(fragment_group_ratio_matrix) <- sp$sp_id


  }


  ### int matrix
  {

    fragment_group_int_sum_matrix <- lapply(1:length(sp),function(i){

      fg <- sp$fragment_group[[i]]
      isoc <- sp$iso_count[[i]]
      fg.r <- sp$fragment_group_int_sum[[i]]
      fg.f <-  paste0(fg,"_M",isoc)

      names(fg.r) <- fg.f
      fg.r <- na.omit(fg.r)
      fg.r <-sapply(split(fg.r,names(fg.r )),sum)

      return(fg.r)


    })

    fragment_group_int_sum_matrix <- do.call(bind_rows,fragment_group_int_sum_matrix)%>%as.matrix()

    rownames(fragment_group_int_sum_matrix) <- sp$sp_id


  }


  ### cda matrix
  {


    fg <- unique(na.omit(unlist(sp$fragment_group)))
    cda <- expand.grid(fragment_group = fg, iso_count = 0:iso_count_max,stringsAsFactors  = F)%>%
      dplyr::mutate(FG_isotopologue =  paste0(fragment_group,"_M",iso_count))%>%
      dplyr::arrange(FG_isotopologue)

    fragment_group_ratio_matrix <- get_matrix_value_fill_with_NA(
      fragment_group_ratio_matrix,
      colnames_vec = cda$FG_isotopologue
    )

    fragment_group_int_sum_matrix <- get_matrix_value_fill_with_NA(
      fragment_group_int_sum_matrix,
      colnames_vec = cda$FG_isotopologue
    )


  }


  ### store in se
  {
    sp.se <- SummarizedExperiment::SummarizedExperiment(
      fragment_group_ratio_matrix,colData = cda
    )

    assay(sp.se,2) <- fragment_group_int_sum_matrix

  }


  return(sp.se)


}



get_Spectra_fg_ratio_se_merge <- function(sp.se){




  fgs <- unique(colData(sp.se)$fragment_group)
  rda <- data.frame(
    fragment_group = fgs,
    int_sum = NA,
    icc = NA,
    cos = NA,
    peaks_count = 0,

    row.names = fgs
  )

  fg_ratio_list <- list()
  for (i.fg in fgs) {

    i.se <- sp.se[,sp.se$fragment_group==i.fg]
    i.rm <- assay(i.se,1)
    i.rm[is.na(i.rm)] <- 0

    i.intsum <- assay(i.se,2)%>%
      apply(1,function(x) mean(x , na.rm = T))


    ### exclude empty sp
    {
      idx <- which(!is.na(i.intsum))
      i.rm <- i.rm[idx,,drop = F]
      i.intsum <- i.intsum[idx]

      }

    ### weighted mearge
    {
      i.rm.weighted <- apply(i.rm,2,weighted.mean,w = i.intsum)
      i.intsum.weighted <- weighted.mean(i.intsum,i.intsum)
      }


    ### consistency
    {
      if(length(i.intsum)>1){

        i.icc <- irr::icc(t(i.rm), model = "twoway",
                          type = "consistency", unit = "single")$value
        i.cos <- lsa::cosine(i.rm.weighted,t(i.rm))
        i.cos.weight <- weighted.mean(i.cos,w = log10(i.intsum))

      }else{

        i.icc <- NA
        i.cos.weight <- 1
      }

    }

    ### data integration
    {
      names(i.rm.weighted) <- paste0("M",i.se$iso_count)
      fg_ratio_list[[i.fg]] <- i.rm.weighted
      rda[i.fg,"int_sum"] <- i.intsum.weighted
      rda[i.fg,"icc"] <- i.icc
      rda[i.fg,"cos"] <- i.cos.weight
      rda[i.fg,"peaks_count"] <- length(i.intsum)
    }


  }


  ### return se
  {

    fg.ratio.matrix <-do.call(bind_rows,fg_ratio_list)%>%
      as.matrix()
    rownames(fg.ratio.matrix) <- fgs
    fg.se <- SummarizedExperiment::SummarizedExperiment(
      fg.ratio.matrix,rowData =rda
    )

    return(fg.se)
  }


}


get_MSIPFragmentMap <- function(sp,
                                cfmd,
                                iso_ele = "[13]C",
                                iso_count_max = 3,
                                ppm= 10){



  fg.map <- new("MSIPFragmentMap")


  ### FG - Ratio
  {

    if (!"fragment_group" %in% spectraVariables(sp.iso)) {
      sp <- Spectra_annotate_cfmd(sp = sp,cfmd = cfmd,iso_ele ,iso_count_max,ppm )
      sp <- Spectra_calculate_fragment_iso_ratio(sp)
    }


    fg.ratio.se <- get_Spectra_fg_ratio_se(sp,iso_count_max = iso_count_max)
    fg.se <- get_Spectra_fg_ratio_se_merge(fg.ratio.se)

  }


  ### FG - atom
  {

    target_atoms <- get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfmd),get_ele_uniso(iso_ele))
    frag.atom.matrix <- matrix(ncol = length(target_atoms),
                               nrow = nrow(fg.se),
                               dimnames = list(rownames(fg.se),
                                               target_atoms))
    for (i.fg in rownames(fg.se)) {

      this.frag.group <- i.fg

      this.frag.c <-cfmd@fragment_group_map[this.frag.group,]
      frag.atom.matrix[this.frag.group,names(this.frag.c)] <- this.frag.c
    }

  }

  ### FG data
  {
    fg.data <- rowData(fg.se)%>%
      as.data.frame()
    fg.data$include <- T

    fg.data <- fg.data[rownames(frag.atom.matrix),]
  }



  fg.map@FG.atom.matrix <- frag.atom.matrix[fg.data$fragment_group,]
  fg.map@FG.ratio.matrix <- assay(fg.se)[fg.data$fragment_group,]
  fg.map@FG.data <- fg.data

  return(fg.map)






}
