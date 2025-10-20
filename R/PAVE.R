#' @describeIn PVAE initialize from MSdev
#' @title PAVE analysis
#' @param MSdev
#'
#' @returns PAVE
#' @export
#'
#' @examples
#' get_PAVE_from_MSdev(msdev)
get_PAVE_from_MSdev <- function(object){

  object@sampleInfo <-
    object@sampleInfo%>%
    dplyr::mutate(
      isotope_tracer = paste0(ifelse(grepl("13C",sample.name),"[13]C","")),
      isotope_tracer = paste0(isotope_tracer,ifelse(grepl("15N",sample.name),"[15]N","")),
      .after = sample.name
    )
  message_with_time("Please check the column ",
                    crayon::red("isotope_tracer"),
                    " to confirm tracer")

  object <- MSdev_checkSampleInfo(object)
  object <- MSdev_update_xcms_pdata(object)
  return(object)


}


PAVE_get_atom_count <- function(object){


  polarity.index <- c("0" = "Negative",
                      "1"="Positive")
  for (i.pol in 0:1) {


    ### get xcms
    {
      polarity.tag <- paste0(polarity.index[as.character(i.pol)],"MS1")
      xcms.xcms <- object@xcmsData[[polarity.tag]]
      if (is.null(xcms.xcms)) next
      #xcms.se <- get_xcms_feature_se(xcms.xcms)
    }

    ### find +C +N from peaks in unlabeled sample
    {

      xcms.xcms <- PAVE_xcms_find_CN(xcms.xcms)
      xcms.xcms -> object@xcmsData[[polarity.tag]]

    }


  }

  return(object)


}

PAVE_junk_remover <- function(object){

  object <- MSdev_get_Stat(object,metabolite = F)


  ### isotope
  {
    iso.diff <- data.frame(
      chemform = c("[13]C1C-1","[13]C1C-1","[18]O1O-1","[18]O1O-1","[15]N1N-1","[34]S1S-1","[37]ClCl-1"),
      charge = c(1,2,1,2,1,1,1)
    )%>%
      dplyr::mutate(mz = MSCC::chemform_mz(chemform)/charge)

    adduct.diff

    polarity.index <- c("0" = "Negative",
                        "1"="Positive")
    for (i.pol in 0:1) {

      polarity.tag <- paste0(polarity.index[as.character(i.pol)],"MS1")
      xcms.xcms <- object@xcmsData[[polarity.tag]]



    }


  }




}

PAVE_xcms_find_CN <- function(xcms.xcms,rt.tol = 20, ppm= 10 ){



  ### prepare data
  {

    CN_mass_diff_df <- get_CN_mass_diff_table(C_max = 99,N_max = 10)

    xcms.fdf <- featureDefinitions(xcms.xcms)%>%
      as.data.frame()%>%
      dplyr::mutate(pave_seed = "",
                    pave_CN = "",
                    pave_cor = NA)
    xcms.pda <- pData(xcms.xcms)
    xcms.pave.sample <- xcms.pda%>%
      dplyr::filter(sample.type %in% c("S12C14N","S13C14N","S12C15N","S13C15N"))
    xcms.val <- featureValues(xcms.xcms, missing  = 0,value = "maxo")
    #pave.sample.val <- apply(xcms.val,1,median_f, f= xcms.pda$sample.type)%>%t
    #pave.sample.val <- pave.sample.val[,c("S12C14N","S13C14N","S12C15N","S13C15N")]
    #xcms.fdf[,c("S12C14N","S13C14N","S12C15N","S13C15N")] <- pave.sample.val
  }


  ### peaks count in sample
  if(F) {
    sample.12C14N <- xcms.pda%>%
      dplyr::filter(sample.type == "S12C14N")
    xcms.fdf$S12C14N <- apply(xcms.val[,sample.12C14N$sampleNames],1,
                         function(x)sum(!is.na(x)))%>%unname()

    sample.13C14N <- xcms.pda%>%
      dplyr::filter(sample.type == "S13C14N")
    xcms.fdf$S13C14N <- apply(xcms.val[,sample.13C14N$sampleNames],1,
                                function(x)sum(!is.na(x)))

    sample.12C15N <- xcms.pda%>%
      dplyr::filter(sample.type == "S12C15N")
    xcms.fdf$S12C15N <- apply(xcms.val[,sample.12C15N$sampleNames],1,
                                function(x)sum(!is.na(x)))
    sample.13C15N <- xcms.pda%>%
      dplyr::filter(sample.type == "S13C15N")
    xcms.fdf$S13C15N <- apply(xcms.val[,sample.13C15N$sampleNames],1,
                                function(x)sum(!is.na(x)))
   feature.S12C14N <- xcms.fdf %>%dplyr::filter(S12C14N >= nrow(sample.12C14N))
   feature.S13C14N <- xcms.fdf %>%dplyr::filter(S13C14N >= nrow(sample.13C14N))
   feature.s12C15N <- xcms.fdf %>%dplyr::filter(S12C15N >= nrow(sample.12C15N))
   feature.S13C15N <- xcms.fdf %>%dplyr::filter(S13C15N >= nrow(sample.13C15N))
  }

  cn.list <- BiocParallel::bplapply(
    seq_along(xcms.fdf$feature_id),
    FUN = function(i.ft,xcms.fdf,rt.tol,ppm,xcms.val,xcms.pave.sample){

      this.fid <- xcms.fdf$feature_id[i.ft]
      #message_with_time(this.fid)
      this.mz <- xcms.fdf$mzmed[i.ft]
      this.rt <- xcms.fdf$rtmed[i.ft]
      this.CN_mass_diff_df <- CN_mass_diff_df%>%
        dplyr::filter( C_count <= this.mz/14)

      mz.pred <- this.mz+this.CN_mass_diff_df$mass_diff
      this.fdf <- xcms.fdf%>%
        dplyr::filter(abs(rtmed - this.rt) < rt.tol,
                      mzmed >= this.mz
        )%>%
        dplyr::mutate(idx = match_mz(mzmed,mz.pred,mz.ppm = ppm),
                      this.CN_mass_diff_df[idx,]  )%>%
        dplyr::filter(!is.na(idx)   )%>%
        dplyr::mutate(
          pave_seed = paste0(pave_seed,this.fid,";"),
          pave_CN = paste0("C",C_count,"N",N_count)
        )
      #this.fdf[,c("S12C14N","S13C14N","S12C15N","S13C15N")] <-
      #  this.fdf[,c("S12C14N","S13C14N","S12C15N","S13C15N")]/
      #  this.fdf[1,c("S12C14N")]


      possible.c.count <- unique(this.fdf$C_count)%>%setdiff(0)
      possible.n.count <- unique(this.fdf$N_count)
      cn.comb <- expand.grid(C = possible.c.count,
                             N = possible.n.count,
                             p.cor = NA)

      ### score possible C and N pattern
      for (i.cn in 1:nrow(cn.comb)) {
        this.c <- cn.comb$C[i.cn]
        this.n <- cn.comb$N[i.cn]
        all.form <- c("C0N0",paste0("C0N",this.n),paste0("C",this.c,"N0"),paste0("C",this.c,"N",this.n))
        if (all(all.form %in% this.fdf$pave_CN) ) {

          cn.ft <-this.fdf%>%
            dplyr::filter(pave_CN %in% all.form)%>%
            dplyr::group_by(pave_CN)%>%
            dplyr::slice_max(peakMaxo)%>%
            dplyr::ungroup()
          m.detected <- xcms.val[cn.ft$feature_id,xcms.pave.sample$sampleNames]
          colnames(m.detected) <- xcms.pave.sample$sample.type
          rownames(m.detected) <- cn.ft$pave_CN
          m.detected <- m.detected/m.detected[1,1]
          m.ideal <- MSdev:::get_ideal_CN_ratio(this.c,this.n)%>%t
          m.ideal <- m.ideal[rownames(m.detected),colnames(m.detected)]

          p.cor <- cor(as.vector(m.detected),as.vector(m.ideal))
          cn.comb$p.cor[i.cn] <- p.cor

        }
      }


      ### save to xcms.fdf
      if (any(cn.comb$p.cor > 0.75,na.rm =T)) {
        cn.comb <- cn.comb%>%
          dplyr::slice_max(p.cor)
        all.form <- c("C0N0",paste0("C0N",cn.comb$N),paste0("C",cn.comb$C,"N0"),paste0("C",cn.comb$C,"N",cn.comb$N))
        message_with_time("Pattern: ",all.form[4],"; Cor = ",cn.comb$p.cor)
        this.fdf <- this.fdf %>%
          dplyr::filter(pave_CN %in% all.form)
        this.fdf[this.fid,]$pave_cor <- max(cn.comb$p.cor,na.rm  =T)
        return(this.fdf)
        #xcms.fdf[this.fdf$feature_id,]$pave_seed <-this.fdf$pave_seed
        #xcms.fdf[this.fdf$feature_id,]$pave_CN <-this.fdf$pave_CN
      }

#


    },
    xcms.fdf = xcms.fdf,rt.tol =rt.tol,ppm = ppm,xcms.val=xcms.val,xcms.pave.sample=xcms.pave.sample,
    BPPARAM = SnowParam(workers = 6,progressbar = T)
  )
  cn.list <- cn.list[!sapply(cn.list,is.null)]
  for (i.cnl in seq_along(cn.list)) {
    this.fdf <- cn.list[[i.cnl]]
    xcms.fdf[this.fdf$feature_id,]$pave_seed <- this.fdf$pave_seed
    xcms.fdf[this.fdf$feature_id,]$pave_CN <- this.fdf$pave_CN
    xcms.fdf[this.fdf$feature_id,]$pave_cor <- this.fdf$pave_cor
  }

  featureDefinitions(xcms.xcms)$pave_seed <- xcms.fdf$pave_seed
  featureDefinitions(xcms.xcms)$pave_CN <- xcms.fdf$pave_CN
  featureDefinitions(xcms.xcms)$pave_cor <- xcms.fdf$pave_cor

  return(xcms.xcms)
}

PAVE_xcms_junk_remover <- function(xcms.xcms){

}

get_CN_mass_diff_table <- function(C_max=100,N_max=50){



  ### mass and max count define
  {
    C13_mass_diff= MSCC::chemform_mz("[13]CC-1")
    N15_mass_diff= MSCC::chemform_mz("[15]NN-1")


  }

  ### mass diff matrix
  if(F){
    C_mass_diff_matrix <-
      matrix(
        rep(C13_mass_diff * (0:C_max), N_max+1),
        nrow = C_max+1
      )
    N_mass_diff_matrix <-
      matrix(
        rep(N15_mass_diff * (0:N_max), C_max+1),
        ncol = N_max+1,byrow = T
      )
    CN_mass_diff_matrix <-
      C_mass_diff_matrix+N_mass_diff_matrix
    rownames(CN_mass_diff_matrix) <- paste0("C",num2str(0:C_max))
    colnames(CN_mass_diff_matrix) <- paste0("N",num2str(0:N_max))

  }

  ### mass diff data.frame
  {

    CN_mass_diff_df <-
      expand.grid(
        C_count = 0:C_max,
        N_count = 0:N_max
      )%>%
      dplyr::mutate(
        mass_diff = C_count * C13_mass_diff + N_count * N15_mass_diff
      )

  }

  return(CN_mass_diff_df)



}

get_ideal_CN_ratio <- function(C = 10 , N = 2){

  m <- diag(rep(1,4))
  colnames(m) <- c("C0N0",paste0("C",C,"N0"),paste0("C0N",N),paste0("C",C,"N",N))
  rownames(m) <-c("S12C14N","S13C14N","S12C15N","S13C15N")
 # m <- m[,!duplicated(colnames(m))]
  return(m)
}

PAVE_report <- function(object,mzr = c(0,Inf)){

  ### CN labeled peaks
  {

    pos.fdf <- featureDefinitions(object@xcmsData$PositiveMS1)%>%
      as.data.frame()%>%
      dplyr::filter(between.range(mzmed,r = mzr))

    pos.cn.m <- sum(pos.fdf$pave_cor>0,na.rm = T)
    pos.cn.pks <- sum(pos.fdf$pave_seed!="",na.rm = T)
    pos.pks <-nrow(pos.fdf)

    neg.fdf <- featureDefinitions(object@xcmsData$NegativeMS1)%>%
      as.data.frame()%>%
      dplyr::filter(between.range(mzmed,r = mzr))
    neg.cn.m <- sum(neg.fdf$pave_cor>0,na.rm = T)
    neg.cn.pks <- sum(neg.fdf$pave_seed!="",na.rm = T)
    neg.pks <-nrow(neg.fdf)

    df <- data.frame(
      Positive = c( pos.pks,
                    paste0(pos.cn.pks," (",str_digit(pos.cn.pks/pos.pks*100,2),"%)"),
                    paste0(pos.cn.m," (",str_digit(pos.cn.m/pos.pks*100,2),"%)")),
      Negative = c( neg.pks,
                    paste0(neg.cn.pks," (",str_digit(neg.cn.pks/neg.pks*100,2),"%)"),
                    paste0(neg.cn.m," (",str_digit(neg.cn.m/neg.pks*100,2),"%)"))
    )
    rownames(df) <- c("Total peaks",
                      "CN labeled peaks",
                      "CN labeled metabolites")
    print(df)
    edit_df_in_excel(df)

  }

  ### SNR
  {

    xcms.xcms <- object@xcmsData$PositiveMS1
    xcms.xcms@.processHistory[[1]]@type <- "Peak detection"
    xcms.xcms@.processHistory[[1]]@param <- object@processingInfo$MSdevParam$findChromPeaks
    p1 <- plot_xcms_peaks_SN_distribution(xcms.xcms)

    xcms.xcms <- object@xcmsData$NegativeMS1
    xcms.xcms@.processHistory[[1]]@type <- "Peak detection"
    xcms.xcms@.processHistory[[1]]@param <- object@processingInfo$MSdevParam$findChromPeaks
    p2 <- plot_xcms_peaks_SN_distribution(xcms.xcms)
    open_plot_win(p1+p2,10,5)
  }

}

xcms_get_feature_isotopologues_multi_tracer <- function(xcms.xcms,
                                           iso_ele = c("[13]C","[15]N"),
                                           max_label = c(100,10),
                                           ppm = 10,
                                           rt.tol = 5,
                                           net.degree.ratio = 0.3){



  fdf.iso.connect <- get_xcms_feature_iso_connection(xcms.xcms,
                                                     iso_ele,max_label,
                                                     ppm,rt.tol )


  ### assign isotopologues
  {

    xcms.fdf <- featureDefinitions(xcms.xcms)%>%
      as.data.frame()
    xcms.fdf[,paste0("iso_seed")] <- NA
    xcms.fdf[,paste0("iso_count")] <- NA
    fdf.iso.igraph <- igraph::graph_from_data_frame(fdf.iso.connect)
    #fdf.iso.igraph <- igraph_filter_vertex(fdf.iso.igraph,degree(fdf.iso.igraph)>2)
    node.group <- igraph::components(fdf.iso.igraph)$membership
    xcms.fdf <- as.data.frame(xcms.fdf)
    rownames(xcms.fdf) <-xcms.fdf$feature_id
    message( length(unique(na.omit(node.group)))," iso-group"  )
    message( (length(node.group))," iso-features"  )

    for (i in seq_along(unique(node.group))) {

      #message(i)
      this.nodes <- names(which(node.group==i))
      this.fdf <- xcms.fdf[this.nodes,]
      this.iso <- fdf.iso.connect %>%
        dplyr::filter(from%in%this.nodes | to %in% this.nodes)

      this.igraph <- igraph::graph_from_data_frame(this.iso,vertices =this.fdf[,1:7] )
      #visNetwork::visIgraph(this.igraph)
      this.iso.assign <- get_iso_net_assign(this.igraph,net.degree.ratio = net.degree.ratio)
      xcms.fdf[names(this.iso.assign$iso.seed),
               "iso_seed"] <- this.iso.assign$iso.seed
      xcms.fdf[names(this.iso.assign$iso_count),
               "iso_count"] <- this.iso.assign$iso_count
      xcms.fdf[this.nodes,"iso_connection_group"] <- i
      this.fdf <- xcms.fdf[this.nodes,]


    }





  }


  ### save to featuredef
  {

    xcms.fdf.temp <- featureDefinitions(xcms.xcms)
    rownames(xcms.fdf) <- xcms.fdf$feature_id
    xcms.fdf.temp[,"iso_seed"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_seed"]
    xcms.fdf.temp[,"iso_count"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_count"]
    xcms.fdf.temp[,"iso_connection_group"] <- xcms.fdf[xcms.fdf.temp$feature_id,"iso_connection_group"]
    xcms.fdf.temp -> featureDefinitions(xcms.xcms)
    message("Get ",
            sum(!is.na(xcms.fdf.temp[,"iso_count"])),
            " isotopologues")

  }

  return(xcms.xcms)

}
