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


PAVE_get_atom_count <- function(object,BPPARAM = SnowParam(workers = 6,progressbar = T)){


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

      cn.list <- PAVE_find_xcms_CN(xcms.xcms,BPPARAM = BPPARAM)
      object@statData$PAVE[[polarity.index[as.character(i.pol)]]] <- cn.list

    }


  }

  return(object)


}


PAVE_junk_remover <- function(object,ppm = 10,rt.tol = 20){



  ### inside polarity
  {
    cn.seed.pol <- list()
    for (i.pol in 0:1) {

      pol <- ifelse(i.pol==0,"Negative","Positive")

      cn.list <- object@statData$PAVE[[pol]]
      cn.seed <- lapply(cn.list,function(x){
        x %>%dplyr::mutate(pave_formula = paste0("C",max(C_count),
                                                 "N",max(N_count)))%>%
          dplyr::filter(feature_id == pave_seed)
      })%>%data.table::rbindlist()%>%
        as.data.frame()

      cn.seed <- cn.seed%>%
        dplyr::filter(pave_cor > 0.75)%>%
        dplyr::mutate(rtg = cluster_rt(rt = rtmed,rt.tol = 20),
                      pave_junkremover = "")


      ### adduct and isotope
      {

        message_with_time("Find isotope and adduct in ",pol)

        iso.diff <- data.frame(
          chemform = c("[13]C1C-1","[13]C1C-1","[18]O1O-1","[18]O1O-1","[15]N1N-1","[34]S1S-1","[37]ClCl-1"),
          charge = c(1,2,1,2,1,1,1)
        )%>%
          dplyr::mutate(
            type = "isotope",
            mass_diff = MSCC::chemform_mz(chemform)/charge )%>%
          dplyr::filter(charge%in% c(1,2))

        data("pave_adduct")
        adducts.diff <- pave_adduct%>%
          dplyr::filter(polarity == pol)%>%
          dplyr::mutate(type = "adduct")

        mass_diff <- bind_rows(iso.diff,adducts.diff)




        cn.seed.adduct.isotope <- lapply(unique(cn.seed$rtg),
                                         function(x){

                                           #message(x)
                                           this.seed <- cn.seed %>%
                                             dplyr::filter(rtg == x)

                                           if (nrow(this.seed) <2) {
                                             return(this.seed)
                                           }

                                           this.seed.diff <-
                                             expand.grid(ft1 = seq_along(this.seed$feature_id),
                                                         ft2 = seq_along(this.seed$feature_id))%>%
                                             dplyr::filter(ft2>ft1)%>%
                                             dplyr::mutate(mz1 = this.seed$mzmed[ft1],
                                                           mz2 = this.seed$mzmed[ft2],
                                                           mzd = abs(mz1-mz2),
                                                           idx = match_mz( mzd, mass_diff$mass_diff , mz.ppm = 10000 ),
                                                           mass_diff[idx,],
                                                           mzd.ppm = abs(mzd-mass_diff)/mean(mz1+mz2)*1e6
                                             )%>%
                                             dplyr::filter(!is.na(idx),mzd.ppm < ppm)

                                           if (nrow(this.seed.diff)) {

                                             this.seed$pave_junkremover[this.seed.diff$ft2] <-
                                               this.seed.diff$type

                                           }

                                           return(this.seed)

                                         })%>%data.table::rbindlist()%>%
          as.data.frame()



        }

      ### dimer
      {

        message_with_time("Find dimer in ",pol)
        cn.seed.dimer <- lapply(unique(cn.seed$rtg),
                                function(x){

                                  #message(x)
                                  this.seed <- cn.seed %>%
                                    dplyr::filter(rtg == x)

                                  if (nrow(this.seed) <2) {
                                    return(this.seed)
                                  }

                                  this.seed.dimer <-
                                    this.seed%>%
                                    dplyr::mutate(dimer.mz = (mzmed + (i.pol-0.5)*2*1.00727) / 2 )

                                  this.seed.is.dimer <- this.seed%>%
                                    dplyr::mutate(
                                      dimer.matched = match_mz(mzmed , this.seed.dimer$dimer.mz,mz.ppm = ppm ),
                                      pave_junkremover = case_when(!is.na(dimer.matched)~ "dimer",
                                                                   T~""))

                                  this.seed$pave_junkremover <-this.seed.is.dimer$pave_junkremover


                                  return(this.seed)

                                })%>%data.table::rbindlist()%>%
          as.data.frame()



      }

      ### ring
      {


        message_with_time("Find ringing in ",pol)
        xcms.xcms <- object@xcmsData[[paste0(pol,'MS1')]]
        xcms.fdf <- featureDefinitions(xcms.xcms)%>%
          as.data.frame()
        xcms.pda <- pData(xcms.xcms)
        xcms.pave.sample <- xcms.pda%>%
          dplyr::filter(sample.type %in% c("S12C14N","S13C14N","S12C15N","S13C15N"))
        xcms.val <- featureValues(xcms.xcms, missing  = 0,value = "maxo")
        pave.sample.val <- apply(xcms.val,1,mean_f, f= xcms.pda$sample.type)%>%t
        pave.sample.val <- pave.sample.val[,c("S12C14N")]
        xcms.fdf$peakMaxo <- pave.sample.val

        cn.seed.ring <- cn.seed
        for (i in 1:nrow(cn.seed.ring)) {

          this.mz <- cn.seed.ring$mzmed[i]
          this.rt <- cn.seed.ring$rtmed[i]
          this.peakmaxo <- cn.seed.ring$peakMaxo[i]
          this.mz.range50 <- mz.range.ppm(this.mz,50)
          this.mz.range500 <- mz.range.ppm(this.mz,500)
          this.ring.maxo <- xcms.fdf %>%
            dplyr::filter(between.range(rtmed , this.rt+ c(-rt.tol,rt.tol)),
                          between.range(mzmed,c(this.mz.range500[1],this.mz.range50[1])) |
                            between.range(mzmed,c(this.mz.range50[2],this.mz.range500[2])) )%>%
            dplyr::pull(peakMaxo)%>%
            max()

          if(this.ring.maxo>this.peakmaxo*100)
            cn.seed.ring$pave_junkremover[i] <- "ringing"

        }

      }


      ### integrate
      {

        table(cn.seed.adduct.isotope$pave_junkremover)
        table(cn.seed.dimer$pave_junkremover)
        table(cn.seed.ring$pave_junkremover)

        rownames(cn.seed.adduct.isotope) <- cn.seed.adduct.isotope$feature_id
        rownames(cn.seed.dimer) <- cn.seed.dimer$feature_id
        rownames(cn.seed.ring) <- cn.seed.ring$feature_id

        cn.seed.adduct.isotope <- cn.seed.adduct.isotope[cn.seed$feature_id,]
        cn.seed.dimer <- cn.seed.dimer[cn.seed$feature_id,]
        cn.seed.ring <- cn.seed.ring[cn.seed$feature_id,]


        cn.seed$pave_junkremover <-
          integrate_anntation(cn.seed$pave_junkremover,
                              cn.seed.adduct.isotope$pave_junkremover)

        cn.seed$pave_junkremover <-
          integrate_anntation(cn.seed$pave_junkremover,
                              cn.seed.dimer$pave_junkremover)

        cn.seed$pave_junkremover <-
          integrate_anntation(cn.seed$pave_junkremover,
                              cn.seed.ring$pave_junkremover)
        cn.seed.pol[[pol]] <- cn.seed


      }




    }
  }

  ### between polarity
  {

    diff_to_neg <- pave_adduct %>%
      dplyr::filter(polarity=="Positive")%>%
      dplyr::mutate(mass_diff = mass_diff + 1.0078250320*2 - 2* 0.00054857990943 )

    diff_to_pos <- pave_adduct %>%
      dplyr::filter(polarity=="Negative")%>%
      dplyr::mutate(mass_diff = mass_diff - 1.0078250320*2 + 2* 0.00054857990943 )


    cn.seed <- data.table::rbindlist(cn.seed.pol)
    cn.seed <- cn.seed%>%
      dplyr::mutate(rtg = cluster_rt(rt = rtmed,rt.tol = 20))

    cn.seed.list <- list()
    for (i in unique(cn.seed$rtg)) {


      this.cn.seed <- cn.seed%>%
        dplyr::filter(rtg == i)

      possible.adduct.neg <- this.cn.seed %>%
        dplyr::filter(#pave_junkremover=="",
          polarity == 1)%>%
        dplyr::pull(mzmed)%>%
        expand.grid(mz = ., mzd = diff_to_pos$mass_diff)%>%
        dplyr::mutate(mz.expected = mz + mzd)


      possible.adduct.pos <- this.cn.seed %>%
        dplyr::filter(#pave_junkremover=="",
          polarity == 0)%>%
        dplyr::pull(mzmed)%>%
        expand.grid(mz = ., mzd = diff_to_neg$mass_diff)%>%
        dplyr::mutate(mz.expected = mz + mzd)

      this.cn.seed <- this.cn.seed%>%
        dplyr::ungroup()%>%
        dplyr::mutate(
          adduct.match = case_when(
            polarity == 0 ~ match_mz(mzmed,possible.adduct.neg$mz.expected,mz.ppm = ppm),
            polarity == 1 ~ match_mz(mzmed,possible.adduct.pos$mz.expected,mz.ppm = ppm)
          ),
          pave_junkremover = case_when(
            is.na(adduct.match) ~ pave_junkremover,
            T ~ integrate_anntation(pave_junkremover,"opposite_adduct")
          )
        )#%>%
      #dplyr::select(-adduct.match)
      #message(sum(!is.na(this.cn.seed$adduct.match)))

      cn.seed.list[[i]] <- this.cn.seed

    }

    cn.seed <- data.table::rbindlist(cn.seed.list)%>%
      as.data.frame()

  }


  ### Low C
  {
    data("PAVE_LowC_cutoff")
    cn.seed <- cn.seed%>%
      dplyr::mutate(
        pave_lowC_cutoff = PAVE_LowC_cutoff$mass_max[
          match(get_formula_ele_count(pave_formula,"C"),PAVE_LowC_cutoff$c.count) ],
        pave_junkremover =
          case_when(
            mzmed > pave_lowC_cutoff~integrate_anntation(pave_junkremover,"LowC"),
            T~pave_junkremover) )



  }




  ### return
  {

    for (i.pol in 0:1) {

      pol <- ifelse(i.pol==0,"Negative","Positive")
      cn.list <- object@statData$PAVE[[pol]]
      cn.seed.pol <- cn.seed%>%
        dplyr::filter(polarity == i.pol)%>%
        dplyr::mutate(tmp = feature_id)%>%
        tibble::column_to_rownames("tmp")

      cn.list.junkremoved <- lapply(cn.list,function(x){

        x %>%
          dplyr::mutate(pave_formula =  cn.seed.pol[pave_seed,"pave_formula"],
                        pave_junkremover = cn.seed.pol[pave_seed,"pave_junkremover"])
      })



      cn.list.junkremoved -> object@statData$PAVE[[pol]]

    }

    return(object)

  }







}

PAVE_formula_assign <- function(object,ppm = 10,rt.tol = 20){


  #object@statData$PAVE
  cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite"
  cpdb <- CompoundDb::CompDb(cpdb_path)


  cn.list.pol <- list()
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    if (is.null(xcms.xcms)) next
    message_with_time("Find MS1 candidate...")
    xcms.xcms <- xcms_get_feature_ms1_candidate(xcms.xcms,
                                                cpdb,
                                                ppm = ppm)
    xcms.fdf <- featureDefinitions(xcms.xcms)

    cn.list <- object@statData$PAVE[[pol]]
    cn.seed <- lapply(cn.list,function(x){
      x %>%dplyr::mutate(pave_formula = paste0("C",max(C_count),
                                               "N",max(N_count)))%>%
        dplyr::filter(feature_id == pave_seed)
    })%>%data.table::rbindlist()%>%
      as.data.frame()%>%
      dplyr::mutate(pave_formula_matched = F)%>%
      dplyr::filter(pave_cor > 0.75)

    for (i.seed in 1:nrow(cn.seed)) {

      this.fid <- cn.seed$feature_id[i.seed]
      this.pave.formula <- cn.seed$pave_formula[i.seed]

      this.adduct.candidate <- xcms.fdf[this.fid,]$candidate.adduct%>%unlist()
      this.chemform.candidate <- xcms.fdf[this.fid,]$candidate.formula%>%unlist()

      candidate.c.count <- get_formula_ele_count(this.chemform.candidate,"C")
      candidate.n.count <- get_formula_ele_count(this.chemform.candidate,"N")
      candidate.cn.formula <- paste0("C",candidate.c.count,"N",candidate.n.count)

      idx <- match(this.pave.formula,candidate.cn.formula)
      if(!is.na(idx))
        cn.seed$pave_formula_matched[i.seed] <- T

    }


    cn.list.pol[[pol]] <- cn.seed


  }


  ### return
  {

    for (i.pol in 0:1) {

      pol <- ifelse(i.pol==0,"Negative","Positive")
      cn.list <- object@statData$PAVE[[pol]]
      cn.seed.pol <- cn.list.pol[[pol]]%>%
        dplyr::filter(polarity == i.pol)%>%
        dplyr::mutate(tmp = feature_id)%>%
        tibble::column_to_rownames("tmp")

      cn.list.formula.assigned <- lapply(cn.list,function(x){

        x %>%
          dplyr::mutate(pave_formula_matched =  cn.seed.pol[pave_seed,"pave_formula_matched"])
      })



      cn.list.formula.assigned -> object@statData$PAVE[[pol]]

    }

    return(object)

  }

}



PAVE_report <- function(object,file = tempfile(fileext = "pdf"),mzr = c(0,Inf)){


  cn.stat.list <- list()
  for (i.pol in 0:1) {

    pol <- ifelse(i.pol==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.fdf <- featureDefinitions(xcms.xcms)%>%
      as.data.frame()%>%
      dplyr::filter(between.range(mzmed,mzr))

    cn.list <- object@statData$PAVE[[pol]]
    cn.peaks <- cn.list%>%
      data.table::rbindlist()%>%
      as.data.frame()%>%
      dplyr::filter(between.range(mzmed,mzr))%>%
      dplyr::group_by(pave_seed)%>%
      dplyr::mutate(pave_cor = na.omit(pave_cor))%>%
      dplyr::ungroup()%>%
      dplyr::group_by(feature_id)%>%
      dplyr::slice_max(pave_cor)%>%
      dplyr::ungroup()%>%
      dplyr::distinct(feature_id,.keep_all = T)

    cn.peaks.high.cor <- cn.peaks%>%
      dplyr::filter(pave_cor >= 0.75)%>%
      dplyr::mutate(
        isotope = grepl("^isotope",pave_junkremover),
        adduct = grepl("^adduct",pave_junkremover),
        LowC = grepl("^LowC",pave_junkremover),
        opposite_adduct  = grepl("^opposite_adduct",pave_junkremover),
        dimer   = grepl("^dimer",pave_junkremover),
        ringing    = grepl("^ringing",pave_junkremover)
      )

    cn.peaks.high.cor.formula <-cn.peaks.high.cor %>%
      dplyr::filter(pave_junkremover == "")




    cn.stat.list[[pol]]$ATOMCOUNT <-
      list(total_peaks = nrow(xcms.fdf),
           peaks_in_blak = sum(xcms.fdf$Blank>0),
           peaks_withou_labeling = length(setdiff(xcms.fdf$feature_id,cn.peaks$feature_id)),
           peaks_low_cor = sum(cn.peaks$pave_cor < 0.75),
           peaks_high_cor = sum(cn.peaks$pave_cor >= 0.75) )


    cn.stat.list[[pol]]$JUNKREMOVER <-
      list( isotopes = sum(cn.peaks.high.cor$isotope),
            adduct = sum(cn.peaks.high.cor$adduct),
            LowC = sum(cn.peaks.high.cor$LowC),
            opposite_adduct = sum(cn.peaks.high.cor$opposite_adduct),
            dimer = sum(cn.peaks.high.cor$dimer),
            ringing = sum(cn.peaks.high.cor$ringing)
            )

    cn.stat.list[[pol]]$`Formula assignment` <-
      list( formula.matched = sum(cn.peaks.high.cor.formula$pave_formula_matched),
            formula.non.matched =sum(!cn.peaks.high.cor.formula$pave_formula_matched)    )


    cn.stat.list[[pol]] <- lapply(cn.stat.list[[pol]],function(x){
      do.call(rbind,x)%>%
        `colnames<-`(pol)%>%
        as.data.frame()%>%
        tibble::rownames_to_column("PAVE annotation")
    })%>%
      rbindlist(idcol = "PAVE FUN")
  }


  ###
  {
    cn.stat.df <- cn.stat.list$Positive%>%
      dplyr::mutate(Negative = cn.stat.list$Negative$Negative)%>%
      dplyr::mutate(tmp.ratio = Positive/Positive[1],
                    tmp.ratio = num2percent(tmp.ratio),
                    #Positive = paste0(Positive,"(",tmp.ratio,")"),

                    tmp.ratio = Negative/Negative[1],
                    tmp.ratio = num2percent(tmp.ratio),
                    #Negative = paste0(Negative,"(",tmp.ratio,")"),
                    `PAVE annotation` = c(
                      "total peaks number",
                      "peaks in procedure blank",
                      "other peaks without labeling",
                      "labeling but ρ<0.75",
                      "logical labeling (i.e., biological)",
                      "isotopes",
                      "dimer or double charge",
                      "adducts(assigned using same polarity mode)",
                      "adducts(assigned only using opposite polarity mode)",
                      "too low C count for mass",
                      "ringing peaks",
                      "formula match to metabolite",
                      "no formula match in database"
                    )
                    )%>%
      dplyr::select(-tmp.ratio)

   # cn.stat.df%>%
   #   gt::gt(rowname_col = "PAVE annotation", groupname_col = "PAVE FUN")

    return(cn.stat.df)

  }

  ### SNR
  if(F) {

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

PAVE_find_xcms_CN <- function(xcms.xcms, rt.tol = 20, ppm= 10 ,
                              BPPARAM = SnowParam(workers = 6,progressbar = T) ){



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


  ### find CN candidate
  {
    cn.list <- BiocParallel::bplapply(
      seq_along(xcms.fdf$feature_id),
      #1:300,
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
            pave_seed = paste0(pave_seed,this.fid,""),
            pave_CN = paste0("C",C_count,"N",N_count,"")
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
        cn.comb.list <- list()
        for (i.cn in 1:nrow(cn.comb)) {
          this.c <- cn.comb$C[i.cn]
          this.n <- cn.comb$N[i.cn]
          all.form <- c("C0N0",paste0("C0N",this.n,""),paste0("C",this.c,"N0"),paste0("C",this.c,"N",this.n,""))
          if (all(all.form %in% this.fdf$pave_CN) ) {

            cn.ft <-this.fdf%>%
              dplyr::filter(pave_CN %in% all.form)%>%
              dplyr::group_by(pave_CN)%>%
              dplyr::slice_max(peakMaxo)%>%
              dplyr::ungroup()
            cn.comb.list[[i.cn]] <- cn.ft
            m.detected <- xcms.val[cn.ft$feature_id,xcms.pave.sample$sampleNames]
            colnames(m.detected) <- xcms.pave.sample$sample.type
            rownames(m.detected) <- cn.ft$pave_CN
            m.detected <- m.detected/m.detected[1,1]
            m.ideal <- get_ideal_CN_ratio(this.c,this.n)%>%t
            m.ideal <- m.ideal[rownames(m.detected),colnames(m.detected)]

            p.cor <- cor(as.vector(m.detected),as.vector(m.ideal))
            cn.comb$p.cor[i.cn] <- p.cor

          }
        }


        ### filter xcms.fdf to return
        {
          if (any(cn.comb$p.cor > 0,na.rm =T)) {
            cn.fdf <- cn.comb.list[[which.max(cn.comb$p.cor)]]
            cn.fdf$pave_cor <- max(cn.comb$p.cor,na.rm  =T)

            message_with_time(this.fid," Pattern: ",cn.fdf$pave_CN[-1],"; Cor = ",cn.fdf$pave_cor[1])
            return(cn.fdf)
          }

        }


        #


      },
      xcms.fdf = xcms.fdf,rt.tol =rt.tol,ppm = ppm,xcms.val=xcms.val,xcms.pave.sample=xcms.pave.sample,
      BPPARAM = BPPARAM
    )

    cn.list <- cn.list[!sapply(cn.list,is.null)]
    names(cn.list) <- sapply(cn.list,function(x) x$feature_id[1])
  }


  ### save to xcms fdf
  if(F){

    for (i.cnl in seq_along(cn.list)) {
      this.fdf <- cn.list[[i.cnl]]
      xcms.fdf[this.fdf$feature_id,]$pave_seed <-
        paste0(xcms.fdf[this.fdf$feature_id,]$pave_seed , this.fdf$pave_seed)
      xcms.fdf[this.fdf$feature_id,]$pave_CN <-
        paste0(xcms.fdf[this.fdf$feature_id,]$pave_CN , this.fdf$pave_CN)
      xcms.fdf[this.fdf$feature_id[1],]$pave_cor <- this.fdf$pave_cor[1]

    }
    featureDefinitions(xcms.xcms)$pave_seed <- xcms.fdf$pave_seed
    featureDefinitions(xcms.xcms)$pave_CN <- xcms.fdf$pave_CN
    featureDefinitions(xcms.xcms)$pave_cor <- xcms.fdf$pave_cor

  }



  return(cn.list)
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
  rownames(m) <- c("S12C14N","S13C14N","S12C15N","S13C15N")

  if(N==0) m <- m[,1:2]+m[,3:4]

  return(m)
}


get_PAVE_LowC_cutoff <- function( c_max = 100){

  hmdb.cp <- MSdb:::get_HMDB_Compound_DF()
  atom.count <- MSCC:::chemform_parse(hmdb.cp$chemform)

  hmdb.pave.stat <- data.frame(
    chemform = hmdb.cp$chemical_formula,
    mass = as.numeric(hmdb.cp$monisotopic_molecular_weight),
    c.count =atom.count[,"C"]
  )%>%
    dplyr::filter(!is.na(mass),!is.na(c.count))


  PAVE_LowC_cutoff <- data.frame(
    c.count = 1:c_max,
    mass_min = NA,
    mass_max = NA
  )
  for (i.c in 1:c_max) {

    c.count.cp <- hmdb.pave.stat%>%
      dplyr::filter(c.count >= i.c-2,c.count <= i.c +2)

    c.c.m.q <- quantile(c.count.cp$mass,c(0.01,0.99))
    PAVE_LowC_cutoff$mass_min[i.c] <- c.c.m.q[1]
    PAVE_LowC_cutoff$mass_max[i.c] <- c.c.m.q[2]
  }

  return(PAVE_LowC_cutoff)

}

integrate_anntation <- function(str1,str2,sep = ";"){

  str_sep <- ifelse(str1==""|is.na(str1)|str2 == "","",sep)

  paste0(str1,str_sep,str2)

}
