




annotation_by_MSdb <- function(ms.ana,
                               ion_mode = "positive",
                               database_to_match = "all"){
  message(Sys.time(), " Annotating...")
  ### param
  {
    mz.ppm <- 20
    rt.tol <- 30

  }

  ### feature spectra
  {
    xcms.xcms <- ms.ana[[paste0("xcms.",ion_mode)]]
    norm_fun <- function(z, ...) {
      z[, "intensity"] <- z[, "intensity"] /
        max(z[, "intensity"], na.rm = TRUE) * 100
      z
    }
    feature.sp <- ms.ana[[paste0("Spectra.",ion_mode)]][["Spectra"]]%>%
      filterEmptySpectra()%>%
      addProcessing(norm_fun)%>%
      filterIntensity(c(3,Inf))%>%
      applyProcessing()

    feature.sp <- split(feature.sp , feature.sp$feature_id)

    feature.definition <-featureDefinitions(xcms.xcms)
    feature.definition$feature_id <- rownames(feature.definition)
    feature.definition$Spectra  <- sapply(feature.definition$feature_id,function(x)feature.sp[[x]])
  }

  ### database to match

  {
    data("Spectra_database",package = "MSdb")

    if (database_to_match == "all") {
      database_to_match <- unique(Spectra_database$database_origin)
    }else if( !database_to_match %in%Spectra_database$database_origin){
      message(database_to_match ," not exesited in MSdb, using all")
      database_to_match <- unique(Spectra_database$database_origin)
    }
    Spectra_database <- Spectra_database %>%
      filterPolarity(ifelse(ion_mode == "positive",1,0))%>%
      filterEmptySpectra%>%
      `[` (.$database_origin %in% database_to_match )
    Spectra_database$rtime[Spectra_database$database_origin != "in_house"] <- NA

  }


  ### candidate lib
  {
    feature.mz_rt <- data.frame(mz = feature.definition$mzmed,
                                rt = feature.definition$rtmed)
    lib.precursormz <- precursorMz(Spectra_database)
    lib.rtime <- rtime(Spectra_database)
    lib.candidate <- apply(feature.mz_rt,1,function(x){

      mz.hit <- abs( lib.precursormz-unlist(x["mz"])) < Spectra::ppm(unlist(x["mz"]),mz.ppm)
      rt.hit <- abs( lib.rtime -unlist(x["rt"])) < rt.tol
      rt.hit[is.na(rt.hit)] <- T
      which(mz.hit & rt.hit  )

    })
    feature.definition$candidate <- lapply(lib.candidate,  function(x){
      sp <- Spectra_database[x]
      if (length(sp)== 0) {
        return(NULL)
      }else{
        return(sp)
      }
    })

    ### annotation
    {
      message(Sys.time()," Matching spectra...")
      feature.annotation <- bplapply(1:nrow(feature.definition),mz.ppm,rt.tol, FUN = function(x,mz.ppm,rt.tol){


        to.return <- data.frame(
          feature.id = feature.definition$feature_id[x],
          mz = feature.definition$mzmed[x],
          rt = feature.definition$rtmed[x],
          ref.mz = NA,
          ref.rt = NA,
          score = 0,
          compound = NA,
          adduct = NA,
          inchikey = NA,
          kegg.id = NA,
          origin = NA,
          sp.exp = NA,
          sp.lib = NA
        )

        #mz.ppm <- 20
        #rt.tol <- 30
        sp.exp <- feature.definition$Spectra[[x]]
        sp.lib <- feature.definition$candidate[[x]]

        mz.error <- abs(to.return$mz - sp.lib$precursorMz)
        mz.score <- 1 - mz.error/Spectra::ppm(to.return$mz , mz.ppm)

        rt.error <- abs(to.return$rt - sp.lib$rtime)
        rt.score <- 2- rt.error/rt.tol
        rt.score[is.na(rt.score)] <- 0

        if(length(sp.exp)>0 & length(sp.lib)>0){
          sp.score <- Spectra::compareSpectra(sp.exp,sp.lib)
          sp.score[is.nan(sp.score)] <- 0
        }else{
          sp.score <- rep(0,length(sp.lib))
        }

        if (length(sp.exp) ==0 ) {
          score <- mz.score + rt.score + sp.score
        }else{
          score <- mz.score * 0 + rt.score + sp.score
          to.return$sp.exp <- sp.exp
        }

        if (length(sp.lib) == 0) {
        }else{

          sp.lib <- sp.lib[which.max(score)]
          score <- sp.score[which.max(score)]
          to.return$ref.mz <- sp.lib$precursorMz
          to.return$ref.rt <- sp.lib$rtime
          to.return$score <- score
          to.return$compound <- sp.lib$Compound_name
          to.return$Lipid_subclass <-sp.lib$Lipid_subclass
          to.return$adduct <- sp.lib$adduct
          to.return$inchikey <- sp.lib$inchikey
          to.return$kegg.id <- sp.lib$KEGG_id
          to.return$origin <- sp.lib$database_origin
          to.return$sp.lib <- sp.lib
        }
        return(to.return)



      },BPPARAM = SerialParam(
        progressbar = T))

      annotation.table <- lapply(feature.annotation,function(x){
        x[1:11]
      }) %>%data.table::rbindlist(fill = T) %>%as.data.frame()


    }


  }
  ### save and return
  {
    ms.ana[[paste0("annotation.",polarity)]][["annotation"]] <- feature.annotation
    ms.ana[[paste0("annotation.",polarity)]][["annotation.table"]] <- annotation.table
    save(ms.ana , file = ms.ana$processing.info$project.info$ms.ana.file)

    return(ms.ana)
  }


}
