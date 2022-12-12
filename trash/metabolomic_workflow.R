#' @title metabolomic_workflow
#'
#' @param project.dir
#' @param raw.data.dir
#'
#' @return
#' @export
#'
#' @examples
metabolomic_workflow <- function(project.dir = "d:/2022_07_05-Lirui/",
                                 raw.data.dir = "d:/2022_07_05-Lirui/Data/") {
  thread <- parallel::detectCores() - 1




  ### MS file convert and information
  {
    ms.ana <-
      export_sample_information_from_wiff(raw.data.dir = raw.data.dir,
                                          project.dir = project.dir)

    #ms.ana$sample.info <- readxl::read_excel(paste0(project.dir,"/sample.info.xlsx"))

    ms.ana <- ms_convert(ms.ana)

  }

  ### xcms processing
  {
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "positive")
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "negative")
    ms.ana <- extract_spectra(ms.ana)

  }
  ### annotation
  {
    ms.ana <- annotation_by_database(ms.ana, polarity = "positive")
    ms.ana <- annotation_by_database(ms.ana, polarity = "negative")

  }

  ### Statistic
  {
    ms.ana <- get_feature(ms.ana )
    ms.ana <- get_unique_compound(ms.ana)
    sta_workflow(ms.ana )

  }
}

lipidomic_workflow <- function(project.dir = "d:/2022_07_05-Lirui/",
                               raw.data.dir = "d:/2022_07_05-Lirui/Data/"){
  thread <- parallel::detectCores() - 1

  ### MS file convert and information
  {
    ms.ana <-
      export_sample_information_from_wiff(raw.data.dir = raw.data.dir,
                                          project.dir = project.dir)

    #ms.ana$sample.info <- readxl::read_excel(paste0(project.dir,"/sample.info.xlsx"))

    ms.ana <- ms_convert(ms.ana)

  }

  ### xcms processing
  {
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "positive")
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "negative")
    ms.ana <- extract_spectra(ms.ana)

  }
  ### annotation
  {

    ms.ana <- annotation_by_MSdb(ms.ana = ms.ana,
                                 ion_mode = "positive",
                                 database_to_match = "LipidBlast_from_MSDIAL")

    ms.ana <- annotation_by_MSdb(ms.ana = ms.ana,
                                 ion_mode = "negative",
                                 database_to_match = "LipidBlast_from_MSDIAL")
  }


}




metabolomic_workflow_single_file <-
  function(xcms.xcms ,
           polarity = "positive"){

  centwave.param <- CentWaveParam(ppm = 20,
                                  peakwidth = c(5,30),
                                  snthresh = 10,
                                  prefilter = c(3,100))
  xcms.xcms <- findChromPeaks(xcms.xcms , param = centwave.param)

  mz_width <- chromPeaks(xcms.xcms)%>%
    as.data.frame()%>%
    mutate(mz_error = mzmax - mzmin )%>%
    pull(mz_error)%>%
    max%>%
    `*`(0.3)

  peak.density.param <- PeakDensityParam(sampleGroups = "a",
                                         minFraction = 0.5,bw = 30,
                                         binSize = 0.05)
  xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)
  xcms.spectra <- featureSpectra(xcms.xcms,return.type = "Spectra")
  Spectra::combineSpectra(xcms.spectra,
                          f = xcms.spectra$feature_id,
                          peaks = "intersect",
                          minProp = 0.5,
                          ppm = 20)->xcms.spectra
  norm_fun <- function(z, ...) {
    z[, "intensity"] <- z[, "intensity"] /
      max(z[, "intensity"], na.rm = TRUE) * 100
    z
  }
  feature.sp <- xcms.spectra%>%
    filterEmptySpectra()%>%
    addProcessing(norm_fun)%>%
    filterIntensity(c(3,Inf))%>%
    applyProcessing()
  feature.sp <- split(feature.sp , feature.sp$feature_id)

  feature.definition <-featureDefinitions(xcms.xcms)
  feature.definition$feature_id <- rownames(feature.definition)
  feature.definition$Spectra  <- sapply(feature.definition$feature_id,function(x)feature.sp[[x]])

  ### database to match

  {
    database.file = "C:\\Users\\91879\\OneDrive\\Documents\\Code\\R\\Projecct\\2022.1.17_Compounds.database\\Spectra.integrated.database.integration.2022_02_12.Rdata"
    load(database.file)

    spectra.database %<>%
      filterPolarity(ifelse(polarity == "positive",1,0))%>%
      filterEmptySpectra%>%
      `[` (.$database %in% c("HMDB","KEGG","MassBnak","MoNA","inHouse"))


  }
  ### candidate lib
  {
    feature.mz_rt <- data.frame(mz = feature.definition$mzmed,
                                rt = feature.definition$rtmed)
    lib.precursormz <- precursorMz(spectra.database)
    lib.rtime <- rtime(spectra.database)
    mz.ppm <- 20
    rt.tol <- 10
    lib.candidate <- apply(feature.mz_rt,1,function(x){

      mz.hit <- abs( lib.precursormz-x["mz"]) < Spectra::ppm(x["mz"],mz.ppm)
      rt.hit <- abs( lib.rtime -x["rt"]) < rt.tol
      rt.hit[is.na(rt.hit)] <-T
      which(mz.hit & rt.hit  )

    })
    feature.definition$candidate <- lapply(lib.candidate,  function(x){
      sp <- spectra.database[x]
      if (length(sp)== 0) {
        return(NULL)
      }else{
        return(sp)
      }
    })


  }
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
        to.return$compound <- sp.lib$name
        to.return$adduct <- sp.lib$adduct
        to.return$inchikey <- sp.lib$inchikey
        to.return$kegg.id <- sp.lib$kegg.id
        to.return$origin <- sp.lib$database
        to.return$sp.lib <- sp.lib
      }
      return(to.return)



    },BPPARAM = SerialParam(
      progressbar = T))

    annotation.table <- lapply(feature.annotation,function(x){
      x[1:11]
    }) %>%data.table::rbindlist(fill = T) %>%as.data.frame()


  }
  annotation.table

  }



metabolomics_workflow_for_QE_fullscan_both <- function(){


  MS_dev_QE <- MSdev(rawDataDir = "d:/2022.9.1_QE.Test/rawData",
                     projectDir = "d:/2022.9.1_QE.Test",
                     experimentInfo = MS_Experiment[5])
  MS_dev_QE <- checkSampleInfo(MS_dev_QE)
  MS_dev_QE <- msConvert(MS_dev_QE)
  MS_dev_QE <- xcmsProcessing_fullscan_DDA(MS_dev_QE)
  MS_dev_QE <- extractSpectra_fullscan_DDA(MS_dev_QE)
  MS_dev_QE <- featureSpectra_fullscan_DDA(MS_dev_QE)
  MS_dev_QE <- featureCandidate(MS_dev_QE,mz.ppm = 20)
  MS_dev_QE <- annotateMSdev(MS_dev_QE)
  MS_dev_QE <- dropSpectra(MS_dev_QE)


}


