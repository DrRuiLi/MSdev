


#' @title MSdev input and output
#' @description
#' save `MSdev` using `qs::qsave()` and  qs::qread()`
#'
#' @describeIn MSdev_IO MSdev_save
#'
#' this function ....
#'
#' @param object MSdev
#'
#' @return MSdev
#' @export
#'
MSdev_save <- function(object,file = object@projectInfo$MSdevFile){
  MSdev <- object
  save.dir <- dirname(file)
  if (!dir.exists(save.dir)) {
    dir.create(save.dir,recursive = T)
  }
  #save(MSdev, file =  object@projectInfo$MSdevFile)
  qs::qsave(MSdev, file =  file)
  invisible(MSdev)
}



#' @title Load an MSdev object from a file
#' @description Load an MSdev object from a file using `qs::qread()`.
#' @describeIn MSdev_IO load
#' @param file_to_load file path
#' @return MSdev
#' @export
#'
MSdev_load <- function(file_to_load){

  qs::qread(file_to_load)


}

.updateProjectInfoFromSampleInfo <- function(object){

  object@sampleInfo -> sampleInfo

  if ("scanType"%in% colnames(sampleInfo)) {
    object@projectInfo$msAcquisition <- unique(sampleInfo$scanType)
  }
  if ("msLevels"%in% colnames(sampleInfo)) {
    object@projectInfo$msLevel <- unique(sampleInfo$msLevels)
  }
  if ("polarity"%in% colnames(sampleInfo)) {
    object@projectInfo$polarity <- unique(sampleInfo$polarity)
  }
  if ("manufacturer"%in% colnames(sampleInfo)) {
    object@projectInfo$msManufacturer <- unique(sampleInfo$manufacturer)
  }
  if ("model"%in% colnames(sampleInfo)) {
    object@projectInfo$msModel <- unique(sampleInfo$model)
  }
  object@projectInfo$sampleCount <- sampleInfo%>%
    dplyr::filter(sample.type=="Sample")%>%
    dplyr::pull(sample.name)%>%unique()%>%length()
  object@projectInfo$rawDataFileCount <- sum(!is.na(sampleInfo$raw.files))
  p.index <- c("-1"="Unknow","0"="Negative","1"="Positive")
  object@projectInfo$rawDatafiles <-table(p.index[sampleInfo$polarity],sampleInfo$sample.type)
  object
}



MSdev_get_MSinfo <- function(object){


  ### define acquire Type
  ### note, these model string are identified by mzR
  {
    HRMS <- c("Q Exactive Plus","TripleTOF 6600","Orbitrap IQ-X","Q Exactive",
              "Orbitrap Exploris 480","Orbitrap Astral",
              "Thermo Electron instrument model")
    TQMS <- c("TSQ Quantis")
    model.df <- data.frame(
      model = c(HRMS,TQMS),
      type = c( rep("HRMS",length(HRMS)),
                rep("TQMS",length(TQMS)))
    )

  }



  ### update SampleInfo
  if (
    !"manufacturer"%in%colnames(object@sampleInfo)|
    any(is.na(object@sampleInfo$manufacturer))
  ) {
    sampleInfo <-object@sampleInfo%>%
      dplyr::mutate(get_MSinfo_mzR(msData.files),
                    msType = model.df$type[match(model,model.df$model)],
                    scanType = case_when(msType=="HRMS"&msLevels=="1" ~ "FS",
                                         msType=="HRMS"&msLevels=="2" ~ "DDA",
                                         msType=="TQMS" ~ "MRM",
                                         T~NA),
                    xcmsProcessing = case_when(scanType=="FS"~ "MS1",
                                               scanType=="DDA" ~ "Both",
                                               scanType=="MRM" ~ "MRM",
                                               T~NA))%>%
      dplyr::group_by(sample.name)%>%
      dplyr::mutate(polarity_paired = case_when(
        "0"%in%polarity&"1"%in%polarity~T,
        "0;1"%in% polarity~T,
        T~F
      ),
      .after = polarity)%>%
      dplyr::ungroup()
    sampleInfo -> object@sampleInfo

  }
  return(object)

}


#' @title Add new sample files to MSdev object
#' @description Add raw data files from a directory to the MSdev object, converting them as needed.
#' @describeIn MSdev_workflow add samples
#' @param object MSdev
#' @param raw.data.dir file path
#' @return MSdev
#' @export
MSdev_add_sample <- function(object,
                             raw.data.dir = object@projectInfo$rawDataDir){
  sample.info <- object@sampleInfo
  sample.info.new <- get_MS_sampleinfo(raw.data.dir,
                                       rawDataFormat = object@projectInfo$rawDataFormat,
                                       verbose = T)
  if (all(sample.info.new$raw.files %in% sample.info$raw.files)) {

    return(object)

  }
  sample.info.new <- sample.info.new%>%
    dplyr::filter(!raw.files%in% sample.info$raw.files)%>%
    dplyr::mutate(sample.name = case_when(
      sample.name %in% sample.info$sample.name ~
        paste0(sample.name,"_",1:n()+nrow(sample.info)),
      T~ sample.name
    ))
  message("Get ",nrow(sample.info.new)," samples")

  object@sampleInfo <- sample.info%>%
    bind_rows(sample.info.new)
  message("MSconvert ",nrow(sample.info.new)," samples...")
  object <- MSdev_msConvert(object)
  object <- .updateProjectInfoFromSampleInfo(object )
  show(object)
  object

}


#' @title Process MS data using xcms
#' @description Perform peak detection and grouping using xcms functions based on acquisition mode (FS/DDA/MRM).
#' @describeIn MSdev_workflow use xcms to Processing data
#' @param object MSdev
#' @param ... additional arguments passed to xcms functions
#' @return MSdev
#' @export
MSdev_xcmsProcessing <- function(object,...){

  MS.mode <- object@projectInfo$msAcquisition

  if ("FS" %in% MS.mode|"DDA" %in% MS.mode) {
    object <- MSdev_get_xcms(object)
    object <- xcmsProcessingMSdev.DDA(object,...)
    return(object)
  }

  if ("MRM" %in% MS.mode) {
    object <-  xcmsProcessingMSdev.MRM(object)
    return(object)
  }

  return(object)



}


MSdev_get_xcms <- function(object){

  polarity.index <- c("0" = "Negative",
                      "1"="Positive")
  for (i in c(0,1)) {
    sample.info.polarity <- object@sampleInfo%>%
      dplyr::filter(grepl(i,polarity))
    if (!nrow(sample.info.polarity)) {
      xcms.xcms <-NA
      next
    }else{
      xcms.xcms <- readMSData(sample.info.polarity$msData.files,
                              mode = "onDisk")
      xcms.xcms <- filterPolarity(xcms.xcms,i)
    }
    Biobase::pData(xcms.xcms ) <-
      cbind(Biobase::pData(xcms.xcms),
            sample.info.polarity)
    polarity.tag <- paste0(polarity.index[as.character(i)],"MS1")
    xcms.xcms -> object@xcmsData[[polarity.tag]]

  }
  return(object)


}

#' @title Process DDA data using xcms
#' @description Perform peak detection on DDA data using xcms, grouping features across samples.
#' @param object MSdev
#' @param ... additional arguments passed to xcms functions
#' @return MSdev
#' @export
xcmsProcessingMSdev.DDA <- function(object,...){

  xcms.param <- get_MSdev_param(object )
  sampleInfo <-dplyr::filter(object@sampleInfo,
                                xcmsProcessing %in% c("MS1","Both")
  )%>%
    dplyr::filter(!is.na(msData.files))

  polarity.index <- c("0" = "Negative",
                      "1"="Positive")
  for (i in c(0,1)) {
    sample.info.polarity <- sampleInfo%>%
      #dplyr::filter(polarity %in% c(i,-1,"0;1"))%>%
      dplyr::filter(grepl(i,polarity))
    polarity.tag <- paste0(polarity.index[as.character(i)],"MS1")
    if (!nrow(sample.info.polarity)) {
      xcms.xcms <-NA
      next
    }
    xcms.xcms <- filterFile(object@xcmsData[[polarity.tag]],
      which(Biobase::pData(object@xcmsData[[polarity.tag]])$sample.name%in% sample.info.polarity$sample.name)
    )

    xcms.xcms <-
      xcmsProcessingMS1(xcms.xcms = xcms.xcms,
                        ion_mode = i,
                        xcms_param  = xcms.param,
                        ...
    )

    xcms.xcms <- xcms_get_feature_stat(xcms.xcms )
    xcms.xcms -> object@xcmsData[[polarity.tag]]



  }


  object




}


xcmsProcessingMSdev.MRM <- function(object){

  xcms.param <- get_MSdev_param(object )
  sampleInfoPos <-object@sampleInfo%>%
    dplyr::filter(!is.na(msData.file.positive))

  object@xcmsData$positiveMS1  <- xcmsProcessingMRM(msDataFiles = sampleInfoPos$msData.file.positive,
                                                    peaksGroup =sampleInfoPos$sample.type,
                                                    centWaveParam = xcms.param$findpeak.param
  )

  Biobase::pData(object@xcmsData$positiveMS1) <- cbind(Biobase::pData(object@xcmsData$positiveMS1),
                                                       sampleInfoPos)



}

extractFeature <- function(object){

  object@xcmsData$positiveFeature <- as.data.frame(featureDefinitions(object@xcmsData$PositiveMS1))
  object@xcmsData$negativeFeature <- as.data.frame(featureDefinitions(object@xcmsData$NegativeMS1))
  object

}

get_MSdev_param <- function(object){


  if (!is.null(object@processingInfo$MSdevParam )) {
    return(object@processingInfo$MSdevParam )
  }

  if (object@experimentInfo@General$Name == "Metabolomics_YLF") {

    fpp <- xcms::MassifquantParam(ppm = 5,
                            peakwidth = c(10,100),
                            mzCenterFun  = "wMeanApex3",
                            snthresh = 10,
                            prefilter = c(5,1000),
                            verboseColumns=T,
                            withWave = T)
    gpp <- xcms::PeakDensityParam("A",bw = 5,
                            minFraction = 0.3,binSize = 0.002)

    ### temp for Astral
    if ("Orbitrap Astral"  %in% object@projectInfo$msModel ) {

      fpp <- xcms::MassifquantParam(ppm = 5,
                              peakwidth = c(5,60),
                              mzCenterFun  = "wMeanApex3",
                              snthresh = 10,
                              prefilter = c(3,100),
                              verboseColumns=T,
                              withWave = T)
    }

    msdev.param <- list(findChromPeaks = fpp,
                        groupChromPeaks = gpp)
    return(msdev.param)
  }

  msdev.param <- get_MSdev_xcms_param_by_exp(object)
  return(msdev.param)




}

MSdev_set_param <- function(object,
                            findChromPeaks = xcms::CentWaveParam(),
                            groupChromPeaks = xcms::PeakDensityParam("A")
                            ){

  msdev.param <- list(findChromPeaks = findChromPeaks,
                      groupChromPeaks = groupChromPeaks)
  object@processingInfo$MSdevParam  <- msdev.param

  return(object)

}



get_MSdev_xcms_param_by_exp <- function(object){



  MS.mode <- object@projectInfo$msAcquisition
  MS.instru <-object@projectInfo$msModel
  MS.LC.rate <- object@experimentInfo@Chroma_gradient[[1]]$Flow_rate%>%mean
  MS.LC.time <- object@experimentInfo@Chroma_gradient[[1]]$time%>%max
  cwp <- xcms::CentWaveParam(fitgauss = F,verboseColumns = T)

  ### ppm
  cwp@ppm <- switch(MS.instru,
                    "Q Exactive Plus" = 10,
                    "TripleTOF 6600" = 25,
                    20)

  cwp@peakwidth <-switch(as.character(MS.LC.rate),
                         "0.5" = c(5,20),
                         "0.3" = c(15,200),
                         c(5,20))

  cwp@snthresh <- switch(MS.instru,
                         "Q Exactive Plus" = 100,
                         "SCIEX TripleTOF 6600" = 100,
                         "Thermo Quantis" = 0,
                         100)
  cwp@prefilter <- switch(MS.instru,
                          "Q Exactive Plus" = c(3,100),
                          "SCIEX TripleTOF 6600" = c(3,100),
                          "Thermo Quantis" = c(3,100),
                          c(3,100))

  ### group peaks param
  {
    gpp <- xcms::PeakDensityParam(sampleGroups = "A",
                            bw = 5,
                            minFraction = 0.6,
                            binSize = 0.015)


  }

  msdev.param <- list(findChromPeaks = cwp,
                      groupChromPeaks = gpp)


  return(msdev.param)


}


#' @title extractSpectra_fullscan_DDA
#' @description extract all MS2 Spectra from `object@sampleInfo$msDataFile` which `sampleInfo$xcmsProcessing` %in% % c("Both","MS2"),
#' return store in `object@spectra$positiveMS2`
#' @param object a `MSdev` object
#'
#' @return  a `MSdev` object
#' @export
#'

extractSpectra_fullscan_DDA <- function(object){

  sampleInfo <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing %in% c("Both","MS2"))
  if (nrow(sampleInfo)==0) {
    spectra.pos <- Spectra::Spectra()
    spectra.neg <- Spectra::Spectra()
  } else {
    spectra.pos <- Spectra::Spectra(na.omit(sampleInfo$msData.file.positive),backend = Spectra::MsBackendDataFrame())%>%
      filterMsLevel(2)

    spectra.neg <- Spectra::Spectra(na.omit(sampleInfo$msData.file.negative),backend = Spectra::MsBackendDataFrame())%>%
      filterMsLevel(2)
  }


  object@spectra$positiveMS2 <- spectra.pos
  object@spectra$negativeMS2 <- spectra.neg
  return(object)

}




#' @title featureSpectra_fullscan_DDA
#' @description extrat spectra from `MSdev@spectra` according to mz and rt of feature,
#' extracted spectra store in `object@spectra$positiveFeatureMS2` and `object@spectra$negativeFeatureMS2`,
#' a list contain `Spectra` object of each feature, empty `Spectra` with precursorMz and rtime
#'
#' @param object a `MSdev` object
#'
#' @return a `MSdev` object
#'
#' @export
#'

featureSpectra_fullscan_DDA <- function(object){

  .matchSP <- function(x,spectras,
                       mz_ppm = 10,
                       rt_tol = 30){
    mz <- x$mzmed
    rt <- x$rtmed
    mzError <- abs((mz - spectras$precursorMz)/mz*1e6)
    rtError <- abs(rt- spectras$rtime)
    matchedspectras_id <- which(mzError < mz_ppm &rtError < rt_tol)
    if (length(matchedspectras_id)==0) {
      matchedspectras <- makeSpectra(mz,rt)
      matchedspectras$feature_id <- rownames(x)
      return(matchedspectras)
    }else{
      matchedspectras <- spectras[matchedspectras_id]
      matchedspectras$feature_id <- rownames(x)
      return(matchedspectras)
    }


  }

  object@spectra$positiveFeatureMS2 <- apply(featureDefinitions(object@xcmsData$positiveMS1), 1, .matchSP , object@spectra$positiveMS2 )
  object@spectra$negativeFeatureMS2 <- apply(featureDefinitions(object@xcmsData$negativeMS1), 1, .matchSP , object@spectra$negativeMS2 )


  object
}


#' @title featureCandidate
#' @description match feature with database by mz, return all spectra matched in a list  splited by feature
#'
#' @param object a `MSdev` object
#' @param mz.ppm mz error
#' @param spectraDatabase databse path
#'
#' @return a `MSdev` object
#' @export
#'

featureCandidate<- function(object,mz.ppm = 20,
                            spectraDatabase =
                              "C:\\Users\\91879\\OneDrive\\Documents\\Code\\R\\Projecct\\2022.1.17_Compounds.database\\Spectra.integrated.database.integration.2022_02_12.Rdata")
{

  ### load spectra database
  Spectra_database <- load_as_var(spectraDatabase)

  ### match exact mz to database, extract all matched record as candidate
  {

    .matchMz <- function(xcmsFeature,spectraDB){
      feature.mz_rt <- data.frame(mz = xcmsFeature$mzmed,
                                  rt = xcmsFeature$rtmed)
      lib.precursormz <- ProtGenerics::precursorMz(spectraDB)
      lib.rtime <- rtime(spectraDB)
      lib.candidate <- apply(feature.mz_rt,1,function(x){

        mz.hit <- abs( lib.precursormz-x[["mz"]]) < x[["mz"]]*mz.ppm /1e6
        which(mz.hit   )

      })
      #sum(sapply(lib.candidate, length)>0)
      featureCandidate.list <- lapply(lib.candidate,  function(x){
        sp <- spectraDB[x]
        if (length(sp)== 0) {
          return(NULL)
        }else{
          return(sp)
        }
      })
      return(featureCandidate.list)

    }
    object@projectInfo$MSDB_path <- spectraDatabase
    object@annotation$positiveCandidate <- .matchMz(object@xcmsData$positiveFeature,
                                                    ProtGenerics::filterPolarity(Spectra_database,1))
    object@annotation$negativeCandidate <- .matchMz(object@xcmsData$negativeFeature,
                                                    ProtGenerics::filterPolarity(Spectra_database,0))
  }

  ### according to isotopes, from all candidate, expand possible isotope
  {
    .expand_iso <- function( featureCandidate.list,xcms.xcms){

      featuredef <- featureDefinitions(xcms.xcms)%>%as.data.frame()
      featureval <- featureValues(xcms.xcms,missing = "rowmin_half")
      .match.featurecandidate <- function(x,featuredef,featureval){
        if (is.null(x)) {
          return(NULL)
        }else{

          iso.all <- list()
          for (i in 1:length(x)) {
            iso.table <- MSCC::chemform_isotopes_pattern_enviPat(
              MSCC::chemform_adduct(x$formula[i],x$adduct[i]) )%>%
              MSCC:::match_isotopes_to_featuredef(featuredef ,rt.tol = 10)%>%
              MSCC:::match_isotopes_to_featureval(featureval)%>%
              dplyr::mutate(MSDB_id = x$MSDB_id[i])

            iso.all[[i]] <- iso.table

          }
          iso.all <- do.call("rbind",iso.all)

        }
        iso.all


      }

      bplapply(featureCandidate.list,
               .match.featurecandidate,BPPARAM = SnowParam(progressbar = T),featuredef,featureval) -> expanded.spectra




    }


  }


  object

}

annotateMSdev <- function(object){

  .annotateMSdev <- function(featureMS2 , candidate){

    BiocParallel::bplapply(1:length(featureMS2),function(i){
      annotateSpectraMSdb(featureMS2[[i]],candidate[[i]])
    },BPPARAM = BiocParallel::SerialParam(
      progressbar = T))
  }
  object@annotation$positiveAnnotation <- .annotateMSdev(featureMS2 = object@spectra$positiveFeatureMS2,
                                                         candidate = object@annotation$positiveCandidate)
  object@annotation$negativeAnnotation <- .annotateMSdev(featureMS2 = object@spectra$negativeFeatureMS2,
                                                         candidate = object@annotation$negativeCandidate)



  return(object)
}

dropSpectra <- function(object){

  object@spectra$positiveMS2 <- NULL
  object@spectra$negativeMS2 <- NULL
  object@annotation$positiveCandidate <- NULL
  object@annotation$negativeCandidate <- NULL
  return(object)
}



adjustFeatureByIS <-function(object,to.adjust = "featureRaw"){

  object <- findISMSdev(object,corr.thred = 0.3)
  features <- object@statData[[to.adjust]]%>%
    dplyr::mutate(internal_standard =object@statData$featureRaw$internal_standard[match(
      feature_id ,object@statData$featureRaw$feature_id
    )],.before = qc_rsd )
  sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("Both","MS1"),
                  sample.type!= "GQC")

  feature.matrix <- features%>%
    column_to_rownames("feature_id")%>%
    dplyr::select(sample.info$sample.name)

  {### pos

    is.norm <- feature.matrix[features%>%
                                dplyr::filter(ion_mode == "positive",
                                              !is.na(internal_standard))%>%
                                pull(feature_id),]%>%
      apply(1,function(x){x/mean(x,na.rm = T)})%>%
      apply(1,mean)
    feature.matrix.pos <-feature.matrix[features%>%
                                          dplyr::filter(ion_mode == "positive")%>%
                                          pull(feature_id),]%>%
      t%>%
      `/`(is.norm)%>%
      t

    }
  {### neg

    is.norm <- feature.matrix[features%>%
                                dplyr::filter(ion_mode == "negative",
                                              !is.na(internal_standard))%>%
                                pull(feature_id),]%>%
      apply(1,function(x){x/mean(x,na.rm = T)})%>%
      apply(1,mean)
    feature.matrix.neg <-feature.matrix[features%>%
                                          dplyr::filter(ion_mode == "negative")%>%
                                          pull(feature_id),]%>%
      t%>%
      `/`(is.norm)%>%
      t

  }
  feature.matrix <- rbind(feature.matrix.neg,feature.matrix.pos)[features$feature_id,sample.info$sample.name]
  features[,sample.info$sample.name] <- feature.matrix


  features -> object@statData[["feature"]]

  return(object)





}

adjustFeatureByGQC <- function(msdev.object,to.adjust = "featureRaw"){

  sampleinfo <- msdev.object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("Both","MS1") )
  GQC.sampleinfo <- msdev.object@sampleInfo%>%
    dplyr::filter(sample.type == "GQC")

  sample.matrix <- msdev.object@statData[[to.adjust]]%>%
    column_to_rownames("feature_id")%>%
    dplyr::select(sampleinfo$sample.name)%>%t

  .adjust.fun <- function(x){
    fit.df <- data.frame(y =GQC.sampleinfo$QC.gradient.concentraion,
                         x = x[GQC.sampleinfo$sample.name])%>%
      dplyr::filter(!is.na(x))

    if (nrow(fit.df)<=1) {
      return(c(x, r2 = 0))

    }
    #fit <- lm(y~x , data = fit.df)
    fit <- e1071::svm(y~x , data = fit.df)
    fit.pred <- predict(fit,newdata = data.frame(x = x ))
    fit.pred <- fit.pred[names(x)]
    names(fit.pred) <- names(x)
    r2 <- cor(fit.df$y,predict(fit))^2
    return(c(fit.pred , r2 = r2))

  }



  adjusted.matrix  <- apply(sample.matrix,2, .adjust.fun)%>%t
  msdev.object@statData[["feature"]] <- msdev.object@statData[[to.adjust]]
  msdev.object@statData[["feature"]][ ,sampleinfo$sample.name] <-adjusted.matrix[,sampleinfo$sample.name]
  msdev.object@statData[["feature"]] <- msdev.object@statData[["feature"]]%>%
    dplyr::mutate(gqc_r2 = adjusted.matrix[,"r2"],.before = qc_rsd)

  msdev.object
}

findFeature <- function(object,
                        exact_mass =100,
                        retention_time = 100,
                        ppm = 10,rt.err = 10, ion_mode = 1 ){

  ion_mz <- exact_mass+ifelse(ion_mode==1 , 1.007825,-1.007825)
  ion_mode_char <- ifelse(ion_mode==1 , "positive","negative")

  rt.err <- ifelse(is.na(retention_time),Inf, rt.err)
  retention_time <- ifelse(is.na(retention_time),0, retention_time)

  feature <- object@statData$featureRaw
  feature_matched <- feature%>%
    dplyr::filter( ion_mode == ion_mode_char,
                   mz > ion_mz-ion_mz*ppm/1e6,
                   mz <  ion_mz+ion_mz*ppm/1e6)%>%
    dplyr::filter(rt > retention_time - rt.err,
                  rt < retention_time + rt.err)

  return(feature_matched)
}

#' @title findISMSdev
#' @description find features of internals standard listed in `object@experimentInfo@Internal_Standard`
#' by `Exact_mass` and `Retention_time` (if provide),
#' only `[M+H]` and `[M-H]` are considered. Correlation and intensity will be plot based on `object@statData[["featureRaw"]]`, please check.
#' A column "internal_standard" will be added in `object@statData[["featureRaw"]]`
#'
#' @param object MSdev
#' @param corr.thred cor
#'
#' @return MSdev
#' @export
#'

findISMSdev <- function(object ,to.adjust = "featureRaw",corr.thred = 0.6){

  internal.standard <- object@experimentInfo@Internal_Standard%>%as.data.frame()
  feature <-  object@statData[[to.adjust]]%>%
    dplyr::mutate( .before = qc_rsd,
                   internal_standard = NA)
  for (i in 1:nrow(internal.standard)) {
    ft.pos <- findFeature(object ,
                          exact_mass = internal.standard$Exact_mass[i],
                          retention_time = internal.standard$Retention_time[i]*60,
                          ppm = 10,rt.err = 10,ion_mode = 1)$feature_id

    ft.neg <- findFeature(object ,
                          exact_mass = internal.standard$Exact_mass[i],
                          retention_time = internal.standard$Retention_time[i]*60,
                          ppm = 10,rt.err = 10,ion_mode = 0)$feature_id

    feature <- feature%>%
      dplyr::mutate( .before = qc_rsd,
                     internal_standard = ifelse(feature_id %in% c(ft.pos,ft.neg),internal.standard$Compound_name[i],internal_standard))
  }



  feature.internal.standard <- feature %>%
    dplyr::filter(!is.na(internal_standard))

  { ### confirm by correlation

    sampleinfo <- object@sampleInfo%>%
      dplyr::filter(xcmsProcessing %in% c("Both","MS1"))%>%
      dplyr::filter(sample.type %in% c("QC","Sample","GQC"))

    feature.matrix <- feature.internal.standard%>%
      column_to_rownames("feature_id")%>%
      dplyr::select(sampleinfo$sample.name)%>%
      as.matrix()

    feature.matrix <- apply(feature.matrix,1,function(x){
      x[is.na(x)] <- min(x,na.rm = T)/10
      return(x)}
    )%>%t
    cor.matrix <- cor(t(feature.matrix))
    high.cor <-apply(cor.matrix, 1, mean)>corr.thred

    feature.internal.standard <- feature.internal.standard[high.cor,]
    feature.matrix <- feature.matrix[high.cor,]
    cor.matrix <- cor(t(feature.matrix))

    corrplot::corrplot(cor.matrix  ,is.corr = F,
                       tl.col = "black",
                       col = colorRampPalette(c("#0A3A70","white","#FF6666"))(100),
                       col.lim = c(min(cor.matrix),1))
    dir.create(paste0(object@projectInfo[["msDataDir"]],"/dataProcessing"))
    export::graph2pdf(file = paste0(object@projectInfo[["msDataDir"]],"/dataProcessing/Internal_Standard_Corr.pdf"),
                      width = nrow(cor.matrix)*0.7,height = nrow(cor.matrix)*0.7)
    openxlsx::write.xlsx(feature.internal.standard,
                         file = paste0(object@projectInfo[["msDataDir"]],"/dataProcessing/Internal_Standard.xlsx"))

    }
  feature <- feature%>%
    dplyr::mutate(internal_standard = ifelse(feature_id %in% feature.internal.standard$feature_id,
                                             internal_standard,NA))
  { # plot internal standard intensity
    p.list <- list()
    for (i in 1:nrow(feature.internal.standard)) {
      if (feature.internal.standard$ion_mode[i] == "positive") {
        plot_xcms_feature_intensity(object@xcmsData$positiveMS1,
                                    sub(feature.internal.standard$feature_id[i],pattern = "_pos",replacement = ""))+
          labs(title =feature.internal.standard$feature_id[i] )->p
        p.list[[i]] <- p
      }else{
        plot_xcms_feature_intensity(object@xcmsData$negativeMS1,
                                    sub(feature.internal.standard$feature_id[i],pattern = "_neg",replacement = ""))+
          labs(title =feature.internal.standard$feature_id[i] )->p
        p.list[[i]] <- p

      }


    }

    (ggplot()+theme_void())/p.list+plot_layout(guides = 'collect')->p.all
    export::graph2pdf(p.all,file = paste0(object@projectInfo[["msDataDir"]],"/dataProcessing/Internal_Standard_Intensity.pdf"),
                      width =5,height = 2*length(p.list))


    }

  object@statData$featureRaw <-feature
  return(object)




}







#' @title Plot MS/MS spectrum for a feature
#' @description Plot experimental and reference MS/MS spectra for a given feature, with annotation details.
#' @param MSdev.obj MSdev object
#' @param feature.id Character string specifying the feature ID
#' @return ggplot object (or NULL if no spectra)
#' @export
#'

plot_MSdev_feature_spectrum <- function(MSdev.obj,feature.id  ){

  feature.data <- MSdev.obj@statData$featureRaw%>%
    dplyr::filter(feature_id == feature.id)
  feature.annotation <- MSdev.obj@annotation[[paste0(feature.data$ion_mode,"Annotation")]][[which(
    rownames(MSdev.obj@xcmsData[[paste0(feature.data$ion_mode,"Feature")]])==
      gsub(pattern = "_|pos|neg",x = feature.id,replacement = ""))]]

  sp.exp <-feature.annotation$expSpec
  sp.ref <-feature.annotation$refSpec

  if (is.null(sp.exp)) {
    return()
  }else{
    sp.exp <- sp.exp%>%
      #Spectra::combineSpectra(
      #  peaks = "intersect",minProp = 0.3,ppm = 50,
      #  intensityFun = median,mzFun = median)%>%
      normalizeSpectra()%>%
      Spectra::filterIntensity(intensity = c(0.05,Inf))%>%
      Spectra::applyProcessing()

    if (is.null(sp.ref)) {

      sp.exp <- sp.exp%>%
        Spectra::combineSpectra(
          peaks = "intersect",minProp = 0.3,ppm = 50,
          intensityFun = median,mzFun = median)

      text.to.show <- paste0("Feature ID :",feature.id,"\n",
                             "Exp Precursor mz :",sp.exp$precursorMz,"\n",
                             "Exp Retention time :",sp.exp$rtime,"\n",
                             "Score: ",feature.data$score,"\n",
                             "Compound: ",feature.data$name,"\n",
                             "Adduct: ",feature.data$adduct,"\n",
                             "Ref Precursor mz: ",sp.ref$precursorMz,"\n",
                             "Ref Retention time: ",sp.ref$rtime,"\n",
                             "INCHIKEY: ",sp.ref$inchikey,"\n",
                             #"KEGG ID: ",sp.ref$kegg.id,"\n",
                             "Reference Source: ",sp.ref$database

      )
      Spectra::plotSpectra(sp.exp,labels = function(z) {
        lbls <- round(mz(z)[[1L]], digits = 4)
        lbls[intensity(z)[[1L]] <= 15] <- ""
        lbls},
        main = text.to.show,
        adj = 0,
        cex.main = 1.5,
        cex.axis = 1,
        labelCex = 1)

    }else{

      sp.score <- Spectra::compareSpectra(sp.exp,sp.ref)
      sp.exp <- sp.exp[which.max(sp.score)]
      text.to.show <- paste0("Feature ID :",feature.id,"\n",
                             "Exp Precursor mz :",feature.data$mz,"\n",
                             "Exp Retention time :",feature.data$rt,"\n",
                             "Score: ",feature.data$score,"\n",
                             "Compound: ",feature.data$name,"\n",
                             "Adduct: ",feature.data$adduct,"\n",
                             "Ref Precursor mz: ",sp.ref$precursorMz,"\n",
                             "Ref Retention time: ",sp.ref$rtime,"\n",
                             "INCHIKEY: ",sp.ref$inchikey,"\n",
                             # "KEGG ID: ",sp.ref$kegg.id,"\n",
                             "Reference Source: ",sp.ref$database

      )
      Spectra::plotSpectraMirror(sp.exp,sp.ref,
                                 ylab = "relative intensity",
                                 labels = function(z) {
                                   lbls <- round(mz(z)[[1L]], digits = 4)
                                   lbls[intensity(z)[[1L]] <= 15] <- ""
                                   lbls},
                                 tolerance = 0.2,
                                 main = text.to.show,
                                 adj = 0,
                                 cex.main = 1.5,
                                 cex.axis = 1,
                                 labelCex = 1

      )


    }


  }



}


#' @title Export MS/MS spectrum and chromatogram for a feature
#' @description Export PNG images of the MS/MS spectrum and chromatogram for a given feature.
#' @param MSdev.obj MSdev object
#' @param feature_id Character string specifying the feature ID
#' @param out.dir Output directory path
#' @return NULL (writes files to disk)
#' @export
#'

export_MSdev_feature_MSMS <- function(MSdev.obj,feature_id,out.dir ){




  png(paste0(out.dir,"/",feature_id,".MSMS.png"),
      res = 100,width = 1000,height = 800)
  par(mar = c(2,2,17,2))
  plot_MSdev_feature_spectrum(MSdev.obj,feature_id)
  dev.off()


  is.pos <-grepl(pattern = "pos",
                 x = feature_id)
  if (is.pos) {
    xcms.xcms <- MSdev.obj@xcmsData$positiveMS1
  }else{
    xcms.xcms <- MSdev.obj@xcmsData$negativeMS1}
  gp <- plot_xcms_feature_chromatogram(xcms.xcms,
                                       feature.id = gsub(x = feature_id,pattern = "[_|pos|neg]",replacement = ""))+
    theme(legend.position = "none")
  export::graph2png(gp , file = paste0(out.dir,"/",feature_id,".Chrom.png"))

}



#' @title Generate sample information table from raw data files
#' @description Read raw MS data files from a directory and generate a sample information data frame.
#' @param raw.data.dir Path to directory containing raw data files
#' @param rawDataFormat File extension of raw data files (default ".raw")
#' @param verbose Logical indicating whether to print messages
#' @return data.frame with sample information
#' @export
#'

get_MS_sampleinfo <- function(raw.data.dir,
                              rawDataFormat=".raw",
                              verbose=T){


  raw.files <- dir(path = raw.data.dir,
                   pattern = paste0(rawDataFormat,"$"),
                   full.names = T,recursive = T)
  raw.files <- normalizePath(raw.files,winslash = "/")
  if (length(raw.files)==0) {
    stop("No ",rawDataFormat," files exist")
  }

  ### generate sampleInfo
  {
    .select_char <- function(char_vector){
      if (all(is.na(char_vector))) {
        return(NA)
      }
      max(char_vector,na.rm = T)
    }
    sample.info <- data.frame(raw.files = raw.files)%>%
      #dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
      dplyr::mutate(ms.labels = gsub(pattern = paste0(rawDataFormat,"$"),
                                     x = basename(raw.files),
                                     ignore.case = T,
                                     replacement = "") ,
                    ms.name = gsub(pattern = "[^0-9A-z]",
                                           x = ms.labels ,
                                           replacement = "_"),
                    polarity = case_when(grepl("pos",x= ms.name, ignore.case =T)~"1",
                                         grepl("neg",x= ms.name, ignore.case =T)~"0",
                                         T~"-1"))%>%
      dplyr::group_by(ms.name)%>%
      dplyr::mutate(raw.files = .select_char(raw.files))%>%
      dplyr::ungroup()%>%
      dplyr::distinct(ms.name,.keep_all = T)%>%
      dplyr::mutate(analysis.time = as.character(file.info(raw.files)$mtime))%>%
      dplyr::mutate(sample.labels = gsub(pattern = paste0("pos|neg|"),
                                         x = ms.labels,
                                         ignore.case = T,
                                         replacement = ""),
                    sample.labels = gsub(pattern = paste0("_$"),
                                         x = sample.labels,
                                         ignore.case = T,
                                         replacement = ""),
                    sample.name =  gsub(pattern = "[^0-9A-z]",
                                        x = sample.labels ,
                                        replacement = "_"),
                    sample.type = case_when(grepl(pattern = "QC",x = sample.labels,ignore.case = T)~ "QC",
                                            grepl(pattern = "blank|blk",x = sample.labels,ignore.case = T)~ "Blank",
                                            T~"Sample"),
                    .before = sample.labels)%>%
      dplyr::mutate(msData.files = case_when(is.na(raw.files)~raw.files,
                                             T~paste0(dirname(raw.files),"/msData/",ms.name)))%>%
      dplyr::arrange(analysis.time)%>%
      dplyr::mutate(no = 1:nrow(.),
                    group = case_when(sample.type=="Sample"~gsub(
                      pattern = "[^A-Za-z]",
                      x= sample.labels,
                      replacement = ""
                    ),
                    T~sample.type),
                    sample.source = group,
                    weight = NA ,
                    xcmsProcessing = "Both",
                    isotope_tracer = NA)%>%
      dplyr::select(no,sample.name,sample.type,sample.labels,sample.source,
                    group, weight,
                    raw.files,
                    polarity,
                    analysis.time,
                    msData.files,
                    ms.name,
                    xcmsProcessing,
                    isotope_tracer)


  }

  if (verbose ) {
    message("Default sample group:")
    show(table(sample.info$group))
  }

  return(sample.info)


}




#' @title Group features across samples using xcms
#' @description Group detected features across samples based on retention time and intensity correlation.
#' @param object MSdev object
#' @param diffRt maximum retention time difference for grouping (seconds)
#' @param intCor minimum intensity correlation threshold
#' @param eicCor minimum EIC correlation threshold
#' @param ... additional arguments passed to xcms grouping functions
#' @return MSdev object with updated feature groups
#' @export
#'
MSdev_xcms_group_features <- function(object,
                                      diffRt = 5,
                                      intCor = 0.5,
                                      eicCor = 0.3,
                                      ...
                                      ){

  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.xcms <- xcms_get_feature_group(xcms.xcms,
                                        diffRt = diffRt,
                                        intCor = intCor,
                                        eicCor = eicCor,
                                        ...)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]


  }

  return(object)



}










#' @title Manually check and edit sample information using Excel
#' @description Open the sample information data frame in Excel for manual editing, then update the MSdev object.
#' @describeIn MSdev_workflow manually check sampleInfo using excel
#' @param object a MSdev object
#' @return MSdev a MSdev object
#' @export
#'

MSdev_checkSampleInfo <- function(object){

  sampleInfo <- object@sampleInfo%>%
    dplyr::group_by(polarity)%>%
    dplyr::mutate(sample.name = str_duplicate_suffix(sample.name))%>%
    dplyr::ungroup()
  sampleInfo <- edit_df_in_excel(sampleInfo,rowname = F)
  ### save
  {
    object@sampleInfo <- sampleInfo
    if(!rlang::is_empty(object@xcmsData)){
      object <- MSdev_update_xcms_pdata(object )}
    object <- .updateProjectInfoFromSampleInfo(object )

  }

  object
}




#' @title Convert raw data files to mzML format
#' @description Convert raw data files to mzML format using MSconvertR, updating sample information.
#' @describeIn MSdev_workflow convert raw files
#' @param object MSdev object
#' @param format.to target format (default "mzML")
#' @return MSdev object with converted files
#' @export

#'

MSdev_msConvert<- function(object,format.to = "mzML"){


  ### filter files
  {
    sample.info <- object@sampleInfo%>%
      dplyr::mutate(
        msData.files = paste0(dirname(raw.files),"/msData/",ms.name,".",format.to),
        raw.exist = file.exists(raw.files),
        ms.exist = file.exists(msData.files))%>%
      dplyr::filter(raw.exist,!ms.exist)
  }

  ### convert
  if (nrow(sample.info)) {
    MSconvertR::msConvert(
      raw.files  = sample.info$raw.files,
      ms.data.names = sample.info$msData.files,
      format.to = format.to,
      BPPARAM = SnowParam(workers =parallel::detectCores()-1,
                                                   progressbar = T))


  }

  ### filter non converted
  {
    sample.info <- object@sampleInfo%>%
      dplyr::mutate(
        msData.files = paste0(dirname(raw.files),"/msData/",ms.name,".",format.to),
        raw.exist = file.exists(raw.files),
        ms.exist = file.exists(msData.files))
    object@sampleInfo <- sample.info[sample.info$ms.exist,]

  }

  object <- MSdev_get_MSinfo(object)
  object <- .updateProjectInfoFromSampleInfo(object)

  ### return
  {
    object@processingInfo$rawDataConvert <- list(
      done = T,
      time = Sys.time(),
      rawFormat =object@projectInfo$rawDataFormat,
      msDataFormat =".mzML"

    )
    MSdev_save(object )
    return(object)

  }
}



#' @title Extract and store MS1 and MS2 spectra from raw data files
#' @description Read raw data files, split spectra by MS level, evaluate noise and purity, and store as on-disk data.
#' @describeIn MSdev_workflow Extract all spectra, split to MS1 and MS2, store as onDiskData
#' @param object MSdev object
#' @param rt.tol retention time tolerance for matching spectra to features (seconds)
#' @param eval.noise logical, whether to evaluate noise in MS2 spectra
#' @param eval.ms1 logical, whether to evaluate purity using MS1 scans
#' @return MSdev object with spectra stored
#' @export
#'
MSdev_extract_Spectra <- function(object,
                                  rt.tol = 10,
                                  eval.noise = F,
                                  eval.ms1 = F){

  sampleInfo <- object@sampleInfo%>%
    dplyr::mutate(msData.files = normalizePath(msData.files))

  sp.list <- list()
  if (nrow(sampleInfo)==0) {
    sp <- Spectra::Spectra()
  } else {
    sp <- Spectra::Spectra(na.omit(sampleInfo$msData.files),
                           backend = Spectra::MsBackendMemory())
    if(1 %in% msLevel(sp)){
      sp.ms1 <- filterMsLevel(sp,1)
      sp.ms1$sp_id <- paste0("MS1_SP",num2str(1:length(sp.ms1)))
      Spectra::spectraNames(sp.ms1) <- sp.ms1$sp_id
    }
    if (2 %in% msLevel(sp)) {
      sp.ms2 <- filterMsLevel(sp,2)
      sp.ms2$sp_id <- paste0("MS2_SP",num2str(1:length(sp.ms2)))
      sp.ms2$precursorMz <- sp.ms2$isolationWindowTargetMz
      if(eval.noise) sp.ms2 <- Spectra_get_noise(sp.ms2)
      if(eval.ms1) sp.ms2 <- Spectra_get_purity(sp.ms2,msLevel = 1,sp.ms1 = sp.ms1)
      Spectra::spectraNames(sp.ms2) <- sp.ms2$sp_id
    }


  }

  ### iso-labeled ms2
  {
    if ("isotope_tracer"%in% colnames(sampleInfo)) {

      sp.ms2$isotope_tracer <- sampleInfo$isotope_tracer[match(sampleNames(sp.ms2),
                                     basename(sampleInfo$msData.files))]
      sp.ms2$from_isotope_tracer <- !is.na(sp.ms2$isotope_tracer)

    }

  }


  ### save on disk
  {
    if(1 %in% msLevel(sp)){
      sp.ms1 <- onDiskData(sp.ms1,
                           path = paste0(object@projectInfo$projectDir,"/MS1_Spectra.rds"))
      object@spectra$MS1_Spectra <- sp.ms1
    }

    if(2 %in% msLevel(sp)){
      sp.ms2 <- onDiskData(sp.ms2,
                           path = paste0(object@projectInfo$projectDir,"/MS2_Spectra.rds"))
      object@spectra$MS2_Spectra <- sp.ms2
    }

  }


  ### assign to feature
  {

    object <-  MSdev_match_Spectra_to_feature(object,rt.tol = rt.tol)
  }


  return(object)

}


#' @title Assign MS2 spectra to features based on precursor m/z and retention time
#' @description Match MS2 spectra to features using precursor m/z and retention time tolerances, updating the MSdev object.
#' @describeIn MSdev_workflow assign Spectra to feature
#' @param object MSdev object
#' @param rt.tol retention time tolerance (seconds)
#' @param ppm m/z tolerance in parts per million
#' @return MSdev object with updated feature-spectra assignments
#' @export
#'
MSdev_match_Spectra_to_feature <- function(object,
                                           rt.tol = 10,
                                           ppm = 20){



  MS2_Spectra <- onDiskData_retrieve(object@spectra$MS2_Spectra)
  MS2_Spectra$feature_id<-NA
  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    sp.ms2 <- MS2_Spectra%>%
      ProtGenerics::filterPolarity(i)
    if (length(sp.ms2)==0) next
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
      as.data.frame()
    sp.ms2.data <- get_Spectra_ms2_feature_id(sp.ms2,
                                              xcms.fdf,
                                              ppm = ppm,
                                              rt.tol = rt.tol)


    ### update MS2_Spectra
    sp.ms2.total <-MS2_Spectra %>%
      Spectra::spectraData()%>%
      as.data.frame()%>%
      dplyr::mutate(feature_id= case_when(
        polarity==i ~ sp.ms2.data[sp_id,]$feature_id,
        T~feature_id
      ))
    Spectra::spectraData(MS2_Spectra ) <- S4Vectors::DataFrame(sp.ms2.total)

    ### update xcms featuredef
    xcms.fdf$ms2_id <- sapply(xcms.fdf$feature_id,
                              function(i){
                                sp_id <- sp.ms2.data$sp_id[which(sp.ms2.data$feature_id==i)]
                                return(sp_id)
                              }
    )
    xcms::featureDefinitions(xcms.xcms) <- xcms.fdf%>%
      S4Vectors::DataFrame()
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]


  }
  object@spectra$MS2_Spectra <- onDiskData(MS2_Spectra,
                                           path = object@spectra$MS2_Spectra@path)

  return(object)


}

#' @title Annotate features using a compound database
#' @description Perform feature annotation using a CompoundDb database, including MS1 candidate search, MS2 scoring, and isotope pattern scoring.
#' @describeIn MSdev_workflow annotation
#' @param object MSdev object
#' @param cpdb_path path to CompoundDb SQLite database
#' @param calc_isopattern_score logical, whether to calculate isotope pattern scores
#' @param ppm m/z tolerance in parts per million
#' @param ... additional arguments passed to annotation functions
#' @return MSdev object with annotation results
#' @export
#'

MSdev_annotation <- function(object,
                             cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite",
                             calc_isopattern_score = F,
                             ppm = 10,
                             ...){

  cpdb <- CompoundDb::CompDb(cpdb_path)

  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    if (is.null(xcms.xcms)) next
    message_with_time("Find MS1 candidate...")
    xcms.xcms <- xcms_get_feature_ms1_candidate(xcms.xcms,
                                                cpdb,
                                                ...)
    message_with_time("Calculate MS2 score...")
    sp.ms2 <- onDiskData_retrieve(object@spectra$MS2_Spectra)
    xcms.xcms <- xcms_get_feature_ms2_score(xcms.xcms ,
                                            cpdb = cpdb,
                                            sp.ms2 = sp.ms2,
                                            ppm = ppm,
                                            ...)
    message_with_time("Calculate isopattern score...")
    xcms.xcms <- xcms_get_feature_isopattern_score(xcms.xcms,
                                                   ppm = 10,
                                                   calc_isopattern_score = calc_isopattern_score,
                                                   BPPARAM = SerialParam(
                                                     #workers = 6,
                                                     progressbar = T
                                                   ))
    xcms.xcms <- xcms_get_feature_annotation(xcms.xcms,
                                             cpdb,
                                             ...)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]

  }

  object@projectInfo$CompoundDB_path <- cpdb_path
  return(object)


}


#'
#' @title Extract and format statistical data from processed MS features
#' @description Extract feature data from xcms, retrieve compound information, filter based on scores, and generate SummarizedExperiment objects.
#' @describeIn MSdev_workflow
#' @param object MSdev object
#' @param keys character vector of compound database keys to retrieve (e.g., "name", "formula")
#' @param score_thresh minimum annotation score threshold
#' @param rt_bin retention time binning width (seconds)
#' @param polarity_paired logical, whether to pair positive and negative polarity features
#' @param candi logical, whether to include all candidates in output
#' @param metabolite logical, whether to filter for metabolite features only
#' @param ... additional arguments
#' @return MSdev object with statData populated
#' @export
#'
MSdev_get_Stat <- function(object,
                           keys = c("name","formula",
                                    "kegg_id",
                                    "inchikey","lipidclass"),
                           score_thresh = 0.5,rt_bin = NA,
                           polarity_paired = T,
                           candi = F,
                           metabolite = T,...){

  ### make se
  {

    sample.info <- object@sampleInfo%>%
      dplyr::filter(polarity_paired|(!polarity_paired)
                    )
    col.order <- sample.info%>%
      dplyr::distinct(sample.name)%>%
      dplyr::pull(sample.name)
    se <- list()
    for (i in 0:1) {
      pol <- ifelse(i==0,"Negative","Positive")
      xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
      if (is.null(xcms.xcms)) {
        se[[pol]] <- SummarizedExperiment::SummarizedExperiment()
        next
      }
      pol.se <- get_xcms_feature_se(xcms.xcms,...)
      se[[pol]] <- pol.se[,intersect(col.order,colnames(pol.se))]
      se[[pol]]$sampleNames<- NULL
      se[[pol]]$no<- NULL
      se[[pol]]$raw.files<- NULL
      se[[pol]]$polarity<- NULL
      se[[pol]]$analysis.time<- NULL
      se[[pol]]$msData.files<- NULL
      se[[pol]]$ms.name<- NULL
      se[[pol]]$files<- NULL
      se[[pol]]$ExpTime<- NULL
    }
    overlap <- do.call("intersect",unname(lapply(se,colnames)))
    feature.se <- do.call("rbind",sapply(se,`[`,,overlap))
    #feature.se <- se[[2]]
  }

  ### formate
  {



    ### sort colname
    rda <- rowData(feature.se)%>%
      as.data.frame()%>%
      dplyr::select(any_of(
        c("feature_id","mzmed","rtmed","compound_id", "adduct","mz_ref","rt_ref","polarity",
          "pave_seed",  "pave_CN",   "pave_cor",
          "score","qc_rsd","sample_rsd","peakMaxo",#ms2_id,
                             "candidate.id","candidate.adduct","candidate.mz","score.ms2")))

    ### retrieve data
    if (!is.null(rda$compound_id)) {
      cpdb <- CompoundDb::CompDb(object@projectInfo$CompoundDB_path)
      db.info <- get_CompDb_info(compound_id = rda$compound_id,
                                 keys = keys,
                                 cpdb = cpdb)
      rda <- rda%>%
        dplyr::mutate(db.info,.after = rtmed,
                      #KEGG_get_cp_linked_gene(kegg_id)
        )
    }


    rowData(feature.se) <- rda


    ### adjust
    feature.se <- se_adjuset_by_weight(feature.se)
    #feature.se <- DEP_impute_mean(feature.se)



    object@statData$feature.se <- feature.se

  }


  ### all candidate
  if (candi) {
    {
      candi.rda <- rda%>%
        dplyr::mutate(candidate.n = sapply(candidate.id,length))
      candi.rda.split <- candi.rda[rep(candi.rda$feature_id,candi.rda$candidate.n),]%>%
        dplyr::group_by(feature_id)%>%
        dplyr::mutate(temp_id = 1:n())%>%
        dplyr::rowwise()%>%
        dplyr::mutate(compound_id = candidate.id[[temp_id]],
                      adduct = candidate.adduct[[temp_id]],
                      mz_ref = candidate.mz[[temp_id]],
                      score = score.ms2[[temp_id]])%>%
        dplyr::ungroup()
      cpdb <- CompoundDb::CompDb(object@projectInfo$CompoundDB_path)
      db.info <- get_CompDb_info(compound_id = candi.rda.split$compound_id,
                                 keys = keys,
                                 cpdb = cpdb)
      candi.rda.split <- candi.rda.split%>%
        dplyr::mutate(db.info,.after = rtmed)
      candi.se <- feature.se[candi.rda.split$feature_id,]
      rowData(candi.se) <- candi.rda.split

    }

    object@statData$candidate.se <- candi.se
  }



  ###  metabolite
  if(metabolite){
    ### select unique feature
    .uniqueFeatures <- function(score,intensity){
      score <- ifelse(score > 0.75 , 10,1)
      unique.score <- score*log10(intensity)
      unique.score
    }

    if (is.null(rda$score)) {
      rda$score <- 0
    }

    rda.filter <- rda%>%
      as.data.frame()%>%
      dplyr::filter(score >= score_thresh,
                    !is.na(compound_id))%>%
      dplyr::mutate( rt_bins = ceiling((rtmed)/rt_bin ),
                     temp_id = paste0( compound_id, "_", rt_bins )
                     )%>%
      dplyr::group_by(temp_id)%>%
      dplyr::slice_max(.uniqueFeatures(score,peakMaxo))%>%
      ungroup()


    metabolite.se <- feature.se[rownames(feature.se)%in%rda.filter$feature_id,]
    object@statData$metabolite.se <- metabolite.se

  }



  object

}









MSdev_update_xcms_pdata <- function(object,
                                    XCMSnExp = T,
                                    XChromatograms = F){


  sample_info <- object@sampleInfo
  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    ###XCMSnExp
    if(XCMSnExp){
      xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
      if (is.null(xcms.xcms)) next
      xcms.pdata <- Biobase::pData(xcms.xcms)%>%
        dplyr::mutate(sample_info[match(msData.files,
                                        sample_info$msData.files),  ])
      xcms.pdata -> Biobase::pData(xcms.xcms)
      xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
    }

    ###XChromatograms
    if(XChromatograms){
      xcms.chrom <- object@xcmsData[[paste0(pol,"_Chromatograms")]]
      if (is.null(xcms.chrom)) next
      xcms.chrom.data <- onDiskData_retrieve(xcms.chrom)
      xcms.pdata <- Biobase::pData(xcms.chrom.data)%>%
        dplyr::mutate(sample_info[match(msData.files,
                                        sample_info$msData.files),  ])
      xcms.pdata -> Biobase::pData(xcms.chrom.data)
      xcms.chrom <- onDiskData_update(xcms.chrom,xcms.chrom.data )
      xcms.chrom -> object@xcmsData[[paste0(pol,"_Chromatograms")]]
    }
  }

  return(object)


}






#' @title Retrieve MS2 spectra from MSdev object
#' @description Retrieve the stored MS2 Spectra object from the MSdev object.
#' @param object MSdev object
#' @return Spectra object
#' @export
#'
get_MSdev_ms2_Spectra <- function(object){

  sp <- object@spectra$MS2_Spectra%>%
    onDiskData_retrieve()
  #sp.total <- do.call(`c`,unname(object@spectra))
  return(sp)
}


#' @title Extract chromatograms for features
#' @description Extract chromatograms for specified features from xcms data, storing them as on-disk data.
#' @param object MSdev object
#' @param BPPARAM BiocParallel backend for parallel processing
#' @param feature.list optional list of feature IDs with names "Positive" and "Negative"
#' @return MSdev object with chromatograms stored
#' @export
#'
MSdev_get_feature_chrom <- function(object,BPPARAM =  SnowParam(
  workers  = parallel::detectCores() -1,
  progressbar = T),feature.list = NULL){

  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    if (is.null(feature.list[[pol]])) {
      fid = featureDefinitions(xcms.xcms)$feature_id
    }else{
      fid = feature.list[[pol]]
    }
    message("Extract ", length(fid) , " features ",pol)
    #xcms.chrom <- get_xcms_feature_chrom(xcms.xcms,
    #                                     feature.id = fid,
    #                                    sample = "all",rt = "all",
    #                                    BPPARAM = SnowParam(progressbar = T))
    xcms.chrom <- featureChromatograms(xcms.xcms,
                                       features = fid,
                                       expandRt = Inf,
                                       filled = T,
                                       BPPARAM = BPPARAM
                                       )
    xcms.chrom <- onDiskData(xcms.chrom,
                             path = paste0(object@projectInfo$projectDir,"/",pol,"_Chromatograms.rds"))
    #featureValues( xcms.chrom ,value = "maxo")
    #rownames(xcms.chrom) <- fid
    object@xcmsData[[paste0(pol,"_Chromatograms")]] <- xcms.chrom
  }

  return(object)

}



get_MSdev_isotopologues_data <- function(object){


 data.list <- list()

  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]

    xcms.fdf.iso <- featureDefinitions(xcms.xcms)%>%
      as.data.frame()%>%
      dplyr::mutate(is_seed = feature_id %in% iso_seed)%>%
      dplyr::filter(!is.na(iso_seed))

    xcms.fdf.iso[!xcms.fdf.iso$is_seed,
                 c("compound_id","name","adduct","score",
                   "mz_ref","rt_ref","formula","smiles")] <- NA

    xcms.fdf.iso <- xcms.fdf.iso%>%
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
        smiles = na.unique(smiles))%>%
      dplyr::filter(!is.na(compound_id))

    natural_prob <- mapply(iso_count = xcms.fdf.iso$iso_count,
           C_count = xcms.fdf.iso$formula_max_iso_count,
           FUN = iso_prob)
    xcms.fdf.iso$natural_prob <- natural_prob

    iso.fraction <- get_xcms_iso_fraction(xcms.xcms)
    iso.fraction <- iso.fraction[xcms.fdf.iso$feature_id,]

    val <- featureValues(xcms.xcms)[xcms.fdf.iso$feature_id,]
    colnames(val) <- paste0("int_",pData(xcms.xcms)$sample.name)

    frac <- iso.fraction
    colnames(frac) <- paste0("fraction_",pData(xcms.xcms)$sample.name)

    adj <- iso.fraction - natural_prob
    adj[adj<0] <- 0
    colnames(adj) <- paste0("fraction_adjusted_",pData(xcms.xcms)$sample.name)

    data <- cbind(xcms.fdf.iso,val,frac,adj)
    data.list[[pol]] <- data

  }

  do.call(rbind,data.list)



}


count_MSdev_peaks <- function(object){

  xcms <- object@xcmsData$NegativeMS1
  pks <- table(chromPeaks(xcms)[,"sample"])
  names(pks) <- pData(xcms)[,1]
  print(pks)


  xcms <- object@xcmsData$PositiveMS1
  pks <- table(chromPeaks(xcms)[,"sample"])
  names(pks) <- pData(xcms)[,1]
  print(pks)
}


get_MSdev_instrument <- function(object){

  manufacturer <- object@projectInfo$msManufacturer
  manufacturer <- ifelse(manufacturer == "instrument model","Unknow manufacturer",manufacturer)
  model <- object@projectInfo$msModel
  model <- ifelse(model == "Applied Biosystems instrument model",
                         "Unknow model",model)

  paste(manufacturer,model)


}


MSdev_get_feature_wmean <- function(object){

  for (i.pol in 0:1) {
    pol <- ifelse(i.pol==0,"Negative","Positive")
    polarity.tag <- paste0(pol,"MS1")
    xcms.xcms <- object@xcmsData[[polarity.tag]]
    xcms.xcms <- xcms_get_feature_wmean(xcms.xcms)
    xcms.xcms -> object@xcmsData[[polarity.tag]]
  }

  return(object)
}
