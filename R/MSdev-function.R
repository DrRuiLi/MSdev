

saveMSdev <- function(object){
  MSdev <- object
  save.dir <- dirname(object@projectInfo$MSdevFile)
  if (!dir.exists(save.dir)) {
    dir.create(save.dir,recursive = T)
  }
  save(MSdev, file =  object@projectInfo$MSdevFile)
  invisible(MSdev)
}

#' @title readInRawData
#' @description read in ms raw data from `object@projectInfo$msDataDir`
#' and generate a table `sampleInfo`
#' note this function read in file according to their file names, those contaion both "pos" and "neg"
#' will be regarded as two files
#' @param object a `MSdev` object
#' @details discriminate sample type and ion mode according char in file names
#'
#' grep "pos" and "neg" for ion mode
#'
#' grep "blank", "blk" and "QC" for sample type, other samples will regard as "Sample"
#' @return a `MSdev` object
#' @export
#'
#' @examples
readInRawData <- function(object){

  projectInfo <- object@projectInfo
  msData.dir <- projectInfo$msDataDir
  raw.files <- dir(path = projectInfo$rawDataDir,
                   pattern = paste0(projectInfo$rawDataFormat,"$"),
                   full.names = T)

  ### generate sampleInfo
  {
    .select_char <- function(char_vector){
      if (all(is.na(char_vector))) {
        return(NA)
      }
      max(char_vector,na.rm = T)
    }
    sample.info <- data.frame(raw.files = raw.files)%>%
      dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
      dplyr::mutate(raw.file.positive = case_when(grepl(pattern = "pos", x = raw.files,ignore.case = T)~raw.files),
                    raw.file.negative = case_when(grepl(pattern = "neg", x = raw.files,ignore.case = T)~raw.files))%>%
      dplyr::mutate(sample.abbreviation= gsub(pattern = paste0("pos|neg|",projectInfo$rawDataFormat,"$"),
                                              x = basename(raw.files) ,
                                              ignore.case = T,
                                              replacement = ""),
                    sample.abbreviation = tolower(sample.abbreviation))%>%
      dplyr::group_by(sample.abbreviation)%>%
      dplyr::mutate(raw.file.positive = .select_char (raw.file.positive),
                    raw.file.negative = .select_char (raw.file.negative))%>%
      dplyr::ungroup()%>%
      dplyr::distinct(sample.abbreviation,.keep_all = T)%>%
      dplyr::mutate(analysis.time.positive =as.character( file.info(raw.file.positive)$mtime),
                    analysis.time.negative = as.character(file.info(raw.file.negative)$mtime))%>%
      dplyr::mutate(sample.type = case_when(grepl(pattern = "QC",x = sample.abbreviation,ignore.case = T)~ "QC",
                                            grepl(pattern = "blank|blk",x = sample.abbreviation,ignore.case = T)~ "Blank",
                                            T~"Sample"),
                    .before = sample.abbreviation)%>%
      dplyr::mutate(sample.name = paste0(sample.type,str_pad(1:nrow(.), ceiling(log10(nrow(.))),pad = "0")),
                    .before = sample.type)%>%
      dplyr::mutate(msData.file.positive = case_when(is.na(raw.file.positive)~raw.file.positive,
                                                     T~paste0(msData.dir,"/pos/",sample.name,".mzML")),
                    msData.file.negative = case_when(is.na(raw.file.negative)~raw.file.negative,
                                                     T~paste0(msData.dir,"/neg/",sample.name,".mzML")))%>%
      dplyr::arrange(analysis.time.positive)%>%
      dplyr::mutate(no = 1:nrow(.),
                    label = sample.abbreviation,
                    group = sample.type,
                    weight = NA ,
                    xcmsProcessing = "Both")%>%
      dplyr::select(no,sample.name,sample.type,group ,label, weight,
                    sample.abbreviation,
                    raw.file.positive,raw.file.negative,
                    analysis.time.positive,analysis.time.negative,
                    msData.file.positive,msData.file.negative,
                    xcmsProcessing)

  }
  ### save
  {
    object@sampleInfo <- sample.info
    object <- .updateProjectInfoFromSampleInfo(object )
    object@processingInfo$readInRawData$done <- T


  }
  object


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



get_MSdev_MSinfo <- function(object){


  ### define acquire Type
  ### note, these model string are identified by mzR
  {
    HRMS <- c("Q Exactive Plus","TripleTOF 6600","Orbitrap Exploris 480")
    TQMS <- c("TSQ Quantis")
    model.df <- data.frame(
      model = c(HRMS,TQMS),
      type = c( rep("HRMS",length(HRMS)),
                rep("TQMS",length(TQMS)))
    )

  }



  ### update SampleInfo
  if (!"manufacturer"%in%colnames(object@sampleInfo) ) {
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


get_MSdev_newSamples <- function(object,
                                 raw.data.dir = object@projectInfo$rawDataDir){
  sample.info <- object@sampleInfo
  sample.info.new <- get_MS_sampleinfo(raw.data.dir,
                                       rawDataFormat = object@projectInfo$rawDataFormat,F)%>%
    dplyr::filter(!raw.files%in% sample.info$raw.files)
  message("Get ",nrow(sample.info.new)," samples")

  object@sampleInfo <- sample.info%>%
    bind_rows(sample.info.new)
  object <- msConvert_MSdev(object)
  show(object)
  object

}


MSdev_xcmsProcessing <- function(object,...){

  MS.mode <- object@projectInfo$msAcquisition

  if ("FS" %in% MS.mode|"DDA" %in% MS.mode) {
    object <- xcmsProcessingMSdev.DDA(object,...)
    return(object)
  }

  if ("MRM" %in% MS.mode) {
    object <-  xcmsProcessingMSdev.MRM(object)
    return(object)
  }





}

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
    xcms.xcms <-
      xcmsProcessingMS1(msDataFiles = sample.info.polarity$msData.files,
                        ion_mode = i,
                        peaksGroup =sample.info.polarity$sample.type,
                        findChromPeaks_param  = xcms.param$findChromPeaks,
                        ...
    )

    Biobase::pData(xcms.xcms ) <-
      cbind(Biobase::pData(xcms.xcms),
            sample.info.polarity)

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


  if (object@experimentInfo@General$Name == "Metabolomics_YLF") {

    fpp <- MassifquantParam(ppm = 5,
                            peakwidth = c(10,100),
                            snthresh = 10,
                            prefilter = c(5,1000),
                            verboseColumns=T,
                            withWave = T)

    msdev.param <- list(findChromPeaks = fpp)
    return(msdev.param)
  }

  msdev.param <- get_MSdev_xcms_param_by_exp(object)
  return(msdev.param)




}


get_MSdev_xcms_param_by_exp <- function(object){



  MS.mode <- object@projectInfo$msAcquisition
  MS.instru <-object@projectInfo$msModel
  MS.LC.rate <- object@experimentInfo@Chroma_gradient[[1]]$Flow_rate%>%mean
  MS.LC.time <- object@experimentInfo@Chroma_gradient[[1]]$time%>%max
  cwp <- CentWaveParam(fitgauss = F,verboseColumns = T)

  ### ppm
  cwp@ppm <- switch(MS.instru,
                    "Q Exactive Plus" = 5,
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
                          "Q Exactive Plus" = c(5,5000),
                          "SCIEX TripleTOF 6600" = c(3,100),
                          "Thermo Quantis" = c(3,100),
                          c(3,100))



  msdev.param <- list(findChromPeaks = cwp)


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
#' @examples
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
#' @examples
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
#' @examples
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


#' @title getStaDataMSdev
#'
#' @param object
#' @param MSDB.keys
#'
#' @return
#' @export
#'
#' @examples
getStaDataMSdev <- function(object,missing = NA,
                            MSDB.keys =c("name","adduct","formula","inchikey" ,"database_origin")
){

  sampleInfo <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing %in% c("Both","MS1"))



  annotationPos <- lapply(object@annotation$positiveAnnotation, function(x){
    x[c( "mz"    ,
         "rt" ,
         "ref.mz" ,
         "ref.rt" ,
         "score" ,
         "MSDB_id" )]
  })%>%data.table::rbindlist()
  annotationNeg <- lapply(object@annotation$negativeAnnotation, function(x){
    x[c( "mz"    ,
         "rt" ,
         "ref.mz" ,
         "ref.rt" ,
         "score" ,
         "MSDB_id" )]
  })%>%data.table::rbindlist()

  featurePos <- get_features_from_xcms(object@xcmsData$positiveMS1,missing = missing)
  featurePos <- SummarizedExperiment::rowData(featurePos)%>%
    as.data.frame()%>%
    rownames_to_column("feature_id")%>%
    dplyr::select(feature_id,qc_rsd,sample_rsd,med_intensity)%>%
    dplyr::mutate(feature_id = paste0(feature_id , "_pos"),
                  ion_mode = "positive")%>%
    cbind(annotationPos,SummarizedExperiment::assay(featurePos))%>%
    dplyr::rename_with( ~sub(pattern = ".mzML",replacement = "",x = .x))%>%
    remove_rownames()

  featureNeg <- get_features_from_xcms(object@xcmsData$negativeMS1,missing = missing)
  featureNeg <- SummarizedExperiment::rowData(featureNeg)%>%
    as.data.frame()%>%
    rownames_to_column("feature_id")%>%
    dplyr::select(feature_id,qc_rsd,sample_rsd,med_intensity)%>%
    dplyr::mutate(feature_id = paste0(feature_id , "_neg"),
                  ion_mode = "negative")%>%
    cbind(annotationNeg,SummarizedExperiment::assay(featureNeg))%>%
    dplyr::rename_with( ~sub(pattern = ".mzML",replacement = "",x = .x))%>%
    remove_rownames()

  featureAll <-  dplyr::bind_rows(featurePos,featureNeg)
  featureAllanno <- MSdb:::getInfoFromMSDB(featureAll$MSDB_id,
                                           msdb_path = object@projectInfo$MSDB_path,
                                           keys =  MSDB.keys)
  featureAll<- add_column(featureAll,featureAllanno[,-1],.after = "feature_id")
  object@statData$featureRaw <-featureAll

  ### adjust
  # object <- adjustFeatureByIS(object)
  #  object <- adjustFeatureByGQC(object,to.adjust = "featureRaw")
  object <- adjustFeatureByweight(object,to.adjust = "featureRaw")

  object@statData$feature <- object@statData$feature%>%
    dplyr::filter(qc_rsd <0.3)%>%
    #dplyr::filter(gqc_r2 >0.8)%>%
    dplyr::filter()
  .uniqueFeatures <- function(score,intensity){
    score <- ifelse(score >0.3 , 10,1)
    unique.score <- score*log10(intensity)
    unique.score

  }
  object@statData$metabolites <- object@statData$feature%>%
    dplyr::filter(!is.na(inchikey))%>%
    dplyr::group_by(inchikey)%>%
    dplyr::slice_max(.uniqueFeatures(score,med_intensity ))%>%
    dplyr::ungroup()

  return(object)
}


adjustFeatureByweight <- function(object,to.adjust = "featureRaw"){

  sampleInfoToAdjust <- object@sampleInfo%>%
    dplyr::filter(!is.na(weight))
  weight <- sampleInfoToAdjust$weight / mean(sampleInfoToAdjust$weight )
  featureMatrix <- object@statData[[to.adjust]]%>%
    column_to_rownames("feature_id")%>%
    dplyr::select(sampleInfoToAdjust$sample.name)%>%
    as.matrix()

  featureMatrixAdjusted <-t( t(featureMatrix)/weight)
  object@statData$feature <-object@statData[[to.adjust]]
  featureMatrixAdjusted -> object@statData$feature[,sampleInfoToAdjust$sample.name]
  object

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
#' only [M+H] and [M-H] are considered. Correlation and intensity will be plot based on `object@statData[["featureRaw"]]`, please check.
#' A column "internal_standard" will be added in `object@statData[["featureRaw"]]`
#'
#' @param object
#' @param corr.thred
#'
#' @return
#' @export
#'
#' @examples
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




#' @title getSEMSdev
#' @description extrat SummarizedExperiment::SummarizedExperiment, which combine coldata(sample info) and rowdata (metabolites/feature)
#' and perform data filter, normalization and imputation, refer to package DEP
#'
#' @param MSdev.obj
#'
#' @return
#' @export
#'
#' @examples
getSEMSdev <- function(MSdev.obj){

  col.info <- MSdev.obj@sampleInfo%>%
    dplyr::filter(sample.type == "Sample")%>%
    dplyr::mutate(label = sample.name,
                  condition = group ,
                  replicate = 1)%>%
    dplyr::group_by(condition)%>%
    dplyr::mutate(replicate = 1:n())

  data.unique <- DEP::make_unique(MSdev.obj@statData$metabolites , names ="feature_id" , ids = "feature_id")
  data.colum <- which(colnames(data.unique )%in% col.info$sample.name)
  data.se <- DEP::make_se(proteins_unique = data.unique,
                          columns = data.colum,
                          expdesign = col.info )
  data_filt <- DEP::filter_missval(data.se, thr = min(table(col.info$group))*0.3)
  data_norm <- DEP::normalize_vsn(data_filt)
  data_imp <- DEP::impute(data_norm, fun = "MinProb")

  MSdev.obj@statData$data.se$data.raw$data.raw <- data_imp
  MSdev.obj
}





#' Title
#'
#' @param MSdev.obj
#' @param feature_id
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param MSdev.obj
#' @param feature_id
#' @param out.dir
#'
#' @return
#' @export
#'
#' @examples
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



get_MS_sampleinfo <- function(raw.data.dir,
                              rawDataFormat=".raw",
                              verbose=T){


  raw.files <- dir(path = raw.data.dir,
                   pattern = paste0(rawDataFormat,"$"),
                   full.names = T)

  ### generate sampleInfo
  {
    .select_char <- function(char_vector){
      if (all(is.na(char_vector))) {
        return(NA)
      }
      max(char_vector,na.rm = T)
    }
    sample.info <- data.frame(raw.files = raw.files)%>%
      dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
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
                                             T~paste0(dirname(raw.files),"/mzML/",ms.name,".mzML")))%>%
      dplyr::arrange(analysis.time)%>%
      dplyr::mutate(no = 1:nrow(.),
                    group = case_when(sample.type=="Sample"~gsub(
                      pattern = "[^A-Za-z]",
                      x= sample.labels,
                      replacement = ""
                    ),
                    T~sample.type),
                    weight = NA ,
                    xcmsProcessing = "Both")%>%
      dplyr::select(no,sample.name,sample.type,sample.labels,group, weight,
                    raw.files,
                    polarity,
                    analysis.time,
                    msData.files,
                    ms.name,
                    xcmsProcessing)


  }

  if (verbose ) {
    message("Default sample group:")
    show(table(sample.info$group))
  }

  return(sample.info)


}


get_MSdev_Inclusion_Queue <- function(object){

  for (i in 0:1) {
    polarity <-ifelse(i==0,"Negative","Positive")
    polarity.tag <- paste0(polarity,"MS1")
    xcms.xcms <- object@xcmsData[[polarity.tag]]
    if (is.na(xcms.xcms)) {
      next
    }
    feature.rsd <- get_features_from_xcms(xcms.xcms)@elementMetadata%>%as.data.frame()
    feature.stat <- featureDefinitions_PeakSta(xcms.xcms)%>%
      cbind(feature.rsd[,c("qc_rsd","sample_rsd","med_intensity")])
    dda.mine.queue <- feature.stat%>%
      dplyr::mutate( qc_rsd.score = (log(0.3)-log(qc_rsd)),
                     MS1.score = qc_rsd.score*log10(peakMaxo)  )%>%
      dplyr::arrange(-MS1.score)%>%
      rownames_to_column("feature.id")%>%
      dplyr::mutate(CE10=10,CE20=20,CE30=30,CE40=40,CE50=50)%>%
      pivot_longer(CE10:CE50,names_to = "CE.tag",values_to = "collisionEnergy")%>%
      dplyr::mutate(DDA.id= paste0(feature.id,"_",CE.tag),
                    acquired =F,
                    acquired.in.list = NA,
                    queued.in.list = NA,
                    queued.time = 0)

    object@statData[[paste0("DDA_mine_queue_",polarity)]] <- dda.mine.queue
    object@statData[[paste0("DDA_mine_list_",polarity)]] <- list()

  }


  object

}

get_MSdev_Inclusion_List <- function(object){

  for (i in 0:1) {
    polarity <-ifelse(i==0,"Negative","Positive")
    DDA.queue <- object@statData[[paste0("DDA_mine_queue_",polarity)]]
    if ( is.null(DDA.queue)) {
      next
    }
    DDA.mine.list <-DDA.queue%>%
      dplyr::ungroup()%>%
      dplyr::filter(!acquired)%>%
      dplyr::mutate(ion.cluster = cluster_ion(mzmed,
                                              rtmed,
                                              rt.tol = 60))%>%
      dplyr::group_by(ion.cluster)%>%
      dplyr::slice_max(MS1.score,n=1,with_ties =F)%>%
      dplyr::ungroup()%>%
      dplyr::slice_max(MS1.score,n=5000,with_ties =F)%>%
      dplyr::mutate(feature.id = paste0(feature.id ,"_", CE.tag))

    ### update list
    queue.list <- object@statData[[paste0("DDA.mine.list.",polarity)]]
    if (length(queue.list)) {

      list.name =str_add( max(names(queue.list)),1)
      queue.list <- append(queue.list,
                            list( DDA.mine.list))
      names(queue.list)[length(queue.list)] <-list.name

    }else{
      list.name <- "DDA.mine.list001"
      queue.list <- list("DDA.mine.list001" = DDA.mine.list)

    }
    object@statData[[paste0("DDA.mine.list.",polarity)]] <-queue.list

    DDA.mine.list.qe <- feature_def_to_QE_inclusion(DDA.mine.list,polarity = polarity)
    write.csv(DDA.mine.list.qe,
              file = paste0(object@projectInfo$projectDir,"/",
                            list.name,".csv"))

    ### update DDA.queue
    DDA.queue <- DDA.queue %>%
      dplyr::mutate(queued.in.list = case_when(
        DDA.id %in% DDA.mine.list$feature.id ~ paste0(queued.in.list,";",list.name),
        T~queued.in.list),
        queued.time = case_when(
          DDA.id %in% DDA.mine.list$feature.id ~ queued.time+1,
          T~queued.time))
    DDA.queue -> object@statData[[paste0("DDA.mine.queue.",polarity)]]

  }


  object

}


get_MSdev_MS2acquisitionStat <- function(object){

  assign_ms2_list <- function(pmz,rt,ce ,il){

    il %>%
      dplyr::filter(abs(mzmed-pmz)/pmz < 1e-5,
                    rtmed < rtmax,
                    rtmed > rtmin,
                    collisionEnergy==ce)%>%
      dplyr::pull( DDA.id)->x
    if (length(x)==0) {
      return(NA)

    }
    return(paste0(x,collapse = ";"))

  }

  for (i in 0:1) {

    polarity <-ifelse(i==0,"Negative","Positive")
    DDA.queue <- object@statData[[paste0("DDA.mine.queue.",polarity)]]
    if ( is.null(DDA.queue)) {
      next
    }
    sample.info <- object@sampleInfo%>%
      dplyr::filter(polarity %in% c(i,-1),
                    msLevels %in% c(2))
    xcms.xcms <- readMSData(sample.info$msData.files,mode = "onDisk")
    xcms.scan <- get_xcms_scan_Stat(xcms.xcms)%>%
      dplyr::filter(msLevel==2)%>%
      dplyr::rowwise()%>%
      dplyr::mutate(assigned.id = assign_ms2_list(pmz = precursorMZ,
                                    rt = retentionTime,
                                    ce = collisionEnergy,
                                    il = DDA.queue))%>%
      dplyr::filter(!is.na(assigned.id))%>%
      dplyr::ungroup()
    ms2.stat <- xcms.scan%>%
      dplyr::ungroup()%>%
      dplyr::group_by(fileIdx)%>%
      dplyr::mutate(total.id = paste0(assigned.id,collapse = ";"))%>%
      dplyr::distinct(fileIdx,total.id)%>%
      dplyr::mutate(files = sampleNames(xcms.xcms)[fileIdx])

    ms2.list <-lapply(ms2.stat$total.id,function(x){
              strsplit(x,";")%>%unlist()
            })
    for (i in 1:length(ms2.list)) {
      ids <- ms2.list[[i]]
      a <-  ids%in% DDA.queue$DDA.id
      DDA.queue <- DDA.queue%>%
        dplyr::ungroup()%>%
        dplyr::mutate(acquired = case_when(DDA.id %in% ids~T,
                                           T~acquired),
                      acquired.time = case_when(
                        DDA.id %in% ids~paste0(acquired.time,
                                               ";",sampleNames(xcms.xcms)[i]),
                                                T~acquired.time))%>%
        dplyr::ungroup()%>%
        dplyr::mutate(fail.time = case_when(acquired~0,
                                            T~queued.time))%>%
        dplyr::group_by(feature.id)%>%
        dplyr::mutate(
          acquired = case_when(any(fail.time >3)~T,
                                           T~acquired))

    }

    DDA.queue -> object@statData[[paste0("DDA.mine.queue.",polarity)]]



  }


  return(object)

}

