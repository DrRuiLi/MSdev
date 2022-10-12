
saveMSdev <- function(object){
  MSdev <- object
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
                  msData.file.negative = case_when(is.na(raw.file.negative)~raw.file.positive,
                                                   T~paste0(msData.dir,"/neg/",sample.name,".mzML")))%>%
    dplyr::mutate(group = sample.type,weight = NA ,
                  xcmsProcessing = "Both")%>%
    dplyr::select(sample.name,sample.type,group , weight,
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

#' @title checkSampleInfo
#' @description manually check sampleInfo using excel
#' @param object a `MSdev` object
#'
#' @return a `MSdev` object
#' @export
#'
#' @examples
checkSampleInfo <- function(object){

  sampleInfo <- object@sampleInfo
  sampleInfo <- edit_df_in_excel(sampleInfo)
  ### save
  {
    object@sampleInfo <- sampleInfo
    object <- .updateProjectInfoFromSampleInfo(object )

  }

  object
}

.updateProjectInfoFromSampleInfo <- function(object){

  object@sampleInfo -> sampleInfo
  object@projectInfo$sampleCount <- sum(sampleInfo$sample.type=="Sample")
  object@projectInfo$rawDataFileCount <- sum(!is.na(sampleInfo$raw.file.positive),
                                             !is.na(sampleInfo$raw.file.negative))
  sampleInfo <- object@sampleInfo%>%
    dplyr::mutate(pos = ifelse(is.na(raw.file.positive),NA,"positive"),
                  neg = ifelse(is.na(raw.file.negative),NA,"negative"))%>%
    pivot_longer(c("pos","neg"),names_to = "ion_mode",values_to = "p")
  object@projectInfo$rawDatafiles <- table(sampleInfo$p,sampleInfo$sample.type)

  object
}


msConvert_MSdev <- function(object){

  ### convert
  if (object@projectInfo$rawDataFormat == ".raw") {
    msconvert_raw2mzML(raw.files  = object@sampleInfo$raw.file.positive,
                       mzML.files = object@sampleInfo$msData.file.positive,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

    msconvert_raw2mzML(raw.files  = object@sampleInfo$raw.file.negative,
                       mzML.files = object@sampleInfo$msData.file.negative,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

  }else if(object@projectInfo$rawDataFormat == ".wiff"){
    msconvert_wiff2mzML(wiff.files  = object@sampleInfo$raw.file.positive,
                       mzML.files = object@sampleInfo$msData.file.positive,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

    msconvert_wiff2mzML(wiff.files  = object@sampleInfo$raw.file.negative,
                       mzML.files = object@sampleInfo$msData.file.negative,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

  }

  ### update SampleInfo
  {
    sampleInfo <-object@sampleInfo
    sampleInfo$analysis.time.positive <- getmsExpTime(sampleInfo$msData.file.positive)$ExpTime
    sampleInfo$analysis.time.negative <- getmsExpTime(sampleInfo$msData.file.negative)$ExpTime
    sampleInfo -> object@sampleInfo

  }
  ### return
  {
    object@processingInfo$rawDataConvert <- list(
      done = T,
      time = Sys.time(),
      rawFormat =object@projectInfo$rawDataFormat,
      msDataFormat =".mzML"

    )
    saveMSdev(object )
    return(object)

  }
}

xcmsProcessing_fullscan_DDA <- function(object){



  sampleInfoPos <-dplyr::filter(object@sampleInfo,
                                xcmsProcessing %in% c("MS1","Both")
                                )%>%
    dplyr::filter(!is.na(msData.file.positive))
  object@xcmsData$positiveMS1  <- xcmsProcessingMS1(msDataFiles = sampleInfoPos$msData.file.positive,
                                                       ion_mode = 1,
                                                       peaksGroup =sampleInfoPos$sample.type,
                                                       centWaveParam =xcms::CentWaveParam(ppm = 10,snthresh = 10,
                                                                                          peakwidth = c(5,50),
                                                                                          prefilter = c(3,1000))
  )

  pData(object@xcmsData$positiveMS1) <- cbind(pData(object@xcmsData$positiveMS1),
                                              sampleInfoPos)

  sampleInfoNeg <-dplyr::filter(object@sampleInfo,
                                xcmsProcessing %in% c("MS1","Both")
  )%>%
    dplyr::filter(!is.na(msData.file.negative))
  object@xcmsData$negativeMS1  <- xcmsProcessingMS1(msDataFiles = sampleInfoNeg$msData.file.negative,
                                                       ion_mode = 0,
                                                       peaksGroup = sampleInfoNeg$sample.type,
                                                       centWaveParam =xcms::CentWaveParam(ppm = 10,snthresh = 10,
                                                                                          peakwidth = c(5,50),
                                                                                          prefilter = c(3,1000))
  )

  pData(object@xcmsData$negativeMS1)$sampleType <-  cbind(pData(object@xcmsData$negativeMS1),
                                                          sampleInfoNeg)

  extractFeature(object)






}

extractFeature <- function(object){

  object@xcmsData$positiveFeature <- as.data.frame(featureDefinitions(object@xcmsData$positiveMS1))
  object@xcmsData$negativeFeature <- as.data.frame(featureDefinitions(object@xcmsData$negativeMS1))
  object

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
    spectra.pos <- Spectra::Spectra(sampleInfo$msData.file.positive,backend = MsBackendDataFrame())%>%
    filterMsLevel(2)

    spectra.neg <- Spectra::Spectra(sampleInfo$msData.file.negative,backend = MsBackendDataFrame())%>%
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
                       rt_tol = 20){
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
featureCandidate<- function(object,mz.ppm = 10,
                            spectraDatabase =
                              "C:\\Users\\91879\\OneDrive\\Documents\\Code\\R\\Projecct\\2022.1.17_Compounds.database\\Spectra.integrated.database.integration.2022_02_12.Rdata")
  {
  Spectra_database <- load_as_var(spectraDatabase)
  .matchMz <- function(xcmsFeature,spectraDB){
    feature.mz_rt <- data.frame(mz = xcmsFeature$mzmed,
                                rt = xcmsFeature$rtmed)
    lib.precursormz <- precursorMz(spectraDB)
    lib.rtime <- rtime(spectraDB)
    lib.candidate <- apply(feature.mz_rt,1,function(x){

      mz.hit <- abs( lib.precursormz-x[["mz"]]) < x[["mz"]]*mz.ppm /1e6
      which(mz.hit   )

    })
    #sum(sapply(lib.candidate, length)>0)
    featureCandidate <- lapply(lib.candidate,  function(x){
      sp <- spectraDB[x]
      if (length(sp)== 0) {
        return(NULL)
      }else{
        return(sp)
      }
    })
    return(featureCandidate)

  }
  object@projectInfo$MSDB_path <- spectraDatabase
  object@annotation$positiveCandidate <- .matchMz(object@xcmsData$positiveFeature,
                                                  filterPolarity(Spectra_database,1))
  object@annotation$negativeCandidate <- .matchMz(object@xcmsData$negativeFeature,
                                                  filterPolarity(Spectra_database,0))
  object

}

annotateMSdev <- function(object){

  .annotateMSdev <- function(featureMS2 , candidate){

    BiocParallel::bplapply(1:length(featureMS2),function(i){
      annotateSpectraMSdb(featureMS2[[i]],candidate[[i]])
    },BPPARAM = SerialParam(
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


getStaData <- function(object,MSDB.keys =c("Compound_name","adduct","formula","inchikey","Lipid_subclass" ,"database_origin")){

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

  featurePos <- get_features_from_xcms(object@xcmsData$positiveMS1)
  featurePos <- rowData(featurePos)%>%
    as.data.frame()%>%
    rownames_to_column("feature_id")%>%
    dplyr::select(feature_id,qc_rsd,sample_rsd,med_intensity)%>%
    dplyr::mutate(feature_id = paste0(feature_id , "_pos"),
                  ion_mode = "positive")%>%
    cbind(annotationPos,assay(featurePos))%>%
    dplyr::rename_with( ~sub(pattern = ".mzML",replacement = "",x = .x))%>%
    remove_rownames()

  featureNeg <- get_features_from_xcms(object@xcmsData$negativeMS1)
  featureNeg <- rowData(featureNeg)%>%
    as.data.frame()%>%
    rownames_to_column("feature_id")%>%
    dplyr::select(feature_id,qc_rsd,sample_rsd,med_intensity)%>%
    dplyr::mutate(feature_id = paste0(feature_id , "_neg"),
                  ion_mode = "negative")%>%
    cbind(annotationNeg,assay(featureNeg))%>%
    dplyr::rename_with( ~sub(pattern = ".mzML",replacement = "",x = .x))%>%
    remove_rownames()

  featureAll <-  rbind(featurePos,featureNeg)
  featureAllanno <- MSdb:::getInfoFromMSDB(featureAll$MSDB_id,
                                       msdb_path = object@projectInfo$MSDB_path,
                                       keys =  MSDB.keys)
  featureAll<- add_column(featureAll,featureAllanno[,-1],.after = "feature_id")
  object@statData$featureRaw <-featureAll

  object <- adjusetFeautreByweight(object)

  object@statData$feature <- featureAll%>%
    dplyr::filter(qc_rsd <0.3)
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


adjusetFeautreByweight <- function(object){

  sampleInfoToAdjust <- object@sampleInfo%>%
    dplyr::filter(!is.na(weight))
  weight <- sampleInfoToAdjust$weight / mean(sampleInfoToAdjust$weight )
  featureMatrix <- object@statData$featureRaw%>%
    column_to_rownames("feature_id")%>%
    dplyr::select(sampleInfoToAdjust$sample.name)%>%
    as.matrix()

  featureMatrixAdjusted <-t( t(featureMatrix)/weight)
  featureMatrixAdjusted-> object@statData$featureRaw[,sampleInfoToAdjust$sample.name]
  object

}



findFeature <- function(object,exact_mass =100,ppm = 10,ion_mode = 1 ){

  ion_mz <- exact_mass+ifelse(ion_mode==1 , 1.007825,-1.007825)
  ion_mode_char <- ifelse(ion_mode==1 , "positive","negative")
  feature <- object@statData$featureRaw
  feature_matched <- feature%>%
    dplyr::filter( ion_mode == ion_mode_char,
                   mz > ion_mz-ion_mz*ppm/1e6,
                   mz <  ion_mz+ion_mz*ppm/1e6)


}





