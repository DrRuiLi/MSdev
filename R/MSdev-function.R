load_demo <- function(){

  var <- load( "C:/Users/91879/OneDrive/Documents/Code/R/Projecct/2022.1.8_MS.demo/Demo3/MSdev_2022_10_15.Rdata")
  if (length(var)!=1) {
    stop("Too many variabls in ",file_to_load)
  }
  eval(str2expression(var))

}
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


#' @title msConvert_MSdev
#'
#' @param object
#'
#' @return
#' @export
#' @importFrom BiocParallel  bplapply
#' @examples
#'

msConvert_MSdev <- function(object){

  ### convert
  if (object@projectInfo$rawDataFormat == ".raw") {
    MSconvertR::msConvert2mzML(raw.files  = object@sampleInfo$raw.file.positive,
                       mzML.files = object@sampleInfo$msData.file.positive,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

    MSconvertR::msConvert2mzML(raw.files  = object@sampleInfo$raw.file.negative,
                       mzML.files = object@sampleInfo$msData.file.negative,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

  }else if(object@projectInfo$rawDataFormat == ".wiff"){
    MSconvertR::msConvert2mzML(wiff.files  = object@sampleInfo$raw.file.positive,
                       mzML.files = object@sampleInfo$msData.file.positive,
                       BPPARAM = SnowParam(workers =parallel::detectCores()-1 ))

    MSconvertR::msConvert2mzML(wiff.files  = object@sampleInfo$raw.file.negative,
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

xcmsProcessingMSdev <- function(object,
                                xcms.findpeak.param = xcms::CentWaveParam(ppm = 10,snthresh = 100,
                                                                   peakwidth = c(5,20),
                                                                   prefilter = c(3,1000))){

  ### determine xcms param
  {
    if (is.na(xcms.findpeak.param)) {

    }

  }

  sampleInfoPos <-dplyr::filter(object@sampleInfo,
                                xcmsProcessing %in% c("MS1","Both")
                                )%>%
    dplyr::filter(!is.na(msData.file.positive))
  object@xcmsData$positiveMS1  <- xcmsProcessingMS1(msDataFiles = sampleInfoPos$msData.file.positive,
                                                       ion_mode = 1,
                                                       peaksGroup =sampleInfoPos$sample.type,
                                                    centWaveParam = xcms.findpeak.param
  )

  Biobase::pData(object@xcmsData$positiveMS1) <- cbind(Biobase::pData(object@xcmsData$positiveMS1),
                                              sampleInfoPos)

  sampleInfoNeg <-dplyr::filter(object@sampleInfo,
                                xcmsProcessing %in% c("MS1","Both")
  )%>%
    dplyr::filter(!is.na(msData.file.negative))
  object@xcmsData$negativeMS1  <- xcmsProcessingMS1(msDataFiles = sampleInfoNeg$msData.file.negative,
                                                       ion_mode = 0,
                                                       peaksGroup = sampleInfoNeg$sample.type,
                                                       centWaveParam = xcms.findpeak.param
  )

  Biobase::pData(object@xcmsData$negativeMS1) <-  cbind(Biobase::pData(object@xcmsData$negativeMS1),
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
                            MSDB.keys =c("Compound_name","adduct","formula","inchikey","Lipid_subclass" ,"database_origin")
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
                             "Compound: ",feature.data$Compound_name,"\n",
                             "Adduct: ",feature.data$adduct,"\n",
                             "Ref Precursor mz: ",sp.ref$precursorMz,"\n",
                             "Ref Retention time: ",sp.ref$rtime,"\n",
                             "INCHIKEY: ",sp.ref$inchikey,"\n",
                             "KEGG ID: ",sp.ref$kegg.id,"\n",
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
                             "Compound: ",feature.data$Compound_name,"\n",
                             "Adduct: ",feature.data$adduct,"\n",
                             "Ref Precursor mz: ",sp.ref$precursorMz,"\n",
                             "Ref Retention time: ",sp.ref$rtime,"\n",
                             "INCHIKEY: ",sp.ref$inchikey,"\n",
                             "KEGG ID: ",sp.ref$kegg.id,"\n",
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
                                       feature.id = gsub(x = feature_id,pattern = "[_|pos|neg]",replacement = ""))
  export::graph2png(gp , file = paste0(out.dir,"/",feature_id,".Crhom.png"))

}





