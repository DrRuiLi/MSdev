#' check SampleInfo in excel
#' @description manually check sampleInfo using excel
#' @param object a `MSdev` object
#'
#' @return MSdev a `MSdev` object
#' @export
#'

MSdev_checkSampleInfo <- function(object){

  sampleInfo <- object@sampleInfo
  sampleInfo <- edit_df_in_excel(sampleInfo)
  ### save
  {
    object@sampleInfo <- sampleInfo
    if(!is_empty(object@xcmsData)){
      object <- MSdev_update_xcms_pdata(object )}
    object <- .updateProjectInfoFromSampleInfo(object )

  }

  object
}




#' @title msConvert_MSdev
#'
#' @param object
#'
#' @return MSdev
#' @export
#' @importFrom BiocParallel  bplapply

#'

MSdev_msConvert<- function(object){


  ### filter files
  {
    sample.info <- object@sampleInfo%>%
      dplyr::mutate(raw.exist = file.exists(raw.files),
                    ms.exist = file.exists(msData.files))%>%
      dplyr::filter(raw.exist,!ms.exist)

  }

  ### convert
  if (nrow(sample.info)) {
    MSconvertR::msConvert2mzML(raw.files  = sample.info$raw.files,
                               mzML.files = sample.info$msData.files,
                               BPPARAM = SnowParam(workers =parallel::detectCores()-1,
                                                   progressbar = T))


  }

  ### filter non converted
  {
    sample.info <- object@sampleInfo%>%
      dplyr::mutate(raw.exist = file.exists(raw.files),
                    ms.exist = file.exists(msData.files))
    object@sampleInfo <- object@sampleInfo[sample.info$ms.exist,]

  }

  object <- get_MSdev_MSinfo(object)
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



#' MSdev_extract_Spectra
#'
#' @param object  MSdev
#'
#' @return MSdev
#' @export
#'
MSdev_extract_Spectra <- function(object){

  sampleInfo <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing %in% c("Both","MS2"))%>%
    dplyr::mutate(msData.files = normalizePath(msData.files))

  if (nrow(sampleInfo)==0) {
    sp <- Spectra::Spectra()
  } else {
    sp <- Spectra::Spectra(na.omit(sampleInfo$msData.files),
                           backend = Spectra::MsBackendMemory())%>%
      filterMsLevel(2)
    sp$sp_id <- paste0("MS2_SP",num2str(1:length(sp)))
    Spectra::spectraNames(sp) <- sp$sp_id

  }

  ### iso-label
  {
    if ("isotope_label"%in% colnames(sampleInfo)) {
      sp.data <- spectraData(sp)%>%
        as.data.frame()%>%
        rownames_to_column("sp.name" )%>%
        dplyr::mutate(isotope_label = sampleInfo$isotope_label[match(
          dataOrigin , sampleInfo$msData.files
        )])
      object@spectra <- split(sp,sp.data$isotope_label)
      object@spectra$MS2_Spectra <- sp[is.na(sp.data$isotope_label)]
      return(object)
    }



  }

  object@spectra$MS2_Spectra <- sp
  return(object)

}



#' MSdev_match_Spectra_to_feature
#'
#' @param object MSdev
#'
#' @return MSdev
#' @export
#'
MSdev_match_Spectra_to_feature <- function(object){


  object@spectra$MS2_Spectra$feature_id<-NA
  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    sp.ms2 <- object@spectra$MS2_Spectra%>%
      ProtGenerics::filterPolarity(i)
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
      as.data.frame()
    sp.ms2.data <- get_Spectra_ms2_feature_id(sp.ms2,xcms.fdf)


    ### update MS2_Spectra
    sp.ms2.total <-object@spectra$MS2_Spectra %>%
      Spectra::spectraData()%>%
      as.data.frame()%>%
      dplyr::mutate(feature_id= case_when(
        polarity==i ~ sp.ms2.data[sp_id,]$feature_id,
        T~feature_id
      ))
    Spectra::spectraData(object@spectra$MS2_Spectra ) <- DataFrame(sp.ms2.total)

    ### update xcms featuredef
    xcms.fdf$ms2_id <- sapply(xcms.fdf$feature_id,
                              function(i){
                                sp_id <- sp.ms2.data$sp_id[which(sp.ms2.data$feature_id==i)]
                                return(sp_id)
                              }
    )
    xcms::featureDefinitions(xcms.xcms) <- xcms.fdf%>%
      DataFrame()
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]


  }

  return(object)


}

#' MSdev annotation
#'
#' @param object MSdev
#' @param db.path CompoundDB
#' @param ...
#'
#' @return MSdev
#' @export
#'

MSdev_annotation <- function(object,
                             cpdb_path,
                             ...){

  cpdb <- CompoundDb::CompDb(cpdb_path)

  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    message(Sys.time()," Find MS1 candidate...")
    xcms.xcms <- xcms_get_feature_ms1_candidate(xcms.xcms,
                                                cpdb,
                                   ...)
    message(Sys.time()," calculate MS2 score...")
    xcms.xcms <- xcms_get_feature_ms2_score(xcms.xcms ,
                                            cpdb = cpdb,
                                            object@spectra$MS2_Spectra,
                                            ...)
    xcms.xcms <- xcms_get_feature_annotation(xcms.xcms,
                                             ...)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]

  }

  object@projectInfo$CompoundDB_path <- cpdb_path
  return(object)


}


MSdev_get_Stat <- function(object,QC_RSD = 0.3){

  sample.info <- object@sampleInfo%>%
    dplyr::filter(polarity_paired)
  col.order <- sample.info%>%
    dplyr::distinct(sample.name)%>%
    dplyr::pull(sample.name)
  se <- list()
  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    if (is.null(xcms.xcms)) {
      se[[pol]] <- SummarizedExperiment()
      next
    }
    pol.se <- get_xcms_feature_se(xcms.xcms)
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
  feature.se <- do.call("rbind",se)

  ### sort colname
  rda <- rowData(feature.se)%>%
    as.data.frame()%>%
    dplyr::select(feature_id,mzmed,rtmed,MSDB_id, adduct,mz_ref,rt_ref,score,qc_rsd,sample_rsd,peakMaxo,
                  candidate,candidate.adduct,candidate.mz,candidate.score)

  ### all candidate
  {
    candi.rda <- rda%>%
      dplyr::mutate(candidate.n = sapply(candidate,length))
    candi.rda.split <- candi.rda[rep(candi.rda$feature_id,candi.rda$candidate.n),]%>%
      dplyr::group_by(feature_id)%>%
      dplyr::mutate(temp_id = 1:n())%>%
      dplyr::rowwise()%>%
      dplyr::mutate(MSDB_id = candidate[[temp_id]],
                    adduct = candidate.adduct[[temp_id]],
                    mz_ref = candidate.mz[[temp_id]],
                    score = candidate.score[[temp_id]])%>%
      dplyr::ungroup()
    db.info <- get_MSDB_info(candi.rda.split$MSDB_id,
                               keys = c("name","formula",
                                        "kegg_id",
                                        "inchikey","Lipid_subclass"),
                               object@projectInfo$MSdbPath)
    candi.rda.split <- candi.rda.split%>%
      dplyr::mutate(db.info,.after = rtmed)
    candi.se <- feature.se[candi.rda.split$feature_id,]
    rowData(candi.se) <- candi.rda.split

   }


  ### retrieve data
  db.info <- get_MSDB_info(rda$MSDB_id,
                                   keys = c("name","formula",
                                            "kegg_id",
                                            "inchikey","Lipid_subclass"),
                                   object@projectInfo$MSdbPath)
  rda <- rda%>%
    dplyr::mutate(db.info,.after = rtmed)
  rowData(feature.se) <- rda


  ### filter
  .uniqueFeatures <- function(score,intensity){
    score <- ifelse(score >0.3 , 10,1)
    unique.score <- score*log10(intensity)
    unique.score
  }
  rda <- rda%>%
    as.data.frame()%>%
    dplyr::filter(qc_rsd < QC_RSD,!is.na(MSDB_id))%>%
    dplyr::group_by(inchikey)%>%
    dplyr::slice_max(.uniqueFeatures(score,peakMaxo))%>%
    ungroup()
  metabolite.se <- feature.se[rda$feature_id,]


  object@statData$feature.se <- feature.se
  object@statData$candidate.se <- candi.se
  object@statData$metabolite.se <- metabolite.se
  object

}




get_MSdev_DEP_se <- function(object,
                             from = c("feature.se","metabolite.se")){

  from <- match.arg(from)
  data.se <- object@statData[[from]]

  sampleinfo <- object@sampleInfo
  ### col
  cda <- colData(data.se)%>%
    as.data.frame()%>%
    dplyr::mutate(group = sampleinfo$group[match(sample.name,sampleinfo$sample.name)],
                  condition = group,
                  sample.labels = sampleinfo$sample.labels[match(sample.name,sampleinfo$sample.name)],
                  label =sample.labels)%>%
    dplyr::group_by(condition)%>%
    dplyr::mutate(replicate = 1:n(),
                  ID = paste0(condition,num2str(1:n())))
  rownames(cda) <- cda$ID
  colData(data.se) <- cda%>%DataFrame()

  ### row
  rda <- rowData(data.se)%>%
    as.data.frame()%>%
    dplyr::mutate( label = name,
                   name = feature_id,
                   ID= feature_id)
  rowData(data.se) <- rda%>%DataFrame()

  assay(data.se) <- log2(assay(data.se))


  return(data.se)
}


MSdev_update_xcms_pdata <- function(object){


  sample_info <- object@sampleInfo
  for (i in 0:1) {
    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.pdata <- pData(xcms.xcms)%>%
      dplyr::mutate(sample_info[match(msData.files,sample_info$msData.files),  ])
    xcms.pdata -> pData(xcms.xcms)
    xcms.xcms -> object@xcmsData[[paste0(pol,"MS1")]]
  }

  return(object)


}

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





