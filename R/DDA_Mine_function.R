MSdev_get_Inclusion_Queue <- function(object){

  for (i in 0:1) {
    polarity <-ifelse(i==0,"Negative","Positive")
    polarity.tag <- paste0(polarity,"MS1")
    xcms.xcms <- object@xcmsData[[polarity.tag]]
    if (is.null(xcms.xcms) ) next
    feature.rsd <- get_features_from_xcms(xcms.xcms)@elementMetadata%>%as.data.frame()
    xcms.xcms <- xcms_get_feature_def_stat(xcms.xcms)
    feature.stat <- get_xcms_feature_definitions(xcms.xcms)
    dda.mine.queue <- feature.stat%>%
      dplyr::mutate( qc_rsd.score = (log(0.3)-log(qc_rsd)),
                     MS1.score = qc_rsd.score*log10(peakMaxo)  )%>%
      dplyr::arrange(-MS1.score)%>%
      tibble::rownames_to_column("feature.id")%>%
      dplyr::mutate(CE10=10,
                    CE20=20,
                    CE30=30,
                    CE40=40,
                    CE50=50
      )%>%
      tidyr::pivot_longer(CE10:CE50,names_to = "CE.tag",values_to = "collisionEnergy")%>%
      dplyr::filter(collisionEnergy %in% c(30))%>%
      dplyr::mutate(DDA.id= paste0(feature.id,"_",CE.tag),
                    acquired =F,
                    acquired.in.list = "",
                    queued.in.list = "",
                    queued.time = 0)

    object@advancedAna[[paste0("DDA_mine_queue_",polarity)]] <- dda.mine.queue
    object@advancedAna[[paste0("DDA_mine_list_",polarity)]] <- list()

  }


  object

}

MSdev_get_Inclusion_List <- function(object){

  for (i in 0:1) {
    polarity <-ifelse(i==0,"Negative","Positive")
    DDA.queue <- object@advancedAna[[paste0("DDA_mine_queue_",polarity)]]
    if ( is.null(DDA.queue)) next

    DDA.mine.list <- DDA.queue%>%
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
    queue.list <- object@advancedAna[[paste0("DDA_mine_list_",polarity)]]
    if (length(queue.list)) {

      list.name =str_add( max(names(queue.list)),1)
      queue.list <- append(queue.list,
                           list( DDA.mine.list))
      names(queue.list)[length(queue.list)] <-list.name

    }else{
      list.name <- "DDA_mine_list001"
      queue.list <- list("DDA_mine_list001" = DDA.mine.list)

    }
    object@advancedAna[[paste0("DDA_mine_list_",polarity)]] <-queue.list

    DDA.mine.list.qe <- QE_list_2feature_def(DDA.mine.list)
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
    DDA.queue -> object@advancedAna[[paste0("DDA_mine_queue_",polarity)]]

  }


  object

}

MSdev_get_MS2acquisitionStat <- function(object){

  assign_ms2_list <- function(pmz,rt,ce ,il){


    idx.mz <- match_mz(pmz,il$mzmed,mz.ppm = 5)
    rt.match <- between(rt, il$peakRtMin[idx.mz]-30,
                        il$peakRtMax[idx.mz]+30)

    rt.match[!rt.match ] <-NA
    idx.mz[rt.match]

    x <- il$DDA.id[idx.mz[rt.match]]
    if (length(x)==0) {
      return(NA)

    }
    return(x)

  }

  for (i in 0:1) {

    polarity <-ifelse(i==0,"Negative","Positive")
    DDA.queue <- object@advancedAna[[paste0("DDA_mine_queue_",polarity)]]
    if ( is.null(DDA.queue)) next
    sample.info <- object@sampleInfo%>%
      dplyr::filter(polarity %in% c(i,-1),
                    msLevels %in% c(2))
    xcms.xcms <- readMSData(sample.info$msData.files,mode = "onDisk")
    xcms.scan <- get_xcms_scan_Stat(xcms.xcms)%>%
      dplyr::filter(msLevel==2)%>%
      #dplyr::rowwise()%>%
      dplyr::mutate(assigned.id = assign_ms2_list(pmz = precursorMZ,
                                                  rt = retentionTime,
                                                  ce = collisionEnergy,
                                                  il = DDA.queue)) %>%
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
                      acquired.in.list = case_when(
                        DDA.id %in% ids ~ paste0( acquired.in.list,
                                                  ";",sampleNames(xcms.xcms)[i]),
                        T~ acquired.in.list))%>%
        dplyr::ungroup()%>%
        dplyr::mutate(fail.time = case_when(acquired~0,
                                            T~queued.time))%>%
        dplyr::group_by(feature.id)%>%
        dplyr::mutate(
          acquired = case_when(any(fail.time >3)~T,
                               T~acquired))

    }

    DDA.queue -> object@advancedAna[[paste0("DDA_mine_queue_",polarity)]]



  }


  return(object)

}

