
edit_sample_info<- function(ms.ana) {

  sample.info <- ms.ana$sample.info
  sample.info <- edit_df_in_excel(sample.info)
  tibble::as_tibble(sample.info) -> ms.ana$sample.info

  {### save

    ms.ana[["sample.info"]] <-sample.info
    pData(ms.ana$xcms.positive) <- cbind(sample.info,sampleNames(ms.ana$xcms.positive))
    pData(ms.ana$xcms.negative) <- cbind(sample.info,sampleNames(ms.ana$xcms.negative))
    save_ms_ana( ms.ana)
  }
  return(ms.ana)

}


save_ms_ana <- function(ms.ana){

    save(ms.ana , file = ms.ana[["processing.info"]][["project.info"]][["ms.ana.file"]])


}

get_feature <- function(ms.ana,qc_rsd_thresh = 0.3){



  {### extract feature data
    feature.data.pos <- get_features_from_xcms(ms.ana$xcms.positive)
    feature.pos <- cbind(SummarizedExperiment::rowData(feature.data.pos),
                         SummarizedExperiment::assay(feature.data.pos))%>%
      tibble::as_tibble()%>%
      select_all(gsub,pattern = ".mzML", replacement = "")%>%
      select(-c(peakidx))%>%
      cbind(ms.ana$annotation.positive$annotation.table,.)%>%
      mutate(feature.id = paste0(feature.id , "_pos"))

    feature.data.neg <- get_features_from_xcms(ms.ana$xcms.negative)
    feature.neg <- cbind(SummarizedExperiment::rowData(feature.data.neg),
                         SummarizedExperiment::assay(feature.data.neg))%>%
      tibble::as_tibble()%>%
      select_all(gsub,pattern = ".mzML", replacement = "")%>%
      select(-c(peakidx))%>%
      cbind(ms.ana$annotation.negative$annotation.table,.)%>%
      mutate(feature.id = paste0(feature.id , "_neg"))
    feature <- rbind(feature.pos,
                     feature.neg)%>%
      filter(qc_rsd < qc_rsd_thresh)
    message("Total ",
            nrow(feature),"/",
            nrow(feature.pos)+nrow(feature.neg),
            " features RSD < ",qc_rsd_thresh*100," %")

  }

  {###save
    ms.ana[["feature"]] <-feature
    save_ms_ana( ms.ana)

  }


  return(ms.ana)
}

get_feature_corrected <- function(ms.ana){

  if (isTRUE(ms.ana$processing.info$extract.feature$weight.correct)) {
    message("Weight Corrected ")
    return(ms.ana)
  }
  if (is.null(ms.ana$sample.info$weight)) {
    message("No weight to correct")
    return(ms.ana)
  }

  {### correct weight
    sample.info <- ms.ana$sample.info
    to.correct <- sample.info%>%
      filter(!is.na(weight))
    feature.matrix <-ms.ana$feature%>%
      select(to.correct$sample.name)
    feature.matrix <- apply(feature.matrix , 1 , function(x , weight){

      x / weight
    },weight = to.correct$weight/mean(to.correct$weight))%>%t

    ms.ana$feature[,to.correct$sample.name] <- feature.matrix
  }


  {### save
    message(Sys.time(),", Correct weight")
    ms.ana[["feature"]] <-feature
    ms.ana$processing.info$extract.feature$weight.correct <- T
    save_ms_ana( ms.ana)


  }
  return(ms.ana)

}
get_unique_compound <- function(ms.ana){

  {### check compound
    if (!is.null(ms.ana[["compound"]]) ) {
      ms.ana[["compound"]] ->compound

      return(ms.ana)
    }
  }
  {### get unique

    compound <- ms.ana$feature%>%
      filter(!is.na(inchikey),
             score > 0.3)%>%
      group_by(inchikey)%>%
      slice_max(score*10+log10(med_intensity))%>%
      ungroup()
  }


  {### save

    ms.ana[["compound"]] <-compound
    save_ms_ana( ms.ana)
  }

  return(ms.ana)
}



plot_feature_intensity_distribution<-function(ms.ana,feature.id.to.plot){

  sample.info <- ms.ana$sample.info%>%
    dplyr::arrange(analysis.time.positive)%>%
    dplyr::mutate(injection.order = 1:nrow(.))
  feature <- ms.ana$feature

  sample.info$value <- feature %>%
    dplyr::filter(feature.id == feature.id.to.plot)%>%
    dplyr::select(sample.info$sample.name)%>%
    as.numeric()
  qc.value <- sample.info$value[sample.info$sample.type == "QC"]
  qc.rsd <- sd(qc.value,na.rm = T)/mean(qc.value,na.rm = T)*100
  sample.value <- sample.info$value[sample.info$sample.type != "Blank"]
  sample.rsd <- sd(sample.value,na.rm = T)/mean(sample.value,na.rm = T)*100
  ggplot(sample.info)+
    geom_point(aes(x =injection.order , y = value ,col = sample.type))+
    scale_color_manual(values = c("QC"= "green","Blank" = "grey","Sample" = "blue"))+
    labs(title = paste0(feature.id.to.plot),
         subtitle = paste0("mz = ",feature$mz[feature$feature.id == feature.id.to.plot]%>%
                             sprintf("%.4f",.)," ; rt = ",
                           feature$rt[feature$feature.id == feature.id.to.plot]%>%
                             sprintf("%.2f",.),"\n",
                           "QC RSD = ",sprintf("%.2f",qc.rsd),"%\n",
                           "Sample RSD = ",sprintf("%.2f",sample.rsd),"%"),
         col = "Sample type")


}
