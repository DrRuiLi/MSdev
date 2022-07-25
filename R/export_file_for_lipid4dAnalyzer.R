
export_file_for_lipid4dAnalyzer <- function(ms.ana){


  sample.info <- ms.ana$sample.info
  ms1.table.pos <- cbind(
    featureDefinitions(ms.ana$xcms.positive),
    featureValues(ms.ana$xcms.positive)
  )%>%
    as.data.frame()%>%
    rename_all(gsub,pattern = ".mzML",replacement = "")%>%
    dplyr::mutate(mz = mzmed,
                  rt = rtmed/60,
                  ccs = NA
    )%>%
    dplyr::select(mz,rt,ccs,sample.info$sample.name)

  file.to.save <- paste0(project.dir,"Lipid4DAnalyzer/positive/ms1.table.csv")
  dir.create(dirname(file.to.save),recursive = T)
  write.csv(ms1.table.pos , file = file.to.save , row.names = F)

  ms1.table.neg <- cbind(
    featureDefinitions(ms.ana$xcms.negative),
    featureValues(ms.ana$xcms.negative)
  )%>%
    as.data.frame()%>%
    rename_all(gsub,pattern = ".mzML",replacement = "")%>%
    dplyr::mutate(mz = mzmed,
                  rt = rtmed/60,
                  ccs = NA
    )%>%
    dplyr::select(mz,rt,ccs,sample.info$sample.name)

  file.to.save <- paste0(project.dir,"Lipid4DAnalyzer/negative/ms1.table.csv")
  dir.create(dirname(file.to.save),recursive = T)
  write.csv(ms1.table.neg , file = file.to.save , row.names = F)




}

