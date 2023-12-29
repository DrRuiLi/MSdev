MSdev_1.2_workflow <- function(
    project.dir,
    MS_instrument = "AB6600",
    LC_condition = "Metabolomics"){


  MS_dev_obj <- MSdev(rawDataDir = paste0(project.dir,"/rawData"))
  MS_dev_obj <- checkSampleInfo(MS_dev_obj)
  MS_dev_obj <- msConvert_MSdev(MS_dev_obj)
  MS_dev_obj
  MS_dev_obj <- xcmsProcessingMSdev(MS_dev_obj,
                                    xcms.findpeak.param = switch(MS_instrument,
                                             "AB6600" = MSdev_param_set$xcms.param$findpeakparam$AB6600,
                                             "QEplus" = MSdev_param_set$xcms.param$findpeakparam$QEplus)
  )

  MS_dev_obj <- extractSpectra_fullscan_DDA(MS_dev_obj)
  MS_dev_obj <- featureSpectra_fullscan_DDA(MS_dev_obj)
  MS_dev_obj <- featureCandidate(MS_dev_obj,
                                 mz.ppm = 20,
                                spectraDatabase = switch(LC_condition,
                                   "Metabolomics" = "d:/MSdb/msdb.HMDB.Rdata",
                                   "Lipidomics" = "d:/MSdb/MSdb_LipidBlast_from_MSDIAL.Rdata",
                                  ))
  MS_dev_obj <- annotateMSdev(MS_dev_obj)
  MS_dev_obj <- getStaDataMSdev(MS_dev_obj,missing = "rowmin_half",
                                MSDB.keys = switch(LC_condition,
                                   "Metabolomics" = c("Compound_name","adduct","formula","inchikey" ,"database_origin"),
                                  "Lipidomics" = c("Compound_name","adduct","formula","inchikey","Lipid_subclass" ,"database_origin")
                                  )
                                )
  saveMSdev(MS_dev_obj)




return(MS_dev_obj)






}



MSdev_Stat_workflow <- function(MS_dev_obj){
  exportMSdev(MS_dev_obj)
  plotMSdevPCA(MS_dev_obj)

  MS_dev_obj <- analyzeMSdevANOVA(MS_dev_obj)
  plotMSdevANOVA(MS_dev_obj)

  MS_dev_obj <- analyzeMSdevDiffMetabolites(MS_dev_obj)
  plotMSdevDiffHeatmap(MS_dev_obj)
  plotMSdevDiffVolcano(MS_dev_obj,p.adjusted = F)

  saveMSdev(MS_dev_obj)

}


MSdev_param <- function(){

  MSdev_param_set <- list(
    xcms.param = list(
    findpeakparam = list(
      AB6600 = xcms::CentWaveParam(
        ppm = 25,
        peakwidth = c(5,50),
        snthresh = 100,
        prefilter = c(3,100)),
      QEplus = xcms::CentWaveParam(
        ppm = 20,
        peakwidth = c(5,50),
        snthresh = 100,
        prefilter = c(3,1000)
      ),
      Default = xcms::CentWaveParam(
        ppm = 20,
        peakwidth = c(5,50),
        snthresh = 100,
        prefilter = c(3,1000)
      )




    )

  )
  )


}


MSdev_1.3_workflow <- function(){

  msdev.demo <- load_demo("MSdev")
  msdev.demo <- MSdev_checkSampleInfo(msdev.demo)
  msdev.demo <- MSdev_msConvert(msdev.demo)
  msdev.demo <- MSdev_xcmsProcessing(msdev.demo)
  msdev.demo <- MSdev_extract_Spectra(msdev.demo)
  msdev.demo <- MSdev_match_Spectra_to_feature(msdev.demo)
  msdev.demo <- MSdev_annotation(msdev.demo)
  msdev.demo <- MSdev_get_Stat(msdev.demo)
  MSdev_export(msdev.demo)


}



DEP_Stat_workflow <- function(){

  proj.dir <- msdev.demo@projectInfo$projectDir
  data.se <- get_MSdev_DEP_se(msdev.demo,from = "metabolite")
  p.pca <- DEP_plot_PCA(data.se)
  export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                   width = 5,height = 5)

  data.se <- DEP_normalization(data.se)
  data.se <- data.se[,data.se$sample.type=="Sample"]
  data.diff <- DEP_test_diff(data.se)
  data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

  p.diff.list <- DEP_plot_volcano(data.diff,"all")
  p.diff <- ggplot_sum_patchwork(p.diff.list)
  export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                   width = 10,height = 7,append = T)

  table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
  xlsx.write.list(table.diff,
                  paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
                  )

  diff.path <- DEP_pathway_enrich(data.diff,
                                  contrast = "all")



}



