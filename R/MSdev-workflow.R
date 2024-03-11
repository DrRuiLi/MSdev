MSdev_1.4_workflow <- function(){

  msdev.demo <- load_demo("MSdev")
  msdev.demo <- MSdev_checkSampleInfo(msdev.demo)
  msdev.demo <- MSdev_msConvert(msdev.demo)
  get_MSdev_param(msdev.demo)
  msdev.demo <- MSdev_xcmsProcessing(msdev.demo)
  msdev.demo <- MSdev_extract_Spectra(msdev.demo)
  msdev.demo <- MSdev_match_Spectra_to_feature(msdev.demo)
  msdev.demo <- MSdev_annotation(msdev.demo)
  msdev.demo <- MSdev_get_Stat(msdev.demo)
  MSdev_export(msdev.demo)
  MSdev_save(msdev.demo)


}



MSdev_Stat_workflow <- function(MS_dev_obj){
  exportMSdev(MS_dev_obj)
  plotMSdevPCA(MS_dev_obj)

  MS_dev_obj <- analyzeMSdevANOVA(MS_dev_obj)
  plotMSdevANOVA(MS_dev_obj)

  MS_dev_obj <- analyzeMSdevDiffMetabolites(MS_dev_obj)
  plotMSdevDiffHeatmap(MS_dev_obj)
  plotMSdevDiffVolcano(MS_dev_obj,p.adjusted = F)

  MSdev_save(MS_dev_obj)

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




DEP_Stat_workflow <- function(){

  proj.dir <- msdev.demo@projectInfo$projectDir
  data.se <- get_MSdev_DEP_se(msdev.demo,from = "metabolite")
  p.pca <- DEP_plot_PCA(data.se)
  export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                   width = 5,height = 5)

  data.se <- DEP_normalization(data.se)
  data.se <- data.se[,!data.se$sample.type %in% c("Blank","QC")]
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



