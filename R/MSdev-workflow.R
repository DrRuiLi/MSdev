MSdev_workflow_for_lipidomic <- function(
    project.dir,
    MS_instrument = "AB6600",
    LC_condition = "Metabolomics"){


  MS_dev_obj <- MSdev(rawDataDir = paste0(project.dir,"/rawData"))
  MS_dev_obj <- checkSampleInfo(MS_dev_obj)
  MS_dev_obj <- msConvert_MSdev(MS_dev_obj)
  MS_dev_obj
  MS_dev_obj <- xcmsProcessingMSdev(MS_dev_obj,
                                    xcms.findpeak.param = dplyr::case_when(
                                      MS_instrument == "AB6600" ~ MSdev_param_set$xcms.param$findpeakparam$AB6600,
                                      MS_instrument == "QEplus" ~ MSdev_param_set$xcms.param$findpeakparam$QEplus,
                                      T ~ MSdev_param_set$xcms.param$findpeakparam$Default))
  MS_dev_obj <- extractSpectra_fullscan_DDA(MS_dev_obj)
  MS_dev_obj <- featureSpectra_fullscan_DDA(MS_dev_obj)
  MS_dev_obj <- featureCandidate(MS_dev_obj,
                                 mz.ppm = 20,
                                spectraDatabase = case_when(
                                  LC_condition == "Metabolomics" ~ "d:/MSdb/HMDB/Spectra/HMDB_spectras_experiment_2022_08_12.Rdata",
                                  LC_condition == "Lipidomics" ~ "d:/MSdb/MSdb_LipidBlast_from_MSDIAL.Rdata",
                                  ))
  MS_dev_obj <- annotateMSdev(MS_dev_obj)
  MS_dev_obj <- getStaDataMSdev(MS_dev_obj,missing = "rowmin_half",
                                MSDB.keys = case_when(
                                  LC_condition == "Metabolomics"~c("Compound_name","adduct","formula","inchikey" ,"database_origin"),
                                  LC_condition == "Lipidomics"~c("Compound_name","adduct","formula","inchikey","Lipid_subclass" ,"database_origin")
                                  )
                                )
  saveMSdev(MS_dev_obj)




return(MS_dev_obj)






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

