# Sat Jan 13 10:43:20 2024 ------------------------------
### DCX Metabolomics
{
  msdev.dcx <- MSdev("d:/2024_01_08-Duchenxi/Data/")
  msdev.dcx <- load_as_var("d:/2024_01_08-Duchenxi/MSdev_2024_01_13.Rdata")
  msdev.dcx <- MSdev_msConvert(msdev.dcx)
  msdev.dcx <- MSdev_checkSampleInfo(msdev.dcx)
  msdev.dcx <- MSdev_xcmsProcessing(msdev.dcx)
  msdev.dcx <- MSdev_extract_Spectra(msdev.dcx)
  msdev.dcx <- MSdev_match_Spectra_to_feature(msdev.dcx)
  msdev.dcx <- MSdev_match_Spectra_to_feature(msdev.dcx)
  msdev.dcx <- MSdev_annotation(msdev.dcx,
                                expand_adduct= T,
                                selected_adduct = c("[M-H]-",
                                                    "[M-H2O-H]-",
                                                    "[2M-H]-" ,
                                                    "[M+FA-H]-" ,
                                                    "[M+H]+" ,
                                                    "[M-H2O+H]+",
                                                    "[M+NH4]+"
                                ) ,
                                db.path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/SpectraDB/HMDB_KEGG_export_human_pathway.rda")
  msdev.dcx <- MSdev_get_Stat(msdev.dcx,QC_RSD = Inf)
  MSdev_save(msdev.dcx)
  MSdev_export(msdev.dcx)


}




# Sun Jan 14 12:31:01 2024 ------------------------------
### WYQ
{
  msdev.wyq <- MSdev("d:/WYQ/2024_01_10-Wangyongqiang/Data/")
  msdev.wyq <- MSdev_msConvert(msdev.wyq)
  msdev.wyq <- MSdev_checkSampleInfo(msdev.wyq)
  msdev.wyq <- MSdev_xcmsProcessing(msdev.wyq)
  msdev.wyq <- MSdev_extract_Spectra(msdev.wyq)
  msdev.wyq <- MSdev_match_Spectra_to_feature(msdev.wyq)
  msdev.wyq <- MSdev_annotation(msdev.wyq,
                                expand_adduct= T,
                                selected_adduct = c("[M-H]-",
                                                    "[M-H2O-H]-",
                                                    "[2M-H]-" ,
                                                    "[M+FA-H]-" ,
                                                    "[M+H]+" ,
                                                    "[M-H2O+H]+",
                                                    "[M+NH4]+"
                                ) ,
                                db.path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/SpectraDB/HMDB_KEGG_export_human_pathway.rda")
  msdev.wyq <- MSdev_get_Stat(msdev.wyq)
  MSdev_save(msdev.wyq)
  MSdev_export(msdev.wyq)



}
