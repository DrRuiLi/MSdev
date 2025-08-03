load_demo <- function(demo = c("MSdev",
                               "XCMSnExp","xcms",
                               "SummarizedExperiment","data.se",
                               "Spectra","sp")){

  demo <-match.arg(demo)
  file_path <- switch(demo,
         "MSdev" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/MSdev_2024_05_30.Rdata",
         "XCMSnExp" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/XCMSnExp_2023_11_17.rda",
         "xcms" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/XCMSnExp_2023_11_17.rda",
         "SummarizedExperiment" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/SummarizedExperiment_2023_11_03.rda",
         "data.se" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/SummarizedExperiment_2023_11_03.rda",
         "Spectra" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/Spectra_2023_11_23.rda",
         "sp" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/Spectra_2023_11_23.rda"
         )

  fun <- switch(demo,
                "MSdev" = MSdev_load,
                "XCMSnExp" = readRDS,
                "xcms" = readRDS,
                "SummarizedExperiment" = readRDS,
                "data.se" = readRDS,
                "Spectra" = readRDS,
                "sp" = readRDS
  )
  demo <- fun(file_path)

  return(demo)
}

make_demo <- function(){


  msdev.demo <- MSdev("C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/raw.data/")
  msdev.demo <- MSdev_msConvert(msdev.demo)
  msdev.demo <- MSdev_xcmsProcessing(msdev.demo)
  msdev.demo <- MSdev_extract_Spectra(msdev.demo)
  msdev.demo <- MSdev_match_Spectra_to_feature(msdev.demo)
  msdev.demo <- MSdev_annotation(msdev.demo,
                                 cpdb_path ="c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb" )
  MSdev_save(msdev.demo)


}


