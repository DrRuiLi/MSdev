#' @title metabolomic_workflow
#'
#' @param project.dir
#' @param raw.data.dir
#'
#' @return
#' @export
#'
#' @examples
metabolomic_workflow <- function(project.dir = "d:/2022_07_05-Lirui/",
                                 raw.data.dir = "d:/2022_07_05-Lirui/Data/") {
  thread <- parallel::detectCores() - 1




  ### MS file convert and information
  {
    ms.ana <-
      export_sample_information_from_wiff(raw.data.dir = raw.data.dir,
                                          project.dir = project.dir)

    #ms.ana$sample.info <- readxl::read_excel(paste0(project.dir,"/sample.info.xlsx"))

    ms.ana <- ms_convert(ms.ana)

  }

  ### xcms processing
  {
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "positive")
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "negative")
    ms.ana <- extract_spectra(ms.ana)

  }
  ### annotation
  {
    ms.ana <- annotation_by_database(ms.ana, polarity = "positive")
    ms.ana <- annotation_by_database(ms.ana, polarity = "negative")

  }
}

lipidomic_workflow <- function(project.dir = "d:/2022_07_05-Lirui/",
                               raw.data.dir = "d:/2022_07_05-Lirui/Data/"){
  thread <- parallel::detectCores() - 1

  ### MS file convert and information
  {
    ms.ana <-
      export_sample_information_from_wiff(raw.data.dir = raw.data.dir,
                                          project.dir = project.dir)

    #ms.ana$sample.info <- readxl::read_excel(paste0(project.dir,"/sample.info.xlsx"))

    ms.ana <- ms_convert(ms.ana)

  }

  ### xcms processing
  {
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "positive")
    ms.ana <- xcms_dda_processing(ms.ana ,  polarity = "negative")
    ms.ana <- extract_spectra(ms.ana)

  }
  ### annotation
  {
    ms.ana <- annotation_by_database(ms.ana,
                                     polarity = "positive",
                                     database.file = "d:/MSdb_temp/LipidBlast/LipidBlast_spectra_2022_08_12.Rdata")
    ms.ana <- annotation_by_database(ms.ana, polarity = "negative")

  }


}











