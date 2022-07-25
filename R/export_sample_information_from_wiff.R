#' Detect raw data and export sample information
#'
#' Detect raw data (*.wiff and *.wiff.scan) from Sciex 6600
#'
#' @param raw.data.dir directory of raw data
#' @param project.dir directory to export
#'
#' @return ms.ana objective
#' @export
#'
#' @examples
#' export_sample_information_from_wiff(raw.data.dir,project.dir)
export_sample_information_from_wiff <- function(
    raw.data.dir,
    project.dir){


  ### format dir
  {

    raw.data.dir <- gsub(raw.data.dir,pattern = "\\",replacement = "/",fixed = T)
    project.dir <- gsub(project.dir,pattern = "\\",replacement = "/",fixed = T)
    mzMl.dir <- paste0(project.dir,"mzML/")
  }

  ### check if analysis existed
  {
    ms.ana.file <- dir(path = project.dir , pattern = "ms.ana",full.names = T)
    if(!is_empty(ms.ana.file)){
      load(ms.ana.file)
      return(ms.ana)
    }



  }

  ### project information
  {
    ms.ana <- list(processing.info = list(
      project.info = list(project.dir = project.dir,
                          raw.data.dir = raw.data.dir,
                          mzMl.dir = mzMl.dir)
    ))


  }


  ### import raw data
  {
    wiff.files <- dir(path = raw.data.dir,
                      pattern = "wiff$",
                      full.names = T)

    sample.info <- data.frame(raw.files = wiff.files)%>%
      dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
      dplyr::mutate(polarity = case_when(grepl(pattern = "pos", x = raw.files,ignore.case = T)~"raw.file.positive",
                                         grepl(pattern = "neg", x = raw.files,ignore.case = T)~"raw.file.negative",
                                         T~"error"),
                    sample.abbreviation= gsub(pattern = "pos|neg|.wiff$",
                                              x = basename(raw.files) ,
                                              ignore.case = T,
                                              replacement = ""),
                    sample.abbreviation = tolower(sample.abbreviation))%>%
      pivot_wider(names_from = polarity , values_from = raw.files)%>%
      dplyr::mutate(analysis.time.positive = file.info(raw.file.positive)$mtime,
                    analysis.time.negative = file.info(raw.file.negative)$mtime)%>%
      dplyr::mutate(sample.type = case_when(grepl(pattern = "QC",x = sample.abbreviation,ignore.case = T)~ "QC",
                                            grepl(pattern = "blank|blk",x = sample.abbreviation,ignore.case = T)~ "Blank",
                                            T~"Sample"),
                    .before = sample.abbreviation)%>%
      dplyr::mutate(sample.name = paste0(sample.type,str_pad(1:nrow(.), ceiling(log10(nrow(.))),pad = "0")),
                    .before = sample.type)%>%
      dplyr::mutate(mzML.file.positive = case_when(is.na(raw.file.positive)~raw.file.positive,
                                                   T~paste0(mzMl.dir,"pos/",sample.name,".mzML")),
                    mzML.file.negative = case_when(is.na(raw.file.negative)~raw.file.positive,
                                                   T~paste0(mzMl.dir,"neg/",sample.name,".mzML")))%>%
      dplyr::select(sample.name,sample.type,sample.abbreviation,
                    raw.file.positive,raw.file.negative,
                    analysis.time.positive,analysis.time.negative,
                    mzML.file.positive,mzML.file.negative)

  }

  ### export sample.info
  {
    ms.ana$processing.info$project.info$ms.expriment.time <- min(sample.info$analysis.time.negative,
                                                                 sample.info$analysis.time.positive,
                                                                 na.rm = T)
    ms.ana$sample.info <- sample.info
    ms.ana$processing.info$raw.data.import = list(sample.no = nrow(sample.info),
                                                  positive = sum(!is.na(sample.info$raw.file.positive)),
                                                  negative = sum(!is.na(sample.info$raw.file.negative)),
                                                  time = Sys.time(),
                                                  done = T)
    sample.info.file <- paste0(project.dir,"/sample.info.xlsx")
    openxlsx::write.xlsx(list(sample.info = sample.info),
                         file = sample.info.file)
  }

  ### check if sample.info should be edited
  {
    cat(paste0("Raw data reading...\nExport sample information to\n",
               sample.info.file,
               "\nPlease check and press 1 to reload sample information,\nor any other key to continue\n"))
    input <- readline()
    if (input == 1) {
      sample.info <- readxl::read_excel(sample.info.file)

      ms.ana$processing.info$project.info$ms.expriment.time <- min(sample.info$analysis.time.negative,
                                                                   sample.info$analysis.time.positive,
                                                                   na.rm = T)
      ms.ana$sample.info <- sample.info
      ms.ana$processing.info$raw.data.import = list(sample.no = nrow(sample.info),
                                                    positive = sum(!is.na(sample.info$raw.file.positive)),
                                                    negative = sum(!is.na(sample.info$raw.file.negative)),
                                                    analysis.time = paste0(Sys.time()))
    }

  }

  ### save data
  {
    file.to.save <- paste0(project.dir,"ms.ana.",Sys.Date(),".Rdata")
    ms.ana$processing.info$project.info$ms.ana.file <- file.to.save
    save(ms.ana,file = file.to.save)
    return(ms.ana)
  }

}

