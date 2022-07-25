#' Extract Spectra
#'
#' Extract feature MSpectra from XCMSnExp and convert to Spectra
#'
#' @param ms.ana
#'
#' @return ms.ana
#' @export
#'
#' @examples
extract_spectra <- function(ms.ana){

  message(Sys.time()," Extract Spectra...")
  ### check done
  if (isTRUE(ms.ana$processing.info$extract.spectra$done)) {
    return(ms.ana)
  }

  ### extract spectra
  {
    bp.param <- bpparam()
    register(SerialParam())
    ms.ana$Spectra.positive$MSpectra <- featureSpectra(ms.ana$xcms.positive,msLevel = 2,return.type = "MSpectra")
    try(ms.ana$Spectra.positive$Spectra <- MSpectra2Spectra(ms.ana$Spectra.positive$MSpectra)%>%
          combineSpectra(f = .$feature_id,
                         peaks = "intersect",
                         minProp = 0.3,
                         ppm = 20,
                         BPPARAM =SnowParam(workers = 31) ))
    if (!is.null(ms.ana$Spectra.positive$Spectra)) {
      ms.ana$Spectra.positive$MSpectra <- NULL
    }

    ms.ana$Spectra.negative$MSpectra <- featureSpectra(ms.ana$xcms.negative,msLevel = 2,return.type = "MSpectra")
    try(ms.ana$Spectra.negative$Spectra <- MSpectra2Spectra(ms.ana$Spectra.negative$MSpectra )%>%
          combineSpectra(f = .$feature_id,
                         peaks = "intersect",
                         minProp = 0.5,
                         ppm = 20))
    if (!is.null(ms.ana$Spectra.negative$Spectra)) {
      ms.ana$Spectra.negative$MSpectra <- NULL
    }
  }
  register(bp.param)
  ### save and return
  {

    ms.ana$processing.info$extract.spectra$done <- T
    ms.ana$processing.info$extract.spectra$time <- Sys.time()
    save(ms.ana , file = ms.ana$processing.info$project.info$ms.ana.file)

    return(ms.ana)

  }


}

