#' get xcms param by Autotuner
#'
#' @param xcms.xcms xcms
#'
#' @return Autotuner
#' @export
#' @import Autotuner
#'

get_xcms_Autotuner <- function(xcms.xcms ){

  xcms.xcms <- filterFile(xcms.xcms,
                          which(Biobase::pData(xcms.xcms)$sample.type!="Blank"))

  autotuner <- createAutotuner(
    data_paths = filepaths(xcms.xcms),
    runfile = Biobase::pData(xcms.xcms),
    file_col = "msData.files",
    factorCol = "sample.type")

  ### 1
  lag <- 25
  threshold<- 3.1
  influence <- 0.1
  signals <- lapply(getAutoIntensity(autotuner),
                    ThresholdingAlgo, lag, threshold, influence)

  autotuner <- isolatePeaks(Autotuner = autotuner,
                            returned_peaks = 10,
                            signals = signals)
  ### 2
  eicParamEsts <- EICparams(Autotuner = autotuner,
                            massThresh = .005,
                            verbose = T,
                            returnPpmPlots = FALSE,
                            useGap = TRUE)
  ###3
  returnParams(eicParamEsts, autotuner)

}






