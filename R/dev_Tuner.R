#' @describeIn xcms_extenstion get xcms parameters via Autotuner
#' @title Get XCMS Parameters via Autotuner
#' @description Uses the Autotuner algorithm to automatically estimate XCMS peak detection parameters from an xcms object.
#'
#' @param xcms.xcms An xcmsSet object (or similar) containing chromatographic data.
#'
#' @return A list of XCMS parameters as returned by Autotuner's `returnParams`.
#' @export

get_xcms_Autotuner <- function(xcms.xcms ){

  xcms.xcms <- xcms::filterFile(xcms.xcms,
                          which(Biobase::pData(xcms.xcms)$sample.type!="Blank"))

  autotuner <- createAutotuner(
    data_paths = xcms::filepaths(xcms.xcms),
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






