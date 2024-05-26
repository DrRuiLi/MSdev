
MSpectra2Spectra <- function(msp){

  msp.DF <- DataFrame( rtime = rtime(msp),
                       precursorMz = precursorMz(msp),
                       precursorCharge = precursorCharge(msp),
                       precScanNum = precScanNum(msp),
                       precursorIntensity = precursorIntensity(msp),
                       acquisitionNum = acquisitionNum(msp),
                       scanIndex = scanIndex(msp),
                       peaksCount = peaksCount(msp),
                       msLevel = msLevel(msp),
                       tic = tic(msp),
                       ionCount = ionCount(msp),
                       collisionEnergy = collisionEnergy(msp),
                       fromFile = fromFile(msp),
                       polarity = polarity(msp),
                       smoothed = smoothed(msp),
                       isEmpty = isEmpty(msp),
                       centroided = centroided(msp),
                       isCentroided = isCentroided(msp)
  )
  msp.DF$mz <- mz(msp)
  msp.DF$intensity <- intensity(msp)

  if (!is.null(msp@elementMetadata@listData[["feature_id"]])) {
    msp.DF$feature_id <- msp@elementMetadata@listData[["feature_id"]]
  }

  if (!is.null(msp@elementMetadata@listData[["peak_id"]])) {
    msp.DF$peak_id <- msp@elementMetadata@listData[["peak_id"]]
  }

  sp <- Spectra(msp.DF)
  return(sp)

}
