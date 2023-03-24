annotateSpectra <- function(expSpec,refSpec){

  to.return <- list(
    mz = median(ProtGenerics::precursorMz(expSpec)),
    rt = median(ProtGenerics::rtime(expSpec)),
    ref.mz = NA,
    ref.rt = NA,
    score = 0,
    compound = NA,
    adduct = NA,
    inchikey = NA,
    kegg.id = NA,
    origin = NA,
    "expSpec" = expSpec
  )
  if (length(refSpec)==0) {
    return(to.return)
  }

  mz.error <- abs(matrixSub(expSpec$precursorMz , refSpec$precursorMz))/to.return$mz*1e6
  mz.score <- 1 - mz.error/20 ###ppm 20

  rt.error <- abs(matrixSub(expSpec$rtime , refSpec$rtime))
  rt.score <- 2- rt.error/20
  rt.score[is.na(rt.score)] <- 0

  sp.score <- Spectra::compareSpectra(expSpec,refSpec)
  sp.score[is.nan(sp.score)] <- 0

  score <- mz.score*0.2 + rt.score*0.2 + sp.score*0.6


    refSpecMatched <- refSpec[ceiling(which.max(score)/length(expSpec))]
    score <- sp.score[which.max(score)]
    to.return$ref.mz <- refSpecMatched$precursorMz
    to.return$ref.rt <- refSpecMatched$rtime
    to.return$score <- score
    to.return$compound <- refSpecMatched$name
    to.return$adduct <- refSpecMatched$adduct
    to.return$inchikey <- refSpecMatched$inchikey
    to.return$kegg.id <- refSpecMatched$kegg.id
    to.return$origin <- refSpecMatched$database
    to.return$refSpec <- refSpecMatched
  return(to.return)




}

annotateSpectraMSdb <- function(expSpec,refSpec){

  to.return <- list(
    mz = median(ProtGenerics::precursorMz(expSpec)),
    rt = median(ProtGenerics::rtime(expSpec)),
    ref.mz = NA,
    ref.rt = NA,
    score = 0,
    MSDB_id =NA,
    "expSpec" = expSpec,
    "refSpec" = refSpec
  )
  if (length(refSpec)==0) {
    return(to.return)
  }

  mz.error <- abs(matrixSub(expSpec$precursorMz , refSpec$precursorMz))/to.return$mz*1e6
  mz.score <- 1 - mz.error/20 ###ppm 20

  rt.error <- abs(matrixSub(expSpec$rtime , refSpec$rtime))
  rt.score <- 2- rt.error/20
  rt.score[is.na(rt.score)] <- 0

  sp.score <- Spectra::compareSpectra(expSpec,refSpec)
  sp.score[is.nan(sp.score)] <- 0

  score <- mz.score*0.2 + rt.score*0.2 + sp.score*0.6


  refSpecMatched <- refSpec[ceiling(which.max(score)/length(expSpec))]
  score <- sp.score[which.max(score)]
  to.return$ref.mz <- refSpecMatched$precursorMz
  to.return$ref.rt <- refSpecMatched$rtime
  to.return$score <- score
  to.return$MSDB_id <-refSpecMatched$MSDB_id
  to.return$refSpec <- refSpecMatched
  return(to.return)




}



makeSpectra <- function(precursorMz ,
                        rtime ,...){

  Spectra::Spectra(S4Vectors::DataFrame(precursorMz = precursorMz,
                             rtime = rtime,
                             ...))

}


filterSpectraIntensity <- function(sp,r){

  Spectra::filterIntensity( sp, intensity =c(max(intensity(sp))*r,Inf   ) )

}



normalizeSpectra <- function(sp){

  nf <-  function(z, ...) {
    #z[,"maxIntensity"] <- max( z[, "intensity"] )
    z[, "intensity"] <- z[, "intensity"] /
      max(z[, "intensity"], na.rm = TRUE) * 100
    z
  }
  sp <- Spectra::addProcessing(sp,nf)

  return(sp)

}
