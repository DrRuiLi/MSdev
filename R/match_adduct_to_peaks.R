#' @title match_adduct_to_peaks
#' @description match adduct candidate in xcms.peaks table
#' @param MS.network
#' @param ppm.thresh
#' @param xcms.peaks
#'
#' @return
#' @export
#'
#' @examples
match_adduct_to_peaks <- function(MS.network , xcms.peaks, ppm.thresh = 10){

match.adduct <-function(x){
  adduct.candidate <- x[["adduct.candidate"]]
  adduct.mz <- adduct.candidate$exact.mz
  peak.mz <- xcms.peaks[,"mz"]
  adduct.matrix <- matrix(rep(adduct.mz,length(peak.mz)) , nrow = length(adduct.mz))
  peak.matrix <- matrix(rep(peak.mz,length(adduct.mz)) ,
                        nrow = length(adduct.mz),
                        byrow = T)
  sub.matrix <- adduct.matrix - peak.matrix
  tol.matrix <- peak.matrix * ppm.thresh*1e-6
  pass.matrix <- abs(sub.matrix) < tol.matrix

  matched.id <- which(pass.matrix,arr.ind = T)
  matched.id

  adduct <- data.frame( adduct.candidate[matched.id[,1],],
                        peak.mz = xcms.peaks[matched.id[,2],"mz"],
                        peak.rt = xcms.peaks[matched.id[,2],"rt"],
                        peak.intb = xcms.peaks[matched.id[,2],"intb"])

  x[["adduct"]] <- adduct
  return(x)
}

MS.network <- lapply(MS.network ,match.adduct )
return(MS.network)

}
