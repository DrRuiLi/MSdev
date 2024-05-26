#' @title match_adduct_to_peaks
#' @description match adduct candidate in xcms.peaks table
#'
#' @param MS.network
#' @param xcms.xcms
#' @param ppm.thresh
#'
#' @return
#' @export
#'
#' @examples
match_adduct_to_peaks <- function(MS.network , xcms.xcms, ppm.thresh = 10,rt.tol = 10){

  xcms.peaks <-chromPeaks(xcms.xcms)
match.adduct <-function(x){
  #x <- MS.network[[2]]
  adduct.candidate <- x[["adduct.candidate"]]
  adduct.mz <- adduct.candidate$exact.mz
  peak.mz <- xcms.peaks[,"mz"]
  peak.id <-rownames(xcms.peaks)
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
                        peak.id = peak.id[matched.id[,2]],
                        peak.rt = xcms.peaks[matched.id[,2],"rt"],
                        peak.intb = xcms.peaks[matched.id[,2],"intb"])%>%
    dplyr::mutate(peak.error = (peak.mz-exact.mz)/exact.mz*1e6, .before = peak.id)
  rt_max <- adduct$peak.rt[which.max(adduct$peak.intb)]
  adduct <- adduct%>%
    dplyr::mutate(rt.filter = case_when(abs(peak.rt - rt_max) <rt.tol ~ T,
                                        T ~F ))

  x[["adduct"]] <- adduct
  return(x)
}

MS.network <- lapply(MS.network ,match.adduct )
return(MS.network)

}
