#' Title
#'
#' @param MS.network
#' @param xcms.features
#' @param ppm.thresh
#'
#' @return
#' @export
#'
#' @examples
match_adduct_to_features <-
  function(MS.network , xcms.xcms, ppm.thresh = 10,rt.tol = 10)
  {

    xcms.features <- featureDefinitions(xcms.xcms)%>%as.data.frame()
    xcms.features.intb <- apply(featureValues(xcms.xcms , missing = "rowmin_half"),
                                1,median)
    match.adduct <-function(x){
     # x <- MS.network[[2]]
      adduct.candidate <- x[["adduct.candidate"]]
      adduct.mz <- adduct.candidate$exact.mz
      peak.mz <- xcms.features$mzmed
      peak.id <-rownames(xcms.features)
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
                            peak.mz = xcms.features[matched.id[,2],"mzmed"],
                            peak.id = peak.id[matched.id[,2]],
                            peak.rt = xcms.features[matched.id[,2],"rtmed"],
                            peak.intb = xcms.features.intb[matched.id[,2]])%>%
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
