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
      feature.mz <- xcms.features$mzmed
      feature.id <-rownames(xcms.features)
      adduct.matrix <- matrix(rep(adduct.mz,length(feature.mz)) , nrow = length(adduct.mz))
      feature.matrix <- matrix(rep(feature.mz,length(adduct.mz)) ,
                            nrow = length(adduct.mz),
                            byrow = T)
      sub.matrix <- adduct.matrix - feature.matrix
      tol.matrix <- feature.matrix * ppm.thresh*1e-6
      pass.matrix <- abs(sub.matrix) < tol.matrix

      matched.id <- which(pass.matrix,arr.ind = T)
      matched.id

      adduct <- data.frame( adduct.candidate[matched.id[,1],],
                            feature.mz = xcms.features[matched.id[,2],"mzmed"],
                            feature.id = feature.id[matched.id[,2]],
                            feature.rt = xcms.features[matched.id[,2],"rtmed"],
                            feature.intb = xcms.features.intb[matched.id[,2]])%>%
        dplyr::mutate(feature.error = (feature.mz-exact.mz)/exact.mz*1e6, .before = feature.id)
      rt_max <- adduct$feature.rt[which.max(adduct$feature.intb)]
      adduct <- adduct%>%
        dplyr::mutate(rt.filter = case_when(abs(feature.rt - rt_max) <rt.tol ~ T,
                                            T ~F ))

      x[["adduct"]] <- adduct
      return(x)
    }

    MS.network <- lapply(MS.network ,match.adduct )
    return(MS.network)










  }
