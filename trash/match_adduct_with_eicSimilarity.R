#' @title match_adduct_with_eicSimilarity
#' @description extract feature's chromatogram and compare their similarity to main peak
#'
#' @param MS.network
#' @param xcms.xcms
#'
#' @return
#' @export
#'
#' @examples
match_adduct_with_eicSimilarity <-
  function(MS.network, xcms.xcms) {

    compare_eic <- function(x , xcms.xcms) {
      #x <- MS.network[[2]]
      compound <- x[["compound"]]
      adduct <- x[["adduct"]]


      adduct
      chrom <-
        featureChromatograms(xcms.xcms , features = adduct$feature.id)
      cor.matrix <- compareChromatograms(chrom)
      cor.to.main.peak <- cor.matrix[which.max(adduct$feature.intb), ]
      cor.to.main.peak[is.na(cor.to.main.peak)] <- 0
      adduct$cor.to.main.peak <- cor.to.main.peak

      x[["chromatogram"]] <- chrom
      x[["adduct"]] <- adduct
      return(x)

    }

    MS.network <- lapply(MS.network, compare_eic, xcms.xcms)
    return(MS.network)


  }
