#' @title  get_features_from_xcms
#' @description extract feature data from xcms::XCMSnExp,
#'  calculate RSD of QC and Sample
#'  ( note this rely on character "QC" and "Sample" in `sampleNames(xcms.xcms)` )
#' @param xcms.xcms
#'
#' @return
#' @export
#'
#' @examples
get_features_from_xcms <- function(xcms.xcms){

  xcms.sum <- quantify(xcms.xcms,missing = "rowmin_half")
  feature.def <- SummarizedExperiment::rowData(xcms.sum)%>%
    tibble::as_tibble()

  feature.matrix <- SummarizedExperiment::assay(xcms.sum)
  rsd <- function(x){sd(x,na.rm =  T)/mean(x , na.rm = T)}
  feature.matrix.qc <- feature.matrix[,which(grepl("QC",colnames(feature.matrix)))]
  feature.matrix.sample <- feature.matrix[,which(grepl("Sample",colnames(feature.matrix)))]
  feature.def$qc_rsd <- apply(feature.matrix.qc, 1, rsd)
  feature.def$sample_rsd <- apply(feature.matrix.sample, 1, rsd)
  feature.def$med_intensity <- apply(feature.matrix , 1 ,median)
  SummarizedExperiment::rowData(xcms.sum) <-feature.def
  return(xcms.sum)
}
