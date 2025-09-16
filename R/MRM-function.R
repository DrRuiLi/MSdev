#' export_MChromatograms_Metabolites
#' @describeIn MRM_data_analysis MRM data analysis
#' @param mchroms XChromatograms or MChromatograms
#' @param file tempfile(fileext = ".pdf")
#'
#' @returns tempfile(fileext = ".pdf")
#' @export
#'
export_MChromatograms_Metabolites <- function(mchroms,file = tempfile(fileext = ".pdf")){

  fda <- fData(mchroms)%>%
    dplyr::mutate(row = 1:n(),
                  label =sub('.*name=([^ ]+).*', '\\1', chromatogramId ),
                  temp_idf = paste0(polarity,"_",precursorIsolationWindowTargetMZ)
                  )
  idfs <- unique(fda$temp_idf)
  for (idf in idfs) {
    idf.fda <- fda%>%
      dplyr::filter(temp_idf == idf)%>%
      dplyr::mutate(color_f = paste0(label,": ",productIsolationWindowTargetMZ) )


    p <- plot_XChromatograms(mchroms[idf.fda$row,],norm = F,move = F,color_by = "row",
                        color_f =idf.fda$color_f)+
      labs(subtitle  = paste0(idf.fda$chromatogramId,collapse = "\n"))
    export_graph2pdf(p,file_path = file,append = T)
  }

  return(file)

}
