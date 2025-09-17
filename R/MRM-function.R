#' export_MChromatograms_Metabolites
#' @describeIn MRM_data_analysis MRM data analysis
#' @param mchroms XChromatograms or MChromatograms
#' @param file tempfile(fileext = ".pdf")
#'
#' @returns tempfile(fileext = ".pdf")
#' @export
#'
export_MChromatograms_Metabolites <- function(mchroms,file = tempfile(fileext = ".pdf"),patched = F){

  fda <- fData(mchroms)%>%
    dplyr::mutate(row = 1:n(),
                  label =sub('.*name=([^ ]+).*', '\\1', chromatogramId ),
                  temp_idf = paste0(polarity,"_",precursorIsolationWindowTargetMZ)
                  #temp_idf = chromatogramId
                  )%>%
    dplyr::arrange(label)
  idfs <- unique(fda$temp_idf)
  pl <- list()
  for (idf in idfs) {
    idf.fda <- fda%>%
      dplyr::filter(temp_idf == idf)%>%
      dplyr::mutate(color_f = paste0(label,": ",productIsolationWindowTargetMZ) )


    p <- plot_XChromatograms(mchroms[idf.fda$row,],norm = F,move = F,color_by = "row",
                        color_f =idf.fda$color_f)+
      labs(subtitle = paste0("mz = ",unique(idf.fda$precursorIsolationWindowTargetMZ)))
    # labs(subtitle  = paste0(idf.fda$chromatogramId,collapse = "\n"))
    pl[[idf]] <- p
    if (!patched)
      export_graph2pdf(p,file_path = file,append = T)

  }
  if (patched) {
    p <- ggplot_sum_patchwork(pl)+
      plot_layout(ncol = 2)
    export_graph2pdf(p,file_path = file,append = F,width = 12, height = 15)
  }

  return(file)

}
