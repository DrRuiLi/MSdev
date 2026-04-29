

setClass(Class = "MSIPIsotopologueData",
         contains = "SummarizedExperiment")


MSIPIsotopologueData <- function(assays,
                                 rowData = S4Vectors::DataFrame(),
                                 colData = S4Vectors::DataFrame(),
                                 metadata = list()) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = rowData,
    colData = colData,
    metadata = metadata
  )
  methods::as(se, "MSIPIsotopologueData")
}


setMethod(f = "show", signature = "MSIPIsotopologueData",
          definition = function(object) {
            rda <- SummarizedExperiment::rowData(object)
            if (!nrow(rda)) {
              cat("MSIPIsotopologueData (empty)\n")
              return(invisible())
            }
            iso_form <- NULL
            if ("isotopologue_form" %in% colnames(rda)) {
              iso_form <- as.character(rda$isotopologue_form)
            } else {
              iso_form <- rownames(rda)
            }
            iso_form <- iso_form[!is.na(iso_form) & nzchar(iso_form)]
            if (!length(iso_form)) {
              cat("MSIPIsotopologueData (", nrow(rda), " rows)\n", sep = "")
              return(invisible())
            }
            l <- paste0(utils::head(iso_form, 10), collapse = "; ")
            if (length(iso_form) > 10) l <- paste0(l, "; ...")
            cat(length(iso_form), " isotopologues: ", l, "\n", sep = "")
          })


setClass(Class = "MSIPMetaboliteData",
         slots = list(
           "CompoundInfo" = "list",
           "Spectra" = "list",
           # legacy container; kept flexible for backward compatibility
           "MSIPIsotopologueDatas" = "ANY"
         ))

MSIPMetaboliteData <- function(CompoundInfo = list(),
                               Spectra = list(),
                               MSIPIsotopologueDatas = MSIPIsotopologueData()
){

  new("MSIPMetaboliteData",
      CompoundInfo = CompoundInfo,
      Spectra = Spectra,
      MSIPIsotopologueDatas = MSIPIsotopologueDatas
  )


}

setMethod(f = "show",signature = "MSIPMetaboliteData",definition = function(object){

  cat("MSIPMetaboliteData")
})
