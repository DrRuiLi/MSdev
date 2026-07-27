#' @title Load bundled demo data objects
#' @description Load pre-built demo objects from the MSdev demo project directory.
#' Dispatches to [MSdev_load()] for `"MSdev"` and `readRDS()` for other demo types.
#' `"xcms"` is an alias for `"XcmsExperiment"` (preferred). Use `"XCMSnExp"` for the
#' legacy on-disk `XCMSnExp` demo.
#'
#' @param demo Character. One of `"MSdev"`, `"XcmsExperiment"`, `"xcms"`,
#'   `"XCMSnExp"`, `"SummarizedExperiment"`, `"data.se"`, `"Spectra"`, or `"sp"`.
#'
#' @return The loaded demo object (type depends on `demo`).
#'
#' @examples
#' \dontrun{
#' msdev.demo <- load_demo("MSdev")
#' xcms.demo <- load_demo("xcms")
#' xe.demo <- load_demo("XcmsExperiment")
#' }
#'
#' @export
load_demo <- function(demo = c("MSdev",
                               "XcmsExperiment", "xcms",
                               "XCMSnExp",
                               "SummarizedExperiment", "data.se",
                               "Spectra", "sp")) {

  demo <- match.arg(demo)
  demo_dir <- "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo"
  file_path <- switch(
    demo,
    "MSdev" = file.path(demo_dir, "MSdev_2025_04_01.Rdata"),
    "XcmsExperiment" = file.path(demo_dir, "XcmsExperiment_2023_11_17.rda"),
    "xcms" = file.path(demo_dir, "XcmsExperiment_2023_11_17.rda"),
    "XCMSnExp" = file.path(demo_dir, "XCMSnExp_2023_11_17.rda"),
    "SummarizedExperiment" = file.path(demo_dir, "SummarizedExperiment_2023_11_03.rda"),
    "data.se" = file.path(demo_dir, "SummarizedExperiment_2023_11_03.rda"),
    "Spectra" = file.path(demo_dir, "Spectra_2023_11_23.rda"),
    "sp" = file.path(demo_dir, "Spectra_2023_11_23.rda")
  )

  fun <- switch(
    demo,
    "MSdev" = MSdev_load,
    "XcmsExperiment" = readRDS,
    "xcms" = readRDS,
    "XCMSnExp" = readRDS,
    "SummarizedExperiment" = readRDS,
    "data.se" = readRDS,
    "Spectra" = readRDS,
    "sp" = readRDS
  )
  fun(file_path)
}

make_demo <- function(){


  msdev.demo <- MSdev("C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/raw.data/")
  msdev.demo <- MSdev_msConvert(msdev.demo)
  msdev.demo <- MSdev_xcmsProcessing(msdev.demo)
  msdev.demo <- MSdev_match_Spectra_to_feature(msdev.demo)
  msdev.demo <- MSdev_annotation(msdev.demo,
                                 cpdb_path ="c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb" )
  MSdev_save(msdev.demo)


}


