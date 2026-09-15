#' @title Load bundled demo data objects
#' @description Load pre-built demo objects from package **MSdemo**.
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
  if (!requireNamespace("MSdemo", quietly = TRUE)) {
    stop(
      "Demo data now lives in package MSdemo. Install MSdemo and try again.",
      call. = FALSE
    )
  }
  MSdemo::load_demo(demo = demo)
}

#' @title Rebuild the MSdev demo project
#' @description Wrapper for `MSdemo::make_demo()`. Requires package **MSdemo**
#'   (raw SCIEX files) plus **MSconvertR**.
#'
#' @param ... Passed to `MSdemo::make_demo()`.
#' @return The processed `MSdev` object (invisibly).
make_demo <- function(...) {
  if (!requireNamespace("MSdemo", quietly = TRUE)) {
    stop(
      "Demo rebuild now lives in package MSdemo. Install MSdemo and try again.",
      call. = FALSE
    )
  }
  MSdemo::make_demo(...)
}
