#' @title setwdActivedFilePath
#' @description see `rstudioapi::getSourceEditorContext()`
#' @return
#' @export
#'
setwdActivedFilePath <- function(){

  path <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(path)

}
