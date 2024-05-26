#' @title setwdActivedFilePath
#' @description Set Working Directory as currently opened file in RStudio editor,
#' see `rstudioapi::getSourceEditorContext()`
#' @return null
#' @export
#'
setwdActivedFilePath <- function(){

  path <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(path)

}


#' setRStudioDir
#'
#' @param path path
#'
#' @return null
#' @export
#'
setRStudioDir <- function(path = rstudioapi::getSourceEditorContext()$path){

  rstudioapi::filesPaneNavigate(path)

}
getRStudioDir <- function(){

  path = rstudioapi::getSourceEditorContext()$path
  return(path)
}


getActivedFilePath <- function(){
  rstudioapi::getSourceEditorContext()$path

}

