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
#' @title Setrstudiodir
#' @description SetRStudioDir.
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


get_dir_ActivedFilePath <- function(){
  rstudioapi::getSourceEditorContext()$path

}

