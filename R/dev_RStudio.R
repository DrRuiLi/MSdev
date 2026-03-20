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

#' Set RStudio Files Pane Directory
#'
#' @description Navigates the RStudio files pane to the specified directory or the directory of the currently active file.
#'
#' @param path Path to a file or directory. Defaults to the path of the currently active file in the RStudio editor.
#'
#' @return NULL (invisibly)
#' @export
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

