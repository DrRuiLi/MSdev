


#' @title open_ggplot_win
#' @description
#' create a temp.png file and open in Windows
#'
#'
#' @param p ggplot objective
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
open_ggplot_win <- function(p,width = NA,height = NA){

  temp.file <- tempfile(fileext = ".png")
  ggplot2::ggsave(filename = temp.file,plot = p,width = width,height= height)
  open_dir(temp.file)

}
