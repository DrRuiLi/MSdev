


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
  ggplot2::ggsave(filename = temp.file,plot = p,
                  width = width,height= height,dpi = 600)
  open_dir(temp.file)

}




#' @title ggplot_sum_patchwork
#' @description
#' add all ggplot by patchwork
#'
#' @param ggplot.list a list with all item as ggplot objective
#'
#' @return
#' @export
#' @import patchwork
#'
#' @examples
ggplot_sum_patchwork <- function(ggplot.list){
  x <- ggplot.list
  x.len <- length(x)
  sum.exp <- 1
  x.exp <- paste0( "x.sum <- ",paste0(paste0("x[[",1:x.len,"]]"),collapse = " + "),
                   "+plot_annotation(tag_levels=\"A\")")%>%
    str2expression()
  eval(x.exp)
  return(x.sum)
}



ggplot_ggsci <- function(pal = "npg"){

  ggsci_db <- ggsci:::ggsci_db
  ggsci_pal <- names(ggsci_db)
  pal <- match.arg(pal,ggsci_pal)
  pal.col <- ggsci_db[[pal]][[1]]
  message("Show color palette from ggsci: ", crayon::magenta(pal))
  scales::show_col(pal.col)

  return(pal.col)

}

colored_text <- function(x , color = "#E64B35"){

  crayon::make_style(color)(x)


}




