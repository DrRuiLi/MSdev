


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
  open_file(temp.file)

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
  message("Show color palette from ggsci: ",
          crayon::magenta(pal))
  scales::show_col(pal.col)

  return(pal.col)

}

colored_text <- function(x , color = "#E64B35"){

  crayon::make_style(color)(x)


}



#' export plot to pdf file
#' pdf() and export::graph2pdf() not support `append` arg
#' using qpdf::pdf_combine() to realize that function
#'
#' @param p
#' @param file_path
#' @param append
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
export_graph2pdf <- function(p ,
                             file_path ,
                             append = F,
                             ...
                             ){

  if("list" %in% class(p)){
    if((!append) & file.exists(file_path)){
      pdf(file_path)
      dev.off()
    }
    for(i in seq_along(p)){

      export_graph2pdf(p[[i]],file_path = file_path ,append = T,...)


    }
    return(invisible())


  }

  file_path <- normalizePath(file_path,mustWork = F)
  tpf1 <- paste0(tempfile(),".pdf")
  tpf2 <- paste0(tempfile(),".pdf")
  suppressMessages(export::graph2pdf(plot_multi_formate(p),
                                     file = tpf1,...))

  if (append & file.exists(file_path)) {
    qpdf::pdf_combine(input = c(file_path,tpf1),
                      output = tpf2)
  }else{
    tpf2 <- tpf1
  }
  message("Exported graph as ",file_path)
  file.copy(tpf2,file_path,overwrite = T)
  suppressWarnings(file.remove(tpf1,tpf2))
  return(invisible())

}

plot_multi_formate <- function(p){


  if ("Heatmap" %in% class(p)) {

    ComplexHeatmap::draw(p)

  }else{

    p
  }


}



ggplot_from_img <- function(img_path,coord = F,...){

#
  #p.jpg <- readJPEG(img_path,native = T)
  #p.frame <- ggplot()+
  #  xlim(c(0,10))+
  #  ylim(c(0,10))+
  #  theme_void()
  #p <- p.frame + inset_element(p = p.jpg,
  #                             position[1],
  #                             position[2],
  #                             position[3],
  #                             position[4]
  #                             )
  #return(p)


  p <- ggplot() +
    #xlim(c(0,10))+
    #ylim(c(0,10))+
    ggimage::geom_image(aes(x=5,y=5,image = img_path),...)+
    theme_void()
  if (coord) {
    p <- p+theme_classic()
  }

  return(p)

}



ggplot_km <- function(km.km,legend_tile = "group",
                      legend_label = names(km.km$strata) ){

  #km.km <- surv_fit(Surv(time, status) ~ 1, data = colon)
  km.sum <- summary(km.km)
  km.pval <- survminer::surv_pvalue(km.km)
  plot.data <- data.frame(
    time = km.sum$time,
    strata= km.sum$strata,
    surv=km.sum$surv,
    lower=km.sum$lower,
    upper=km.sum$upper
  )

  ggplot(plot.data)+
    geom_line(aes(x = time , y = surv,
                  col = strata))+
    geom_ribbon(aes(x = time ,ymin = lower,
                    ymax = upper,fill = strata),alpha = 0.3)+
    annotate(geom="text",x = quantile(range(plot.data$time),0.2),
             y = 0.2,
             label =paste0("p = ", format(km.pval[2],scientific = T,digits =2)))+
    scale_color_manual(values = ggsci::pal_npg()(8),
                       labels = legend_label)+
    scale_fill_manual(values = ggsci::pal_npg()(8),
                      labels = legend_label)+
    ylim(c(0,1))+
    labs(x = "Time",y = "Survival Probability",
         col = legend_tile,
         fill = legend_tile)+
    theme_bw()

}

