#' @title capture_base_plot
#' @description
#' Evaluate base graphics in a null device and return a \code{recordedplot}
#' object for replay or export (e.g. with \code{\link{open_plot_win}}).
#'
#' @param expr Base graphics commands, typically passed as a braced expression.
#' @param envir Environment in which \code{expr} is evaluated.
#'
#' @return A \code{recordedplot} object (from \code{\link[grDevices]{recordPlot}}).
#' @export
#'

capture_base_plot <- function(expr, envir = parent.frame()) {
  grDevices::pdf(nullfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  grDevices::dev.control(displaylist = "enable")
  eval(substitute(expr), envir = envir)
  grDevices::recordPlot()
}


#' @title open_plot_win
#' @description
#' create a temp.png file and open in Windows
#'
#'
#' @param p ggplot/Complexheatmap
#' @param width num
#' @param height num
#'
#' @return null
#' @export
#'

open_plot_win <- function(p,width = 5,height = 4){

  temp.file <- tempfile(fileext = ".png")
  if (any(c("Heatmap","HeatmapList")%in%class(p))) {
    export::graph2png(ComplexHeatmap::draw(p,padding = unit(c(1, 1, 1, 3), "mm")),
                     file =temp.file,
                     width = width,height= height
                     )
  }else if(any(c("ggplot")%in%class(p))){
    ggplot2::ggsave(filename = temp.file,plot = p,
                    width = width,height= height,dpi = 600)
  }else{
    export::graph2png(p,
                      file =temp.file,
                      width = width,height= height
    )
  }
  open_file(temp.file)

}


open_plot_pdf <- function(p,width = 5,height = 4){

  temp.file <- tempfile(fileext = ".pdf")
  if (any(c("Heatmap","HeatmapList")%in%class(p))) {
    export::graph2png(ComplexHeatmap::draw(p),
                      file =temp.file,
                      width = width,height= height
    )
  }else if(any(c("ggplot")%in%class(p))){
    ggplot2::ggsave(filename = temp.file,plot = p,
                    width = width,height= height,dpi = 600)
  }else{
    export::graph2pdf(p,
                      file =temp.file,
                      width = width,height= height
    )
  }
  open_file(temp.file)

}


open_plot_ppt <- function(p,width = 5,height = 4){

  temp.file <- tempfile(fileext = ".pptx")
  if (any(c("Heatmap","HeatmapList")%in%class(p))) {
    export::graph2ppt(ComplexHeatmap::draw(p,padding = unit(c(1, 1, 1, 3), "mm")),
                      file =temp.file,
                      width = width,height= height
    )
  }else if(any(c("ggplot")%in%class(p))){
    export::graph2ppt(p,
                      file =temp.file,
                      width = width,height= height
    )
  }else{
    export::graph2ppt(p,
                      file =temp.file,
                      width = width,height= height
    )
  }
  open_file(temp.file)

}



#' @title ggplot_sum_patchwork
#' @description
#' add all ggplot by patchwork
#'
#' @param ggplot.list a list with all item as ggplot objective
#'
#' @return null
#' @export
#'

ggplot_sum_patchwork <- function(ggplot.list){
  x <- ggplot.list
  x.len <- length(x)
  sum.exp <- 1
  x.exp <- paste0( "x.sum <- ",paste0(paste0("x[[",1:x.len,"]]"),collapse = " + "),
                   "+patchwork::plot_annotation(tag_levels=\"A\")")%>%
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
#' @title Export Graph2pdf
#' @description Graph2pdf.
#' @param p ggplot
#' @param file_path file path
#' @param append logic
#' @param ... additional arguments passed to export functions
#'
#' @return null
#' @export
#'

export_graph2pdf <- function(p ,
                             file_path ,
                             append = F,
                             ...
                             ){
  dir.exists(dirname(file_path))
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
    ggimage::geom_image(aes(x=0,y=0,image = img_path),...)+
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
             label =paste0("p = ", format(km.pval[2],scientific = F,digits =3)))+
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


colramp<- function(breaks = c(0,0.5,1),
                   colors = c("white","#F7844F","#B20C26"),
                   na.col = "#AAAAAA",
                   ...){

  circlize::colorRamp2(breaks =breaks,
                       colors = colors,
                       ...)
}


make_group_color <- function(x,palette = "random",verbose=F){

  palette <- match.arg(palette,c(names(ggsci:::ggsci_db),"random"))
  x <- unique(x)%>%sort()
  if (palette == "random") {
    col  <- randomcoloR::distinctColorPalette(length(x))
  }else{
    col <- ggsci:::ggsci_db[[palette]][[1]][1:length(x)]

  }
  if (verbose) {
    message(paste0(mapply(x,col,FUN = colored_text),collapse = " "))

  }
  names(col)<-x
  return(col)
}


scale_color_random <- function(...) {
  ggplot2::discrete_scale("colour", "random", function(n) randomcoloR::distinctColorPalette(n), ...)
}



scale_fill_random <- function(...) {
  ggplot2::discrete_scale("fill", "random", function(n) randomcoloR::distinctColorPalette(n), ...)
}

# Function to format legend labels with subscripts, with custom color and labels
scale_subscript_legend <- function(values = NULL, labels = NULL) {

  # Helper function to format labels with subscripts
  format_labels <- function(labels) {
    lapply(labels, function(label) {
      # Check if label contains a number
      if (grepl("\\d", label)) {
        # If number found, split into text and number
        split_label <- strsplit(label, "(?<=\\D)(?=\\d)", perl = TRUE)[[1]]
        # Create an expression with subscript
        expression_text <- bquote(.(split_label[1])[.(split_label[2])])
      } else {
        # If no number, return label as is
        expression_text <- label
      }
      return(expression_text)
    })
  }

  # Create a custom scale for color and fill
  scale_custom_labels <- list(
    scale_color_manual(values = values, labels = format_labels(labels)),
    scale_fill_manual(values = values, labels = format_labels(labels))
  )

  # Return the list of scale components to add to the plot
  return(scale_custom_labels)
}



heatmap_set_size <- function(hm,width  = 5,height= 5){

  hm@heatmap_param$width <- unit(width,"cm")
  hm@heatmap_param$height <- unit(height,"cm")

  return(hm)
}


get_ggplot_from_heatmap <- function(hm){
  wrap_elements(grid::grid.grabExpr(draw(hm)))
}

ggplot_irange <- function(IR, scale = 1e-6){

  ir.data <- as.data.frame(IR)
  ir.data <- (ir.data *  scale)%>%
    dplyr::mutate(no = 1:n())

  ggplot(ir.data)+
    geom_segment( aes(x = start,xend = end, y = no,yend = no) )

}
