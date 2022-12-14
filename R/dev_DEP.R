#' @title list_DEP_contrast
#' @description return contrasts in data.se from DEP::test_diff
#' @param data.se
#'
#' @return
#'
#' @examples
list_DEP_contrast <- function(data.se){

  contrast <- grep(names(data.se@elementMetadata@listData),
                   pattern = "vs.+diff",value = T)%>%
    sub(pattern = "_diff",
        replacement = "")
  contrast
}



#' @title add_rejections_no_p.adj
#' @description reference to DEP::add_rejections(),which not support significant with out p adjust,
#' this function as supplymentary
#' @param data.se
#' @param alpha
#' @param lfc
#'
#' @return
#'
#' @examples
add_rejections_no_p.adj <- function(data.se , alpha = 0.05,lfc = 1){

  data.aaa <- DEP::add_rejections(data.se , alpha ,lfc )
  sum(data.aaa@elementMetadata@listData$significant)
  contrasts <-list_DEP_contrast(data.aaa )

  for (i in contrasts) {

    n.sig <- grep(paste0(i,"_significant"),names(data.aaa@elementMetadata@listData))


    data.aaa@elementMetadata@listData[["significant"]]<-
      data.aaa@elementMetadata@listData[[n.sig]] <-
        (data.aaa@elementMetadata@listData[[paste0(i,"_p.val")]] < alpha&
           abs(data.aaa@elementMetadata@listData[[paste0(i,"_diff")]]) > lfc)

  }
  data.aaa

}



#' @title DEP.test.diff
#' @description warpper of DEP::test_diff
#' @param data.se
#' @param p.adj
#'
#' @return
#' @export
#'
#' @examples
DEP.test.diff <- function(data.se){

  groups <- data.se$condition%>%
    groupStringFactor()
  data.diff<- DEP::test_diff(data.se,type = "all",control = levels(groups)[1])
  #data.diff<- add_rejections(data.diff,alpha = 0.05,lfc =1)
#
  #if (!p.adj) {
#
  #  data.diff <- add_rejections_no_p.adj(data.diff)
  #  #sum(data.diff@elementMetadata@listData$significant)
#
  #}
  data.diff
}

DEP.get.diff.table <- function(data.se,
                               contrast = list_DEP_contrast(data.se )[1],
                               p.adujst = T
                               ){
  if (p.adujst) {
    data.se <- DEP::add_rejections(data.se)

  } else{
    data.se <- add_rejections_no_p.adj(data.se)
  }

  diff.table.adj <- DEP::plot_volcano(data.se,
                                      contrast = list_DEP_contrast(data.se )[1],
                                      plot = F,adjusted = T)
  diff.table <- DEP::plot_volcano(data.se,
                                  contrast = list_DEP_contrast(data.se )[1],
                                  plot = F,adjusted = F)
  diff.table<- diff.table%>%
    dplyr::mutate(`adjusted_p_value_-log10` = diff.table.adj$`adjusted_p_value_-log10`,
                  .after = `p_value_-log10`)

  diff.table
}

DEP.plot.volcano <- function(data.se,
                             contrast = list_DEP_contrast(data.se )[1],
                             p.adjust = F,
                             show.label = T){
  if (p.adujst) {
    data.se <- DEP::add_rejections(data.se)

  } else{
    data.se <- add_rejections_no_p.adj(data.se)
  }

  volcano.data <- DEP::plot_volcano(data.se,contrast,plot = F,adjusted = p.adjust )%>%
    dplyr::mutate(diff = case_when(log2_fold_change > 0 & significant~ "up",
                                   log2_fold_change < 0 & significant~ "down",
                                   T~ "no"))

  if (p.adjust) {
    volcano.data$p<-volcano.data$`adjusted_p_value_-log10`

  }else{
    volcano.data$p <- volcano.data$`p_value_-log10`
  }


  x.max <- max(abs(volcano.data$log2_fold_change))
  y.max <- max(volcano.data$p)

  ggplot(volcano.data, aes(log2_fold_change, p)) +
    geom_point(aes(col = diff),size = 0.1) +
    labs( x = expression(log[2] ~ "Fold change")) +
    labs(y = expression(-log[10] ~ "P-value"))+
    geom_vline(xintercept = c(-1,1),lty = 2,size = 0.2) +
    geom_hline(yintercept = c(-log10(0.05)),lty = 2,size = 0.2) +
    xlim(c(-x.max,x.max))+
    ylim(c(0,y.max))+
    coord_fixed(ratio = 2*x.max/y.max)+
    theme_DEP1() +
    scale_color_manual(values = c("up" = "#DD001B",
                                  "down" = "#2A77D3",
                                  "no" = "grey"))+
    geom_text(data = data.frame(),
              aes(x = c(-x.max*0.8, x.max*0.8),
                  y = c(y.max*0.9, y.max*0.9),
                  hjust = c(1, 0),
                  vjust = c(-1, -1),
                  label = c(sum(volcano.data$diff=="down"),
                            sum(volcano.data$diff=="up")),
                  fontface = "bold"),size = 2)+
    theme(legend.position = "none",
          text = element_text(size = 4),
          panel.border = element_rect(size = 0.2),
          plot.title = element_text(size = 6),
          axis.text = element_text(size = 4),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6)
    ) ->p

  if (show.label) {
    row.data <- SummarizedExperiment::rowData(data.se)%>%
      as.data.frame()
    volcano.data <- volcano.data%>%
      dplyr::mutate(Compound_name = row.data$Compound_name[match(protein , row.data$feature_id)])

    p <- p+
      ggrepel::geom_text_repel(data = filter(volcano.data, significant),
                               aes(label = Compound_name), size = 1,
                               box.padding = unit(0.1, "lines"),
                               point.padding = unit(0.1, "lines"),
                               segment.size = 0.5)
  }
  return(p)




}

DEP.plot.volcano.lipidomic <- function(){





}





