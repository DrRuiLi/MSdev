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

    n.sig <- grep(paste0("^",i,"_significant"),names(data.aaa@elementMetadata@listData))


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
                               lfc = 0.5,
                               contrast = list_DEP_contrast(data.se )[1],
                               p.adjust = T
                               ){
  if (length(grep("_significant", colnames(rowData(data.se))))<1) {
    if (p.adjust) {
      data.se <- DEP::add_rejections(data.se)

    } else{
      data.se <- add_rejections_no_p.adj(data.se,lfc = lfc)
    }
  }

  diff.table.adj <- DEP::plot_volcano(data.se,
                                      contrast = contrast,
                                      plot = F,adjusted = T)
  diff.table <- DEP::plot_volcano(data.se,
                                  contrast =contrast,
                                  plot = F,adjusted = F)
  diff.table<- diff.table%>%
    dplyr::mutate(`adjusted_p_value_-log10` = diff.table.adj$`adjusted_p_value_-log10`,
                  .after = `p_value_-log10`)

  data.row <- rowData(data.se)%>%as.data.frame()
  diff.table<-diff.table%>%
    dplyr::mutate(data.row[protein,1:15])
  diff.table
}

DEP.plot.volcano <- function(data.se,
                             lfc = 0.5,
                             contrast = list_DEP_contrast(data.se )[1],
                             p.adjust = F,
                             show.label = T){
  if (length(grep("_significant", colnames(rowData(data.se))))<1) {
    if (p.adjust) {
      data.se <- DEP::add_rejections(data.se)

    } else{
      data.se <- add_rejections_no_p.adj(data.se,lfc = lfc)
    }
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
    labs(title = contrast ,
         x = expression(log[2] ~ "Fold change")) +
    labs(y = expression(-log[10] ~ "P-value"))+
    geom_vline(xintercept = c(-lfc,lfc),lty = 2,size = 0.2) +
    geom_hline(yintercept = c(-log10(0.05)),lty = 2,size = 0.2) +
    xlim(c(-x.max,x.max))+
    ylim(c(0,y.max))+
    coord_fixed(ratio = 2*x.max/y.max)+
    DEP::theme_DEP1() +
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
          plot.title = element_text(size = 6,hjust = 0),
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
                               segment.size = 0.1)
  }
  return(p)




}

DEP.plot.volcano.lipidomic <- function(data.se,
                                       contrast = list_DEP_contrast(data.se )[1],
                                       p.adjust = F,
                                       show.label = T){
  if (p.adjust) {
    data.se <- DEP::add_rejections(data.se)

  } else{
    data.se <- add_rejections_no_p.adj(data.se)
  }

  data("lipid.classification",package = "MSdb" )
  volcano.data <- DEP::plot_volcano(data.se,contrast,plot = F,adjusted = p.adjust )%>%
    dplyr::mutate(diff = case_when(log2_fold_change > 0 & significant~ "up",
                                   log2_fold_change < 0 & significant~ "down",
                                   T~ "no"),
                  SummarizedExperiment::rowData(data.se[protein,])%>%
                    as.data.frame()
                  )%>%
    dplyr::mutate(
      Lipid_class=lipid.classification$Lipid_class[match(Lipid_subclass,lipid.classification$Lipid_subclass)],
      col = case_when(log2_fold_change > 0 & significant~ Lipid_class,
                       log2_fold_change < 0 & significant~ Lipid_class,
                       T~ "no"),
      .after = Lipid_subclass
    )
  volcano.data$col[volcano.data$col=="no"]<-NA

  if (p.adjust) {
    volcano.data$p<-volcano.data$`adjusted_p_value_-log10`

  }else{
    volcano.data$p <- volcano.data$`p_value_-log10`
  }


  x.max <- max(abs(volcano.data$log2_fold_change))
  y.max <- max(volcano.data$p)

  col.list <-ggsci::pal_npg()(length(unique(volcano.data$col)))%>%
    `names<-`(unique(volcano.data$col))
  ggplot(volcano.data, aes(log2_fold_change, p)) +
    geom_point(aes(col = col),size = 0.3) +
    labs( x = expression(log[2] ~ "Fold change")) +
    labs(y = expression(-log[10] ~ "P-value"),
         col = "Lipid class")+
    geom_vline(xintercept = c(-1,1),lty = 2,size = 0.2) +
    geom_hline(yintercept = c(-log10(0.05)),lty = 2,size = 0.2) +
    scale_color_manual(values = col.list,breaks = na.omit(names(col.list)))+
    xlim(c(-x.max,x.max))+
    ylim(c(0,y.max))+
    coord_fixed(ratio = 2*x.max/y.max)+
    DEP::theme_DEP1() +
    geom_text(data = data.frame(),
              aes(x = c(-x.max*0.8, x.max*0.8),
                  y = c(y.max*0.9, y.max*0.9),
                  hjust = c(1, 0),
                  vjust = c(-1, -1),
                  label = c(sum(volcano.data$diff=="down"),
                            sum(volcano.data$diff=="up")),
                  fontface = "bold"),size = 2)+
    theme(legend.position = "right",
          legend.text =element_text(size = 4),
          legend.title =element_text(size = 6),
          legend.key.size = unit(0.1,"inch"),
          text = element_text(size = 4),
          panel.border = element_rect(size = 0.2),
          plot.title = element_text(size = 6),
          axis.text = element_text(size = 4),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6)
    ) ->p
  p
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

DEP.plot.lfc.lipid.class <- function(data.se,
                                     contrast = list_DEP_contrast(data.se )[1],
                                     p.adjust = F){
  if (p.adjust) {
    data.se <- DEP::add_rejections(data.se)

  } else{
    data.se <- add_rejections_no_p.adj(data.se)
  }

  data("lipid.classification",package = "MSdb" )
  plot.data <- DEP::plot_volcano(data.se,contrast,plot = F,adjusted = p.adjust )%>%
    dplyr::mutate(diff = case_when(log2_fold_change > 0 & significant~ "up",
                                   log2_fold_change < 0 & significant~ "down",
                                   T~ "no"),
                  SummarizedExperiment::rowData(data.se[protein,])%>%
                    as.data.frame()
    )

  if (p.adjust) {
    plot.data$p<-plot.data$`adjusted_p_value_-log10`

  }else{
    plot.data$p <- plot.data$`p_value_-log10`
  }

  plot.data <- plot.data%>%
    dplyr::mutate(
      Lipid_class=lipid.classification$Lipid_class[match(Lipid_subclass,lipid.classification$Lipid_subclass)],
      col = case_when(log2_fold_change > 0 & significant~ Lipid_class,
                      log2_fold_change < 0 & significant~ Lipid_class,
                      T~ "no"),
      size = case_when(significant ~ p,
                       T~ 1),
      .after = Lipid_subclass
    )

  plot.data$col[plot.data$col=="no"]<-NA

  plot.lfc.data <- plot.data[sample(1:nrow(plot.data),nrow(plot.data)),]%>%
    dplyr::arrange(Lipid_class,Lipid_subclass)%>%
    dplyr::group_by(Lipid_class)%>%
    dplyr::mutate(protein = factor(protein,levels = protein),
                  col.x = case_when(abs(log2_fold_change) > 1&p > 1.30103~ Lipid_class,
                                    T~"no"),
                  size.x = case_when(abs(log2_fold_change) > 1&p > 1.30103~ p,
                                     T~ 1),
                  x = cur_group_rows()+cur_group_id()*100-1)%>%
    dplyr::ungroup()

  x.labels <- plot.lfc.data%>%
    dplyr::group_by(Lipid_class)%>%
    dplyr::mutate(x.mean = mean(x))%>%
    dplyr::distinct(Lipid_class,x.mean)

  col.list <-ggsci::pal_npg()(length(unique(plot.data$col)))%>%
    `names<-`(unique(plot.data$col))
  ggplot(plot.lfc.data)+
    geom_hline(yintercept = c(-1,1),lty = "dashed",col = "black")+
    geom_point(aes(x = x , y = log2_fold_change ,
                   col = col.x ,size= size.x),
               alpha = 0.6)+
    scale_x_continuous(breaks = x.labels$x.mean,
                       labels = x.labels$Lipid_class)+
    scale_size(range = c(0,3))+
    scale_color_manual(values = col.list)+
    labs(x = NULL, y = expression(log[2] ~ "Fold change"),
         size =  expression(-log[10] ~ "P-value"),
         col = "Lipid Class")+
    DEP::theme_DEP1()+
    theme(legend.position = "right",
          text = element_text(size = 6),
          plot.title = element_text(size = 6),
          strip.text = element_text(size = 5,
                                    margin = margin(0,0.03,0,0.03,"inch")),
          strip.background = element_rect(fill = "grey",
                                          size = 0.2
          ),
          axis.text = element_text(size = 4),
          axis.ticks = element_line(size =0.2),
          axis.ticks.length = unit(0.01,"inch"),
          axis.title.x = element_text(size = 4),
          axis.title.y = element_text(size = 4),
          legend.title = element_text(size = 4),
          legend.text = element_text(size = 4),
          legend.key.height  = unit(0.1,"inch"),
          legend.key.width   = unit(0.05,"inch"),
          panel.border = element_rect(size = 0.2),
          panel.grid =element_blank()
    )->p
  p

  return(p)




}

DEP.plot.heatmap <- function(data.se,
                             contrast = list_DEP_contrast(data.se )[1],
                             p.adjust = F){

  if (length(grep("_significant", colnames(rowData(data.se))))<1) {
    if (p.adjust) {
      data.se <- DEP::add_rejections(data.se)

    } else{
      data.se <- add_rejections_no_p.adj(data.se,lfc = lfc)
    }
  }


  col.info <-SummarizedExperiment::colData(data.se)%>%
    as.data.frame()%>%
    dplyr::mutate(col.group = condition,
                  col.label = sample.label)
  row.info <-DEP::plot_volcano(data.se,contrast,plot = F,adjusted = p.adjust )%>%
    dplyr::mutate(SummarizedExperiment::rowData(data.se[protein,])%>%
                    as.data.frame()
    )%>%
    dplyr::mutate(row.group = "",
                  row.label = Compound_name)%>%
    dplyr::filter(significant)

  heatmap.matrix <-SummarizedExperiment::assay(data.se[row.info$ID,col.info$ID])%>%
    `^`(2,.)%>%t%>%scale%>%t



  plotHeatmap(heatmap.matrix,col.info,row.info)



}



