#' @title DEP_list_contrast
#' @description return contrasts in data.se from DEP::test_diff
#' @param data.se
#'
#' @return
#'
#' @examples
DEP_list_contrast <- function(data.se){

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
DEP_add_rejections <- function(data.se , p.adjust = T,lfc = 0.5){

  if (any(grepl(pattern = "significant",names(data.se@elementMetadata@listData)))) {
    message("significant existed, skip")
    data.diff <- data.se
  }else{
    data.diff <- DEP::add_rejections(data.se , alpha = 0.05 ,lfc )
  }

  data.diff <- DEP_p_adjust(data.diff)
  if (p.adjust) {
    contrasts <-DEP_list_contrast(data.diff )
    for (i in contrasts) {
      n.sig <- grep(paste0("^",i,"_significant$"),names(data.diff@elementMetadata@listData))
      data.diff@elementMetadata@listData[["significant"]]<-
        data.diff@elementMetadata@listData[[n.sig]] <-
        (data.diff@elementMetadata@listData[[paste0(i,"_p.adj")]] < 0.05 &
           abs(data.diff@elementMetadata@listData[[paste0(i,"_diff")]]) > lfc)

    }


  }else{
    contrasts <-DEP_list_contrast(data.diff )
    for (i in contrasts) {
      n.sig <- grep(paste0("^",i,"_significant$"),names(data.diff@elementMetadata@listData))
      data.diff@elementMetadata@listData[["significant"]]<-
        data.diff@elementMetadata@listData[[n.sig]] <-
        (data.diff@elementMetadata@listData[[paste0(i,"_p.val")]] < 0.05 &
           abs(data.diff@elementMetadata@listData[[paste0(i,"_diff")]]) > lfc)

    }

  }


  data.diff

}

DEP_p_adjust <- function(data.se , p.adjust.method = "fdr"){


  for (i in DEP_list_contrast(data.se)) {

    p.val <- data.se@elementMetadata[[paste0(i,"_p.val")]]
    p.adj <- data.se@elementMetadata[[paste0(i,"_p.adj")]]
    data.se@elementMetadata[[paste0(i,"_p.adj")]] <- p.adjust(p.val ,
                                                              method =p.adjust.method )

  }
  return(data.se)

}

DEP_check_sig <- function(data.se){

  if (length(grep("_significant", colnames(rowData(data.se))))<1) {

    stop("non signifcant, please run DEP_add_rejections")
  }


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
DEP_test_diff <- function(data.se,type = "all",...){

  groups <- data.se$condition%>%
    groupStringFactor()
  data.diff<- DEP::test_diff(data.se,
                             #control = levels(groups)[1],
                             ...)
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

DEP_get_diff_table <- function(data.se,
                               contrast = DEP_list_contrast(data.se )[1],
                               keep.all =F
                               ){
  DEP_check_sig(data.se)


  if (contrast == "all") {
    table.list <- list()
    for (i in 1:length( DEP_list_contrast(data.se ))) {
      i_contrast <-  DEP_list_contrast(data.se )[i]
      table.list[[i_contrast]] <- DEP_get_diff_table(data.se ,
                                                   contrast = i_contrast,
                                                   keep.all = keep.all)
    }
    return(table.list)
  }


  diff.table.adj <- DEP::plot_volcano(data.se,
                                      contrast = contrast,
                                      plot = F,adjusted = T)
  diff.table <- DEP::plot_volcano(data.se,
                                  contrast =contrast,
                                  plot = F,adjusted = F)
  row.data <- rowData(data.se)%>%as.data.frame()
  diff.table<- diff.table%>%
    dplyr::mutate(`adjusted_p_value_-log10` = diff.table.adj$`adjusted_p_value_-log10`,
                  .after = `p_value_-log10`)
  if (keep.all) {
    diff.table<- diff.table%>%
      dplyr::mutate(
        row.data[match(protein ,row.data$name),
                 !colnames(row.data)%in% "significant"])

  }

  diff.table
}

DEP_plot_volcano <- function(data.se,
                             contrast = DEP_list_contrast(data.se )[1],
                             show.label = T){
  if (length(grep("_significant", colnames(rowData(data.se))))<1) {

    stop("non signifcant, please run DEP_add_rejections")
  }

  if (contrast == "all") {
    p.vol.list <- list()
    for (i in 1:length( DEP_list_contrast(data.se ))) {
      i_contrast <-  DEP_list_contrast(data.se )[i]
      p.vol.list[[i_contrast]] <- DEP_plot_volcano(data.se ,
                                                   contrast = i_contrast,
                                                   show.label = show.label)
    }
    return(p.vol.list)
  }


  volcano.data <- DEP_get_diff_table(data.se,contrast = contrast)

  volcano.sig <- volcano.data%>%
    dplyr::filter(significant)
  lfc <- min(abs(volcano.sig$log2_fold_change))
  lfc <- ifelse(is.infinite(lfc),0,lfc)
  p.sig <- volcano.data[(abs(volcano.data$log2_fold_change)>lfc  )&
                         (volcano.data$`p_value_-log10` > -log10(0.05)),  ]

  p.adjust <- ifelse(any(!p.sig$protein%in% volcano.sig$protein),
                     T,F)

  if (p.adjust) {
    volcano.data$y <-volcano.data$`adjusted_p_value_-log10`

  }else{
    volcano.data$y <- volcano.data$`p_value_-log10`
  }

  volcano.data <- volcano.data %>%
    dplyr::mutate(diff = case_when(log2_fold_change > 0 & significant~ "up",
                                   log2_fold_change < 0 & significant~ "down",
                                   T~ "no") )



  x.max <- max(abs(volcano.data$log2_fold_change))
  y.max <- max(volcano.data$y)
  if (y.max==0 ) y.max <- 1

  ggplot(volcano.data, aes(log2_fold_change, y)) +
    geom_point(aes(col = diff),size = 0.1) +
    labs(title = contrast ,
         x = expression(log[2] ~ "Fold change")) +
    labs(y = ifelse(p.adjust,expression(-log[10] ~ "Adjusted P"),
                    expression(-log[10] ~ "P-value")))+
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
      dplyr::mutate(label = row.data$label[match(protein , row.data$feature_id)])

    p <- p+
      ggrepel::geom_text_repel(data = filter(volcano.data, significant),
                               aes(label = label), size = 1,
                               box.padding = unit(0.1, "lines"),
                               point.padding = unit(0.1, "lines"),
                               segment.size = 0.1)
  }
  return(p)




}

DEP.plot.volcano.lipidomic <- function(data.se,
                                       contrast = DEP_list_contrast(data.se )[1],
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
      dplyr::mutate(name = row.data$name[match(protein , row.data$feature_id)])

    p <- p+
      ggrepel::geom_text_repel(data = filter(volcano.data, significant),
                               aes(label = name), size = 1,
                               box.padding = unit(0.1, "lines"),
                               point.padding = unit(0.1, "lines"),
                               segment.size = 0.5)
  }
  return(p)




}

DEP.plot.lfc.lipid.class <- function(data.se,
                                     contrast = DEP_list_contrast(data.se )[1],
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
                             contrast = DEP_list_contrast(data.se )[1],
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
                  row.label = name)%>%
    dplyr::filter(significant)

  heatmap.matrix <-SummarizedExperiment::assay(data.se[row.info$ID,col.info$ID])%>%
    `^`(2,.)%>%t%>%scale%>%t



  plotHeatmap(heatmap.matrix,col.info,row.info)



}


#' DEP_export_data
#' wirte coldata and rowdata to excel
#'
#' @param data.se
#' @param file_path
#' @import SummarizedExperiment
#' @return
#' @export
#'
#' @examples
DEP_export_data <- function(data.se,file_path){

  col.data <- colData(data.se)%>%as.data.frame()
  row.data <- rowData(data.se)%>%as.data.frame()%>%
    cbind(assay(data.se))

  data.to.write <-list(sample_info = col.data,
                       compound = row.data)

  xlsx.write.list(data.to.write,file = file_path)
  return(invisible())

}


DEP_plot_PCA <- function(data.se,
                         col.group =NULL,
                         showlabel = F){

  se.coldata <- colData(data.se)%>%
    as.data.frame()
  pca.matrix <- assay(data.se)%>%
    `colnames<-`(se.coldata$label)%>%t

  if (is.null(col.group)) {
    col.group <- ggsci::pal_npg()(length(unique(se.coldata$condition)))
  }
 plotPCA(pca.matrix,
                  pca.group = se.coldata$condition,
                  showlabel = showlabel)+
    scale_color_manual(values = col.group)+
    scale_fill_manual(values = col.group)->p.pca

return(p.pca)

}

DEP_pathway_enrich <- function(data.se,
                               contrast ,
                               method = c("HyperTest","GlobalTest")){


  if (contrast == "all") {
    enrich.list <- list()
    for (i in 1:length( DEP_list_contrast(data.se ))) {
      i_contrast <-  DEP_list_contrast(data.se )[i]
      enrich.list[[i_contrast]] <- DEP_pathway_enrich(data.se ,
                                                   contrast = i_contrast,
                                                   method = method)
    }
    return(enrich.list)
  }
  method <- match.arg(method)
  if (method == "HyperTest") {
    pathway.table <- DEP_get_diff_table(data.se,
                                        contrast = contrast,
                                        keep.all = T)%>%
        dplyr::filter(significant)%>%
        dplyr::pull(kegg.id)%>%
      analyzePathwayHyperTest()

  }
  if (method == "GlobalTest") {

    data.se <- data.se[!is.na(rowData(data.se)$kegg.id),
                       data.se$condition%in% strsplit(contrast,"_vs_")[[1]]]
    pathway.matrix <- assay(data.se)
    rownames(pathway.matrix) <- rowData(data.se)$kegg.id
    pathway.table <-  analyzePathwayGlobalTest(t(pathway.matrix) ,data.se$condition)


  }

  return(pathway.table)



}

DEP_normalization <- function(data.se){

  data_filt <- DEP::filter_missval(data.se,
                                   thr = min(table(data.se$group))*0.3)
  data_norm <- DEP::normalize_vsn(data_filt)
  data_imp <- DEP::impute(data_norm, fun = "MinProb")
  data_imp
}
