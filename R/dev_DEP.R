#' @describeIn DEP_Style_se get DEP style \code{\link[SummarizedExperiment]{SummarizedExperiment}} from MSdev
#' @title \code{\link[DEP]{DEP}} styled \code{\link[SummarizedExperiment]{SummarizedExperiment}} and related analysis
#' @description DEP styled SummarizedExperiment and related analysis
#' @export
get_MSdev_DEP_se <- function(object,
                             from = c("metabolite.se",
                                      "feature.se"),
                             preprocess = T,...){

  from <- match.arg(from)
  data.se <- object@statData[[from]]

  ### format to DEP
  {


    sampleinfo <- object@sampleInfo
    ### col
    {
      cda <- colData(data.se)%>%
        as.data.frame()%>%
        dplyr::mutate(group = sampleinfo$group[match(sample.name,sampleinfo$sample.name)],
                      group = groupStringFactor(group),
                      condition = group,
                      sample.labels = sampleinfo$sample.labels[match(sample.name,sampleinfo$sample.name)],
                      label =sample.labels)%>%
        dplyr::group_by(condition)%>%
        dplyr::mutate(replicate = 1:n(),
                      ID = paste0(condition,num2str(1:n())))%>%
        dplyr::ungroup()%>%
        as.data.frame()
      rownames(cda) <- cda$ID



      colData(data.se) <- cda%>%S4Vectors::DataFrame()
    }

    ### row
    {
      rda <- rowData(data.se)%>%
        as.data.frame()%>%
        dplyr::mutate( label = name,
                       name = feature_id,
                       ID= feature_id)
      rowData(data.se) <- rda%>%S4Vectors::DataFrame()
    }

    assay(data.se) <- log2(assay(data.se))
  }


  ### pre process
  {
    data.se <- DEP_get_QC_RSD(data.se)
    data.se <- DEP_get_group_color(data.se)
    if (preprocess) {
      data.se <- DEP_preprocess(data.se,...)

    }

  }

  return(data.se)
}





#' @describeIn DEP_Style_se list all contrast in SummarizedExperiment
#' @export
DEP_list_contrast <- function(data.se){

  contrast <- grep(names(data.se@elementMetadata@listData),
                   pattern = "vs.+diff",value = T)%>%
    sub(pattern = "_diff",
        replacement = "")
  contrast
}





#' @param data.se SE
#' @param alpha p threshold
#' @param lfc fc threshold
#'
#' @return se
#'
#' @describeIn DEP_Style_se Add significant,Reference to \code{\link[DEP]{add_rejections}},which not support significant with out p adjust,
#' this function as supplymentary
#' @export
DEP_add_rejections <- function(data.se , p.adjust = T, p = 0.05,lfc = 0.5){

  if (any(grepl(pattern = "significant",names(data.se@elementMetadata@listData)))) {
    message("significant existed, skip")
    data.diff <- data.se
  }else{
    data.diff <- DEP::add_rejections(data.se , alpha = p ,lfc )
  }

  data.diff <- DEP_p_adjust(data.diff)
  if (p.adjust) {
    contrasts <-DEP_list_contrast(data.diff )
    for (i in contrasts) {
      n.sig <- grep(paste0("^",i,"_significant$"),names(data.diff@elementMetadata@listData))
      data.diff@elementMetadata@listData[["significant"]]<-
        data.diff@elementMetadata@listData[[n.sig]] <-
        (data.diff@elementMetadata@listData[[paste0(i,"_p.adj")]] < p  &
           abs(data.diff@elementMetadata@listData[[paste0(i,"_diff")]]) > lfc)

    }


  }else{
    contrasts <-DEP_list_contrast(data.diff )
    for (i in contrasts) {
      n.sig <- grep(paste0("^",i,"_significant$"),names(data.diff@elementMetadata@listData))
      data.diff@elementMetadata@listData[["significant"]]<-
        data.diff@elementMetadata@listData[[n.sig]] <-
        (data.diff@elementMetadata@listData[[paste0(i,"_p.val")]] < p  &
           abs(data.diff@elementMetadata@listData[[paste0(i,"_diff")]]) > lfc)

    }

  }
  data.diff@metadata$DEP_add_rejections$p <- p
  data.diff@metadata$DEP_add_rejections$p.adjust <- p.adjust
  data.diff@metadata$DEP_add_rejections$lfc <- lfc

  data.diff

}

#' @describeIn DEP_Style_se multiple test
#' @export
DEP_p_adjust <- function(data.se , p.adjust.method = "fdr"){


  for (i in DEP_list_contrast(data.se)) {

    p.val <- data.se@elementMetadata[[paste0(i,"_p.val")]]
    p.adj <- data.se@elementMetadata[[paste0(i,"_p.adj")]]
    data.se@elementMetadata[[paste0(i,"_p.adj")]] <- p.adjust(p.val ,
                                                              method =p.adjust.method )

  }
  return(data.se)

}

#' @describeIn DEP_Style_se check if \code{\link[MSdev]{DEP_add_rejections}} performed
DEP_check_sig <- function(data.se){

  if (length(grep("_significant", colnames(rowData(data.se))))<1) {

    stop("non signifcant, please run DEP_add_rejections")
  }


}




#' @param data.se SE
#'
#' @return SE
#'

#' @describeIn DEP_Style_se warpper of DEP::test_diff
#' @export
DEP_test_diff <- function(se,type,...){


  mc <- match.call()
  arg.list <- as.list(mc[-1])
  groups <- se$condition%>%
    groupStringFactor()
  if (is.null(arg.list$control)) arg.list$control<- levels(groups)[1]
  if (is.null(arg.list$type)) arg.list$type<- "all"
  do.call(DEP::test_diff,arg.list)
}


DEP_t_test_diff <- function(se,...){

  diff <- DEP_test_diff(se,...)
  contras <- DEP_list_contrast(diff)
  for (i.contra in contras) {

    groups <- strsplit(i.contra,"_vs_")[[1]]
    groups.se <- se[,se$condition%in%groups]
    groups.assay <- assay(groups.se)
    p.val <- apply(groups.assay,1,function(x){t.test_dev(x~groups.se$condition)})
    p.val <-  unname(p.val[rownames(diff)])
    p.adj <- p.adjust(p.val,method = "BH")
    diff@elementMetadata[[paste0(i.contra,"_p.val")]] <-p.val

  }
  return(diff)

}


#' @describeIn DEP_Style_se filter significant feature
#' @export
DEP_filter_significant <- function(data.se,
                                   contrast = DEP_list_contrast(data.se )[1],
                                   top = Inf){


  rd <- DEP_get_diff_table(data.se,contrast = contrast,keep.all = T)%>%
    dplyr::slice_max(`p_value_-log10`,n = top,with_ties = F)%>%
    dplyr::filter(significant)

  data.se[rd$feature_id,]




}


#' @describeIn DEP_Style_se get differential table
#' @export
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



DEP_get_group_color <- function(data.se,col.group = NULL){


  if (is.null(col.group)) {

    groups <- levels(data.se$group)%>%setdiff(c("QC","Blank"))
    col.group <- setNames(ggsci::pal_npg()(length(groups)),groups)
    col.group[c("QC","Blank")] <- c("#FF7F0E","#808180")


  }else{

  }


  data.se$group.color <- col.group[data.se$group]

  return(data.se)

}

#' @describeIn DEP_Style_se plot volcano
#' @export
DEP_plot_volcano <- function(data.se,
                             contrast = DEP_list_contrast(data.se )[1],
                             show.label = T,
                             label.top = 10,
                             label.max.char = 15){
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
    p <- ggplot_sum_patchwork(p.vol.list)
    return(p)
  }


  volcano.data <- DEP_get_diff_table(data.se,contrast = contrast)

  volcano.sig <- volcano.data%>%
    dplyr::filter(significant)
  lfc <- data.se@metadata$DEP_add_rejections$lfc


  p.adjust <- data.se@metadata$DEP_add_rejections$p.adjust

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
 # x.max <- 2
 # y.max <- 3

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
    volcano.data2 <- volcano.data%>%
      dplyr::mutate(label = row.data$label[match(protein , row.data$feature_id)],
                    label = str_short(label,label.max.char))%>%
      dplyr::slice_max(`p_value_-log10`,n = label.top)

    p <- p+
      ggrepel::geom_text_repel(data = filter(volcano.data2, significant),
                               aes(label = label), size = 2.4,
                               box.padding = unit(0.1, "lines"),
                               point.padding = unit(0.1, "lines"),
                               segment.size = 0.1)
  }
  return(p)




}

#' @describeIn DEP_Style_se plot volcano with lipid class
#' @export
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

#' @describeIn DEP_Style_se plot lfc-class
#' @export
DEP.plot.lfc.lipid.class <- function(data.se,
                                     contrast = DEP_list_contrast(data.se )[1],
                                     p.adjust = F){
  if (p.adjust) {
    data.se <- DEP::add_rejections(data.se)

  } else{
    data.se <- add_rejections_no_p.adj(data.se)
  }

  data("lipid.classification",package = "MSdb" )
  plot.data <- DEP::plot_volcano(data.se,contrast,plot = F )%>%
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


DEP_plot_lipid_change_ratio <- function(data.se,
                                        contrast = DEP_list_contrast(data.se )[1]){


  data("lipid.classification",package = "MSdb" )

  plot.data <- DEP_get_diff_table(data.se,contrast = contrast,keep.all = T)%>%
    dplyr::left_join(lipid.classification,by = "Lipid_subclass" )%>%
    dplyr::mutate(
      change = case_when(
        significant&log2_fold_change >0 ~ "UP",
        significant&log2_fold_change <0 ~ "DOWN",
        T ~ "No Change"
      )
    )%>%
    dplyr::group_by(Lipid_subclass)%>%
    dplyr::mutate(sig.ratio = sum(significant)/n())%>%
    dplyr::ungroup()%>%
    dplyr::count(Lipid_subclass,change,sig.ratio)%>%
    dplyr::group_by(Lipid_subclass)%>%
    dplyr::mutate(ratio = n/sum(n),
                  total = sum(n)
                    )%>%
    dplyr::ungroup()%>%
    dplyr::arrange(-total)%>%
    dplyr::mutate(Lipid_subclass = factor(Lipid_subclass,levels = rev(unique(Lipid_subclass) )))

  l <- plot.data%>%
    dplyr::distinct(Lipid_subclass,.keep_all = T)%>%
    dplyr::arrange(sig.ratio,total)%>%
    dplyr::slice_tail(n = 20)%>%
    dplyr::pull(Lipid_subclass)
  plot.data <- plot.data%>%
    dplyr::filter(Lipid_subclass %in% l)

  ggplot(plot.data)+
    geom_bar(aes(y = Lipid_subclass , x = n,fill = change),
             stat = "identity", color = "black")+
    scale_fill_manual(values = c(
      "DOWN" = "#2C7BB6",      # blue
      "No Change" = "#BDBDBD", # grey
      "UP" = "#D7191C"         # red
    ))+
    labs(title = contrast,x = "Count", y = "Lipid Subclass",fill = "Change")+
    theme_bw(base_size = 6)


}

#' @describeIn DEP_Style_se plot heatmap
#' @export
DEP_plot_heatmap <- function(data.se,
                             feature_id = NULL,
                             ...){



  col.info <-SummarizedExperiment::colData(data.se)%>%
    as.data.frame()%>%
    dplyr::mutate(col_group = condition,
                  column_labels = sample.labels,
                  column_split = group)
  row.info <-SummarizedExperiment::rowData(data.se)%>%
                    as.data.frame()%>%
    dplyr::mutate(row_group = "",
                  row_labels = label)
  rownames(data.se) <- row.info$name
  if (!is.null(feature_id)) row.info <- row.info[feature_id,,drop = F]

  heatmap.matrix <-SummarizedExperiment::assay(data.se[row.info$name,col.info$ID])%>%
    `^`(2,.)%>%t%>%scale%>%t


  ComplexHeatmap::ht_opt(TITLE_PADDING = unit(c(4, 4), "points"))
  #ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
  hm <- plotHeatmap(
    heatmap.matrix,
    col.info,row.info,
    column_title_gp = gpar(fill = get_DEP_se_group_color(data.se),fontsize = 12),
    ...
    )
  hm



}



#'
#' @param data.se SE
#' @param file_path path
#' @return null

#' @describeIn DEP_Style_se export data, wirte coldata and rowdata to excel
#' @export
DEP_export_data <- function(data.se,file_path){

  col.data <- colData(data.se)%>%as.data.frame()
  row.data <- rowData(data.se)%>%as.data.frame()%>%
    cbind(assay(data.se))

  data.to.write <-list(sample_info = col.data,
                       compound = row.data)

  xlsx.write.list(data.to.write,file = file_path)
  return(invisible())

}


#' @describeIn DEP_Style_se plot PCA
#' @export
DEP_plot_PCA <- function(data.se,
                         col.group = get_DEP_se_group_color(data.se),
                         showlabel = F,
                         ...){
  se.coldata <- colData(data.se)%>%
    as.data.frame()
  pca.matrix <- assay(data.se)%>%
    `colnames<-`(se.coldata$label)%>%t

  if (is.null(col.group)) {
    col.group <- ggsci::pal_npg()(length(unique(se.coldata$condition)))
  }
  plot_PCA(pca.matrix,
                  pca.group = se.coldata$condition,
                  showlabel = showlabel,
           ...)+
    scale_color_manual(values = col.group)+
    scale_fill_manual(values = col.group)->p.pca

return(p.pca)

}

#' @describeIn DEP_Style_se plot pathway enrich
#'
#' @param filter_Metabolism only output pathway of Metabolism
#' @export
DEP_pathway_enrich <- function(data.se,
                               contrast ,
                               method = c("HyperTest","GlobalTest"),
                               filter_Metabolism =  F){


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
        dplyr::pull(kegg_id)%>%
      analyzePathwayHyperTest(filter_Metabolism= filter_Metabolism)

  }
  if (method == "GlobalTest") {

    data.se <- data.se[!is.na(rowData(data.se)$kegg_id),
                       data.se$condition%in% strsplit(contrast,"_vs_")[[1]]]
    pathway.matrix <- assay(data.se)
    rownames(pathway.matrix) <- rowData(data.se)$kegg_id
    pathway.table <-  analyzePathwayGlobalTest(t(pathway.matrix) ,data.se$condition,filter_Metabolism= filter_Metabolism)


  }

  return(pathway.table)



}


#' @describeIn DEP_Style_se plot gene pathway enrich
#' @export
DEP_pathway_enrich_gene <- function(data.se,
                               contrast ,
                               database = c("KEGG_2021_Human")){


  if (contrast == "all") {
    enrich.list <- list()
    for (i in 1:length( DEP_list_contrast(data.se ))) {
      i_contrast <-  DEP_list_contrast(data.se )[i]
      enrich.list[[i_contrast]] <- DEP_pathway_enrich_gene(data.se ,
                                                      contrast = i_contrast,
                                                      database = database)
    }
    return(enrich.list)
  }
  database <- match.arg(database)
  {
    pathway.table <- DEP_get_diff_table(data.se,
                                        contrast = contrast,
                                        keep.all = T)%>%
      dplyr::filter(significant)%>%
      dplyr::pull(genes)%>%
      enrichR::enrichr(databases = database)
    pathway.table <- pathway.table[[1]]
  }

  return(pathway.table)



}

#' @describeIn DEP_Style_se adjust by weight
#' @export
se_adjuset_by_weight <- function(data.se){

  weight <- data.se$weight
  if (is.null(weight)|all(is.na(weight))) {
    message("no weight input")
    return(data.se)
  }
  weight <- weight/mean(weight,na.rm =T)
  data.matrix <- assay(data.se)
  to.weight <- which(!is.na(weight))
  data.matrix[,to.weight] <- t(t(data.matrix[,to.weight])/weight[to.weight])
  data.matrix -> assay(data.se)

  return(data.se)
}


#' @describeIn DEP_Style_se ANOVA test
#' @export
DEP_test_ANOVA <- function(data.se){


  data.matrix <- assay(data.se)
  p.kruskal <- apply(data.matrix, 1, function(x){
    kruskal.test(x~data.se$condition)$p.value
  })
  p.kruskal.fdr <- p.adjust(p.kruskal)

  rowData(data.se)$p.kruskal <- p.kruskal
  rowData(data.se)$p.kruskal.fdr <- p.kruskal.fdr
  return(data.se)
}



#' @describeIn DEP_Style_se plot bar plot for feature
#' @export
DEP_plot_single_bar <- function(data.se,
                                id){

  plot.data <- data.frame(
    group = data.se$condition,
    val = assay(data.se)[id,]
  )%>%
    dplyr::filter(val >-2)

  ggplot(plot.data ,
         aes(x = group , y = val,col = group) )+
    #geom_bar(stat = "summary",fun  = mean)+
    geom_boxplot()+
    geom_jitter()+
    #ggsignif::geom_signif(
    #  comparisons = apply(combn(names(col.groups),2),2
    #                      ,c,simplify = F),
    #  test  = "t.test",
    #  y_position = seq(8000,12000,length.out = 6))+
    theme_bw()


}



#' @describeIn DEP_Style_se impute with mean
#' @export
DEP_impute_mean <- function(data.se){


  se.data <- assay(data.se)
  se.data <- apply(se.data,2,function(x){
    x[is.na(x)] <- mean(x,na.rm = T)
    return(x)
  })
  se.data -> assay(data.se)
  return(data.se)
}

#' @describeIn DEP_Style_se filter feature with miss value
#' @export
DEP_filter_miss <- function(data.se,group.miss.ratio = 0.3 ){


  datam <- assay(data.se)

  group <- data.se$group

  split_mat <- split(seq_len(ncol(datam)), as.character(group))
  split_list <- lapply(split_mat, function(cols) datam[, cols, drop = FALSE])
  miss.ratio <- lapply(split_list,function(x){  apply(x,1,function(y) sum(is.na(y))/length(y)  )})

  miss.ratio <- do.call(cbind,miss.ratio)
  miss.ratio.min <- apply(miss.ratio,1,function(x) min(x))


  data.se <- data.se[which(miss.ratio.min < group.miss.ratio),]

  return(data.se)

}

#' @describeIn DEP_Style_se filter feature with QC RSD
#' @export
DEP_filter_QC_RSD <- function(data.se,QC_RSD = 0.3){

  if (is.infinite(QC_RSD)) return(data.se)
  rda <- rowData(data.se)
  se <- data.se[which(rda$qc_rsd<QC_RSD),]

  return(se)
}


#' @describeIn DEP_Style_se calculate RSD of QC
#' @export
DEP_get_QC_RSD <- function(data.se){

  qc.se <- data.se[,data.se$group=="QC"]
  rowData(data.se)$qc_rsd <- apply(2^assay(qc.se),1,rsd)
  return(data.se)

}



#' @describeIn DEP_Style_se filter miss, filter QC rsd, normalization, imputation
#' @export
DEP_preprocess <- function(data.se,
                           group.miss.ratio =0.3,
                           QC_RSD = 0.3,
                           keep_before_norm = F){

  #data_filt <- DEP::filter_missval(data.se, thr = min(table(cda$group))*0.3)
  data_filt <- DEP_filter_miss(data.se, group.miss.ratio = group.miss.ratio)
  if ("QC" %in%  data.se$condition){
    data_filt <- DEP_filter_QC_RSD(data_filt, QC_RSD = QC_RSD)
  }else
    message_with_time("No QC, skip QC RSD filter")

  assay.before <- assay(data_filt)
  data_norm <- DEP::normalize_vsn(data_filt)
  if (keep_before_norm) {
    SummarizedExperiment::assay(data_norm,2) <- assay.before
    SummarizedExperiment::assayNames(data_norm) <- c("data","before_norm")
  }
  data_imp <- DEP::impute(data_norm, fun = "MinProb")

  return(data_imp)

}



#' @describeIn DEP_Style_se update `plot_normalization`, add group color
#' @export
DEP_plot_normalization <- function (se, ...)
{
  call <- match.call()
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x, "SummarizedExperiment"),
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "ID", "condition", "replicate") %in%
            colnames(colData(x)))) {
      stop("'label', 'ID', 'condition' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })
  gather_join <- function(se) {
    assay(se) %>% data.frame() %>% tidyr::gather(ID, val) %>% left_join(.,
                                                                 data.frame(colData(se)), by = "ID")
  }
  df <- purrr::map_df(arglist, gather_join, .id = "var") %>% mutate(var = factor(var,
                                                                          levels = names(arglist)))
  ggplot(df, aes(x = ID, y = val, fill = condition)) +
    geom_boxplot(notch = TRUE,
                 na.rm = TRUE) +
    scale_fill_manual(values = get_DEP_se_group_color(se))+
    coord_flip() +
    facet_wrap(~var, ncol = 1) +
    labs(x = "", y = expression(log[2] ~ "Intensity")) +
    DEP::theme_DEP1()
}


#' @describeIn DEP_Style_se get a vector of group color
#' @export
get_DEP_se_group_color <- function(se){

  if(is.null(se$group.color)){
    se$group.color <- ggsci::pal_aaas()(10)[groupStringFactor(se$group)]
  }
  x <- setNames(se$group.color,se$group)
  x <- x[!duplicated(x)]
  #x <- x[levels(se$group)]
  na.omit(x)
}


#' @describeIn DEP_Style_se get significant feature
#' @export
get_DEP_se_sig_feature <- function(data.diff,contrast = DEP_list_contrast(data.diff)[1] ){


  if("all" %in% contrast){

    x <- DEP_get_diff_table(data.diff,"all" )
    fid <- sapply(x,function(y)
      y%>%
        dplyr::filter(significant)%>%
        dplyr::pull(protein)
    )



  }else{
    if (is.numeric(contrast)) {
      contrast <- DEP_list_contrast(data.diff)[contrast]
    }
    fid <- DEP_get_diff_table(data.diff,contrast = contrast )%>%
      dplyr::filter(significant)%>%
      dplyr::pull(protein)
  }

  return(unique(unlist(fid)))


}




#' @describeIn DEP_Style_se import data from ME result
#' @export
#'
get_DEP_se_from_ME_result <- function(ME_file ){

  data.file <- ME_file
  sample.info <- readxl::read_excel(data.file,sheet = "sample.info")
  data <- readxl::read_excel(data.file,sheet = "data")



  col.data <- sample.info%>%
    dplyr::mutate(group = groupStringFactor(group),
                  label = sample.name,
                  sample.labels = label,
                  condition = group,
                  replicate = 1:n())

  row.data <- data%>%
    dplyr::mutate(
      feature_id = paste0("FT",  MSdev:::num2str(1:n())),# num2str
      id = feature_id,
      name = feature_id,
      label = Name,
      kegg_id =`KEGG ID`)


  data_unique <- DEP::make_unique(row.data, "id", "id", delim = ";")
  data.se <- DEP::make_se(data_unique,
                     match(col.data$sample.name,colnames(data_unique)), col.data)


  col.data <- colData(data.se)%>%
    as.data.frame()%>%
    dplyr::mutate(
                  label = sample.name,
                  sample.labels = label,
                  sample.type = case_when(
                    grepl("QC",ignore.case = T,x = label)~"QC",
                    grepl("Blank|BLK",ignore.case = T,x = label)~"Blank",
                    T~"Sample"
                  ),
                  group = groupStringFactor(group),
                  condition = group,
                  replicate = 1:n())

  rownames(col.data) <- colnames(data.se)
  colData(data.se) <- col.data %>%
    S4Vectors::DataFrame()

  data.se <- DEP_get_group_color(data.se)
  return(data.se)

}

#' @describeIn DEP_Style_se remove QC and Blank
#' @export
#'
DEP_remove_QC <- function(data.se,remove_QC = T, remove_Blank = T){

  if (remove_QC) data.se <- data.se[, !data.se$sample.type %in% "QC"]
  if (remove_QC) data.se <- data.se[, !data.se$sample.type %in% "Blank"]

  return(data.se)


}


DEP_test_gsea <-function (dep, databases = c("GO_Molecular_Function_2017b",
                                             "GO_Cellular_Component_2017b", "GO_Biological_Process_2017b"),
                          contrasts = TRUE)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(databases), is.logical(contrasts), length(contrasts) ==
                            1)
  if (!"enrichR" %in% rownames(installed.packages())) {
    stop("test_enrichR() requires the 'enrichR' package",
         "\nTo install the package run: install.packages('enrichR')")
  }
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  libraries <- enrichR::listEnrichrDbs()$libraryName
  if (all(!databases %in% libraries)) {
    stop("Please run `test_gsea()` with valid databases as argument",
         "\nSee http://amp.pharm.mssm.edu/Enrichr/ for available databases")
  }
  if (any(!databases %in% libraries)) {
    databases <- databases[databases %in% libraries]
    message("Not all databases found", "\nSearching the following databases: '",
            paste0(databases, collapse = "', '"), "'")
  }
  message("Background")
  background <- gsub("[.].*", "", row_data$name)
  background_enriched <- enrichR::enrichr(background, databases)
  df_background <- NULL
  for (database in databases) {
    temp <- background_enriched[database][[1]] %>% mutate(var = database)
    df_background <- rbind(df_background, temp)
  }
  df_background$contrast <- "background"
  df_background$n <- length(background)
  OUT <- df_background %>% mutate(bg_IN = as.numeric(gsub("/.*",
                                                          "", Overlap)), bg_OUT = n - bg_IN) %>% select(Term,
                                                                                                        bg_IN, bg_OUT)
  if (contrasts) {
    df <- row_data %>% as.data.frame() %>% select(name,
                                                  ends_with("_significant")) %>% mutate(name = gsub("[.].*",
                                                                                                    "", name))
    df_enrich <- NULL
    for (contrast in colnames(df[2:ncol(df)])) {
      message(gsub("_significant", "", contrast))
      significant <- df[df[[contrast]], ]
      genes <- significant$name
      enriched <- enrichR::enrichr(genes, databases)
      contrast_enrich <- NULL
      for (database in databases) {
        temp <- enriched[database][[1]] %>% mutate(var = database)
        contrast_enrich <- rbind(contrast_enrich, temp)
      }
      contrast_enrich$contrast <- contrast
      contrast_enrich$n <- length(genes)
      cat("Background correction... ")
      contrast_enrich <- contrast_enrich %>% mutate(IN = as.numeric(gsub("/.*",
                                                                         "", Overlap)), OUT = n - IN) %>% select(-n) %>%
        left_join(OUT, by = "Term") %>% mutate(log_odds = log2((IN *
                                                                  bg_OUT)/(OUT * bg_IN)))
      cat("Done.")
      df_enrich <- rbind(df_enrich, contrast_enrich) %>%
        mutate(contrast = gsub("_significant", "", contrast))
    }
  }
  else {
    significant <- row_data %>% as.data.frame() %>% select(name,
                                                           significant) %>% filter(significant) %>% mutate(name = gsub("[.].*",
                                                                                                                       "", name))
    genes <- significant$name
    enriched <- enrichR::enrichr(genes, databases)
    df_enrich <- NULL
    for (database in databases) {
      temp <- enriched[database][[1]] %>% mutate(var = database)
      df_enrich <- rbind(df_enrich, temp)
    }
    df_enrich$contrast <- "significant"
    df_enrich$n <- length(genes)
    cat("Background correction... ")
    df_enrich <- df_enrich %>% mutate(IN = as.numeric(gsub("/.*",
                                                           "", Overlap)), OUT = n - IN) %>% select(-n) %>%
      left_join(OUT, by = "Term") %>% mutate(log_odds = log2((IN *
                                                                bg_OUT)/(OUT * bg_IN)))
    cat("Done.")
  }
  return(df_enrich)
}
