sta_workflow <- function(ms.ana){

  ms.ana <-edit_sample_info(ms.ana)
  ms.ana <- get_feature(ms.ana)
  ms.ana <- get_unique_compound(ms.ana)

  output.dir <- paste0(ms.ana$processing.info$project.info$project.dir,
                       "/Statistic")
  dir.create(output.dir)
  p.pca <- sta_pca(ms.ana )
  export::graph2ppt(p.pca ,
                    file = paste0(output.dir,"/PCA.pptx"),
                    width = 4,height = 4)
  ms.ana<-sta_diff(ms.ana)
  for (i in 1:length(ms.ana$Statistic$Difference)) {
    diff.table <- ms.ana$Statistic$Difference[[i]]$table
    plot_volcano(diff.table,F)->vp
    plot_diff_heatmap(ms.ana , i) ->hp
    dir.create( paste0(output.dir,"/",
                       names(ms.ana$Statistic$Difference)[i]))
    openxlsx::write.xlsx(diff.table , file = paste0(output.dir,"/",
                                                    names(ms.ana$Statistic$Difference)[i],
                                                    "/Diff_",names(ms.ana$Statistic$Difference)[i],".xlsx"))
    export::graph2ppt(vp ,
                      file = paste0(output.dir,"/",
                                     names(ms.ana$Statistic$Difference)[i],
                                    "/Volcano_",names(ms.ana$Statistic$Difference)[i],".pptx"))

    export::graph2ppt(hp ,
                      file = paste0(output.dir,"/",
                                    names(ms.ana$Statistic$Difference)[i],
                                    "/Heatmap",names(ms.ana$Statistic$Difference)[i],".pptx"),
                      width = 8,
                      height = 8)

  }
  ms.ana <- sta_anova(ms.ana)
  hp <- plot_anova_heatmap(ms.ana )

  export::graph2ppt(hp ,
                    file = paste0(output.dir,"/",
                                  "/Heatmap_ANOVA.pptx"),
                    width = 8,height = 8)

  }



sta_pca <- function(ms.ana) {
  sample.info <- ms.ana$sample.info %>%
    dplyr::filter(sample.type != "Blank")
  compound <- ms.ana$compound
  compound.matrix <- compound %>%
    dplyr::ungroup() %>%
    dplyr::select(sample.info$sample.name) %>%
    t %>%
    as.matrix()
  pca.pca <- ropls::opls(x = compound.matrix,
                         predI = 2)
  pca.data <- data.frame(sample.info,
                         pca.pca@scoreMN) %>%
    dplyr::mutate(pca.group = case_when(
      sample.type == "Sample" ~ group,
      sample.type != "Sample" ~ sample.type
    ))
  col.list <-
    ggsci::pal_lancet()(length(unique(pca.data$pca.group)))
  names(col.list) <- unique(pca.data$pca.group)
  col.list["QC"] <- "grey"
  #col.list["Blank"] <- "grey"
  ggplot(pca.data) +
    geom_point(
      aes(x = p1, y = p2 , col = pca.group),
      size = 1,
      alpha = 0.9,
      pch = 16
    ) +
    stat_ellipse(aes(x = p1, y = p2 , fill = pca.group),
                 alpha = 0.2,
                 geom = "polygon") +
    scale_color_manual(values = col.list, labels = c(N = "Normal", P = "HCC")) +
    scale_fill_manual(values = col.list, labels = c(N = "Normal", P = "HCC")) +
    #xlim(-40,40)+
    #ylim(-20,20)+
    xlab(paste0("PC1 (", round(pca.pca@modelDF[["R2X"]][1] * 100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pca.pca@modelDF[["R2X"]][2] * 100, 2), "%)")) +
    labs(
      title = "PCA",
      col = "Group",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 8),
      aspect.ratio = 1,
      legend.key.size = unit(0.1, "inch"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      panel.border = element_rect(fill = NA, size = 0.1)
    ) -> p
  p

}


sta_diff <- function(ms.ana) {
  sample.info <- ms.ana$sample.info
  compound <- ms.ana$compound
  groups.pair <- sample.info$group %>%
    unique() %>%
    setdiff(c("QC", "Blank", NA)) %>%
    combn(m = 2)
  diff.list <- list()
  for (i in 1:ncol(groups.pair)) {
    if (any(grepl(
      pattern = "wt|Con",
      x = groups.pair[, i],
      ignore.case = T
    ))) {
      group.c <-
        grep(
          pattern = "wt|Con",
          x = groups.pair[, i],
          ignore.case = T,
          value = T
        ) %>%
        ifelse(length(.) > 1 , .[1], .)
      group.t <- setdiff(groups.pair[, i] , group.c)

    } else{
      group.c <- groups.pair[1, i]
      group.t <- groups.pair[2, i]
    }

    diff.sample.info <- sample.info %>%
      dplyr::filter(group %in% c(group.c, group.t))
    group.c.id <- which(diff.sample.info$group==group.c)
    group.t.id <- which(diff.sample.info$group==group.t)
    diff.matrix <- compound%>%
      column_to_rownames("feature.id")%>%
      ungroup()%>%
      dplyr::select(diff.sample.info$sample.name)%>%
      as.matrix()
    diff.compound <- compound%>%
      ungroup()%>%
      dplyr::select( !sample.info$sample.name)%>%
      dplyr::mutate(.after = feature.id,
                    foldchange =apply(diff.matrix,1,
                                      function(x)
                                        mean(x[group.t.id])/mean(x[group.c.id])),
                    p.value = apply(diff.matrix,1,
                                    function(x)
                                      t.test_dev(x[group.c.id],x[group.t.id])),
                    logp = -log10(p.value),
                    FDR = p.adjust(p.value),
                    logFDR = -log10(FDR))
    diff.list[[paste0(group.t,"_",group.c)]][["table"]] <- diff.compound
    diff.list[[paste0(group.t,"_",group.c)]][["matrix"]] <- diff.matrix

  }

  {### save

    ms.ana[["Statistic"]][["Difference"]] <-diff.list
    save_ms_ana( ms.ana)


  }
  return(ms.ana)

}


sta_anova <- function(ms.ana){

  sample.info <- ms.ana$sample.info%>%
    dplyr::filter(!sample.type %in% c("QC","Blank"))
  compound <- ms.ana$compound
  anova.matrix <- compound%>%
    ungroup()%>%
    column_to_rownames("feature.id")%>%
    dplyr::select(sample.info$sample.name)%>%
    as.matrix()%>%
    t
  p.values <- rep(1,ncol(anova.matrix))
  anova.data <- data.frame(group = sample.info$group,
                           anova.matrix)
  for (i  in 2:ncol(anova.data)) {

    ft <- colnames(anova.data)[i]
    aov.formula <- paste0(ft,"~ group") %>%as.formula()
    p.values[i-1]<- summary(aov(aov.formula, anova.data) )[[1]][["Pr(>F)"]][1]

  }
  anova.table <- compound%>%
    ungroup()%>%
    dplyr::select(!ms.ana$sample.info$sample.name)%>%
    dplyr::mutate(p.value = p.values,
                  FDR = p.adjust(p.value),
                  logp = -log10(p.value),
                  logFDR = -log10(FDR),.after = feature.id)


  {### save

    ms.ana[["Statistic"]][["ANOVA"]][["table"]] <-anova.table
    ms.ana[["Statistic"]][["ANOVA"]][["matrix"]] <- anova.matrix
    save_ms_ana( ms.ana)


    }
  return(ms.ana)

}






plot_volcano <- function(diff.table,adjusted.p = T){
  diff.table <- diff.table%>%
    dplyr::mutate(logfc = log2(foldchange) ,
                  diff = case_when(logfc > 0.4150375 & FDR <0.05 ~ "up",
                                   logfc < -0.4150375 & FDR <0.05 ~ "down",
                                   T~ "no")
    )
  logfc.max <- max(abs(diff.table$logfc))
  logp.max <-max(abs(diff.table$logp))
  logfdr.max <- max(abs(diff.table$logFDR))

  if (adjusted.p) {
    ggplot(diff.table)+
      geom_point(aes(x = logfc , y = logFDR,col = diff),size = 1,alpha = 0.9,pch = 16)+
      scale_color_manual(values = c(up = "#DC0000" , no = "#BEBEBE",down = "#3C5488"),
                         labels = c("Up","No change","Down"))+
      geom_abline(slope = 0,intercept = -log10(0.05),lty = "dashed" , col = "#E9B574",size = 0.5)+
      geom_vline(xintercept = -0.4150375,lty = "dashed" , col = "#E9B574",size = 0.5)+
      geom_vline(xintercept = 0.4150375,lty = "dashed" , col = "#E9B574",size = 0.5)+
      annotate("text",x = logfc.max*0.8, y = logfdr.max*0.8,label = sum(diff.table$diff == "up"),size= 2.8)+
      annotate("text",x = -logfc.max*0.8, y = logfdr.max*0.8,label = sum(diff.table$diff == "down"),size= 2.8)+
      labs(col = "Significant" , x = "Log2(Foldchange)",y = "-Log10(FDR)")+
      xlim(c(-logfc.max,logfc.max))+
      ylim(c(0,logfdr.max))+
      theme_bw()+
      theme(legend.key.size = unit(0.1,"inch"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.position = "none",
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 8),
            panel.border = element_rect(fill= NA,size = 0.1),
            text = element_text(size=8)) ->vp
  }else{

    diff.table <- diff.table%>%
      dplyr::mutate(logfc = log2(foldchange) ,
                    diff = case_when(logfc > 0.4150375 & p.value <0.05 ~ "up",
                                     logfc < -0.4150375 & p.value <0.05 ~ "down",
                                     T~ "no")
      )
    ggplot(diff.table)+
      geom_point(aes(x = logfc , y = logp,col = diff),size = 1,alpha = 0.9,pch = 16)+
      scale_color_manual(values = c(up = "#DC0000" , no = "#BEBEBE",down = "#3C5488"),
                         labels = c("Up","No change","Down"))+
      geom_abline(slope = 0,intercept = -log10(0.05),lty = "dashed" , col = "#E9B574",size = 0.5)+
      geom_vline(xintercept = -0.4150375,lty = "dashed" , col = "#E9B574",size = 0.5)+
      geom_vline(xintercept = 0.4150375,lty = "dashed" , col = "#E9B574",size = 0.5)+
      annotate("text",x = logfc.max*0.8, y = logp.max*0.8,label = sum(diff.table$diff == "up"),size= 2.8)+
      annotate("text",x = -logfc.max*0.8, y = logp.max*0.8,label = sum(diff.table$diff == "down"),size= 2.8)+
      labs(col = "Significant" , x = "Log2(Foldchange)",y = "-Log10(P value)")+
      xlim(c(-logfc.max,logfc.max))+
      ylim(c(0,logp.max))+
      theme_bw()+
      theme(legend.key.size = unit(0.1,"inch"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.position = "none",
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 8),
            panel.border = element_rect(fill= NA,size = 0.1),
            text = element_text(size=8)) ->vp
  }

  return(vp)


}

plot_anova_heatmap <- function(ms.ana){

  anova.table <- ms.ana$Statistic$ANOVA$table%>%
    dplyr::filter(p.value <0.05)
  anova.matrix <- ms.ana$Statistic$ANOVA$matrix%>%
    scale
  sample.info <- ms.ana$sample.info%>%
    dplyr::filter(!sample.type %in% c("QC","Blank"))%>%
    dplyr::mutate(group = factor(group))
  heatmap.matrix <- anova.matrix[,anova.table[p.value <0.05,]$feature.id]%>%
    t
  col.list <- randomcoloR::randomColor(length(levels(sample.info$group)))
  names(col.list) <- levels(sample.info$group)
  ComplexHeatmap::Heatmap(heatmap.matrix,
                          show_row_names = F,
                          row_names_side = "left",
                          show_column_names = F,
                          show_column_dend = F,
                          column_split= sample.info$group,
                          cluster_column_slices = F,
                          rect_gp = grid::gpar(col = "white",lwd = 0.001),
                          column_gap = unit(2,"mm"),
                         # row_dend_width = unit(40,"mm"),
                         # row_dend_gp = grid::gpar(lwd = 5),
                          row_dend_side = "right",
                          column_title = NULL,
                          top_annotation = columnAnnotation(Diagnose = anno_block(gp = grid::gpar(fill = col.list,
                                                                                            lty = 0),
                                                                                  labels = names(col.list),
                                                                                  labels_gp = grid::gpar(cex = 1))),
                          heatmap_legend_param = list(title = "z-score")
  ) ->hp
  hp

}

plot_diff_heatmap <- function(ms.ana,i ){

  diff.table <- ms.ana$Statistic$Difference[[i]]$table
  diff.matrix <- ms.ana$Statistic$Difference[[i]]$matrix
  diff.sample.info <- ms.ana$sample.info%>%
    dplyr::filter(sample.name %in% colnames(diff.matrix))%>%
    dplyr::mutate(group = factor(group))
  diff.matrix <- diff.matrix[diff.table$feature.id[diff.table$p.value <0.05],]%>%
    t%>%scale%>%t
  col.list <- randomcoloR::randomColor(length(levels(diff.sample.info$group)))
  names(col.list) <- levels(diff.sample.info$group)

  ComplexHeatmap::Heatmap(diff.matrix,
                          show_row_names = F,
                          row_names_side = "left",
                          show_column_names = F,
                          show_column_dend = F,
                          column_split= diff.sample.info$group,
                          cluster_column_slices = F,
                          rect_gp = grid::gpar(col = "white",lwd = 0.001),
                          column_gap = unit(2,"mm"),
                          # row_dend_width = unit(40,"mm"),
                          # row_dend_gp = grid::gpar(lwd = 5),
                          row_dend_side = "right",
                          column_title = NULL,
                          top_annotation = columnAnnotation(Diagnose = anno_block(gp = grid::gpar(fill = col.list,
                                                                                            lty = 0),
                                                                                  labels = names(col.list),
                                                                                  labels_gp = grid::gpar(cex = 1))),
                          heatmap_legend_param = list(title = "z-score")
  ) ->hp
  hp
}


