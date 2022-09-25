analyzeDiff <- function(diff.matrix,diff.group){

  groups.pair <- unique(diff.group)
  if(any(grepl(pattern = "CON|WT",x = groups.pair,ignore.case = T))){
    group.con <- groups.pair[grepl(pattern = "CON|WT",x = groups.pair,ignore.case = T)][1]
    group.case <- setdiff(groups.pair,group.con)

  }else {
    group.con <- sort(groups.pair)[1]
    group.case <- sort(groups.pair)[2]
  }
  idx.con <- which(diff.group==group.con)
  idx.case <- which(diff.group==group.case)

  diff.table <- data.frame(feature_id = colnames(diff.matrix),
                           foldchange =apply(diff.matrix,2,
                                             function(x)
                                               mean(x[idx.case])/mean(x[idx.con])),
                           p.value = apply(diff.matrix,2,function(x){
                             t.test_dev(x[idx.case],x[idx.con])
                           }))%>%
    dplyr::mutate(log2foldchange = log2(foldchange),
                  p.fdr = p.adjust(p.value),
                  log10p = -log10(p.value),
                  log10fdr = -log10(p.fdr))%>%
    cbind(t(diff.matrix))
  diff.table

}

analyzeANOVA <- function(anova.matrix , anova.group){

  p.values <- rep(1,ncol(anova.matrix))
  anova.data <- data.frame(group = anova.group,
                           anova.matrix)
  for (i  in 2:ncol(anova.data)) {

    ft <- colnames(anova.data)[i]
    aov.formula <- paste0(ft,"~ group") %>%as.formula()
    p.values[i-1]<- summary(aov(aov.formula, anova.data) )[[1]][["Pr(>F)"]][1]

  }
  anova.table <- data.frame(feature_id = colnames(anova.matrix),
                            p.value = p.values )%>%
    dplyr::mutate(p.fdr = p.adjust(p.value),
                  log10p = -log10(p.value),
                  log10fdr = -log10(p.fdr))
  anova.table

}

plotPCA <- function(pca.matrix,pca.group){

  pca.pca <- ropls::opls(x = pca.matrix,
                         predI = 5)
  pca.data <- data.frame(pca.group,
                         pca.pca@scoreMN)
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
    scale_color_manual(values = col.list) +
    scale_fill_manual(values = col.list) +
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

plotVolcano <- function(diff.table,p.adjusted = T){

  diff.table <- diff.table%>%
    dplyr::mutate(p = case_when(p.adjusted~p.fdr,
                                T~p.value),
                  log10p = -log10(p),
                  diff =  case_when(log2foldchange > 0.4150375 & p <0.05 ~ "up",
                                    log2foldchange < -0.4150375 & p <0.05 ~ "down",
                                    T~ "no")
                  )
  logfc.max <- max(abs(diff.table$log2foldchange))
  log10p.max <-max(abs(diff.table$log10p))

  ggplot(diff.table)+
    geom_point(aes(x = log2foldchange  , y = log10p,col = diff),size = 1,alpha = 0.9,pch = 16)+
    scale_color_manual(values = c(up = "#DC0000" , no = "#BEBEBE",down = "#3C5488"),
                       labels = c("Up","No change","Down"))+
    geom_abline(slope = 0,intercept = -log10(0.05),lty = "dashed" , col = "#E9B574",size = 0.5)+
    geom_vline(xintercept = -0.4150375,lty = "dashed" , col = "#E9B574",size = 0.5)+
    geom_vline(xintercept = 0.4150375,lty = "dashed" , col = "#E9B574",size = 0.5)+
    annotate("text",x = logfc.max*0.8, y = log10p.max*0.8,label = sum(diff.table$diff == "up"),size= 2.8)+
    annotate("text",x = -logfc.max*0.8, y = log10p.max*0.8,label = sum(diff.table$diff == "down"),size= 2.8)+
    labs(col = "Significant" , x = "Log2(Foldchange)",y = paste0("-Log10(",ifelse( p.adjusted,"FDR","P" ),")"))+
    xlim(c(-logfc.max,logfc.max))+
    ylim(c(0,log10p.max))+
    theme_bw()+
    theme(legend.key.size = unit(0.1,"inch"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.position = "none",
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          panel.border = element_rect(fill= NA,size = 0.1),
          text = element_text(size=8)) ->vp
  vp


}

plotHeatmap <- function(heatmap.matrix,col.info,row.info){

  ComplexHeatmap::Heatmap(heatmap.matrix,
                          name = "Z score",
                          cluster_column_slices = F,
                          column_split = col.info$col.group,
                          row_split = row.info$row.group,


                          column_labels = col.info$col.label,
                          column_names_rot =- 45,
                          row_labels = row.info$row.label,
                          row_names_side = "left",
                          row_names_gp = grid::gpar(fontsize= 6),
                          column_names_gp  = grid::gpar(fontsize= 6),

                          show_column_dend = F,
                          show_row_dend = F
                          )->p
  p


}
