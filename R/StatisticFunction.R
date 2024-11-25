analyzeDiff <- function(diff.matrix,diff.group){

  groups.pair <- unique(diff.group)%>%
    groupStringFactor()

  group.con <- levels(groups.pair)[1]
  group.case <- levels(groups.pair)[2]
  idx.con <- which(diff.group==group.con)
  idx.case <- which(diff.group==group.case)

  diff.table <- data.frame(feature_id = colnames(diff.matrix),
                           foldchange =apply(diff.matrix,2,
                                             function(x)
                                               mean(x[idx.case],na.rm =T)/mean(x[idx.con],na.rm =T)),
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

analyzeOddRate <- function(or.matrix,or.group,cov.matrix ){

  or.table <- apply(or.matrix,2 , odd.rate.test , y = or.group, cov.matrix )%>%
    data.table::rbindlist()%>%
    dplyr::mutate(var = colnames(or.matrix))

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

#' @title analyzePathwayGlobalTest
#' @describeIn stat-pathway analyzePathwayGlobalTest
#' @description pathway enrichment by global test
#' @param pathway.matrix should be a matrix, sample as rowname, kegg id as colname
#' @param pathway.group should be a vector with length same as nrow(pathway.matrix), indicate group of sample
#'
#' @return global.test.result
#' @export
#'

analyzePathwayGlobalTest <- function(pathway.matrix,pathway.group ){

  kegg.pathway <- load_as_var("C:/Users/91879/OneDrive/Code/R/data/database.for.compounds.identification/kegg.pathway.database.2022.04.10.Rda")
  pathway.matrix <-data.frame( group = groupStringFactor(pathway.group),
                               scale(pathway.matrix))
  global.test.result <- data.frame(pathway.name =0,
                                   pathway.id = 0,
                                   pathway.class = 0,
                                   Hit=0,
                                   Total=0,
                                   p.value = rep(1,length(kegg.pathway)),
                                   Statistic = 1,
                                   Expected = 1,
                                   Std.dev = 0,
                                   Cov = 0,
                                   compounds = "")
  for (i in 1:length(kegg.pathway)) {

    kegg.pathway.name <- kegg.pathway[[i]][["NAME"]]%>%sub(pattern = " - Homo sapiens (human)",replacement = "",fixed = T)
    kegg.pathway.id <- kegg.pathway[[i]][["ENTRY"]]
    kegg.pathway.class <- kegg.pathway[[i]][["CLASS"]]
    kegg.pathway.compounds <- names(kegg.pathway[[i]][["COMPOUND"]])
    kegg.pathway.hits <- colnames(pathway.matrix)[ colnames(pathway.matrix)%in%kegg.pathway.compounds]
    global.test.result$pathway.name[i] <- kegg.pathway.name
    global.test.result$pathway.id[i] <- kegg.pathway.id
    if (!rlang::is_empty(kegg.pathway.class)) {
      global.test.result$pathway.class[i] <- kegg.pathway.class
    }
    global.test.result$Hit[i] <- length(kegg.pathway.hits)
    global.test.result$Total[i] <- length(kegg.pathway.compounds)
    global.test.result$compounds[i] <- paste(kegg.pathway.hits,collapse = ";")
    if (length(kegg.pathway.hits) ==0) {
      next
    }
    gt.data <-   pathway.matrix[,c("group",kegg.pathway.hits)]

    gt.gt <-globaltest::gt( group~. , data = gt.data,model = "logistic")
    global.test.result[i,6:10] <- gt.gt@result
  }
  global.test.result



}

#' @title analyzePathwayHypertest
#' @describeIn stat-pathway analyzePathwayHypertest
#'
#' @param kegg.id kegg.id
#'
#' @return DF
#' @export
#'

analyzePathwayHyperTest <- function(kegg.id){

  kegg.pathway <- MSdb::get_KEGG_pathway()
  kegg.pathway.compound <- MSdb::get_KEGG_compound_pathway_df()
  kegg.id <- intersect(kegg.id,kegg.pathway.compound$COMPOUND.ID)
  N <- length(unique(kegg.pathway.compound$COMPOUND.ID))  #number of compounds in kegg.pathway.database
  n <- length(kegg.id)
  pathway.hyper.test <- data.frame( pathway.id = 1:length(kegg.pathway) )
  for (i in 1: length(kegg.pathway)) {

    pathway <- kegg.pathway[[i]]

    pathway.hyper.test$pathway.name[i] <- pathway$NAME%>%sub(pattern = " - Homo sapiens (human)",replacement = "",fixed = T)
    pathway.hyper.test$pathway.id[i] <- pathway$ENTRY
    pathway.hyper.test$pathway.class[i] <- ifelse(rlang::is_empty(pathway$CLASS),NA,pathway$CLASS)

    M <- pathway$COMPOUND %>% length
    k <- sum(names(pathway$COMPOUND) %in%  kegg.id )


    pathway.hyper.test$Hit[i] <- k
    pathway.hyper.test$Total[i] <- M
    pathway.hyper.test$p.value[i] <- phyper(k-1,M,N-M,n,lower.tail = F)
    pathway.hyper.test$Statistic[i] <- NA
    pathway.hyper.test$Expected[i] <- NA
    pathway.hyper.test$Std.dev[i] <- NA
    pathway.hyper.test$Cov[i] <- k
    #d<- data.frame(a=c(M-k,N-M-n+k),b=c(k,n-k))
    #pathway.hyper.test$p[i] <- fisher.test(d)$p.value
    compounds.richment <- names(pathway$COMPOUND)[names(pathway$COMPOUND) %in%  kegg.id ]
    pathway.hyper.test$compounds[i] <- ifelse(length(compounds.richment)!=0,stringr::str_c(compounds.richment,collapse = ";"),NA)

  }
  pathway.hyper.test
}

#' plotPCA
#'
#' @param pca.matrix pca.matrix
#' @param pca.group pca.group
#' @param showlabel showlabel
#'
#' @return ggplot
#' @export
#'

plot_PCA <- function(pca.matrix,pca.group,showlabel = F){

  pca.pca <- ropls::opls(x = pca.matrix,
                         crossvalI = ifelse(nrow(pca.matrix)>6,7,nrow(pca.matrix)),
                         predI = 5)
  pca.data <- data.frame(pca.group,
                         pca.label = rownames(pca.matrix),
                         pca.pca@scoreMN)
  col.list <-
    ggsci::pal_lancet()(length(unique(pca.data$pca.group)))
  names(col.list) <- unique(pca.data$pca.group)
  col.list["QC"] <- "grey"
  #col.list["Blank"] <- "grey"

  if (showlabel) {
    p<-  ggplot(pca.data) +
      geom_text(aes(x = p1, y = p2 , label = pca.label,col = pca.group),
                size = 2,show.legend = F)

  }else{
    p <- ggplot(pca.data) +
      geom_point(
        aes(x = p1, y = p2 , col = pca.group),
        size = 1,
        alpha = 0.9,
        pch = 16
      )
  }
    p+
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



plotPLSDA<- function(plsda.matrix,plsda.group,showlabel = F){

  plsda.plsda <- ropls::opls(x = plsda.matrix,
                             y =plsda.group,
                         predI = 5)
  plsda.data <- data.frame(plsda.group,
                         plsda.label = rownames(plsda.matrix),
                         plsda.plsda@scoreMN)
  col.list <-
    ggsci::pal_lancet()(length(unique(plsda.data$plsda.group)))
  names(col.list) <- unique(plsda.data$plsda.group)
  col.list["QC"] <- "grey"
  #col.list["Blank"] <- "grey"

  if (showlabel) {
    p<-  ggplot(plsda.data) +
      geom_text(aes(x = p1, y = p2 , label = plsda.label,col = plsda.group),
                size = 2,show.legend = F)

  }else{
    p <- ggplot(plsda.data) +
      geom_point(
        aes(x = p1, y = p2 , col = plsda.group),
        size = 1,
        alpha = 0.9,
        pch = 16
      )
  }
  p+
    stat_ellipse(aes(x = p1, y = p2 , fill = plsda.group),
                 alpha = 0.2,
                 geom = "polygon") +
    scale_color_manual(values = col.list) +
    scale_fill_manual(values = col.list) +
    #xlim(-40,40)+
    #ylim(-20,20)+
    xlab(paste0("PC1 (", round(plsda.plsda@modelDF[["R2X"]][1] * 100, 2), "%)")) +
    ylab(paste0("PC2 (", round(plsda.plsda@modelDF[["R2X"]][2] * 100, 2), "%)")) +
    labs(
      title = "PLSDA",
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

plotOPLSDA<- function(oplsda.matrix,oplsda.group,showlabel = F){

  oplsda.oplsda <- ropls::opls(x = oplsda.matrix,
                             y =oplsda.group,
                             orthoI= 1,
                             predI = 5)
  oplsda.data <- data.frame(oplsda.group,
                           oplsda.label = rownames(oplsda.matrix),
                           oplsda.oplsda@scoreMN,
                           oplsda.oplsda@orthoScoreMN)
  col.list <-
    ggsci::pal_lancet()(length(unique(oplsda.data$oplsda.group)))
  names(col.list) <- unique(oplsda.data$oplsda.group)
  col.list["QC"] <- "grey"
  #col.list["Blank"] <- "grey"

  if (showlabel) {
    p<-  ggplot(oplsda.data) +
      geom_text(aes(x = p1, y = p2 , label = oplsda.label,col = oplsda.group),
                size = 2,show.legend = F)

  }else{
    p <- ggplot(oplsda.data) +
      geom_point(
        aes(x = p1, y = o1 , col = oplsda.group),
        size = 1,
        alpha = 0.9,
        pch = 16
      )
  }
  xmax <- max(abs(oplsda.data$p1))*1.05
  ymax <- max(abs(oplsda.data$o1))*1.05
  p+
    stat_ellipse(aes(x = p1, y = o1 , fill = oplsda.group),
                 alpha = 0.2,
                 geom = "polygon") +
    geom_hline(yintercept = 0 , linetype = "dashed")+
    geom_vline(xintercept = 0 , linetype = "dashed")+
    geom_text(aes(x = xmax*0.1,y=ymax*0.8 ,
                  label = paste0("R2X = ",oplsda.oplsda@summaryDF$`R2X(cum)`,"\n",
                                 "R2Y = ",oplsda.oplsda@summaryDF$`R2Y(cum)`)),
              hjust = "left",check_overlap = T)+
    scale_color_manual(values = col.list) +
    scale_fill_manual(values = col.list) +
    #xlim(-40,40)+
    #ylim(-20,20)+
    xlab(paste0("PC1 (", round(oplsda.oplsda@modelDF[["R2X"]][1] * 100, 2), "%)")) +
    ylab(paste0("PC2 (", round(oplsda.oplsda@modelDF[["R2X"]][2] * 100, 2), "%)")) +
    labs(
      title = "PLSDA",
      col = "Group",
      fill = "Group"
    ) +
    xlim(c(-xmax,xmax))+
    ylim(c(-ymax,ymax))+
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
    )->p
  p


}

plotVolcano <- function(diff.table,p.adjusted = T,point.label =F){

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

  if(point.label){
    label.df <-diff.table%>%
      dplyr::filter(diff != "no")
    vp <- vp+
      ggrepel::geom_text_repel(data = label.df,
                                aes(x = log2foldchange  , y = log10p,
                                    label = Compound_name),
                               size = 1,
                               segment.size = 0.1,
                               max.overlaps = 30)




  }

  vp


}

plotHeatmap <- function(heatmap.matrix,col.info,row.info){

  ComplexHeatmap::Heatmap(heatmap.matrix,
                          name = "Z score",
                          cluster_column_slices = F,
                          column_split = col.info$column_split,
                          row_split = row.info$column_split,


                          column_labels = col.info$column_labels,
                          column_names_rot =- 45,
                          row_labels = str_short(row.info$row_labels,20),
                          row_names_side = "left",
                          row_names_gp = grid::gpar(fontsize= 6),
                          column_names_gp  = grid::gpar(fontsize= 6),

                          show_column_dend = F,
                          show_row_dend = F,
                          show_row_names = T
                          )->p
  p


}

plotPathwayEnrichment <- function(pathway.table,top = 20 , method = "set1"){

if (!"diff"%in%colnames(pathway.table)) {
  pathway.table$diff <- "NA"
}
  pathway.table <- pathway.table%>%
    dplyr::group_by(diff)%>%
    dplyr::arrange(-p.value)%>%
    dplyr::mutate(p.fdr = p.adjust(p.value),
                  enrich.ratio = Hit/Total,
                  log10p = -log10(p.value),
                  pathway.name = factor(pathway.name , levels = pathway.name)
    )%>%
   # dplyr::filter(grepl(x = pathway.class , pattern = "Metabolism"),
   #               !grepl(x = pathway.class , pattern = "Glycan biosynthesis "),
   #               !grepl(x = pathway.class , pattern = "terpenoids "),
   #               !grepl(x = pathway.class , pattern = "Xenobiotics "),
   #               !grepl(x = pathway.class , pattern = "secondary metabolites")
   #               )%>%
    dplyr::mutate(pathway.class = factor(pathway.class))%>%
    dplyr::slice_tail(n = top)%>%
    dplyr::ungroup()

 if (method == "set1") {
   ### for global test plot
   ggplot(pathway.table)+
     geom_bar(aes(x = pathway.name , y = enrich.ratio , fill = log10p),
              col ="black",stat = "identity",size = 0.1)+
     scale_fill_gradient2(low = "white" ,mid = "yellow",high = "#DC0000",
                          midpoint = median(-log10(x = pathway.table$p.value)))+
     labs(x = NULL ,y = "Enrich Ratio", fill = "-Log10(P)")+
     coord_flip()+
     theme_classic()+
     theme(text = element_text(size = 8),
           legend.key.size = unit(0.1,"inch"),
           legend.text = element_text(size = 8),
           legend.title = element_text(size = 8),
           axis.text = element_text(size = 5),
           axis.title = element_text(size = 5),
           axis.line.x = element_line(size = 0.1),
           axis.line.y = element_line(colour = NA),
           axis.ticks.x = element_line(size = 0.1),
           axis.ticks.y = element_line(colour = NA),
           panel.background = element_blank())->p
   return(p)

 }
  if (method == "set2") {
    ### for hyper test, and distinct up and down
    pathway.table <- pathway.table%>%
      dplyr::filter(diff %in% c("up","down")  )%>%
      dplyr::arrange(abs(log10p))%>%
      dplyr::mutate(log10p = ifelse(diff == "up",1,-1)*log10p,
                    pathway.name = factor(pathway.name , levels = unique(pathway.name)))

    limit.abs <-max(abs(pathway.table$log10p))*1.5
    ggplot(pathway.table)+
      geom_bar(aes(x = pathway.name , y = log10p , fill = log10p),
               col ="black",stat = "identity",size = 0.01)+
      scale_fill_gradient2(low = "#0271B6",high = "#992307",
                           midpoint = 0,
                           limits = c(-limit.abs,limit.abs),
                           breaks = c(-floor(limit.abs),0,floor(limit.abs)),
                           labels = c(floor(limit.abs),0,floor(limit.abs)))+
      scale_y_continuous(limits = c(-limit.abs,limit.abs))+
      labs(x = NULL ,y = "Log10(P)", fill = "-Log10(P)")+
      coord_flip()+
      theme_classic()+
      theme(text = element_text(size = 8),
            legend.key.size = unit(0.1,"inch"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 5),
            axis.line.x = element_line(size = 0.1),
            axis.line.y = element_line(colour = NA),
            axis.ticks.x = element_line(size = 0.1),
            axis.ticks.y = element_line(colour = NA),
            panel.background = element_blank())->p
    p
    return(p)


  }



}


odd.rate.test <- function(x,y,cov.matrix = NULL){

  if (is.factor(y)) {
    or.data <- data.frame(y=y , x= x ,cov.matrix)
    glm.glm <- glm(data = or.data,formula = "y ~ x + .",family = binomial())
    or.result <- epiDisplay::logistic.display(glm.glm,simplified = T)$table
    return(data.frame(
      OR = or.result["x","OR"],
      lower95ci = or.result["x","lower95ci"],
      upper95ci = or.result["x","upper95ci"],
      p.value = or.result["x","Pr(>|Z|)"]
    ))
  }


}


plot_cor_density <- function(x,y,xlab = "x",ylab = "y"){

  plot.data <- data.frame(
    x,y
  )%>%
    dplyr::mutate(x = as.numeric(x),
                  y = as.numeric(y),
                  var = "a")%>%
    dplyr::filter(!is.na(x),!is.na(y))

  ggplot(plot.data,aes(x = x,y = y))+
    geom_point()+
    stat_smooth(method = "lm")+
    ggpubr::stat_cor(label.y.npc = 1)+
    scale_x_continuous(expand = c(0.1,0))+
    scale_y_continuous(expand = c(0.1,0))+
    labs(x = xlab,
         y = ylab)+
    theme_bw()+
    theme(plot.margin = unit(c(0,0,0,0),"inch"))->p.cor
  p.cor
  ggplot(plot.data)+
    geom_density(aes(x = x),fill = "grey")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_void()->p.x
  ggplot(plot.data)+
    geom_density(aes(y=y),fill = "grey")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_void()->p.y



}


