
plotMSdevPCA <- function(object){

  sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("Both","MS1"),
                  !sample.type%in%   c("Blank"))
  pca.matrix <- object@statData$feature%>%
    column_to_rownames("feature_id")%>%
    dplyr::select(sample.info$sample.name)%>%
    t

  plotPCA(pca.matrix,sample.info$group)->p
  dir.create(paste0(object@projectInfo$projectDir,"/Statistic"),recursive = T)
  export::graph2ppt(p,
                    file= paste0(object@projectInfo$projectDir,"/Statistic/PCA.pptx"),
                    width = 4,height = 4)
  p



}

plotMSdevDiffVolcano <- function(object,p.adjusted = T){

  n.diff <- length(object@statData$DifferentialMetabolites)
  metabolites.table <- object@statData$metabolites%>%
    dplyr::select(1:17)

  for (i in 1:n.diff) {

    diff.title <- names(object@statData$DifferentialMetabolites)[i]
    diff.dir <- paste0(object@projectInfo$projectDir,"/Statistic/",diff.title)
    dir.create(diff.dir,recursive = T,showWarnings = F)
    diff.table <- object@statData$DifferentialMetabolites[[i]]%>%
      dplyr::mutate(metabolites.table,.after = log10fdr)
    diff.volcano <- plotVolcano(diff.table,p.adjusted )+
      labs(title = diff.title)
    diff.volcano
    openxlsx::write.xlsx(diff.table,
                         file = paste0(diff.dir,"/",diff.title,".xlsx"))
    export::graph2ppt(diff.volcano,
                      file= paste0(diff.dir,"/",diff.title,".pptx"),
                      width = 3,height = 3)



  }


}

plotMSdevDiffHeatmap <- function(object){


  sample.info <- object@sampleInfo
  n.diff <- length(object@statData$DifferentialMetabolites)
  metabolites.table <- object@statData$metabolites
  for (i in 1:n.diff) {

    diff.title <- names(object@statData$DifferentialMetabolites)[i]
    diff.dir <- paste0(object@projectInfo$projectDir,"/Statistic/",diff.title)
    dir.create(diff.dir,recursive = T,showWarnings = F)

    diff.table <- object@statData$DifferentialMetabolites[[i]]
    diff.col.info <- sample.info%>%
      dplyr::filter(sample.name %in% colnames(diff.table))%>%
      dplyr::mutate(col.group = groupStringFactor(group),
                    col.label = label)
    diff.row.info <- metabolites.table[diff.table$p.fdr < 0.05,]%>%
      dplyr::mutate(row.label = fixStringLength(Compound_name , 20),
                    row.group = " ")
    diff.matrix <- object@statData$metabolites%>%
      dplyr::filter(feature_id %in% diff.row.info$feature_id)%>%
      dplyr::select(diff.col.info$sample.name)%>%
      t%>%scale%>%t

    diff.heatmap<- plotHeatmap(heatmap.matrix = diff.matrix ,
                               col.info = diff.col.info,
                               row.info = diff.row.info)
    diff.heatmap

    export::graph2pdf(diff.heatmap,
                      file= paste0(diff.dir,"/Heatmap.",diff.title,".pptx"),
                      width = 1*nrow(diff.col.info),height = 0.08*nrow(diff.row.info))



  }



}

plotMSdevANOVA <- function(object){

  sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("Both","MS1"),
                  sample.type != "Blank",
                  sample.type != "QC")
  n.anova <- length(object@statData$ANOVA)
  metabolites.table <- object@statData$metabolites

  for (i in 1:n.anova) {

    anova.title <- names(object@statData$ANOVA)[i]
    anova.dir <- paste0(object@projectInfo$projectDir,"/Statistic/",anova.title)
    dir.create(anova.dir,recursive = T,showWarnings = F)

    anova.table <- object@statData$ANOVA[[i]]%>%
      cbind(metabolites.table)

    if (anova.title == "All Group") {
      anova.group <- sample.info$group%>%
        unique()
    }else{
      anova.group <- sample.info$group[str_detect(pattern = sample.info$group,string = anova.title)]%>%
        unique()

    }

    anova.col.info <- sample.info%>%
      dplyr::filter(group %in% anova.group)%>%
      dplyr::mutate(col.group = groupStringFactor(group),
                    col.label = label)
    anova.row.info <- metabolites.table[which(anova.table$p.fdr < 0.05),]%>%
      dplyr::mutate(row.label = fixStringLength(Compound_name , 20),
                    row.group = " ")
    anova.matrix <- object@statData$metabolites%>%
      dplyr::filter(feature_id %in% anova.row.info$feature_id)%>%
      dplyr::select(anova.col.info$sample.name)%>%
      t%>%scale%>%t

    anova.heatmap<- plotHeatmap(heatmap.matrix = anova.matrix ,
                               col.info = anova.col.info,
                               row.info = anova.row.info)
    anova.heatmap

    export::graph2pdf(anova.heatmap,
                      file= paste0(anova.dir,"/Heatmap.",anova.title,".pptx"),
                      width = 0.2*nrow(anova.col.info),height = 0.08*nrow(anova.row.info)+3)
    openxlsx::write.xlsx(anova.table,
                         file = paste0(anova.dir,"/ANOVA.",anova.title,".xlsx"))


  }


}


plotMSdevDiffVennDiagram <- function(object,
                                     diff.select = "manual",
                                     change = c("both","up","down"),
                                     p = "p.value"){

  diff.all <-data.frame(idx = 1:length(object@statData$DifferentialMetabolites),
                        selected = names(object@statData$DifferentialMetabolites))
  metabolite.table <- object@statData$metabolites[,1:18]
  if (diff.select == "manual") {
    diff.select <- edit_df_in_excel(diff.all)%>%
      dplyr::filter(selected%in%diff.all$selected)

    message("selected: ", paste0(diff.select$idx,collapse = ", "),"\n",
            paste0(diff.select$selected,collapse = "\n"))
  }else{
    diff.select <- diff.all[diff.select,]
  }

  if ("both" %in% change) {
    lapply(object@statData$DifferentialMetabolites, function(x){
      x%>%
        dplyr::filter(eval(str2expression(p)) <0.05)%>%
        pull(feature_id)
    })->venn.list
    venn.list <- venn.list[diff.select$idx]
    #dev.off()
    VennDiagram::venn.diagram(venn.list,
                              width = 10,height = 10,
                              filename = NULL,
                              fill =pal_aaas()(nrow(diff.select)),
                              alpha = 0.7,
                              col = "white",
                              cex= 0.67)%>%
      grid.draw()%>%
      export::graph2ppt(file = paste0(object@projectInfo$projectDir,"/Statistic/VennDiagram_All_Significant_Metabolites.pptx"),
                        width = 4,height = 4)
    venn.overlap <-VennDiagram::get.venn.partitions (venn.list)
    venn.overlap.list <- lapply(venn.overlap$..values.., function(x){
        metabolite.table%>%
        dplyr::filter(feature_id %in% x)

    })
    names(venn.overlap.list) <- venn.overlap$..set..
    for (i in 1:length(venn.overlap.list)) {
      openxlsx::write.xlsx(venn.overlap.list[[i]],
                           file =paste0(object@projectInfo$projectDir,"/Statistic/VennDiagram_All_Significant_Metabolites.Overlap.",
                                        names(venn.overlap.list)[i],".xlsx") )
    }

  }

  if ("up" %in% change) {
    lapply(object@statData$DifferentialMetabolites, function(x){
      x%>%
        dplyr::filter(eval(str2expression(p)) <0.05,foldchange >1)%>%
        pull(feature_id)
    })->venn.list
    venn.list <- venn.list[diff.select$idx]
    #dev.off()
    VennDiagram::venn.diagram(venn.list,
                              width = 10,height = 10,
                              filename = NULL,
                              fill =pal_aaas()(nrow(diff.select)),
                              alpha = 0.7,
                              col = "white",
                              cex= 0.67)%>%
      grid.draw()%>%
      export::graph2ppt(file = paste0(object@projectInfo$projectDir,"/Statistic/VennDiagram_Up_Metabolites.pptx"),
                        width = 4,height = 4)
    venn.overlap <-VennDiagram::get.venn.partitions (venn.list)
    venn.overlap.list <- lapply(venn.overlap$..values.., function(x){
      metabolite.table%>%
        dplyr::filter(feature_id %in% x)

    })
    names(venn.overlap.list) <- venn.overlap$..set..
    for (i in 1:length(venn.overlap.list)) {
      openxlsx::write.xlsx(venn.overlap.list[[i]],
                           file =paste0(object@projectInfo$projectDir,"/Statistic/VennDiagram_Up_Metabolites.Overlap.",
                                        names(venn.overlap.list)[i],".xlsx") )
    }

  }

  if ("down" %in% change) {
    lapply(object@statData$DifferentialMetabolites, function(x){
      x%>%
        dplyr::filter(eval(str2expression(p)) <0.05,foldchange <1)%>%
        pull(feature_id)
    })->venn.list
    venn.list <- venn.list[diff.select$idx]
    #dev.off()
    VennDiagram::venn.diagram(venn.list,
                              width = 10,height = 10,
                              filename = NULL,
                              fill =pal_aaas()(nrow(diff.select)),
                              alpha = 0.7,
                              col = "white",
                              cex= 0.67)%>%
      grid.draw()%>%
      export::graph2ppt(file = paste0(object@projectInfo$projectDir,"/Statistic/VennDiagram_Down_Metabolites.pptx"),
                        width = 4,height = 4)
    venn.overlap <-VennDiagram::get.venn.partitions (venn.list)
    venn.overlap.list <- lapply(venn.overlap$..values.., function(x){
      metabolite.table%>%
        dplyr::filter(feature_id %in% x)

    })
    names(venn.overlap.list) <- venn.overlap$..set..
    for (i in 1:length(venn.overlap.list)) {
      openxlsx::write.xlsx(venn.overlap.list[[i]],
                           file =paste0(object@projectInfo$projectDir,"/Statistic/VennDiagram_Down_Metabolites.Overlap.",
                                        names(venn.overlap.list)[i],".xlsx") )

  }



  }
}

plotMSdevPathway <- function(object,method = "set1",topN=20){

  n.pathway <- length(object@statData$PathwayEnrichment)
  metabolites.table <- object@statData$metabolites%>%
    dplyr::select(1:17)

  for (i in 1:n.pathway) {

    pathway.title <- names(object@statData$PathwayEnrichment)[i]
    pathway.dir <- paste0(object@projectInfo$projectDir,"/Statistic/",pathway.title)
    dir.create(pathway.dir,recursive = T,showWarnings = F)
    pathway.table <- object@statData$PathwayEnrichment[[i]]%>%
      dplyr::arrange(p.value)

    pathway.plot <- plotPathwayEnrichment(pathway.table,
                                          method = method,
                                          top = topN)+
      labs(title = pathway.title)
    pathway.plot
    openxlsx::write.xlsx(pathway.table,
                         file = paste0(pathway.dir,"/PathwayEnrichment.",pathway.title,".xlsx"))
    export::graph2ppt(pathway.plot,
                      file= paste0(pathway.dir,"/PathwayEnrichment.",pathway.title,".pptx"),
                      width = 4,height = 0.8+0.1*topN)



  }



}


analyzeMSdevDiffMetabolites <- function(object){

  sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("Both","MS1"),
                  sample.type != "Blank",
                  sample.type != "QC")
  sample.groups <- unique(sample.info$group)
  groups.comb <- combn(sample.groups,2)


  for (i in 1:ncol(groups.comb)) {
    groups.pair <- groups.comb[,i]%>%
      groupStringFactor()


    group.con <- levels(groups.pair)[1]
    group.case <- levels(groups.pair)[2]


    diff.sample.info <- sample.info%>%
      dplyr::filter(group %in% groups.pair)
    diff.matrix <- object@statData$metabolites %>%
      column_to_rownames("feature_id")%>%
      dplyr::select(diff.sample.info$sample.name)%>%
      t


    diff.table <- analyzeDiff(diff.matrix,diff.sample.info$group)
    object@statData$DifferentialMetabolites[[paste0(group.case," vs ",group.con)]] <- diff.table

  }

  object

}


analyzeMSdevANOVA <- function(object,groupANOVA = "All Group"){


  anova.sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("Both","MS1"),
                  sample.type != "Blank",
                  sample.type != "QC")
  if (!"All Group"%in%groupANOVA) {
    anova.sample.info <- dplyr::filter(anova.sample.info,group %in% groupANOVA)
  }


  anova.matrix <- object@statData$metabolites %>%
    column_to_rownames("feature_id")%>%
    dplyr::select(anova.sample.info$sample.name)%>%
    t

  anova.table <-analyzeANOVA(anova.matrix,anova.sample.info$group)

  object@statData$ANOVA[[paste0(unique(groupANOVA),collapse = "_and_")]] <- anova.table
  object

}


analyzeMSdevPathway <- function(object,method = "global.test",p.adjusted=F){

  sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("Both","MS1"),
                  sample.type != "Blank",
                  sample.type != "QC")
  sample.groups <- unique(sample.info$group)
  groups.comb <- combn(sample.groups,2)

  metabolites.table <- object@statData$metabolites%>%
    dplyr::select(1:17)

  if (method == "global.test") {
    for (i in 1:ncol(groups.comb)) {
      groups.pair <- groups.comb[,i]%>%
        groupStringFactor()


      group.con <- levels(groups.pair)[1]
      group.case <- levels(groups.pair)[2]


      pathway.sample.info <- sample.info%>%
        dplyr::filter(group %in% groups.pair)
      pathway.matrix <- object@statData$metabolites %>%
        dplyr::filter(!is.na(kegg.id))%>%
        dplyr::distinct(kegg.id,.keep_all = T)%>%
        column_to_rownames("kegg.id")%>%
        dplyr::select(pathway.sample.info$sample.name)%>%
        t
      pathway.table <- analyzePathwayGlobalTest(pathway.matrix,pathway.sample.info$group)%>%
        dplyr::mutate(diff = "all")
      object@statData$PathwayEnrichment[[paste0(group.case," vs ",group.con)]] <- pathway.table

    }

  }
  if (method == "hyper.test") {
    for (i in 1:length(object@statData$DifferentialMetabolites)) {


      diff.table <- object@statData$DifferentialMetabolites[[i]]%>%
        dplyr::mutate(p = case_when(p.adjusted~p.fdr,
                                    T~p.value),
                      log10p = -log10(p),
                      diff =  case_when(log2foldchange > 0.4150375 & p <0.05 ~ "up",
                                        log2foldchange < -0.4150375 & p <0.05 ~ "down",
                                        T~ "no"),
                      kegg.id = metabolites.table$kegg.id[match(.$feature_id,metabolites.table$feature_id)]
        )%>%data.table::as.data.table()
      diff.title <- names(object@statData$DifferentialMetabolites)[i]

      pathway.table <- rbind(analyzePathwayHypertest(diff.table$kegg.id[diff.table$diff=="down"])%>%
                               dplyr::mutate(diff = "down"),
                             analyzePathwayHypertest(diff.table$kegg.id[diff.table$diff=="up"])%>%
                               dplyr::mutate(diff = "up"),
                             analyzePathwayHypertest(diff.table$kegg.id[diff.table$diff%in% c("up","down")])%>%
                               dplyr::mutate(diff = "all"))

      object@statData$PathwayEnrichment[[diff.title]] <- pathway.table

    }

  }
  return(object)


}


exportMSdev <- function(object){

  dir.create(paste0(object@projectInfo$projectDir,"/Statistic"),recursive = T)
  openxlsx::write.xlsx(object@sampleInfo,
                       file = paste0(object@projectInfo$projectDir,"/Statistic/Sample.info.xlsx"))

  openxlsx::write.xlsx(object@statData$metabolites,
                       file = paste0(object@projectInfo$projectDir,"/Statistic/Metabolites.xlsx"))

}






