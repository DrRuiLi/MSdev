
plotMSdevPCA <- function(object){

  sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("both","MS1"),
                  sample.type != "Blank")
  pca.matrix <- object@statData$feature%>%
    column_to_rownames("feature_id")%>%
    dplyr::select(sample.info$sample.name)%>%
    t%>%scale

  plotPCA(pca.matrix,sample.info$group)->p
  export::graph2ppt(p,
                    file= paste0(object@projectInfo$projectDir,"/Statistic/PCA.pptx"),
                    width = 4,height = 4)
  p



}

plotMSdevDiffVolcano <- function(object){

  n.diff <- length(object@statData$DifferentialMetabolites)
  metabolites.table <- object@statData$metabolites%>%
    dplyr::select(1:17)

  for (i in 1:n.diff) {

    diff.title <- names(object@statData$DifferentialMetabolites)[i]
    diff.dir <- paste0(object@projectInfo$projectDir,"/Statistic/",diff.title)
    dir.create(diff.dir,recursive = T,showWarnings = F)
    diff.table <- object@statData$DifferentialMetabolites[[i]]%>%
      dplyr::mutate(metabolites.table,.after = log10fdr)
    diff.volcano <- plotVolcano(diff.table,p.adjusted = F)+
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
    diff.row.info <- metabolites.table[diff.table$p.value < 0.05,]%>%
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
    dplyr::filter(xcmsProcessing%in% c("both","MS1"),
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
    anova.row.info <- metabolites.table[which(anova.table$p.value < 0.05),]%>%
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
                      width = 1*nrow(anova.col.info),height = 0.08*nrow(anova.row.info))
    openxlsx::write.xlsx(anova.table,
                         file = paste0(anova.dir,"/ANOVA.",anova.title,".xlsx"))


  }


}

analyzeMSdevDiffMetabolites <- function(object){

  sample.info <- object@sampleInfo%>%
    dplyr::filter(xcmsProcessing%in% c("both","MS1"),
                  sample.type != "Blank",
                  sample.type != "QC")
  sample.groups <- unique(sample.info$group)
  groups.comb <- combn(sample.groups,2)


  for (i in 1:ncol(groups.comb)) {
    groups.pair <- groups.comb[,i]
    if(any(grepl(pattern = "CON|WT",x = groups.pair,ignore.case = T))){
      group.con <- groups.pair[grepl(pattern = "CON|WT",x = groups.pair,ignore.case = T)][1]
      group.case <- setdiff(groups.pair,group.con)

    }else {
      group.con <- sort(groups.pair)[1]
      group.case <- sort(groups.pair)[2]
    }

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
    dplyr::filter(xcmsProcessing%in% c("both","MS1"),
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


exportMSdevSampleInfo <- function(object){
  openxlsx::write.xlsx(object@sampleInfo,
                       file = paste0(object@projectInfo$projectDir,"/Statistic/Sample.info.xlsx"))

  openxlsx::write.xlsx(object@statData$metabolites,
                       file = paste0(object@projectInfo$projectDir,"/Statistic/Metabolites.xlsx"))

}






