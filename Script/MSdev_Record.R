# Sat Jan 13 10:43:20 2024 ------------------------------
### DCX Metabolomics
{
  msdev.dcx <- MSdev("d:/2024_01_08-Duchenxi/Data/")
  msdev.dcx <- load_as_var("d:/2024_01_08-Duchenxi/MSdev_2024_01_13.Rdata")
  msdev.dcx <- MSdev_msConvert(msdev.dcx)
  msdev.dcx <- MSdev_checkSampleInfo(msdev.dcx)
  msdev.dcx <- MSdev_xcmsProcessing(msdev.dcx)
  msdev.dcx <- MSdev_extract_Spectra(msdev.dcx)
  msdev.dcx <- MSdev_match_Spectra_to_feature(msdev.dcx)
  msdev.dcx <- MSdev_annotation(msdev.dcx,
                                expand_adduct= T,
                                cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.dcx <- MSdev_get_Stat(msdev.dcx,QC_RSD = Inf)
  MSdev_save(msdev.dcx)
  MSdev_export(msdev.dcx)


}




# Sun Jan 14 12:31:01 2024 ------------------------------
### WYQ
{
  msdev.wyq <- MSdev("d:/WYQ/2024_01_10-Wangyongqiang/Data/")
  msdev.wyq <- MSdev_msConvert(msdev.wyq)
  msdev.wyq <- MSdev_checkSampleInfo(msdev.wyq)
  msdev.wyq <- MSdev_xcmsProcessing(msdev.wyq)
  msdev.wyq <- MSdev_extract_Spectra(msdev.wyq)
  msdev.wyq <- MSdev_match_Spectra_to_feature(msdev.wyq)
  msdev.wyq <- MSdev_annotation(msdev.wyq,
                                expand_adduct= T,
                                selected_adduct = c("[M-H]-",
                                                    "[M-H2O-H]-",
                                                    "[2M-H]-" ,
                                                    "[M+FA-H]-" ,
                                                    "[M+H]+" ,
                                                    "[M-H2O+H]+",
                                                    "[M+NH4]+"
                                ) ,
                                db.path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/SpectraDB/HMDB_KEGG_export_human_pathway.rda")
  msdev.wyq <- MSdev_get_Stat(msdev.wyq)
  MSdev_save(msdev.wyq)
  MSdev_export(msdev.wyq)

  ### STAT
  {
    msdev.wyq <- load_as_var("d:/WYQ/2024_01_10-Wangyongqiang/MSdev_2024_01_14.Rdata")
    msdev.wyq <- MSdev_checkSampleInfo(msdev.wyq)
    msdev.wyq <- MSdev_update_xcms_pdata(msdev.wyq)
    {

      proj.dir <- msdev.wyq@projectInfo$projectDir
      data.se <- get_MSdev_DEP_se(msdev.wyq,from = "metabolite")
      data.se$condition <- data.se$group <- data.se$condition%>%
        gsub(pattern = "_",replacement = "")
      p.pca <- DEP_plot_PCA(data.se)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5)


      ### Sample_P
      data.se.Sample_P <- DEP_normalization(data.se)
      data.se.Sample_P <- data.se.Sample_P[,data.se.Sample_P$sample.type== "Sample_P"]
      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 7,append = T)

      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)
      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 7,append = T)
      table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(table.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
      )

      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      rowData(data.diff)$kegg.id <- rowData(data.diff)$kegg_id
      diff.path <- DEP_pathway_enrich(data.diff,
                                      contrast = "all")
      p.diff.path <- lapply(diff.path,
                  plotPathwayEnrichment)%>%
        ggplot_sum_patchwork()+
        plot_layout(ncol = 2)
      export_graph2pdf(p.diff.path , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 15,append = T)
      xlsx.write.list(diff.path,
                      paste0(proj.dir,"/Statistic/Pathway.enrich.xlsx"))


      ### Sample_P
      data.se.Sample_P <- DEP_normalization(data.se)
      data.se.Sample_P <- data.se.Sample_P[,data.se.Sample_P$sample.type== "Sample_ACNH"]
      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 7,append = T)

      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)
      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 7,append = T)
      table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(table.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites2.xlsx")
      )

      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      rowData(data.diff)$kegg.id <- rowData(data.diff)$kegg_id
      diff.path <- DEP_pathway_enrich(data.diff,
                                      contrast = "all")
      p.diff.path <- lapply(diff.path,
                            plotPathwayEnrichment)%>%
        ggplot_sum_patchwork()+
        plot_layout(ncol = 2)
      export_graph2pdf(p.diff.path , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 15,append = T)
      xlsx.write.list(diff.path,
                      paste0(proj.dir,"/Statistic/Pathway.enrich2.xlsx"))

    }

  }



}


### WYQ
{
  msdev.wyq <- MSdev("d:/WYQ/2024_01_11-Wangyongqiang/Data/")
  msdev.wyq <- MSdev_msConvert(msdev.wyq)
  msdev.wyq <- MSdev_checkSampleInfo(msdev.wyq)
  msdev.wyq <- MSdev_xcmsProcessing(msdev.wyq)
  msdev.wyq <- MSdev_extract_Spectra(msdev.wyq)
  msdev.wyq <- MSdev_match_Spectra_to_feature(msdev.wyq)
  msdev.wyq <- MSdev_annotation(msdev.wyq,
                                db.path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/SpectraDB/LipidBlast_export.rda")
  msdev.wyq <- MSdev_get_Stat(msdev.wyq, QC_RSD = Inf)
  MSdev_save(msdev.wyq)
  MSdev_export(msdev.wyq)


  ### STAT
  {
    msdev.wyq <- load_as_var("d:/WYQ/2024_01_11-Wangyongqiang/MSdev_2024_01_23.Rdata")
    msdev.wyq <- MSdev_checkSampleInfo(msdev.wyq)
    msdev.wyq <- MSdev_update_xcms_pdata(msdev.wyq)
    msdev.wyq <- MSdev_get_Stat(msdev.wyq)
    {

      proj.dir <- msdev.wyq@projectInfo$projectDir
      data.se <- get_MSdev_DEP_se(msdev.wyq,from = "metabolite")
      data.se$condition <- data.se$group <- data.se$condition%>%
        gsub(pattern = "_",replacement = "")
      p.pca <- DEP_plot_PCA(data.se)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5)


      ### sample_p
      data.se.Sample_P <- DEP_normalization(data.se)
      data.se.Sample_P <- data.se.Sample_P[,data.se.Sample_P$sample.type== "Sample_P"]
      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 30,append = T)
      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 30,append = T)
      table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(table.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
      )

      ### sample_p
      data.se.Sample_P <- DEP_normalization(data.se)
      data.se.Sample_P <- data.se.Sample_P[,data.se.Sample_P$sample.type== "Sample_A"]
      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 30,append = T)
      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 30,append = T)
      table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(table.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
      )




    }

  }



}



### 2024.3.11 HGH lipidomic
{
  msdev.hgh  <- MSdev("D:/2024_03_07-HGH/Data/")
  msdev.hgh <- load_as_var("d:/2024_03_07-HGH/MSdev_2024_03_11.Rdata")
  msdev.hgh <- MSdev_msConvert(msdev.hgh)
  msdev.hgh <- MSdev_checkSampleInfo(msdev.hgh)
  msdev.hgh <- MSdev_xcmsProcessing(msdev.hgh)
  msdev.hgh <- MSdev_extract_Spectra(msdev.hgh)
  msdev.hgh <- MSdev_match_Spectra_to_feature(msdev.hgh)
  msdev.hgh <- MSdev_annotation(msdev.hgh,cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast3.sqlite",
                                selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.hgh <- MSdev_get_Stat(msdev.hgh)
  MSdev_export(msdev.hgh)
  MSdev_save(msdev.hgh)

  {

    proj.dir <- msdev.hgh@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.hgh,from = "metabolite")

    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )






  }



}


### 2024.3.11 fudan lipidomic
{
  msdev.fudan  <- MSdev("D:/2024_03_07-Fudan/Data/")
  msdev.fudan <- load_as_var("D:/2024_03_07-Fudan/MSdev_2024_03_11.Rdata")
  msdev.fudan <- MSdev_msConvert(msdev.fudan)
  msdev.fudan <- MSdev_checkSampleInfo(msdev.fudan)
  msdev.fudan <- MSdev_xcmsProcessing(msdev.fudan)
  msdev.fudan <- MSdev_extract_Spectra(msdev.fudan)
  msdev.fudan <- MSdev_match_Spectra_to_feature(msdev.fudan)
  msdev.fudan <- MSdev_annotation(msdev.fudan,cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast3.sqlite",
                                selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.fudan <- MSdev_get_Stat(msdev.fudan)
  MSdev_export(msdev.fudan)
  MSdev_save(msdev.fudan)

  {

    proj.dir <- msdev.fudan@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.fudan,from = "metabolite")
    #data.se <- data.se[,!colnames(data.se)%in% c( "WT_C1",   "WT_C2" ,  "WT_C3") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )






  }



}


# Wed Apr 24 13:50:45 2024 ------------------------------
{
  msdev.CX  <- MSdev("D:/2024_04_25-Chenxin/Data/")
  #msdev.CX <- load_as_var("")
  msdev.CX <- MSdev_msConvert(msdev.CX)
  msdev.CX <- MSdev_checkSampleInfo(msdev.CX)
  msdev.CX <- MSdev_xcmsProcessing(msdev.CX)
  msdev.CX <- MSdev_extract_Spectra(msdev.CX)
  msdev.CX <- MSdev_match_Spectra_to_feature(msdev.CX)
  msdev.CX <- MSdev_annotation(msdev.CX,cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                                  selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.CX <- MSdev_get_Stat(msdev.CX,candi = T,QC_RSD = 10)
  MSdev_export(msdev.CX,candi = F)
  MSdev_save(msdev.CX)

  {

    proj.dir <- msdev.CX@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.CX,from = "metabolite")
    #data.se <- data.se[,!colnames(data.se)%in% c( "WT_C1",   "WT_C2" ,  "WT_C3") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca,
                     paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}
# Thu May 23 13:45:06 2024 cx------------------------------
{
  msdev.CX  <- MSdev("D:/2024_05_20-Chenxin/Data/")
  msdev.CX <- load_as_var("D:/2024_05_20-Chenxin/MSdev_2024_05_23.Rdata")
  msdev.CX <- MSdev_msConvert(msdev.CX)
  msdev.CX <- MSdev_checkSampleInfo(msdev.CX)
  msdev.CX <- MSdev_xcmsProcessing(msdev.CX)
  msdev.CX <- MSdev_extract_Spectra(msdev.CX)
  msdev.CX <- MSdev_match_Spectra_to_feature(msdev.CX)
  msdev.CX <- MSdev_annotation(msdev.CX,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.CX <- MSdev_get_Stat(msdev.CX,candi = T,QC_RSD = 10)
  MSdev_export(msdev.CX,candi = F)
  MSdev_save(msdev.CX)

  {

    proj.dir <- msdev.CX@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.CX,from = "metabolite")
    data.se <- data.se[,!colnames(data.se)%in% c( "WT_C1",   "WT_C2" ,  "WT_C3") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}




# Fri Jun  7 09:14:22 2024 WYQ------------------------------
{
  msdev.wyq <- MSdev("d:/2024_06_04-Wangyongqiang/Data/")
  msdev.wyq <- load_as_var("d:/2024_06_04-Wangyongqiang/MSdev_2024_06_07.Rdata")
  msdev.wyq <- MSdev_msConvert(msdev.wyq)
  msdev.wyq <- MSdev_checkSampleInfo(msdev.wyq)
  msdev.wyq <- MSdev_xcmsProcessing(msdev.wyq)
  msdev.wyq <- MSdev_extract_Spectra(msdev.wyq)
  msdev.wyq <- MSdev_match_Spectra_to_feature(msdev.wyq)
  msdev.wyq <- MSdev_annotation(msdev.wyq,
                                expand_adduct= T,
                                #selected_adduct = c("[M-H]-",
                                #                    "[M-H2O-H]-",
                                #                    "[2M-H]-" ,
                                #                    "[M+FA-H]-" ,
                                #                    "[M+H]+" ,
                                #                    "[M-H2O+H]+",
                                #                    "[M+NH4]+"
                                #) ,
                                cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.wyq <- MSdev_get_Stat(msdev.wyq,QC_RSD = 0.3)
  MSdev_save(msdev.wyq)
  MSdev_export(msdev.wyq,candi = F)
  ### STAT
  {

      proj.dir <- msdev.wyq@projectInfo$projectDir
      data.se <- get_MSdev_DEP_se(msdev.wyq,from = "metabolite")
      data.replace <- readxl::read_excel("d:/2024_06_04-Wangyongqiang/Statistic/Metabolites2.xlsx",sheet = 2)
      rowData(data.se)$label <- data.replace$name[match(rownames(data.se),
                                                       data.replace$feature_id)]
      data.se$condition <- data.se$group <- data.se$condition%>%
        gsub(pattern = "_",replacement = "")


      data.se  <- data.se[,!data.se$sample.type=="Blank"]




      p.pca <- DEP_plot_PCA(data.se)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5)


      ### Sample_P
      data.se.Sample_P <- DEP_normalization(data.se)
      data.se.Sample_P <- data.se.Sample_P[,data.se.Sample_P$sample.type== "Sample"]
      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 7,append = T)

      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)
      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff ,
                       paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 7,append = T)
      table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(table.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
      )

      data.diff <- DEP_test_diff(data.se.Sample_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      rowData(data.diff)$kegg.id <- rowData(data.diff)$kegg_id
      diff.path <- DEP_pathway_enrich(data.diff,
                                      contrast = "all")
      p.diff.path <- lapply(diff.path,
                            plotPathwayEnrichment)%>%
        ggplot_sum_patchwork()+
        plot_layout(ncol = 2)
      export_graph2pdf(p.diff.path , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 15,append = T)
      xlsx.write.list(diff.path,
                      paste0(proj.dir,"/Statistic/Pathway.enrich.xlsx"))


      data.diff$sample.label <- data.diff$ID
      fid <- table.diff$KOF_vs_KOCON%>%
        dplyr::filter(significant)%>%
        dplyr::pull(protein)
      data.plot <- data.diff[fid,]
      rda <- rowData(data.plot)
      cda <- colData(data.plot)
      hm <- assay(data.plot)%>%
        t%>%scale()%>%t

      ComplexHeatmap::Heatmap(
        hm,name = "Z score",
        show_row_names = F,
        show_column_names = F,
        cluster_columns = F,
        column_split = cda$group)
      export::graph2pdf(file = "d:/temp/heatmap.pdf",
                        width = 10,height = 10)

    }




}

# Wed Jun 12 13:50:53 2024 LJW------------------------------
{
  msdev.ljw <- MSdev("d:/LJW_MS data20240611/rawdata/",
                     experimentInfo =MS_Experiment[10] )
  msdev.ljw <- load_as_var("d:/LJW_MS data20240611/MSdev_2024_06_12.Rdata")
  msdev.ljw <- MSdev_msConvert(msdev.ljw)
  msdev.ljw <- MSdev_checkSampleInfo(msdev.ljw)
  msdev.ljw <- MSdev_xcmsProcessing(msdev.ljw)
  msdev.ljw <- MSdev_annotation(msdev.ljw,
                                expand_adduct= T,
                                cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.ljw <- MSdev_get_Stat(msdev.ljw,QC_RSD = 0.3)
  MSdev_save(msdev.ljw)
  MSdev_export(msdev.ljw,candi = F)
  ### STAT
  {

    proj.dir <- msdev.ljw@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.ljw,from = "metabolite")


    data.se  <- data.se[,!data.se$sample.type=="Blank"]




    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### Sample_P
    data.se.Sample_P <- DEP_normalization(data.se)
    data.se.Sample_P <- data.se.Sample_P[,data.se.Sample_P$sample.type== "Sample"]
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 7,append = T)

    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)
    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff ,
                     paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 7,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )

    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    rowData(data.diff)$kegg.id <- rowData(data.diff)$kegg_id
    diff.path <- DEP_pathway_enrich(data.diff,
                                    contrast = "all")
    p.diff.path <- lapply(diff.path,
                          plotPathwayEnrichment)%>%
      ggplot_sum_patchwork()+
      plot_layout(ncol = 2)
    export_graph2pdf(p.diff.path , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 15,append = T)
    xlsx.write.list(diff.path,
                    paste0(proj.dir,"/Statistic/Pathway.enrich.xlsx"))
    p.hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(draw(p.hm) ,
                     paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 15,append = T)



    ### Heatmap
    data.selected <- lapply(1:3,function(x){
      x.df <- readxl::read_excel("d:/LJW_MS data20240611/Statistic/LJW.Labeled.xlsx",
                         sheet = x)
      x.df%>%
        dplyr::mutate(show = as.logical(show))%>%
        dplyr::filter(show)%>%
        dplyr::pull(protein)
    })%>%
      unlist()%>%
      unique()

    data.se <- data.se[data.selected,]
    rda <- rowData(data.se)
    heatmap.matrix <- assay(data.se)%>%t%>%scale%>%t

    Heatmap(heatmap.matrix,
            name = "Z score",
            column_split = data.se$condition,
            cluster_column_slices = F,
            cluster_columns = F,
            show_row_dend = F,
            row_names_side = "left",
            row_labels = rda$label)->p.hm

    export_graph2pdf(draw(p.hm,padding = unit(c(1,8,1,1),"cm")) ,
                     paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 15,height = 15,append = F)



  }




}





# Sat Jun 15 13:06:45 2024 CHB------------------------------
{
  msdev.CHB <- MSdev("d:/2024_06_13-Lirui/Data/" )
  #msdev.CHB <- load_as_var("d:/CHB_MS data20240611/MSdev_2024_06_12.Rdata")
  msdev.CHB <- MSdev_msConvert(msdev.CHB)
  msdev.CHB <- MSdev_checkSampleInfo(msdev.CHB)
  msdev.CHB <- MSdev_xcmsProcessing(msdev.CHB)
  msdev.CHB <- MSdev_annotation(msdev.CHB,
                                expand_adduct= T,
                                cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.CHB <- MSdev_get_Stat(msdev.CHB,QC_RSD = 0.3)
  MSdev_save(msdev.CHB)
  MSdev_export(msdev.CHB,candi = F)


}

# Thu Jun 20 10:20:17 2024 XCH------------------------------

{
  MSdev.XCH  <- MSdev("D:/2023_12_19-Xinchenhao/Data/")
  #MSdev.XCH <- load_as_var("D:/2024_05_20-Chenxin/MSdev_2024_05_23.Rdata")
  MSdev.XCH <- MSdev_msConvert(MSdev.XCH)
  MSdev.XCH <- MSdev_checkSampleInfo(MSdev.XCH)
  MSdev.XCH <- MSdev_xcmsProcessing(MSdev.XCH)
  MSdev.XCH <- MSdev_extract_Spectra(MSdev.XCH)
  MSdev.XCH <- MSdev_match_Spectra_to_feature(MSdev.XCH)
  MSdev.XCH <- MSdev_annotation(MSdev.XCH,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.XCH <- MSdev_get_Stat(MSdev.XCH,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.XCH,candi = F)
  MSdev_save(MSdev.XCH)

  {

    proj.dir <- MSdev.XCH@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.XCH,from = "metabolite")
    data.se <- data.se[,!colnames(data.se)%in% c( "WT_C1",   "WT_C2" ,  "WT_C3") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}




# Tue Jul  9 15:08:14 2024 CHB------------------------------
{
  msdev.CHB <- MSdev("d:/2024_07_08-Lirui/Data/" )
  #msdev.CHB <- load_as_var("d:/CHB_MS data20240611/MSdev_2024_06_12.Rdata")
  msdev.CHB <- MSdev_msConvert(msdev.CHB)
  msdev.CHB <- MSdev_checkSampleInfo(msdev.CHB)
  msdev.CHB <- MSdev_xcmsProcessing(msdev.CHB)
  msdev.CHB <- MSdev_annotation(msdev.CHB,
                                expand_adduct= T,
                                cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.CHB <- MSdev_get_Stat(msdev.CHB,QC_RSD = 0.3)
  MSdev_save(msdev.CHB)
  MSdev_export(msdev.CHB,candi = F)


}

# Wed Jul 31 15:46:54 2024 XCH ------------------------------
{
  MSdev.XCH  <- MSdev("D:/2024_07_25-Xinchenhao/Data/")
  MSdev.XCH <- load_as_var("D:/2024_07_25-Xinchenhao/MSdev_2024_07_31.Rdata")
  MSdev.XCH <- MSdev_msConvert(MSdev.XCH)
  MSdev.XCH <- MSdev_checkSampleInfo(MSdev.XCH)
  MSdev.XCH <- MSdev_xcmsProcessing(MSdev.XCH)
  MSdev.XCH <- MSdev_extract_Spectra(MSdev.XCH)
  MSdev.XCH <- MSdev_match_Spectra_to_feature(MSdev.XCH)
  MSdev.XCH <- MSdev_annotation(MSdev.XCH,
                                cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                                selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.XCH <- MSdev_get_Stat(MSdev.XCH,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.XCH,candi = F)
  MSdev_save(MSdev.XCH)

  {

    proj.dir <- MSdev.XCH@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.XCH,from = "metabolite")
    data.se <- data.se[,!colnames(data.se)%in% c( "WT_C1",   "WT_C2" ,  "WT_C3") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}


# Mon Aug  5 15:51:08 2024 CX ------------------------------
{
  MSdev.CX  <- MSdev("D:/2024_08_01-Chenxin/Data/")
  MSdev.CX <- load_as_var("D:/2024_08_01-Chenxin/MSdev_2024_08_05.Rdata")
  MSdev.CX <- MSdev_msConvert(MSdev.CX)
  MSdev.CX <- MSdev_checkSampleInfo(MSdev.CX)
  MSdev.CX <- MSdev_xcmsProcessing(MSdev.CX)
  MSdev.CX <- MSdev_extract_Spectra(MSdev.CX)
  MSdev.CX <- MSdev_match_Spectra_to_feature(MSdev.CX)
  MSdev.CX <- MSdev_annotation(MSdev.CX,
                                cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                                selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.CX <- MSdev_get_Stat(MSdev.CX,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.CX,candi = F)
  MSdev_save(MSdev.CX)

  {

    proj.dir <- MSdev.CX@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.CX,from = "metabolite")
    data.se <- data.se[,!colnames(data.se)%in% c( "G1",   "G2" ,  "G3","G4") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}

# Sat Aug 10 15:36:57 2024 CX------------------------------
{
  MSdev.CX  <- MSdev("D:/2024_08_07-Chenxin/Data/")
  MSdev.CX <- load_as_var("D:/2024_08_07-Chenxin/MSdev_2024_08_10.Rdata")
  MSdev.CX <- MSdev_msConvert(MSdev.CX)
  MSdev.CX <- MSdev_checkSampleInfo(MSdev.CX)
  MSdev.CX <- MSdev_update_xcms_pdata(MSdev.CX)
  MSdev.CX <- MSdev_xcmsProcessing(MSdev.CX)
  MSdev.CX <- MSdev_extract_Spectra(MSdev.CX)
  MSdev.CX <- MSdev_match_Spectra_to_feature(MSdev.CX)
  MSdev.CX <- MSdev_annotation(MSdev.CX,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.CX <- MSdev_get_Stat(MSdev.CX,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.CX,candi = F)
  MSdev_save(MSdev.CX)

  {

    proj.dir <- MSdev.CX@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.CX,from = "metabolite")
    data.se <- data.se[,!colnames(data.se)%in% c( "G1",   "G2" ,  "G3","G4") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}

# Mon Aug 12 21:55:13 2024 ZQ------------------------------
{
  MSdev.ZQ  <- MSdev("d:/2024_08_08-Zhaoqiang/Data/")
  MSdev.ZQ <- load_as_var("D:/2024_08_08-Zhaoqiang/MSdev_2024_08_12.Rdata")
  MSdev.ZQ <- MSdev_msConvert(MSdev.ZQ)
  MSdev.ZQ <- MSdev_checkSampleInfo(MSdev.ZQ)
  MSdev.ZQ <- MSdev_update_xcms_pdata(MSdev.ZQ)
  MSdev.ZQ <- MSdev_xcmsProcessing(MSdev.ZQ)
  MSdev.ZQ <- MSdev_extract_Spectra(MSdev.ZQ)
  MSdev.ZQ <- MSdev_match_Spectra_to_feature(MSdev.ZQ)
  MSdev.ZQ <- MSdev_annotation(MSdev.ZQ,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+",
                                                   "[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.ZQ <- MSdev_get_Stat(MSdev.ZQ,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.ZQ,candi = F)
  MSdev_save(MSdev.ZQ)



  ### ZQ
  {

    proj.dir <- MSdev.ZQ@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.ZQ,from = "metabolite")
    data.se <- data.se[,data.se$group%in% c( "G1",   "G2" ,  "G3","G4","QC","Blank") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



  ### ESCC
  {

    proj.dir <- MSdev.ZQ@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.ZQ,from = "metabolite")
    data.se <- data.se[,data.se$group%in% c( "ESCC_CON",   "ESCC_KO" ,"QC","Blank" ) ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }


}


# Thu Sep 19 13:54:07 2024 CHB------------------------------
{
  #msdev.CHB <- MSdev("d:/2024_07_08-Lirui/Data/" )
  msdev.CHB <- load_as_var("d:/2024.09.19.CHB/MSdev_2024_09_19.Rdata")
  msdev.CHB <- MSdev_msConvert(msdev.CHB)
  msdev.CHB <- MSdev_checkSampleInfo(msdev.CHB)
  msdev.CHB <- MSdev_xcmsProcessing(msdev.CHB)
  msdev.CHB <- MSdev_annotation(msdev.CHB,
                                expand_adduct= T,
                                selected_adduct = c("[M-H]-","[M+H]+"),
                                cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.CHB <- MSdev_get_Stat(msdev.CHB,QC_RSD = 0.3)
  MSdev_save(msdev.CHB)
  MSdev_export(msdev.CHB,candi = F)


}

# Wed Jul 31 15:46:54 2024 XCH ------------------------------
{
  MSdev.XCH  <- MSdev("D:/2024.09.25.XCH/Data/")
  MSdev.XCH <- load_as_var("D:/2024.09.25.XCH/MSdev_2024_09_27.Rdata")
  MSdev.XCH <- MSdev_msConvert(MSdev.XCH)
  MSdev.XCH <- MSdev_checkSampleInfo(MSdev.XCH)
  MSdev.XCH <- MSdev_xcmsProcessing(MSdev.XCH)
  MSdev.XCH <- MSdev_extract_Spectra(MSdev.XCH)
  MSdev.XCH <- MSdev_match_Spectra_to_feature(MSdev.XCH)
  MSdev.XCH <- MSdev_annotation(MSdev.XCH,
                                cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                                selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.XCH <- MSdev_get_Stat(MSdev.XCH,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.XCH,candi = F)
  MSdev_save(MSdev.XCH)

  {

    proj.dir <- MSdev.XCH@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.XCH,from = "metabolite")
    data.se <- data.se[,data.se$group%in% c( "G10",   "G50" ,  "G_B","QC") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}


# Mon Sep 30 13:49:14 2024 ------------------------------
{
  MSdev.CX  <- MSdev("D:/2024_09_27-Chenxin/Data/")
  MSdev.CX <- load_as_var("D:/2024_09_27-Chenxin/MSdev_2024_09_30.Rdata")
  MSdev.CX <- MSdev_msConvert(MSdev.CX)
  MSdev.CX <- MSdev_checkSampleInfo(MSdev.CX)
  MSdev.CX <- MSdev_update_xcms_pdata(MSdev.CX)
  MSdev.CX <- MSdev_xcmsProcessing(MSdev.CX)
  MSdev.CX <- MSdev_extract_Spectra(MSdev.CX)
  MSdev.CX <- MSdev_match_Spectra_to_feature(MSdev.CX)
  MSdev.CX <- MSdev_annotation(MSdev.CX,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.CX <- MSdev_get_Stat(MSdev.CX,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.CX,candi = F)
  MSdev_save(MSdev.CX)

  {

    proj.dir <- MSdev.CX@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.CX,from = "metabolite")
    #data.se <- data.se[,!colnames(data.se)%in% c( "G1",   "G2" ,  "G3","G4") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )



  }



}

# Fri Nov 22 13:08:48 2024 ZJ lipidomic------------------------------
{

  MSdev.ZJ  <- MSdev("D:/2024_11_19-Zhujing/Data/")
  MSdev.ZJ <- load_as_var("D:/2024_11_19-Zhujing/MSdev_2024_11_23.Rdata")
  MSdev.ZJ <- MSdev_msConvert(MSdev.ZJ)
  MSdev.ZJ <- MSdev_checkSampleInfo(MSdev.ZJ)
  MSdev.ZJ <- MSdev_update_xcms_pdata(MSdev.ZJ)
  MSdev.ZJ <- MSdev_xcmsProcessing(MSdev.ZJ)
  MSdev.ZJ <- MSdev_extract_Spectra(MSdev.ZJ)
  MSdev.ZJ <- MSdev_match_Spectra_to_feature(MSdev.ZJ)
  MSdev.ZJ <- MSdev_annotation(MSdev.ZJ,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.ZJ <- MSdev_get_Stat(MSdev.ZJ,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.ZJ,candi = F)
  MSdev_save(MSdev.ZJ)

  {

    proj.dir <- MSdev.ZJ@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.ZJ,from = "metabolite")
    #data.se <- data.se[,!colnames(data.se)%in% c( "G1",   "G2" ,  "G3","G4") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height =8,append = T)


  }



}
# Wed Dec  4 14:30:02 2024 ZJ lipidomic------------------------------
{

  MSdev.ZJ  <- MSdev("D:/2024_12_04-Zhujing/Data/")
  MSdev.ZJ <- load_as_var("D:/2024_12_04-Zhujing/MSdev_2024_12_04.Rdata")
  MSdev.ZJ <- MSdev_msConvert(MSdev.ZJ)
  MSdev.ZJ <- MSdev_checkSampleInfo(MSdev.ZJ)
  MSdev.ZJ <- MSdev_update_xcms_pdata(MSdev.ZJ)
  MSdev.ZJ <- MSdev_xcmsProcessing(MSdev.ZJ)
  MSdev.ZJ <- MSdev_extract_Spectra(MSdev.ZJ)
  MSdev.ZJ <- MSdev_match_Spectra_to_feature(MSdev.ZJ)
  MSdev.ZJ <- MSdev_annotation(MSdev.ZJ,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.ZJ <- MSdev_get_Stat(MSdev.ZJ,candi = F,QC_RSD = 10)
  MSdev_export(MSdev.ZJ,candi = F)
  MSdev_save(MSdev.ZJ)

  {

    proj.dir <- MSdev.ZJ@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.ZJ,from = "metabolite")
    #data.se <- data.se[,!colnames(data.se)%in% c( "G1",   "G2" ,  "G3","G4") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height =8,append = T)


  }



}

# Fri Dec 27 12:22:34 2024 ZQ lipidomic------------------------------
{

  MSdev.ZJ  <- MSdev("D:/ZQ/data/")
  MSdev.ZJ <- load_as_var("D:/ZQ/MSdev_2024_12_27.Rdata")
  MSdev.ZJ <- MSdev_msConvert(MSdev.ZJ)
  MSdev.ZJ <- MSdev_checkSampleInfo(MSdev.ZJ)
  MSdev.ZJ <- MSdev_update_xcms_pdata(MSdev.ZJ)
  MSdev.ZJ <- MSdev_xcmsProcessing(MSdev.ZJ)
  MSdev.ZJ <- MSdev_extract_Spectra(MSdev.ZJ)
  MSdev.ZJ <- MSdev_match_Spectra_to_feature(MSdev.ZJ)
  MSdev.ZJ <- MSdev_annotation(MSdev.ZJ,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  MSdev.ZJ <- MSdev_get_Stat(MSdev.ZJ,candi = F,polarity_paired = F,QC_RSD = 10)
  MSdev_export(MSdev.ZJ,candi = F)
  MSdev_save(MSdev.ZJ)

  ### Tumor
  {

    proj.dir <- MSdev.ZJ@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(MSdev.ZJ,from = "metabolite")
    data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 15,append = T)


  }



}

# Sat Dec 28 14:13:10 2024 ZJ Metabolomics------------------------------
{
  msdev.zj <- MSdev("d:/2024_12_28-Zhujing/data/")
  msdev.zj <- load_as_var("d:/2024_12_28-Zhujing/MSdev_2024_12_28.Rdata")
  msdev.zj <- MSdev_msConvert(msdev.zj)
  msdev.zj <- MSdev_checkSampleInfo(msdev.zj)
  msdev.zj <- MSdev_xcmsProcessing(msdev.zj)
  msdev.zj <- MSdev_extract_Spectra(msdev.zj)
  msdev.zj <- MSdev_match_Spectra_to_feature(msdev.zj)
  msdev.zj <- MSdev_annotation(msdev.zj,
                                expand_adduct= T,
                                cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.zj <- MSdev_get_Stat(msdev.zj,QC_RSD = Inf)
  MSdev_save(msdev.zj)
  MSdev_export(msdev.zj)



  {

    proj.dir <- msdev.zj@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.zj,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### sample_p
   # data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 15,append = T)

    data.path <- DEP_pathway_enrich(data.diff,contrast = "all")
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}
# Sun Feb  9 21:50:46 2025 LE------------------------------
{
  msdev.LE <- MSdev("d:/20250208.LE/2250125_WYJ_Metabolomics_15min_Folch//")
  msdev.LE <- load_as_var("d:/20250208.LE/MSdev_2025_02_08.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)



  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### sample_p
    # data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 15,append = T)

    data.path <- DEP_pathway_enrich(data.diff,contrast = "all")
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}

# Tue Feb 11 10:08:20 2025 LE Lipidomics------------------------------
{
  msdev.LE <- MSdev("d:/20250211.LE/20240122_WYJ_Lipidomics_17min_2-inject/data/")
  msdev.LE <- load_as_var("d:/20250211.LE/20240122_WYJ_Lipidomics_17min_2-inject/MSdev_2025_02_11.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)




  ### figure
  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff,p.adjust = T)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 120,append = T)


  }


}

# Tue Feb 11 10:08:20 2025 LE Lipidomics 2------------------------------
{
  msdev.LE <- MSdev("d:/20250211.LE/20250125_WYJ_Lipidomics_30min_P&N-Switch/data/")
  msdev.LE <- load_as_var("d:/20250211.LE/20250125_WYJ_Lipidomics_30min_P&N-Switch/MSdev_2025_02_11.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)



  ### figure
  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff,p.adjust = T)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 80,append = T)


  }


}
# Tue Feb 11 16:04:38 2025 LE Metabolomics -----------------------------
{
  msdev.LE <- MSdev("d:/20250211.LE/20250125_WYJ_Metabolomics_15min_Folch/data/")
  msdev.LE <- load_as_var("d:/20250211.LE/20250125_WYJ_Metabolomics_15min_Folch/MSdev_2025_02_11.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)



  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)

    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 60,append = T)


    data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
    p <- plotPathwayEnrichment(data.path[[1]],method = "bubble")
    export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 4,append = T)
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}
# Tue Feb 11 16:04:38 2025 LE Metabolomics 2 -----------------------------
{
  msdev.LE <- MSdev("d:/20250211.LE/20250125_WYJ_Metabolomics_15min_Routine/data/")
  msdev.LE <- load_as_var("d:/20250211.LE/20250125_WYJ_Metabolomics_15min_Routine/MSdev_2025_02_11.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)



  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)
    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 60,append = T)

    data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
    p <- plotPathwayEnrichment(data.path[[1]],method = "bubble")
    export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 4,append = T)
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}

# Tue Feb 25 23:45:20 2025 LE Metabolomics -----------------------------
{
  msdev.LE <- MSdev("d:/2025.02.25.LE/20250222_LHQ_Metabolite/20250222_LHQ_Metabolite/")
  msdev.LE <- load_as_var("d:/2025.02.25.LE/20250222_LHQ_Metabolite/MSdev_2025_02_25.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = 0.3)
  MSdev_save(msdev.LE)



  MSdev_export(msdev.LE)
  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    data.se <- data.se[,!data.se$label%in%  c("Meta__Con26")]
    data.se <- data.se[rowData(data.se)$score>0.5,]
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)

    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 60,append = T)


    data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
    p <- plotPathwayEnrichment(data.path[[1]],method = "bubble")
    export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 4,append = T)
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}

# Wed Feb 26 02:55:00 2025 LE Lipidomics------------------------------
{
  msdev.LE <- MSdev("d:/2025.02.25.LE/20250215_LHQ_Lipid/data/")
  msdev.LE <- load_as_var("d:/2025.02.25.LE/20250215_LHQ_Lipid/MSdev_2025_02_26.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_save(msdev.LE)




  MSdev_export(msdev.LE)
  ### figure
  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    data.se <- data.se[rowData(data.se)$score>-1,]
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff,p.adjust = T)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 120,append = T)


  }


}


# Fri Mar  7 20:34:07 2025 Lipidomics------------------------------
{
  msdev.LE <- MSdev("d:/20250306.LE/LXM_LIP/LXM_Lipid/")
  msdev.LE <- load_as_var("d:/20250306.LE/LXM_LIP/MSdev_2025_03_07.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)




  ### figure
  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff,p.adjust = T)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 120,append = T)


  }


}

# Fri Mar  7 20:34:07 2025  LE Lipidomics 2------------------------------
{
  msdev.LE <- MSdev("d:/20250306.LE/Zafar_Lipid/Zafar_Lipid/")
  msdev.LE <- load_as_var("d:/20250306.LE/Zafar_Lipid/MSdev_2025_03_07.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/Lipidblast.sqlite",
                               selected_adduct = c("[M]+","[M+NH4]+","[M+H]+","[M+Na]+","[M-H]-","[M+HCOO]-","[M+CH3COO]-"))
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)



  ### figure
  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }


    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff,p.adjust = T)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 80,append = T)


  }


}
# Fri Mar  7 20:34:07 2025 LE Metabolomics -----------------------------
{
  msdev.LE <- MSdev("d:/20250306.LE/LXM_META/LXM_Meta/")
  msdev.LE <- load_as_var("d:/20250306.LE/LXM_META/MSdev_2025_03_08.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)



  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)

    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 60,append = T)


    data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
    p.list <- lapply(names(data.path),
                     function(x){
                       p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                     })
    p <- ggplot_sum_patchwork(p.list)
    export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 20,height = 10,append = T)
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}
# Fri Mar  7 20:34:07 2025 LE Metabolomics 2 -----------------------------
{
  msdev.LE <- MSdev("d:/20250306.LE/Zafar_Meta/Zafar_Meta/")
  msdev.LE <- load_as_var("d:/20250306.LE/Zafar_Meta/MSdev_2025_03_08.Rdata")
  msdev.LE <- MSdev_msConvert(msdev.LE)
  msdev.LE <- MSdev_checkSampleInfo(msdev.LE)
  msdev.LE <- MSdev_xcmsProcessing(msdev.LE)
  msdev.LE <- MSdev_extract_Spectra(msdev.LE)
  msdev.LE <- MSdev_match_Spectra_to_feature(msdev.LE)
  msdev.LE <- MSdev_annotation(msdev.LE,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.LE <- MSdev_get_Stat(msdev.LE,QC_RSD = Inf)
  MSdev_export(msdev.LE)
  MSdev_save(msdev.LE)



  {

    proj.dir <- msdev.LE@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.LE,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)
    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.LE@xcmsData$PositiveMS1)+
        labs(title = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.LE@xcmsData$NegativeMS1)+
        labs(title = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### sample_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.Sample_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.Sample_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.Sample_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 60,append = T)

    data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
    p.list <- lapply(names(data.path),
                     function(x){
                       p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                     })
    p <- ggplot_sum_patchwork(p.list)
    export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 20,height = 10,append = T)
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}
# Thu Apr  3 14:25:40 2025 CityU NS------------------------------
{
  msdev.NS <- MSdev("d:/2025.04.02.CityU.NS/data/")
  msdev.NS <- MSdev_load("d:/2025.04.02.CityU.NS/MSdev_2025_04_04.Rdata")
  msdev.NS <- MSdev_msConvert(msdev.NS)
  msdev.NS <- MSdev_checkSampleInfo(msdev.NS)
  msdev.NS <- MSdev_update_xcms_pdata(msdev.NS)
  msdev.NS <- MSdev_set_param(msdev.NS)
  msdev.NS <- MSdev_xcmsProcessing(msdev.NS)
  msdev.NS <- MSdev_extract_Spectra(msdev.NS)
  msdev.NS <- MSdev_annotation(msdev.NS,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.NS <- MSdev_get_Stat(msdev.NS,QC_RSD = 0.3)
  MSdev_export(msdev.NS)
  MSdev_save(msdev.NS)



  {

    proj.dir <- msdev.NS@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.NS,from = "metabolite")
    data.se$condition <- data.se$sample.type
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 5,height = 5)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)
    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.NS@xcmsData$PositiveMS1)+
        labs(titNS = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.NS@xcmsData$NegativeMS1)+
        labs(titNS = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### APC
    {
      stp <- "APC"
      data.se <- get_MSdev_DEP_se(msdev.NS,from = "metabolite")
      data.se <- data.se[,data.se$sample.type%in%c(stp)]
      data.se.SampNS_P <- DEP_normalization(data.se)

      p.pca <- DEP_plot_PCA(data.se.SampNS_P)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5,append = T)


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      data.diff <- DEP_test_diff(data.se.SampNS_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      tabNS.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(tabNS.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.",stp,".xlsx")
      )


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
      hm <- DEP.plot.heatmap(data.diff)
      export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 6,height = 60,append = T)

      data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
      p.list <- lapply(names(data.path),
                       function(x){
                         p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                       })
      p <- ggplot_sum_patchwork(p.list)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 8,height = 8,append = T)
      xlsx.write.list(
        data.path,
        file =paste0(proj.dir,"/Statistic/pathway.",stp,".xlsx")
      )


    }

    ###   HT
    {
      stp <- "HT"
      data.se <- get_MSdev_DEP_se(msdev.NS,from = "metabolite")
      data.se <- data.se[,data.se$sample.type%in%c(stp)]
      data.se.SampNS_P <- DEP_normalization(data.se)

      p.pca <- DEP_plot_PCA(data.se.SampNS_P)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5,append = T)


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      data.diff <- DEP_test_diff(data.se.SampNS_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      tabNS.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(tabNS.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.",stp,".xlsx")
      )


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
      hm <- DEP.plot.heatmap(data.diff)
      export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 6,height = 60,append = T)

      data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
      p.list <- lapply(names(data.path),
                       function(x){
                         p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                       })
      p <- ggplot_sum_patchwork(p.list)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 8,height = 8,append = T)
      xlsx.write.list(
        data.path,
        file =paste0(proj.dir,"/Statistic/pathway.",stp,".xlsx")
      )


    }

    ###    Liver
    {
      stp <- "Liver"
      data.se <- get_MSdev_DEP_se(msdev.NS,from = "metabolite")
      data.se <- data.se[,data.se$sample.type%in%c(stp)]
      data.se.SampNS_P <- DEP_normalization(data.se)

      p.pca <- DEP_plot_PCA(data.se.SampNS_P)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5,append = T)


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      data.diff <- DEP_test_diff(data.se.SampNS_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      tabNS.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(tabNS.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.",stp,".xlsx")
      )


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
      hm <- DEP.plot.heatmap(data.diff)
      export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 6,height = 60,append = T)

      data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
      p.list <- lapply(names(data.path),
                       function(x){
                         p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                       })
      p <- ggplot_sum_patchwork(p.list)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 8,height = 8,append = T)
      xlsx.write.list(
        data.path,
        file =paste0(proj.dir,"/Statistic/pathway.",stp,".xlsx")
      )


    }

    ###   M
    {
      stp <- "M"
      data.se <- get_MSdev_DEP_se(msdev.NS,from = "metabolite")
      data.se <- data.se[,data.se$sample.type%in%c(stp)]
      data.se.SampNS_P <- DEP_normalization(data.se)

      p.pca <- DEP_plot_PCA(data.se.SampNS_P)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5,append = T)


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      data.diff <- DEP_test_diff(data.se.SampNS_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      tabNS.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(tabNS.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.",stp,".xlsx")
      )


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
      hm <- DEP.plot.heatmap(data.diff)
      export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 6,height = 60,append = T)

      data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
      p.list <- lapply(names(data.path),
                       function(x){
                         p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                       })
      p <- ggplot_sum_patchwork(p.list)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 8,height = 8,append = T)
      xlsx.write.list(
        data.path,
        file =paste0(proj.dir,"/Statistic/pathway.",stp,".xlsx")
      )


    }

    ### P
    {
      stp <- "P"
      data.se <- get_MSdev_DEP_se(msdev.NS,from = "metabolite")
      data.se <- data.se[,data.se$sample.type%in%c(stp)]
      data.se.SampNS_P <- DEP_normalization(data.se)

      p.pca <- DEP_plot_PCA(data.se.SampNS_P)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5,append = T)


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      data.diff <- DEP_test_diff(data.se.SampNS_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      tabNS.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(tabNS.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.",stp,".xlsx")
      )


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
      hm <- DEP.plot.heatmap(data.diff)
      export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 6,height = 60,append = T)

      data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
      p.list <- lapply(names(data.path),
                       function(x){
                         p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                       })
      p <- ggplot_sum_patchwork(p.list)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 8,height = 8,append = T)
      xlsx.write.list(
        data.path,
        file =paste0(proj.dir,"/Statistic/pathway.",stp,".xlsx")
      )


    }

    ###  SCF
    {
      stp <- "SCF"
      data.se <- get_MSdev_DEP_se(msdev.NS,from = "metabolite")
      data.se <- data.se[,data.se$sample.type%in%c(stp)]
      data.se.SampNS_P <- DEP_normalization(data.se)

      p.pca <- DEP_plot_PCA(data.se.SampNS_P)
      export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 5,height = 5,append = T)


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      data.diff <- DEP_test_diff(data.se.SampNS_P)
      data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

      p.diff.list <- DEP_plot_volcano(data.diff,"all")
      p.diff <- ggplot_sum_patchwork(p.diff.list)
      export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 10,append = T)
      tabNS.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
      xlsx.write.list(tabNS.diff,
                      paste0(proj.dir,"/Statistic/diff.metabolites.",stp,".xlsx")
      )


      data.diff <- DEP_test_diff(data.se.SampNS_P,type = "all")
      data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
      #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
      hm <- DEP.plot.heatmap(data.diff)
      export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 6,height = 60,append = T)

      data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
      p.list <- lapply(names(data.path),
                       function(x){
                         p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                       })
      p <- ggplot_sum_patchwork(p.list)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 8,height = 8,append = T)
      xlsx.write.list(
        data.path,
        file =paste0(proj.dir,"/Statistic/pathway.",stp,".xlsx")
      )


    }




  }


}

# Thu Jun  5 00:04:30 2025 XHL Metabolomics  -----------------------------
{
  msdev.XHL <- MSdev("d:/20250501/rawdata/")
  #msdev.XHL <- load_as_var("d:/20250306.XHL/Zafar_Meta/MSdev_2025_03_08.Rdata")
  msdev.XHL <- MSdev_msConvert(msdev.XHL,format.to = "mzXML")
  msdev.XHL <- MSdev_checkSampleInfo(msdev.XHL)
  msdev.XHL <- MSdev_xcmsProcessing(msdev.XHL)
  msdev.XHL <- MSdev_extract_Spectra(msdev.XHL)
  msdev.XHL <- MSdev_match_Spectra_to_feature(msdev.XHL)
  msdev.XHL <- MSdev_annotation(msdev.XHL,
                               expand_adduct= T,
                               cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.XHL <- MSdev_get_Stat(msdev.XHL,QC_RSD = Inf)
  MSdev_export(msdev.XHL)
  MSdev_save(msdev.XHL)



  {

    proj.dir <- msdev.XHL@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.XHL,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)
    ### Plot TIC
    {
     # p1 <- plot_xcms_TIC(msdev.XHL@xcmsData$PositiveMS1)+
     #   labs(titXHL = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.XHL@xcmsData$NegativeMS1)+
        labs(titXHL = "Negative TIC")
      p <- p2+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### sampXHL_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.SampXHL_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.SampXHL_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.SampXHL_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.SampXHL_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 60,append = T)

    data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
    p.list <- lapply(names(data.path),
                     function(x){
                       p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title = x)
                     })
    p <- ggplot_sum_patchwork(p.list)
    export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 20,height = 10,append = T)
    xlsx.write.list(
      data.path,
      file =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}

# Thu Jun  5 00:04:30 2025 XHL Metabolomics2  -----------------------------
{
  msdev.XHL <- MSdev("d:/20250519/rawdata/")
  msdev.XHL <- load_as_var("d:/20250519/MSdev_2025_06_05.Rdata")
  msdev.XHL <- MSdev_msConvert(msdev.XHL)
  msdev.XHL <- MSdev_checkSampleInfo(msdev.XHL)
  msdev.XHL <- MSdev_xcmsProcessing(msdev.XHL)
  msdev.XHL <- MSdev_extract_Spectra(msdev.XHL)
  msdev.XHL <- MSdev_match_Spectra_to_feature(msdev.XHL)
  msdev.XHL <- MSdev_annotation(
    msdev.XHL,
    expand_adduct= T,
    cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")
  msdev.XHL <- MSdev_get_Stat(msdev.XHL,
                              QC_RSD = Inf,rt_bin = 30,
                              score = 0.5)
  MSdev_export(msdev.XHL)
  MSdev_save(msdev.XHL)



  {

    proj.dir <- msdev.XHL@projectInfo$projectDir
    data.se <- get_MSdev_DEP_se(msdev.XHL,from = "metabolite")
    #data.se <- data.se[,grepl(data.se$group,pattern = "Tumor") ]
    p.pca <- DEP_plot_PCA(data.se)
    export_graph2pdf(p.pca , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 3,height = 3)


    #maped.genes <- KEGG_get_cp_linked_gene(rowData(data.se)$kegg_id)
    #edit_df_in_excel(maped.genes)
    ### Plot TIC
    {
      p1 <- plot_xcms_TIC(msdev.XHL@xcmsData$PositiveMS1)+
        labs(titXHL = "Positive TIC")
      p2 <- plot_xcms_TIC(msdev.XHL@xcmsData$NegativeMS1)+
        labs(titXHL = "Negative TIC")
      p <- p1+p2+plot_layout(ncol = 1)
      export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                       width = 10,height = 8,append = T)
    }



    ### sampXHL_p
    data.se <- data.se[,!data.se$condition%in%c("QC","Blank","A0")]
    data.se.SampXHL_P <- DEP_normalization(data.se)
    data.diff <- DEP_test_diff(data.se.SampXHL_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    data.diff <- DEP_test_diff(data.se.SampXHL_P)
    data.diff <- DEP_add_rejections(data.diff,p.adjust = T)

    p.diff.list <- DEP_plot_volcano(data.diff,"all")
    p.diff <- ggplot_sum_patchwork(p.diff.list)
    export_graph2pdf(p.diff , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 10,height = 10,append = T)
    table.diff <- DEP_get_diff_table(data.diff,contrast = "all",keep.all = T)
    xlsx.write.list(table.diff,
                    paste0(proj.dir,"/Statistic/diff.metabolites.xlsx")
    )


    data.diff <- DEP_test_diff(data.se.SampXHL_P,type = "all")
    data.diff <- DEP_add_rejections(data.diff,p.adjust = F)
    #data.diff <- data.diff[,grepl("LB|V|F",data.diff$group )]
    hm <- DEP.plot.heatmap(data.diff)
    export_graph2pdf(hm , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 6,height = 60,append = T)

    data.path <- DEP_pathway_enrich(data.diff,contrast = "all",method = "GlobalTest")
    p.list <- lapply(names(data.path),
                     function(x){
                       p <- plotPathwayEnrichment(data.path[[x]],method = "bubble",title  = x)
                     })
    p <- ggplot_sum_patchwork(p.list)
    export_graph2pdf(p , paste0(proj.dir,"/Statistic/Figures.pdf"),
                     width = 20,height = 10,append = T)
    xlsx.write.list(
      data.path,
      file  =paste0(proj.dir,"/Statistic/pathway.xlsx")
    )
  }


}

# Thu Jun 19 15:16:10 2025 DDA Mine pos------------------------------
{



  ### run after QC FS, once
  data.dir <- "d:/DDAmine/pos/data/"
  msdev.qe <- MSdev(rawDataDir = data.dir)
  msdev.qe <- MSdev_load("d:/DDAmine/pos/MSdev_2025_06_19.Rdata")
  msdev.qe <- MSdev_msConvert(msdev.qe)
  msdev.qe <- MSdev_xcmsProcessing(msdev.qe)

  msdev.qe@statData <- list()
  msdev.qe <- MSdev_get_Inclusion_Queue(msdev.qe)

  ### run after every time DDA acquired
  msdev.qe <- MSdev_get_Inclusion_List(msdev.qe)
  msdev.qe <- MSdev_add_sample(msdev.qe,
                               raw.data.dir = "d:/DDAmine/pos/data/")
  msdev.qe <- MSdev_get_MS2acquisitionStat(msdev.qe)

  table(msdev.qe@statData$DDA_mine_queue_Positive$acquired)




}


# Thu Jun 19 15:16:10 2025 DDA Mine neg------------------------------
{



  ### run after QC FS, once
  data.dir <- "d:/DDAmine/neg/data/"
  msdev.qe <- MSdev(rawDataDir = data.dir)
  msdev.qe <- MSdev_load("d:/DDAmine/neg/MSdev_2025_06_19.Rdata")
  msdev.qe <- MSdev_msConvert(msdev.qe)
  msdev.qe <- MSdev_xcmsProcessing(msdev.qe)

  msdev.qe@statData <- list()
  msdev.qe <- MSdev_get_Inclusion_Queue(msdev.qe)

  ### run after every time DDA acquired
  msdev.qe <- MSdev_get_Inclusion_List(msdev.qe)
  msdev.qe <- MSdev_add_sample(msdev.qe,raw.data.dir = "d:/DDAmine/neg/data/")
  msdev.qe <- MSdev_get_MS2acquisitionStat(msdev.qe)

  table(msdev.qe@statData$DDA_mine_queue_Negative$acquired)




}

# Thu Jun 19 22:25:32 2025 ------------------------------
{

  msdev.ddamine <- MSdev("d:/20250619_LR/data/")
  msdev.ddamine <- MSdev_load("d:/20250619_LR/MSdev_2025_06_19.Rdata")
  msdev.ddamine <- MSdev_msConvert(msdev.ddamine)
  msdev.ddamine <- MSdev_checkSampleInfo(msdev.ddamine)
  msdev.ddamine <- MSdev_xcmsProcessing(msdev.ddamine)
  msdev.ddamine <- MSdev_extract_Spectra(msdev.ddamine)
  msdev.ddamine <- MSdev_match_Spectra_to_feature(msdev.ddamine)
  msdev.ddamine <- MSdev_annotation(
    msdev.ddamine,
    expand_adduct= T,
    cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")

  msdev.ddamine <- MSdev_get_Stat(msdev.ddamine)
  MSdev_export(msdev.ddamine)
  MSdev_save(msdev.ddamine)



  msdev.demo <- MSdev_get_Stat(msdev.demo)
  plot_MSdev_normalization(msdev.demo)
  plot_MSdev_QC_RSD_CDF(msdev.demo)
  plot_MSdev_QC_RSD_hist(msdev.demo)




}
