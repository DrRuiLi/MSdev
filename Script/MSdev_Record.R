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
                                #selected_adduct = c("[M-H]-",
                                #                    "[M-H2O-H]-",
                                #                    "[2M-H]-" ,
                                #                    "[M+FA-H]-" ,
                                #                    "[M+H]+" ,
                                #                    "[M-H2O+H]+",
                                #                    "[M+NH4]+"
                                #) ,
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
  #msdev.ljw <- load_as_var("d:/2024_06_04-Wangyongqiang/MSdev_2024_06_07.Rdata")
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

  }




}
