# Fri Jun 28 21:17:17 2024 cfmd atom map prob------------------------------
cfmd <- CFM_annotate_by_predict()

get_CFM_data_trans_igraph(cfmd)%>%
  visIgraph()

cfmd <- CFM_data_get_igraph(cfmd)


vis_cfm_data_atom_map(cfmd,2)%>%
  open_visNet()

kegg.cp <- MSdb::get_CompoundDB_Compound()
kegg.pathway <- MSdb:::get_KEGG_compound_pathway_df()

kegg.pathway <- kegg.pathway%>%
  dplyr::filter(ENTRY=="hsa00510")

edit_df_in_excel(kegg.pathway)

# Wed Jul  3 21:30:25 2024 ------------------------------
smiles <- "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4[13CH2][C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)C"



cfmd <- CFM_annotate_by_predict(smiles)
cfmd <- CFM_data_get_igraph(cfmd)

load("chole.temp.rda")
get_CFM_data_trans_igraph(cfmd)%>%
  visIgraph()%>%
  open_visNet()


sim.sp <- get_CFM_data_Spectra(cfmd)
sim.sp.pd <- peaksData(sim.sp)

sim.sp.pd <- lapply(sim.sp.pd,function(z){
  z[,"mz"] <- z[,"mz"] +C13_mass_diff
  return(z)
})


sim.sp.pd -> Spectra::peaksData(sim.sp@backend)

cfmd <- cfm_data_get_fragment_group(cfmd)
msip.core <- get_MSIPCoreData(sim.sp,
                              cfmd,0,ppm = 10)

fg.map <- msip.core@FG_map %>%
    MSIPFragmentMap_reduce_fragment()
heatmap_MSIPFragmentMap(fg.map)->p
open_plot_win(p,10,10)


msdev.13C1@xcmsData$Negative_Chromatograms



xcms.xcms <- msdev.13C1@xcmsData$NegativeMS1
xcms.xcms <- xcms_remove_feature_var(xcms.xcms,
                                     grep(pattern = "Ratio_to_seed",value = T,
                                          x = colnames(featureDefinitions(xcms.xcms))))
xcms.xcms <- xcms_remove_feature_var(xcms.xcms,"ms1_purity")
colnames(featureDefinitions(xcms.xcms))
xcms.xcms -> msdev.13C1@xcmsData$NegativeMS1


chrom <- msdev.13C1@xcmsData$Positive_Chromatograms
chrom.data <- onDiskData_retrieve(chrom)
xcms.fdf<- featureDefinitions(chrom.data)
var <- grep(pattern = "Ratio_to_seed",value = T,
     x = colnames(xcms.fdf))
var.selected <- setdiff(colnames(xcms.fdf),var)
xcms.fdf <- xcms.fdf[,var.selected]
xcms.fdf -> chrom.data@featureDefinitions
msdev.13C1@xcmsData$Negative_Chromatograms <-
  onDiskData(chrom.data,chrom@path)



# Mon Jul  8 11:25:24 2024 ------------------------------
A <- MSdb::get_KEGG_compound_pathway_df()%>%
  dplyr::filter(ENTRY == 'hsa00020')
edit_df_in_excel(A)
cp <- MSdb::get_CompoundDB_Compound()
cp.tca <-cp %>% dplyr::filter(
  kegg_id %in% A$COMPOUND.ID)
edit_df_in_excel(cp.tca)

# Mon Jul  8 16:18:04 2024 ------------------------------


###  Get metabolites
{
  kegg.pathway <- MSdb:::get_KEGG_compound_pathway_df()%>%
    dplyr::filter(ENTRY == "hsa00020")
  kegg.compound <- MSdb:::get_KEGG_compound_df()
  compound.tca <- kegg.compound[kegg.compound$KEGG_id %in% kegg.pathway$COMPOUND.ID,]%>%
    as.data.frame()%>%
    dplyr::filter(!is.na(Exact_mass))



}



### pos
{
  xcms.pos <- msdev.purity@xcmsData$PositiveMS1
  xcms.fdf.pos <- get_xcms_feature_definitions(xcms.pos)
  tca.cp.pos <- chemform_adduct(compound.tca$Formula,"[M+H]+",value = "all")%>%
    dplyr::mutate(compound.tca[,c("KEGG_id","Name")])
  peak.match.pos <- match_mz_rt(mz1 = tca.cp.pos$chemform.adduct.mz,
                                mz2 = xcms.fdf.pos$mzmed,
                                mz.ppm = 10)%>%
    dplyr::mutate(tca.cp.pos[ion1,c("KEGG_id","Name")])


  p.list.pos <- list()
  for (i in seq_along(unique(peak.match.pos$ion1))) {
    ion <-unique(peak.match.pos$ion1)[i]
    plot_xcms_peaks_Chromatogram(
      xcms.pos,
      peak_id =peak.match.pos$ion2[peak.match.pos$ion1==ion]
    )->p
    p.list.pos[[i]] <- p+
      labs(title = paste0(tca.cp.pos$Name[ion],
                          "; ",tca.cp.pos$adduct[ion],
                          "; mz = ",tca.cp.pos$chemform.adduct.mz[ion]%>%
                            sprintf("%.4f",.)))

  }
  p.pos <- ggplot_sum_patchwork(p.list.pos)
  export_graph2pdf(p.pos,
                   file_path = paste0(proj.dir,"/Chrom.pos.pdf"),
                   width = 15,
                   height = 8)


}


### neg
{
  xcms.neg <- msdev.glu@xcmsData$NegativeMS1
  xcms.peaks.neg <- get_xcms_peaks_stat(xcms.neg)
  tca.cp.neg <- chemform_adduct(compound.tca$Formula,"[M-H]-")%>%
    dplyr::mutate(compound.tca[,c("KEGG_id","Name")])
  peak.match.neg <- match_mz_rt(mz1 = tca.cp.neg$chemform.adduct.mz,
                                mz2 = xcms.peaks.neg$mz,
                                mz.ppm = 5)%>%
    dplyr::mutate(tca.cp.neg[ion1,c("KEGG_id","Name")])


  p.list.neg <- list()
  for (i in seq_along(unique(peak.match.neg$ion1))) {
    ion <-unique(peak.match.neg$ion1)[i]
    plot_xcms_peaks_Chromatogram(
      xcms.neg,
      peak_id =peak.match.neg$ion2[peak.match.neg$ion1==ion]
    )->p
    p.list.neg[[i]] <- p+
      labs(title = paste0(tca.cp.neg$Name[ion],
                          "; ",tca.cp.neg$adduct[ion],
                          "; mz = ",tca.cp.neg$chemform.adduct.mz[ion]%>%
                            sprintf("%.4f",.)))

  }
  p.neg <- ggplot_sum_patchwork(p.list.neg)+
    plot_layout(ncol = 2)
  export_graph2pdf(p.neg,
                   file_path = paste0(proj.dir,"/Chrom.neg.pdf"),
                   width = 15,
                   height = 20)


}


### TCA adduct.mz
{

  tca.cp <- list("Positive"=tca.cp.pos,
                 "Negative"=tca.cp.neg)

  xlsx.write.list(tca.cp,
                  file = paste0(proj.dir,
                                "/TCA.metabolites.mz.xlsx"))

}

{
  roi <- readxl::read_excel("C:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/MSIP.interest.list.xlsx")
  roi <- roi %>%
    dplyr::mutate(id = match(name,a$name),
                  compound_id = case_when(!is.na(id)~a$compound_id[id],
                                          T~ compound_id),
                  formula = case_when(!is.na(id)~a$formula[id],
                                          T~ formula),
                  )
  edit_df_in_excel(roi)

}

MSIP_update_compoundDB_from_interest_list()





# Mon Jul  8 20:41:22 2024 CFM MAP------------------------------
msdev.purity <- load_as_var("C:/Users/91879/OneDrive/Code/R/data/MSIP_data/240601_ThreeGroup/MSdev_2024_06_04.Rdata")

iso.data <-msdev.purity@statData$MSIP$isotopologues_data[[21]]

msip.core <- get_MSIPCoreData(sp.iso = iso.data$Spectra$M3$U,
                              cfmd = iso.data$CFM_annotation,
                              iso_count = 3)
a <- MSIPCore_solve(msip.core)

# Thu Jul 11 13:35:04 2024 ESCC ------------------------------
gene <- "SLAMF6"
plot_single_feature(#data.se = analysis.data$diff.se$tissue.transcriptome$T1_vs_T0.based.on.All ,
                    data.se = tissue.transcriptome.se,
                    feature_id = gene,
                    diff.table = analysis.data$diff.feature$tissue.transcriptome$R_vs_NR.based.on.T1_T0
)->p
p
MSdev:::open_plot_win(p,3,4)

plot_single_feature(data.se = tissue.proteome.se,
                    feature_id = gene,
                    diff.table = analysis.data$diff.feature$tissue.proteomic$R_vs_NR.based.on.T1_T0
) -> p
p
MSdev:::open_plot_win(p,3,4)

# Thu Jul 11 20:15:41 2024 ------------------------------
load("chole.temp.rda")
get_cfm_data_sdf_igraph(cfmd)%>%
  vis_sdf_igraph()


vis_cfm_data_atom_map(cfmd,100)%>%
  open_visNet()
#错误: BiocParallel errors
#1 remote errors, element index: 34
#34 unevaluated and other errors
#first remote error:
#  Error in quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 0): 外接函数调用时不能有NA/NaN/Inf(arg1)



