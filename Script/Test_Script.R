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


CFM_data_get_igraph()


smiles2sdf("C1=NC2=C(N1)N=CN=C2N")%>%
  plot

vis_smiles("C1=NC2=C(N1)N=CN=C2N")%>%
  open_visNet()

cfmd <- CFM_annotate_by_predict("C1=NC2=C(N1)N=CN=C2N")
cfmd <- cfm_data_get_fragment_group(cfmd)
cfmd <- CFM_data_get_igraph(cfmd)


vis_cfm_data_atom_map(cfmd ,59,show.id = T)%>%
  open_visNet()

vis_cfm_data_fragment(cfmd,3)%>%
  open_visNet()

get_cfm_data_fragment_group_atom_map(cfmd,"2")

### debug for fragment59


# Sat Jul 13 14:22:20 2024 adenine------------------------------
smiles <- "C1=NC2=C(N1)N=CN=C2N"
sdf <- get_smiles_sdf(smiles)
{
  sdf.igraph <- get_sdf_igraph(sdf)[[1]]

  sdf.igraph%>%
    sdf_igraph_add_border_color(value = 0.99)%>%
    sdf_igraph_add_background_color(value = 0.8)%>%
    vis_sdf_igraph()%>%
    open_visNet()

}

{
  sdf.igraph1 <- sdf.igraph%>%
    sdf_igraph_add_background_color(0.8)
  sdf.igraph2 <- sdf.igraph%>%
    sdf_igraph_add_border_color(0.6)
  sdf_igraph_merge(sdf.igraph1,
                   sdf.igraph2)%>%
    vis_sdf_igraph()%>%
    open_visNet()

}


cfmd <- CFM_data_get_igraph(cfmd)
atom.map <- get_atom_map(sdf.parent ,sdf.product ,ig.parent ,ig.product )
vis_sdf_igraph(ig.parent,T)%>%open_visNet()
vis_sdf_igraph(ig.product,T)%>%open_visNet()

apply(atom.map,2,sum)

vis_cfm_data_trans_map(cfmd,32 ,show_id = T)%>%
  visOptions(width = "150%")%>%
  open_visNet()

cfmd <- CFM_annotate_by_predict()


trans.maps <- sapply(1:115,get_CFM_data_trans_map,cfmd = cfmd)

atoms <- sapply(trans.maps,function(x){
 sum( colSums(x)==0)
})


which(atoms>0)

# Sat Jul 13 20:14:51 2024 CFM map for NAD------------------------------
smiles <- "NC(=O)C1=C[N+](=CC=C1)[C@@H]1O[C@H](COP([O-])(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)N2C=NC3=C2N=CN=C3N)[C@@H](O)[C@H]1O"
cfmd <- CFM_annotate_by_predict(smiles)
cfmd <- cfm_data_get_fragment_group(cfmd)
cfmd <- CFM_data_get_igraph(cfmd)

trans.maps <- sapply(1:10,
                     get_CFM_data_trans_map,
                     cfmd = cfmd)

vis_cfm_data_trans_map(cfmd,1)%>%
  open_visNet()




# Sat Jul 13 21:19:20 2024 adenine ------------------------------
smiles <- "C1=NC2=C(N1)N=CN=C2N"
cfmd <- CFM_annotate_by_predict(smiles)
cfmd <- CFM_data_get_igraph(cfmd)
cfmd <- cfm_data_get_fragment_group(cfmd)
cfmd.trans <- check_CFM_data_trans_map(cfmd)

vis_cfm_data_trans_map(cfmd,1,show_id = F)%>%
  open_visNet()

debug(get_atom_map)
get_CFM_data_trans_map(cfmd,1)


# Sun Jul 14 15:12:51 2024 Glu------------------------------
{
  smiles <- "C(C1C(C(C(C(O1)O)O)O)O)O"
  cfmd <- CFM_annotate_by_predict(smiles)
  cfmd <- CFM_data_get_igraph(cfmd)
  cfmd <- cfm_data_get_fragment_group(cfmd)
  cfmd.trans <- check_CFM_data_trans_map(cfmd)

  vis_cfm_data_trans_map(object,90,show_id = F)%>%
    open_visNet()

  debug(get_atom_map)
  a <- get_CFM_data_trans_map(object,5)

  rowSums(a)
  colSums(a)

  a <- CFM_data_get_atom_map(cfmd,BPPARAM = SnowParam(workers = 4,
                                                      progressbar = T))

  vis_cfm_data_fragment_atom_map(a,"Fragment204")%>%
    open_visNet()
  vis_cfm_data_trans_map(a,70,show_id = F)%>%
    open_visNet()
}

# Mon Jul 15 11:30:58 2024 ------------------------------
{
  iso.data <- msdev.purity@statData$MSIP$isotopologues_data
  a <- sapply(iso.data,function(x){
    nrow(x$CFM_annotation@fragment_transition)
  })
  sapply(a,function(x){
    sum(is.na(x))
  })/ lengths(a)

  which(lengths(a)==0)%>%names()%>%
    vector2str()%>%cat()
}


a <- distances(ig.trans,1,mode = "out")

# Tue Jul 16 20:37:11 2024 ------------------------------
chemform_adduct("C6[13]C0H12O5",adduct = "[M-H]-")%>%format(digit=10)


# Wed Jul 17 14:16:47 2024 CP_DB for CFM_data------------------------------
cpdbt <- MSdb::get_CompoundDB_Compound()
cpdbt <- cpdbt[1:5,]%>%
  DataFrame()
cfmds <- sapply(cpdbt$smiles,CFM_annotate_by_predict)
cpdbt$cfmd <- unname(cfmds)

cpdb.new <- insertCompound(cpdb.new,cpdbt,
                           addColumns   = T)

b <- compounds(cpdb.new)
compoundVariables(cpdb.new)



get_CFM_data_from_smiles(smiles = iso.data$FT11035_Positive $compound_info$smiles,
                         temp_dir = "d:/ttttmp")

msdev.purity <- MSIP_get_isotopologues_CFM_annotation(msdev.purity,
                                                      BPPARAM = SerialParam(progressbar = T))


file.name <- rep(letters,10)%>%paste0(collapse = "")

saveRDS(a,file = file.name)


cpdbt <- MSdb::get_CompoundDB_Compound()

cpdbt$smiles%>%
  sapply(nchar)->a

get_CFM_data_from_smiles(smiles = iso.data$FT00684_Negative$compound_info$smiles,
                         temp_dir = "d:/ttttmp")
# Fri Jul 19 16:21:44 2024 ------------------------------
{
  # Define the probabilities
  probs <- sample(seq(0,1,0.1),30,replace = T)


}



system.time(get_iso_prob_chatgpt(fc,ifc)[1:7])
system.time(get_iso_prob(fc,ifc)[1:7])
# Fri Jul 19 18:28:25 2024 ------------------------------
iso.data <- msdev.purity@statData$MSIP$isotopologues_data
x <- iso.data$FT11035_Positive
sapply(iso.data,function(x){

  n_c <- chemform_parse(x$compound_info$formula)[,"C"]
  n_iso <- x$Spectra%>%names%>%str_isotope2_num()%>%max
  message_with_time(x$compound_info$compound_id)
  choose(n_c,n_iso)

  #x$compound_info$compound_id
}) ->a



# Wed Jul 24 11:06:13 2024 MSIP Solve------------------------------
a <- MSIP_solve_test(msdev.purity,
                     BPPARAM = SnowParam(workers = 6,
                                         progressbar = T))





# Wed Jul 24 15:09:44 2024 ------------------------------
{

  my_long_running_function <- function(data) {
    Sys.sleep(40)  # Simulating a function that takes 40 seconds to complete
    return(paste("Result for", data))
  }

  library(BiocParallel)
  library(parallel)

  # Example large list-like data
  input_list <- list(
    element1 = "data1",
    element2 = "data2",
    element3 = "data3"
    # Add more elements as needed...
  )

  # Initialize BiocParallel with SnowParam
  bp_param <- SnowParam(workers = 2, timeout = 30)  # Set timeout to 30 seconds

  # Wrapper function to handle timeout
  timeout_wrapper <- function(data) {
    result <- tryCatch({
      my_long_running_function(data)
    }, error = function(e) {
      message("Timeout or error occurred for", data, ":", e$message)
      return("Timeout or error occurred")
    })
    return(result)
  }

  # Use BiocParallel's bplapply to apply the timeout wrapper function
  timeout_results <- bplapply(input_list, timeout_wrapper, BPPARAM = bp_param)

  # Print the results
  print(timeout_results)

}


rm(a)
a <- R.utils::withTimeout(my_long_running_function(1) ,
                     timeout = 3,
                     onTimeout = "silent")

my_long_running_function <- function(data) {
  Sys.sleep(data)  # Simulating a function that takes 40 seconds to complete
  return(paste("Result for", data))
}

a <- bplapply(1:10,
         BPPARAM = SerialParam(progressbar = T),
         FUN = function(i){
           R.utils::withTimeout(my_long_running_function(i) ,
                       timeout = 3,
                       onTimeout = "silent")
})



