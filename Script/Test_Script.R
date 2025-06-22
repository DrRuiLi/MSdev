# Fri Dec 27 13:15:52 2024 ------------------------------
gsh.smiles <- "C(CC(=O)N[C@@H](CS)C(=O)NCC(=O)O)[C@@H](C(=O)O)N"
Molecule_igraph <- get_Molecule_igraph_from_smiles(gsh.smiles)

vis_Molecule_igraph(Molecule_igraph,show_id = T)


pubchem.retrive <- webchem::pc_prop(vdata(ig)$PubChem )
vdata(ig)$smiles <- pubchem.retrive$IsomericSMILES
saveRDS(ig,file = paste0("temp/",str_date(),"Glucose_glutamate.ig.rds"))


ig <- readRDS("temp/20241227Glucose_glutamate.ig.rds")
mfn <- new("Metabolic_flux_network",
    metabolic_network = igraph::simplify(ig))

a <- get_Molecule_igraph_from_smiles(vdata(mfn)$smiles)

vdata(mfn)$mig <- a


# Sat Dec 28 14:01:26 2024 ------------------------------
f.neg <- dir("d:/2024_12_28-Zhujing/data/neg/",full.names = T)
f.new <- paste0(dirname(f.neg),"/",paste0("neg_",basename(f.neg)))
file.rename(f.neg,f.new)


# Sat Dec 28 23:00:41 2024 ------------------------------
gsh.smiles <- "C(CC(=O)N[C@@H](CS)C(=O)NCC(=O)O)[C@@H](C(=O)O)N"
Molecule_igraph <- get_Molecule_igraph_from_smiles(gsh.smiles)

mol.igs <- get_Molecule_igraph_from_smiles(vdata(mfn)$smiles)
vdata(mfn)$Molecule_igraph <- mol.igs


edata(mfn)[1,]
mol.ig.from <- V(mfn)[["C00022"]]$Molecule_igraph
mol.ig.to <- V(mfn)[["C00036"]]$Molecule_igraph

a <- get_atom_map(sdf.parent = mol.ig.from@sdf,sdf.product = mol.ig.to@sdf,
             ig.parent = mol.ig.from@igraph,
             ig.product = mol.ig.to@igraph)

get_Molecule_atom_transfer_by_map(mol.ig.from,
                                  mol.ig.to)

mfn.ig <- readRDS("temp/20241227Glucose_glutamate.ig.rds")

Metabolic_flux_network <- new("Metabolic_flux_network")

Metabolic_flux_network@metabolic_network <- mfn.ig


mfn <- Metabolic_flux_network_get_atom_transfer(Metabolic_flux_network)

edata(Metabolic_flux_network)$x <- mfn.transfer
#edata(Metabolic_flux_network)$x <- 1:68
edata(Metabolic_flux_network)$xxxxxxxx <- mfn.transfer


a <- edge.attributes(Metabolic_flux_network@metabolic_network)

# Mon Dec 30 20:40:53 2024 ------------------------------
mfn.ig <- readRDS("temp/20241227Glucose_glutamate.ig.rds")
mol.igs <- get_Molecule_igraph_from_smiles(vdata(mfn.ig)$smiles)
vdata(mfn.ig)$Molecule_igraph <- mol.igs
#mfn.ig <- igraph::simplify(mfn.ig)

Metabolic_flux_network <- new("Metabolic_flux_network")
Metabolic_flux_network@metabolic_network <- mfn.ig
Metabolic_flux_network <- Metabolic_flux_network_get_atom_transfer(Metabolic_flux_network)


edata(Metabolic_flux_network)$atom_transfer

saveRDS(Metabolic_flux_network,
        file = "temp/20250106Glucose_glutamate.MFN.rds")

# Mon Jan  6 10:53:47 2025 ------------------------------
Metabolic_flux_network <- readRDS(
        file = "temp/20250106Glucose_glutamate.MFN.rds")


vda <- vdata(Metabolic_flux_network)
glu.id <- "C00267"
glu.smiles <- "CC1=NCCC2=C1NC3=C2C=CC(=C3)O"
glu.mi <- get_Molecule_igraph_from_smiles(glu.smiles)

glu.mi <- Molecule_igraph_add_isotopomer(glu.mi,
                                         isotopomer = "tracer",
                                         iso_vec = "all_C")
glu.mi@isotopomer

# Tue Jan  7 17:27:54 2025 ------------------------------
glu.gsh.net <- "temp/20241227Glucose_glutamate.ig.xlsx"
glu.gsh.net <- igraph_import(glu.gsh.net)
cid <- webchem::get_cid(vdata(glu.gsh.net)$PubChem,
                        from = "sid",domain = "substance")
pubchem.retrive <- webchem::pc_prop(cid$cid)
vdata(glu.gsh.net)$smiles <- pubchem.retrive$IsomericSMILES

a <- vdata(glu.gsh.net)%>%
  dplyr::mutate(smiles.formula = get_smile_formula(smiles))

saveRDS(glu.gsh.net,file = "temp/20241227Glucose_glutamate.ig.rds")


mfn.ig <- readRDS("temp/20241227Glucose_glutamate.ig.rds")
mol.igs <- get_Molecule_igraph_from_smiles(vdata(mfn.ig)$smiles)
vdata(mfn.ig)$Molecule_igraph <- mol.igs
#mfn.ig <- igraph::simplify(mfn.ig)

Metabolic_flux_network <- new("Metabolic_flux_network")
Metabolic_flux_network@metabolic_network <- mfn.ig
Metabolic_flux_network <- Metabolic_flux_network_get_atom_transfer(Metabolic_flux_network)


edata(Metabolic_flux_network)$atom_transfer

saveRDS(Metabolic_flux_network,
        file = "temp/20250106Glucose_glutamate.MFN.rds")


vda <- vdata(Metabolic_flux_network)
vda <- vdata(Metabolic_flux_network)
glu.id <- "C00267"
glu.smiles <- "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
glu.mi <- get_Molecule_igraph_from_smiles(glu.smiles)

glu.mi <- Molecule_igraph_add_isotopomer(glu.mi,
                                         isotopomer = "tracer",
                                         iso_vec = "all_C")
vis_Molecule_igraph(glu.mi)
glu.mi@isotopomer

# Wed Jan  8 13:45:59 2025 ------------------------------
Metabolic_flux_network <- readRDS("temp/20250106Glucose_glutamate.MFN.rds")


test.xlsx <- tempfile(fileext = ".xlsx")
test.data <- rep(list(adducts),1000)
names(test.data) <- 1:1000
xlsx.write.list(test.data,  test.xlsx)
open_file(test.xlsx)

# Thu Jan  9 22:32:37 2025 ------------------------------
a <- c(glu.mi,glu.mi,glu.mi,glu.mi)
am <- matrix(a,nrow = 2)
Glutamate <- Molecule_igraph_matrix[17,]
Glutamate <- Molecule_igraph_matrix[16,]

p.pie.list <- list()
data.list <- list()
for (i.sample  in names(Glutamate)) {

  mol.ig <- Glutamate[[i.sample]]
  isotopomer.data <- mol.ig@isotopomer%>%
    dplyr::arrange(isotopologue)%>%
    dplyr::mutate(abundance = abundance/sum(abundance),
                  label = factor(label,levels = label))
  data.list[[i.sample]] <- isotopomer.data


  plot.data <- isotopomer.data%>%
    dplyr::slice_max(abundance,n = 10,with_ties = F)

  p <- ggplot(plot.data)+
    geom_bar(aes( x= 1, y = abundance,fill = label,group = isotopologue),stat = "identity")+
    ggsci::scale_fill_npg()+
    labs( title = sub("FS_",i.sample,replacement = ""))+
    coord_polar(theta = "y")+
    theme_void()+
    theme(plot.title = element_text(hjust = 0.5))
  p
  p.pie.list[[i.sample]] <- p

}

p <- ggplot_sum_patchwork(p.pie.list)
open_plot_win(p,16,9)


# Wed Jan 15 10:07:51 2025 Atom transfer shiny------------------------------
{

  Metabolic_flux_network <- readRDS("temp/20250106Glucose_glutamate.MFN.rds")
  MFN_manul_Shiny(Metabolic_flux_network)
  Metabolic_flux_network <- Metabolic_flux_network_get_atom_transfer(Metabolic_flux_network)


  #edata(Metabolic_flux_network)$id<- edata(Metabolic_flux_network)$name
  MFN_manul_Shiny(Metabolic_flux_network)

}


# Thu Jan 16 00:09:58 2025 ------------------------------
Metabolic_flux_network <- load_MFN(name = "GLN")
MFN_manul_Shiny(Metabolic_flux_network)


visGetEdges <- readRDS("d:/temp/edges_20250121150706.rds")



# Thu Jan 23 01:27:19 2025 	5. Initialize a network from KEGG  ------------------------------
{
  kegg.net <- get_KEGG_Reaction_network()
  kegg.pathway.df <- MSdb:::get_KEGG_compound_pathway_df()
  kegg.selected <- kegg.pathway.df %>%
    dplyr::filter(ENTRY %in% c("hsa00010",### glycolysis
                               "hsa00020",### TCA
                               "hsa00250",### glutamate
                               "hsa00480" ### GSH
    ))%>%
    dplyr::pull(COMPOUND.ID)%>%
    unique()%>%
    setdiff(c("C00005","C00006","C00014"))
  kegg.net.selected <- igraph_filter_vertex(
    kegg.net,names(V(kegg.net)) %in% kegg.selected|V(kegg.net)$node.type=="Reaction"
  )

  kegg.net.selected.vda <- vdata(kegg.net.selected)
  cp.node <- which(kegg.net.selected.vda$node.type=="Compound")
  kegg.net.selected.cid <- webchem::get_cid(kegg.net.selected.vda$PubChem[cp.node],
                                            from = "sid",domain = "substance")
  pubchem.retrive <- webchem::pc_prop(kegg.net.selected.cid$cid)
  vdata(kegg.net.selected)$smiles <- NA
  vdata(kegg.net.selected)$smiles[cp.node] <- pubchem.retrive$IsomericSMILES

 # kegg.net.selected.filter <- igraph_filter_vertex(
 #   kegg.net.selected,get_formula_ele_count(vdata(kegg.net.selected)$Formula,"C")>0
 # )
  kegg.net.selected.filter <- kegg.net.selected
  vdata(kegg.net.selected.filter)$Molecule_igraph  <- NA
  vdata(kegg.net.selected.filter)$Molecule_igraph[cp.node] <- get_Molecule_igraph_from_smiles(
    pubchem.retrive$IsomericSMILES  )

  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net.selected.filter)
  mfn <- Metabolic_flux_network_get_atom_transfer(mfn)
  edata(mfn)$id <- edata(mfn)$name
  MFN_manul_Shiny(mfn)
}


# Tue Feb  4 13:53:59 2025 start shiny------------------------------
Metabolic_flux_network <- load_MFN(name = "GSH")
MFN_manul_Shiny(Metabolic_flux_network)


# Wed Feb  5 13:12:22 2025 set tracer ------------------------------
{
  mfn.v <- vdata(Metabolic_flux_network)
  glucose.smiles <- "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
  glucose.mig <- get_Molecule_igraph_from_smiles(glucose.smiles)

  Glu_1_2.mig <- Molecule_igraph_add_isotopomer(Molecule_igraph = glucose.mig,
                                 isotopomer = "Tracer",
                                 iso_vec = c("C_6" = "[13]C","C_10" = "[13]C") ,
                                 abundance = 1)
  Glu_1_2.mig <- Molecule_igraph_remove_isotopomer(Glu_1_2.mig,"base")

  Metabolic_flux_network <- Metabolic_flux_network_set_tracer(Metabolic_flux_network,
                            "C00267",Glu_1_2.mig)



  MFN_manul_Shiny(Metabolic_flux_network)
}

# Wed Feb  5 15:33:51 2025 reaction tracing ------------------------------
{


  edata <- edata(Metabolic_flux_network)
  edge.selected <- 78
  mat <- edata$atom_transfer[[edge.selected]]
  mol.ig.from <-  V(Metabolic_flux_network@metabolic_network)[[edata$from[edge.selected]]]$Molecule_igraph
  mol.ig.to <-  V(Metabolic_flux_network@metabolic_network)[[edata$to[edge.selected]]]$Molecule_igraph



}




# Wed Feb  5 20:29:27 2025 Initialize a network from KEGG, reaction as node------------------------------
{
  kegg.net <- get_KEGG_Reaction_network()
  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net)
  kegg.pathway.df <- MSdb:::get_KEGG_compound_pathway_df()
  kegg.selected <- kegg.pathway.df %>%
    dplyr::filter(ENTRY %in% c("hsa00010",### glycolysis
                               "hsa00020",### TCA
                               "hsa00250",### glutamate
                               "hsa00480" ### GSH
    ))%>%
    dplyr::pull(COMPOUND.ID)%>%
    unique()%>%
    setdiff(c("C00005","C00006","C00014","C00080"))%>%
    intersect(names(V(kegg.net)))
  mfn <- Metabolic_flux_network_select_compound(mfn,kegg.selected)
  mfn <-  Metabolic_flux_network_get_compound_data_from_cid(mfn)
  mfn <- Metabolic_flux_network_clean_reactions(mfn)
  mfn <- Metabolic_flux_network_get_Reaction_atom_transfer(mfn)
  edata(mfn)$id <- edata(mfn)$name
  MFN_manul_Shiny(mfn)

}
# Mon Feb 10 10:33:09 2025 RXNMAPPER ------------------------------
{
  library()
  rxnmp <- import("rxnmapper")
  rxn_mapper = rxnmp$RXNMapper()
  rxns = c('CC(C)S.CN(C)C=O.Fc1cccnc1F.O=C([O-])[O-].[K+].[K+]>>CC(C)Sc1ncccc1F',
           'C>>C')
  x <- rxn_mapper$get_attention_guided_atom_maps(rxns,detailed_output=T)

  'CC(C)S.CN(C)C=O.Fc1cccnc1F.O=C([O-])[O-].[K+].[K+]>>CC(C)Sc1ncccc1F'%>%
    get_Molecule_igraph_from_smiles()%>%
    vis_Molecule_igraph()->a



  mfn <- load_MFN()
  edata(mfn@metabolic_network)[1,c("from","to","equation")]

  x <- rep(from,times = equation.coef[from])
  from.smiles <- V(mfn@metabolic_network)[x]$smiles
  from.smiles <- paste0(from.smiles,collapse = ".")


  x <- rep(to,times = equation.coef[to])
  to.smiles <- V(mfn@metabolic_network)[x]$smiles
  to.smiles <- paste0(to.smiles,collapse = ".")

  req <- paste0(from.smiles, ">>",to.smiles)
  rxns <- c(req,"C>>C")
  x <- rxn_mapper$get_attention_guided_atom_maps(rxns)
  x

}

# Mon Feb 10 15:41:02 2025 demo------------------------------
{
  rid <- "R00258"
  #kegg.net <- get_KEGG_Reaction_network()
  kegg.net <- readRDS("temp/kegg.net.rds")
  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net)
  kegg.selected <- edata(mfn) %>%
    dplyr::filter(REACTION_id == rid
    )
  kegg.selected <- c(kegg.selected$from,kegg.selected$to)%>%
    unique()%>%
    setdiff(c("C00005","C00006","C00014","C00080",rid))%>%
    intersect(names(V(kegg.net)))
  mfn <- Metabolic_flux_network_select_compound(mfn,kegg.selected)
  mfn <-  Metabolic_flux_network_get_compound_data_from_cid(mfn)


  vda <- vdata(mfn)%>%
    dplyr::filter(id == rid
    )

  from.id <- c("C00041","C00026")
  to.id <- c("C00022","C00025")

  from.smiles <- V(mfn@metabolic_network)[from.id]$smiles
  to.smiles <- V(mfn@metabolic_network)[to.id]$smiles


  from.string <- paste0(from.smiles,collapse = ".")
  to.string <- paste0(to.smiles,collapse = ".")

  req <- paste0(from.string, ">>",to.string)
  rxns <- c(req,"C>>C")


  rxn.result <- RXNMapper(rxns)[[1]]
  rxn.result.split <- str_split(rxn.result$mapped_rxn,">>")[[1]]%>%
    sapply(function(x){
      str_split(x,"\\.")
    })


  smi.rxn <- rxn.result.split[[2]][2]
  smi.rxn
  #smi.rxn%>%get_smiles_sdf(canonicalize = F)%>%sdf2smiles()


  can = T

  to.smiles[1]%>%
    get_smiles_sdf(canonicalize = can)%>%
    get_Molecule_igraph_from_sdf()%>%
    `[[`(1)%>%
    atom()

  smi.rxn%>%
    get_smiles_sdf(canonicalize = can)%>%
    get_Molecule_igraph_from_sdf()%>%
    `[[`(1)%>%
    atom()


}
# Mon Feb 10 19:37:01 2025 DEMO------------------------------
{
  from.smiles <- c("C[C@@H](C(=O)O)N","C(CC(=O)O)C(=O)C(=O)O")
  to.smiles <- c("CC(=O)C(=O)O","C(CC(=O)O)[C@@H](C(=O)O)N")

  rxn.map <- RXNMapper_map(from.smiles = from.smiles,
                           to.smiles = to.smiles)

}
# Tue Feb 11 10:57:07 2025 implement RXN into Metabolic flux network------------------------------
{


  kegg.net <- get_KEGG_Reaction_network()
  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net)
  kegg.pathway.df <- MSdb:::get_KEGG_compound_pathway_df()
  kegg.selected <- kegg.pathway.df %>%
    dplyr::filter(ENTRY %in% c("hsa00010",### glycolysis
                               "hsa00020",### TCA
                               "hsa00250",### glutamate
                               "hsa00480" ### GSH
    ))%>%
    dplyr::pull(COMPOUND.ID)%>%
    unique()%>%
    setdiff(c("C00005","C00006","C00014","C00080"))%>%
    intersect(names(V(kegg.net)))
  mfn <- Metabolic_flux_network_select_compound(mfn,kegg.selected)
  mfn <-  Metabolic_flux_network_get_compound_data_from_cid(mfn)
  mfn <- Metabolic_flux_network_clean_reactions(mfn)
  mfn <- Metabolic_flux_network_get_Reaction_atom_transfer(mfn)


  MFN_manul_Shiny(mfn)


}

# Tue Feb 11 19:06:33 2025 ------------------------------
{
  from.smiles[1]%>%
    vis_Molecule_igraph_smiles()

  from.smiles[1]%>%
    sub(pattern = "\\.","-",.)%>%
    canonicalize_smiles()%>%
    vis_Molecule_igraph_smiles()


}

# Wed Feb 12 14:55:01 2025 shiny update------------------------------
{

  MFN_manul_Shiny(mfn)

}


# Wed Feb 12 15:58:50 2025 tracing------------------------------
{
  glucose.smiles <- "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
  glucose.mig <- get_Molecule_igraph_from_smiles(glucose.smiles)

  Glu_1_2.mig <- Molecule_igraph_add_isotopomer(Molecule_igraph = glucose.mig,
                                                isotopomer = "Tracer",
                                                iso_vec = c("C_6" = "[13]C","C_10" = "[13]C") ,
                                                abundance = 1)
  Glu_1_2.mig <- Molecule_igraph_remove_isotopomer(Glu_1_2.mig,"base")


  #mfn <- load_MFN()
  mfn <- Metabolic_flux_network_set_tracer(mfn,
                                           "C00267",Glu_1_2.mig)

  mfn <- Metabolic_flux_tracing(mfn)
  MFN_manul_Shiny(mfn)


}

# Thu Feb 13 21:30:13 2025 a simple MFN for test ------------------------------
{

  kegg.net <- get_KEGG_Reaction_network()
  #kegg.net <- readRDS("temp/kegg.net.rds")
  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net)
  kegg.pathway.df <- MSdb:::get_KEGG_compound_pathway_df()
  kegg.selected <- kegg.pathway.df %>%
    dplyr::filter(ENTRY %in% c("hsa00010",### glycolysis
                               "hsa00020"### TCA
    ))%>%
    dplyr::pull(COMPOUND.ID)%>%
    unique()%>%
    intersect(names(V(kegg.net)))

  mfn <- Metabolic_flux_network_select_compound(mfn,kegg.selected)
  mfn <- Metabolic_flux_network_clean_reactions(mfn)
  mfn <-  Metabolic_flux_network_get_compound_data_from_cid(mfn)
  mfn <- Metabolic_flux_network_get_Reaction_atom_transfer(mfn)
  #mfn <- Metabolic_flux_tracing(mfn)
  #"R01063" %in% vdata(mfn)$id

  MFN_manul_Shiny(mfn)

}
# Fri Feb 14 02:05:56 2025 ------------------------------

rxns.bak -> rxns
rxns <- gsub(pattern = "#","",rxns)
rxns
rxn.result <- RXNMapper(rxns,detailed_output = T)[[1]]

# Fri Feb 14 10:36:33 2025 TCA MFN TEST------------------------------
{
  glucose.smiles <- "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
  glucose.mig <- get_Molecule_igraph_from_smiles(glucose.smiles)

  Glu_1_2.mig <- Molecule_igraph_add_isotopomer(Molecule_igraph = glucose.mig,
                                                isotopomer = "Tracer",
                                                iso_vec = c("C_6" = "[13]C","C_10" = "[13]C") ,
                                                abundance = 1)
  Glu_1_2.mig <- Molecule_igraph_remove_isotopomer(Glu_1_2.mig,"base")


  #mfn <- load_MFN()
  mfn <- Metabolic_flux_network_set_tracer(mfn,
                                           "C00267",Glu_1_2.mig)

  mfn <- Metabolic_flux_tracing(mfn)
  MFN_manul_Shiny(mfn)


}
# Fri Feb 14 13:41:57 2025 selected reaction ------------------------------
{
 # kegg.net <- get_KEGG_Reaction_network()
  kegg.net <- readRDS("temp/kegg.net.no.filter.hsa.rds")
  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net)
  rid <- link.pathway.reaction[names(link.pathway.reaction )%in% c("path:map00010","path:map00020")]%>%
    sub("rn:","",x=.)
  mfn <- Metabolic_flux_network_filter_reactions(mfn,rid)
  mfn <- Metabolic_flux_network_get_compound_data_from_cid(mfn)
  mfn <- Metabolic_flux_network_get_Reaction_atom_transfer(mfn)
  mfn <- Metabolic_flux_network_set_tracer(mfn,
                                           "C00267",Glu_1_2.mig)
  mfn <- Metabolic_flux_tracing(mfn)

  MFN_manul_Shiny(mfn)

}

# Wed Feb 19 19:16:03 2025 ------------------------------
{

  plot.data <- msdev.13C1@statData[["MSIP"]][["isotopologues_table"]][["Positive"]]
  plot.data <- plot.data%>%
    dplyr::mutate(peakwidth = rtmax-rtmin,
                  peakwidth = case_when(peakwidth>50~50,
                                        T~peakwidth))%>%
    dplyr::filter(is_labeled)

  ggplot(plot.data)+
    geom_point(aes(x = rtmed,y= mzmed,
                   size = peakwidth),pch= 21,
               fill = "#4DBBD5",col = "black",alpha = 0.5)+
    labs(x = "Retention time",y = "m/z",size = "Peak width")+
    theme_bw()->p1


  p2 <-  ggplot(plot.data[rep(1:nrow(plot.data),5),])+
    geom_histogram(aes(x = rtmed),binwidth = 5)+
    labs(x ="Retention time", y = "Scan count in 5s window")+
    theme_bw()

  p <- p1/p2+
    plot_layout(guides = "collect")
  open_plot_win(p,8,6)

}

# Thu Feb 20 14:45:34 2025 combined isotopomers for GSH------------------------------
{


  GSH <- Molecule_igraph_matrix[40,]

  p.pie.list <- list()
  data.list <- list()
  for (i.sample  in names(GSH)[2:3]) {

    mol.ig <- GSH[[i.sample]]
    isotopomer.data <- mol.ig@isotopomer%>%
      dplyr::arrange(isotopologue)%>%
      dplyr::mutate(abundance = abundance/sum(abundance),
                    label = factor(label,levels = label))
    data.list[[i.sample]] <- isotopomer.data


    plot.data <- isotopomer.data%>%
      dplyr::slice_max(abundance,n = 10,with_ties = F)

    p <- ggplot(plot.data)+
      geom_bar(aes( x= 1, y = abundance,fill = label,group = isotopologue),stat = "identity")+
      ggsci::scale_fill_npg()+
      labs( title = gsub("[FS_]",i.sample,replacement = " "),fill = "Isotopomers")+
      coord_polar(theta = "y")+
      theme_void()+
      theme(plot.title = element_text(hjust = 0.5))
    p
    p.pie.list[[i.sample]] <- p

  }

  p <- ggplot_sum_patchwork(p.pie.list)
  open_plot_win(p,16,9)
  export_graph2pdf(p,"d:/temp/temp.pdf",
                   width = 16,height= 9)


}

# Thu Feb 20 19:08:47 2025 MSIP STAT TABLE------------------------------
{
  EVA <- MSIP_solve_computation_evaluate(msdev.M1)%>%
    dplyr::filter(FSIS.count>0)


  EVA%>%
    dplyr::filter(samples=="Con")%>%
    dplyr::pull(feature_id)%>%
    unique()%>%
    length()

  EVA%>%
    dplyr::filter(samples=="Liver")%>%
    nrow()

  EVA%>%
    dplyr::filter(samples=="Liver")%>%
    dplyr::pull(FSIS.count)%>%
    sum(na.rm = T)


  cp.table <- get_MSIP_compound_info(msdev.M1,vars= "all")
  cp.df <- MSdb:::get_CompoundDB_Compound()
  cp.table$kegg.id <-  cp.df$kegg_id[ match(cp.table$compound_id, cp.df$compound_id)]

  path.table <- analyzePathwayHyperTest(cp.table$kegg.id)
  plot.data <- path.table%>%
    dplyr::mutate(ratio = Hit/Total,
                  label = case_when(ratio >0.23|Hit>10~pathway.name))%>%
    dplyr::filter(grepl(x = pathway.class,pattern = "Metabolism"))

  ggplot(plot.data)+
    geom_point(aes(x = Hit, y = ratio ,
                   fill = ratio,size = ratio),pch = 21)+
    ggrepel::geom_text_repel(aes(x = Hit, y = ratio ,label = label),
                             nudge_y  = 0.1 , size =3) +
    scico::scale_fill_scico(direction = 1,palette = "vikO")+
    scale_size(range = c(0,10))+
    labs(x = "Number of labled metabolites", y = "Ratio")+
    theme_classic()+
    theme(legend.position = "none")->p
  p
  open_plot_win(p,5,4)
}



# Thu Feb 20 22:34:55 2025 MFN  reaction------------------------------
{

  mfn <- load_MFN(name = "TCA_20250214153408")
  MFN_manul_Shiny(mfn)
  mfn.filter <- Metabolic_flux_network_filter_reactions(mfn,"R00014")
  MFN_manul_Shiny(mfn)


  vis_Reaction_atom_transfer(rat)%>%
    #visOptions(width = "200%")%>%
    open_visNet()


  mfn <- load_MFN(name = "TCA_20250214153408")
  mfn <- Metabolic_flux_remove_tracing(mfn)
  MFN_manul_Shiny(mfn)
  Metabolic_flux_tracing(mfn)


  plot.data <- labeled.stat%>%
    dplyr::filter(labeled.cp>0)
  ggplot(plot.data)+
    geom_bar(aes(x = round, y = labeled.cp ),stat = "identity",alpha = 0.8,
             color = "black",fill = "#4DBBD5",width = 0.8)+
    geom_point(aes(x = round, y = isotopomers.count/35),alpha = 1,
               col = "#E64B35",size = 3)+
    geom_line(aes(x = round, y = isotopomers.count/35),alpha = 0.8,
              col = "#E64B35",size = 1)+
    scale_y_continuous(name = "Count of labeled metabolites",sec.axis = sec_axis(~.*35,name = "Count of isotopomers"))+
    labs(x = "Round")+
    theme_bw()->p1
  p1
  open_plot_win(p1,6,3)


  mfn.c <- vdata(mfn)%>%
    dplyr::filter(node.type == "Compound")

  cp.stat <- data.frame(
    id = mfn.c$id,
    isotopomer.count = NA,
    total.count = NA
  )
  for (i in 1:nrow(mfn.c)) {

    formula <- mfn.c$Formula[i]
    c.count <- get_formula_ele_count(formula)

    cp.stat$isotopomer.count[i] <- nrow(mfn.c$Molecule_igraph[[i]]@isotopomer)
    cp.stat$total.count[i] <- sum(choose(c.count,0:c.count))

  }



  plot.data <- cp.stat%>%
    dplyr::mutate(
      isotopomer.count = case_when(isotopomer.count > total.count~total.count,
                                   T~isotopomer.count),
      ratio = isotopomer.count/total.count)

  ggplot(plot.data)+
    geom_jitter(aes(x =isotopomer.count, y =  log10(total.count) ,
                   size = ratio,fill = ratio),stroke = 1,
               alpha = 0.5,pch = 21)+
    scale_size(range = c(3,10))+
    scico::scale_fill_scico(direction = -1,palette = "batlow")+
    guides(
      fill = guide_legend(title = "Isotopomers\nRatio"),
      size = guide_legend(title = "Isotopomers\nRatio")
    )+
    labs(x = "Count of detected isotopmers",
         y = "Log10 count of theoritical isotopomers")+
    theme_bw()->p

  open_plot_win(p)


}


# Fri Feb 21 19:03:41 2025 ------------------------------
{


  gln.mol.ig <- Molecule_igraph_matrix[17,]




}

# Fri Feb 21 19:43:00 2025 MFN to glutamate------------------------------
{
  kegg.net <- readRDS("temp/kegg.net.rds")
  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net)
  link.pathway.reaction <- KEGGREST::keggLink("reaction","pathway")
  rid <- link.pathway.reaction[names(link.pathway.reaction )%in% c("path:map00010","path:map00020","path:map00250")]%>%
    sub("rn:","",x=.)
  mfn <- Metabolic_flux_network_filter_reactions(mfn,rid)
  mfn <- Metabolic_flux_network_get_compound_data_from_cid(mfn)
  mfn <- Metabolic_flux_network_get_Reaction_atom_transfer(mfn)



  ###
  {
    glucose.smiles <- "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
    glucose.mig <- get_Molecule_igraph_from_smiles(glucose.smiles)

    Glu_1_2.mig <- Molecule_igraph_add_isotopomer(Molecule_igraph = glucose.mig,
                                                  isotopomer = "Tracer",
                                                  iso_vec = c("C_6" = "[13]C","C_10" = "[13]C") ,
                                                  abundance = 1)
    Glu_1_2.mig <- Molecule_igraph_remove_isotopomer(Glu_1_2.mig,"base")


    #mfn <- load_MFN()
    mfn <- Metabolic_flux_remove_tracing(mfn)
    mfn <- Metabolic_flux_network_set_tracer(mfn,
                                             "C00267",Glu_1_2.mig)

    mfn <- Metabolic_flux_tracing(mfn)
  }

  MFN_manul_Shiny(mfn)


}

# Sun Feb 23 15:22:20 2025 ------------------------------
{
  time.df <- a%>%
    pivot_longer(2:4)%>%
    dplyr::mutate(times = as.character(times))

  ggplot(time.df)+
    geom_bar(aes(x = times, y = value,fill = name),
             position = "dodge",
             stat = "identity")+
    ggsci::scale_fill_nejm()+
    labs(x = "Count of calulation",
         y = "Time consume",
         fill = "Calculation\ntype")+
    theme_bw()->p

  open_plot_win(p,4,3)


}
# Tue Feb 25 11:42:20 2025 ------------------------------
{
  citation("MSnbase")


}

# Thu Mar  6 16:49:44 2025 ------------------------------
{
  mfn <- load_MFN()


}


# Wed Mar 12 14:51:41 2025 shiny------------------------------
{
  library(shiny)

  ui <- fluidPage(
    titlePanel("My Shiny App"),
    mainPanel(
      tags$iframe(
        src = "https://drruili.github.io/MSIP/index.html",
        width = "100%",
        height = "600px",
        style = "border:none;"
      )
    )
  )

  server <- function(input, output, session) {}

  shinyApp(ui, server)

}

{
  library(shiny)

  ui <- fluidPage(
    titlePanel("My Shiny App"),
    mainPanel(
      actionButton("show_modal", "Open GitHub Page")
    )
  )

  server <- function(input, output, session) {
    observeEvent(input$show_modal, {
      showModal(
        modalDialog(
          title = "GitHub Page",
          HTML('<iframe src="https://drruili.github.io/MSIP/index.html" width="100%" height="500px" style="border:none;"></iframe>'),
          easyClose = TRUE,
          footer = modalButton("Close"),
          size = "l"  # Large modal for better viewing
        )
      )
    })
  }

  shinyApp(ui, server)

}
{
  library(shiny)

  ui <- fluidPage(
    titlePanel("My Shiny App"),
    mainPanel(
      p("Welcome to my Shiny app!")
    ),

    # JavaScript to open a new pop-up window
    tags$script(HTML("
    Shiny.addCustomMessageHandler('openPopup', function(url) {
      window.open(url, '_blank', 'width=1000,height=600,scrollbars=yes,resizable=yes');
    });
  "))
  )

  server <- function(input, output, session) {
    # Trigger the pop-up when the app starts
    session$sendCustomMessage("openPopup", "https://drruili.github.io/MSIP/index.html")
  }

  shinyApp(ui, server)

}

# Wed Mar 12 16:38:20 2025 Glu structure------------------------------
{

  glucose.smiles <- "C(C1C(C(C(C(O1)O)O)O)O)O"
  glucose.mig <- get_Molecule_igraph_from_smiles(glucose.smiles)
  glucose.mig <- Molecule_igraph_add_isotopomer(glucose.mig,
                                 iso_vec = c("C_6"="[13]C"))
  glucose.mig <- Molecule_igraph_add_isotopomer(glucose.mig,
                                                iso_vec = c("C_10"="[13]C"))
  glucose.mig <- Molecule_igraph_add_isotopomer(glucose.mig,
                                                iso_vec = c("C_8"="[13]C"))
  vis_Molecule_igraph_isotopomer(glucose.mig,isotopomer = 4,show_id = F)



}

# Wed Mar 12 22:28:08 2025 ------------------------------
{

  h2 <- ComplexHeatmap::Heatmap(frag.ratio.matrix,
                                na_col  ="#999999",
                                cell_fun = cellfun,
                                # width = unit( ncol(frag.ratio.matrix)*length.unit,cell.unit),
                                # height =unit( nrow(frag.atom.matrix)*length.unit,cell.unit),
                                name = "Isotope labeled\nratio",
                                col = circlize::colorRamp2(breaks = c(0,0.5,1),
                                                           c("white","#F7844F","#B20C26")),
                                right_annotation  = rowAnnotation(
                                  intensity = ComplexHeatmap::anno_numeric(c(6,5),
                                                                           bg_gp = gpar(fill = "#AFAFAF", col = "black")),
                                  width  = unit(0.8,"inch"),
                                  annotation_label = list(intensity = "Log10\nIntensity"),
                                  annotation_name_rot  = 0,
                                  annotation_name_side  = "top"),
                                cluster_columns = F,
                                row_names_side  = "left",
                                column_names_side = "top",
                                column_names_rot = 0.5,
                                column_names_centered = T,
                                rect_gp =  grid::gpar(lwd=2,col = "black"),
                                cluster_rows = F)

  open_plot_win(h2)
}

# Fri Mar 14 14:29:57 2025 ------------------------------
a <- r_bg(func = function(){

  msdev.13C1 <- MSdev::load_as_var(
    "C:/Users/91879/OneDrive/Code/R/data/MSIP_data/240701_FS_ONE_POSITION/MSdev_2024_07_04.Rdata"
  )
  MSIP_shiny_start(msdev.13C1,port = 8303)
})


# Run an R script using system2() function


# Tue Mar 18 12:55:21 2025 RXNMapper------------------------------
{


  mfn <- load_MFN()
  MFN_manul_Shiny(mfn)


}
# Wed Mar 19 00:05:45 2025 Diffusion button------------------------------
{
  kegg.net <- get_KEGG_Reaction_network()
  mfn <- new("Metabolic_flux_network",metabolic_network = kegg.net)
  kegg.pathway.df <- MSdb:::get_KEGG_compound_pathway_df()
  kegg.selected <- kegg.pathway.df %>%
    dplyr::filter(ENTRY %in% c("hsa00010",### glycolysis
                               "hsa00020",### TCA
                               "hsa00250",### glutamate
                               "hsa00480" ### GSH
    ))%>%
    dplyr::pull(COMPOUND.ID)%>%
    unique()%>%
    setdiff(c("C00005","C00006","C00014","C00080"))%>%
    intersect(names(V(kegg.net)))
  mfn <- Metabolic_flux_network_select_compound(mfn,kegg.selected)
  mfn <-  Metabolic_flux_network_get_compound_data_from_cid(mfn)

  mfn <- Metabolic_flux_network_clean_reactions(mfn)
  mfn <- Metabolic_flux_network_get_Reaction_atom_transfer(mfn)
  edata(mfn)$id <- edata(mfn)$name
  MFN_manul_Shiny(mfn)




}





# Tue Apr  1 01:21:44 2025 ------------------------------
{

  migs <- list()
  for (i.smi in na.omit( cp.node.data$smiles)) {

    migs[[i.smi]] <- get_Molecule_igraph_from_smiles(i.smi)


  }

}



# Thu Apr 10 23:04:40 2025 ------------------------------
{
  library(ComplexHeatmap)
  library(tibble)
  data <- read.csv("d:/tmp/Heatmap-c57-HNF4a-metabolites.csv")%>%
    column_to_rownames("X")%>%
    t%>%scale%>%t
  Heatmap(data,
          name = "Z score",
          rect_gp = gpar(color = "white"),
          cluster_rows = F,
          row_names_side = "left",
          row_dend_side = "right")
  export::graph2pdf(file = "d:/tmp/heatmap.c57.metabolites.pdf",
                    width = 6,height = 5)


  data <- read.csv("d:/tmp/Heatmap-data-HNF4a_cre-Linifanib.csv")%>%
    column_to_rownames("X")%>%
    t%>%scale%>%t
  Heatmap(data,
          name = "Z score",
          cluster_rows = F,
          rect_gp = gpar(color = "white"),
          row_names_side = "left",
          row_dend_side = "right")
  export::graph2pdf(file = "d:/tmp/heatmap.Linifanib.pdf",
                    width = 6,height = 5)


}

# Mon Apr 14 16:21:06 2025 2025 combined isotopomers for GSH------------------------------
{


  GSH <- Molecule_igraph_matrix[40,]

  p.pie.list <- list()
  data.list <- list()
  for (i.sample  in names(GSH)[2:3]) {

    mol.ig <- GSH[[i.sample]]
    isotopomer.data <- mol.ig@isotopomer%>%
      dplyr::arrange(isotopologue)%>%
      dplyr::mutate(abundance = abundance/sum(abundance),
                    label = factor(label,levels = label))
    data.list[[i.sample]] <- isotopomer.data


    split.col <- function(col,idx){

      if(length(idx)==1) return(col)
      col <- col[1]
      col_fun <- circlize::colorRamp2(breaks = c(-1,max(idx)),
                           colors = c("white",col))
      col_fun(idx)
    }

    col.isotopologues <- make_vector(
      ggsci::pal_npg()(10)[seq_along(unique(isotopomer.data$isotopologue))],
      unique(isotopomer.data$isotopologue)
    )
    col.isotopologues["M0"] <- "#7E6148"
    plot.data <- isotopomer.data%>%
      dplyr::select(isotopomer,isotopologue,label,abundance)%>%
      dplyr::slice_max(abundance,n = 10,with_ties = F)%>%
      dplyr::mutate(abundance = abundance/sum(abundance))%>%
      dplyr::group_by(isotopologue)%>%
      dplyr::mutate(id = paste0(isotopologue,num2str(1:n())),
                    col = col.isotopologues[isotopologue],
                    col = split.col(col,1:n()))

    isotpologues.data <-      plot.data%>%
      dplyr::group_by(isotopologue)%>%
      dplyr::summarise(abundance = sum(abundance))
  #$data.split <- data.frame(
  #$  isotopologue = unique(plot.data$isotopologue),
  #$  id = paste0(unique(plot.data$isotopologue),"0"),
  #$  abundance = 0.02,
  #$  label = " "
  #$)
  #$plot.data <- bind_rows(plot.data,data.split)%>%
  #$  dplyr::arrange(id)
#$
    p <- ggplot()+
      geom_bar(aes( x= 1, y = abundance,
                    #col = isotopologue,
                    fill = label,group = isotopologue),
               data = plot.data,
               col = "white",
               linewidth = 1,width = 0.95,
               stat = "identity")+
      geom_bar(aes( x= 1.7, y = abundance,
                    #col = isotopologue,
                    fill = isotopologue,group = isotopologue),
               col = "black",
               data = isotpologues.data,
               linewidth = 1,width = 0.2,
               stat = "identity")+
     # ggsci::scale_fill_npg()+
      scale_fill_manual(values = c(
        make_vector(plot.data$col,plot.data$label),
        col.isotopologues))+
     # scale_color_manual(values = c("#EEEEEE","#888888","#222222"))+
      labs( title = gsub("[FS_]",i.sample,replacement = " "),fill = "Isotopomers")+
      xlim(c(0.3,2))+
      coord_polar(theta = "y")+
      theme_void()+
      theme(plot.title = element_text(hjust = 0.5))
    p
    p.pie.list[[i.sample]] <- p

  }

  p <- ggplot_sum_patchwork(p.pie.list)
  open_plot_win(p,16,9)
  export::graph2ppt(p,"d:/temp.pdf",
                   width = 16,height= 9)


}




# Wed Apr 16 21:26:06 2025 ------------------------------
{

  wyq.file <- "d:/temp/20230228-3NPH-非靶向-ZScore(1).xlsx"
  data <- readxl::read_excel(wyq.file)
  hmdb.df <- MSdb:::get_HMDB_Compound_DF()
  kegg.df <- MSdb:::get_KEGG_compound_df()
  hmdb.df$chemical_formula[match(data$`Database ID`,hmdb.df$kegg_id)]
  data.formula<- data %>%
    dplyr::mutate(
      Formula = case_when(
        grepl("HMDB",`Database ID`) ~ hmdb.df$chemical_formula[match(`Database ID`,hmdb.df$accession)],
        grepl("C",`Database ID`) ~ kegg.df$Formula[match(`Database ID`,kegg.df$KEGG_id)]
      ),
      Formula.dev = MSCC::chemform_calc( Formula,"C6H5N3O1" ,return = "chemform"),
      polarity = str_extract(data$ID...1,"^[NP]"),
      rt = sapply(str_split(data$ID...1,"[\\|_]"),`[`,2),
      .after = Formula
    )

  openxlsx::write.xlsx(data.formula,file = wyq.file)


}





# Tue Apr 22 15:14:31 2025 word cloud------------------------------
{

  library(pdftools)

  text <- pdf_text("c:/Users/91879/OneDrive/Documents/ShanghaiTech/2025毕业/明审与答辩/明审材料/盲审意见修改/2.毕业论文_基于质谱的代谢物同位素异构体解析算法与代谢组学应用.pdf")
  text <- paste(text, collapse = " ")


  library(jiebaR)
  whitel.list <- c("同位素异构体","同位素同系物","代谢组学",
                   "代谢流","原子映射","天然同位素","MSdev",
                   "isotopomer","isotopologue","","","","","","","","",
                   "MSIP","PD-1",
                   "C13","13C")
  tpf <- tempfile()
  writeLines(whitel.list, tpf)

  cutter <- worker(user = tpf)  # 创建分词器
  words <- segment(text, cutter)

  # 过滤掉无用字符、停用词（如有）
  words <- words[nchar(words) > 1]  # 只保留长度 >1 的词（防止很多"的"、"了"）
  # 将所有英文统一为小写
  words <- textstem::lemmatize_words(words)

  # 去除复数（简单处理 -s 结尾）
  words <- gsub("s$", "", words)
  ### black list
  black.list <- c("进行","过程","可能","提供","能够","不同","具有",
                  "the","to","Figure","","","","","","","","","","","","",
                  "at","of","et","for","al","and","in")
  words <- words[!words %in% black.list]
  # 统计词频
  stop_words <- readLines("https://raw.githubusercontent.com/stopwords-iso/stopwords-zh/master/stopwords-zh.txt", encoding = "UTF-8")
  words <- words[!words %in% stop_words]


  df <- as.data.frame(table(words))
  df <- df[order(df$Freq, decreasing = TRUE), ]

  library(wordcloud2)


  # 随机浅色
 #wordcloud2(df, color = "random-light")
  # 固定某种颜色（红色）
  wordcloud2(df, size = 0.6, #minSize = 1,
             shape = "circle",ellipticity = 0.3,
             fontFamily = "微软雅黑", color = "#A40006")

  wordcloud2(df, size =0.5, minSize = 2,gridSize = 1,
             fontFamily = "微软雅黑", color = "random-dark")


}

# Sun Apr 27 16:39:54 2025 overlap of biobank------------------------------
{
  data <- readxl::read_excel("d:/temp/data.xlsx")

  data.ubk <- data %>%
    dplyr::mutate(
      Lipid_protein =
        case_when(grepl(pattern = "HDL",x = Metabolite)~"HDL",
                  grepl(pattern = "LDL",x = Metabolite)~"LDL",
                  grepl(pattern = "IDL",x = Metabolite)~"IDL")
    )

  library(ComplexHeatmap)
  library(circlize)

  ###col
  {
    col1 <- make_vector(
      x= rand_color(length(unique(data.ubk$Group))),
      name = unique(data.ubk$Group)
    )
    col2 <- make_vector(
      x= rand_color(length(unique(data.ubk$Lipid_protein))),
      name = unique(data.ubk$Lipid_protein)
    )

    col3 <- make_vector(
      x= c("#E64340","grey"),
      name = unique(data.ubk$LCMS)
    )

  }

  {
    circos.clear()
    #circos.initialize(factors = rownames(hm), xlim = c(1, 88))
    circos.heatmap(data.ubk$Group,
                   col = col1,
                   #split = data.l.selected$lclass,
                   cluster = F, track.height = 0.1,
                   rownames.side ="outside",
                   track.margin = c(0.01, 0.01)
    )

    circos.heatmap(data.ubk$Lipid_protein,
                   col = col2,
                   #split = data.l.selected$lclass,
                   cluster = F, track.height = 0.1,
                   rownames.side ="outside",
                   track.margin = c(0.01, 0.01)
    )

    circos.heatmap(data.ubk$LCMS,
                   col = col3,
                   #split = data.l.selected$lclass,
                   cluster = F, track.height = 0.1,
                   rownames.side ="outside",
                   track.margin = c(0.01, 0.01)
    )



    l1 <- Legend(labels  = names(col1),title = "Metabolite Class",
                 legend_gp = gpar(fill = col1))
    l2 <- Legend(labels  = names(col2),title = "Lipid Protein",
                 legend_gp = gpar(fill = col2))
    l3 <- Legend(labels  = names(col3),title = "LCMS detected",
                 legend_gp = gpar(fill = col3))
    lp <- ComplexHeatmap::packLegend(l1,l2,l3)
    draw(lp,x = unit(25, "cm"), y = unit(7, "cm"))
  }

  export::graph2pdf(file = "d:/temp/f.pdf",width = 10,height = 6)


}
# Tue Apr 29 14:49:09 2025 list all mz data------------------------------
{

  files.msraw <- dir("e:/",pattern = "raw$|wiff$",full.names = T,recursive = T)
  files.msraw.info <- file.info(files.msraw)
  files.zip <- dir("e:/",pattern = "zip$",full.names = T,recursive = T)
  fzl <- list()
  for (i.zip in files.zip) {
    message_with_time(i.zip)
    file.in.zip <- unzip(i.zip,list = T)%>%
      dplyr::mutate(zip = i.zip)
    fzl[[i.zip]] <- file.in.zip

  }

  db.files <- list(files.msraw = files.msraw.info,files.in.zip = fzl)
  #save(db.files,file = "temp/20250429.ms.files.rda")


  {
    load("temp/20250429.ms.files.rda")

    ### raw data
    {

      msraw <-  db.files$files.msraw%>%
        rownames_to_column("file.path")%>%
        dplyr::mutate(file.name = basename(file.path),
                      file.dir = dirname(file.path))%>%
        dplyr::filter(!grepl(pattern = "^Cal",x = file.name),
                      !grepl(pattern = "^E_Raw",x = file.name),
                      !grepl(pattern = "^Batch",x = file.name),
                      !grepl(pattern = "BLK",x = file.name),
                      !grepl(pattern = "blank",ignore.case = T,x = file.name),
                      !grepl(pattern = "condition",x = file.name),
                      !grepl(pattern = "BLK",x = file.name),
                      !grepl(pattern = "BLK",x = file.name))%>%
        dplyr::mutate(project =
                        gsub(pattern = "e:/|YHY.Cloud/|raw data|Raw.data.from.6600|Data|Analyst |Projects|Example|Results",
                             replacement = "",x = file.dir),
                      project =
                        gsub(pattern = "//*",
                             replacement = "_",x = project)
                      )%>%
        dplyr::group_by(project)%>%
        dplyr::filter(n()>3)%>%
        dplyr::mutate(size = size,
                      date = mtime)%>%
        dplyr::select(project,size,date)

      table(msraw$project)

    }

    ### zip file
    {
      library(tools)
      library(stringi)

      zip.files <- db.files$files.in.zip
      zip.files <- zip.files[grepl(pattern = "e:/Raw.data.from",x=names(zip.files))]

      zip.ms.files <- zip.files%>%
        rbindlist()%>%
        dplyr::filter(Length>0)%>%
        dplyr::mutate(Name = stri_replace_all_regex(str = Name,"[^[:print:]]", "_"),
                      ext = file_ext(Name))%>%
        dplyr::filter(ext %in% c("wiff","lcd","raw"))%>%
        dplyr::mutate(file.name = basename(Name),
                      file.dir = paste0(zip,"/",dirname(Name)))%>%
        dplyr::filter(!grepl(pattern = "^Cal",x = file.name),
                      !grepl(pattern = "^E_Raw",x = file.name),
                      !grepl(pattern = "^Batch",x = file.name),
                      !grepl(pattern = "BLK",x = file.name),
                      !grepl(pattern = "blank",ignore.case = T,x = file.name),
                      !grepl(pattern = "condition",x = file.name),
                      !grepl(pattern = "BLK",x = file.name),
                      !grepl(pattern = "BLK",x = file.name))%>%
        dplyr::mutate(project =
                        gsub(pattern = "e:/|YHY.Cloud/|raw data|Raw.data.from.|Data|Analyst |Projects|Example|Results|.zip|Positive|Negative|^_|^$|8050_2021_2021|6600_2022_01_10",
                             replacement = "",x = file.dir,ignore.case = T),
                      project =
                        gsub(pattern = "//*",
                             replacement = "_",x = project),
                      project =
                        gsub(pattern = "8050_2021_2021|6600_2022_01_10|8050_2021.12.10.gout.volidation_|6600_2023_09_18-|6600_2021_12_28_Liangningning_|6600_2021_10_09-Liangningning_|QEPlus_2022.10.04.ESCC.Tissue.Metabolomics_",
                             replacement = "",x = project,ignore.case = T),
                      project = str_short(project,100)
        )%>%
        dplyr::group_by(project)%>%
        dplyr::filter(n()>3)%>%
        dplyr::ungroup()%>%
        dplyr::mutate(size = Length,
                      date = Date)%>%
        dplyr::select(Name,project,size,date)



      table(zip.ms.files$project)


    }


    ### total
    {
      ms.files <- bind_rows(msraw,zip.ms.files)%>%
        dplyr::distinct(project,.keep_all = T)%>%
        dplyr::mutate(year = format(date,"%Y"))%>%
        dplyr::filter(year>2020)


      ggplot(ms.files,aes(x=1,group = year,fill = year))+
        geom_bar(col = "black",position = "dodge",width = 1)+
        geom_text(stat = "count",
                                 aes(x=1,label = paste0(..count..,"\n(", (2021:2024)[group],")")),
                                 position = position_dodge(width = 1),#force_pull  = 2,
                                 color = "black", size = 5)+
        #ggsci::scale_fill_npg()+
        scale_fill_manual(values = colramp(colors = c("white","#FE9D71","#A40006"))(seq(0.25,1,0.25)))+
        #xlim(c(0.55,2))+
        ylim(c(-2,80))+
        #coord_polar(theta = "y",direction = -1)+
        theme_void()+
        theme(legend.position = "none")->p
      p
      open_plot_win(p,4,2.5)



      ggplot(ms.files,aes(x=1,group = year,fill = year))+
        geom_bar(col = "white")+
        ggrepel::geom_text_repel(stat = "count",
                  aes(x=1.5,label = paste0(..count..,"\n(", (2020:2024)[group],")")),
                  position = position_stack(vjust = 0.5),#force_pull  = 2,
                  color = "black", size = 5)+
        #ggsci::scale_fill_npg()+
        scale_fill_manual(values = colramp(colors = c("white","#FE9D71","#A40006"))(seq(0,1,0.25)))+
        xlim(c(0.55,2))+
        #coord_polar(theta = "y",direction = -1)+
        theme_void()+
        theme(legend.position = "none")->p
    p
    open_plot_win(p,3,3)




    }

    ### proj-data
    {
      ms.files <- bind_rows(msraw,zip.ms.files)%>%
        dplyr::mutate(year = format(date,"%Y"))%>%
        dplyr::filter(year>2019)%>%
        dplyr::count(project,name = "count")%>%
        dplyr::arrange(count)%>%
        dplyr::ungroup()%>%
        dplyr::slice_max(count,n=20)%>%
        dplyr::mutate(
          project = sub(x = project,pattern = "^_|__$|_$",replacement = ""),
                      project = factor(project,levels = rev(project)))

      p <- ggplot(ms.files,aes(x = count, y = project,fill = count))+
        geom_bar(
                 col = "black",
                 stat = "identity")+
        geom_text(aes(x = count + 30,label = count))+
        scale_fill_gradient(low = "white",high = "#A40006")+
        labs(x = "Data count", y = "")+
        theme_classic()+
        theme(legend.position = "none",
              axis.text.y  = element_text(face = "bold",colour = "black",family = "sans") )
      p
      open_plot_win(p,6,6)


    }



  }
}
# Wed Apr 30 02:19:44 2025 MSdev stat------------------------------
{


  ### parse by files
  {
    # 1. 创建一个新的临时环境
    temp_env <- new.env()

    # 2. 使用 source() 将 .R 文件加载到临时环境中
    # 请替换为实际的文件路径
    file_path <- "R/dev_CFM.R"
    source(file_path, local = temp_env)

    # 3. 获取临时环境中的所有对象
    all_objects <- ls(envir = temp_env)

    # 4. 过滤出函数
    functions <- all_objects[sapply(all_objects, function(x) is.function(get(x, envir = temp_env)))]

    # 5. 统计函数的数量
    num_functions <- length(functions)

    # 6. 打印函数数量
    print(num_functions)
  }

  ### parse by name
  {
    library(MSCC)
    library(MSdb)
    msdev.obj <- c(ls(envir = environment(MSdev)),
                   ls(envir = environment(chemform_adduct)),
                   ls(envir = environment(get_KEGG_pathway)))

    msdev.obj.info <- data.frame(
      name = msdev.obj
    )%>%
      dplyr::mutate(
        class = case_when(
          grepl("xcms|Feature|feature|Chromatogram|chrom|Peaks",name)~"dev-xcms",
          grepl("MSdev|MS_Exp",name)~"MSdev",
          grepl("^analyze|DEP",name)~"Bioinfo",
          grepl("Spectra|spectra|sp",name)~"dev-Spectra",
          grepl("CFM|cfm",name)~"MSCC",
          grepl("MSCC|atom|bond|smile|sdf|chemform|formula",name)~"MSCC",
          grepl("dplyr|str",name)~"Base",
          grepl("MSIP|Metabolic_flux_network|MFN|iso|frag|Molecule|shiny|flux|mfn",name)~"MSIP",
          grepl("graph|node|edge|vis",name)~"MSgraph",
          grepl("KEGG|HMDB|MSDB|MSdb|msdb|Compound",name)~"MSdb",
          grepl("ggplot|plot|fella",name)~"Bioinfo",
          T~"Unclassified",
          T~"aaaa"
        )
      )%>%
      dplyr::group_by(class)%>%
      dplyr::mutate(count = n())%>%
      dplyr::ungroup()%>%
      dplyr::arrange(count)%>%
      dplyr::mutate(class = factor(class,levels = unique(class)))

    table(msdev.obj.info$class)

    ggplot(msdev.obj.info)+
      geom_bar(aes(x = 1,fill = class),position = "stack",col = "white")+
      ggsci::scale_fill_npg()+
      coord_polar(theta = "y",direction = -1)+
      labs(fill = "Function")+
      xlim(c(-0.5,2))+
      theme_void()+
      theme(legend.position = "bottom")->p

    open_plot_win(p,5.5,4)
  }


  ### code line stat
  {
    r.files <- c(dir("R/",pattern = "R$",full.names = T),
      dir("../MSCC/R/",pattern = "R$",full.names = T),
      dir("../MSconvertR/R/",pattern = "R$",full.names = T),
      dir("../MSdb/R/",pattern = "R$",full.names = T))

    r.stat.list <- list()
    for (i in r.files) {
      r.txt <- readLines(i)
      data.frame(
        r.file = basename(i),
        nchar = sum(nchar(r.txt)  ),
        nline = length(r.txt)
      )->r.stat.list[[i]]
    }
    r.stat <- rbindlist(r.stat.list)


  }

  ### plot
  {

    d <- tibble(value   = c(1,     2,     3,     5,     6,     7,     8,     9),
                s1 = c(TRUE,  FALSE, TRUE,  TRUE,  FALSE, TRUE,  FALSE, TRUE),
                s2= c(TRUE,  FALSE, FALSE, TRUE,  FALSE, FALSE, FALSE, TRUE),
                s3 = c(TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE),
                s4 = c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE))

    library(ggvenn)
    ggplot(d)+
      geom_venn(aes(A = s1,B=s2,C = s3),
                text_size = 0,
                set_name_size = 0,
                fill_color  = ggsci::pal_aaas(alpha = 0.1)(3),
                fill_alpha  = 0.2,
                stroke_color  = "white")+
      annotate("text",x = 0,y=-1,label = "26730\nCode Lines")+
      annotate("text",x = -1,y=1,label = "713\nFunction")+
      annotate("text",x = 1,y=1,label = "13\nClass")+
      theme_void()->p
    open_plot_win(p,1.5,1.5)







  }




  ### git commit
  {

    # 按天统计提交次数
    daily_commits <- df %>%
      mutate(date = as.Date(commit_time)) %>%
      group_by(date) %>%
      summarise(commit_count = n())

    # 绘制按天统计的提交密度图
    ggplot(daily_commits, aes(x = date, y = commit_count)) +
      geom_bar(stat = "identity",color = "#B02B1C") +
      labs(title = "Commit Frequency by Day",
           x = "Date", y = "Number of Commits") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # 按月分组并统计每月提交次数
    monthly_commits <- df %>%
      mutate(month = floor_date(commit_time, "month")) %>%
      group_by(month) %>%
      summarise(commit_count = n())

    # 检查按月统计的结果
    head(monthly_commits)

    ggplot(monthly_commits, aes(x = month, y = commit_count)) +
      geom_bar(stat = "identity",fill = "#B02B1C") +
      #labs(title = "Update Frequency of MSdev",
      #     x = "Date", y = "Number of Commits") +
      labs(x = NULL,,y=NULL)+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            )->p

    ggplot(monthly_commits, aes(x = month, y = commit_count)) +
      geom_bar(stat = "identity",fill = "#B02B1C") +
      labs(title = "Update Frequency of MSdev",
           x = "Date", y = "Number of Commits") +
      labs(x = NULL,,y=NULL)+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
      )->p

    open_plot_win(p,4,3)

  }



}



# Wed May 21 16:16:39 2025 CFMD fix------------------------------
{

  a <- get_CFM_data_from_smiles()
  6999 * (1-0.02-0.02-0.05)
  6598 * (1-0.02)

}

# Sat Jun  7 15:24:27 2025 MSIP ------------------------------
{

  msip.test <- MSdev("d:/MSIP_Test/FS/",
                     experimentInfo = MS_Experiment[10] )
  msip.test <- MSdev_msConvert(msip.test)
  msip.test <- MSdev_checkSampleInfo(msip.test)
  msip.test <- MSdev_xcmsProcessing(msip.test)


  MSdev_save(msip.test)
}
# Sun Jun  8 23:54:41 2025 plot_MSIPCore_spectra_consistency------------------------------
{
  cfmd <- msdev.13C1@statData$MSIP$isotopologues_data$FT01592_Positive$CFM_annotation
  MSIPCoreData <- msdev.13C1@statData$MSIP$isotopologues_data$FT01592_Positive$MSIP_result$M3$FS_1_13C
  sp.data <- MSIPCoreData@Spectra_data

  p <- plot_MSIPCore_spectra_consistency(MSIPCoreData)
  open_plot_win(p,12,8)
  p <- plot_MSIPCore_spectra_consistency_hm(MSIPCoreData)
  open_plot_win(p,17,24)
  ratio.annotated.peaks <- sp.data %>%
    dplyr::filter(!merged)%>%
    dplyr::mutate(
      annotated = !is.na(fragment_group)
    )%>%
    dplyr::group_by(annotated)%>%
    dplyr::summarise(int_sum = sum(intensity))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(ratio = int_sum/sum(int_sum))%>%
    dplyr::filter(annotated )%>%
    dplyr::pull(ratio)


  process.info <- MSIP_solve_computation_evaluate(
    msdev.13C1,
    show_message = F)
  process.info <- process.info%>%
    dplyr::filter(iso_count > 0)
  for (i in 1:nrow(process.info)) {

    msip.core <- msdev.13C1@statData$MSIP$isotopologues_data[[
      process.info$feature_id[i]
    ]]$MSIP_result[[str_isotope2_num(process.info$iso_count[i])]][[process.info$samples[i]]]
    if (nrow(msip.core@FG_map@FG.ratio.matrix)<3) {
      next
    }
    p <- plot_MSIPCore_spectra_consistency_hm(msip.core,
                                              title = paste0(process.info$name[i],"\n",
                                                             process.info$samples[i],"\n",
                                                             str_isotope2_num(process.info$iso_count[i])
                                                             ))

    export_graph2pdf(p,
                     file_path = "d:/temp/sp.cons.pdf",
                     width = 15,height = 1+0.5*length(unique(msip.core@Spectra_data$sp.id)),
                     append = T)
  }


}

# Fri Jun 13 14:54:34 2025 MFNA Simulate PDHi------------------------------
{


  mfn <- load_MFN()

  ### WT
  {
    glucose.smiles <- "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
    glucose.mig <- get_Molecule_igraph_from_smiles(glucose.smiles)

    Glu_full.mig <- Molecule_igraph_add_isotopomer(Molecule_igraph = glucose.mig,
                                                  isotopomer = "Tracer",
                                                  iso_vec = make_vector("[13]C",atom(glucose.mig,"C")) ,
                                                  abundance = 1)
    Glu_full.mig <- Molecule_igraph_remove_isotopomer(Glu_full.mig,"base")


    #
    mfn.wt <- Metabolic_flux_network_set_tracer(mfn,
                                             "C00267",Glu_full.mig)
    mfn.wt <- Metabolic_flux_remove_tracing(mfn.wt)
    mfn.wt <- Metabolic_flux_tracing(mfn.wt)
    MFN_manul_Shiny(mfn.wt,port = 9991)


  }

  ### PDHi
  {
    glucose.smiles <- "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
    glucose.mig <- get_Molecule_igraph_from_smiles(glucose.smiles)

    Glu_full.mig <- Molecule_igraph_add_isotopomer(Molecule_igraph = glucose.mig,
                                                   isotopomer = "Tracer",
                                                   iso_vec = make_vector("[13]C",atom(glucose.mig,"C")) ,
                                                   abundance = 1)
    Glu_full.mig <- Molecule_igraph_remove_isotopomer(Glu_full.mig,"base")


    #
    mfn.pdhi <- Metabolic_flux_network_set_tracer(mfn,
                                                "C00267",Glu_full.mig)
    mfn.pdhi <- Metabolic_flux_network_filter_reactions(mfn.pdhi,
                                                        rid = setdiff(vdata(mfn.pdhi)$id,"R00014"))

    mfn.pdhi <- Metabolic_flux_remove_tracing(mfn.pdhi)
    mfn.pdhi <- Metabolic_flux_tracing(mfn.pdhi)
    MFN_manul_Shiny(mfn.pdhi,port=9992)

  }

}




# Sun Jun 15 00:31:29 2025 MSIP Merge spectra------------------------------
{

  msdev.13C1 <- MSdev_load(
    "C:/Users/91879/OneDrive/Code/R/data/MSIP_data/240701_FS_ONE_POSITION/MSdev_2024_07_04.Rdata"
  )
  cfmd <- msdev.13C1@statData$MSIP$isotopologues_data$FT01592_Positive$CFM_annotation
  MSIPCoreData <- msdev.13C1@statData$MSIP$isotopologues_data$FT01592_Positive$MSIP_result$M3$FS_1_13C

  sp.iso <- msdev.13C1@statData$MSIP$isotopologues_data$FT01592_Positive$Spectra$M3$FS_1_13C

  MSIPCoreData <- get_MSIPCoreData(sp.iso = sp.iso,cfmd = cfmd,iso_count_max  = 3,ppm = 10)


  plot_MSIPCore_spectra_consistency_hm(MSIPCoreData)



}

# Tue Jun 17 22:36:26 2025 BiocParallel test ------------------------------
{

  bplapply(1:4,FUN = function(x){

    MSdev::atom(a)

  },BPPARAM = SnowParam(workers = 1,progressbar = T))


  atom(cfmd)

}



# Wed Jun 18 22:56:52 2025 DDA mine test------------------------------
{

  ### run after QC FS, once
  data.dir <- "d:/DDA.mine.test/pos"
  msdev.qe <- MSdev(rawDataDir = data.dir)
  msdev.qe <- MSdev_load("d:/DDA.mine.test/MSdev_2025_06_18.Rdata")
  msdev.qe <- MSdev_msConvert(msdev.qe)
  msdev.qe <- MSdev_xcmsProcessing(msdev.qe)

  msdev.qe@statData <- list()
  msdev.qe <- MSdev_get_Inclusion_Queue(msdev.qe)

  ### run after every time DDA acquired
  msdev.qe <- MSdev_get_Inclusion_List(msdev.qe)
  msdev.qe <- MSdev_add_sample(msdev.qe,raw.data.dir = "d:/DDA.mine.test/pos")
  msdev.qe <- MSdev_get_MS2acquisitionStat(msdev.qe)

  table(msdev.qe@statData$DDA_mine_queue_Positive$acquired)

}
# Sat Jun 21 17:05:31 2025 ------------------------------
{


  msdev.demo <- MSdev_extract_Spectra(msdev.demo)
  msdev.demo <- MSdev_match_Spectra_to_feature(msdev.demo)
  msdev.demo <- MSdev_annotation(
    msdev.demo,
    expand_adduct= T,
    cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite")

  msdev.demo <- MSdev_get_Stat(msdev.demo)
  MSdev_export(msdev.demo)
  MSdev_save(msdev.demo)



}

# Sun Jun 22 14:04:02 2025 Compare ------------------------------
{

  msdev.dda <- MSdev_load("d:/20250619_LR/MSdev_2025_06_19.Rdata")
  msdev.dda <- MSdev_checkSampleInfo(msdev.dda)
  msdev.dda <- MSdev_update_xcms_pdata(msdev.dda)
  msdev.dda <- MSdev_get_Stat(msdev.dda,candi = F)


  msdev.demo <- load_demo()
  msdev.demo <- MSdev_checkSampleInfo(msdev.demo)
  msdev.demo <- MSdev_update_xcms_pdata(msdev.demo)
  msdev.demo <- MSdev_get_Stat(msdev.demo,candi = F)




  ### tic
  {



    p1 <- plot_xcms_TIC(msdev.dda@xcmsData$PositiveMS1,title = "Positive TOF6600")+theme(legend.position = "none")
    p2 <- plot_xcms_TIC(msdev.demo@xcmsData$PositiveMS1,title = "Positive QE")+theme(legend.position = "none")

    p <- p1/p2
    open_plot_win(p)

    p1 <- plot_xcms_TIC(msdev.dda@xcmsData$NegativeMS1,title = "Negative TOF6600")+theme(legend.position = "none")
    p2 <- plot_xcms_TIC(msdev.demo@xcmsData$NegativeMS1,title = "Negative QE")+theme(legend.position = "none")
    p <- p1/p2
    open_plot_win(p)
  }


  ### FEATURE
  {

    cp.qe <- msdev.dda@statData$feature.se %>%rowData()
    cp.qe <- cp.qe[which(cp.qe$qc_rsd < 0.3),]
    cp.6600 <- msdev.demo@statData$feature.se %>%rowData()
    cp.6600 <- cp.6600[which(cp.6600$qc_rsd < 0.3&cp.6600$peakMaxo > 1e3),]

    match.df <-  match_mz_rt(mz1 = cp.qe$mzmed,
                             rt1 = cp.qe$rtmed,
                             mz2 = cp.6600$mzmed,
                             rt2 = cp.6600$rtmed ,mz.ppm = 20,rt.tol = Inf )

    overlap.id <- paste0("Overlap",num2str(1:nrow(match.df)))
    qe.uid <- paste0("QE",setdiff(1:nrow(cp.qe),match.df$ion1))
    tof.uid <- paste0("TOF",setdiff(1:nrow(cp.6600),match.df$ion2))

    library(ggvenn)

    p <- ggvenn(data = list("QE" =  c(qe.uid,overlap.id),
                            "6600" =   c(tof.uid,overlap.id)),
                fill_color = c("#E64B35",
                               "#4DBBD5"),
                show_outside =  "none",
                stroke_color = "white")

    open_plot_win(p,6,6)


    cp.overlap <- match.df %>%
      dplyr::mutate(log_int_QE = log10(cp.qe$peakMaxo)[ion1],
                    log_int_6600 = log10(cp.6600$peakMaxo)[ion1])


    p <- plot_cor_density(cp.overlap$log_int_QE,
                     cp.overlap$log_int_6600,
                     xlab = "Log10 intensity QE",
                     ylab = "Log10 intensity 6600")


    open_plot_win(p,5,5)


    ### QE
    {
      cp.qe <- cp.qe%>%
        dplyr::mutate(no = 1:n(),
                      overlap = no%in% match.df$ion1,
                      log_int = log10(peakMaxo))

      p1 <- ggplot(cp.qe)+
        geom_histogram(aes(x = log_int,colour = overlap),position = "dodge")

      p2 <- ggplot(cp.qe)+
        geom_histogram(aes(x = qc_rsd,colour = overlap),bins = 30,position = "dodge")

      open_plot_win(p1/p2,8,8)



    }



  }


  ### MS2
  {

    cp.qe <- msdev.dda@statData$feature.se %>%rowData()
    cp.qe <- cp.qe[which(cp.qe$qc_rsd < 0.3),]
    cp.6600 <- msdev.demo@statData$feature.se %>%rowData()
    cp.6600 <- cp.6600[which(cp.6600$qc_rsd < 0.3&cp.6600$peakMaxo > 1e3),]

    sum(lengths(cp.6600$score.ms2)==0)

    data.frame(
      MS2 = c(T,F),
      count = c( sum(lengths(cp.6600$ms2_id)!=0),
                 sum(lengths(cp.6600$ms2_id)==0))
    )%>%
      ggplot()+
      geom_bar(aes(x = 1 , y = count,fill = MS2),stat = "identity",position = "stack")+
      geom_text(aes(x = 1 , y = count * 0.8,label = count))+
      coord_polar(theta = "y")+
      scale_fill_npg()+
      #labs(title = "6600")+
      theme_void()->p1


    data.frame(
      MS2 = c(T,F),
      count = c( sum(lengths(cp.qe$ms2_id)!=0),
                 sum(lengths(cp.qe$ms2_id)==0))
    )%>%
      ggplot()+
      geom_bar(aes(x = 1 , y = count,fill = MS2),stat = "identity",position = "stack")+
      geom_text(aes(x = 1 , y = count * 0.8,label = count))+
      coord_polar(theta = "y")+
      scale_fill_npg()+
      #labs(subtitle = "QE")+
      theme_void()->p2

    open_plot_win(p2/p1,5,10)


    p1 <- plot_xcms_ms2_distribution(msdev.dda@xcmsData$NegativeMS1)
    p2 <- plot_xcms_ms2_distribution(msdev.demo@xcmsData$NegativeMS1)
    open_plot_win(p1/p2,10,10)
  }


  ### MS2 score
  {

    hist(x = cp.6600$score[which(cp.6600$score>0)])
    hist(x = cp.qe$score[which(cp.qe$score>0)])




  }

  cp.qe <- msdev.dda@statData$feature.se %>%rowData()%>%as.data.frame()%>%
    dplyr::filter(qc_rsd < 0.3,score > 0.6 )
  cp.6600 <- msdev.demo@statData$feature.se %>%rowData()%>%as.data.frame()%>%
    dplyr::filter(peakMaxo > 1e3,qc_rsd < 0.3,score > 0.6 )

  length(intersect(cp.qe$compound_id,cp.6600$compound_id))

  match.df <-  match_mz_rt(mz1 = cp.qe$mzmed,
                           rt1 = cp.qe$rtmed,
                           mz2 = cp.6600$mzmed,
                           rt2 = cp.6600$rtmed,
                           mz.ppm = 25,rt.tol = Inf )

  overlap.id <- paste0("Overlap",num2str(1:nrow(match.df)))
  qe.uid <- paste0("QE",setdiff(1:nrow(cp.qe),match.df$ion1))
  tof.uid <- paste0("TOF",setdiff(1:nrow(cp.6600),match.df$ion2))

  library(ggvenn)

  p <- ggvenn(data = list("QE" =  c(qe.uid,overlap.id),
                          "6600" =   c(tof.uid,overlap.id)),
              fill_color = c("#E64B35",
                             "#4DBBD5"),
              show_outside =  "none",
              stroke_color = "white")

  p
  open_plot_win(p,6,6)

  cp.overlap <- cp.6600[match.df$ion2,]
  cp.not <- cp.6600[-match.df$ion2,]


  cp.overlap <- cp.qe[match.df$ion1,]%>%
    dplyr::mutate(
      kegg_id = case_when(
        is.na(kegg_id)~cp.6600$kegg_id[match.df$ion2],
        T~kegg_id)
    )

  cp.not <- cp.qe[-match.df$ion1,]

  #cp.path <- MSdb::get_KEGG_compound_pathway_df()
  cp.pathm <- cp.path%>%
    dplyr::filter(grepl("Metabolism",CLASS))

  sum(cp.overlap$kegg_id %in% cp.path$COMPOUND.ID)/nrow(cp.overlap)
  sum(cp.not$kegg_id%in% cp.pathm$COMPOUND.ID)/nrow(cp.not)


  {


    data.frame(
      include.in.KEGG = c(T,F),
      count = c( sum(cp.overlap$kegg_id %in% cp.path$COMPOUND.ID),
                 sum(!cp.overlap$kegg_id %in% cp.path$COMPOUND.ID))
    )%>%
      ggplot()+
      geom_bar(aes(x = 1 , y = count,fill = include.in.KEGG),stat = "identity",position = "stack")+
      geom_text(aes(x = 1 , y = count * 0.8,label = count))+
      coord_polar(theta = "y")+
      scale_fill_npg()+
      #labs(title = "6600")+
      theme_void()->p1


    data.frame(
      include.in.KEGG = c(T,F),
      count = c( sum(cp.not$kegg_id%in% cp.pathm$COMPOUND.ID),
                 sum(!cp.not$kegg_id%in% cp.pathm$COMPOUND.ID))
    )%>%
      ggplot()+
      geom_bar(aes(x = 1 , y = count,fill = include.in.KEGG),stat = "identity",position = "stack")+
      geom_text(aes(x = 1 , y = count * 0.8,label = count))+
      coord_polar(theta = "y")+
      scale_fill_npg()+
      #labs(subtitle = "QE")+
      theme_void()->p2

    open_plot_win(p2/p1,5,7)
  }


  score.df <-
    data.frame(
      score =c(cp.overlap$score,cp.not$score),
      overlap = c(rep(T,length(cp.overlap$score)),rep(F,length(cp.not$score)))
    )
  score.df$score <- score.df$score+rnorm(710,0,0.1)

  ggplot(score.df)+
    geom_histogram(aes( x = score,col = overlap),
                   position = "dodge")




}
