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
from.smiles[1]%>%
  vis_Molecule_igraph_smiles()

from.smiles[1]%>%
  sub(pattern = "\\.","-",.)%>%
  canonicalize_smiles()%>%
  vis_Molecule_igraph_smiles()


