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
Metabolic_flux_network <- load_MFN()
MFN_manul_Shiny(Metabolic_flux_network)



