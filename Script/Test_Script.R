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


# Wed Jul 24 20:08:40 2024 vis------------------------------
cfmd <- msdev.purity@statData$MSIP$isotopologues_data$FT00210_Negative$CFM_annotation

get_cfm_data_sdf_igraph(cfmd)%>%
  vis_sdf_igraph() %>%
  htmlwidgets::onRender("
    function(el, x) {
      var container = document.getElementById(el.id);
      container.style.border = '2px solid black';  // Add border
      container.style.padding = '10px';  // Optional: Add padding
      container.style.margin = '0 auto';
      container.style.display = 'block';
    }
  ")%>%
  open_visNet()


# Thu Jul 25 13:49:01 2024 ------------------------------
{
  library(devtools)
  load_all()
  msdev.purity <- load_as_var("C:/Users/91879/OneDrive/Code/R/data/MSIP_data/240601_ThreeGroup/MSdev_2024_07_19.Rdata")
  process.info <- MSIP_solve_computation_evaluate(msdev.purity)

}

# Fri Jul 26 14:57:41 2024 ------------------------------
{
  msdev.purity <- load_as_var("C:/Users/91879/OneDrive/Code/R/data/MSIP_data/240601_ThreeGroup/MSdev_2024_07_19.Rdata")
  iso.data <- msdev.purity@statData$MSIP$isotopologues_data$FT12300_Positive
  msip.core <- iso.data$MSIP_result$M4$U
  heatmap_MSIPFragmentMap(msip.core@FG_map)
  sp.data <- shiny_get_sp_data(iso.data,sample = "U","M4")
  shiny_plotly_iso_data_spectra(sp.data)%>%
    open_visNet()

  heatmap_MSIPIsoformMap(msip.core@solve$MSIPIsoformMap)

}

# Sat Jul 27 16:03:41 2024 MSIP solve------------------------------
{
  library(nloptr)
  library(pracma)

  find_all_solutions <- function(Q, c, A, b, lb, ub, x0 = NULL) {
    # Ensure x0 is provided, if not set to a reasonable default
    if (is.null(x0)) {
      x0 <- rep(0.5, ncol(A))
    }

    # Define the objective function
    objective_function <- function(x) {
      0.5 * sum(x * (Q %*% x)) + sum(c * x)
    }

    # Define the gradient of the objective function
    gradient_function <- function(x) {
      Q %*% x + c
    }

    # Define the equality constraint function
    constraint_function <- function(x) {
      A %*% x - b
    }

    # Define the constraint gradient function
    constraint_gradient <- function(x) {
      A
    }

    # Solve the quadratic programming problem
    result <- nloptr(
      x0 = x0,
      eval_f = objective_function,
      eval_grad_f = gradient_function,
      lb = lb,
      ub = ub,
      eval_g_eq = constraint_function,
      eval_jac_g_eq = constraint_gradient,
      opts = list(
        algorithm = "NLOPT_LD_AUGLAG_EQ",
        xtol_rel = 1.0e-8,
        print_level = 2,
        local_opts = list(
          algorithm = "NLOPT_LD_SLSQP",
          xtol_rel = 1.0e-8
        )
      )
    )

    # Find the particular solution
    x_p <- result$solution

    # Compute the null space of A
    null_space_A <- null(A)

    # Define a function to get general solutions
    get_general_solution <- function(z) {
      x_p + null_space_A %*% z
    }

    return(list(
      particular_solution = x_p,
      null_space_basis = null_space_A,
      get_general_solution = get_general_solution
    ))
  }

  # Define matrices and vectors for the problem
  Q <- matrix(c(2, 0, 0,
                0, 2, 0,
                0, 0, 2), nrow = 3, byrow = TRUE)
  c <- c(-2, -5, -3)
  A <- matrix(c(1, 1, 0,
                0, 1, 1), nrow = 2, byrow = TRUE)
  b <- c(1, 2)
  lb <- c(0, 0, 0)
  ub <- c(1, 1, 1)

  # Find all solutions
  solutions <- find_all_solutions(Q, c, A, b, lb, ub)

  # Example usage: Compute some specific general solutions
  # The null space dimension
  null_dim <- ncol(solutions$null_space_basis)

  # Define free variable vectors based on null space dimension
  z1 <- rep(1, null_dim)  # Example vector in the null space dimension
  z2 <- rep(0, null_dim)  # Another example vector

  x1 <- solutions$get_general_solution(z1)
  x2 <- solutions$get_general_solution(z2)

  print("Particular Solution:")
  print(solutions$particular_solution)

  print("Null Space Basis:")
  print(solutions$null_space_basis)

  print("General Solution with z1:")
  print(x1)

  print("General Solution with z2:")
  print(x2)

}
# Sat Jul 27 19:10:34 2024 ------------------------------
{
  object <- msdev.purity
  iso.data <- object@statData$MSIP$isotopologues_data
  iso_ele <- get_MSdev_iso_ele(object)
  target_ele <-get_ele_uniso(iso_ele)
  all.sample <- .get_MSIP_tracer(object)
  traced.sample <- names(na.omit(all.sample))


  comp.eval.list <- list()
  for (i in seq_along(iso.data)) {

    cfmd <- iso.data[[i]]$CFM_annotation
    cfmd.ig <- get_cfm_data_sdf_igraph(cfmd)
    this.atom <- get_sdf_igraph_atom(cfmd.ig,ele = target_ele)
    this.ele.count <-length(this.atom)
    iso_count <- names(iso.data[[i]]$Spectra)%>%
      str_isotope2_num()%>%
      setdiff(0)



    ms2_count.matrix <- iso.data[[i]]$compound_info$ms2_count

    for (j in 1:nrow(comp.eval)) {
      msip.core <- iso.data[[i]][["MSIP_result"]][[str_isotope2_num(comp.eval$iso_count[j])]][[
        comp.eval$samples[j]]]
      if (is.null(msip.core)) {
        next
      }

      isoform.map <- (msip.core@solve$MSIPIsoformMap@isoform.map)
      mes <- paste0("constraints: ",nrow(isoform.map),", variables: ",
                    ncol((isoform.map)))%>%
        paste0(collapse = "\n")
      message(mes)
      message("")
    }





    comp.eval.list[[i]] <-comp.eval

  }

  comp.eval <- do.call(rbind,comp.eval.list)

  return(invisible( comp.eval ))

}

{

  compound.df <- get_MSIP_compound_info(msdev.purity@statData$MSIP$isotopologues_data)
  iso.data <- msdev.purity@statData$MSIP$isotopologues_data[[10]]
  natural.ratio <- 0.8
  cfmd <- iso.data$CFM_annotation
  sp.iso <-iso.data$Spectra$M1$Con
  plotSpec(sp.iso)
  ppm = 10
  iso_count <- 1

  sp.raw.data <- get_Spectra_data(sp.iso)
  sp.data <- CFM_annotate_isotopologues(sp.iso,
                                        cfmd  = cfmd,
                                        ppm = ppm,
                                        iso_count = iso_count)


  msip.core <- get_MSIPCoreData(sp.iso = sp.iso,
                                cfmd = cfmd,
                                iso_count = iso_count,
                                ppm = ppm)
  msip.core <- MSIPCore_solve(msip.core,int_thresh = 10^3.68)
  #msip.core
  heatmap_MSIPFragmentMap(msip.core@FG_map)
  heatmap_MSIPIsotopomerMap(msip.core@solve$MSIPIsotopomerMap)
  im <- msip.core@solve$MSIPIsotopomerMap
  plot(lengths(im@solve$isoform.set)/length(im@isoform.defination),
       im@solve$isoform.set.prob,
       xlim = c(0,max( im@solve$isoform.set.prob)*1.2),
       ylim = c(0,max( im@solve$isoform.set.prob)*1.2))
  abline(a=0,b=1)



}

# Tue Jul 30 10:10:21 2024 SWY biomarkers------------------------------
msdev.gout.biomarkers <- MSdev("d:/2024.07.29.Gout.biomarkers/Result/")
msdev.gout.biomarkers <- MSdev_msConvert(msdev.gout.biomarkers)
msdev.gout.biomarkers <- MSdev_checkSampleInfo(msdev.gout.biomarkers)

sp <- Spectra(
  msdev.gout.biomarkers@sampleInfo$msData.files)
sp.ms2 <- filterMsLevel(sp,2)
sp.data <- spectraData(sp.ms2)%>%
  as.data.frame()

prm.file <- "c:/Users/91879/OneDrive/Code/R/Projecct/2023.09.11.Gout.SWY/Figure/biomarkers.to.PRM.xlsx"
prm.list <- rbind( readxl::read_excel(prm.file,sheet = 1),
                   readxl::read_excel(prm.file,sheet = 2))
umz <- unique(precursorMz(sp.ms2))
na.idx <- match_mz(umz, prm.list$mz)%>%is.na()
umz[na.idx]

sp.ms2 <- Spectra_set_MEM_backend(sp.ms2)

sp.ms2.list  <- split(sp.ms2,precursorMz(sp.ms2))


sp.ms2.list[[6]]%>%
  combineSpectra_groupby_ce()%>%
  plotSpec()


sp.ms2.list[[2]]%>%
  plot_Spectra_product_CE_curve()

{
this.sp <- sp.ms2.list[[3]]
df <- data.frame(rt = rtime(this.sp),
           int = (this.sp$totIonCurrent),
           ce = collisionEnergy(this.sp))
ggplot(df,aes(x = rt,y=int,col = factor(ce)))+
  geom_point()

rt.max <- df$rt[which.max(df$int)]
this.sp <- this.sp%>%
  filterRt(c(rt.max-5,rt.max+5))
plot_Spectra_product_CE_curve(this.sp)+
  scale_y_log10()

  }

# Tue Jul 30 15:35:18 2024 ------------------------------

idx <- match_mz(precursorMz(sp.ms2),prm.list$mz  )
sp.ms2$id <- idx
sp.ms2.list  <- split(sp.ms2,sp.ms2$id )
trans.list <-list()
for (i in 1:nrow(prm.list)) {

  this.sp <- sp.ms2.list[[i]]
  df <- data.frame(rt = rtime(this.sp),
                   int = (this.sp$totIonCurrent),
                   ce = collisionEnergy(this.sp))
  p1 <-ggplot(df,aes(x = rt,y=int,col = factor(ce)))+
    geom_point()+
    scale_color_npg()+
    labs(x = "Retention time",y = "Intensity", col = "CE")

  rt.max <- df$rt[which.max(df$int)]
  this.sp <- this.sp%>%
    filterRt(c(rt.max-5,rt.max+5))
  p2 <- plot_Spectra_product_CE_curve(this.sp)+
    #scale_y_log10()+
    guides(col = guide_legend(ncol= 1))
  p2
  p3 <- this.sp%>%
    combineSpectra_groupby_ce(minProp = 0.2)%>%
    plot_Spectra_CE()
  p3
  p.merged <- p1/p2/p3+plot_annotation(title = paste0(
    "Compound: ",prm.list$Compound_name[i],"\n",
    "mz = ",format(prm.list$mz[i],digit = 4)," ; ",
    "rt = ",format(prm.list$rt[i],digit = 2)
  ))
  #p.merged
  #export_graph2pdf(p.merged,file_path = "a.pdf",
  #              width = 7,height = 10,append = T)
  message(i)
   trans.df <- get_Spectra_transition(this.sp)%>%
     dplyr::mutate(fid = prm.list$name[i],
                   compund = prm.list$Compound_name[i])

   trans.list[[i]] <-trans.df
}

trans.df.merged <-  do.call(rbind, trans.list)
write.xlsx(trans.df.merged,file.dir = "b.xlsx")


# Wed Jul 31 19:09:09 2024 ------------------------------
qe.list <- list()
for (i in 1:2) {

  fdf <- msdev.purity@statData$MSIP$isotopologues_table[[i]]
  fdf.m1 <- fdf%>%
    dplyr::filter(int_mean_nontracer>1e5,
                  ms1_purity > 0.8,
                  !is.na(iso_seed),
                  !is.na(compound_id),
                  iso_count %in% c(0,1))%>%
    dplyr::group_by(iso_seed)%>%
    dplyr::filter(all(c(0,1) %in%iso_count ))

  qe.list[[i]] <- QE_list_2feature_def(fdf.m1)


}


write.xlsx(qe.list , file.dir = "inclusion.list.xlsx")





shiny_plotly_natural_ratio(0.123456789)%>%
  open_visNet()


# Fri Aug  2 15:25:08 2024 ------------------------------
iso.data <- msdev.purity@statData$MSIP$isotopologues_data[[9]]

m1.idx <- 1
for (i in 1:63) {


  iso.data <- msdev.purity@statData$MSIP$isotopologues_data[[i]]
  if ("M1"%in%names(iso.data$MSIP_result)) {
    m1.idx <- c(m1.idx,i)
  }
}
m1.idx <- unique(m1.idx)

i=6
{

  iso.data <- msdev.purity@statData$MSIP$isotopologues_data[[m1.idx[i]]]
  message_with_time(iso.data$compound_info$name)
  natural.ratio <- 0.8
  cfmd <- iso.data$CFM_annotation
  sp.iso <-iso.data$Spectra$M1$Con
  plotSpec(sp.iso)
  ppm = 5
  iso_count <- 1

  sp.raw.data <- get_Spectra_data(sp.iso)
  sp.data <- CFM_annotate_isotopologues(sp.iso,
                                        cfmd  = cfmd,
                                        ppm = ppm,
                                        iso_count = iso_count)


  msip.core <- get_MSIPCoreData(sp.iso = sp.iso,
                                cfmd = cfmd,
                                iso_count = iso_count,
                                ppm = ppm)
  msip.core <- MSIPCore_solve(msip.core,
                              int_thresh = 10^3.6,
                              certainty_thresh = 0.8)
  #msip.core
  heatmap_MSIPFragmentMap(msip.core@FG_map)
  #heatmap_MSIPIsotopomerMap(msip.core@solve$MSIPIsotopomerMap)
  im <- msip.core@solve$MSIPIsotopomerMap

  df <- data.frame(
    natural.prob =lengths(im@solve$isotopomer.set)/length(im@isotopomer.defination),
    predict.prob = im@solve$isotopomer.set.prob
  )

  ggplot(df)+
    geom_point(aes(x = natural.prob,y = predict.prob),
               color ="#C43E1C",size = 5)+
    geom_abline(slope = 1,intercept = 0)+
    labs(title = iso.data$compound_info$name)+
    xlim(c(0,1))+
    ylim(c(0,1))+
    theme_bw()->p
  p



}
open_plot_win(p)


# Fri Aug  2 16:27:54 2024 ------------------------------
{

  load_all()
  df <- data.frame(intensity = 10^runif(200,1,8))%>%
    dplyr::mutate(weight = .intensity_weight(intensity))


  ggplot(df)+
    geom_point(aes(x = log10(intensity),y=weight))

}


# Sat Aug  3 19:45:47 2024 ------------------------------
cfmd <- CFM_annotate_by_predict("N[C@@H](CCC(O)=O)C(O)=O",param_adduct = "[M-H]-")
cfmd <- CFM_data_get_igraph(cfmd)%>%
  cfm_data_get_fragment_group()%>%
  CFM_data_get_atom_map()

{
  load_all()
  shiny_vis_cfmd(cfmd)

  }



a <- msip.core@FG_map
b <- a %>%
  MSIPFragmentMap_filter_intensity()%>%
  MSIPFragmentMap_filter_certainty()

# Tue Aug  6 16:07:14 2024 ------------------------------
i=9
{

  iso.data <- msdev.purity@statData$MSIP$isotopologues_data[[i]]
  message_with_time(iso.data$compound_info$name)
  natural.ratio <- 0.8
  cfmd <- iso.data$CFM_annotation
  sp.iso <-iso.data$Spectra$M1$Con
  plotSpec(sp.iso)
  ppm = 5
  iso_count <- 1

  sp.raw.data <- get_Spectra_data(sp.iso)
  sp.data <- CFM_annotate_isotopologues(sp.iso,
                                        cfmd  = cfmd,
                                        ppm = ppm,
                                        iso_count = iso_count)


  msip.core <- get_MSIPCoreData(sp.iso = sp.iso,
                                cfmd = cfmd,
                                iso_count = iso_count,
                                ppm = ppm)
  msip.core <- MSIPCore_solve(msip.core,
                              int_thresh = 10^6,
                              certainty_thresh = 0.5)

  #msip.core
  heatmap_MSIPFragmentMap(msip.core@FG_map)
  #heatmap_MSIPIsotopomerMap(msip.core@solve$MSIPIsotopomerMap)
  im <- msip.core@solve$MSIPIsotopomerMap

  df <- data.frame(
    natural.prob =lengths(im@solve$isotopomer.set)/length(im@isotopomer.defination),
    predict.prob = im@solve$isotopomer.set.prob
  )

  ggplot(df)+
    geom_point(aes(x = natural.prob,y = predict.prob),
               color ="#C43E1C",size = 5)+
    geom_abline(slope = 1,intercept = 0)+
    labs(title = iso.data$compound_info$name)+
    xlim(c(0,1))+
    ylim(c(0,1))+
    theme_bw()->p
  p



}



# Fri Aug  9 14:04:43 2024 shiny spectra update------------------------------


iso.data <- msdev.purity@statData$MSIP$isotopologues_data[[8]]

shiny_get_sp_data(iso.data,"U","M2")%>%
  shiny_plotly_iso_data_spectra(show.rawdata = F)%>%
  open_visNet()



iso_data <- iso.data
sample <- "U"
iso_count = "M1"



sp.m0.frag.data <- CFM_annotate_isotopologues(sp.m0,
                                              cfmd  = iso_data$CFM_annotation,
                                              ppm = 10,
                                              iso_count = 0)

sp.frag.data <- CFM_spectra_data_merge(sp.frag.data,iso_count)


shiny_plotly_iso_distribution()%>%
  open_visNet()


plot_ly(x=1,y=1) %>%
  layout(
    xaxis = list(showline = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(showline = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
  )%>%open_visNet()


# Sat Aug 10 19:17:52 2024 ------------------------------
ggplot(a.df)+
  geom_point(aes(x = rtime , y = precursorMz))


plot_ly(a.df)%>%
  add_markers(x = ~rtime , y = ~precursorMz )%>%
  open_visNet()





# Mon Aug 12 23:44:20 2024 ------------------------------
sp.df <- spectraData(sp.ms2)%>%
  as.data.frame()%>%
  dplyr::distinct(polarity,precursorMz)


object <- msdev.M1
ppm = 5


sp.ms2 <- onDiskData_retrieve(object@spectra$MS2_Spectra)

sp.pol <- filterPolarity(sp.ms2,1)
sample.info.pol <- object@sampleInfo%>%
  dplyr::filter(polarity == 1)%>%
  dplyr::mutate(sample = as.numeric(factor(msData.files)))

{
  sp.peaks.df <- data.frame(
    mz = precursorMz(sp.pol),
    rt = rtime(sp.pol),
    sample.files = dataOrigin(sp.pol),
    into = sp.pol$totIonCurrent
  )%>%
    dplyr::mutate(sample = match_path(sample.files,
                                      sample.info.pol$msData.files),
                  sample = sample.info.pol$sample[sample])

  sp.peaks.matrix <- sp.peaks.df%>%
    dplyr::mutate(mzmin = mz,
                  mzmax = mz,
                  rtmin = rt-30,
                  rtmax = rt+30)%>%
    dplyr::select(any_of(c("mz","mzmin","mzmax",
                           "rt","rtmin","rtmax",
                           "into","intb","maxo",
                           "sn","sample")))%>%
    as.matrix()
  sp.peaks.data <- sp.peaks.df%>%
    dplyr::mutate(ms_level = 1,
                  ms_level = as.integer(ms_level),
                  is_filled = F)%>%
    S4Vectors::DataFrame()

  ion_df <- do_groupChromPeaks_density(sp.peaks.df,
                                       bw = 30,
                                       sampleGroups = sample.info.pol$sample.source,
                                       binSize = 0.001,
                                       ppm = ppm)

  ion_table <-ion_df %>%
    dplyr::mutate(feature_id = paste0("FTS",num2str(1:n())),
                  .before = mzmed)

  }

MsFeatureData <- new("MsFeatureData",
                     chromPeaks = sp.peaks.matrix,
                     chromPeakData = sp.peaks.data,
                     featureDefinitions =  S4Vectors::DataFrame(ion_table))
XCMSnExp <- new("XCMSnExp")
XCMSnExp@msFeatureData <- MsFeatureData


b <- featureDefinitions(a)%>%
  as.data.frame()
edit_df_in_excel(b)



plot_ly(b)%>%
  add_markers(x = ~mzmed,
              y = ~rtmed,
              color = ~iso_seed)%>%
  open_visNet()



xcms.xcms <- xcms_from_ms2_spectra(sp.pol,
                           sample.info,ppm = 5,
                           peak_width = 30
                           )
a <- xcms_get_feature_isotopologues(xcms.xcms,
                                    ppm = 10,
                                    max_label = 5,
                                    rt.tol = 30)
nrow(featureDefinitions(a))


chemform_adduct("C5[13]C1H12O5",adduct = "[M+H]+")%>%format(digit=10)

chemform_adduct("C3[13]C2H10N2O3",adduct = "[M+H]+")%>%format(digit=10)


