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


# Sat Aug 17 11:55:31 2024 ------------------------------
iso.list <- msdev.M1@statData$MSIP$isotopologues_data
iso.data <- iso.list[[68]]







# Mon Aug 19 12:15:25 2024 ------------------------------

fp <- "c:/Users/91879/OneDrive/Code/R/data/MSIP_data/240601_ThreeGroup/Negative_Chromatograms.rds"
chrom.neg <- readRDS(fp)
chrom.neg <- onDiskData(chrom.neg,fp)
object@xcmsData$Negative_Chromatograms <- chrom.neg


MSdev_save(object)



# Mon Aug 19 16:12:47 2024 PRM pos+neg------------------------------
sample.info <- msdev.STD@sampleInfo
prm.both.file <- sample.info%>%
  dplyr::filter(sample.info$ExpTime>"2024-08-15")%>%
  dplyr::pull(msData.files)

prm.both.xcms <- readMSData(prm.both.file,mode = "onDisk")
prm.both.fdata <- get_xcms_scan_Stat(prm.both.xcms)%>%
  dplyr::filter(msLevel==2)
hist(log10(prm.both.fdata$totIonCurrent))
#edit_df_in_excel(prm.both.fdata)

table(prm.both.fdata$precursorMZ)
prm.both.fdata <- prm.both.fdata%>%
  dplyr::mutate(polarity = c("1"="pos",
                             "0"="neg")[as.character(polarity)])%>%
  dplyr::select()


ggplot(prm.both.fdata) +
  geom_point(aes(x = retentionTime ,
                 y = precursorMZ,
                 col = polarity))


plotSpec(a[100])
# Mon Aug 19 18:44:52 2024 add param record to msip solve------------------------------
iso.list <- msdev.M1@statData$MSIP$isotopologues_data


msip.core <- iso.list$FT13758_Positive$MSIP_result$M1$Con

a <- MSIPCore_solve(msip.core)


sp <- msdev.M1@spectra$MS2_Spectra%>%
  onDiskData_retrieve()

size_of(sp)
sp.list <- split(sp,(sp$sp_id))

size_of(sp.list)

# Wed Aug 21 13:36:41 2024 ------------------------------
iso.list <- msdev.M1@statData$MSIP$isotopologues_data

a <- get_MSIP_compound_info(msdev.M1,vars = c("all"))%>%
  dplyr::arrange(compound_id)

p <- heatmap_MSIPFragmentMap(MSIPCoreData1@FG_map)
open_plot_win(p,8,5)


object <- msdev.M1
MSIP_merge()


msip.core <- MSIPCore_solve(
  MSIPCoreData,
  int_thresh = 10^4,
  certainty_thresh = 0.99
)

msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set%>%
  lengths()
msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob

heatmap_MSIPFragmentMap(msip.core@FG_map)

# Thu Aug 22 09:54:45 2024 ------------------------------
###  working on MSIPCore_merge
msip.core.merge.data <- list(MSIPCoreData1,MSIPCoreData2)
save(msip.core.merge.data,
     file = "msip.core.merge.rda")
# Fri Aug 23 13:13:45 2024 ------------------------------
load("msip.core.merge.rda")
MSIPCoreData1 <- msip.core.merge.data[[1]]
MSIPCoreData2 <- msip.core.merge.data[[2]]

MSIPCore_merge()

MSIPCoreData <- MSIPCore_merge(MSIPCoreData1,
                               MSIPCoreData2)
MSIPCoreData <- MSIPCore_solve(MSIPCoreData,
                               int_thresh = 10^4,
                               certainty_thresh = 0.8
                                 )


MSIPCoreData <- msdev.M1@statData$MSIP$isotopologues_data$FT13478_Positive$MSIP_result$M1$Con

heatmap_MSIPFragmentMap(MSIPCoreData@FG_map)


sp.df <- spectraData(sp.ms2)%>%
  as.data.frame()
sp.df <- sp.df%>%
  dplyr::filter(grepl(pattern = "STD",dataOrigin),
                sample.source=="STD")%>%
  dplyr::mutate(groupMz(precursorMz,return.type = "data.frame"))


sp.file<- sp.df$dataOrigin%>%unique()

p.list <- list()
for (i in seq_along(sp.file)) {

  sp.file[[i]]
  this.sp.df <- sp.df%>%
    dplyr::filter(dataOrigin%in% sp.file[[i]])

  p <- this.sp.df%>%
    ggplot()+
    geom_point(aes(x = rtime,
                  y = precursorMz))+
    labs(title = basename(sp.file[[i]]))
  p.list[[i]] <- p

}


gp <- ggplot_sum_patchwork(p.list)

open_plot_win(gp,20,20)
# Sat Aug 24 15:09:14 2024 ------------------------------


msip.core <- msdev.M1@statData[[
  "MSIP"]][[
  "isotopologues_data"]][[
    "FT02170_Negative"]][[
      "MSIP_result"]][["M1"]][["STD"]]
msip.core <- msdev.M1@statData$MSIP$
    isotopologues_data$HMDB0000902_merged$MSIP_result$M1$Con

msip.core@solve <- list()
msip.core <- MSIPCore_solve(msip.core,
                            int_thresh = 10^3.8,
                            certainty_thresh = 0.8)
#heatmap_MSIPFragmentMap(msip.core@FG_map)

plotly_MSIPCore_pred_nature_prob(msip.core)%>%
  ggplotly()

cp.table <- get_MSIP_compound_info(msdev.M1)

cp.vars <- lapply(iso.list,FUN = function(x){
  names(x$compound_info)
})


cp.vars.df <- do.call(rbind,cp.vars)


# Sat Aug 24 22:24:50 2024 ------------------------------
iso.list <- msdev.M1@statData$MSIP$isotopologues_data
iso.data <- iso.list$FT01159_Negative

a <- shiny_get_sp_data(iso.data,"Con","M0")

mz_values <- mz(spectrum)[[1]]
intensity_values <- intensity(spectrum)[[1]]


# Define a noise region, e.g., between m/z 1000 and 1100
noise_region <- intensity_values

# Calculate RMS noise
rms_noise <- sd(noise_region)

# Calculate peak-to-peak noise
peak_to_peak_noise <- max(noise_region) - min(noise_region)

# Average noise level
average_noise <- mean(noise_region)

# Output the results
cat("RMS Noise:", rms_noise, "\n")
cat("Peak-to-Peak Noise:", peak_to_peak_noise, "\n")
cat("Average Noise Level:", average_noise, "\n")

plotly_Spectra(spectrum)


spm <- msdev.M1@xcmsData$NegativeMS1[[1]]
estimateNoise(spm)
spm->spectrum
mz_values <- mz(spectrum)
intensity_values <- intensity(spectrum)

# Choose a noise region (m/z range where you expect no peaks)
noise_region <- intensity_values[mz_values >= 500 & mz_values <= 600]

# Calculate basic noise statistics
rms_noise <- sd(noise_region)
mean_noise <- mean(noise_region)

cat("RMS Noise:", rms_noise, "\n")
cat("Mean Noise:", mean_noise, "\n")




hist(log10(intensity(sp)[[1]]))


plotly_Spectra(sp)





shiny_get_sp_data(iso.data,"Con","M1")%>%
  shiny_plotly_iso_data_spectra()



iso.list <- msdev.M1@statData$MSIP$isotopologues_data
is.merged <- sapply(iso.list,function(x){
  x$compound_info$merged
})
msdev.M1@statData$MSIP$isotopologues_data <- iso.list[!is.merged]


iso.list[[211]]$MSIP_result$M0$Con@FG_map@fragment.atom.matrix


# Mon Aug 26 20:41:09 2024 ------------------------------
df <- data.frame(
  raw = 1.1^(1:200)
)%>%
  dplyr::mutate(weight = .intensity_weight(raw))

ggplot(df)+
  geom_point(aes(x = raw, y = weight))+
  #geom_vline(xintercept = c(1e3,1e4))+
  scale_x_log10()+
  labs(x = "intensity",
       title = "Fragment weight curve")+
  theme_bw()->p
open_plot_win(p)




# Thu Aug 29 10:21:02 2024 MSIP TUNE------------------------------
iso.list <- msdev.M1@statData$MSIP$isotopologues_data
iso.data <- iso.list$FT00684_Negative
msip.core <- iso.data$MSIP_result$M1$Con



# Thu Aug 29 13:15:26 2024 get_MSIPCoreData------------------------------
iso.list <- msdev.M1@statData$MSIP$isotopologues_data
iso.data <- iso.list$FT11987_Negative
msip.core <- iso.data$MSIP_result$M3$Con
msip.core
msip.core <- get_MSIPCoreData(sp.iso = iso.data$Spectra$M1$Con,
                              cfmd = iso.data$CFM_annotation,
                              iso_count = 1,
                              ppm = 10)
msip.core <- MSIPCore_solve(msip.core)
MSIPFragmentMap <- msip.core@FG_map

MSIPFragmentMap <- MSIPFragmentMap_add_constraint(MSIPFragmentMap)
heatmap_MSIPFragmentMap(MSIPFragmentMap)

MSIPCore_solve(msip.core)



# Thu Aug 29 15:06:35 2024 isotopomers set vis------------------------------
msip.frag.map <- msip.core@solve$MSIPIsotopomerMap
possible.set <- msip.frag.map@solve$isotopomer.set[msip.frag.map@solve$isotopomer.set.prob > 1e-5]

for (i_set  in seq_along(possible.set)) {
  msip.frag.map@isotopomer.defination[possible.set[[i_set]]]%>%
    unlist()%>%
    table()%>%
    `/`(sum(.))%>%
    `*`(3)%>%
    print()
}

# Thu Aug 29 17:33:55 2024 eval cor and RMSE------------------------------
process.info <- MSIP_solve_computation_evaluate(msdev.M1)
process.infoL <- process.info%>%
  dplyr::filter(iso_count == 1&samples == "Liver")
#process.info <- process.info[1:10,]
#msdev.M1 <- MSIP_solve_isotopologues(msdev.M1,
#                                     process.info,
#                                     ppm = 20,
#                                     timeout = 15)
process.infoL$cor <- NA
process.infoL$count <- NA
iso.list <- msdev.M1@statData$MSIP$isotopologues_data
for (i in 1:nrow(process.infoL)) {

  msip.core <- iso.list[[process.infoL$feature_id[i]]]$
    MSIP_result[[str_isotope2_num(process.infoL$iso_count[i])]][[process.infoL$samples[i]]]
  if (is.null(msip.core)|isEmpty(msip.core)) {
    next
  }
  msip.ifmap <- msip.core@solve$MSIPIsotopomerMap
  x <- lengths(msip.ifmap@solve$isotopomer.set)

  y <- msip.ifmap@solve$isotopomer.set.prob
  if (length(x)) {
    x <- x/sum(x)
    process.infoL$count[i] <- length(x)
    process.infoL$cor[i]<-  cor(x,y)^2
    process.infoL$rmse[i] <- sqrt(mean((x - y) ^ 2))
  }

}


plot.data <- process.infoL%>%
  dplyr::filter(count > 0)%>%
  dplyr::mutate(meaningful = count>1)

ggplot(plot.data)+
  geom_histogram(aes(x  = rmse,fill = meaningful ),position = "dodge")->p.cor

open_plot_win(p)

# Sat Aug 31 15:10:26 2024 ------------------------------

iso.list <- msdev.M1@statData$MSIP$isotopologues_data
iso.data <- iso.list[[196]]
shiny_vis_cfmd(iso.data$CFM_annotation0)

cfmd <- iso.data$CFM_annotation0
msip.core <- iso.data$MSIP_result$M1$Liver
heatmap_MSIPFragmentMap(msip.core@FG_map)
msip.core <- MSIPCore_solve(msip.core,
                            int_thresh = 8000,
                            certainty_thresh = 0)
plotly_MSIPCore_pred_nature_prob(msip.core)%>%
  open_visNet()



cfmd <- CFM_data_get_atom_map(cfmd)
# Sat Aug 31 16:54:37 2024 MSIP single solve pipline ------------------------------
object <- msdev.M1
iso.list <- object@statData$MSIP$isotopologues_data
iso.data <- iso.list[[9]]
### CFM
{
  cfmd <- get_CFM_data_from_smiles(
    smiles = iso.data$compound_info$smiles,
    compound_id = iso.data$compound_info$compound_id,
    ppm =  10,
    adduct = switch(as.character(iso.data$compound_info$polarity),
                    "0"="[M-H]-",
                    "1"="[M+H]+"),
    check_temp = F,
    temp_dir = paste0(object@projectInfo$CompoundDB_path,"_cfmd"))
  iso.data$CFM_annotation <- cfmd


}

get_CFM_data_MSIPFragmentMap(cfmd)%>%
  heatmap_MSIPFragmentMap()


ggplot(cfmd@fragment_define)+
  geom_point(aes(x = ratio ,y= bond.score))


# Tue Sep  3 18:21:37 2024 evaluate Cor and rmse ------------------------------
### Liver
{
  process.info <- MSIP_solve_computation_evaluate(msdev.M1)
  process.infoL <- process.info%>%
    dplyr::filter(iso_count == 1&samples == "Liver")
  #process.info <- process.info[1:10,]
  #msdev.M1 <- MSIP_solve_isotopologues(msdev.M1,
  #                                     process.info,
  #                                     ppm = 20,
  #                                     timeout = 15)
  process.infoL$cor <- NA
  process.infoL$count <- NA
  iso.list <- msdev.M1@statData$MSIP$isotopologues_data
  for (i in 1:nrow(process.infoL)) {

    msip.core <- iso.list[[process.infoL$feature_id[i]]]$
      MSIP_result[[str_isotope2_num(process.infoL$iso_count[i])]][[process.infoL$samples[i]]]
    if (is.null(msip.core)|isEmpty(msip.core)) {
      next
    }
    msip.ifmap <- msip.core@solve$MSIPIsotopomerMap
    x <- lengths(msip.ifmap@solve$isotopomer.set)

    y <- msip.ifmap@solve$isotopomer.set.prob
    if (length(x)) {
      x <- x/sum(x)
      process.infoL$count[i] <- length(x)
      process.infoL$cor[i]<-  cor(x,y)^2
      process.infoL$rmse[i] <- sqrt(mean((x - y) ^ 2))
      process.infoL$r2[i]<-  summary(lm(x~y))$r.squared



    }

  }


  plot.data <- process.infoL%>%
    dplyr::filter(count > 0)%>%
    dplyr::mutate(meaningful = count>1,
                  solve.ratio = count/target_ele_count)

  p.rmse <- ggplot(plot.data)+
    geom_histogram(aes(x  = rmse,fill = meaningful ),position = "dodge")+
    scale_fill_manual(values = c("TRUE" = "#00BFC4","FALSE" = "#F8766D"))+
    labs(x = "RMSE",y = "Count")+
    theme_bw()

  p.cor <- ggplot(plot.data)+
    geom_histogram(aes(x  = r2,fill = meaningful ),position = "dodge",show.legend = F)+
    scale_fill_manual(values = c("TRUE" = "#00BFC4","FALSE" = "#F8766D"))+
    labs(x = "R²",y = "Count")+
    theme_bw()


  p.stat <- ggplot(plot.data)+
    geom_point(aes(x = solve.ratio, y = count,col = meaningful),show.legend = F)+
    scale_fill_manual(values = c("TRUE" = "#00BFC4","FALSE" = "#F8766D"))+
    labs(x = "Solve ratio",y = "Isotopomers set count")+
    theme_bw()

  p <- p.rmse+p.cor+p.stat+plot_layout(guides="collect")+
    plot_annotation(title = "Natural isotopomers M1 validation in Liver")
  open_plot_win(p, 10, 3)

}


### Cell
{
  process.info <- MSIP_solve_computation_evaluate(msdev.M1)
  process.infoL <- process.info%>%
    dplyr::filter(iso_count == 1&samples == "Con")
  #process.info <- process.info[1:10,]
  #msdev.M1 <- MSIP_solve_isotopologues(msdev.M1,
  #                                     process.info,
  #                                     ppm = 20,
  #                                     timeout = 15)
  process.infoL$cor <- NA
  process.infoL$count <- NA
  iso.list <- msdev.M1@statData$MSIP$isotopologues_data
  for (i in 1:nrow(process.infoL)) {

    msip.core <- iso.list[[process.infoL$feature_id[i]]]$
      MSIP_result[[str_isotope2_num(process.infoL$iso_count[i])]][[process.infoL$samples[i]]]
    if (is.null(msip.core)|isEmpty(msip.core)) {
      next
    }
    msip.ifmap <- msip.core@solve$MSIPIsotopomerMap
    x <- lengths(msip.ifmap@solve$isotopomer.set)

    y <- msip.ifmap@solve$isotopomer.set.prob
    if (length(x)) {
      x <- x/sum(x)
      process.infoL$count[i] <- length(x)
      process.infoL$cor[i]<-  cor(x,y)^2
      process.infoL$rmse[i] <- sqrt(mean((x - y) ^ 2))
      process.infoL$r2[i]<-  summary(lm(x~y))$r.squared
    }

  }


  plot.data <- process.infoL%>%
    dplyr::filter(count > 0)%>%
    dplyr::mutate(meaningful = count>1,
                  solve.ratio = count/target_ele_count)

  p.rmse <- ggplot(plot.data)+
    geom_histogram(aes(x  = rmse,fill = meaningful ),position = "dodge")+
    scale_fill_manual(values = c("TRUE" = "#00BFC4","FALSE" = "#F8766D"))+
    labs(x = "RMSE",y = "Count")+
    theme_bw()

  p.cor <- ggplot(plot.data)+
    geom_histogram(aes(x  = r2,fill = meaningful ),position = "dodge",show.legend = F)+
    scale_fill_manual(values = c("TRUE" = "#00BFC4","FALSE" = "#F8766D"))+
    labs(x = "R²",y = "Count")+
    theme_bw()

  p.stat <- ggplot(plot.data)+
    geom_point(aes(x = solve.ratio, y = count,col = meaningful),show.legend = F)+
    scale_fill_manual(values = c("TRUE" = "#00BFC4","FALSE" = "#F8766D"))+
    labs(x = "Solve ratio",y = "Isotopomers set count")+
    theme_bw()

  p <- p.rmse+p.cor+p.stat+plot_layout(guides="collect")+
    plot_annotation(title = "Natural isotopomers M1 validation in Cell")
  open_plot_win(p, 10, 3)

}


# Tue Sep  3 22:04:40 2024 NMR data ------------------------------
gln.iso.data.mix <- msdev.M1@statData$MSIP$isotopologues_data$FT02220_Negative

ms1.int.matirx <- gln.iso.data.mix$compound_info$ratio_matrix
Heatmap(ms1.int.matirx,
        name = "Ratio",
        col = colramp(),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", ms1.int.matirx[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        row_names_side = "left",
        column_names_rot = 0,
        cluster_columns = F,
        cluster_rows = F
        ) -> p
open_plot_win( p )

nmr.data <- gln.iso.data.mix$MSIP_result$M1$NMR
p <- heatmap_MSIPFragmentMap(nmr.data@FG_map)
open_plot_win(p, 6, 8)
nmr.data@solve$MSIPIsotopomerMap@isotopomer.defination
nmr.data@solve$MSIPIsotopomerMap@solve$isotopomer.set
nmr.data@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob

nmr.data <- gln.iso.data.mix$MSIP_result$M2$NMR
p <- heatmap_MSIPFragmentMap(nmr.data@FG_map)
open_plot_win(p,7,8)

nmr.data <- gln.iso.data.mix$MSIP_result$M3$NMR
p <- heatmap_MSIPFragmentMap(nmr.data@FG_map)
open_plot_win(p,8,8)

vis_sdf_igraph(gln.iso.data.mix$CFM_annotation@fragment_igraph$Fragment01,show_id = T)%>%
  open_visNet()

plot.data = data.frame(
  x = c(0.216,0.2875,0.2875,0.0753,0.066),
  y = c(0.275,0.1718,0.378,0,0.174)
)

ggplot(plot.data)+
  geom_point(aes(x = x,y=y))+
  geom_abline(slope = 1,intercept = 0,col = "grey",linewidth = 2)+
  xlim(c(0,0.5))+
  ylim(c(0,0.5))+
  theme_bw()+
  labs(x = "MSIP", y = "NMR")->p

open_plot_win(p,3.3,3)






# Sat Sep  7 19:39:21 2024 CFM network------------------------------
{
  nodes <- data.frame(id = 1:4, label = c("Node 1", "Node 2", "Node 3", "Node 4"))

  # Create edge data frame with a 'dis' attribute for distance
  edges <- data.frame(
    from = c(1, 2, 3),
    to = c(2, 3, 4),
    length = c(100, 200, 50)  # Example distances
  )

  # Visualize the network and set the edge lengths based on the 'dis' column
  visNetwork(nodes, edges) %>%
    visEdges(smooth = FALSE, scaling = list(min = 10, max = 200)) %>%  # Control scaling of edges if necessary
    #visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -200)) %>%
    #visEdges(length = edges$dis)%>%
    open_visNet()

}

# Tue Sep 10 15:43:19 2024 ------------------------------
{


}

# Tue Sep 10 19:08:57 2024 ------------------------------
smile <-  "CCCC(CCC)C(O)=O"
sdf <- get_smiles_sdf(smile)
sdf.ig <- get_sdf_igraph(sdf)[[1]]

vis_sdf_igraph(sdf.ig)%>%
  open_visNet()







# Wed Sep 11 20:04:58 2024 ------------------------------


nodes <- data.frame(id = 1:3,
                    shape = 'image',
                    image = paste('data:image/png;base64', img.txt, sep = ','),
                    stringsAsFactors = F)
edges <- data.frame(from = c(1,1), to = c(2,3))

visNetwork(nodes, edges) %>%
  visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>%
  visHierarchicalLayout()

# Fri Sep 13 14:46:17 2024 ------------------------------
m1.stat <- get_MSIP_M1_Statistic(msdev.M1)
fragment.stat <- get_MSIP_fragment_Statistic(msdev.M1)

fragment.stat <- fragment.stat%>%
  dplyr::mutate(solve.ratio = is.count/isotopomer)
ggplot(fragment.stat)+
  geom_point(aes(x = target_ele_count,
                 y = fragment.count,
                 size = iso_count,
                 col =  solve.ratio))+
  scale_color_gsea()

fragment.stat$is.count[fragment.stat$is.count>10]<-10
ggplot(fragment.stat)+
  geom_histogram(aes(x = solve.ratio) )

cell.stat <- fragment.stat%>%
  dplyr::filter(samples%in% c("Con" ,   "Glu" ,  "U"  ))
edit_df_in_excel(cell.stat)



this.path <- all_simple_paths(ig.trans,from = 1,to = this.frag,
                              mode = "out",
                              cutoff = 3)
# Fri Sep 20 19:57:15 2024 GLN STD------------------------------
iso.data.list <- msdev.M1@statData$MSIP$isotopologues_data
iso.data <- iso.data.list$FT02170_Negative
iso.data$Spectra$M1$STD$dataOrigin%>%table()



cfmd <- get_CFM_data_from_smiles(iso.data$compound_info$smiles,
                                 adduct = "[M-H]-",
                                 check_temp = F)
msip.core <- get_MSIPCoreData(sp.iso = iso.data$Spectra$M1$STD,
                              cfmd = cfmd,
                              iso_count = 1,
                              ppm = 20)
msip.core <- MSIPCore_solve(msip.core,
                            int_thresh = 10^4,
                            re_split_isotopomers = T)
plotly_MSIPCore_pred_nature_prob(msip.core)
heatmap_MSIPFragmentMap(msip.core@FG_map,T)


# Sat Sep 21 14:15:18 2024 cfmd fragmetn0------------------------------
gln.smiles <- "N[C@@H](CCC(N)=O)C(O)=O"
cfmd0 <- get_CFM_data_from_smiles(gln.smiles,
                                  adduct = "[M-H]-",
                                  check_temp = F)

cfmd1<- get_CFM_data_from_smiles(gln.smiles,
                                 adduct = "[M+H]+",
                                 check_temp = F)

vis_sdf_igraph(cfmd0@fragment_igraph$Fragment00,show_id = T)






# Mon Sep 23 08:16:27 2024 M1 select for Liver------------------------------
a <- featureDefinitions(xcms.chrom.data$Positive_Chromatograms)
xcms.chrom.data$Positive_Chromatograms@featureDefinitions <-  featureDefinitions(xcms.chrom.data$Positive_Chromatograms)[,1:39]
xcms.chrom.data$Negative_Chromatograms@featureDefinitions <-  featureDefinitions(xcms.chrom.data$Negative_Chromatograms)[,1:39]
msdev.liver@xcmsData$Positive_Chromatograms <- xcms.chrom.data$Positive_Chromatograms
msdev.liver@xcmsData$Negative_Chromatograms <- xcms.chrom.data$Negative_Chromatograms



xcms.chrom.data <- object@xcmsData[c("Positive_Chromatograms",
                                       "Negative_Chromatograms")]

xcms.chrom.data <- lapply(xcms.chrom.data,
                          onDiskData_retrieve  )

coln <- xcms.chrom.data$Positive_Chromatograms@featureDefinitions%>%colnames()
idx <- grepl("Ratio_to_seed_FT",coln)
xcms.chrom.data$Positive_Chromatograms@featureDefinitions[,idx] <- NULL
xcms.chrom.data$Positive_Chromatograms[1,1,drop = F]%>%featureDefinitions()


coln <- xcms.chrom.data$Negative_Chromatograms@featureDefinitions%>%colnames()
idx <- grepl("Ratio_to_seed_FT",coln)
xcms.chrom.data$Negative_Chromatograms@featureDefinitions[,idx] <- NULL
xcms.chrom.data$Negative_Chromatograms[1,1,drop = F]%>%featureDefinitions()




# Tue Sep 24 15:19:39 2024 ------------------------------
MSconvertR::msConvertDir("d:/temp/")
xcms.xcms <- readMSData("d:/temp/MCE-pos.mzML",mode = "onDisk")
xcms.scans <- get_xcms_scan_Stat(xcms.xcms)
xcms.sp <- get_xcms_Spectra(xcms.xcms)

xcms.sp.ms2 <- xcms.sp[msLevel(xcms.sp)==2]
plotly_Spectra(xcms.sp.ms2[100])

plot.data <- xcms.scans%>%
  dplyr::distinct(ms1_no,ms2_count,cycle_time,.keep_all = T)

p <- ggplot(plot.data)+
  geom_point(aes(x = retentionTime,
                 y = ms2_count,
                 color = cycle_time))+
  scale_color_gradient(low = "yellow",high = "red")+
  theme_bw()
open_plot_win(p,5,3)

xcms.sp.ms2 <- Spectra_get_noise(xcms.sp.ms2)


# Wed Sep 25 11:48:49 2024 ------------------------------
{

  x <- msdev.liver@statData$MSIP$isotopologues_table$Positive
  id <- grepl("Ratio_to_seed_FT",colnames(x))
  msdev.liver@statData$MSIP$isotopologues_table$Positive <- x[,!id]

  x <- msdev.liver@statData$MSIP$isotopologues_table$Negative
  id <- grepl("Ratio_to_seed_FT",colnames(x))
  msdev.liver@statData$MSIP$isotopologues_table$Negative <- x[,!id]
}

{

  x <- msdev.liver@xcmsData$PositiveMS1
  id <- grepl("Ratio_to_seed_FT",colnames(featureDefinitions(x)))
  featureDefinitions(x) <- featureDefinitions(x)[,!id]
  x -> msdev.liver@xcmsData$PositiveMS1

  x <- msdev.liver@xcmsData$NegativeMS1
  id <- grepl("Ratio_to_seed_FT",colnames(featureDefinitions(x)))
  featureDefinitions(x) <- featureDefinitions(x)[,!id]
  x -> msdev.liver@xcmsData$NegativeMS1
}


x <- msdev.liver@xcmsData$Positive_Chromatograms
x.data <- onDiskData_retrieve(x)
id <- grepl("Ratio_to_seed_FT",colnames(x.data@featureDefinitions))
x.data@featureDefinitions <- x.data@featureDefinitions[,!id]
x <- onDiskData_update(x,x.data)


# Thu Sep 26 11:01:26 2024 ------------------------------
a <- apply(xcms.ratio.to.seed ,1,
      function(x){  mean_f(x , f = c("A","A","B"),
                           simplify = F,na.rm=T)})%>%
  do.call(rbind,.)

b <- apply(xcms.ratio.to.seed ,1,
           function(x){  mean_f(x , f = c("A","A","A"),simplify =F,na.rm=T)})%>%
  do.call(rbind,.)
dim(a)
dim(b)

# Thu Sep 26 15:22:47 2024 GLN STD------------------------------
##neg
{
  iso.data.list <- msdev.M1@statData$MSIP$isotopologues_data
  iso.data <- iso.data.list$FT02170_Negative
  iso.data$Spectra$M1$STD$dataOrigin%>%table()



  cfmd <- get_CFM_data_from_smiles(iso.data$compound_info$smiles,
                                   adduct = "[M-H]-",
                                   check_temp = F)
  msip.core <- get_MSIPCoreData(sp.iso = iso.data$Spectra$M1$STD,
                                cfmd = cfmd,
                                iso_count = 1,
                                ppm = 20)
  msip.core <- MSIPCore_solve(msip.core,
                              int_thresh = 10^3.8,
                              certainty_thresh = 0.5,
                              re_split_isotopomers = T)
  sprintf("%2f",msip.core@solve$Atom_prob)
  fg <- MSIPFragmentMap_filter_intensity(msip.core@FG_map,10^3)%>%
    MSIPFragmentMap_include_fragment
  heatmap_MSIPFragmentMap(fg,T)

}
##pos
{
  {
    iso.data.list <- msdev.M1@statData$MSIP$isotopologues_data
    iso.data <- iso.data.list$FT03168_Positive
    iso.data$Spectra$M1$STD$dataOrigin%>%table()



    cfmd <- get_CFM_data_from_smiles(iso.data$compound_info$smiles,
                                     adduct = "[M+H]+",
                                     check_temp = F)

    msip.core <- get_MSIPCoreData(sp.iso = iso.data$Spectra$M1$STD,
                                  cfmd = cfmd,
                                  iso_count = 1,
                                  ppm = 20)
    msip.core <- MSIPCore_solve(msip.core,
                                int_thresh = 10^3.8,
                                certainty_thresh = 0,
                                re_split_isotopomers = T)
    sprintf("%2f",msip.core@solve$Atom_prob)
    fg <- MSIPFragmentMap_filter_intensity(msip.core@FG_map,10^3.8)%>%
      MSIPFragmentMap_include_fragment
    heatmap_MSIPFragmentMap(fg,T)
  }
}


### Check for FG map
{
  vis_cfm_data_fragment_atom_map(cfmd,
                                 "Fragment05")%>%
    open_visNet()
  vis_cfm_data_trans_map(cfmd,37)
  vis_cfm_data_fragment(cfmd,1)
}


### merged map
{

  iso.data.list <- msdev.M1@statData$MSIP$isotopologues_data
  iso.data <- iso.data.list$HMDB0000641_merged
  msip.core <- iso.data$MSIP_result$M1$STD

  fg <- msip.core@FG_map%>%
    MSIPFragmentMap_filter_intensity(10^4)%>%
    MSIPFragmentMap_include_fragment()
  p <- heatmap_MSIPFragmentMap(fg,show_ratio = T)
  open_plot_win(p,10,10)
}

# Thu Sep 26 18:10:12 2024 colramp------------------------------

library(wesanderson)

ggplot(mtcars, aes(x = wt, y = mpg, color = hp)) +
  geom_point(size = 4) +
  scale_color_gradientn(colors = wes_palette("Moonrise2", 100,
                                             type = "continuous"))


wes_palette("GrandBudapest2", 100,
            type = "continuous")

library(scico)

ggplot(mtcars, aes(x = wt, y = mpg, color = hp)) +
  geom_point(size = 4) +
  scale_color_scico(palette = "acton")





# Thu Sep 26 18:27:55 2024 Glutamate compare to NMR------------------------------
cp.id <- "HMDB0000148"

fdf <- msdev.13C1@statData$MSIP$isotopologues_table$Negative%>%
  dplyr::filter(iso_seed =="FT01908")
ratio.matrix <- msdev.13C1@statData$MSIP$isotopologues_matrix$ratio_to_seed$Negative[fdf$feature_id,]


ratio.matrix <- apply(ratio.matrix,1,
                      mean_f,f = paste0(rep(0:3,each = 3),"-Glu"))%>%t
rownames(ratio.matrix) <- paste0("M",fdf$iso_count)
Heatmap(ratio.matrix,
        col = colramp(),
        cell_fun =  function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", ratio.matrix[i, j]), x, y, gp = gpar(col = "black", fontsize = 10))
        },show_heatmap_legend = F,
        cluster_columns = F,
        cluster_rows = F,
        row_names_side = "left"
        )->p
#open_plot_win(p,5,5)


{

  nmr.data <- readxl::read_excel("C:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/2024.09.26.GLU_13C1_NMR.xlsx")%>%
    column_to_rownames("...1")%>%
    apply(1,function(x)x/sum(x))%>%t
  id <- colnames(nmr.data)


  iso.data <- msdev.13C1@statData$MSIP$isotopologues_data$HMDB0000148_merged
  samples <- names(iso.data$MSIP_result$M0)
  prob.list <- list()
  atom.prob.matrix.list <- list()
  int <- 10^3
  cer = 0.6
  for (i in 1:3) {
    this.sample <- samples[i]

    msip.core.m1 <- iso.data$MSIP_result$M1[[this.sample]]
    msip.core.m1 <- MSIPCore_solve(msip.core.m1,int_thresh = int,
                                   certainty_thresh = cer,
                                   weight_fun = .intensity_weight )
    msip.core.m2 <- iso.data$MSIP_result$M2[[this.sample]]
    msip.core.m2 <- MSIPCore_solve(msip.core.m2,int_thresh = int,
                                   certainty_thresh = cer,
                                   weight_fun = .intensity_weight )
    msip.core.m3 <- iso.data$MSIP_result$M3[[this.sample]]
    msip.core.m3 <- MSIPCore_solve(msip.core.m3,int_thresh = int,
                                   certainty_thresh = cer,
                                   weight_fun = .intensity_weight )

    atom.prob.matrix <- data.frame(
      M1 = msip.core.m1@solve$Atom_prob,
      M2= msip.core.m2@solve$Atom_prob,
      M3= msip.core.m3@solve$Atom_prob)%>%
    #  apply(2,function(x){
    #  x/sum(x)
   # })%>%
      t
    atom.prob.matrix.list[[i]] <- atom.prob.matrix%>%
        apply(1,function(x){
        x/sum(x)
       })%>%t
    #atom.prob.matrix <- atom.prob.matrix.list[[i]][,id]
    x <- apply(atom.prob.matrix,2,
          weighted.mean,
          w = ratio.matrix[paste0("M",1:3),paste0(i,"-Glu")] )
    x <- x[id]/sum(x)
    y <- nmr.data[i,id]
    prob.list[[i]] <-
      data.frame(tracer = paste0(i,"-Glu"),
                 value = c(x,y),
                 source = rep(c("MSIP","NMR"),each = 5))

    prob.list[[i]] <-
      data.frame(tracer = paste0(i,"-Glu"),
                 label = id,
                 MSIP = x,
                 NMR = y)
  }

  prob.df <- prob.list%>%
  do.call(rbind,.)%>%
    dplyr::filter(tracer!= "3-Glu")
  R2 <- summary(lm(prob.df$MSIP~prob.df$NMR))$r.squared

  ggplot(prob.df)+
    geom_point(aes(x = NMR, y =  MSIP,colour = tracer),size = 5)+
    geom_abline(slope = 1)+
    ggrepel::geom_text_repel(aes(x = NMR, y =  MSIP,label = label))+
    ggsci::scale_color_npg()+
    labs(title = R2)+
    xlim(c(0,1))+
    ylim(c(0,1))+
    theme_bw()->p
  print(p)
  open_plot_win(p,6,5)


}

### Heatmap atom prob
{

  rownames(nmr.data) <- paste0(rownames(nmr.data),"-NMR\nTotal Isotopomers")
  glu1.matrix <- atom.prob.matrix.list[[1]][,id]%>%
    `rownames<-`(paste0("GLU1-",rownames(.)))
  glu2.matrix <- atom.prob.matrix.list[[2]][,id]%>%
    `rownames<-`(paste0("GLU2-",rownames(.)))
  glu3.matrix <- atom.prob.matrix.list[[3]][,id]%>%
    `rownames<-`(paste0("GLU3-",rownames(.)))
  total.matrix <- rbind(
    nmr.data[1,id,drop = F],
    glu1.matrix[,id,drop = F],
    nmr.data[2,id,drop = F],
    glu2.matrix[,id,drop = F],
    nmr.data[3,id,drop = F],
    glu3.matrix[,id,drop = F]
  )
  Heatmap(total.matrix,
          col = colramp(),
          row_split = c(10,11,11,11,20,22,22,22,30,33,33,33),
          cluster_row_slices = F,
          cluster_rows = F,
          cluster_columns = F,
          row_title = NULL,
          rect_gp = gpar(color = "grey"),
          row_names_side = "left")->p
  open_plot_win(p,5,7)

}

# Fri Sep 27 13:13:18 2024 Astra MS2 evaluation------------------------------
{
  xcms.qe <- readMSData("d:/temp/GLU1_PRM_NEG_CE10.mzML",mode = "onDisk")
  qe.scans <- get_xcms_scan_Stat(xcms.qe)
  qe.sp <- get_xcms_Spectra(xcms.qe)
  qe.sp.ms2 <- qe.sp[msLevel(qe.sp)==2]
  qe.sp.ms2 <- Spectra_get_noise(qe.sp.ms2)

  plotly_Spectra(qe.sp.ms2[100])

  xcms.astral <- readMSData("d:/temp/MCE-pos.mzML", mode = "onDisk" )
  astral.scans <- get_xcms_scan_Stat(xcms.astral)
  astral.sp <- get_xcms_Spectra(xcms.astral)
  astral.sp.ms2 <- astral.sp[msLevel(astral.sp)==2]
  astral.sp.ms2 <- Spectra_get_noise(astral.sp.ms2)

  plotly_Spectra(astral.sp.ms2[100])

  plot.data <- astral.scans%>%
    dplyr::distinct(ms1_no,ms2_count,cycle_time,.keep_all = T)

  p <- ggplot(plot.data)+
    geom_point(aes(x = retentionTime,
                   y = ms2_count,
                   color = cycle_time))+
    #scale_color_gradient(low = "yellow",high = "red")+
    scico::scale_color_scico(direction = -1)+
    theme_bw()
  open_plot_win(p,5,3)


  ### MS2 scan time compare to QE
  plot.data <- rbind(
    qe.scans%>%
      dplyr::mutate(instrument = "QE"),
    astral.scans%>%
      dplyr::mutate(instrument = "Astra")%>%
      dplyr::filter(msLevel==2)
  )
  ggplot(plot.data)+
    geom_boxplot(aes(x = instrument , y = scan_time,colour = instrument),show.legend = F)+
    ylim(c(0,0.5))+
    ggsci::scale_color_aaas()+
    theme_bw()->p.scan.time
  p.scan.time

  ggplot(plot.data)+
    geom_bar(aes(x = instrument ,fill = instrument))+
    #ylim(c(0,0.5))+
    ggsci::scale_fill_aaas()+
    labs(y = "MS2 count")+
    theme_bw()->p.ms2.count
  p.ms2.count

  p <- p.scan.time+p.ms2.count+
    plot_layout(guides = "collect")
  open_plot_win(p)



   ### noise
  sp.data <- rbind(
    spectraData(qe.sp.ms2)%>%
      as.data.frame()%>%
      dplyr::mutate(instrument = "QE"),
    spectraData(astral.sp.ms2)%>%
      as.data.frame()%>%
      dplyr::mutate(instrument = "Astra")%>%
      dplyr::filter(msLevel==2)
  )%>%
    dplyr::mutate(TNR = totIonCurrent/noise)

  ggplot(sp.data)+
    geom_point(aes(x = log10(totIonCurrent),
                   y = log10(noise),
                   color = instrument),alpha = 0.2)+
    ggsci::scale_color_aaas()+
    #labs(y = "MS2 count")+
    theme_bw()->p.noise.tic


  open_plot_win(p.noise.tic,5,4)

  ggplot(sp.data)+
    geom_violin(aes(x = instrument,
                    y = log10(totIonCurrent),
                    color = instrument),alpha = 0.2)+
    ggsci::scale_color_aaas()+
    #labs(y = "MS2 count")+
    theme_bw()->p.tic
  p.tic
  ggplot(sp.data)+
    geom_violin(aes(x = instrument,
                    y = log10(noise),
                    color = instrument),alpha = 0.2)+
    ggsci::scale_color_aaas()+
    #labs(y = "MS2 count")+
    theme_bw()->p.noise
  p.noise
  ggplot(sp.data)+
    geom_violin(aes(x = instrument,
                   y = log10(TNR),
                   color = instrument),alpha = 0.2)+
    ggsci::scale_color_aaas()+
    #labs(y = "MS2 count")+
    theme_bw()->p.tnr
  p.tnr
  p <- p.noise+p.tic+p.tnr+
    plot_layout(guides = "collect")
  open_plot_win(p,6,4)




}
# Sat Sep 28 13:47:49 2024 ------------------------------
xcms.xcms <- msdev.Astral@sampleInfo$msData.files[1]%>%
  readMSData(mode = "onDisk")
xcms.scans <- get_xcms_scan_Stat(xcms.xcms)
xcms.scans.ms1 <- xcms.scans%>%
  dplyr::filter(msLevel==1)

ggplot(xcms.scans.ms1)+
  geom_point(aes(x = retentionTime, y = cycle_time))


# Sun Sep 29 18:43:16 2024 ------------------------------
all.path <- lapply(1:199,
                   function(x){
                     all_shortest_paths(ig.trans,1,x,mode = "out")$vpath
                   })

all.path[54]

cfmd@fragment_group%>%
  ggplot()+
  geom_jitter(
    aes(x = fragment_count,
                 y= certainty),
    color = "red",
    size = 3,
    alpha = 0.2)


get_CFM_data_from_smiles(smiles = "OC(=O)[C@@H]1CCC(=O)N1",
                         adduct = "[M-H]-",
                         check_temp = T)

msdev.M1@statData$MSIP$isotopologues_data <- list()
# Wed Oct  2 13:39:19 2024 debug for atm error------------------------------
get_CFM_data_from_smiles()

cfmd.files <- dir(
  "c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb_cfmd/",
  pattern = "cfmd",full.names = T
)

cfmd <- readRDS(cfmd.files[1])
cfmd <- CFM_data_get_igraph(cfmd)
cfmd <- CFM_data_get_atom_map(cfmd,iso_ele = iso_ele)
cfmd <- cfm_data_get_FG_map(cfmd,iso_ele = iso_ele)



sdf_igraph_merge(sdf.igraphA =ig.parent,sdf.igraphB =ig.product
)%>%vis_sdf_igraph(show_id = T)%>%
  open_visNet()



system.time(mcs <- fmcsR::fmcs(sdf.parent,sdf.product,bu = 10))
system.time(mcs <- fmcsR::fmcs(sdf.parent,sdf.product,bu = 10,fast = T))


### this compound leat to cfm annotation error
hmdb.id <- "HMDB0001206"
hmdb.df <- MSdb:::get_HMDB_Compound_DF()
idx <- match(hmdb.id,hmdb.df$accession)
smiles <- hmdb.df$smiles[idx]

a <- get_CFM_data_from_smiles(smiles ,adduct = "[M+H]+")


# Sat Oct  5 10:22:13 2024 ATOM debug------------------------------
{
  cfmd <- readRDS("c:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb_cfmd/cfmd.temp.HMDB0001206.[M+H]+.rds")
  cfmd <- CFM_data_get_igraph(cfmd)
  cfmd <- CFM_data_get_atom_map(cfmd)


}

# Sat Oct  5 10:05:19 2024 Astral noise------------------------------
{
  iso.list <- msdev.Astral@statData$MSIP$isotopologues_data
  i <- 10
  iso.data <- iso.list[[i]]
  message(iso.data$compound_info$name," ",iso.data$compound_info$adduct)

  plotly_Spectra(iso.data$Spectra$M0$Glu2)%>%
    open_visNet()

}




# Thu Sep 26 18:27:55 2024 Glutamate compare to NMR Astral------------------------------
cp.id <- "HMDB0000148"

fdf <- msdev.Astral@statData$MSIP$isotopologues_table$Negative%>%
  dplyr::filter(iso_seed =="FT0845")
ratio.matrix <- msdev.Astral@statData$MSIP$isotopologues_matrix$ratio_to_seed$Negative[fdf$feature_id,]


ratio.matrix <- apply(ratio.matrix,1,
                      mean_f,f = paste0(rep(1:3,each = 3),"-Glu"))%>%t
rownames(ratio.matrix) <- paste0("M",fdf$iso_count)
Heatmap(ratio.matrix,
        col = colramp(),
        cell_fun =  function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", ratio.matrix[i, j]), x, y, gp = gpar(col = "black", fontsize = 10))
        },show_heatmap_legend = F,
        cluster_columns = F,
        cluster_rows = F,
        row_names_side = "left"
)->p
open_plot_win(p,5,5)


{

  nmr.data <- readxl::read_excel("C:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/2024.09.26.GLU_13C1_NMR.xlsx")%>%
    column_to_rownames("...1")%>%
    apply(1,function(x)x/sum(x))%>%t
  id <- colnames(nmr.data)


  iso.data <- msdev.Astral@statData$MSIP$isotopologues_data$FT0845_Negative
  vis_sdf_igraph(iso.data$CFM_annotation@fragment_igraph[[1]],show_id = T)
  samples <- names(iso.data$MSIP_result$M0)
  prob.list <- list()
  atom.prob.matrix.list <- list()
  int = 10^2
  cer = 0.8
  #plot_MSIPCore_solve_weight_fun(.intensity_weight_astral)
  {

    for (i in 1:3) {
    this.sample <- samples[i]
    msip.core.m1 <- iso.data$MSIP_result$M1[[this.sample]]
    msip.core.m1 <- MSIPCore_solve(msip.core.m1,int_thresh = int,
                                   certainty_thresh = cer,
                                   weight_fun = .intensity_weight_astral )
    msip.core.m2 <- iso.data$MSIP_result$M2[[this.sample]]
    msip.core.m2 <- MSIPCore_solve(msip.core.m2,int_thresh = int,
                                   certainty_thresh = cer,
                                   weight_fun = .intensity_weight_astral )


    atom.prob.matrix <- data.frame(M1 = msip.core.m1@solve$Atom_prob,
                                   M2 = msip.core.m2@solve$Atom_prob
                                   #M3= iso.data$MSIP_result$M3[[this.sample]]@solve$Atom_prob
                                   )%>%
      #  apply(2,function(x){
      #  x/sum(x)
      # })%>%
      t
    atom.prob.matrix.list[[i]] <- atom.prob.matrix%>%
      apply(1,function(x){
        x/sum(x)
      })%>%t
    #atom.prob.matrix <- atom.prob.matrix.list[[i]]
    x <- apply(atom.prob.matrix,2,
               weighted.mean,
               w = ratio.matrix[paste0("M",1:2),paste0(i,"-Glu")] )
    x <- x[id]
    y <- nmr.data[i,id]
    prob.list[[i]] <-
      data.frame(tracer = paste0(i,"-Glu"),
                 value = c(x,y),
                 source = rep(c("MSIP","NMR"),each = 5))

    prob.list[[i]] <-
      data.frame(tracer = paste0(i,"-Glu"),
                 label = id,
                 MSIP = x,
                 NMR = y)
  }

  prob.df <- prob.list%>%
    do.call(rbind,.)%>%
   dplyr::filter(tracer!= "3-Glu",
                 !label %in% c("C_5","C_2")
                 )%>%
    dplyr::group_by(tracer)%>%
    dplyr::mutate(MSIP = MSIP/sum(MSIP),
                  NMR = NMR/sum(NMR))
  R2 <- summary(lm(prob.df$MSIP~prob.df$NMR))$r.squared
  ggplot(prob.df)+
    geom_abline(slope = 1,lty = 2, col = "grey")+
    geom_point(aes(x =NMR , y =  MSIP,colour = tracer),size = 5)+
    ggrepel::geom_text_repel(aes(x = NMR, y =  MSIP,colour = tracer,label = label))+
    ggsci::scale_color_npg()+
    xlim(c(0,1))+
    ylim(c(0,1))+
    labs(title = R2)+
    theme_bw()->p
  print(p)
  open_plot_win(p,6,5)
  }


}

### Heatmap atom prob
{

  rownames(nmr.data) <- paste0(rownames(nmr.data),"-NMR\nTotal Isotopomers")
  glu1.matrix <- atom.prob.matrix.list[[1]][,id]%>%
    `rownames<-`(paste0("GLU1-",rownames(.)))
  glu2.matrix <- atom.prob.matrix.list[[2]][,id]%>%
    `rownames<-`(paste0("GLU2-",rownames(.)))
  glu3.matrix <- atom.prob.matrix.list[[3]][,id]%>%
    `rownames<-`(paste0("GLU3-",rownames(.)))
  total.matrix <- rbind(
    nmr.data[1,id,drop = F],
    glu1.matrix[,id,drop = F],
    nmr.data[2,id,drop = F],
    glu2.matrix[,id,drop = F],
    nmr.data[3,id,drop = F],
    glu3.matrix[,id,drop = F]
  )
  Heatmap(total.matrix,
          col = colramp(),
          row_split = c(10,11,11,11,20,22,22,22,30,33,33,33),
          cluster_row_slices = F,
          cluster_rows = F,
          cluster_columns = F,
          row_title = NULL,
          rect_gp = gpar(color = "grey"),
          row_names_side = "left")->p
  open_plot_win(p,5,7)

}

# Mon Oct  7 23:42:06 2024 check QE glutamic data------------------------------
{





  m1.sp <- iso.data$Spectra$M1$FS_1_13C%>%
    filterPolarity(1)
  m1.sp.df <- spectraData(m1.sp)%>%
    as.data.frame()%>%
    dplyr::filter(collisionEnergy == 10)

  plotly_Spectra_mirror(m1.sp[m1.sp.df$sp_id[1]],
                    m1.sp[m1.sp.df$sp_id[2]])


  ggplot(m1.sp.df)+
    geom_point(aes(x = rtime,
                   y = log10(totIonCurrent),
                   col = collisionEnergy))



}

# Thu Oct 10 14:20:01 2024 13C1 isotopomers data NMR and QE------------------------------
{

  msdev.13C1 <- load_as_var(
    "C:/Users/91879/OneDrive/Code/R/data/MSIP_data/240701_FS_ONE_POSITION/MSdev_2024_07_04.Rdata"
  )

  iso.data <- msdev.13C1@statData$MSIP$isotopologues_data$HMDB0000148
  fdf <- msdev.13C1@statData$MSIP$isotopologues_table$Negative%>%
    dplyr::filter(iso_seed =="FT01908")
  ratio.matrix <- msdev.13C1@statData$MSIP$isotopologues_matrix$ratio_to_seed$Negative[fdf$feature_id,]


  ratio.matrix <- apply(ratio.matrix,1,
                        mean_f,f = paste0(rep(0:3,each = 3),"-Glu"))%>%t
  rownames(ratio.matrix) <- paste0("M",fdf$iso_count)

  {

    ### NMR data
    {
      nmr.data <- readxl::read_excel("c:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/2024.09.26.GLU_13C1_NMR.xlsx",sheet = 4)
      nmr.data <- nmr.data%>%
        dplyr::group_by(isotopomer)%>%
        dplyr::mutate(Glu1 = sum(Glu1),
                      Glu2 = sum(Glu2))%>%
        dplyr::ungroup()%>%
        dplyr::distinct(isotopomer,Glu1,.keep_all = T)%>%
        dplyr::group_by(atom)%>%
        dplyr::mutate(ratio_glu1 = Glu1/sum(Glu1),
                      ratio_glu2 = Glu2/sum(Glu2))%>%
        dplyr::ungroup()

    }

    ### MSIP formate
    {
      samples <- c("FS_1_13C","FS_2_13C","FS_3_13C")
      samples.name <- make_vector(c("1-Glu","2-Glu","3-Glu"),samples)
      samples.list <- list()
      for (i.sample in samples) {

        lapply(c("M1","M2","M3"), function(i.m) {
          i.atoms <- iso.data$MSIP_result[[i.m]][[i.sample]]@solve$MSIPIsotopomerMap@isotopomer.defination
          i.names <- i.atoms%>%
            sapply(function(x){
              x%>%
                sub("_","",.)%>%
                sub("8","1",.)%>%
                paste0(collapse = ",")
            })
          i.code <- i.atoms%>%
            sapply(function(x){
              x <- x%>%
                sub("_","",.)%>%
                sub("8","1",.)
              z <- make_vector(0,paste0("C",1:5))
              z[x]<-1
              paste0(z,collapse = "")
            })
          i.matrix <- i.atoms%>%
            sapply(function(x){
              x <- x%>%
                sub("_","",.)%>%
                sub("8","1",.)
              z <- make_vector(0,paste0("C",1:5))
              z[x]<-1
              z
            })%>%t
          i.prob <- iso.data$MSIP_result[[i.m]][[i.sample]]@solve$MSIPIsotopomerMap@isotopomer.probability
          i.ratio <- ratio.matrix[i.m,samples.name[i.sample]]
          data.frame(isotopomer = i.names,
                     isotopologue = i.m,
                     code = i.code,
                     prob = i.prob,
                     value = i.prob*i.ratio,
                     i.matrix)%>%
            remove_rownames()
        })%>%
          data.table::rbindlist()->samples.list[[i.sample]]


      }



    }


    ### get data in NMR formate
    {
      nmr.code.table <-readxl::read_excel("c:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/2024.09.26.GLU_13C1_NMR.xlsx",sheet = 5)

      msip.nmr.list <- lapply(samples.list,function(msip.data){
        sum.list <- list()
        for (i in 1:nrow(nmr.code.table)) {

          idx <- which(str_sub(msip.data$code,start = nmr.code.table$char.start[i],
                               end = nmr.code.table$char.end[i])==nmr.code.table$code[i])
          val <- msip.data$value[idx]
          if (is.null(val))  val <- 0
          data.frame(
            atom = nmr.code.table$atom[i],
            isotopomer = nmr.code.table$isotopomer[i],
            code = nmr.code.table$code[i],
            val = sum(val)
          )->sum.list[[i]]

        }

        sum.df <- do.call(rbind,sum.list)%>%
          dplyr::group_by(isotopomer)%>%
          dplyr::mutate(val = sum(val))%>%
          dplyr::ungroup()%>%
          dplyr::distinct(isotopomer,val,.keep_all = T)%>%
          dplyr::group_by(atom)%>%
          dplyr::mutate(ratio = val/sum(val))%>%
          dplyr::ungroup()
        return(sum.df)

      })


      nmr.data <- nmr.data%>%
        dplyr::mutate(
          msip_glu1 =msip.nmr.list[[1]]$ratio[
            match(isotopomer,msip.nmr.list[[1]]$isotopomer)
          ],
          msip_glu2 =msip.nmr.list[[2]]$ratio[
            match(isotopomer,msip.nmr.list[[2]]$isotopomer)
          ]
        )


    }


    ### plot
    {

      ### glu1
      {
        plot.data <- nmr.data%>%
          tidyr::pivot_longer(c(ratio_glu1,msip_glu1))%>%
          dplyr::mutate(label = paste0(isotopomer,"/",code),
                        method = case_when(name=="ratio_glu1"~"NMR",
                                           name =="msip_glu1"~"MSIP"),
                        method = factor(method,level = c("NMR","MSIP")))

        ggplot(plot.data)+
          geom_bar(aes(x = label, y = value,col = method),
                   fill = "white",width =0.6,linewidth = 1,
                   stat = "identity",position =  position_dodge(width = 0.8))+
          scale_color_manual(values = c("MSIP"="#FF9D00", "NMR" ="#2871FF"))+
          scale_y_continuous(expand = expansion())+
          theme_classic()+
          theme(axis.text.x = element_text(angle = 30,hjust = 1))+
          labs(title = "Glu1",x = "",y = "Fraction")->p1
      }

      ### glu2
      {
        plot.data <- nmr.data%>%
          tidyr::pivot_longer(c(ratio_glu2,msip_glu2))%>%
          dplyr::mutate(label = paste0(isotopomer,"/",code),
                        method = case_when(name=="ratio_glu2"~"NMR",
                                           name =="msip_glu2"~"MSIP"),
                        method = factor(method,level = c("NMR","MSIP")))

        ggplot(plot.data)+
          geom_bar(aes(x = label, y = value,col = method),
                   fill = "white",width =0.6,linewidth = 1,
                   stat = "identity",position =  position_dodge(width = 0.8))+
          scale_color_manual(values = c("MSIP"="#FF9D00", "NMR" ="#2871FF"))+
          scale_y_continuous(expand = expansion())+
          theme_classic()+
          theme(axis.text.x = element_text(angle = 30,hjust = 1))+
          labs(title = "Glu2",x = "",y = "Fraction")->p2
      }
      p<-p1+p2+plot_layout(guides = "collect")
      open_plot_win(p,width = 10,height = 3)
    }


  }

}

# Fri Oct 11 19:44:06 2024 ------------------------------
{
  vis_sdf_igraph(sdf.ig[[7]])
}
# Sat Oct 12 09:06:32 2024 ------------------------------
ms2.sp <- get_MSdev_ms2_Spectra(msdev.Astral)
ms2.info <- spectraData(ms2.sp)%>%
  as.data.frame()


ms2.info.pos <- ms2.info%>%
  dplyr::filter(polarity==1)
xcms.fdf.pos <- featureDefinitions(msdev.Astral@xcmsData$PositiveMS1)%>%
  as.data.frame()

sum(length(unlist(xcms.fdf.pos$ms2_id)))
nrow(ms2.info.pos)

# Sun Oct 13 13:59:59 2024 ------------------------------
MSIP_get_isotopologues_table()
a <- msdev.Astral@statData$MSIP$isotopologues_matrix
b <- msdev.Astral@statData$MSIP$isotopologues_table
c <- msdev.Astral@statData$MSIP$isotopologues_data


MSIP_get_isotopologues


# Sun Oct 13 16:03:36 2024 ------------------------------

cfmd <- msdev.13C1@statData$MSIP$isotopologues_data[[1]]$CFM_annotation
vis_cfm_data_fragment(cfmd,1,show_id = T)

#rcdk
{

  # Get all atoms in the molecule
  atoms <- get.atoms(molecule)

  # Loop through each atom in the molecule
  for (i in seq_along(atoms)) {
    atom <- atoms[[i]]

    # Get the symbol for this atom
    atom_symbol <- get.symbol(atom)
    # Check if the atom is Carbon (C)
    if (atom_symbol == "C") {
      # Print the atom index and symbol
      print(paste("Atom index:", i, "Symbol:", atom_symbol))

      # Get the bonds of the atom from the molecule

    }
  }

}

img <- view.image.2d(parse.smiles(smiles)[[1]])
plot(1:10, 1:10, pch=19)
rasterImage(img, 1,6, 5,10)




smiles <- msdev.13C1@statData$MSIP$isotopologues_data[[1]]$compound_info$smiles
sdf.chemmine <- smiles2sdf(smiles)

sdf.ig <- get_sdf_igraph(sdf.chemmine)[[1]]
vis_sdf_igraph(sdf.ig,show_id = T)

plot(a$x,a$y)

canonicalNumbering(sdf.chemmine)
canonicalNumbering(sdf.chemmine)
sdfs <- sdf.chemmine
sdf<-sdfs[[1]]

canonicalNumbering_OB(obmol(sdf))

n <- 2^(1:20)

for (i in n) {

  sdfs <- smiles2sdf(rep(smiles,i))
  message_with_time(i)
  system.time(a <- canonicalNumbering(sdfs))%>%print()
  system.time(a <- sapply(1:i,function(x){
    canonicalNumbering_OB(obmol(sdfs[[1]]))
  }))%>%print()

}


i <- 17
smiles <- msdev.13C1@statData$MSIP$isotopologues_data[[i]]$compound_info$smiles
msdev.13C1@statData$MSIP$isotopologues_data[[i]]$compound_info$name
sdf.chemmine <- smiles2sdf(smiles)

sdf.ig <- get_sdf_igraph(sdf.chemmine)[[1]]
vis_sdf_igraph(sdf.ig,show_id = T)


canonicalNumbering()

iso.list <- msdev.13C1@statData$MSIP$isotopologues_data
iso.data <- iso.list$FT00366_Negative
msip.core <- iso.data$MSIP_result$M1$FS_3_13C


### assign set to isotopomer
{
  isotopomer.set.prob <- p_estimated/sum(p_estimated)
  isotopomer.set <- MSIPIsotopomerMap@solve$isotopomer.set
  isotopomer.prob <- mapply(x = isotopomer.set.prob ,
                            y = isotopomer.set,function(x,y){
                              make_vector(x/length(y),num2str(y,10))
                            },SIMPLIFY = F)
  isotopomer.prob <- unlist(isotopomer.prob)
  isotopomer.prob <- isotopomer.prob[order(names(isotopomer.prob))]
  isotopomer.prob <- unname(isotopomer.prob)
}


msip.core <- msdev.13C1@statData$MSIP$isotopologues_data[[17]]$MSIP_result$M1$FS_1_13C
msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob%>%str_digit(2)
a <- MSIPCore_correct_natural(msip.core,0.5)
a@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob

msip.core@solve$Atom_prob
a@solve$Atom_prob

MSIPCore_solve(MSIPCoreData)

# Mon Oct 14 14:48:32 2024 show sp merge------------------------------
{

  iso.list <- msdev.Astral@statData$MSIP$isotopologues_data
  iso.data <- iso.list$FT12669_Positive

  MSIPCoreData <- get_MSIPCoreData(sp.iso = iso.data$Spectra$M4$GLUCE1,
                   cfmd = iso.data$CFM_annotation,
                   iso_count = 4)
  a <-  MSIPCore_solve(MSIPCoreData)
  peak.df <- CFM_annotate_isotopologues(sp = iso.data$Spectra$M4$GLUCE1,
                             cfmd = iso.data$CFM_annotation,
                             iso_count = 4)
  sp.frag.data <- CFM_spectra_data_merge(peak.df,4)
  fg.map <- get_MSIPFragmentMap(sp.frag.data,
                                cfmd,
                                iso_count = iso_count)

  dist_matrix <- dist((x.ratio))

}

ggplot(frag.df)+
  geom_point(aes(x = log10(int_sum),
                 y = cos,
                 size = peaks_count))


plot.data%>%
  dplyr::group_by(atom,name)%>%
  dplyr::summarise(sum(value))


p <- ggplot(pi)+
  geom_boxplot(aes(y = name_m,x = solve.ratio))+
  geom_point(aes(y = name_m,x = solve.ratio))+
  labs(x = "Solve Ratio", y = "")+
  theme_bw()

open_plot_win(p,5,10)
export::graph2pdf(p,file = "a.pdf",width = 5,height = 10)

{

  process.info.as <- MSIP_solve_computation_evaluate(msdev.Astral, show_message = F)
  pi.as <- process.info.as%>%
    dplyr::mutate(solve.ratio = FSIS.count/isotopomer,
                  name_m = paste0(name,"-M+",iso_count))%>%
    dplyr::group_by(name_m)%>%
    dplyr::mutate(mean_ratio = mean(solve.ratio,na.rm =T))%>%
    dplyr::ungroup()%>%
    dplyr::arrange(mean_ratio)%>%
    dplyr::mutate(name_m = factor(name_m,level = unique(name_m)),
                  instru = "Astral")%>%
    dplyr::filter((mean_ratio >=0.3))

}

{
  process.info.qe <- MSIP_solve_computation_evaluate(msdev.13C1, show_message = F)
  pi.qe <- process.info.qe%>%
    dplyr::mutate(solve.ratio = FSIS.count/isotopomer,
                  name_m = paste0(name,"-M+",iso_count))%>%
    dplyr::group_by(name_m)%>%
    dplyr::mutate(mean_ratio = mean(solve.ratio,na.rm =T))%>%
    dplyr::ungroup()%>%
    dplyr::arrange(mean_ratio)%>%
    dplyr::mutate(name_m = factor(name_m,level = unique(name_m)),
                  instru = "QE")%>%
    dplyr::filter((mean_ratio >=0.3))
}
plot.data <- rbind(pi.as,pi.qe)

ggplot(plot.data)+
  geom_histogram(aes(x = solve.ratio, after_stat(density),fill = instru),
                 bins = 10,
                 position = "dodge")


# Thu Oct 17 16:14:18 2024 ------------------------------
{

  data.list <- list()
  for (i in 1:nrow(process.info)) {

    msip.core <- object@statData$MSIP$isotopologues_data[[process.info$feature_id[i]]]$MSIP_result[[str_isotope2_num(process.info$isotopomer[i])]][[process.info$samples[i]]]
    if (is.null(msip.core)) next
    plot.data <- msip.core@Spectra_data%>%
      dplyr::filter(merged)->data.list[[i]]
    ggplot(plot.data)+
      geom_point(aes(x = log10(int_sum),
                     y = icc,
                     size = peaks_count ))
  }
  data.list <- data.list[!sapply(data.list,is.null)]
  plot.data <- data.table::rbindlist(data.list, fill=TRUE)%>%
    dplyr::filter(peaks_count > 5,
                  peaks_count < 50)

  ggplot(plot.data)+
    geom_point(aes(x = log10(int_sum),
                   y = icc,
                  size = peaks_count),alpha = 0.1)+
    labs(y = "ICC(Intraclass correlation coefficient)\nof FG ratio")+
    theme_bw()->p
  open_plot_win(p,6,5)


}


# Wed Nov 13 10:48:54 2024 ------------------------------
kegg.ig <- FELLA::buildGraphFromKEGGREST('hsa')
vda <- vdata(kegg.ig)
eda <- edata(kegg.ig)

dis <- distances(kegg.ig,v = "R01224")
ig.sub <- igraph_filter_distance(kegg.ig,c("C01189"),1)

visIgraph(ig.sub)%>%
  open_visNet()

paths <- all_simple_paths(kegg.ig,
                          "C00031","C00051",
                          mode   = "all")


kegg.rc <- KEGGREST::keggList("reaction")

kegg.rc.split <- split(kegg.rc,
                       ceiling(seq_along(kegg.rc) / 10))

kegg.rc.data <- plyr::llply(kegg.rc.split,
                    .fun = function(x){
                      success <- FALSE
                      result <- NULL
                      while (!success) {
                        Sys.sleep(1)
                        try({
                          rc.da <- KEGGREST::keggGet(names(x))
                        success <- TRUE },
                        silent = TRUE)
                      }
                      return(rc.da)
                    },
                    .progress = "text")

kegg.rc.data <- unlist(kegg.rc.data,recursive = F)
names(kegg.rc.data) <- sapply(kegg.rc.data,`[[`,"ENTRY")

a <- lapply( kegg.rc.data, KEGG_reaction_parse)



a <-do.call(rbind,kegg.rc.data[1:3])



vars <- c("ENTRY","NAME","DEFINITION","EQUATION","ENZYME",
          "BRITE","DBLINKS","COMMENT","RCLASS","PATHWAY",
          "ORTHOLOGY","MODULE","REMARK","REFERENCE")
var <- "RCLASS"

sapply(kegg.rc.data,function(rc.data){
 rc.data[[var]]%>%class
  })%>%unlist%>%table()

sapply(kegg.rc.data,function(rc.data){
  rc.data[[var]]%>%length
})%>%unlist%>%table()



eqs <- sapply(kegg.rc.data,`[`,"EQUATION")

a <- lapply(eqs,KEGG_reaction_EQUATION_parse)

kegg.rc.data.p <-
 plyr::llply(kegg.rc.data,
             KEGG_reaction_parse,
             .progress = "text")


rc.ig <- MSdb::get_KEGG_Reaction_network()
eda <- edata(rc.ig)
vda <- vdata(rc.ig)

rc.ig <- igraph_filter_vertex(rc.ig,vda$name!="C00001")
eda <- edata(rc.ig)
vda <- vdata(rc.ig)

paths <- all_simple_paths(kegg.rig,
                          from = "C00031",
                          to = "C00051",
                          mode = "all",
                          cutoff = 9)

igraph_filter_path(rc.ig,
                   paths)%>%
  vis_igraph()%>%
  open_visNet()



kegg.raw.data <- MSdb:::get_KEGG_rawdata()

# Mon Nov 18 12:48:27 2024 ------------------------------
kegg.rig <- MSdb::get_KEGG_Reaction_network()
kegg.rig <- KEGG_Reaction_network_add_label(kegg.rig)


path.stat <- data.frame(
  cutoff = 8:12,
  path.count = NA
)

for (i in 1:nrow(path.stat)) {
  paths <- all_simple_paths(kegg.rig,
                            from = "C00031",
                            to = "C00051",
                            mode = "all",
                            cutoff = path.stat$cutoff[i])
  path.stat$path.count[i] <- length(paths)
}

p1 <- ggplot(path.stat)+
  geom_bar(aes(x = cutoff,y = path.count),stat = "identity")

kegg.rig <- KEGG_Reaction_network_merge_path(kegg.rig)
path.stat <- data.frame(
  cutoff = 8:12,
  path.count = NA
)
for (i in 1:nrow(path.stat)) {
  paths <- all_simple_paths(kegg.rig,
                            from = "C00031",
                            to = "C00051",
                            mode = "all",
                            cutoff = path.stat$cutoff[i])
  path.stat$path.count[i] <- length(paths)
}

p2 <- ggplot(path.stat)+
  geom_bar(aes(x = cutoff,y = path.count),stat = "identity")




kegg.rig <- KEGG_Reaction_network_merge_path(kegg.rig)
paths <- igraph::all_simple_paths(kegg.rig,
                          from = "C00031",
                          to = "C00051",
                          mode = "all",
                          cutoff = 10)

ig.path <- igraph_filter_path(kegg.rig,
                   paths)

vda <- vdata(ig.path)%>%
  dplyr::mutate(
    ff = MSCC:::chemform_formate(Formula)
  )



kegg.rig <- KEGG_Reaction_network_remove_nonformat_node(kegg.rig)
paths <- igraph::all_simple_paths(kegg.rig,
                                  from = "C00031",
                                  to = "C00051",
                                  mode = "all",
                                  cutoff = 10)

ig.path <- igraph_filter_path(kegg.rig,
                              paths)
vis_igraph(ig.path)%>%
  open_visNet()

eda <- edata(ig.path)

kegg.rawdata <- MSdb:::get_KEGG_rawdata()
kegg.rdata <- kegg.rawdata$Reaction_rawdata
R00149 <- kegg.rdata[["R00149"]]


# Mon Nov 18 15:45:32 2024 ------------------------------

kegg.rig <- MSdb::get_KEGG_Reaction_network()
kegg.rig <- KEGG_Reaction_network_add_label(kegg.rig)
kegg.rig <- KEGG_Reaction_network_remove_nonformat_node(kegg.rig)



ig <- kegg.rig
ig <- igraph_add_reverse_edges(ig)

paths <- igraph::all_simple_paths(kegg.rig,
                                  from = "C00031",
                                  to = "C00051",
                                  mode = "all",
                                  cutoff = 10)
path.ig <- igraph_filter_path(ig,paths)
vis_igraph(path.ig)%>%
  open_visNet()
epaths <- igraph_vpath_to_epath(ig,paths)

eda <- edata(ig)

epaths.dis <- plyr::laply(epaths,function(x){
  sum(eda$direction[x])
},.progress = "text")

epaths.len <- lengths(epaths)

plot.data <- data.frame(
  x = epaths.dis,
  y = epaths.len)

p <- ggplot2::ggplot(plot.data)+
  geom_jitter(aes(x,y))

open_plot_win(p)


kegg.ig <- MSdb::get_KEGG_Reaction_network()

eda <- edata(ig.krn)%>%
  dplyr::filter(REACTION_id == "R08575")

# Fri Nov 22 14:18:26 2024 ------------------------------
name = "aaa"
path = "c:/aaa.d"

cli::cli_inform(c(x = "Reading {.path {path}}"))
# Load the cli package
library(cli)

# Define the path and the text to display
path <- "C:/Users/91879/OneDrive/Code/R/Package/MSdev/Script/Test_Script.R"
text <- "Click here to read the file"

# Use cli::cli_inform to display the link with custom text
cli::cli_inform("Reading {.href  [AAA](file://{path}:10)}")
cli::cli_inform("Reading {.file  file://{path}:10}")


all.m.reaction <- lapply(kegg.mdata,
                        function(mdata){
                          mdata$REACTION
                        })%>%
  unname()%>%unlist()

grepl(" <- ",all.m.reaction)%>%table
all.m.reaction[grepl(" <- ",all.m.reaction)]



a <- data.frame(
  n = names(all.m.reaction),
  v = all.m.reaction
)

edge_paths <- lapply(this.path, function(path) {
  # Convert path to numeric vector
  path <- as.numeric(path)

  # Find edge IDs for consecutive vertex pairs
  edge_ids <- sapply(seq_along(path)[-length(path)], function(i) {
    get_edge_ids(this.ig, c(path[i], path[i + 1]))
  })
  edge_ids
})


# Thu Nov 28 15:39:58 2024 ------------------------------
kegg.mdata.df <- kegg.mdata.df%>%
  dplyr::mutate(str_syn =
                  paste0(REACTION_id,"_",
                         from,"_",
                         to),
                x = case_when(
                  str_syn%in% eda$str_syn~1,
                  str_syn%in% eda$str_syn_rev~-1

                )      )


ig.krn.dired<- igraph_filter_edge(ig.krn,
                                  which(!is.na(edata(ig.krn)$direction))

                                  )

paths <- igraph::all_simple_paths(ig.krn,
                                  from = "C00031",
                                  to = "C00051",
                                  cutoff = 8,
                                  mode = "all")



# Mon Dec  2 13:05:16 2024 ------------------------------
ig.krn <- get_KEGG_Reaction_network()
eda <- edata(ig.krn)
table(eda$direction)
ig.krn.hsa <- KEGG_Reaction_network_filter_by_emzyme(
  ig.krn,"hsa"
)
eda <- edata(ig.krn.hsa)
vda <- vdata(ig.krn.hsa)
table(eda$direction)


paths <- igraph::all_simple_paths(ig.krn.hsa,
                          "C00031","C00051",
                          mode   = "all",cutoff = 10)


# Mon Dec  2 15:39:03 2024 MFNA figures------------------------------
### prob shared by isotopomers
{
  msip.core <- iso.data$MSIP_result$M2$U
  msip.core <- MSIPCore_solve(msip.core,
                              int_thresh = 10^4)
  solve.data <- data.frame(
    is = names(msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set),
    Raw =msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob,
    Natural = lengths(msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set)
  )%>%
    dplyr::mutate(Natural =Natural/sum(Natural)*natural.ratio,
                  Adjusted = Raw-Natural,
                  Adjusted = case_when(Adjusted<0~0,
                                       T~Adjusted),
                  Adjusted = Adjusted/sum(Adjusted))%>%
    tidyr::pivot_longer(2:4)

  is.data <- solve.data%>%
    dplyr::filter(name == "Adjusted" )%>%
    dplyr::filter(value > 0.01)
  isotopomers.idx <- msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set[is.data$is]
  plot.data <- lapply(isotopomers.idx,function(idx){
    isotopomer.name <- sapply(msip.core@solve$MSIPIsotopomerMap@isotopomer.defination[idx],
                              function(x){
                                paste0(x,collapse = ",")
                              })
    data.frame(
      isotopomer.name = isotopomer.name
    )
  })%>%
    data.table::rbindlist(idcol = "is")%>%
    dplyr::group_by(is)%>%
    dplyr::mutate(prob = is.data$value[match(is,is.data$is)] ,
                  w = case_when(n()==1~0.8,
                                T~1))%>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(prob))%>%
    dplyr::mutate(isotopomer.name=
                    factor(isotopomer.name,level = isotopomer.name))


  p <- ggplot(plot.data)+
    geom_bar(aes(x = isotopomer.name,
                 y = prob,
                 fill = is,
                 width = w),
             stat = "identity")+
    ggsci::scale_fill_npg()+
    labs(title = "Glutathione M+2",y = "Probability",x = NULL,fill = "FSIS")+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = -45,hjust = 0),
          legend.position.inside = c(0.7,0.8))
  p
  open_plot_win(p,5,3)

  atom.prob <- c("C_8"=0.6,"C_9"=0.6)
  get_cfm_data_sdf_igraph(iso.data$CFM_annotation)%>%
    sdf_igraph_add_background_color(atom.prob,color.ramp = cf)%>%
    sdf_igraph_add_border_color(atom.prob,
                                color.ramp = cf)%>%
    vis_sdf_igraph()%>%
    visNodes(shadow= T)%>%
    open_visNet()


}

###
{
  df <- expand.grid(
    x = 1:20,
    y = 1:20
  )%>%
    dplyr::filter(y <= x) %>%
    dplyr::mutate(n = choose(x,y),
                  isotopologues = paste0("M+",y))


  ggplot(df)+
    geom_point(aes(x = x,y = n ,
                   fill = y, size = y),
               alpha = 0.5,
               shape  = 21)+
    scale_y_log10()+
    scico::scale_fill_scico(direction = -1)+
    labs(x = "Number of C", y = "Number of Isotopomers" )+
    theme_bw() +  # Ensure same name for size
    guides(
      fill = guide_legend(title = "Isotopologues"),
      size = guide_legend(title = "Isotopologues")
    )->p
  open_plot_win(p,4,3)

  plot.data<-df%>%
    dplyr::group_by(x)%>%
    dplyr::mutate(total.isotopomers = sum(n),
                  total.isotopologues =x+1
                  )%>%
    dplyr::distinct(x,total.isotopologues,.keep_all = T)%>%
    pivot_longer(total.isotopomers:total.isotopologues)

  ggplot(plot.data)+
    geom_bar(aes(x = x,y=value,fill = name),
             stat = "identity",position = "dodge")+
    scale_y_log10()+
    ggsci::scale_fill_aaas()+
    coord_flip()+
    labs(x = "Number of C", y = "Data size",
         fill = "Data Type")+
    theme_bw()->p
  open_plot_win(p,4,3)
}

# Tue Dec  3 16:23:33 2024 ------------------------------
{
  ig.krn <- get_KEGG_Reaction_network()
  ig.krn.hsa <- KEGG_Reaction_network_filter_by_emzyme(
    ig.krn,"hsa"
  )
  paths <- igraph::all_simple_paths(ig.krn.hsa,
                                    "C00031","C00051",
                                    mode   = "all",cutoff = 8)
  ep <- igraph_vpath_to_epath(ig.krn.hsa,paths)
  ig <- igraph_filter_path(ig.krn.hsa,paths)
  igraph_add_vcolor(ig,v =  c("C00031","C00051"),color = "#F36482")%>%
    igraph_add_vfill(v =  c("C00031","C00051"),color = "#F36482")%>%
    #igraph_add_vfill(v =  c("C00024"),color = "#F36482")%>%
    igraph_add_ecolor(e = which(edata(ig )$direction==1),color = "#498FED")%>%
    igraph_add_earrow(e = which(edata(ig )$direction==1),"to")%>%
    vis_igraph()%>%
    visNodes(shadow = T)%>%
    visEdges(width = 5,color= "#DDDDDD")%>%
    open_visNet()

}

