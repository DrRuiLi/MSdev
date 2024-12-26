get_Molecule_igraph_from_sdf <- function(sdf = sdfsample) {
  .f <- function(sdf) {


    atom.block <- atomblock(sdf) %>%
      as.data.frame() %>%
      rownames_to_column("id") %>%
      dplyr::mutate(name = id,
                    bonds(sdf),
                    element = atom,
                    x = C1,
                    y = C2,
                    .before = id)

    bond.block  <- bondblock(sdf) %>%
      as.data.frame() %>%
      dplyr::mutate(from = atom.block$id[C1],
                    to = atom.block$id[C2],
                    bond_type = C3,
                    .before = C1)

    sdf.igraph <- igraph::graph_from_data_frame(bond.block, vertices = atom.block)

    isotopomer.df <-  make_vector(atom.block$element,atom.block$id)%>%
      as.list()%>%
      as.data.frame()%>%
      dplyr::mutate(isotopomer = "natural",
                    isotopologue = "M0",
                    label = "natural",
                    abundance = 1,
                    .before = everything()
                    )
    rownames(isotopomer.df) <- isotopomer.df$isotopomer
    Molecule_igraph <- new("Molecule_igraph",
                           sdf = sdf, igraph = sdf.igraph,isotopomer = isotopomer.df)
    Molecule_igraph@molecule_info$smiles <- unname(as.character(get_sdf_smiles(sdf)))
    return(Molecule_igraph)
  }

  if (class(sdf) == "SDF")
    Molecule_igraphs <- .f(sdf)


  if (class(sdf) == "SDFset") {
    Molecule_igraphs <- list()
    sdf.valid <- validSDF(sdf)
    for (i in 1:length(sdf)) {
      Molecule_igraphs[[i]] <- .f(sdf[[i]])
    }
    names(Molecule_igraphs) <- cid(sdf)
  }

  return(Molecule_igraphs)

}


get_Molecule_igraph_from_smiles <- function(smiles) {
  sdf <- get_smiles_sdf(smiles)
  if (length(sdf)==1) sdf <- sdf[[1]]
  get_Molecule_igraph_from_sdf(sdf)
}



setMethod(
  "show",
  "Molecule_igraph",
  definition = function(object) {
    print(paste0("Molecule_igraph: ", unname(MF(object@sdf, addH =  T))))
  }
)

setGeneric(
  "vdata",
  def = function(object) {
    igraph::as_data_frame(object, "vertices")
  }
)
setMethod(
  "vdata",
  "Molecule_igraph",
  definition = function(object) {
    vdata(object@igraph)
  }
)



setGeneric(
  "vdata<-",
  def = function(object, value) {
    igraph::vertex.attributes(object) <- as.list(value)
    object
  }
)

setMethod(
  "vdata<-",
  "Molecule_igraph",
  definition = function(object, value) {
    vdata(object@igraph) <- value
    object
  }
)


setGeneric(
  "edata",
  def = function(object) {
    igraph::as_data_frame(object, "edges")
  }
)
setMethod(
  "edata",
  "Molecule_igraph",
  definition = function(object) {
    edata(object@igraph)
  }
)



setGeneric(
  "edata<-",
  def = function(object, value) {
    value <- value[, !grepl("from|to", colnames(value))]
    igraph::edge.attributes(object) <- as.list(value)
    object
  }
)
setMethod(
  "edata<-",
  "Molecule_igraph",
  definition = function(object, value) {
    edata(object@igraph) <- value
    object
  }
)





setMethod(
  f = "plot",
  signature = "Molecule_igraph",
  definition = function(object,x,y) {
    plot(object@sdf)
    invisible()
  }
)


setGeneric("atom",
           def = function(object,
                          element = element_table$element){
    rownames(atomblock(object))[bonds(object)$atom%in%element ]
  }
)

setMethod("atom",
          "Molecule_igraph",
          definition = function(object,
                                element = element_table$element){
            vdata(object)%>%
              dplyr::filter(element %in% !!element)%>%
              dplyr::pull(name)
          })

setGeneric("element",
           def = function(object,...){
             bonds(object)$atom
           }
)

setMethod("element",
          "Molecule_igraph",
          definition = function(object,...){
            vdata(object)%>%
              dplyr::pull(element)
          })


setMethod("formula","Molecule_igraph",
          definition = function(x,...){
            unname(MF(x@sdf,addH = T))
          })

add_Molecule_igraph_isotopomer <- function(
    Molecule_igraph , isotopomer = NULL,iso_vec = NULL,abundance=NA){


  iso_label <- split(names(iso_vec),iso_vec)
  iso_label <- sapply(seq_along(iso_label),function(x){
    paste0("(",paste0(iso_label[[x]],collapse = ","),")",names(iso_label)[x])
  })%>%paste0(collapse = ";")
  ele_vec <-make_vector(vdata(Molecule_igraph)$element,
                        atom(Molecule_igraph))
  ele_vec[names(iso_vec)] <- iso_vec
  isotopologue <- sum(is.isotope(ele_vec),na.rm = T)
  isotopologue <- paste0("M",isotopologue)
  isotopomer.df <- Molecule_igraph@isotopomer
  if (is.null(isotopomer)){
    i <- 1
    while(paste0(isotopologue,"_",i) %in% isotopomer.df$isotopomer){
      i <- i+1
    }
    isotopomer <- paste0(isotopologue,"_",i)
  }
  to.add <- data.frame(isotopomer = isotopomer,
                       isotopologue = isotopologue,
                       label = iso_label ,
                       abundance = abundance
                       )
  rownames(to.add) <- to.add$isotopomer
  isotopomer.df <- bind_rows(isotopomer.df,to.add)
  isotopomer.df[isotopomer,atom(Molecule_igraph)] <- ele_vec

  isotopomer.df -> Molecule_igraph@isotopomer
  return(Molecule_igraph)


}


get_Molecule_igraph_MS1 <- function(Molecule_igraph,polarity=1,adduct = NULL){

  if (is.null(adduct)  ) {
    adduct <- ifelse(polarity==0,"[M-H]-","[M+H]+")
  }


  ### ms spectra
  {

    mz.m0 <- MSCC::chemform_adduct(formula(Molecule_igraph),adduct = adduct)
    isotopomers.eles <- as.matrix(Molecule_igraph@isotopomer[,atom(Molecule_igraph)])
    ele.iso.diff <- make_vector(element_table$Mass_Dif,element_table$symbol)
    isotopomers.mz.diff <- ele.iso.diff[isotopomers.eles]
    dim(isotopomers.mz.diff) <- dim(isotopomers.eles)
    dimnames(isotopomers.mz.diff) <- dimnames(isotopomers.eles)
    mz.diff <- apply(isotopomers.mz.diff,1,sum)

    ms1.data <- Molecule_igraph@isotopomer%>%
      dplyr::mutate(mz = mz.m0+mz.diff)%>%
      dplyr::group_by(isotopologue)%>%
      dplyr::mutate(intensity = sum(abundance) )%>%
      dplyr::ungroup()%>%
      dplyr::select(any_of(c("isotopomer","isotopologue","label","abundance","mz","intensity")))

  nor.to <- ms1.data$intensity[which.min(ms1.data$mz)]
  ms1.data$intensity <-  ms1.data$intensity/nor.to*100
  }

  return(ms1.data)



}

get_Molecule_igraph_MS2 <- function(Molecule_igraph,cfmd){


  isotopomers <- Molecule_igraph@isotopomer
  FG.map <- cfmd@fragment_group_map

  #### prob
  {
    FG.map[FG.map<0.5] <- 0
    FG.map[FG.map>=0.5] <- 1
  }

  isotopomers.eles <- as.matrix(Molecule_igraph@isotopomer[,atom(Molecule_igraph)])
  ele.iso.diff <- make_vector(element_table$Mass_Dif,element_table$symbol)
  isotopomers.mz.diff <- ele.iso.diff[isotopomers.eles]
  dim(isotopomers.mz.diff) <- dim(isotopomers.eles)
  dimnames(isotopomers.mz.diff) <- dimnames(isotopomers.eles)

  isotopomers.mz.diff <- isotopomers.mz.diff[,colnames(FG.map)]
  fg.isotopomers.mz.diff <- isotopomers.mz.diff %*% t(FG.map)


  ms2.data <- cbind( cfmd@fragment_group,t(fg.isotopomers.mz.diff))%>%
    dplyr::filter(fragment_mz > 0)%>%
    tidyr::pivot_longer(rownames(isotopomers.eles),names_to = "isotopomer",values_to = "mzdiff")%>%
    dplyr::mutate( isotopomers[match(isotopomer,isotopomers$isotopomer),
                               c("isotopologue","label","abundance")],
                   mz = fragment_mz +mzdiff)%>%
    dplyr::group_by(mz)%>%
    dplyr::mutate(intensity = sum(abundance))



  return(ms2.data)


}


Molecule_igraph_get_C_order <- function(Molecule_igraph){

  dism <- igraph::distances(Molecule_igraph@igraph)
  ele <- element(Molecule_igraph)
  ele.dis.R <- apply(dism,1,function(x){
    sum(x[ele!="C"])
  })
  ele.dis.C <- apply(dism,1,function(x){
    sum(x[ele=="C"])
  })
  ele.dis.R[ele=="C"]
  ele.dis.C[ele=="C"]

}


vis_Molecule_igraph <- function(Molecule_igraph,show_id = F){

  Molecule_igraph.formated <- Molecule_igraph%>%
    Molecule_igraph_vis_format()%>%
    sdf_igraph_show_id(show_id)



  vis_igraph(Molecule_igraph.formated) %>%
    visPhysics(enabled = F) %>%
    visOptions(width = "100%", height = "100%")

}


Molecule_igraph_vis_format <- function(Molecule_igraph){

  vdata(Molecule_igraph) <-  vdata(Molecule_igraph)%>%
   # dplyr::group_by(atom)%>%
    dplyr::mutate(
      label = paste0(" ",atom," "),
      x =  ( x-mean(x))*100,
      y = ( y-mean(y)) *100,
      font.size = 30,
      borderWidth = 20,
      font.vadjust = 5,
      font.align = "center",
      color.border = "#AAAAAA",
      color.background = "#FFFFFF",
      borderWidth = 5,
      shape = "circle",
      physics = F)%>%
    dplyr::ungroup()


  edata(Molecule_igraph)  <- edata(Molecule_igraph)%>%
    dplyr::mutate(smooth = FALSE,
                  width = 10*bond_type)

  return(Molecule_igraph)

}


