get_Molecule_igraph_from_sdf <- function(sdf = sdfsample) {
  .f <- function(sdf) {


    atom.block <- atomblock(sdf) %>%
      as.data.frame() %>%
      rownames_to_column("id") %>%
      dplyr::mutate(bonds(sdf),
                    element = atom,
                    x = C1,
                    y = C2,
                    .after = id)

    bond.block  <- bondblock(sdf) %>%
      as.data.frame() %>%
      dplyr::mutate(from = atom.block$id[C1],
                    to = atom.block$id[C2],
                    .before = C1)

    sdf.igraph <- igraph::graph_from_data_frame(bond.block, vertices = atom.block)

    isotopomer.df <-  make_vector(atom.block$element,atom.block$id)%>%
      as.list()%>%
      as.data.frame()%>%
      dplyr::mutate(isotopomer = "natural",
                    isotopologue = "M0",
                    abundance = 1,
                    .before = everything()
                    )
    rownames(isotopomer.df) <- isotopomer.df$isotopomer
    new("Molecule_igraph", sdf = sdf, igraph = sdf.igraph,isotopomer = isotopomer.df)

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
    rownames(atomblock(object))[bonds(sdf)$atom%in%element ]
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



add_Molecule_igraph_isotopomer <- function(
    Molecule_igraph , isotopomer = NULL,iso_vec = NULL,abundance=NA){

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
                       abundance = abundance
                       )
  rownames(to.add) <- to.add$isotopomer
  isotopomer.df <- bind_rows(isotopomer.df,to.add)
  isotopomer.df[isotopomer,atom(Molecule_igraph)] <- ele_vec

  isotopomer.df -> Molecule_igraph@isotopomer
  return(Molecule_igraph)


}

