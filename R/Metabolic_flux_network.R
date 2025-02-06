setClass("Metabolic_flux_network",
         slots = list(
           metabolic_network = "ANY",
           Molecule_igraphs = "list" ))


setMethod("show",signature = "Metabolic_flux_network",
          definition =  function(object){
            message("Metabolic_flux_network with ", nrow(vdata(object)), " nodes")
          })



setMethod(
  "vdata",
  signature = "Metabolic_flux_network",
  definition = function(object) {
    vdata(object@metabolic_network)
  }
)

setMethod(
  "vdata<-",
  "Metabolic_flux_network",
  definition = function(object, value) {
    vdata(object@metabolic_network) <- value
    object
  }
)


setMethod(
  "edata",
  "Metabolic_flux_network",
  definition = function(object) {
    edata(object@metabolic_network)
  }
)

setMethod(
  "edata<-",
  "Metabolic_flux_network",
  definition = function(object, value) {
    edata(object@metabolic_network) <- value
    object
  }
)




setClass("Molecule_atom_transfer",
         slots = list(
           "transfer_def" = "data.frame",
           "transfer_matrix" = "matrix"
         ))


setMethod("show",
          "Molecule_atom_transfer",definition = function(object){
            message(paste0(nrow(object@transfer_matrix)," Molecule atom transfer "))
          })

setMethod("as.character",
          "Molecule_atom_transfer",
          definition =  function(x,...){
            paste0(paste0(nrow(x@transfer_matrix)," Molecule atom transfer "))
          })

setMethod("length",
          "Molecule_atom_transfer",
          definition =  function(x){
            nrow(x@transfer_matrix)
          })



Metabolic_flux_network_get_Molecule_atom_transfer <- function(Metabolic_flux_network){

  mfn.e <- edata(Metabolic_flux_network)
  mfn.transfer <- plyr::mlply(mfn.e,.fun = function(from,to,...){
    get_Molecule_atom_transfer_by_atom_map(
      V(Metabolic_flux_network@metabolic_network)[[from]]$Molecule_igraph,
      V(Metabolic_flux_network@metabolic_network)[[to]]$Molecule_igraph,
      target_ele = "C"
    )
  },
  .progress = "text")
  names(mfn.transfer) <- mfn.e$name
  attributes(mfn.transfer)$split_labels <- NULL
  attributes(mfn.transfer)$split_type <- NULL
  Metabolic_flux_network@metabolic_network <-
    igraph::set_edge_attr(Metabolic_flux_network@metabolic_network,
                        name = "atom_transfer",value =mfn.transfer
                        )
  return(Metabolic_flux_network)

}

Metabolic_flux_network_get_Reaction_atom_transfer <- function(mfn){

  mfn.v <- vdata(mfn)
  mfn.e <- edata(mfn)
  mfn.v.r <- mfn.v %>%
    dplyr::filter(node.type == "Reaction")

  mfn.transfer <- plyr::l_ply( mfn.v.r$id ,.fun = function(rid,...){

    from <- mfn.e%>%
      dplyr::filter(to == rid)%>%
      dplyr::pull(from)
    to <- mfn.e%>%
      dplyr::filter(from == rid)%>%
      dplyr::pull(to)

    reaction.data <- mfn.v.r %>%
      dplyr::filter(id == rid)

    equation.coef <- str_split(reaction.data$equation,"\\+|<=>",simplify = T)%>%
      sapply(function(x){
        x.coef <- str_extract(x,"^[:digit:]")
        x.coef <- ifelse(is.na(x.coef),1,x.coef)
        x.coef <- as.numeric(x.coef)
        x.id <- str_extract(x,"C[:digit:]*")
        make_vector(x.coef,x.id)
      },USE.NAMES = F)

    x <- rep(from,times = equation.coef[from])
    mol.ig.from <- V(mfn@metabolic_network)[x]$Molecule_igraph
    names(mol.ig.from) <- paste0(x,"_",sapply(equation.coef[from], function(x) seq(1, x))%>%unlist())


    x <- rep(to,times = equation.coef[to])
    mol.ig.to <- V(mfn@metabolic_network)[to]$Molecule_igraph
    names(mol.ig.to) <- paste0(x,"_",sapply(equation.coef[to], function(x) seq(1, x))%>%unlist())

    get_Reaction_atom_transfer_by_atom_map(
      mol.ig.from,
      mol.ig.to,
      target_ele = "C"
    )
  },
  .progress = "text")
  names(mfn.transfer) <- mfn.e$name
  attributes(mfn.transfer)$split_labels <- NULL
  attributes(mfn.transfer)$split_type <- NULL
  Metabolic_flux_network@metabolic_network <-
    igraph::set_edge_attr(Metabolic_flux_network@metabolic_network,
                          name = "atom_transfer",value =mfn.transfer
    )
  return(Metabolic_flux_network)

}



load_MFN <- function(path = "C:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/Metabolic_flux_network/",
                     name = "."){

  mfn.files <- dir(path,full.names = T)%>%
    file.info()%>%
    dplyr::mutate(basename = basename(rownames(.)))%>%
    dplyr::filter(grepl(name,basename))%>%
    dplyr::slice_max(mtime)
  readRDS(rownames(mfn.files))
}



vis_Metabolic_flux_network <- function(mfn){


  ### vis formate
  {

    vda <- vdata(mfn)%>%
      dplyr::mutate(#shape = "circle",
                    color.background = "white")
    eda <- edata(mfn)%>%
      dplyr::mutate(color.color = "rgba(84,126,158,0.5)",
                    color.highlight = "rgba(84,126,158,1)",
                    width = 8,
                    selectionWidth  = 12,
                    #arrows.to = T,
                    arrows.middle = T,
                    smooth = T)

  }

  visNetwork(nodes = vda,edges = eda)

}


Metabolic_flux_network_reverse <- function(mfn,edge.id){

  mfn@metabolic_network <- igraph::reverse_edges(mfn@metabolic_network,edge.id)
  mat <- E(mfn@metabolic_network)[[edge.id]]$atom_transfer%>%
    Molecule_atom_transfer_reverse()
  mat -> E(mfn@metabolic_network)[[edge.id]]$atom_transfer
  return(mfn)
}

Metabolic_flux_network_select_compound <- function(mfn,vid){

  ### find reaction between these cp
  eda <- edata(mfn)%>%
    dplyr::filter(from%in%vid|to %in% vid )
  reaction.selected <- setdiff(intersect(eda$from,eda$to),vid)

  ### select all cp in these reaction
  eda <- edata(mfn)%>%
    dplyr::filter(from%in%reaction.selected|to %in% reaction.selected )

  mfn@metabolic_network <- igraph_filter_edge(
    mfn@metabolic_network, which( E(mfn@metabolic_network)$name %in% eda$name)
  )

  return(mfn)

}

Metabolic_flux_network_get_compound_data_from_cid <- function(mfn){


  cp.node.data <- vdata(mfn)%>%
    dplyr::mutate(C.count = get_formula_ele_count(Formula ,"C"))%>%
    dplyr::filter(node.type=="Compound",
                  C.count > 0   )
  cp.cid <- webchem::get_cid(cp.node.data$PubChem,
                                            from = "sid",domain = "substance")
  pubchem.retrive <- webchem::pc_prop(cp.cid$cid)


  V(mfn@metabolic_network)$smiles <- NA
  V(mfn@metabolic_network)[cp.node.data$name]$smiles <- pubchem.retrive$IsomericSMILES



  V(mfn@metabolic_network)$Molecule_igraph  <- NA
  V(mfn@metabolic_network)[cp.node.data$name]$Molecule_igraph <- get_Molecule_igraph_from_smiles(
    pubchem.retrive$IsomericSMILES ,id = cp.node.data$name )



  mfn@metabolic_network <- igraph_filter_vertex(mfn@metabolic_network,
                       !is.na(V(mfn@metabolic_network)$smiles)|V(mfn@metabolic_network)$node.type=="Reaction"
                       )

  return(mfn)

}

Metabolic_flux_network_clean_reactions <- function(mfn){

  vda <- vdata(mfn)
  eda <- edata(mfn)

  x <- (vda$name %in% eda$from & vda$name %in% eda$to)|vda$node.type=="Compound"
  mfn@metabolic_network <- igraph_filter_vertex(mfn@metabolic_network,x)

  return(mfn)
}

Metabolic_flux_network_update_from_visGetEdges <- function(mfn,visGetEdges){


  ### edges assigment
  {
    eda <- edata(mfn)
    edge.delete <- setdiff(eda$id,names(visGetEdges))
    edge.add <- setdiff(names(visGetEdges),eda$id)

  }


  ### delete
  {

    mfn@metabolic_network <-
      igraph::delete_edges(mfn@metabolic_network,edge.delete)


  }


  ### add
  {
    edge.add.data <- visGetEdges[edge.add][[1]]
    #names(edge.add.data) <- paste0("AR",str_time(),num2str(1:length(edge.add.data)))
    #edges.to.add <- lapply(edge.add.data, function(x) c(x$from,x$to))%>%do.call(c,.)
    edges.to.add <- c(edge.add.data$from,edge.add.data$to)
    attr.to.add <- list(id = edge.add.data$id ,
                        name = edge.add.data$id ,
                        atom_transfer = list(
                          get_Molecule_atom_transfer_by_atom_map(
                            mol.ig.from = V(mfn@metabolic_network)[[edge.add.data$from]]$Molecule_igraph,
                            mol.ig.to = V(mfn@metabolic_network)[[edge.add.data$to]]$Molecule_igraph,
                            target_ele = "C"
                          )
                        )
                        )
    mfn@metabolic_network <- igraph::add_edges(mfn@metabolic_network ,
                           edges = edges.to.add,
                           attr  = attr.to.add
                           )
  }




  return(mfn)



}

Metabolic_flux_network_simplify <- function(mfn){

  mfn@metabolic_network <- igraph::simplify(mfn@metabolic_network ,remove.loops = F,edge.attr.comb = "first")
  mfn

}


Metabolic_flux_network_set_tracer <-
  function(Metabolic_flux_network,
           vid = 1,
           Molecule_igraph ){

    V(Metabolic_flux_network@metabolic_network)[[vid]]$Molecule_igraph <-Molecule_igraph


    return(Metabolic_flux_network)

  }


Metabolic_flux_atom_transfer <- function(mat,
                                         mol.ig.from,
                                         mol.ig.to
                                         ){


}

