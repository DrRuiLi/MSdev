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



Metabolic_flux_network_get_atom_transfer <- function(Metabolic_flux_network){

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



Metabolic_flux_set_tracer <-
  function(Metabolic_flux_network){




}



load_MFN <- function(path = "C:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/Metabolic_flux_network/"){

  mfn.files <- dir(path,full.names = T)%>%
    file.info()%>%
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



