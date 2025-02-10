get_RXNMapper <- function(){

  rxnmp <- reticulate::import("rxnmapper")
  rxn_mapper = rxnmp$RXNMapper()
  return( rxn_mapper$get_attention_guided_atom_maps)

}


RXNMapper_map  <- function(RXNMapper = get_RXNMapper(),
                           from.smiles ,
                           to.smiles
                           ){



  from.string <- paste0(from.smiles,collapse = ".")
  to.string <- paste0(to.smiles,collapse = ".")

  req <- paste0(from.string, ">>",to.string)
  rxns <- c(req,"C>>C")


  rxn.result <- RXNMapper(rxns)[[1]]
  rxn.result$from <- from.smiles
  rxn.result$to <- to.smiles

  return(rxn.result)

}


RXNMapper_mapped_rxn_parse <- function(mapped_rxn){


  rxn.result.split <- str_split(mapped_rxn,">>")[[1]]%>%
    sapply(function(x){
      str_split(x,"\\.")
    })
  names(rxn.result.split ) <- c("from","to")
  rxn.canonical.id <- lapply(rxn.result.split,function(x){
    lapply(x, function(y){
      y.sdf <- get_smiles_sdf(y,canonicalize = F)
      y.mig <- get_Molecule_igraph_from_sdf(y.sdf[[1]])
      y.rxn.id <- str_extract_all(y, "(?<=:)(\\d+)")[[1]]%>%
        as.numeric()%>%
        `names<-`(atom(y.mig))
      y.canon.id <- canonicalNumbering(y.sdf)[[1]]
      y.rxn.id[order(y.canon.id)]%>%unname()
    })
  })

  return(rxn.canonical.id)


}


