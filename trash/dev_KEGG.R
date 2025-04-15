
KEGG_parse_REACTION_RCLASS <- function(RCLASS){



  stringr::str_extract_all(pattern = "[R]{0,1}C[0-9]{5}",RCLASS)%>%
    lapply(function(x){

      data.frame(RCLASS_id = x[1],
                 from = x[2],
                 to = x[3])
    })%>%
    data.table::rbindlist()

}

KEGG_reaction_to_network_RCLASS <- function(kegg.rdata){



  ### filter_org
  {
    if(F){

      kegg.link <- KEGGREST::keggLink("enzyme",filter_org)%>%
        FELLA:::sanitise("enzyme",filter_org)
      is.org <- sapply(kegg.rdata,function(reaction.data){
        any(reaction.data[["ENZYME"]] %in% kegg.link)
      })
      kegg.rdata.filtered <- kegg.rdata[is.org]
    }else{
      kegg.rdata.filtered <- kegg.rdata
    }

  }

  ### filter empty rclass
  {
    rclass.empty <- sapply(kegg.rdata.filtered,function(reaction.data){
      nrow(reaction.data$RCLASS_data)==0
    })
    kegg.rdata.filtered <- kegg.rdata.filtered[!rclass.empty]

  }


  ### edge
  {
    rcn.edges <- plyr::llply(
      kegg.rdata.filtered,
      function(reaction.data){

        if (!"RCLASS_data" %in% names(reaction.data)) {
          invisible()
        }
        enz <- reaction.data[["ENZYME"]]
        rda <- reaction.data[["RCLASS_data"]]%>%
          dplyr::mutate(REACTION_id = reaction.data$ENTRY,
                        REACTION_name = reaction.data$NAME,
                        equation = reaction.data$EQUATION)
        enzn <- rep(enz,nrow(rda))
        rdan <- rda[rep(1:nrow(rda),each = length(enz)),]
        rcda <- data.frame(rdan,
                           ENZYME = enzn)
        return(rcda)

      },.progress = "text"
    )
    rcn.edges.df <-do.call(bind_rows,rcn.edges)%>%
      dplyr::select(from,to,everything())%>%
      dplyr::mutate(name = paste0("UR",num2str(1:n())))

    if (F) {
      rcn.edges.df <- rcn.edges.df%>%
        dplyr::filter(ENZYME%in% kegg.link)
    }


  }

  ### node
  {

    kegg.cp <- MSdb:::get_KEGG_compound_df()%>%
      as.data.frame()%>%
      dplyr::filter(KEGG_id %in% rcn.edges.df$from|
                      KEGG_id %in% rcn.edges.df$to)

    rcn.nodes.df <- kegg.cp%>%
      dplyr::mutate(id= KEGG_id,
                    name = KEGG_id,
                    label = Name)



  }


  ### ig
  {
    kegg.reaction.ig <- igraph::graph_from_data_frame(
      rcn.edges.df,
      vertices = rcn.nodes.df)

  }

  return(kegg.reaction.ig)

}
