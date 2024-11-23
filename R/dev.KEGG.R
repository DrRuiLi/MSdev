get_KEGG_Reaction_network <- function(kegg.rdata ){


  kegg.rawdata <- MSdb:::get_KEGG_rawdata()
  kegg.rdata <-  kegg.rawdata$Reaction_rawdata
  kegg.rdata <- lapply(kegg.rdata,KEGG_parse_REACTION)
  ig.krn <- KEGG_reaction_to_network(kegg.rdata,"hsa")

  kegg.mdata <- kegg.rawdata$Module_rawdata
  ig.krn <- KEGG_reaction_get_direction_from_module(ig.krn,kegg.mdata)

  kegg.mdata <- kegg.rawdata$Module_rawdata




}
#1' @import igraph
KEGG_parse_REACTION <- function(reaction.data,
                                parse_by = "RCLASS" ){

  vars <- c("ENTRY",
            "NAME",
            "DEFINITION","EQUATION",
            "ENZYME","RCLASS",
            "PATHWAY","MODULE"
            # "BRITE","DBLINKS","COMMENT","ORTHOLOGY","REMARK","REFERENCE"
  )



  ### edit
  {
    #reaction.data <- reaction.data[vars]
    #reaction.data[sapply(reaction.data,is.null)]<- NA
    # names(reaction.data) <- vars


    reaction.data[["ENTRY"]] <- reaction.data[["ENTRY"]]%>%unname()
    reaction.data[["NAME"]]<- reaction.data[["NAME"]][1]


  }

  ### RCLASS
  {

    if(parse_by=="RCLASS"){

      #reaction.data <- kegg.rdata[[1]]
      RCLASS_data <- KEGG_parse_REACTION_RCLASS(reaction.data[["RCLASS"]])
      reaction.data$RCLASS_data <- RCLASS_data
    }


  }



  return(reaction.data)

}

KEGG_parse_REACTION_EQUATION <- function(EQUATION){


  eq.syn <- str_extract(pattern = "[<=>]+",EQUATION)
  dir<- switch (eq.syn,
                "<=>" = 2,
                "=>" = 1,
                "<=" = -1
  )

  eq.split <- str_split(EQUATION,eq.syn)[[1]]
  eq.l <-  stringr::str_extract_all(eq.split[1],"C[0-9]{5}")[[1]]
  eq.r <-  stringr::str_extract_all(eq.split[2],"C[0-9]{5}")[[1]]


  return(list(
    from = eq.l,
    to = eq.r,
    direction = dir
  ))


}


KEGG_parse_MODULE_REACTION <- function(REACTION){


}

KEGG_parse_REACTION_RCLASS <- function(RCLASS){



  stringr::str_extract_all(pattern = "[R]{0,1}C[0-9]{5}",RCLASS)%>%
    lapply(function(x){

      data.frame(RCLASS_id = x[1],
                 from = x[2],
                 to = x[3])
    })%>%
    data.table::rbindlist()

}

KEGG_reaction_to_network <- function(kegg.rdata,
                                     filter_org = "hsa"){



  ### filter_org
  {
    if(!is.null(filter_org)){

      kegg.link <- KEGGREST::keggLink("enzyme",filter_org)%>%
        FELLA:::sanitise("enzyme",filter_org)
      is.org <- sapply(kegg.rdata,function(reaction.data){
        any(reaction.data[["ENZYME"]] %in% kegg.link)
      })
      kegg.rdata.filtered <- kegg.rdata[is.org]
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
                        name = reaction.data$NAME,
                        equation = reaction.data$EQUATION)
        enzn <- rep(enz,nrow(rda))
        rdan <- rda[rep(1:nrow(rda),each = length(enz)),]
        rcda <- data.frame(rdan,
                   ENZYME = enzn)
        return(rcda)

      },.progress = "text"
    )
    rcn.edges.df <-do.call(bind_rows,rcn.edges)%>%
      dplyr::select(from,to,everything())

    if (!is.null(filter_org)) {
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

KEGG_reaction_get_direction_from_module <-
  function(ig.krn,kegg.mdata){



}

KEGG_Reaction_network_add_label <- function(kegg.rig){

  eda <- edata(kegg.rig)%>%
    dplyr::mutate(label = kegg_id)


  eda -> edata(kegg.rig)
  return(kegg.rig)
}

KEGG_Reaction_network_merge_path <- function(kegg.rig){

  eda <- edata(kegg.rig)%>%
    dplyr::mutate(path_str = paste0(from,"_",to))

  edge.to.remove <- which(duplicated(eda$path_str))
  kegg.rig <- igraph::delete_edges(kegg.rig,edge.to.remove)
  return(kegg.rig)
}

KEGG_Reaction_network_remove_nonformat_node <- function(kegg.rig){

  vda <- vdata(kegg.rig)%>%
    dplyr::mutate(formate_formula = MSCC::chemform_formate(Formula))
  kegg.rig<- igraph_filter_vertex(kegg.rig,!is.na(vda$formate_formula))


  return(kegg.rig)
}



get_KRN_edge <- function(KRN){


}

