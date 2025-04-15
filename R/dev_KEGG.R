get_KEGG_Reaction_network <- function(kegg.rdata ){


  kegg.rawdata <- MSdb:::get_KEGG_rawdata()
  kegg.rdata <-  kegg.rawdata$Reaction_rawdata
  kegg.rdata <- KEGG_parse_REACTION(kegg.rdata,parse_by = "EQUATION")
  #ig.krn <- KEGG_reaction_to_network_RCLASS(kegg.rdata)
  #kegg.reaction.check <- KEGG_EQUATION_data_check(kegg.rdata)

  ig.krn <- KEGG_reaction_to_network_EQUATION(kegg.rdata)

  #kegg.mdata <- kegg.rawdata$Module_rawdata
  #ig.krn <- KEGG_reaction_get_direction_from_module(ig.krn,kegg.mdata)
  #ig.krn <- igraph_sort_direction(ig.krn)
  #ig.krn.filter <- KEGG_Reaction_network_filter_by_emzyme(ig.krn,"hsa")
  #ig.krn.filter <- KEGG_Reaction_network_remove_nonformat_node(ig.krn)
  ig.krn.filter <- ig.krn
  return(ig.krn.filter)


}



KEGG_parse_REACTION <- function(kegg.rdata,
                                parse_by = "EQUATION" ){

  vars <- c("ENTRY",
            "NAME",
            "DEFINITION","EQUATION",
            "ENZYME","RCLASS",
            "PATHWAY","MODULE"
            # "BRITE","DBLINKS","COMMENT","ORTHOLOGY","REMARK","REFERENCE"
  )


  ### EQUATION
  {

    if(parse_by=="EQUATION"){

      #reaction.data <- kegg.rdata[[1]]
      message_with_time("Parsing EQUATION...")
      ### MAP
      {
        kegg.gly.cp.map <- KEGGREST::keggLink("compound","glycan")
        names(kegg.gly.cp.map) <- sub(x = names(kegg.gly.cp.map),pattern = "gl:",replacement = "")
        kegg.gly.cp.map <- sub(x = (kegg.gly.cp.map),pattern = "cpd:",replacement = "")

        kegg.cp <- MSdb:::get_KEGG_compound_df()
        kegg.formula <- make_vector(kegg.cp$Formula,kegg.cp$KEGG_id)
      }
      kegg.rdata <- plyr::llply(kegg.rdata,function(reaction.data){

        ### data process
        {
          ENTRY <- unname(reaction.data$ENTRY)
          EQUATION_data <- KEGG_parse_REACTION_EQUATION(reaction.data[["EQUATION"]])
          if (is.null(EQUATION_data)) return(reaction.data)

          EQUATION_data$id <- ifelse(EQUATION_data$id%in% names(kegg.gly.cp.map),
                                     kegg.gly.cp.map[EQUATION_data$id],
                                     EQUATION_data$id)

        }

        ### edge df
        {
          from.df <- EQUATION_data%>%
            dplyr::filter(side =="from")
          to.df<- EQUATION_data%>%
            dplyr::filter(side =="to")
          rda1 <- data.frame(
            from = from.df$id,
            to = ENTRY
          )
          rda2 <- data.frame(
            from = ENTRY ,
            to = to.df$id
          )
          edge.df <- rbind(rda1,rda2)%>%
            dplyr::mutate(REACTION_id = reaction.data$ENTRY)

        }

        ### Reaction node
        {
          node.reaction.data <- list(
            KEGG_id = ENTRY,
            Name = reaction.data$NAME[1],
            EQUATION = reaction.data$EQUATION,
            equation.format = paste0(
              paste0(from.df$coef," ",from.df$id,collapse = " + ")," <=> ",
              paste0(to.df$coef," ",to.df$id,collapse = " + ")
            ),
            equation.formula = paste0(
              paste0(from.df$coef," ",kegg.formula[from.df$id],collapse = " + ")," <=> ",
              paste0(to.df$coef," ",kegg.formula[to.df$id],collapse = " + ")
            ),
            definition = reaction.data$DEFINITION,
            enzyme = paste0(reaction.data$ENZYME,collapse   = ";")
          )

        }


        reaction.data$EQUATION_data <- EQUATION_data
        reaction.data$NODE_data <- node.reaction.data
        reaction.data$EDGE_data <- edge.df
        reaction.data



      },.progress = "text")

    }


  }



  ### RCLASS
  {

    if(parse_by=="RCLASS"){

      #reaction.data <- kegg.rdata[[1]]
      RCLASS_data <- KEGG_parse_REACTION_RCLASS(reaction.data[["RCLASS"]])
      reaction.data$RCLASS_data <- RCLASS_data

    }


  }



  return(kegg.rdata)

}



KEGG_EQUATION_data_check <- function(kegg.rdata){


  kegg.stat <- lapply(kegg.rdata,function(reaction.data){

    if (anyNA(reaction.data)) {
      return(NULL)
    }

    data.frame(
      ENTRY = unname(reaction.data$ENTRY),
      EQUATION = reaction.data$EQUATION,
      from.count = sum(reaction.data$EQUATION_data$side=="from",na.rm = T),
      to.count = sum(reaction.data$EQUATION_data$side=="to",na.rm = T),
      na.count = sum(is.na(reaction.data$EQUATION_data$id),na.rm = T)
    )%>%
      dplyr::mutate(
        error = case_when(
          from.count==0|to.count ==0~"From or to missing",
          na.count > 0 ~ "Char no recongnized",
          grepl( "m|n",EQUATION) ~ "m or n coef",
          T ~ NA
        )
      )
  })
  kegg.stat <- do.call(rbind,kegg.stat)
  return(kegg.stat)
}

KEGG_parse_REACTION_EQUATION <- function(EQUATION){


  if (grepl("m|n",EQUATION)) {
    return(NULL)
  }

  eq.syn <- str_extract(pattern = "[<=>]+",EQUATION)
 # dir <- switch (eq.syn,
 #               "<=>" = 0,
 #               "=>" = 1,
 #               "<=" = -1
 # )

  eq.split <- str_split(EQUATION,eq.syn)[[1]]%>%
    lapply(function(y){
      str_split(y,"\\+",simplify = T)%>%
        gsub(" ","",.)%>%
        lapply(function(x){
          x.coef <- str_extract(x,"^[:digit:]+")
          x.coef <- ifelse(is.na(x.coef),1,x.coef)
          x.coef <- as.numeric(x.coef)
          x.id <- str_extract(x,"[CG][0-9]{5}")
          #make_vector(x.coef,x.id)
          data.frame(id = x.id,
                     coef = x.coef)
        })%>%do.call(rbind,.)
    })
  names(eq.split) <- c("from","to")
  equation.df <- data.table::rbindlist(eq.split,idcol = "side")



  return(equation.df)


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


KEGG_reaction_to_network_EQUATION <- function(kegg.rdata){



  ### edge
  {
    message_with_time("Generate Edges...")
    rcn.edges <- plyr::llply(
      kegg.rdata,
      function(reaction.data){

        if (!"EQUATION_data" %in% names(reaction.data)) {
          return(invisible())
        }
        return(reaction.data$EDGE_data )

      },.progress = "text"
    )
    rcn.edges.df <-do.call(bind_rows,rcn.edges)%>%
      dplyr::select(from,to,everything())%>%
      dplyr::mutate(name = paste0("UR",num2str(1:n())))

  }



  ### node
  {
    # compound
    {
      message_with_time("Generate Nodes")
      kegg.cp <- MSdb:::get_KEGG_compound_df()%>%
        as.data.frame()%>%
        dplyr::filter(KEGG_id %in% rcn.edges.df$from|
                        KEGG_id %in% rcn.edges.df$to)%>%
        dplyr::mutate(node.type = "Compound",
                      name = KEGG_id,
                      formula.format = MSCC::chemform_formate(Formula))
    }


    # reaction
    {
      kegg.reaction <- plyr::llply(
        kegg.rdata,
        function(reaction.data){

          if (is.null(reaction.data[["EQUATION_data"]])) return(NULL)
          return(reaction.data$NODE_data)

        },.progress = "text")

      kegg.reaction.df <- do.call(bind_rows,kegg.reaction)%>%
        dplyr::filter(KEGG_id %in% rcn.edges.df$from|
                        KEGG_id %in% rcn.edges.df$to)%>%
        dplyr::mutate(node.type = "Reaction"  )


    }


    rcn.nodes.df <- bind_rows(kegg.cp,kegg.reaction.df)%>%
      dplyr::mutate(id= KEGG_id,
                    name = KEGG_id,
                    label = Name)
   #rcn.nodes.df.miss <- data.frame(
   #  name =unique( setdiff(c(rcn.edges.df$from,rcn.edges.df$to),rcn.nodes.df$name))
   #)%>%
   #  dplyr::mutate(id = name)
   #rcn.nodes.df <- bind_rows(rcn.nodes.df,rcn.nodes.df.miss)
    kegg.rdata.error <- kegg.rdata[which(!names(kegg.rdata)%in%kegg.reaction.df$KEGG_id)]

  }


  ### ig
  {
    rcn.edges.df <- rcn.edges.df%>%
      dplyr::filter( from %in% rcn.nodes.df$id  & to %in% rcn.nodes.df$id)
    kegg.reaction.ig <- igraph::graph_from_data_frame(
      rcn.edges.df,
      vertices = rcn.nodes.df)

  }

  return(kegg.reaction.ig)

}

KEGG_reaction_get_direction_from_module <-
  function(ig.krn,kegg.mdata){

    reaction.exist <- sapply(kegg.mdata,
                             function(x){
                               "REACTION" %in% names(x)
                             })
    kegg.mdata <- kegg.mdata[reaction.exist]
    kegg.mdata.df <- plyr::llply(
      kegg.mdata,
      function(mdata){

        #message(mdata$ENTRY)
      mreaction <- mdata[["REACTION"]]
       lapply(seq_along(mreaction),function(i){
            KEGG_parse_MODULE_REACTION(mreaction[i])
        })%>%
          do.call(rbind,.)
      },.progress = "text"
    )%>%
      data.table::rbindlist(,idcol = "MODULE_id")


    ### process multiple reaction
    {

      mr.id <- grepl("\\+",kegg.mdata.df$REACTION_id)
      kegg.mdata.df.split <- split.data.frame(kegg.mdata.df,mr.id)
      mr.df <- kegg.mdata.df.split$`TRUE`
      mr.reaction <- list()
      for (i in seq_len(nrow(mr.df))) {
        this.rid <- mr.df$REACTION_id[i] %>%
          stringr::str_split("\\+")%>%unlist()
        this.eda <- edata(ig.krn)%>%
          dplyr::filter(REACTION_id %in% this.rid)
        this.eda$direction <- NA
        this.from <- mr.df$from[i]
        this.to <- mr.df$to[i]

        this.ig <- igraph_filter_edge(ig.krn,which(edata(ig.krn)$REACTION_id %in% this.rid))
        this.vpath <- igraph::all_simple_paths(this.ig,this.from,this.to,mode = "all")
        this.dir <- lapply(this.vpath,function(vpath){
          # Convert path to numeric vector
          path <- as.numeric(vpath)

          # Find edge IDs for consecutive vertex pairs
          edge_rev <- sapply(seq_along(path)[-length(path)], function(i) {
            igraph::get_edge_ids(this.ig, c(path[i], path[i + 1]),directed = F)
          })
          edges <- sapply(seq_along(path)[-length(path)], function(i) {
            igraph::get_edge_ids(this.ig, c(path[i], path[i + 1]),directed = T)
          })
          edges <- edges[edges!=0]
          edge_rev <- setdiff(edge_rev,edges)
          this.eda$direction[edges] <- 1
          this.eda$direction[edge_rev] <- -1
          return(this.eda$direction)
        })
        this.eda <- this.eda%>%
          dplyr::mutate(direction = this.dir[[1]],
                        MODULE_id = mr.df$MODULE_id[i])%>%
          dplyr::filter(!is.na(direction)
                        )%>%
          dplyr::select(c("MODULE_id","from","to","REACTION_id","direction"))

        mr.reaction[[i]] <- this.eda
      }

      mr.reaction <- data.table::rbindlist(mr.reaction)
      kegg.mdata.df <- rbind(
        kegg.mdata.df.split$`FALSE`,
        mr.reaction
      )
      }

    ### modify krn
    {

      kegg.mdata.df <- kegg.mdata.df%>%
        dplyr::mutate(str_syn =
                        paste0(REACTION_id,"_",
                               from,"_",
                               to))
      eda <- edata(ig.krn)%>%
        dplyr::mutate(str_syn =
                        paste0(REACTION_id,"_",
                               from,"_",
                               to),
                      str_syn_rev =
                        paste0(REACTION_id,"_",
                               to,"_",
                               from),
                      direction = case_when(
                        str_syn%in%kegg.mdata.df$str_syn~kegg.mdata.df$direction[
                          match(str_syn,kegg.mdata.df$str_syn)],
                        str_syn_rev%in%kegg.mdata.df$str_syn~0 - kegg.mdata.df$direction[
                          match(str_syn_rev,kegg.mdata.df$str_syn)
                        ]
                      ) )%>%
        dplyr::select(-str_syn,-str_syn_rev)
      eda -> edata(ig.krn)

    }


    return(ig.krn)




  }


KEGG_parse_MODULE_REACTION <- function(REACTION){


  ### Reaction id
  {
    rid <- names(REACTION)
    ###fix
    {

      ### names(REACTION) miss
      if(length(rid)==0||rid == ""){
        x <- sub(" ","&",REACTION)%>%
          stringr::str_split("&")
        REACTION <- setNames(nm = x[[1]][1],x[[1]][2])
        rid <- names(REACTION)
      }

      ### names(REACTION) contain reaction
      if(grepl("R\\d{5}",rid)&grepl("[CG]\\d{5}",rid)){
        x <- sub(" ","&",rid)%>%
          stringr::str_split("&")
        REACTION <- setNames(nm = x[[1]][1],x[[1]][2])
        rid <- names(REACTION)
      }
    }
    rid.split <- stringr::str_split(rid,",")[[1]]

    ### check
    {

      is.mr <- grepl(pattern = "+",x = rid.split,fixed = T)
      nonformat <- rid.split[!grepl("^R\\d{5}$",rid.split)]
      nonformat <- stringr::str_split(nonformat,"\\+")%>%unlist()
      if(!all(grepl("^R\\d{5}$",nonformat))){
        cli::cli_inform(
          "KEGG_parse_MODULE_REACTION error {cli::col_red(names(REACTION))} : {cli::col_blue(REACTION)}"
        )
      }
    }
  }


  ### reaction
  {

    dir.symbol <- stringr::str_extract(REACTION," [-<>]+ ")
    direction <- switch (dir.symbol,
      " <- " = -1,
      " <-> " = 0,
      " -> " = 1,
    )
    r.split <- stringr::str_split(REACTION,dir.symbol)[[1]]
    r.split <- stringr::str_extract_all(
      r.split,"[GC]\\d{5}"
    )
    if (length(r.split)!=2) {
      cli::cli_inform("  KEGG_parse_MODULE_REACTION error for {REACTION}")
    }
    names(r.split) <- c("from","to")

    r.split$REACTION_id <- rid.split
    r.split$stringsAsFactors <- F
    r.df <- do.call(expand.grid,r.split)
    r.df$direction <- direction
  }


  return(r.df)


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
    dplyr::mutate(formate_formula = MSCC::chemform_formate(Formula),
                  remain = !is.na(formate_formula)|node.type == "Reaction")
  kegg.rig<- igraph_filter_vertex(kegg.rig,vda$remain)


  return(kegg.rig)
}



KEGG_Reaction_network_filter_by_emzyme <- function(ig.krn ,org = "hsa"){

  enz.link <- KEGGREST::keggLink("enzyme",org)%>%
    sub(pattern = 'ec:',x = . , replacement = "")
  vda <- vdata(ig.krn)
  ig.krn.filter <- igraph_filter_vertex(ig.krn,
                       which(vda$node.type== "Compound" | vda$enzyme %in% enz.link))

  ig.krn.filter <- igraph_filter_vertex(ig.krn.filter,
                                        degree(ig.krn.filter)>0)
  return(ig.krn.filter)
}

get_KRN_edge <- function(KRN){


}

get_KEGG_MODULE_reaction_direction_stat <- function(kegg.mdata){

  all.m.reaction <- lapply(kegg.mdata,
                           function(mdata){
                             mdata$REACTION
                           })%>%
    unname()%>%
    unlist()

}



KEGG_get_cp_linked_gene <- function(cp.id){

  cp.link <- KEGGREST::keggLink("compound","enzyme")
  cp.link <- sub(x = cp.link,pattern = "cpd:",replacement = "")
  names(cp.link) <- sub(x = names(cp.link),pattern = "ec:",replacement = "")


  cp.link <- cp.link[cp.link%in% cp.id]

  map <- clusterProfiler::bitr(names(cp.link),fromType = "ENZYME",
                        toType = "PMID",OrgDb = org.Hs.eg.db)
  map.matched <- map%>%
    dplyr::mutate(cp.id = cp.link[ENZYME])


  symbol <- sapply(cp.id,function(x){
    if (is.na(x)) {
      return(NA)
    }
    paste0(map.matched$SYMBOL[map.matched$cp.id==x],collapse  = ";")

  })

  enzyme <- sapply(cp.id,function(x){
    if (is.na(x)) {
      return(NA)
    }
    paste0( names(cp.link)[cp.link==x],collapse  = ";")

  })

  data.frame(SYMBOL = symbol,
             ENZYME = enzyme) %>%
    return()



}
