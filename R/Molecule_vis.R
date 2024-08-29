#' get_sdf_igraph
#'
#' @param sdf sdf
#' @param addH T or F
#'
#' @return igraph
#' @export
#' @import dplyr tibble ChemmineR igraph
get_sdf_igraph <- function(sdf,addH = F){

  .f <- function(sdf,addH){
    atom.data <- atomblock(sdf)[,1:2]%>%
      `colnames<-`(c("x","y"))%>%
      as.data.frame()%>%
      rownames_to_column("Atom_id" )
    atom.data <- cbind(atom.data,bonds(sdf))%>%
      dplyr::group_by(atom)%>%
      dplyr::mutate(id = Atom_id,
                    label = paste0(" ",atom," "),
                    x = x*100,y = y *100,
                    font.size = 30,
                    borderWidth = 20,
                    font.vadjust = 5,
                    font.align = "center",
                    color.border = "#AAAAAA",
                    color.background = "#FFFFFF",
                    borderWidth = 5,
                    shape = "circle",
                    physics = F)

    bond.data <- bondblock(sdf)[,1:3,drop =F]%>%
      `colnames<-`(c("from","to","bond_type"))%>%
      as.data.frame()%>%
      dplyr::mutate(from = atom.data$Atom_id[from],
                    to = atom.data$Atom_id[to],
                    width = 10*bond_type,
                    color = c("#888888","#444444","#000000")[bond_type],
                    smooth = FALSE
                    )

    sdf.igraph <- graph_from_data_frame(
      bond.data,vertices = atom.data
    )
    if (!addH) {
      sdf.igraph <- delete.vertices(sdf.igraph,
                                    V(sdf.igraph)$atom=="H")
    }

    return(sdf.igraph)
  }

  if (class(sdf)=="SDF")
    sdf.igraph <-.f(sdf,addH)


  if (class(sdf)=="SDFset") {
    sdf.igraph <- list()
    sdf.valid <- validSDF(sdf)
    for (i in 1:length(sdf)) {

      sdf.igraph[[i]] <- .f(sdf[[i]],addH  )

    }
  }

  return(sdf.igraph)


}

#' vis_sdf_igraph
#'
#' @param sdf.igraph igraph
#' @param show_id logic
#' @param highlight vector
#'
#' @return vis html
#' @export
#' @import igraph visNetwork
vis_sdf_igraph <- function(sdf.igraph,show_id = F,...){

  sdf.igraph <- sdf.igraph%>%
    sdf_igraph_show_id(show_id)
  vda <- vdata(sdf.igraph)%>%
    dplyr::mutate(x= x-mean(x),
                  y = y-mean(y))
  eda <- edata(sdf.igraph)

  visNetwork(nodes = vda,edges = eda) %>%
    visLayout( hierarchical = FALSE)

}

vis_cfm_data_fragment <- function(cfmd,fragment_id,
                                  show_id= T,...){

  ig <- cfmd@fragment_igraph[[fragment_id]]%>%
    sdf_igraph_show_id(show_id )
  vis_sdf_igraph(ig,show_id = show_id,...)

}

vis_cfm_data_fragment_atom_map<- function(cfmd,fragment_id,
                                          show_id= T){

  maps <- cfmd@fragment_atom_map[[fragment_id]]

  sdf.igraphA <- get_cfm_data_sdf_igraph(cfmd,1)%>%
    sdf_igraph_add_border_color(value = rowSums(maps))%>%
    sdf_igraph_show_id(show_id)
  sdf.igraphB <- get_cfm_data_sdf_igraph(cfmd,fragment_id)%>%
    sdf_igraph_add_border_color(value = colSums(maps))%>%
    sdf_igraph_show_id(show_id)

  sdf_igraph_merge(sdf.igraphA,sdf.igraphB)%>%
    vis_sdf_igraph(show_id = show_id)

}

vis_cfm_data_trans_map<- function(cfmd,trans_id,
                                          show_id= T){

  maps <-get_CFM_data_trans_map(cfmd,trans_id )

  fragment.trans <- cfmd@fragment_transition
  sdf.igraphA <- get_cfm_data_sdf_igraph(cfmd,fragment.trans$from[trans_id])%>%
    sdf_igraph_add_border_color(value = rowSums(maps))%>%
    sdf_igraph_show_id(show_id)
  sdf.igraphB <- get_cfm_data_sdf_igraph(cfmd,fragment.trans$to[trans_id])%>%
    sdf_igraph_add_border_color(value = colSums(maps))%>%
    sdf_igraph_show_id(show_id)

  sdf_igraph_merge(sdf.igraphA,sdf.igraphB)%>%
    vis_sdf_igraph(show_id = show_id)

}

sdf_igraph_add_border_color <- function(sdf.igraph,
                                        value,
                                        color.ramp = colramp(breaks = c(0,Inf,1),
                                                             colors = c("#aaaaaa","#97C2FC","#2B7CE9"))
                                        ){

  sdf.igraph.temp <-sdf.igraph
  ele <-  get_sdf_igraph_atom(sdf.igraph)
  col.border <- .get_vis_col(sdf.igraph,value,color.ramp)
  vdata(sdf.igraph.temp)$color.border <- col.border
  return(sdf.igraph.temp)


}
sdf_igraph_add_background_color <- function(sdf.igraph,
                                        value,
                                        color.ramp = colramp(breaks = c(0,Inf,1),
                                                             colors = c("#FFFFFF","#F7844F","#B20C26"))
){

  sdf.igraph.temp <-sdf.igraph
  ele <-  get_sdf_igraph_atom(sdf.igraph)
  col.background<- .get_vis_col(sdf.igraph,value,color.ramp)
  vdata(sdf.igraph.temp)$color.background <- col.background
  return(sdf.igraph.temp)


}


sdf_igraph_show_id <- function(sdf.igraph,show_id){

  message_with_time(show_id)
  if (show_id) {
    vdata(sdf.igraph)$label <-   vdata(sdf.igraph)$id%>%
      str_format_len()
    vdata(sdf.igraph)$font.size <- 20

  }else{
    vdata(sdf.igraph)$label <-  paste0(" ",vdata(sdf.igraph)$atom," ")
  }

  return(sdf.igraph)
}

sdf_igraph_merge <- function(sdf.igraphA,sdf.igraphB){


  xa <-vdata(sdf.igraphA)$x
  xb <- vdata(sdf.igraphB)$x

  vdata(sdf.igraphA)$x <- xa-max(xa)-diff(range(xa))/5
  vdata(sdf.igraphB)$x <- xb- min(xb)+diff(range(xb))/5

  vdata(sdf.igraphA)$name <- paste0("A_",vdata(sdf.igraphA)$name)
  vdata(sdf.igraphB)$name <- paste0("B_",vdata(sdf.igraphB)$name)

  vdata(sdf.igraphA)$id <- paste0("A_",vdata(sdf.igraphA)$id)
  vdata(sdf.igraphB)$id <- paste0("B_",vdata(sdf.igraphB)$id)

  nodes <- bind_rows(
    vdata(sdf.igraphA),
    vdata(sdf.igraphB)
  )
  #nodes$size <- 20
  edges <- rbind(
    edata(sdf.igraphA),
    edata(sdf.igraphB)
  )
  new.igraph <- graph_from_data_frame(edges,
                                      vertices = nodes)
  return(new.igraph)
}



.get_vis_col <- function(sdf.igraph,prob,colramp = colramp(),na.col = "#AAAAAA"){

  ele <- get_sdf_igraph_atom(sdf.igraph,ele = "all")
  if (is.null(prob)) prob <- 0
  if (is.null(names(prob))&length(prob)==1){
    prob <- rep(prob,length(ele))
    names(prob)<-ele
  }else{
    x <- setdiff(ele,names(prob))
    xx <- rep(0,length(x))
    names(xx) <-x
    prob <- c(prob,xx)[ele]
  }
  col <- colramp(prob)
  names(col)<-names(prob)
  col[is.na(col)] <- na.col
  col
}
