ggplot_sdf <- function(sdf,
                       cex = 1,
                       show_ele = F){


  sdf.formula <- MF(sdf,addH=T)
  sdf.mz <- chemform_mz(sdf.formula)%>%round(digits = 4)
  atom.data <- atomblock(sdf)[,1:2]%>%
    `colnames<-`(c("x","y"))%>%
    as.data.frame()%>%
    rownames_to_column("Atom_id" )%>%
    dplyr::mutate(element = str_extract(Atom_id,
                                        "[:alpha:]*"))
  bond.length.short <- ifelse(show_ele,0.1,0)
  bond.data <- bondblock(sdf)[,1:3]%>%
    `colnames<-`(c("from","to","bond_type"))%>%
    as.data.frame()%>%
    dplyr::mutate(
      bond_id = 1:n(),
      x = atom.data$x[from],
      xend = atom.data$x[to],
      y = atom.data$y[from],
      yend = atom.data$y[to]
    )%>%
    dplyr::mutate(xl = (xend-x),
                  yl = (yend - y),
                  x = x + bond.length.short*xl,
                  xend = xend - bond.length.short*xl,
                  y = y+bond.length.short*yl,
                  yend = yend - bond.length.short*yl)
  for (i in 1:nrow(bond.data)) {
    bond.data <- dplyr_copy_row(bond.data,
                                i,
                                bond.data$bond_type[i]-1)
  }

  lw <- 0.5*cex
  sw <- 0.5*cex
  col.bond <- "#666666"
  ggplot()+
    ### 3 bond
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw+2*sw+2*lw)+
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = "white",linewidth = lw+2*sw)+
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw)+
    ### 2 bond
    geom_segment(data = filter(bond.data,bond_type == 2),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = 2*lw+sw)+
    geom_segment(data = filter(bond.data,bond_type == 2),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = "white",linewidth = sw)+
    ### 1 bond
    geom_segment(data = filter(bond.data,bond_type == 1),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw)+
    geom_text(aes(x = median(range(atom.data$x)),
                  y = max(atom.data$y)+diff(range(atom.data$y))*0.5,
                  label = paste(sdf.formula,"\n",sdf.mz)),
              size = 2)+
    ylim(c(min(atom.data$y),max(atom.data$y)+diff(range(atom.data$y))*0.8))+
    xlim(expand_range(range(atom.data$x),multi = 0.2))+
    theme_void()->p


  if (show_ele) {
    p <- p+geom_text(data = atom.data,
                aes(x = x, y = y ,label = element),
                size = 2 *cex)
  }else{
    p <- p+geom_point(data = atom.data,
                     aes(x = x, y = y ),
                     size = 0.5*cex)
  }
  p
  return(p)
}

check_sdf <- function(sdf){

  atom.matrix <- atomcountMA(sdf)
  atom.matrix <- atom.matrix[,setdiff(colnames(atom.matrix),"0"),drop =F]
  apply(atom.matrix,1,sum) >0

}

check_smile <- function(smile){

  smile.sdf <- get_smile_sdf(smile)
  check_sdf(smile.sdf)

}

get_sdf_formula <- function(sdf){

  sdf.checked <- check_sdf(sdf)
  sdf.formula <- character()
  sdf.formula[sdf.checked] <- MF(sdf[sdf.checked],addH=T)
  sdf.formula <- MSCC::chemform_formate(sdf.formula)
  return(sdf.formula)
}

#' get_smile_sdf
#'
#' @param smiles  smiles
#' @param smiles.id NULL
#'
#' @return sdf
#' @export
#' @import ChemmineR
get_smiles_sdf <- function(smiles,smiles.id = NULL){

  data("smiles_map")
  if (is.null(names(smiles))) {
    if (is.null(smiles.id)) {
      names(smiles) <- paste0("CMP",num2str(seq_along(smiles)))
    }else{
      names(smiles) <- smiles.id
    }

  }
  smiles.sdf <- suppressWarnings(
    ChemmineR::smiles2sdf(smiles)
  )
  for (id in ChemmineR::cid(smiles_map)) {
    which(smiles==id)
    suppressWarnings(smiles.sdf[smiles==id] <- smiles_map[[id]])
  }

  return(smiles.sdf)

}

get_smile_formula <- function(smile){

  smile.sdf <- get_smile_sdf(smile)
  smile.formula <- get_sdf_formula(smile.sdf)
  smile.formula <- case_when(smile=="O"~"H2O1",
                 smile=="[HH]"~"H2",
                 T~smile.formula)

  return(smile.formula)
}

#' get_sdf_igraph
#'
#' @param sdf sdf
#' @param addH T or F
#'
#' @return igraph
#' @export
#' @import dplyr tibble ChemmineR igraph
get_sdf_igraph <- function(sdf,addH = F){

  atom.data <- atomblock(sdf)[,1:2]%>%
    `colnames<-`(c("x","y"))%>%
    as.data.frame()%>%
    rownames_to_column("Atom_id" )
  atom.data <- cbind(atom.data,bonds(sdf))%>%
    dplyr::group_by(atom)%>%
    dplyr::mutate(id = Atom_id,
                  label = paste0(atom,1:n()),
                  shape = "text")

  bond.data <- bondblock(sdf)[,1:3,drop =F]%>%
    `colnames<-`(c("from","to","bond_type"))%>%
    as.data.frame()%>%
    dplyr::mutate(from = atom.data$Atom_id[from],
                  to = atom.data$Atom_id[to],
                  width = 10*bond_type,
                  color = c("#F8C959","#FF8F6B","#D35230")[bond_type])

  sdf.igraph <- graph_from_data_frame(
    bond.data,vertices =atom.data
  )
  if (!addH) {
    sdf.igraph <- delete.vertices(sdf.igraph,
                                  V(sdf.igraph)$atom=="H")
  }

  return(sdf.igraph)

}

#' add_sdf_igraph_highlight
#'
#' @param sdf.igraph igraph
#' @param highlight vector
#'
#' @return igraph
#' @export
#' @import  igraph
add_sdf_igraph_highlight <- function(sdf.igraph,
                                     highlight =NULL){

  sdf.igraph.highlight <- sdf.igraph
  if (!is.numeric(highlight)) {
    highlight <- match(highlight  ,V(sdf.igraph.highlight)$id)
  }
  edge.data <- igraph::as_data_frame(sdf.igraph.highlight)%>%
    dplyr::mutate(highlight = from %in% V(sdf.igraph.highlight)$name[highlight] &
                    to %in% V(sdf.igraph.highlight)$name[highlight])

  E(sdf.igraph.highlight)$width[edge.data$highlight] <-
    E(sdf.igraph.highlight)$width[edge.data$highlight]*2
  E(sdf.igraph.highlight)$color[edge.data$highlight] <- "#2B7CE9"
  V(sdf.igraph.highlight)$shadow <- F
  V(sdf.igraph.highlight)$background <- F
  V(sdf.igraph.highlight)$shape[highlight] <- "circle"
  V(sdf.igraph.highlight)$shadow[highlight] <- T

  return(sdf.igraph.highlight)
}

#' vis_sdf_igraph
#'
#' @param sdf.igraph igraph
#' @param show.label logic
#' @param highlight vector
#'
#' @return vis html
#' @export
#' @import igraph visNetwork
vis_sdf_igraph <- function(sdf.igraph ,
                           show.label = T,
                           highlight = NULL){



  sdf.igraph.highlight <- add_sdf_igraph_highlight(sdf.igraph,
                                                   highlight )

  sdf.igraph.highlight%>%
    visIgraph(idToLabel = !show.label)%>%
    visNodes(font = list(size = 40,
                         strokeWidth = 10),
             size = 60,
             borderWidth  = 5,
             color = list(background = "transparent",
                          border = "#2B7CE9"))%>%
    visEdges(arrows = list(to = F),
             length = 0.8)

}

vis_sdf_igraph_compare <- function(sdf.igraph1,
                                   sdf.igraph2,
                                   atom_id1=NULL,
                                   atom_id2=NULL){

  sdf.igraph1 <- add_sdf_igraph_highlight(sdf.igraph1,atom_id1)
  sdf.igraph2 <- add_sdf_igraph_highlight(sdf.igraph2,atom_id2)
  hw1 <- diff(range(V(sdf.igraph1)$x))
  hw2 <- diff(range(V(sdf.igraph2)$x))
  V(sdf.igraph1)$x <- V(sdf.igraph1)$x-hw1/1.5
  V(sdf.igraph2)$x <- V(sdf.igraph2)$x+hw2/1.5
  V(sdf.igraph1)$y <- V(sdf.igraph1)$y/2
  V(sdf.igraph2)$y <- V(sdf.igraph2)$y/2
  V(sdf.igraph1)$name <- paste0(V(sdf.igraph1)$name,"_1")
  V(sdf.igraph2)$name <- paste0(V(sdf.igraph2)$name,"_2")
  V(sdf.igraph1)$id <- paste0(V(sdf.igraph1)$id,"_1")
  V(sdf.igraph2)$id <- paste0(V(sdf.igraph2)$id,"_2")

  nodes <- rbind(
    as_data_frame(sdf.igraph1,"V"),
    as_data_frame(sdf.igraph2,"V")
  )
  edges <- rbind(
    as_data_frame(sdf.igraph1,"E"),
    as_data_frame(sdf.igraph2,"E")
  )
  new.igraph <- graph_from_data_frame(edges,
                                      vertices = nodes)
  vis_sdf_igraph(new.igraph)

}

get_atom_id_from_parent <- function(parent.sdf.graph,
                                    product.sdf.graph){

  #parent.sdf.graph <- fragment.igraph[[1]]
  #product.sdf.graph <- fragment.igraph[[3]]

  old.root.id <- V(product.sdf.graph)$root_atom_id
  V(product.sdf.graph)$root_atom_id <- "unknown"
  ig <- intersection(parent.sdf.graph,
                     product.sdf.graph,
                     byname = F, keep.all.vertices = F)

  ig <- ig - V(ig)[atom_1!=atom_2]

  ig.vd <- as_data_frame(ig,"vertices")
  rownames(ig.vd) <- (ig.vd$name_2)

  #vis_sdf_igraph(parent.sdf.graph,show.label = F)
  #vis_sdf_igraph(product.sdf.graph,show.label = F)
  new.id <- ig.vd[V(product.sdf.graph)$name,]$root_atom_id_1
  if (!is.null(old.root.id)) {
    new.id <- case_when(is.na(new.id)~old.root.id,
              new.id == "unknown"~old.root.id,
              T~new.id)
  }

  V(product.sdf.graph)$root_atom_id <- new.id

  return(product.sdf.graph)

}

