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

  atom.map.matrix <- atomcountMA(sdf)
  atom.map.matrix <- atom.map.matrix[,setdiff(colnames(atom.map.matrix),"0"),drop =F]
  apply(atom.map.matrix,1,sum) >0

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

  .f <- function(sdf,addH){
    atom.data <- atomblock(sdf)[,1:2]%>%
      `colnames<-`(c("x","y"))%>%
      as.data.frame()%>%
      rownames_to_column("Atom_id" )
    atom.data <- cbind(atom.data,bonds(sdf))%>%
      dplyr::group_by(atom)%>%
      dplyr::mutate(id = Atom_id,
                    label = paste0(atom),
                    shape = "text")

    bond.data <- bondblock(sdf)[,1:3,drop =F]%>%
      `colnames<-`(c("from","to","bond_type"))%>%
      as.data.frame()%>%
      dplyr::mutate(from = atom.data$Atom_id[from],
                    to = atom.data$Atom_id[to],
                    width = 10*bond_type,
                    color = c("#F8C959","#FF8F6B","#D35230")[bond_type])

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
    for (i in 1:length(sdf)) {
      sdf.igraph[[i]] <- .f(sdf[[i]],addH  )
    }
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
  if (!length(highlight)) return(sdf.igraph.highlight)
  if (!is.numeric(highlight)) {
    highlight <- match(highlight  ,names(V(sdf.igraph.highlight)))
  }
  edge.data <- igraph::as_data_frame(sdf.igraph.highlight)%>%
    dplyr::mutate(highlight = from %in% V(sdf.igraph.highlight)$name[highlight] &
                    to %in% V(sdf.igraph.highlight)$name[highlight])

 # E(sdf.igraph.highlight)$width[edge.data$highlight] <-
 #   E(sdf.igraph.highlight)$width[edge.data$highlight]*2
 # E(sdf.igraph.highlight)$color[edge.data$highlight] <-
 #   E(sdf.igraph.highlight)$color[edge.data$highlight] %>%
 #   adjustcolor(blue.f = 1.7,green.f = 0.8,red.f = 0.3)
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

vis_sdf <- function(sdf,...){

  sdf.igraph <- get_sdf_igraph(sdf)
  vis_sdf_igraph(sdf.igraph,...)


}

vis_sdf_igraph_compare <- function(sdf.igraph1,
                                   sdf.igraph2,
                                   atom_id1=NULL,
                                   atom_id2=NULL,
                                   show.label= T){

  sdf.igraph1 <- add_sdf_igraph_highlight(sdf.igraph1,atom_id1)
  sdf.igraph2 <- add_sdf_igraph_highlight(sdf.igraph2,atom_id2)
  #hw1 <- diff(range(V(sdf.igraph1)$x))
  #hw2 <- diff(range(V(sdf.igraph2)$x))
  #V(sdf.igraph1)$x <- V(sdf.igraph1)$x-hw1/1.5
  #V(sdf.igraph2)$x <- V(sdf.igraph2)$x+hw2/1.5
  hw1 <- max(range(V(sdf.igraph1)$x))
  hw2 <- min(range(V(sdf.igraph2)$x))
  V(sdf.igraph1)$x <- -1-hw1+V(sdf.igraph1)$x
  V(sdf.igraph2)$x <- 1-hw2+V(sdf.igraph2)$x
  V(sdf.igraph1)$y <- V(sdf.igraph1)$y
  V(sdf.igraph2)$y <- V(sdf.igraph2)$y
  V(sdf.igraph1)$name <- paste0(V(sdf.igraph1)$name,"_1")
  V(sdf.igraph2)$name <- paste0(V(sdf.igraph2)$name,"_2")
  V(sdf.igraph1)$id <- paste0(V(sdf.igraph1)$id,"_1")
  V(sdf.igraph2)$id <- paste0(V(sdf.igraph2)$id,"_2")

  nodes <- rbind(
    vdata(sdf.igraph1),
    vdata(sdf.igraph2)
  )
  edges <- rbind(
    edata(sdf.igraph1),
    edata(sdf.igraph2)
  )
  new.igraph <- graph_from_data_frame(edges,
                                      vertices = nodes)
  vis_sdf_igraph(new.igraph,show.label = show.label)

}


vis_cfm_data_atom_map <- function(cfmd,fragment_id,
                                  show.label= T){

  maps <- cfmd@fragment_atom_map[[fragment_id]]
  maps <- maps[apply(maps, 1, function(x)any(x!=0)),
               apply(maps, 2, function(x)any(x!=0))
               ]
  vis_sdf_igraph_compare(
    cfmd@fragment_igraph[[1]],
    cfmd@fragment_igraph[[fragment_id]],
    atom_id1 = match( rownames(maps),names(V( cfmd@fragment_igraph[[1]]))),
    atom_id2 = match( colnames(maps),names(V( cfmd@fragment_igraph[[fragment_id]]))),
    show.label = T

  )

}


### to be removed
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

  ig.vd <- vdata(ig)
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


#' get_atom_map
#'
#' map atom of two molecular
#'
#'
#' @param sdf.parent sdf
#' @param sdf.product sdf
#' @param ig.parent ig
#' @param ig.product ig
#' @param return.type string
#'
#' @return if prob_matrix, return a matrix, col: atom of product, row: atom of parent
#' if most_prob, return the vector of most likely map
#' @export
#'
get_atom_map <- function(sdf.parent,
                         sdf.product,
                         ig.parent = get_sdf_igraph(sdf.parent),
                         ig.product = get_sdf_igraph(sdf.product),
                         return.type = c("prob_matrix","most_prob")){
  return.type <- match.arg(return.type)
  mcs <- fmcsR::fmcs(sdf.parent,sdf.product,bu = 10000)
  mcs.map <- get_mcs_atom_map(mcs)
  atom.map.matrix <- matrix(nrow = length(mcs.map),
                        ncol = length(atom(sdf.product)),
                        dimnames = list(seq_along(mcs.map),
                                        atom(sdf.product)))
  ring.diff <- length(rings(sdf.parent))- length(rings(sdf.product))
  bond.score <- rep(0,length(mcs.map))
  for (j in seq_along(mcs.map)) {
    this.map <- mcs.map[[j]]
    this.mapv <-this.map$mc1.atom
    names(this.mapv) <- this.map$mc2.atom
    this.mapv <- this.mapv[atom(sdf.product)]
    names(this.mapv) <- atom(sdf.product)


    ### ring re-assign
    {
      if (ring.diff) {
        ring.atom <- unname(unlist(rings(sdf.parent)))
        ring.atom.to.assign <- ring.atom[!ring.atom%in% this.mapv]
        ring.atom.to.assign <- unique(ring.atom.to.assign)
        adj <- sapply(ring.atom.to.assign,function(x){
          #x <- ring.atom.to.assign
          x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
          x.adj <- x.adj[x.adj%in%this.mapv]
          y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
          y.candi <-apply(distances(ig.product,y.adj),1,function(z){
            zz <- names(z)[which(z==1)]
             zz[!zz%in% names(na.omit(this.mapv))&
                  str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
            })%>%unlist()
          unname(y.candi[1])
          })

        adj <- na.omit(unlist(adj))
        this.mapv[adj] <- names(adj)

      }

    }

    ### bond diff
    {
      temp.map <- na.omit(this.mapv)
      m1 <- as_adj(ig.parent,
                   attr = "bond_type")[(temp.map),(temp.map)]
      m1 <- m1+t(m1)
      m2 <- as_adj(ig.product,
                   attr = "bond_type")[names(temp.map),names(temp.map)]
      m2 <- m2 + t(m2)
      bond.score[j] <- 1- sum((m1-m2)!=0)/sum(m1!=0)


    }
    #p<-vis_sdf_igraph_compare(ig.parent,ig.product,temp.map,names(temp.map),show.label = T)

    atom.map.matrix[j,] <- this.mapv


  }
  atom.count <- apply(atom.map.matrix,1,function(x)sum(!is.na(x)))
  selected <- which.max(atom.count+bond.score)
  if (return.type == "most_prob")
    return(atom.map.matrix[selected,])
  if (return.type== "prob_matrix"){
    apply(atom.map.matrix,2,function(s){
      sp <- table(s)/length(s)
      sp <- sp[atom(sdf.parent)]
      names(sp) <- atom(sdf.parent)
      sp[is.na(sp)] <- 0
      return(sp)
    })

  }

}


#' vis_smiles
#'
#' @param smiles smiles
#' @param show.formula logic
#' @param show.label logic
#' @param highlight v
#'
#' @return vis
#' @export
#' @import ChemmineR igraph visNetwork
vis_smiles <- function(smiles,
                       show.formula = T,
                       show.label =T,
                       highlight =NULL){

  smiles.sdf <- get_smiles_sdf(smiles)[[1]]
  smiles.igraph <- get_sdf_igraph(smiles.sdf)
  smiles.vis <- vis_sdf_igraph(smiles.igraph,
                               show.label = show.label,
                               highlight = highlight)

  if (show.formula) {
    smiles.vis$x$main<- list(text = unname(ChemmineR::MF(smiles.sdf,addH = T)),
                             style = "text-align:center")
  }
  smiles.vis

}


atom <- function(sdf){

  rownames(atomblock(sdf))


}
