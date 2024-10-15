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
  id.atom <- apply(atom.map.matrix,1,sum)>1

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

  smile.sdf <- get_smiles_sdf(smile)
  smile.formula <- get_sdf_formula(smile.sdf)
  smile.formula <- case_when(smile=="O"~"H2O1",
                 smile=="[HH]"~"H2",
                 T~smile.formula)

  return(smile.formula)
}


vis_sdf_igraph_old <- function(sdf.igraph ,
                               show_id = F,
                               prob.border = NULL,
                               prob.fill = NULL,
                               highlight = NULL){

  sdf.igraph.temp <-sdf.igraph

  ### map prob to color and hight
  {

    ele <-  get_sdf_igraph_atom(sdf.igraph)
    #prob.fill <- prob.border <- runif(10,0,1)%>%`names<-`(sample(get_sdf_igraph_atom(sdf.igraph),10))
    if (is.numeric(highlight)|is.logical(highlight)) highlight <- ele[highlight]
    prob.border[highlight] <- 1
    col.border <- .get_vis_col(sdf.igraph,prob.border,
                               colramp(breaks = c(0,Inf,1),
                                       colors = c("#aaaaaa","#97C2FC","#2B7CE9")))
    col.fill <- .get_vis_col(sdf.igraph,prob.fill,
                             na.col = "#DDDDDD",
                             colramp(breaks = c(0,Inf,1),
                                     colors = c("#FFFFFF","#F7844F","#B20C26")))
    ele <- get_sdf_igraph_atom(sdf.igraph)

  }

  vda <- vdata(sdf.igraph.temp)<- vdata(sdf.igraph.temp)%>%
    dplyr::mutate(label = case_when(show_id~id,
                                    T~paste0(" ",atom," ")),
                  label = str_format_len(label),
                  font.size = case_when(show_id~20,T~40),
                  # font.multi= T,
                  # font.bold = T,
                  # font.bold.mod = "bold",
                  # font.bold.size = 500,
                  font.vadjust = 5,
                  # font.strokeWidth = 2,
                  #  font.strokeColor = "black",
                  font.align = "left",
                  borderWidth = 3,
                  color.background = col.fill[name],
                  color.border = col.border[name],
                  shape = "circle"
    )

  sdf.igraph.temp%>%
    visIgraph(idToLabel = F,
              type = "square")%>%
    visEdges(arrows = list(to = F),
             length = 2)

}




get_sdf_igraph_atom <- function(ig,ele = "all"){

  vdf <- vdata(ig)
  if (ele== "all") {
    return(vdf$name)
  }else{
    vdf <- vdf %>%
      dplyr::filter(atom %in% ele)
    return(vdf$name)
  }

}



.get_highlight <- function(sdf.igraph,highlight){

}


vis_sdf <- function(sdf,show_id = F,...){

  sdf.igraph <- get_sdf_igraph(sdf)
  vis_sdf_igraph(sdf.igraph,show_id = show_id,...)


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
                         iso_ele = "[13]C",
                         return.type = c("most_prob","prob_matrix")){
  return.type <- match.arg(return.type)
  mcs <- fmcsR::fmcs(sdf.parent,sdf.product,bu = 10)
  mcs.map <- get_mcs_atom_map(mcs)
  mcs.map <- mcs.map.filter.duplicate(mcs.map,target_ele = get_ele_uniso(iso_ele))
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
      ring.solved <- F
      if (ring.diff&sum(is.na(this.mapv))) {
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
            })%>%unlist()%>%unique()
          #message("ring candi: ",length(y.candi))
          unname(y.candi[1])
          })

        adj <- na.omit(unlist(adj))
        this.mapv[adj] <- names(adj)
        ring.solved <- ifelse(length(adj),T,F)
      }

    }

    ### ring nearest
    {
      if (ring.diff&sum(is.na(this.mapv))&ring.solved){
#
        ring.nearest.to.assign <-apply(distances(ig.parent,ring.atom),
                        1,function(z){
          zz <- names(z)[which(z==1)]
          zz[!zz%in% (na.omit(this.mapv))]
        })%>%unlist()%>%unique()

        adj <- sapply(ring.nearest.to.assign,function(x){
          #x <- ring.nearest.to.assign
          x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
          x.adj <- x.adj[x.adj%in%this.mapv & x.adj%in%ring.atom]
          y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
          y.candi <-apply(distances(ig.product,y.adj),1,function(z){
            zz <- names(z)[which(z==1)]
            zz[!zz%in% names(na.omit(this.mapv))&
                 str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
          })%>%unlist()%>%unique()
          #message("ring nearest candi: ",length(y.candi))
          unname(y.candi[1])
        })

        adj <- na.omit(unlist(adj))
        this.mapv[adj] <- names(adj)


      }

    }

    ### non-match nearest
    {

     if (sum(is.na(this.mapv))) {
       non.match.to.assign <- setdiff(atom(sdf.parent),this.mapv)
       adj <- sapply(non.match.to.assign,function(x){
         #x <- non.match.to.assign
         x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
         x.adj <- x.adj[x.adj%in%this.mapv]
         y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
         y.candi <-apply(distances(ig.product,y.adj),1,function(z){
           zz <- names(z)[which(z==1)]
           zz[!zz%in% names(na.omit(this.mapv))&
                str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
         })%>%unlist()%>%unique()
         #message("non match nearest candi: ",length(y.candi))
         unname(y.candi[1])
       })

       adj <- na.omit(unlist(adj))
       this.mapv[adj] <- names(adj)

     }

    }

    ### bond diff
    {
      temp.map <- na.omit(this.mapv)
      temp.map.t <- make_vector(names(temp.map),temp.map)
      ig.sub <- igraph_filter_distance(ig.parent,from = temp.map,dis = 1)

      m1 <- as_adj(ig.sub,
                   attr = "bond_type")
      m1 <- m1+t(m1)
      m2 <- as_adj(ig.product,
                   attr = "bond_type")[names(temp.map),names(temp.map)]
      m2 <- m2 + t(m2)
      m2 <- get_matrix_value_fill_with_NA(
        m2,temp.map.t[rownames(m1)],temp.map.t[colnames(m1)])
      m2[is.na(m2)] <- 0
      bond.score[j] <- 1- sum((m1-m2)!=0)/sum(m1!=0)


    }
    #p<-vis_sdf_igraph_compare(ig.parent,ig.product,temp.map,names(temp.map),show.label = T)

    atom.map.matrix[j,] <- this.mapv


  }


  ### select map
  {
    atom.ele <- vdata(ig.product)$atom
    iso.atom.map.matrix <- atom.map.matrix[,atom.ele==get_ele_uniso(iso_ele ),drop = F]
    atom.count <- apply(iso.atom.map.matrix,1,function(x)sum(!is.na(x)))
    full.mapped <- apply(iso.atom.map.matrix,1,function(x)sum(is.na(x))==0)
    selected <- which((atom.count+bond.score)==max(atom.count+bond.score,na.rm = T))

  }
  if (return.type == "most_prob"){
    map <-  apply(atom.map.matrix[selected,,drop=F],2,function(x){
      x <- na.omit(x)
      xp <- table(x)/length(x)
      xp <- xp[atom(sdf.parent)]
      names(xp) <- atom(sdf.parent)
      xp[is.na(xp)] <- 0
      return(xp)
    })

    attributes(map)$bond.score <- mean(bond.score[selected])
  }


  if (return.type== "prob_matrix"){
    #apply(atom.map.matrix[atom.count==max(atom.count),,drop=F],2,function(x){
    map <-  apply(atom.map.matrix[full.mapped,,drop=F],2,function(x){
      x <- na.omit(x)
      xp <- table(x)/length(x)
      xp <- xp[atom(sdf.parent)]
      names(xp) <- atom(sdf.parent)
      xp[is.na(xp)] <- 0
      return(xp)
    })
    attributes(map)$bond.score <- mean(bond.score)

  }



  return(map)

}

mcs.map.filter.duplicate <- function(mcs.map,target_ele = "C"){


  mcs.list <- lapply(mcs.map,function(x){

    x <- x %>%
      dplyr::filter(
        grepl(target_ele,mc1.atom)
      )%>%
      dplyr::arrange(mc1.idx)

    make_vector(x$mc2.atom,name = x$mc1.atom)

  })

  return(mcs.map[!duplicated(mcs.list)])
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
                       show_id =T,
                       highlight =NULL){

  smiles.sdf <- get_smiles_sdf(smiles)[[1]]
  smiles.igraph <- get_sdf_igraph(smiles.sdf)
  smiles.vis <- vis_sdf_igraph(smiles.igraph,
                               show_id = show_id,
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


get_isopattern_score <- function(formula,
                                 mzs,
                                 int_matrix,
                                 ppm = 10){


  if (!length(formula)) return(NULL)
  formula.f <- factor(formula)
  iso_patterns <- lapply(levels(formula.f),
                        MSCC::chemform_isotopes_pattern_enviPat )
  iso_pattern <- iso_patterns[[1]]
  ip.score  <- lapply(iso_patterns,
         function(iso_pattern){
           if (nrow(iso_pattern)<=1) return(NA)
           iso_patterng <-iso_pattern %>%
             dplyr::ungroup()%>%
             dplyr::mutate(groupMz(x =m.z, ppm=ppm,return.type = "d"))%>%
             dplyr::group_by(mz.center)%>%
             dplyr::mutate(abundance=sum(abundance))%>%
             dplyr::distinct(mz.center,abundance)%>%
             dplyr::ungroup()

           id <- match_mz(mz1 = iso_patterng$mz.center,
                             mz2 = mzs,
                             mz.ppm  = ppm)
           iso.valm <- int_matrix[id,,drop = F]
           iso.ratio <- t(t(iso.valm)/iso.valm[1,])*100
           apply(iso.ratio , 2, function(iso.ratio.x){
             x <- iso_patterng$abundance[-1]
             y <- iso.ratio.x[-1]
             if (!length(x)) return(NA)
             if (all(is.na(y)))  return(0)
             y[is.na(y)] <- 0
             id.na <- is.na(x)|is.na(y)
             x <- x[!id.na]
             y <- y[!id.na]
             sum(x*y)^2/(sum(x^2L) *
                           sum(y^2L))
             1/exp(weighted.mean((abs(x-y)/x),w = x))
           })


         })

  ip.score <- sapply(ip.score,mean)
  ip.score <- ip.score[as.numeric(formula.f)]
  return(ip.score)
}



