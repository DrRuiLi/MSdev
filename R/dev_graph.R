get_nodes_between_selected <- function(ig, selected.node){

  dism <- distances(ig)[,selected.node]
  node.con <- apply(dism,1,function(x){
    sum(x<2&x>0)
  })
  nodes <- node.con[node.con>=2]%>%names()
  c(nodes,selected.node)%>%unique()%>%return()
}

get_edges_from_epath <- function(ig,v,directed= F){

  vp <- rep(names(v),each=2)
  vp <- vp[-c(1,length(vp))]
  E(ig)[get.edge.ids(ig, vp,directed =directed )]

}

igraph_filter_vertex <- function(ig,v){

  if (is.numeric(v)|is.logical(v)|is.character(v))
    v <- igraph::V(ig)[v]
  igraph::delete.vertices(ig,setdiff(names(V(ig)),names(v)))

}

igraph_filter_edge <- function(ig,e){

  if (is.numeric(e))
    e <- igraph::E(ig)[e]


  ig <- igraph::delete.edges(ig,setdiff(igraph::E(ig),e))
  igraph_filter_vertex(ig,igraph::degree(ig)!=0)


}

igraph_remove_edge <- function(ig,e){

  if (is.numeric(e))
    e <- igraph::E(ig)[e]


  ig <- igraph::delete.edges(ig,e)
  #igraph_filter_vertex(ig,igraph::degree(ig)!=0)
  return(ig)

}

igraph_remove_vertex <- function(ig,v){

  if (is.numeric(v)|is.logical(v)|is.character(v))
    v <- igraph::V(ig)[v]
  igraph::delete.vertices(ig,v)

}

igraph_filter_distance <- function(ig, from , dis = 1,...){

  dis.matrix <- distances(ig,from)
  dis.pass <- dis.matrix <= dis
  id <- apply(dis.pass,2,function(x){
    any(x)
  })

  igraph_filter_vertex(ig,id)

}


igraph_filter_shortest_path <- function(ig,from,to){

  sp <- shortest_paths(ig,from,to)
  ids <- sapply(sp$vpath,as_ids)
  ids <- unique(unlist(ids))
  igraph_filter_vertex(ig,ids)
}


igraph_filter_path<- function(ig,paths){


  ids <- sapply(paths,igraph::as_ids)
  ids <- unique(unlist(ids))
  igraph_filter_vertex(ig,ids)
}

get_path_direction <- function(ig,vpath,epath){

  vs <- names(vpath)
  vs <- c(vs,vs[1])
  es <- igraph::ends(ig,epath)

  ep.dir <- rep(1,length(epath))
  for (i in 1:nrow(es)) {
    ep.dir[i] <- ifelse( es[i,1] == vs[i] &es[i,2] == vs[i+1],1,-1)
  }
  return(ep.dir)

}


igraph_add_reverse_edges <- function(ig){

  eda <- edata(ig)%>%
    dplyr::mutate(direction = 1)
  edata(ig) <- eda
  eda.rev <- eda%>%
    dplyr::mutate(
      tmp = from,
      from = to,
      to = tmp,
      direction = - direction
    )%>%
    dplyr::select(-tmp,-from,-to)

  ig <- igraph::add_edges(ig,
                          as.vector(rbind(eda$to,eda$from)),
                           attr = eda.rev)

  ig
}

show_vis_icon <- function(icon_code = paste0("f",num2str(1:900)),
                          type = c("FontAwesome","Ionicons")){

  type <- match.arg(type)
  n.row <- sqrt(length(icon_code))%>%ceiling()

  vis.v <- igraph::vertex(icon_code)
  vis.v$x <- rep(x = 1:n.row,
               times = n.row)[1:length(icon_code)]
  vis.v$y <- rep(x = 1:n.row,
               each = n.row)[1:length(icon_code)]
  vis.v$shape = "icon"
  if (type=="Ionicons")
    vis.v$icon.face = 'Ionicons'
  vis.v$icon.code  = icon_code
  vis.v$code  = icon_code
  vis.v$label  =length(icon_code)
  vis <-  igraph::make_empty_graph()+ vis.v
  vis <- vis+edge(c(1,2))
  vis <- visIgraph(vis)

  if (type =="Ionicons" ) {
    return(addIonicons(vis))
  }
  if(type =="FontAwesome"){
    return(addFontAwesome(vis))
    }


}





open_visNet <- function(x){

  tpf <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(x,file = tpf)
  open_file(tpf)

}


igraph_vpath_to_epath <- function(ig,vpath){

  epaths <- lapply(vpath, function(path) {
    edges <- c()
    for (i in seq_along(path)[-length(path)]) {
      edges <- c(edges, igraph::get_edge_ids(ig, c(path[i], path[i + 1]),
                                     directed  =F))
    }
    edges
  })

  return(epaths)
}

igraph_sort_direction <- function(ig){

  eda <- edata(ig)
  rev.id <- which(eda$direction==-1)
  igd <- igraph::delete_edges(ig,rev.id)
  eda.to.add <- edata(ig)%>%
    dplyr::filter(direction == -1)%>%
    dplyr::mutate(tmp = from,
                  from = to,
                  to = tmp,
                  direction = -direction)
  eda.attr <- eda.to.add%>%
    dplyr::select(-from,-to,-tmp)

  igd <- igraph::add_edges(igd,
                          as.vector(rbind(eda.to.add$to,eda.to.add$from)),
                          attr =eda.attr)
  return(igd)
}



vis_add_text <- function(vis,text, font.size =10 ){

  to.add.df <- data.frame(
    name = "text",
    shape = "box",
    font.size = font.size,
    label = text,
    color.border = "#888888",
    color.background = "#eeeeee"

  )

  vis$x$nodes <- bind_rows(vis$x$nodes,to.add.df)
  vis
}

vis_add_image <- function(vis,img.file ,size = 500){

  img.txt <- RCurl::base64Encode(readBin(img.file,
                                     'raw',
                                     file.info(img.file)[1, 'size']),
                             'txt')

  to.add.df <- data.frame(
    x = 0,y = 0,
    name = "image",
    size= size,
    shape = 'image',
    color.border = "rgba(256,256,256,1)",
    image = paste('data:image/png;base64', img.txt, sep = ','),
    stringsAsFactors = F
  )

  vis$x$nodes <- bind_rows(to.add.df,vis$x$nodes)
  vis

}


vis_add_arrow_icon <- function(vis,
                               size = 60,
                               color = "#222222"){

  to.add.df <- nodes <- data.frame(x = 0,
                                   y = 0,
                                   name = "ARROW",
                                   shape = "icon",  # Specify that the node shape is an icon
                                   icon.face = 'FontAwesome',  # Use FontAwesome icons
                                   icon.code = "f178",  # Unicode for an arrow (FontAwesome)
                                   icon.size = size,  # Set the size of the icon
                                   icon.color =color)  # Set the color of the icon


  vis$x$nodes <- bind_rows(vis$x$nodes,to.add.df)
  vis%>%
  addFontAwesome(version  = "4.7.0")
}


vis_igraph <- function(ig){


  vda <- vdata(ig)%>%
    dplyr::mutate(id = name, .before = name)
  visNetwork::visNetwork(nodes = vda,edges = edata(ig))%>%
    visEdges(arrows = "to")%>%
    visNetwork::visOptions(width = "100%",
                           height = "100%")

}


vis_pave_igraph <- function(ig){

  vis_igraph(ig)%>%
    visNetwork::visEdges(
      smooth = F,
      length = 300,
      font = list(
       align = "top"
      ))
}

igraph_add_vcolor<- function(ig,v,color){

  vda <- vdata(ig)
  if (!"color.border" %in% colnames(vda)) {
    vda[v,"color.border"] <- NA
  }
  vda[v,"color.border"] <- color
  vda -> vdata(ig)
  return(ig)


}

igraph_add_vfill<- function(ig,v,color){

  vda <- vdata(ig)
  if (!"color.background" %in% colnames(vda)) {
    vda[v,"color.background"] <- NA
  }
  vda[v,"color.background"] <- color
  vda -> vdata(ig)
  return(ig)
}

igraph_add_ecolor<- function(ig,e,color){

  eda <- edata(ig)
  if (!"color" %in% colnames(eda)) {
    eda[e,"color"] <- NA
  }
  eda[e,"color"] <- color
  eda -> edata(ig)
  return(ig)
}

igraph_add_earrow<- function(ig,e,arrows = "to"){

  eda <- edata(ig)
  if (!"arrows" %in% colnames(eda)) {
    eda[e,"arrows"] <- NA
  }
  eda[e,"arrows"] <- arrows
  eda -> edata(ig)
  return(ig)
}



igraph_export <- function(ig,file){

  ig.list <- list(
    edata = edata(ig),
    vdata = vdata(ig)
  )

  xlsx.write.list(ig.list,file )
}

igraph_import <- function(file){

  vdata <- readxl::read_excel(file,sheet = "vdata")
  edata <- readxl::read_excel(file,sheet = "edata")
  igraph::graph.data.frame(edata,vertices = vdata)
}


igraph_get_nodes_distance <- function(ig,from,dis){

  dist <- distances(ig,from)

  id.dis <- which(apply(dist,2,function(x){any(x<=dis)}))
  setdiff(names(id.dis),from)

}


get_igraph_membership <- function(ig){


  node.group <- igraph::components(ig)$membership
  return(node.group)
}
