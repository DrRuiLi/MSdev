get_nodes_between_selected <- function(ig,selected.node){

  dism <- distances(ig)[,selected.node]
  node.con <- apply(dism,1,function(x){
    sum(x<2&x>0)
  })
  nodes <- node.con[node.con>=2]%>%names()
  c(nodes,selected.node)%>%unique()%>%return()
}

get_edges_from_path <- function(ig,v){

  vp <- rep(names(v),each=2)
  vp <- vp[-c(1,length(vp))]
  E(ig)[get.edge.ids(ig, vp)]

}

igraph_filter_vertex <- function(ig,v){

  if (is.numeric(v)|is.logical(v)|is.character(v))
    v <- V(ig)[v]
  delete.vertices(ig,setdiff(V(ig),v))

}

igraph_filter_shortest_path <- function(ig,from,to){

  sp <- shortest_paths(ig,from,to)
  ids <- sapply(sp$vpath,as_ids)
  ids <- unique(unlist(ids))
  igraph_filter_vertex(ig,ids)
}


show_vis_icon <- function(icon_code = paste0("f",num2str(1:900)),
                          type = c("FontAwesome","Ionicons")){

  type <- match.arg(type)
  n.row <- sqrt(length(icon_code))%>%ceiling()

  vis.v <- vertex(icon_code)
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



edata <- function(ig){

  igraph::as_data_frame(ig,"e")

}

vdata <-  function(ig){

  igraph::as_data_frame(ig,"v")

}

`edata<-` <- function(ig,value){

  edge.attributes(ig) <- as.list(value)
  ig
}

`vdata<-` <- function(ig,value){

  vertex.attributes(ig) <- as.list(value)
  ig
}

open_visNet <- function(x){

  tpf <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(x,file = tpf)
  open_file(tpf)

}


