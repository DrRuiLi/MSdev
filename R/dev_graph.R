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


  ids <- sapply(paths,as_ids)
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

  igraph::as_data_frame(ig,"edges")

}

vdata <-  function(ig){

  igraph::as_data_frame(ig,"vertices")

}

`edata<-` <- function(ig,value){

  value <- value[,!grepl("from|to",colnames(value))]
  edge.attributes(ig) <- as.list( value )
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


igraph_node_path_to_edge_path <- function(graph, node_path) {
  edge_paths <- lapply(node_path, function(path) {
    edges <- unlist(lapply(seq_along(path)[-length(path)], function(i) {
      E(graph, path = c(path[i], path[i + 1]))
    }))
    return(edges)
  })
  return(edge_paths)
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

  visNetwork::visNetwork(nodes = vdata(ig),edges = edata(ig))%>%
    visNetwork::visOptions(width = "200%",
                           height = "200%")

}
