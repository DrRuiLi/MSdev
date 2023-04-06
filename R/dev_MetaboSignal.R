#' @title MetaboSignalNetworkTable2df
#'
#' @param network_table matrix with 3 column: source, target, interaction_type
#'
#'
#' @return
#' @export
#'
#' @examples
MetaboSignalNetworkTable2df <- function (network_table)
{
  network = unique(network_table)
  rev_ind = grep("compound:reversible", network[, 3])

  if (length(rev_ind) > 0) {
    rev_edges = paste(network[rev_ind, 1], network[rev_ind,
                                                   2], network[rev_ind, 3], sep = ">")
    new_edges = do.call(rbind, lapply(rev_edges, MetaboSignal:::cyto_directionality))
    network = rbind(network[-rev_ind, ], unique(new_edges))
    rownames(network) = NULL
  }


  cytoscape_net = network

  all_net_edges = paste(cytoscape_net[, 1], cytoscape_net[,
                                                          2], cytoscape_net[, 3], sep = "_")

  sources = as.character(sapply(all_net_edges, MetaboSignal:::get_source))
  cytoscape_netDF = unique(as.data.frame(cbind(cytoscape_net,
                                               database = sources), rownames = NULL))

  all_nodes = unique(as.vector(network[, 1:2]))
  all_nodes_info <- get_node_info(all_nodes)

  net.data <- list(net.df =cytoscape_netDF,
                   node.df = all_nodes_info)
  return(net.data)
}


#' @title MetaboSignalNetworkTableStandardNodeName
#' @description add or remove "hsa:" from network_table
#'
#' @param network_table
#' @param add_hsa
#'
#' @return
#' @export
#'
#' @examples
MetaboSignalNetworkTableStandardNodeName <- function(network_table,add_hsa = T){

  if (add_hsa) {
    network_table <- network_table%>%
      as.data.frame()%>%
      dplyr::mutate(source = case_when(!is.na(as.numeric(source)) ~ paste0("hsa:",source),
                                       T~source),
                    target = case_when(!is.na(as.numeric(target))  ~ paste0("hsa:",target),
                                       T ~ target))%>%
      as.matrix()

  }else{
    network_table <- network_table%>%
      as.data.frame()%>%
      dplyr::mutate(source = case_when(grepl(pattern = "hsa",x = source )~ gsub(pattern = "hsa:",x = source ,replacement = ""),
                                       T~source),
                    target = case_when(grepl(pattern = "hsa",x = target )~ gsub(pattern = "hsa:",x = target ,replacement = ""),
                                       T ~ target))%>%
      as.matrix()




  }
  network_table

}


#' @title get_node_info
#'
#' @param all_nodes kegg id, currently support cpd(compound) and hsa(gene)
#'
#' @return
#' @export
#'
#' @examples
get_node_info <- function(all_nodes){

  if (!exists("KEGG.database")) {

    MSdb::load_KEGG_database()
  }

  fff <- function(node){
    node.name <-NA
    node.type <- "other"

    if (!is.na(as.numeric(node))) {
      node <- paste0("hsa:",node)
    }



    if (grepl(pattern = "cpd",node)) {
      node.data <- KEGG.database$compound.list$data[[sub(pattern = "cpd:",x = node,replacement = "")]]
      node.name = unlist(strsplit(node.data$NAME, "[;]"))[1]
      node.type <- "compound"


    }

    if (grepl(pattern = "hsa",node)) {
      node.data <- KEGG.database$gene.list$data[[sub(pattern = "hsa:",x = node,replacement = "")]]
      symbol <- node.data$SYMBOL
      if (is.null(symbol)) {
        node.name <- NA
      }else{

        name = unlist(strsplit(symbol, "[,]"))
        name = unlist(strsplit(name, "[;]"))
        if (grepl("E", name[1]) == TRUE) {
          name = sort(name, decreasing = FALSE)
        }
        name = name[1]
        node.name = gsub(" ", "", name)
      }

      lines <- as.character(node.data)

        enzyme_lines = grep("EC:", lines[1:5])
        metabo_lines = grep("Metabolism", lines)
        if (length(enzyme_lines) >= 1 & length(metabo_lines) >0) {
          node_type = "metabolic-gene"
        }else if (sub(pattern = "hsa:",x = node,replacement = "") %in% MetaboSignal::regulatory_interactions[, c(1, 3)]) {
          node_type = "signaling-gene"
        }
        else {
          node_type = "other"
        }
        node.type <- node_type

    }

    return(list(node = node ,
                node.name = node.name,
                node.type = node.type))




  }

 #if (length(all_nodes) > 500) {

 #  bpp <- BiocParallel::SnowParam(progressbar = T)

 #}else{

 #  bpp <- BiocParallel::SerialParam(progressbar = T)
 #}
  bpp <- BiocParallel::SerialParam(progressbar = T)
  node.df  <- BiocParallel::bplapply(all_nodes,fff,
                                     BPPARAM = bpp)%>%
    data.table::rbindlist()


return(node.df)




}




