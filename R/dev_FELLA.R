fella_enrich <- function(cpd_kegg,
                         method = c("hypergeom", "diffusion", "pagerank" ) ){

  fella.data <<- loadKEGGdata(internalDir = T,loadMatrix = c("diffusion", "pagerank"))
  fella.fella <- enrich(cpd_kegg,
                        data = fella.data,
                        methods  = method)


  return(fella.fella)
}

#' Title
#'
#' @param fella.fella
#'
#' @return
#' @export
#' @import igraph
#'
#' @examples
fella_igraph <- function(fella.fella,
                         p = 0.05,
                         node = 1000){

  fella.graph <- generateResultsGraph(fella.fella,
                                      method = "diffusion",
                                      plimit = 5,
                                      nlimit = node ,
                                      threshold = p,
                                      thresholdConnectedComponent = 0.05,
                                      LabelLengthAtPlot = 100,
                                      data = fella.data)
  fella.pscore <- FELLA::getPscores(fella.fella,
                                    method = "diffusion")
  com <- V(fella.graph)$com
  v.name <- V(fella.graph)$name
  V(fella.graph)$p.enrich <- fella.pscore[v.name]
  cf <-  circlize::colorRamp2(breaks = c(0,1.30103,5),
                       colors = c("white","white","#D21E31"))
  V(fella.graph)$col.enrich <- cf(-log10(fella.pscore[v.name]))%>%str_sub(.,1,7)
  V(fella.graph)$type <- dplyr::case_when(
    com==1 ~"Pathway",
    com==2~"Module",
    com==3~"Enzyme",
    com==4~"Reaction",
    com==5~"Metabolite")
  V(fella.graph)$type.color <-dplyr::case_when(
    com==1 ~"#CD0000",
    com==2~"#CD96CD",
    com==3~"#FFA200",
    com==4~"#8DB6CD",
    com==5~"#548B54")

  return(fella.graph)


}
