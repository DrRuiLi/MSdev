fella_enrich <- function(cpd_kegg,
                         method = c("hypergeom", "diffusion", "pagerank" ) ){

  fella.data <<- loadKEGGdata(internalDir = T,loadMatrix = c("diffusion", "pagerank"))
  fella.fella <- enrich(cpd_kegg,
                        data = fella.data,
                        methods  = method)


  return(fella.fella)
}

fella_igraph <- function(fella.fella){

  fella.graph <- generateResultsGraph(fella.fella,
                                      method = "diffusion",
                                      plimit = 5,
                                      nlimit = 1000 ,
                                      threshold = 0.05,
                                      thresholdConnectedComponent = 0.05,
                                      LabelLengthAtPlot = 100,
                                      data = fella.data)
  com <- V(fella.graph)$com
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
