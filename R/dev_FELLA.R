fella_enrich <- function(cpd_kegg,
                         method = c("hypergeom", "diffusion", "pagerank")) {
  fella.data <<- FELLA::loadKEGGdata(internalDir = T,
                              loadMatrix = c("diffusion", "pagerank"))
  fella.fella <- FELLA::enrich(cpd_kegg, data = fella.data, methods  = method)


  return(fella.fella)
}

#' fella_igraph
#'
#' @param fella.fella
#'
#' @return igraph
#' @export
#'

fella_igraph <- function(fella.fella,
                         p = 0.05,
                         node = 1000) {
  fella.graph <- FELLA::generateResultsGraph(
    fella.fella,
    method = "diffusion",
    plimit = 5,
    nlimit = node ,
    threshold = p,
    thresholdConnectedComponent = 0.05,
    LabelLengthAtPlot = 100,
    data = fella.data
  )
  fella.pscore <- FELLA::getPscores(fella.fella, method = "diffusion")
  com <- V(fella.graph)$com
  v.name <- V(fella.graph)$name
  V(fella.graph)$p.enrich <- fella.pscore[v.name]
  cf <-  circlize::colorRamp2(breaks = c(0, 1.30103, 5),
                              colors = c("white", "white", "#D21E31"))
  V(fella.graph)$col.enrich <- cf(-log10(fella.pscore[v.name])) %>% str_sub(., 1, 7)
  V(fella.graph)$type <- dplyr::case_when(
    com == 1 ~ "Pathway",
    com == 2 ~ "Module",
    com == 3 ~ "Enzyme",
    com == 4 ~ "Reaction",
    com == 5 ~ "Metabolite"
  )
  V(fella.graph)$type.color <- dplyr::case_when(
    com == 1 ~ "#CD0000",
    com == 2 ~ "#CD96CD",
    com == 3 ~ "#FFA200",
    com == 4 ~ "#8DB6CD",
    com == 5 ~ "#548B54"
  )

  return(fella.graph)


}



fella_get_igraph_for_vis <- function(fella.fella,
                                     p = 0.05,
                                     node = 1000) {
  fella.igraph <- FELLA::generateResultsGraph(
    fella.fella,
    method = "diffusion",
    plimit = 5,
    nlimit = node ,
    threshold = p,
    thresholdConnectedComponent = 0.05,
    LabelLengthAtPlot = 100,
    data = fella.data
  )
  fella.pscore <- FELLA::getPscores(fella.fella, method = "diffusion")
  cf <-  circlize::colorRamp2(breaks = c(0, 1.30103, 5),
                              colors = c("white", "white", "#D21E31"))
  vdata(fella.igraph) <- vdata(fella.igraph) %>%
    dplyr::mutate(
      id = name,
      p.enrich = fella.pscore[name],
      col.enrich = cf(-log10(p.enrich)) %>%
        str_sub(., 1, 7),
      type = case_when(
        com == 1 ~ "Pathway",
        com == 2 ~ "Module",
        com == 3 ~ "Enzyme",
        com == 4 ~ "Reaction",
        com == 5 ~ "Metabolite"
      ),
      label = case_when(
        com == 1 ~ label,
        com == 2 ~ "",
        com == 3 ~ "",
        com == 4 ~ "",
        com == 5 ~ label
      ),
      size = case_when(
        com == 1 ~ 50,
        com == 2 ~ 40,
        com == 3 ~ 20,
        com == 4 ~ 20,
        com == 5 ~ 30
      ),
     #shape = case_when(
     #  com == 1 ~ "circle",
     #  com == 2 ~ "circle",
     #  com == 3 ~ "square",
     #  com == 4 ~ "square",
     #  com == 5 ~ "square"
     #),
      shape = case_when(
        input ~ "square",
        T~ "dot"
      ),
      color.background = case_when(
        com == 1 ~ "#CD0000",
        com == 2 ~ "#CD96CD",
        com == 3 ~ "#FFA200",
        com == 4 ~ "#8DB6CD",
        com == 5 ~ "#548B54"
      ),
     font.color = color.background,
     font.size = 20
    )


  return(fella.igraph)


}

