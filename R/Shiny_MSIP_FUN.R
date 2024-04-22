
get_iso_cfm_compound_info <- function(iso.list){

  compound_info <-  lapply(iso.list, function(x){x$compound_info})%>%
    data.table::rbindlist()%>%
    dplyr::mutate(mz = format(mz,digits = 4),
                  rt =  format(rt,digits = 2),
                  score =format(score,digits = 2) )%>%
    dplyr::select(name,mz,rt,formula,adduct)

  return(compound_info)

}



shiny_plotly_iso_msip_spectra <- function(sp.data){

  ### if
  {
    if (all(is.na(sp.data))) {
      p <- shiny_plotly_empty(source = "plotly_ms2_sp")
      return(p)
    }
  }



  ### plot
  {

    ymax <- max(abs(sp.data$y))*1.2

    sp.data%>%
      plot_ly(source = "plotly_ms2_sp")%>%
      highlight_key(~fragment_group
                    )%>%
      add_segments(x = ~x , xend = ~x,
                   y = I(0),yend = ~y,
                   alpha = 0.5,
                   colors = c(" "= "#BBBBBB",
                              "CE=10"="#4DBBD5FF",
                              "CE=20" ="#FF7F0E",
                              "CE=40" ="#E64B35FF" ),
                   showlegend = F,
                   color = ~ collisionEnergy)%>%
      dplyr::filter(annotated)%>%
      add_markers(x = ~x,
                  y = ~y,
                  colors = c(" "= "#BBBBBB",
                             "10"="#4DBBD5FF",
                             "20" ="#FF7F0E",
                             "40" ="#E64B35FF" ),
                  color = ~ collisionEnergy,
                  text = ~ hover_label,
                  hovertemplate = "%{text}<extra></extra>",
                  hoverinfo = "text",
                  showlegend = T)%>%
      highlight(on = "plotly_click",
                off = NULL,
                selected = attrs_selected(
                  showlegend = F
                ),
                defaultValues=filter(sp.data,!is.na(fragment_group))$fragment_group[1])%>%
    #  dplyr::filter(fragment_group %in% fg)%>%
    #  add_segments(x = ~x , xend = ~x,
    #               y = I(0),yend = ~y,
    #               size = I(5),
    #               colors = c(" "= "#BBBBBB",
    #                          "10"="#4DBBD5FF",
    #                          "20" ="#FF7F0E",
    #                          "40" ="#E64B35FF" ),
    #               showlegend = F,
    #               color = ~ collisionEnergy)%>%
    #  add_markers(x = ~x,
    #              y = ~y,
    #              size = I(50),
    #              colors = c(" "= "#BBBBBB",
    #                         "10"="#4DBBD5FF",
    #                         "20" ="#FF7F0E",
    #                         "40" ="#E64B35FF" ),
    #              color = ~ collisionEnergy,
    #              text = ~ hover_label,
    #              hovertemplate = "%{text}<extra></extra>",
    #              hoverinfo = "text",
    #              showlegend = F)%>%
      layout(dragmode = "zoom",
             xaxis = list(range = c(0, max(sp.data$mz)*1.1),
                          title = "mz",
                          showgrid = T),
             yaxis = list(showgrid = T,
                          range = c(-ymax,ymax),
                          title = "Relative intensity"))%>%
      event_register("plotly_click")

  }
}


shiny_get_sp_data <- function(iso_msip,
                              iso_count = names(iso_msip$msip_data)[1]){


  ### if
  {
    if (!iso_count %in% names(iso_msip$Spectra)) {
      return(NA)
    }
  }

  ### sp data
  {
    sp.m0 <-iso_msip$Spectra$M0%>%
      normalizeSpectra()%>%
      combineSpectra_groupby_ce()
    sp.iso <- iso_msip$Spectra[[iso_count]]%>%
      normalizeSpectra()%>%
      combineSpectra_groupby_ce()

    sp.m0.frag.data <- CFM_annotate_isotopologues(sp.m0,
                                                  cfmd  = iso_msip$CFM_annotation,
                                                  ppm = 10,
                                                  iso.count = 0)
    sp.iso.frag.data <- CFM_annotate_isotopologues(sp.iso,
                                                   cfmd  = iso_msip$CFM_annotation,
                                                   ppm = 10,
                                                   iso.count = str_extract_num(iso_count))

    sp.m0.frag.data <- sp.m0.frag.data%>%
      dplyr::mutate(sp.id=1,x = mz,y = intensity)
    sp.iso.frag.data <- sp.iso.frag.data%>%
      dplyr::mutate(sp.id=2,
                    x = mz,
                    y = - intensity)

    sp.data <- rbind(sp.m0.frag.data,sp.iso.frag.data)%>%
      dplyr::mutate(annotated = !is.na(fragment_group),
                    collisionEnergy = as.character(collisionEnergy),
                    collisionEnergy = case_when(annotated~paste0("CE=",collisionEnergy),
                                                T~ " " ),
                    hover_label = paste0("Frag group: ",fragment_group,"\n",
                                         "mz: ",round(mz,4),"\n",
                                         "intensity: ",round(intensity,4)),
                    point.size = case_when(annotated~0.5,
                                           T~0))

  }

  return(sp.data)


}


shiny_plotly_empty <- function(...){
  plot_ly(...)%>%
    add_markers(x = 0,y = 0,color = I("transparent"))%>%
    layout(xaxis = list(range = c(0,100)),
           yaxis = list(range = c(-100,100)))
}

shiny_get_fg <- function(sp.data,
                         x){

  if (all(is.na(sp.data))) {
    return(NA)
  }
  if (length(x)) {
    idx <- match_mz(x,sp.data$x,mz.ppm = 5)
    return(sp.data$fragment_group[idx])
  }
  return(na.omit(sp.data$fragment_group)[1])
}

shiny_get_fid <- function(iso_msip,
                          fg){

  if (all(is.na(iso_msip))) {
    return(NA)
  }
  if (length(fg)) {
    fid <-  iso_msip$CFM_annotation@fragment_define%>%
      dplyr::filter(fragment_group ==fg)%>%
      dplyr::pull(fragment_id )
    return(fid)
  }
  return(NA)
}

shiny_get_fg_ig <- function(iso_msip,fid = 1){

  if (all(is.na(iso_msip))) {
    return(NA)
  }
  if (all(is.na(fid))) {
    return(NA)
  }
  iso_msip$CFM_annotation@fragment_igraph[[fid]]

}

shiny_get_frag_formula <- function(iso_msip,fid){

  frag.def <- iso_msip$CFM_annotation@fragment_define
  if (fid %in% frag.def$fragment_id) {
    formula <- frag.def%>%
      dplyr::filter(fragment_id == fid)%>%
      dplyr::pull(formula)
  }else{
    formula <- " "
  }
  return(formula)

}

shiny_get_atom_map <- function(iso_msip,
                               fid,
                               prob = T){
  x <- list(NA,NA)
  if (is.null(iso_msip)) {
    return(x)
  }
  if (all(is.na(iso_msip))) {
    return(x)
  }

  if (!fid%in%names( iso_msip$CFM_annotation@fragment_atom_map) ) {
    return(x)
  }

  map.matrix <- iso_msip$CFM_annotation@fragment_atom_map[[fid]]
  if (nrow(map.matrix)) {

    if (prob) {
      atom1 <- apply(map.matrix,1,sum)
      atom2 <- apply(map.matrix,2,sum)
    }else{
      atom1 <- apply(map.matrix,1,sum)
      atom1 <- names(atom1)[atom1!=0]
      atom2 <- apply(map.matrix,2,sum)
      atom2 <- names(atom2)[atom2!=0]

    }
    atom.map <- list(
      atom1,atom2
    )
    return(atom.map)
  }
  return(x)
}


shiny_col_ramp_atom_prob <- function(){

  colramp(breaks = c(0,0.5,1),
          colors = c("white","#F7844F","#B20C26")
  )
}

shiny_vis_ig <- function(ig,
                         prob.border= NA,
                         prob.fill =NA){



  ### if
  {
    if (class(ig)!="igraph"|is.null(ig)) {
      return(visIgraph(make_empty_graph()))
    }
    if(is.null(prob.border))
      prob.border <-0
    if(is.null(prob.fill))
      prob.fill <-0
  }

  ### add color
  {

    atom <- get_atom_from_igraph(ig)
    prob.border <- prob.border[atom]
    names(prob.border) <- atom
    prob.border[is.na(prob.border)] <- 0
    prob.fill <- prob.fill[atom]
    names(prob.fill) <- atom
    prob.fill[is.na(prob.fill)] <- 0

    atom.border.color <- colramp(breaks = c(0,0.00001,0.5,1),
                          colors = c("grey","white","#97C2FC","#2B7CE9")
    )(prob.border)
    atom.fill.color <- shiny_col_ramp_atom_prob()(prob.fill)


    vdata(ig)$color.border <- atom.border.color
    vdata(ig)$color.background  <- atom.fill.color
    vdata(ig)$shape <- "circle"
  }

  ### plot
  {
    ig%>%
      visIgraph(idToLabel = F)%>%
      visNodes(font = list(size = 40,
                           align = "left",
                           vadjust = 3,
                           hadjust = 0.8,
                           strokeWidth = 2),
               size = 40,
               borderWidth  = 5)%>%
      visEdges(arrows = list(to = F),
               length = 0.8)

  }

}

shiny_change_ig_label <- function(ig,
                                  id = F
                                  ){


  if (id) {
    vdata(ig)$label <-
      vdata(ig)$id

  }else{
    vdata(ig)$label <-
      vdata(ig)$atom

  }

  vdata(ig) %>%
    dplyr::select(id,label)


}

