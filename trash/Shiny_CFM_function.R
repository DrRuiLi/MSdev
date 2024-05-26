ggplot_cfm_spectra <- function(peak_data){

  peak_data <- cfm.data$peak.assignment



}

get_cfm_data_shiny_formate <- function(cfm.data = read_CFM_annotate_result()){


  cfm.data <- get_cfm_data_igraph(cfm.data)
  ### data preprocess
  nodes <- cfm.data$fragment_define2%>%
    dplyr::mutate(id = fragment_id,
                  checked = check_smile(smile),
                  formula = get_smile_formula(smile),
                  label = formula
    )%>%
    dplyr::select(-root_atom_id )
  nodes.sdf <- smiles2sdf(nodes$smile%>%
                            `names<-`(nodes$fragment_id))
  edges <- cfm.data$fragment_transition%>%
    dplyr::mutate(arrows = "to",
                  formula = get_smile_formula(smile))

  cfm.data$fragment_sdf <- nodes.sdf
  cfm.data$transition_graph <- graph_from_data_frame(edges,
                                                     vertices = nodes)

  ### tracing stat
  tracing_stat <- cfm.data$peak.assignment%>%
    dplyr::group_by(mz,enery)%>%
    dplyr::mutate(peak_id = cur_group_id(),
                  root_atom = cfm.data$fragment_define2$root_atom_id[fragment_id])%>%
    dplyr::ungroup()%>%
    dplyr::mutate(peak_id = paste0("peak",num2str(peak_id)))%>%
    dplyr::arrange(mz)%>%
    dplyr::select(mz,intensity,
                  fragment_id,
                  fragment_ratio = fragment_score,
                  root_atom)%>%
    dplyr::filter(!is.na(fragment_id))
  cfm.data$tracing_stat <- tracing_stat

  return(cfm.data)

}

plotly_cfm_spectra <- function(peak.data,
                               mz.highlight =NULL){

  peak.data%>%
    dplyr::mutate(annotated = !is.na(fragment_id))%>%
    plot_ly(source = "sp_plot")%>%
    add_segments(x = ~mz , xend = ~mz,
                 y = I(0),yend = ~intensity,
                 #color = I("grey"),
                 color = ~ annotated,
                 colors = c("#EEEEEE","#999999"))%>%
    add_markers(x = ~mz,y = ~intensity,
                color = I("#FF7F0E"),
                hovertemplate = "mz:%{x}\nintensity:%{y}<extra></extra>",
                hoverinfo = "text",
                showlegend = F)%>%
    dplyr::filter(mz %in% mz.highlight)%>%
    add_segments(x = ~mz , xend = ~mz,
                 y = I(0),yend = ~intensity,
                 color = I("#5384CB"))%>%
    add_markers(x = ~mz,y = ~intensity,
                size = 3,color = I("#FF7F0E"),
                alpha = 1,
                #hovertemplate = "mz:%{x}\nintensity:%{y}<extra></extra>",
                hoverinfo = "text",
                showlegend = F)%>%
    layout(dragmode = "zoom",
           showlegend = F,
           xaxis = list(range = c(0, max(peak.data$mz)*1.1),
                        showgrid = FALSE),
           yaxis = list(showgrid = FALSE))%>%
    event_register('plotly_click')

}

vis_fragment_trans <- function(trans.graph ,
                               highligted = NULL,
                               frag.selected = NULL){
  highligted <- c(highligted,frag.selected)
  V(trans.graph)$color <- "#D9D9D9"
  V(trans.graph)[id%in%highligted ]$color <- "#97C2FC"
  vi <- visIgraph(trans.graph)%>%
    visNodes()
  if (is.null(frag.selected)|length(frag.selected)==0) {
    vi <- vi%>%visOptions(
      highlightNearest = list(
        enabled = F,
        degree = 0
      ),
      manipulation = TRUE,
      nodesIdSelection = list(
        enabled=T,
        useLabels = F
      ))
  }else{
    vi <- vi%>%visOptions(
      highlightNearest = list(
        enabled = F,
        degree = 0
      ),
      manipulation = TRUE,
      nodesIdSelection = list(
        enabled=T,
        selected = frag.selected,
        useLabels = F
      ))
  }
  vi

}
#' @import visNetwork MSCC DT
get_cfm_shiny_app <- function(shiny_sever_data){

  ### shiny
  ui <- get_cfm_shiny_ui()
  sever <- get_cfm_shiny_server()
  shinyApp(ui,sever,options = list(host = "0.0.0.0",
                                   port = 9187))

}




