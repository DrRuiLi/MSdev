
get_MSIP_compound_info <- function(object,
                                      vars =c("name","mz","rt","formula","adduct")){

  avalible.vars <- c("name","compound_id","mz","rt","formula","smiles","score","adduct","polarity","merged")

  if ("all" %in% vars) {
    vars <- avalible.vars
    }
  iso.list <- object@statData$MSIP$isotopologues_data

  suppressWarnings(
  compound_info <-  lapply(iso.list, function(x){
    y <- x$compound_info[avalible.vars]
    names(y) <- avalible.vars
    y
    })%>%
    data.table::rbindlist(fill = T)%>%
    dplyr::mutate(mz = format(mz,digits = 4),
                  rt =  format(rt,digits = 2),
                  score =format(score,digits = 2) )%>%
    dplyr::select(all_of(vars))
  )

  return(compound_info)

}



shiny_plotly_iso_data_spectra <- function(sp.data,
                                          show.rawdata = T){

  ### if
  {
    if (all(is.na(sp.data))) {
      p <- shiny_plotly_empty(source = "plotly_ms2_sp")
      return(p)
    }

    if (show.rawdata) {
      sp.data <- sp.data

    }else{
      sp.data <- sp.data %>%
        dplyr::filter(merged)%>%
        dplyr::group_by(sp.id)%>%
        dplyr::mutate(y = 100*y/max(abs(y)))%>%
        dplyr::ungroup()
    }
  }



  ### plot
  {

    ymax <- max(abs(sp.data$y))*1.2
    ymin <- min((sp.data$y))*1.2

    p <- sp.data%>%
      dplyr::mutate(annotated = case_when(merged~T,
                                               T~F))%>%
      plot_ly(source = "plotly_ms2_sp")%>%
      highlight_key(~fragment_group
                    )%>%
      add_segments(x = ~x , xend = ~x,
                   y = I(0),yend = ~y,
                   alpha = 0.5,
                   colors = c("FALSE"= "#BBBBBB",
                              "TRUE"="#FF7F0E",
                              "CE=10"="#4DBBD5FF",
                              "CE=20" ="#FF7F0E",
                              "CE=40" ="#E64B35FF" ),
                   showlegend = F,
                   color = ~ merged)%>%
      dplyr::filter(annotated)%>%
      add_markers(x = ~x,
                  y = ~y,
                  colors = c("FALSE"= "#BBBBBB",
                             "TRUE"="#FF7F0E",
                             "CE=10"="#4DBBD5FF",
                             "CE=20" ="#FF7F0E",
                             "CE=40" ="#E64B35FF"),
                  color = ~ merged,
                  text = ~ hover_label,
                  #hovertemplate = "%{text}<extra></extra>",
                  hoverinfo = "text",
                  showlegend = F)%>%
      highlight(on = "plotly_click",
                off = NULL,
                selected = attrs_selected(
                  showlegend = F
                ),
                defaultValues=filter(sp.data,!is.na(fragment_group))$fragment_group[1])%>%
      layout(dragmode = "zoom",
             xaxis = list(range = c(0, max(sp.data$mz)*1.1),
                          title = "mz",
                          showgrid = T),
             yaxis = list(showgrid = T,
                          tickformat = ifelse(show.rawdata,
                                              ".1e",
                                              "0"),
                          range = c(-ymax,ymax),
                          title = "Relative intensity"))%>%
      event_register("plotly_click")

  }

  ### add message
  {
    fg <- sp.data%>%
      dplyr::filter(y<0,!is.na(fragment_group))%>%
      dplyr::pull(fragment_group)
    if (length(fg)==0) {
      p <- p%>%
        add_text(x = mean(range(sp.data$x)),
                 y = -ymax/2,text = "No fragment")
    }

  }

  return(p)
}


shiny_isotopologues_info <- function(iso_data,
                                     sample,
                                     iso_count){

  if (is.null(iso_data)|is.null(sample)|is.null(iso_count)) {
    return(NA)
  }
  ms1_purity <- ms2_count <- NULL
  purity_matrix <- iso_data$compound_info$purity_matrix
  if (iso_count%in% rownames(purity_matrix)&sample%in% colnames(purity_matrix))
    ms1_purity <- purity_matrix[iso_count,sample]

  ms2_count_matrix <- iso_data$compound_info$ms2_count
  if (iso_count%in% rownames(ms2_count_matrix)&sample%in% colnames(ms2_count_matrix))
    ms2_count <- ms2_count_matrix[iso_count,sample]


  paste0("Purity: ",format(ms1_purity,digits=2),"\n",
         "MS2 count: ",ms2_count)
}


shiny_get_sp_data <- function(iso_data,
                              sample ,
                              iso_count ){


  ### if
  {
    if (length(iso_data)==0|length(sample)==0|length(iso_count)==0) {
      return(NA)
    }

    if (!iso_count %in% names(iso_data$MSIP_result)) {
      return(NA)
    }


  }
  ### debug
  {
    #message(names(iso_data))
    #message(sample)
    #message(iso_count)
  }


  ### SP M0
  if (F)  {
    sp.m0 <-iso_data$Spectra$M0%>%
      unname()%>%
      do.call(c,.)

    sp.m0.frag.data <- CFM_annotate_isotopologues(sp.m0,
                                                  cfmd  = iso_data$CFM_annotation,
                                                  ppm = 10,
                                                  iso_count = 0)
    sp.m0.frag.data <- CFM_spectra_data_merge(sp.m0.frag.data,0)


    sp.m0.frag.data <- sp.m0.frag.data%>%
      dplyr::mutate(sp.id="M0",x = mz,
                    y =intensity)

  }else{
    msip.core.m0 <-iso_data$MSIP_result[["M0"]][[sample]]
    if (!is.null(msip.core.m0)) {
      sp.m0.frag.data <- msip.core.m0@Spectra_data%>%
        dplyr::mutate(sp.id="M0",x = mz,
                      y =intensity)
    }else{
      sp.m0.frag.data <- data.frame()
    }
  }

  ### sp data
  {
    msip.core <- iso_data$MSIP_result[[iso_count]][[sample]]
    if (is.null(msip.core)){
      sp.iso.frag.data<- sp.m0.frag.data[0,]
    }else{
      sp.iso.frag.data <-msip.core@Spectra_data%>%
        dplyr::filter(sp.id == "combined_sp"|is.na(fragment_group  ))%>%
        dplyr::arrange(sp.id)
      suppressWarnings(
        sp.iso.frag.data <- sp.iso.frag.data%>%
          dplyr::mutate(sp.id="Mn",
                        x = mz,
                        y = -intensity)
      )
    }




  }

  sp.data <- bind_rows(sp.m0.frag.data,sp.iso.frag.data)
  if (!nrow(sp.data)) {
    return(NULL)
  }
  sp.data <- sp.data %>%
    dplyr::mutate(annotated = !is.na(fragment_group),
                  collisionEnergy = as.character(collisionEnergy),
                  collisionEnergy = case_when(annotated~paste0("CE=",collisionEnergy),
                                              T~ " " ),
                  hover_label = paste0("Frag group: ",fragment_group,"\n",
                                       "mz: ",round(mz,4),"\n",
                                       "intensity: ",round(intensity,4)),
                  point.size = case_when(annotated~0.5,
                                         T~0))


  return(sp.data)


}


shiny_plotly_empty <- function(...){
  plot_ly(...)%>%
    add_markers(x = 0,y = 0,color = I("transparent"))%>%
    layout(xaxis = list(range = c(0,100)),
           yaxis = list(range = c(-100,100)))
}


shiny_vis_empty <- function(){

  nodes <- data.frame(
    id = 1,               # Unique ID for the node
    label = "No data",    # The text to display
    shape = "text",       # Node shape as text (no border)
    font.size = 20        # Adjust the font size as needed
  )
  visNetwork(nodes)

}

shiny_plotly_void <- function(...){
  plot_ly(...)%>%
    add_markers(x = 0,y = 0,color = I("transparent"))%>%
    layout(xaxis = list(range = c(0,100),
                        showline = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(range = c(-100,100),
                        showline = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>%
    layout(plot_bgcolor = 'rgba(0,0,0,0)',
           paper_bgcolor = 'rgba(0,0,0,0)',
           margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0),
           showlegend = FALSE)%>%
    config(displayModeBar = FALSE)
}


shiny_get_fg <- function(sp.data,
                         x){

  if (all(is.na(sp.data))) {
    return(NULL)
  }
  if (length(x)) {
    idx <- match_mz(x,sp.data$x,mz.ppm = 5)
    return(sp.data$fragment_group[idx])
  }
  return(na.omit(sp.data$fragment_group)[1])
}

shiny_get_fid <- function(iso_data,
                          fg){

  if (all(is.na(iso_data))) {
    return(NULL)
  }


  if (iso_data$compound_info$merged) {

    fg.split <- strsplit(fg,"_")[[1]]

    if (length(fg)) {
      fid <- iso_data[[paste0("CFM_annotation",fg.split[2])]]@fragment_define%>%
        dplyr::filter(fragment_group ==fg.split[1])%>%
        dplyr::pull(fragment_id )

      return(fid)
    }

  }else{
    if (length(fg)) {
      fid <-  iso_data$CFM_annotation@fragment_define%>%
        dplyr::filter(fragment_group ==fg)%>%
        dplyr::pull(fragment_id )

      return(fid)
    }
  }




  return(NA)
}

shiny_get_fg_ig <- function(iso_data,
                            fg = NULL,
                            fid = 1){

  if (all(is.na(iso_data))) {
    return(NULL)
  }
  if (all(is.na(fid))) {
    return(NULL)
  }
  if(is.null(fg)|is.na(fg)){
    return(NULL)
  }


  fg.split <- strsplit(fg,"_")[[1]]
  cfmd <- "CFM_annotation"
  if (iso_data$compound_info$merged) {
    cfmd <- paste0("CFM_annotation",fg.split[2])
  }
  if (fid == "seed") {
    fid <- 1
  }
  iso_data[[cfmd]]@fragment_igraph[[fid]]

}

shiny_get_frag_formula <- function(iso_data, fg = NULL,fid){


  if (all(is.na(iso_data))) {
    return(NULL)
  }
  if (all(is.na(fid))) {
    return(NULL)
  }
  if(is.null(fg)){
    return(NULL)
  }
  fg.split <- strsplit(fg,"_")[[1]]
  cfmd <- "CFM_annotation"
  if (iso_data$compound_info$merged) {
    cfmd <- paste0("CFM_annotation",fg.split[2])
  }
  if (fid == "seed") {
    fid <- 1
  }

  frag.def <- iso_data[[cfmd]]@fragment_define
  if (fid %in% frag.def$fragment_id) {
    formula <- frag.def%>%
      dplyr::filter(fragment_id == fid)%>%
      dplyr::pull(formula)
  }else{
    formula <- " "
  }
  return(formula)

}

shiny_get_atom_map <- function(iso_data,
                               fid,
                               fg = NULL,
                               prob = T){
  x <- list(NA,NA)
  #message("A",names(iso_data))
  #message("B",fid)
  if (is.null(iso_data)|is.null(fid)) {
    return(x)
  }
  if (all(is.na(iso_data))) {
    return(x)
  }

  fg.split <- strsplit(fg,"_")[[1]]
  cfmd <- "CFM_annotation"
  if (iso_data$compound_info$merged) {
    cfmd <- paste0("CFM_annotation",fg.split[2])
  }

  if (!fid%in%names(iso_data[[cfmd]]@fragment_atom_map) ) {
    return(x)
  }

  map.matrix <- iso_data[[cfmd]]@fragment_atom_map[[fid]]
  if (is.null(map.matrix)) return(x)
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

shiny_get_fg_map <- function(iso_data,
                             sample ,
                             iso_count){
  fg.map <- iso_data$MSIP_result[[iso_count]][[sample]]$fg.map
  return(fg.map)
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

    atom <- get_sdf_igraph_atom(ig)
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


shiny_vis_sdf_igraph <- function(ig,
                                 show_id = F,
                                 prob.border= NA,
                                 prob.fill =NA){
  if (is.null(ig)) {
    return(shiny_vis_empty())
  }

  message_with_time(show_id)
  ig%>%
    sdf_igraph_add_background_color(value = prob.fill)%>%
    sdf_igraph_add_border_color(prob.border)%>%
    vis_sdf_igraph(show_id = show_id)

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


as.HTML.checkbox.checked <- function(x){

  x.code <- x
  x.code[x] <- '<input type="checkbox" checked />'
  x.code[!x] <- '<input type="checkbox" />'
  return(x.code)

}

as.HTML.checkbox.logical <- function(x){

  x.logic <- x=='<input type="checkbox" checked />'
  x.logic
}
shiny_format_acq <- function(acq.list.table,acq.selected ){

  tb <- acq.list.table%>%
    dplyr::mutate(mz = round(mzmed,4),
                  rt = round(rtmed,0),
                  ms1_purity = round(ms1_purity,2),
                  selected_to_acq = acq.selected[feature_id] ,
                  selected_to_acq = as.HTML.checkbox.checked(selected_to_acq),
                  group_label = paste0(iso_seed,": ",str_short(name,50),", ",adduct),
                  isotope = paste0("M",iso_count)
                  )%>%
    dplyr::filter(!is.na(iso_seed),!is.na(compound_id),!is.infinite(iso_count))%>%
    dplyr::select(feature_id,mz,rt,isotope,ms1_purity,selected_to_acq,group_label)
  return(tb)

}


shiny_plotly_chrom_int_merged <- function(xchrom ,rtr = range(rtime(xchrom[1,1]))){

  chrom.pda <- pData(xchrom)
  groups <- unique(chrom.pda$group)
  groups.col <- ggsci::pal_npg()(length(groups))

  chrom.data <- get_chroms_data(xchrom)%>%
    dplyr::mutate(group = chrom.pda$group[col],
                  group = as.factor(group))%>%
    dplyr::filter(!is.na(intensity))
  chrom.def <- featureDefinitions(xchrom)
  names(groups.col) <- groups
  plot_ly(chrom.data)%>%
    layout(shapes = list(
      list( type = "line",
        x0 = chrom.def$rtmed, x1 = chrom.def$rtmed,
        y0 = 0,  y1 = max(chrom.data$intensity)*1.1,
        line = list(color = "#666666")
      ),
      list(type = "rect",
           line = list(color = "#88888833"),
           fillcolor= "88888833",
           x0 = chrom.def$peakRtMin, x1 = chrom.def$peakRtMax,
           y0 = 0,y1 = max(chrom.data$intensity))),
      plot_bgcolor = "#00000000")%>%
    add_lines(x = ~rt,
              y = ~intensity,
              split  = ~ col,
              color = ~ group,
              colors = groups.col,
              showlegend = F)%>%
    layout(xaxis = list(title =" retention time",showgrid = F),
           yaxis = list(showgrid = F,tickformat = ".1"))->p
    p

    chrom.val <- data.frame( chrom.pda,
                        val = as.vector(featureValues(xchrom,method="max",value = "maxo")))%>%
      dplyr::arrange(desc(group))%>%
      dplyr::mutate(group = as.factor(group),
                    sample.name = factor(sample.name,level = sample.name))
    chrom.val$val[is.na(chrom.val$val)] <- 0
    p%>%
      add_bars(data = chrom.val,
               y = ~sample.name,
               x = ~val ,
               colors = groups.col,
               color= ~ group,
               xaxis = "x2",
               yaxis = "y2")%>%
      layout(legend = list(x = 0.8,y=0.2,bgcolor = "#00000000"),
             xaxis2 = list(domain = c(0.8, 0.95),tickformat = ".1", anchor='y2'),
             yaxis2 = list(categoryorder = "2",domain = c(0.6, 0.95), anchor='x2', showticklabels = F))




}

shiny_plotly_chrom <- function(xchrom ,
                               col.map ,
                               rtr = range(rtime(xchrom[1,1]))){

  chrom.pda <- pData(xchrom)
  chrom.data <- get_chroms_data(xchrom)%>%
    dplyr::mutate(group = chrom.pda$sample.source[col],
                  group = as.factor(group))%>%
    dplyr::filter(!is.na(intensity))
  chrom.def <- featureDefinitions(xchrom)
  plot_ly(chrom.data)%>%
    layout(shapes = list(
      list( type = "line",
            x0 = chrom.def$rtmed, x1 = chrom.def$rtmed,
            y0 = 0,  y1 = max(chrom.data$intensity)*1.1,
            line = list(color = "#666666")
      ),
      list(type = "rect",
           line = list(color = "#88888833"),
           fillcolor= "#88888833",
           x0 = chrom.def$peakRtMin, x1 = chrom.def$peakRtMax,
           y0 = 0,y1 = max(chrom.data$intensity))),
      plot_bgcolor = "#00000000")%>%
    add_lines(x = ~rt,
              y = ~intensity,
              split  = ~ col,
              color = ~ group,
              colors = col.map,
              showlegend = F)%>%
    layout(xaxis = list(title =" retention time",
                        range = rtr,
                        showgrid = F),
           yaxis = list(showgrid = F,tickformat = ".1"))->p
  p




}

shiny_plotly_feature_int <- function(xchrom,
                                     col.map){

  chrom.val <- data.frame( pData(xchrom),
                           val = as.vector(featureValues(xchrom,method="max",value = "maxo")))%>%
    dplyr::arrange(desc(sample.source),desc(sample.name))%>%
    dplyr::mutate(group = as.factor(sample.source),
                  sample.name = factor(sample.name,level = sample.name))
  chrom.val$val[is.na(chrom.val$val)] <- 0
  plot_ly()%>%
    add_bars(data = chrom.val,
             y = ~sample.name,
             x = ~val ,
             colors = col.map,
             color= ~ group,
             showlegend = F)%>%
    layout(xaxis = list(tickformat = ".1",title = "Intensity",tickangle = 45),
           yaxis = list(title = " "))



}


shiny_plotly_feature_ratio <- function(xchrom,
                                     col.map,
                                     ratio_to_seed
                                     ){

  chrom.val <- data.frame( pData(xchrom))%>%
    dplyr::arrange(desc(sample.source),desc(sample.name))%>%
    dplyr::mutate(group = as.factor(sample.source),
                  val = ratio_to_seed[sampleNames],
                  sample.name = factor(sample.name,level = sample.name))
  chrom.val$val[is.na(chrom.val$val)] <- 0
  plot_ly()%>%
    add_bars(data = chrom.val,
             y = ~sample.name,
             x = ~val ,
             colors = col.map,
             color= ~ group,
             showlegend = F)%>%
    layout(xaxis = list(tickformat = ".1",title = "Ratio",tickangle = 45),
           yaxis = list(title = " "))



}

shiny_plotly_feature_purity <- function(xchrom,
                                       col.map,
                                       ms1_purity
){

  chrom.val <- data.frame( pData(xchrom))%>%
    dplyr::arrange(desc(sample.source),desc(sample.name))%>%
    dplyr::mutate(group = as.factor(sample.source),
                  val = ms1_purity[sampleNames],
                  sample.name = factor(sample.name,level = sample.name))
  chrom.val$val[is.na(chrom.val$val)] <- 0
  plot_ly()%>%
    add_bars(data = chrom.val,
             y = ~sample.name,
             x = ~val ,
             colors = col.map,
             color= ~ group,
             showlegend = F)%>%
    layout(xaxis = list(tickformat = ".1",range= c(0,1),title = "Purity",tickangle = 45),
           yaxis = list(title = " "))



}


shiny_get_C_prob <- function(msip.core.data){

 # ele <- get_sdf_igraph_atom(
 #   get_cfm_data_sdf_igraph(iso_data$CFM_annotation),"C")
 # x <- iso_data$MSIP_result[[iso_count]][[sample]]
 # if (all(is.na(x))|all(is.null(x))) {
 #   prob <- rep(NA,length(ele))
 #   names(prob) <- ele
 # }else{
 #   prob <- x$c.prob
 # }
  if (is.null(msip.core.data)){
    return(NULL)
  }
  prob <- msip.core.data@solve$Atom_prob
  return(prob)
}

shiny_DT_fg_include <- function(msip.core.data){


  if (is.null(msip.core.data)) {
    return(data.frame())
  }
  #if (all(is.null(msip.core.data$fg.map$frag.include))) {
  #  frag.include <- rep(T,nrow(msip.core.data$fg.map$frag.c.matrix))
  #  names(frag.include) <- rownames(msip.core.data$fg.map$frag.c.matrix)
  #}else{
  #  frag.include <- msip.core.data$fg.map$frag.include
  #}
  frag.include <- msip.core.data@FG_map@fragment.include
  data.frame(#fragment.group = names(frag.include),
             include = as.HTML.checkbox.checked(frag.include)
             )

}

shiny_get_fg_include <- function(msip.core.data){

  if (isEmpty(msip.core.data@FG_map)) {
    return(NULL)
  }
  #if (all(is.null(msip.core.data$fg.map$frag.include))) {
  #  frag.include <- rep(T,nrow(msip.core.data$fg.map$frag.c.matrix))
  #  names(frag.include) <- rownames(msip.core.data$fg.map$frag.c.matrix)
  #}else{
  #  frag.include <- msip.core.data$fg.map$frag.include
  #}
  frag.include <- msip.core.data@FG_map@fragment.include
  return(frag.include)
}

shiny_update_fg_include <- function(frag.include,x){


  frag.include[x] <- !frag.include[x]
  return(frag.include)

}

shiny_update_msip_core_data <- function(msip.core.data,frag.include){

  msip.core.data@FG_map@fragment.include <- frag.include
  return(msip.core.data)
}


shiny_update_iso_data_list <- function(iso.data.list,feature_id,
                                       sample,iso_count,msip.core.data){
  iso.data.list[[feature_id]][["MSIP_result"]][[iso_count]][[sample]] <- msip.core.data
  return(iso.data.list)
}

shiny_heatmap_fgmap <- function(msip.core.data){

  if (is.null( msip.core.data)|isEmpty(msip.core.data)) {
    plot.new()
    title("no data")
  }else{
    fg.map <- msip.core.data@FG_map
    heatmap_MSIPFragmentMap(fg.map )

  }

}


shiny_plotly_natural_ratio <- function(natural.ratio = 0.6 ){


  natural.ratio <- ifelse(natural.ratio>1,1,natural.ratio)

  df <- data.frame(
    name = c("labeled","natural"),
    value =c(1-natural.ratio,natural.ratio)
  )%>%
    dplyr::mutate(y = cumsum(value)-0.1,
                  label = paste0(name,"\n",
                                 format(value*100,digits = 2),"%"),
                  label = case_when(value < 0.1~"",
                                    T~label))

  plot_ly(df)%>%
    add_pie(label  = ~name, value  = ~value,
            textinfo = 'label+percent',
            marker = list(colors =c("labeled" = "#AA3310","natural" = "#eeeeee"))
          )%>%
    #open_visNet()
   # add_text(x = I(1), y = ~y,text = ~label)%>%
    layout(plot_bgcolor = 'rgba(0,0,0,0)',
           paper_bgcolor = 'rgba(0,0,0,0)',
           margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0),
           showlegend = FALSE)%>%
    config(displayModeBar = FALSE)

}


shiny_get_frag_iso_distribution<- function(msip.core,fg.id){


  if (is.null(msip.core)) return(NULL)
  msip.fgmap <- msip.core@FG_map

  if (is.null(msip.fgmap)|isEmpty(msip.fgmap)) return(NULL)

  fgm <- msip.fgmap@fragment.ratio.matrix
  if (!fg.id%in% rownames(fgm)) return(NULL)

  iso.dis <- fgm[fg.id,]
  return(iso.dis)
}



shiny_plotly_iso_distribution <- function(frag.iso.distribution){

  #frag.iso.distribution <- make_vector(runif(5),paste0("M",0:4))
  if ( is.null(frag.iso.distribution)) {
    return(shiny_plotly_void())
  }

  df <- data.frame(
    name = names(frag.iso.distribution),
    value = frag.iso.distribution
  )%>%
    dplyr::mutate(iso_count = str_isotope2_num(name),
                  iso_count = normalize_max_min(iso_count),
                  name = factor(name),
                  col = colramp()(iso_count),
                  col = substr(col,1,7),
                  col = case_when(iso_count==0~"#eeeeee",
                                  T~col))


  plot_ly(df)%>%
    dplyr::arrange(name)%>%
    add_pie(label = ~name,value = ~value,
            textinfo = 'label+percent',sort = F,
            marker = list(colors = ~I(col)))%>%
    layout(plot_bgcolor = 'rgba(0,0,0,0)',
           paper_bgcolor = 'rgba(0,0,0,0)',
           margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0),
           showlegend = FALSE)%>%
    config(displayModeBar = FALSE)

}

