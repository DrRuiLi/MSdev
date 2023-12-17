CFM <- function(docker_run_param = "--rm",
                  result_dir = tempdir(),
                  cfmid_cmd = "cfm-id"
){
  paste0("docker run ",docker_run_param," ",
         " -v ",result_dir,":/cfmid",
         "  wishartlab/cfmid ",
         cfmid_cmd)->cmd
  return(cmd)

}


CFM_get_param_config <- function(adduct = c("[M+H]+","[M-H]-"),
                                   param = T,
                                   config = T){
  adduct = match.arg(adduct)

  param_file <- switch (adduct,
                        "[M+H]+" = " /trained_models_cfmid4.0/[M+H]+/param_output.log ",
                        "[M-H]-" =  " /trained_models_cfmid4.0/[M-H]-/param_output.log ",
                        " none "
  )
  config_file <- switch (adduct,
                        "[M+H]+" = " /trained_models_cfmid4.0/[M+H]+/param_config.txt ",
                        "[M-H]-" = " /trained_models_cfmid4.0/[M-H]-/param_config.txt ",
                        " none "
  )
  param_file <- ifelse(param,param_file," none ")
  config_file <- ifelse(config,config_file," none ")
  paste0(param_file,config_file)
}

#' CFM predict
#'
#' @param smiles_or_inchi_or_file
#' @param prob_thresh
#' @param param_adduct
#' @param annotate_fragments
#' @param output_file_or_dir
#' @param apply_postproc
#' @param suppress_exceptions
#'
#' @return
#' @export
#' @import crayon
#' @examples
CFM_predict <- function(smiles_or_inchi_or_file = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                          prob_thresh = 0.001,
                          param_adduct = "[M+H]+",
                          annotate_fragments = 1,
                          output_file_or_dir = NULL,
                          apply_postproc = 0,
                          suppress_exceptions = 1){

  out.file <- ifelse(is.null(output_file_or_dir),
                               tempfile(),
                               normalizePath(output_file_or_dir,mustWork = F))

  model_param = CFM_get_param_config(param_adduct)
  CFM(result_dir = dirname(out.file),
        cfmid_cmd = paste0(
          "cfm-predict ",
          smiles_or_inchi_or_file," ",
          prob_thresh," ",
          model_param,
          " 1 ",
          basename(out.file), " ",
          apply_postproc," ",
          suppress_exceptions

        ))->cmd

  message("Make Command to docker: ")
  message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file_or_dir)) {
    return(invisible(read_CFM_predict_result(out.file)))
  }else{
    message("Result export to: \n",crayon::red(output_file_or_dir))
    return(invisible(out.file))
  }
}

CFM_annotate<- function(smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                          spectrum_file = NULL,
                          id = "AN_ID",
                          ppm_mass_tol = 5.0,
                          abs_mass_tol = 0.01,
                          param_adduct = "[M+H]+",
                          output_file = NULL){

  out.file <- ifelse(is.null(output_file),
                     tempfile(),
                     normalizePath(output_file,mustWork = F))
  win.dir<- dirname(out.file)
  file.copy(spectrum_file,
            paste0(win.dir,"/",basename(spectrum_file)))
  model_param = CFM_get_param_config(param_adduct,param = T,config = T)

  CFM(result_dir = win.dir,
        cfmid_cmd = paste0(
          "cfm-annotate ",
          smiles_or_inchi," ",
          basename(spectrum_file)," ",
          id," ",
          ppm_mass_tol, " ",
          abs_mass_tol," ",
          model_param,
          basename(out.file)
        ))->cmd
  cmd
  message("Make Command to docker: ")
  message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file)) {
    return(invisible(read_CFM_annotate_result(out.file)))
  }else{
    message("Result export to: \n",crayon::red(output_file))
    return(invisible(out.file))
  }
}

CFM_fraggen <- function(smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                        max_depth = 2,
                        param_adduct = "[M+H]+",
                        output_file_or_dir = NULL){

  out.file <- ifelse(is.null(output_file_or_dir),
                     tempfile(),
                     normalizePath(output_file_or_dir,mustWork = F))
  CFM(result_dir = dirname(out.file),
      cfmid_cmd = paste0(
        "fraggraph-gen ",
        smiles_or_inchi," ",
        max_depth," ",
        switch(param_adduct,
               "[M+H]+"="+",
               "[M-H]-"="-"),
        " fullgraph ",
        basename(out.file), " "

      ))->cmd

  message("Make Command to docker: ")
  message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file_or_dir)) {
    return(invisible(read_CFM_fraggen_result(out.file)))
  }else{
    message("Result export to: \n",crayon::red(output_file_or_dir))
    return(invisible(out.file))
  }


}

read_CFM_annotate_result <- function(result_path = "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_annotate_result.txt"){

  cfm.data <- read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm.data),
                       line.data = cfm.data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep)
                  )%>%
    dplyr::mutate(session.energy.sep = grepl(pattern = "energy",x = line.data ),
                  session.energy = cumsum(session.energy.sep),
                  session.energy = case_when(session.type!=0~NA,
                                              T~session.energy)
                  )

  ### session 1
  {
    energy <- unique(cfm.df$session.energy)%>%na.omit()
    peak.assignment <- data.frame()
    for (i in seq_along(energy)  ) {
      this.energy <-cfm.df%>%
        dplyr::filter(session.energy==energy[i],session.energy.sep)%>%
        dplyr::pull(line.data)
      energy.data <- cfm.df%>%
        dplyr::filter(session.energy==energy[i],!session.energy.sep)%>%
        dplyr::mutate(
          assigned = grepl(pattern = "\\(",x = line.data),
          d1 =case_when(
            assigned~str_extract(line.data,".*(?= \\()"),
            T~line.data
          ),
          fragment_score = str_extract(line.data,"(?<= \\()[^\\)]*"))%>%
        dplyr::mutate(
          mz = str_extract(d1,"^([^\\s]*\\s){1}[^\\s]*"),
          intensity = str_extract(mz,"(?<=[^\\s]{1,10}\\s).*"),
          mz = str_extract(mz,".*(?=[^\\s]{1,10}\\s)"),
          fragment_id = str_extract(d1,"(?<=^[^\\s]{1,10}\\s[^\\s]{1,10}\\s)[^\\s].*"),
          fragment_id = strsplit(x = fragment_id,split = "\\s"),
          fragment_score = strsplit(x = fragment_score,split = "\\s")
        )

      if (nrow(energy.data) == 0) {
        next
      }
      for (j in 1:nrow(energy.data)) {
        this.peak.assign <- data.frame(energy = this.energy,
                        mz = energy.data$mz[j],
                        intensity = energy.data$intensity[j],
                        fragment_id = energy.data$fragment_id[[j]],
                        fragment_score = energy.data$fragment_score[[j]])

        peak.assignment <- bind_rows(peak.assignment,this.peak.assign)

      }


    }
    peak.assignment <- peak.assignment%>%
      dplyr::mutate(mz = as.numeric(mz),
                    intensity = as.numeric(intensity),
                    fragment_score = as.numeric(fragment_score))
    }

  ### session 2
  {
    session.data.2 <- cfm.df%>%
      dplyr::filter(session.type==1,!session.sep)%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smile =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
                    )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smile)
  }

  ### session 3
  {

    session.data.3 <- cfm.df%>%
      dplyr::filter(session.type==2,!session.sep)%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smile =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smile)

  }


  ### session 4 fragment transition
  {

    session.data.4 <- cfm.df%>%
      dplyr::filter(session.type==3,!session.sep)%>%
      dplyr::mutate(from  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    to = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    smile =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(from))%>%
      dplyr::select(from,to,smile)

  }

  ### re id
  {
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.3)))%>%
      `names<-`(session.data.3$fragment_id)
    peak.assignment <-peak.assignment%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id]))
    session.data.3 <- session.data.3%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)
    session.data.4 <- session.data.4%>%
      dplyr::mutate(from = unname(fragment.id.new[from]),
                    to = unname(fragment.id.new[to]))

  }

  cfm.data <- list(
    peak_assignment = peak.assignment,
    #fragment_define1 = session.data.2,
    fragment_define = session.data.3,
    fragment_transition = session.data.4
  )

    return(invisible(cfm.data))
}

read_CFM_predict_result <- function(result_path){

  cfm.data <- read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm.data),
                       line.data = cfm.data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep),
                  session.type = case_when(grepl(pattern = "^#",
                                                 x = line.data)~-1,
                                           T~session.type)
    )%>%
    dplyr::mutate(session.energy.sep = grepl(pattern = "energy",x = line.data ),
                  session.energy = cumsum(session.energy.sep),
                  session.energy = case_when(session.type!=0~NA,
                                             T~session.energy)
    )
  ### header
  {


    }

  ### session 1
  {
    energy <- unique(cfm.df$session.energy)%>%na.omit()
    peak.assignment <- data.frame()
    for (i in seq_along(energy)  ) {
      this.energy <-cfm.df%>%
        dplyr::filter(session.energy==energy[i],session.energy.sep)%>%
        dplyr::pull(line.data)
      energy.data <- cfm.df%>%
        dplyr::filter(session.energy==energy[i],!session.energy.sep)%>%
        dplyr::mutate(
          assigned = grepl(pattern = "\\(",x = line.data),
          d1 =case_when(
            assigned~str_extract(line.data,".*(?= \\()"),
            T~line.data
          ),
          fragment_score = str_extract(line.data,"(?<= \\()[^\\)]*"))%>%
        dplyr::mutate(
          mz = str_extract(d1,"^([^\\s]*\\s){1}[^\\s]*"),
          intensity = str_extract(mz,"(?<=[^\\s]{1,10}\\s).*"),
          mz = str_extract(mz,".*(?=[^\\s]{1,10}\\s)"),
          fragment_id = str_extract(d1,"(?<=^[^\\s]{1,10}\\s[^\\s]{1,10}\\s)[^\\s].*"),
          fragment_id = strsplit(x = fragment_id,split = "\\s"),
          fragment_score = strsplit(x = fragment_score,split = "\\s")
        )

      if (nrow(energy.data) == 0) {
        next
      }
      for (j in 1:nrow(energy.data)) {
        this.peak.assign <- data.frame(energy = this.energy,
                                       mz = energy.data$mz[j],
                                       intensity = energy.data$intensity[j],
                                       fragment_id = energy.data$fragment_id[[j]],
                                       fragment_score = energy.data$fragment_score[[j]])

        peak.assignment <- bind_rows(peak.assignment,this.peak.assign)

      }


    }
    peak.assignment <- peak.assignment%>%
      dplyr::mutate(mz = as.numeric(mz),
                    intensity = as.numeric(intensity),
                    fragment_score = as.numeric(fragment_score))
  }

  ### session 2
  {

    session.data.2 <- cfm.df%>%
      dplyr::filter(session.type==1,!session.sep)%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smile =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smile)

  }

  ### re id
  {
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.2)))%>%
      `names<-`(session.data.2$fragment_id)
    peak.assignment <-peak.assignment%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id]))
    session.data.2 <- session.data.2%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)

  }

  cfm.data <- list(
    peak_assignment = peak.assignment,
    #fragment_define1 = session.data.2,
    fragment_define = session.data.2
  )

  return(invisible(cfm.data))

}

read_CFM_fraggen_result <- function(result_path){

  cfm.data <- read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm.data),
                       line.data = cfm.data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep),
                  session.type = case_when(grepl(pattern = "^#",
                                                 x = line.data)~-1,
                                           T~session.type)
    )

  ### session 1 fragment def
  {

    session.data.1 <- cfm.df%>%
      dplyr::filter(session.type==0,!session.sep)%>%
      dplyr::mutate(intermediate = grepl(pattern =" Intermediate Fragment",x = line.data),
                    line.data= gsub(pattern =" Intermediate Fragment",x = line.data,replacement= ""))%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smile =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smile,intermediate)

    }


  ### session 2 fragment transition
  {

    session.data.2 <- cfm.df%>%
      dplyr::filter(session.type==1,!session.sep)%>%
      dplyr::mutate(from  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    to = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    neutral_loss = str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s]*(?=\\s)"),
                    smile =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(from))%>%
      dplyr::select(from,to,neutral_loss,smile)

  }

  ### re id
  {
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.1)))%>%
      `names<-`(session.data.1$fragment_id)
    session.data.1 <- session.data.1%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)
    session.data.2 <- session.data.2%>%
      dplyr::mutate(from = unname(fragment.id.new[from]),
                    to = unname(fragment.id.new[to]))

  }

  cfm.data <- list(
    fragment_define = session.data.1,
    fragment_transition = session.data.2
  )


  return(invisible(cfm.data))

}



plot_CFM_annotated_Spectra <- function(cfm_annoate_result){

  peak.assign <- cfm_annoate_result$peak.assignment
  fragment.define <- cfm_annoate_result$fragment_define2
  fragment.sdf <- smiles2sdf(fragment.define$smile)


  peak.assign.unique <- peak.assign%>%
    dplyr::group_by(mz)%>%
    dplyr::slice_max(fragment_score)%>%
    dplyr::mutate(fragment.define[match(fragment_id,fragment.define$fragment_id),])%>%
    dplyr::mutate(x = mz,
                  xend = mz,
                  y = 0,
                  yend = intensity,
                  annotated = !is.na(fragment_score))
  x.range <-c(min(peak.assign.unique$mz),
              max(peak.assign.unique$mz))
  y.range <-c(0,
              max(peak.assign.unique$intensity)*1.2)

  ggplot(peak.assign.unique)+
   # geom_hline(yintercept = 0 , size = 0.2,col = "grey")+
    geom_segment(aes(x = x,y =y,xend = xend,
                     yend = yend ,col = annotated),
                 size = 0.2,
                 show.legend = F)+
    geom_point(aes(x = x, y = yend ,
                   alpha = annotated,col = annotated),
               show.legend = F,size = 0.5)+
    scale_color_manual(values = c(`FALSE` = "grey",`TRUE` = "#80B1D3"))+
    scale_alpha_manual(values = c(`FALSE` = 0,`TRUE` =1))+
    scale_y_continuous(expand = expansion(0,0),
                       limits = y.range)+
    labs(x = "Mz",y = "Intensity")+
    theme_classic()+
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"inch"),
          axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))->p.sp

  p.annotated <- p.sp
  x.width = 10
  y.height = 10
  for (i in 1:nrow(peak.assign.unique)) {
    smile <- peak.assign.unique$smile[i]
    if (is.na(smile)) {
      next
    }
    sdf <- smiles2sdf(smile)[[1]]
    if (is.null(sdf)) {
      next
    }
    p.sdf <- ggplot_sdf(sdf,show_ele = T,cex = 0.3)
    p.annotated+annotation_custom(
      ggplotGrob(p.sdf),
      xmin = peak.assign.unique$x[i]-x.width/2,
      xmax = peak.assign.unique$x[i]+x.width/2,
      ymin = peak.assign.unique$yend[i]-y.height/2 ,
      ymax= peak.assign.unique$yend[i] + y.height
    )->p.annotated
  }

  return(p.annotated)
}





get_cfm_data_igraph <- function(cfm.data ){

  #cfm.data <- read_CFM_annotate_result()

  ### fragment
  fragment.data <- cfm.data$fragment_define2
  fragment.sdf <- suppressWarnings(smiles2sdf(fragment.data$smile))
  fragment.data$formula <- MF(fragment.sdf,addH = T)
  fragment.data$atom.count <- atomcountMA(fragment.sdf)%>%
    apply(1,sum)
  fragment.igraph <- list()
  for (i in 1:nrow(fragment.data)) {
    fragment.igraph[[fragment.data$fragment_id[i]]] <-
      get_sdf_igraph(fragment.sdf[[i]])
  }

  ### transition
  frag.trans.df <- cfm.data$fragment_transition%>%
    dplyr::mutate(from_formula = fragment.data[from,]$formula ,
                  to_formula = fragment.data[to,]$formula,
                  to_atom_count = fragment.data[to,]$atom.count)
  frag.trans.graph <- graph_from_data_frame(cfm.data$fragment_transition,
                                            vertices =fragment.data )
  frag.assigned <- cfm.data$peak.assignment$fragment_id%>%
    unique()%>%
    na.omit()


  ### graph calc
  frag.trans.df$intersection.atom <- 0
  for (i in 1:nrow(frag.trans.df)) {

    graph.to   <- fragment.igraph[[frag.trans.df$to[i]]]
    graph.from <- fragment.igraph[[frag.trans.df$from[i]]]
    graph.inter <- intersection(graph.to,graph.from,
                                byname =F,
                                keep.all.vertices = F)
    graph.inter <- graph.inter - V(graph.inter)[atom_1!=atom_2]
    frag.trans.df$intersection.atom[i] <- length(graph.inter)
  }
  frag.trans.df <- frag.trans.df%>%
    dplyr::mutate(retrieve.ratio = intersection.atom/to_atom_count)
  frag.trans.graph <- graph_from_data_frame(frag.trans.df,
                                            vertices =fragment.data )

  #visIgraph(frag.trans.graph)


  ### retrieve atom
  root.frag.graph <-fragment.igraph[[1]]
  fragment.data$dis.to.root <- 0
  fragment.data$retrieve.ratio <- 1
  for (i in 2:nrow(fragment.data)) {


      this.frag.id <- fragment.data$fragment_id[i]
      this.path <- all_simple_paths(frag.trans.graph,
                                    1,this.frag.id)
      if (!length(this.path)==0) {
        retrieve.ratio <- sapply(this.path, function(x){
          get_edges_from_path(frag.trans.graph,x)$retrieve.ratio%>%
            prod()
        })
        this.path <- this.path[[which.max(retrieve.ratio)]]
      }else{
        retrieve.ratio <- 0
      }

      fragment.data$retrieve.ratio[i] <- mean(retrieve.ratio)
      fragment.data$dis.to.root[i] <- length(this.path) -1

    }

  fragment.data.tmp <- fragment.data%>%
    dplyr::arrange(dis.to.root)
  V(fragment.igraph[[1]])$root_atom_id <- V(fragment.igraph[[1]])$id
  for (i in 2:nrow(fragment.data.tmp)) {

    this.frag.id <- fragment.data.tmp$fragment_id[[i]]
    message("Tracing atom for ",this.frag.id)
    fragment.igraph[[this.frag.id]] <- get_atom_id_from_parent(fragment.igraph[[1]],
                                                               fragment.igraph[[this.frag.id]])
  }


  fragment.data$root_atom_id  <- lapply(fragment.igraph,function(x){  V(x)$root_atom_id})


  cfm.data$fragment_define2 <- fragment.data
  cfm.data$fragment_transition <- frag.trans.df
  cfm.data$fragment_igraph <- fragment.igraph

  return(cfm.data)
}










