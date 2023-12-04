


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
CFM_predict <- function(smiles_or_inchi_or_file = "CCCNNNC(O)O",
                          prob_thresh = 0.001,
                          param_adduct = "[M+H]+",
                          annotate_fragments = 1,
                          output_file_or_dir = tempfile(),
                          apply_postproc = 1,
                          suppress_exceptions = 1){

  model_param = CFM_get_param_config(param_adduct)
  CFM(result_dir = dirname(output_file_or_dir),
        cfmid_cmd = paste0(
          "cfm-predict ",
          smiles_or_inchi_or_file," ",
          prob_thresh," ",
          model_param,
          " 1 ",
          basename(output_file_or_dir), " ",
          apply_postproc," ",
          suppress_exceptions

        ))->cmd

  message("Make Command to docker: ")
  message(crayon::red(cmd)%>%crayon::reset() )
  #shell(cmd)
  message("Result export to: \n",crayon::red(output_file_or_dir))
  return(invisible(output_file_or_dir))
}







CFM_annotate<- function(smiles_or_inchi = "CCCNNNC(O)O",
                          spectrum_file = NULL,
                          id = "AN_ID",
                          ppm_mass_tol = 5.0,
                          abs_mass_tol = 0.01,
                          param_adduct = "[M+H]+",
                          output_file = tempfile()){
  ### check sp ID
  sp <- load_Spectra(spectrum_file)
  if (!"ID"%in%spectraVariables(sp)) {
    stop("Spectra variable ID not exist")
  }
  if (!id%in%sp$ID) {
    warning(id," not exist in ",spectrum_file, ", using id ",sp$ID[1])
    id <- sp$ID[1]
  }
  win.dir<- dirname(output_file)
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
          basename(output_file)
        ))->cmd
  cmd
  message("Make Command to docker: ")
  message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  message("Result export to: \n",crayon::red(output_file))
  return(invisible(output_file))
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
        this.peak.assign <- data.frame(enery = this.energy,
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

  cfm.data <- list(
    peak.assignment = peak.assignment,
    fragment_define1 = session.data.2,
    fragment_define2 = session.data.3,
    fragment_transition = session.data.4
  )%>%

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





