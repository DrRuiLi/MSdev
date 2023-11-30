


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
  model_param = CFM_get_param_config(param_adduct,param = F,config = T)

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
  #shell(cmd)
  message("Result export to: \n",crayon::red(output_file))
  return(invisible(output_file))
}


