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







setClass("CFM_data",
         slots = list(
           peak_assignment = "data.frame",
           fragment_define = "data.frame",
           fragment_transition = "data.frame",
           fragment_igraph = "list",
           fragment_sdf = "SDFset",
           fragment_atom_map = "list",
           fragment_group = "data.frame",
           fragment_group_map = "matrix"
         ))


setMethod("show",signature ="CFM_data",definition =
            function(object){
              message("CFM_data with ",nrow(object@fragment_define)," fragment")
            } )


#' CFM predict
#' @description
#' see CFM
#' @describeIn CFM predict
#' @param smiles_or_inchi_or_file smiles
#' @param prob_thresh prob
#' @param param_adduct adduct
#' @param annotate_fragments logic
#' @param output_file_or_dir path
#' @param apply_postproc logic
#' @param suppress_exceptions logic
#'
#' @return list of cfm data
#' @export

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

  #message("Make Command to docker: ")
  #message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file_or_dir)) {
    cfm.result <- read_CFM_predict_result(out.file)
    return(invisible(cfm.result))
  }else{
    #message("Result export to: \n",crayon::red(output_file_or_dir))
    return(invisible(out.file))
  }
}

#' CFM_annotate
#' @describeIn CFM anntate
#'
#' @param smiles_or_inchi smiles
#' @param spectrum_file spectra
#' @param id id
#' @param ppm_mass_tol ppm
#' @param abs_mass_tol mz tol
#' @param param_adduct adduct
#' @param output_file file path
#'
#' @return cfm data
#'
#' @export
#'

CFM_annotate<- function(smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                          spectrum_file = NULL,
                          id = "AN_ID",
                          ppm_mass_tol = 5.0,
                          abs_mass_tol = 0,
                          param_adduct = "[M+H]+",
                          output_file = NULL,...){

  out.file <- ifelse(is.null(output_file),
                     tempfile(),
                     normalizePath(output_file,mustWork = F))
  win.dir<- dirname(out.file)
  if ("Spectra" %in% class(spectrum_file) ) {
    spectrum_file <- export_Spectra_peak_list_for_cfm(spectrum_file,tempfile())
  }
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
  #message("Make Command to docker: ")
  #message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file)) {
    return(invisible(read_CFM_annotate_result(out.file)))
  }else{
    message("Result export to: \n",crayon::red(output_file))
    return(invisible(out.file))
  }
}


#' CFM_annotate_by_fraggen
#' @describeIn CFM fraggen and annotate
#' @description
#' this function generate CFM_data by CFM_fraggen,
#' not rely on Spectra input
#'
#' @note
#' filter fragment of mz > 30 to avoid error when smiles transform to sdf
#'
#'
#' @param smiles_or_inchi smiles
#' @param spectrum_file NULL
#' @param max_depth 1
#' @param id ID
#' @param ppm_mass_tol 20
#' @param abs_mass_tol 0
#' @param param_adduct adduct
#' @param output_file NULL
#'
#'
#' @return CFM_data
#' @export
#'
CFM_annotate_by_fraggen <- function(
    smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
    spectrum_file = NULL,
    max_depth = 1,
    id = "AN_ID",
    ppm_mass_tol = 5.0,
    abs_mass_tol = 0,
    param_adduct = "[M+H]+",
    output_file = NULL,...){

  .Deprecated("CFM_annotate_by_predict")
  cfm.fragen <- CFM_fraggen(smiles_or_inchi,
                            max_depth = max_depth,
                            param_adduct = param_adduct)


  peak.data <- cfm.fragen@fragment_define%>%
    dplyr::filter(!intermediate)%>%
    dplyr::mutate(energy = "energy0",
                  intensity = NA,
                  fragment_score= NA)%>%
    dplyr::add_row(energy = c(rep("energy1",nrow(.)),
                              rep("energy2",nrow(.))),
                   fragment_mz=rep(.$fragment_mz,2),
                   fragment_id = rep(.$fragment_id,2))%>%
    dplyr::select(energy,
                  mz =fragment_mz,
                  intensity,
                  fragment_id,
                  fragment_score)%>%
    remove_rownames()

  cfm.fragen@peak_assignment <- peak.data
  return(cfm.fragen)


}


#' CFM_annotate_by_predict
#'
#' @describeIn CFM predict and annotate
#' @description
#' [CFM_predict] and then [CFM_annotate]
#'
#' @export
#' @inheritParams  CFM_annotate
CFM_annotate_by_predict <- function(
    smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
    id = "AN_ID",
    ppm_mass_tol = 5.0,
    abs_mass_tol = 0,
    param_adduct = "[M+H]+",
    output_file = NULL,...){

  cfm.pred <- CFM_predict(smiles_or_inchi,
                            param_adduct = param_adduct)
  cfm.pred.sp <- get_CFM_data_Spectra(cfm.pred)
  cfm.anno <- CFM_annotate(smiles_or_inchi = smiles_or_inchi,id = id,
               spectrum_file = cfm.pred.sp,
               ppm_mass_tol=ppm_mass_tol,
               abs_mass_tol=abs_mass_tol,
               param_adduct=param_adduct   )
  return(cfm.anno)


}


#' @describeIn CFM fraggen
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

  #message("Make Command to docker: ")
  #message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file_or_dir)) {
    return(invisible(read_CFM_fraggen_result(out.file)))
  }else{
    #message("Result export to: \n",crayon::red(output_file_or_dir))
    return(invisible(out.file))
  }


}

read_CFM_annotate_result <- function(result_path = "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_annotate_result.txt")
  {

  cfm_data <- readr::read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm_data),
                       line.data = cfm_data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep)
                  )%>%
    dplyr::mutate(session.energy.sep = grepl(pattern = "energy",x = line.data ),
                  session.energy = cumsum(session.energy.sep),
                  session.energy = dplyr::case_when(session.type!=0~NA,
                                              T~session.energy)
                  )

  ### session 1
  {
    energy <- unique(cfm.df$session.energy)%>%na.omit()
    peak_assignment <- data.frame()
    for (i in seq_along(energy)  ) {
      this.energy <-cfm.df%>%
        dplyr::filter(session.energy==energy[i],session.energy.sep)%>%
        dplyr::pull(line.data)
      energy.data <- cfm.df%>%
        dplyr::filter(session.energy==energy[i],!session.energy.sep)%>%
        dplyr::mutate(
          assigned = grepl(pattern = "\\(",x = line.data),
          d1 =dplyr::case_when(
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

        peak_assignment <- bind_rows(peak_assignment,this.peak.assign)

      }


    }
    peak_assignment <- peak_assignment%>%
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
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
                    )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles)
  }

  ### session 3
  {

    session.data.3 <- cfm.df%>%
      dplyr::filter(session.type==2,!session.sep)%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles)

  }


  ### session 4 fragment transition
  {

    session.data.4 <- cfm.df%>%
      dplyr::filter(session.type==3,!session.sep)%>%
      dplyr::mutate(from  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    to = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(from))%>%
      dplyr::select(from,to,smiles)%>%
      dplyr::distinct(from,to,smiles)

  }

  ### re id
  {
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.3)))%>%
      `names<-`(session.data.3$fragment_id)
    peak_assignment <-peak_assignment%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id]))
    session.data.3 <- session.data.3%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)
    session.data.4 <- session.data.4%>%
      dplyr::mutate(from = unname(fragment.id.new[from]),
                    to = unname(fragment.id.new[to]))

  }

  cfm_data <- list(
    peak_assignment = peak_assignment,
    #fragment_define1 = session.data.2,
    fragment_define = session.data.3,
    fragment_transition = session.data.4
  )

  ### create cfm_data
  {
    cfm_data <- new("CFM_data")
    cfm_data@peak_assignment <- peak_assignment
    cfm_data@fragment_define <- session.data.3
    cfm_data@fragment_transition <- session.data.4

  }

  return(cfm_data)
}

#' read_CFM_predict_result
#' @describeIn CFM read_CFM_predict_result
#' @param result_path path
#'
#' @return list of cfm data
#' @export
#'

read_CFM_predict_result <- function(result_path =  "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_predict_result.txt"){

  cfm_data <- readr::read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm_data),
                       line.data = cfm_data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep))
  cfm.df <- cfm.df[!grepl(pattern = "^#",  x = cfm.df$line.data),]
    #dplyr::mutate(
    #  session.type = dplyr::case_when(grepl(pattern = "^#",
    #                                             x = line.data)~-1,
    #                                       T~session.type)
    #)%>%
  cfm.df <- cfm.df%>%
    dplyr::mutate(session.energy.sep = grepl(pattern = "energy",x = line.data ),
                  session.energy = cumsum(session.energy.sep),
                  session.energy = dplyr::case_when(session.type!=0~NA,
                                             T~session.energy)
    )
  ### header
  {


    }

  ### session 1
  {
    energy <- unique(cfm.df$session.energy)%>%na.omit()
    peak_assignment <- data.frame()
    for (i in seq_along(energy)  ) {
      this.energy <-cfm.df%>%
        dplyr::filter(session.energy==energy[i],session.energy.sep)%>%
        dplyr::pull(line.data)
      energy.data <- cfm.df%>%
        dplyr::filter(session.energy==energy[i],!session.energy.sep)%>%
        dplyr::mutate(
          assigned = grepl(pattern = "\\(",x = line.data))%>%
        dplyr::mutate(
          d1 = case_when(
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

        peak_assignment <- bind_rows(peak_assignment,this.peak.assign)

      }


    }
    peak_assignment <- peak_assignment%>%
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
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles)

  }

  ### re id
  {
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.2)))%>%
      `names<-`(session.data.2$fragment_id)
    peak_assignment <-peak_assignment%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id]))
    session.data.2 <- session.data.2%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)

  }

  cfm_data <- list(
    peak_assignment = peak_assignment,
    #fragment_define1 = session.data.2,
    fragment_define = session.data.2
  )


  ### create cfm_data
  {
    cfm_data <- new("CFM_data")
    cfm_data@peak_assignment <- peak_assignment
    cfm_data@fragment_define <- session.data.2


  }
  return(invisible(cfm_data))

}

read_CFM_fraggen_result <- function(result_path){

  cfm_data <- readr::read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm_data),
                       line.data = cfm_data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep),
                  session.type = dplyr::case_when(grepl(pattern = "^#",
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
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles,intermediate)

    }


  ### session 2 fragment transition
  {

    session.data.2 <- cfm.df%>%
      dplyr::filter(session.type==1,!session.sep)%>%
      dplyr::mutate(from  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    to = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    neutral_loss = str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s]*(?=\\s)"),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(from))%>%
      dplyr::select(from,to,neutral_loss,smiles)

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


  ### create cfm_data
  {
    cfm_data <- new("CFM_data")
    cfm_data@fragment_define <- session.data.1
    cfm_data@fragment_transition <- session.data.2

  }


  return(invisible(cfm_data))

}



plot_CFM_annotated_Spectra <- function(cfmd){

  peak.assign <- cfmd@peak_assignment
  fragment.define <- cfmd@fragment_define

  get_CFM_data_Spectra(cfmd)%>%
    combineSpectra()%>%
    plot_Spectra(label.top = 0)

}



#' CFM_annotate_isotopologues
#' @describeIn MSIP CFM_annotate_isotopologues
#'
#' @param sp Spectra
#' @param cfmd cfmd_data
#' @param isotope `[13]C`
#' @param iso_count num
#' @param ppm 20
#'
#' @return null
#' @export
#'
CFM_annotate_isotopologues <- function(sp,
                                 cfmd,
                                 iso_ele = "[13]C",
                                 iso_count = 0 ,
                                 ppm = 20){

  if (!"fragment_group"%in% colnames(cfmd@peak_assignment)) {
    cfmd <- cfm_data_get_fragment_group(cfmd)
  }

  cfm.peaks.data <- cfmd@peak_assignment%>%
    dplyr::filter(!is.na(fragment_id))%>%
    dplyr::mutate(collisionEnergy = case_when(energy == "energy0"~10,
                                              energy == "energy1"~20,
                                              energy == "energy2"~40,
    ))
  if (!nrow(cfm.peaks.data))
    cfm.peaks.data <-  cfmd@peak_assignment%>%
    dplyr::mutate(mz = 0)

  diff.formula <- paste0(iso_ele,get_ele_uniso(iso_ele),"-1")
  iso.mz.diff <- (0:iso_count)*MSCC::chemform_mz(diff.formula)
  mz.labeled.m <- matrixSub(cfm.peaks.data$mz,-iso.mz.diff)%>%
    `colnames<-`(paste0("M",0:iso_count))%>%
    as.data.frame()%>%
    dplyr::mutate(fragment_group=cfm.peaks.data$fragment_group)%>%
    dplyr::distinct(M0,.keep_all = T)%>%
    tidyr::pivot_longer(paste0("M",0:iso_count),
                        names_to = "iso_count",
                        values_to = "mz")



  sp.data <- get_Spectra_data(sp,var = "collisionEnergy")%>%
    dplyr::mutate(fragment_group = NA,iso_count=NA)
  if (nrow(sp.data)) {


  sp.data <- sp.data%>%
    dplyr::mutate(idx = match_mz(mz1=mz,
                                mz2 = mz.labeled.m$mz,
                                mz.ppm = ppm),
                  fragment_group = mz.labeled.m$fragment_group[idx],
                  iso_count = mz.labeled.m$iso_count[idx])%>%
    dplyr::select(-idx)%>%
    dplyr::group_by(fragment_group,sp.id)%>%
    dplyr::mutate(
      ratio = case_when(
      !is.na(fragment_group)~intensity/sum(intensity),
      T~NA),
      int_sum =  case_when(
        !is.na(fragment_group)~sum(intensity),
        T~NA))%>%
    dplyr::ungroup()
  }


  return(sp.data)
}


#' CFM_data_get_igraph
#'
#' @describeIn CFM_data CFM_data_get_igraph
#' @param object cfmd
#'
#' @return cfmd
#' @export
CFM_data_get_igraph <- function(object){

  ### fragment def
  {
    fragment.data <- object@fragment_define
    fragment.sdf <- suppressWarnings(get_smiles_sdf(fragment.data$smiles))
    fragment.data$formula <- get_sdf_formula(fragment.sdf)

  }

  ### filter valid smiles and sdf
  {
    fragment.data <- fragment.data%>%
      dplyr::filter(!is.na(formula))
    object@fragment_define <- fragment.data
    object@fragment_transition <-
      object@fragment_transition %>%
      dplyr::filter(from %in% object@fragment_define $fragment_id,
                    to %in% object@fragment_define $fragment_id
      )

  }

  ### igraph
  {
    fragment.sdf <- suppressWarnings(get_smiles_sdf(fragment.data$smiles))
    cid(fragment.sdf) <- fragment.data$fragment_id
    fragment.igraph <- get_sdf_igraph(fragment.sdf)
    names(fragment.igraph) <- fragment.data$fragment_id
  }




  ### save to object
  {

    fragment.data -> object@fragment_define
    object@fragment_igraph <- fragment.igraph
    object@fragment_sdf <- fragment.sdf
  }

  return(object)
}



#' Atom tracing map
#' @describeIn Atom_tracing_map get atom map
#'
#'
#' @param object cfmd
#' @param iso_ele `[13]C`
#' @param BPPARAM [BiocParallel::SerialParam]
#'
#' @note CFM_data_get_atom_map take too much times when construct, major in get_CFM_data_trans_map>>>fmcs
#'
#'
#'
#' @return cfmd
#' @export
#'
CFM_data_get_atom_map <- function(object,
                                  iso_ele = "[13]C",
                                  BPPARAM = SerialParam()){




  ### trans.map
  {
    #object <- cfmd
    message_with_time("trans map")
    fragment.trans <- object@fragment_transition
    trans.maps <- bplapply(1:nrow(fragment.trans),
                           get_CFM_data_trans_map,cfmd = object,
                           iso_ele=iso_ele,
                           BPPARAM = BPPARAM)
    names(trans.maps) <- paste0(fragment.trans$from,fragment.trans$to)
    trans.maps.stat <- check_CFM_data_trans_map(object,trans.maps = trans.maps,iso_ele = iso_ele)
    fragment.trans[,colnames(trans.maps.stat)] <- trans.maps.stat
    object@fragment_transition <- fragment.trans

  }

  ### filter involid trans
  {
    object_bak <- object
    object <- CFM_data_remove_trans(object)

  }


  ### fragment map
  {
    message_with_time("fragment map")
    fragment.atom.map <-list()
    ### weight for path selection
    fragment.trans <- object@fragment_transition%>%
      dplyr::mutate(id =  paste0(from,to),
                    weight = loss.distance  )
    fragment.data <- object@fragment_define
    fragment.data$ratio <- NA
    fragment.data$bond.score <- NA
    fragment.data$cumsum.loss.distance <- NA
    fragment.igraph <- object@fragment_igraph
    if (nrow(fragment.trans)) {
      ig.trans <- igraph::graph_from_data_frame(fragment.trans)
      for (i in 1:nrow(fragment.data)) {

       # message_with_time(i)

        this.frag <- fragment.data$fragment_id[i]
        fragment.atom.map[[i]] <- NA
        if(i==1){
          ele <- get_sdf_igraph_atom(fragment.igraph[[1]])
          maps <- diag(nrow = length(ele))
          rownames(maps)<-colnames(maps)<-ele
          fragment.atom.map[[i]] <- maps
          fragment.data$ratio[i] <-  1
          next
        }

        ### find path
        {
          #this.path <- all_simple_paths(ig.trans,from = 1,to = this.frag,
          #                              mode = "out",
          #                              cutoff = distances(ig.trans,1,this.frag,mode = "out")+5)
          #this.paths <- shortest_paths(ig.trans,from = 1,
          #                            to = this.frag,
          #                            mode = "out",output  = "vpath")$vpath
          this.paths <- igraph::all_shortest_paths(ig.trans,
                                          from = 1,
                                      to = this.frag,
                                      mode = "out")$vpath
          this.epaths <- lapply(this.paths,function(path){
            epath <- paste0(names(path),names(path)[-1])
            epath[1:length(epath)-1]
          })
          this.path.ratio <- sapply(this.epaths,function(this.epath){
            sapply(this.epath,function(epath){
              idx.path <- match(epath,fragment.trans$id)
              r <- fragment.trans$ratio[idx.path]
              return(r)
            })%>%prod()
          })
          this.path.bond.score <- sapply(this.epaths,function(this.epath){
            sapply(this.epath,function(epath){
            idx.path <- match(epath,fragment.trans$id)
            r <- fragment.trans$bond.score[idx.path]
            return(r)
            })%>%prod()
          })
          this.path.cumsum.loss.distance <- sapply(this.epaths,function(this.epath){
            sapply(this.epath,function(epath){
              idx.path <- match(epath,fragment.trans$id)
              r <- fragment.trans$loss.distance[idx.path]
              return(r)
            })%>%sum()
          })
        }

        ### prod maps
        {
          maps.list <- lapply(this.epaths,function(this.epath){
            if (!length(this.epath)==0) {
              maps <- trans.maps[this.epath]
              #trans.atom.map[trans.idx] <- maps
              while(length(maps)>1){
                maps[[2]] <- maps[[1]]%*% maps[[2]]
                maps[[1]] <-NULL
              }
              return(maps[[1]])
            }
          })
          map <- do.call(sum_matrix,maps.list)/length(maps.list)

        }

        ### return
        {

          fragment.atom.map[[i]] <- map
          fragment.data$ratio[i] <-  mean(this.path.ratio)
          fragment.data$bond.score[i] <-  mean(this.path.bond.score)
          fragment.data$cumsum.loss.distance[i] <-  mean(this.path.cumsum.loss.distance)
        }
      }
      names(fragment.atom.map) <-fragment.data$fragment_id
    }


  }



  object@fragment_atom_map <- fragment.atom.map
  object@fragment_define <- fragment.data

  return(object)

}


CFM_data_remove_trans <- function(object){


  .f <- function(object){
    fragment.trans <-object@fragment_transition
    fragment.trans <-fragment.trans[fragment.trans$volid,]

    trans.ig <- igraph::graph_from_data_frame(fragment.trans)
    dis.to.fragment1 <- distances(trans.ig,mode  = "out",
                                  v = object@fragment_define$fragment_id[1])
    reachable <- colnames(dis.to.fragment1)[!is.infinite(dis.to.fragment1)]
    to.remove <- setdiff(object@fragment_define$fragment_id,reachable)
    ### remove
    {
      #object <- cfmd
      object@peak_assignment <-object@peak_assignment %>%
        dplyr::filter(!fragment_id%in% to.remove)

      object@fragment_define <-object@fragment_define %>%
        dplyr::filter(!fragment_id%in% to.remove)

      object@fragment_transition <-fragment.trans %>%
        dplyr::filter(!(from%in% to.remove|to %in% to.remove))

      object@fragment_igraph <-object@fragment_igraph[
        !names(object@fragment_igraph)%in%to.remove]

      object@fragment_sdf <-object@fragment_sdf[
        !cid(object@fragment_sdf)%in%to.remove]

      object@fragment_atom_map <-object@fragment_atom_map[
        !names(object@fragment_atom_map)%in%to.remove]

    }

    return(object)
  }

  ### iteration
  i <- 1
  object<- .f(object)
  while(any(is.infinite(distances(get_CFM_data_trans_igraph(object),1,mode = "out")))&i<=5){
    object <- .f(object)
  }
  if (i==5)
    warning("CFM_data_remove_trans abnormal")

  return(object)

}


get_CFM_data_trans_map <- function(cfmd,trans_id,iso_ele="[13]C"){

  #message("trans: ",trans_id)
  fragment.trans <- cfmd@fragment_transition
  fragment.igraph <- cfmd@fragment_igraph
  fragment.sdf <- cfmd@fragment_sdf
  ig.parent <- fragment.igraph[[fragment.trans$from[trans_id]]]
  ig.product <- fragment.igraph[[fragment.trans$to[trans_id]]]
  sdf.parent <- fragment.sdf[[fragment.trans$from[trans_id]]]
  sdf.product <- fragment.sdf[[fragment.trans$to[trans_id]]]

  maps <- get_atom_map(sdf.parent ,sdf.product ,
                       ig.parent ,ig.product ,iso_ele=iso_ele)
  return(maps)
}


check_CFM_data_trans_map <- function(cfmd,
                                     iso_ele = "[13]C",
                                     trans.maps =NULL){

  if (is.null(trans.maps)) {
    trans.maps <- sapply(1:nrow(cfmd@fragment_transition),
                         get_CFM_data_trans_map,cfmd = cfmd)
  }
  bond.score <- sapply(trans.maps,function(x) attributes(x)$bond.score)
  #x <- trans.maps[[30]]
  atom.ele <- get_ele_uniso(iso_ele )
  trans.maps <- lapply(trans.maps,function(x){
    x[grepl(atom.ele,rownames(x)),grepl(atom.ele,colnames(x)),drop = F]
  })
  n_parent <- sapply(trans.maps,function(x){
    nrow(x)
  })
  n_atoms <-sapply(trans.maps,function(x){
    ncol(x)
  })
  n_atoms_compose_map <- sapply(trans.maps,function(x){

    sum(rowSums(x)==1)
  })
  n_atoms_certain_map <- sapply(trans.maps,function(x){
    max.prob <- apply(x,2,max)
    sum( max.prob == 1 )
  })
  n_atoms_noncertain_map <- sapply(trans.maps,function(x){
    max.prob <- apply(x,2,max)
    sum( max.prob <1&max.prob>0 )
  })
  n_atoms_non_map <- sapply(trans.maps,function(x){
    max.prob <- apply(x,2,max)
    sum( max.prob==0 )
  })

  map.stat <- data.frame(n_atoms,
                         n_parent,
                  n_atoms_compose_map,
                  n_atoms_certain_map,
                  n_atoms_noncertain_map,
                  n_atoms_non_map,
                  bond.score)%>%
    dplyr::mutate(volid = n_atoms_non_map ==0,
                  ratio = n_atoms_compose_map/n_atoms,
                  atom.loss = n_atoms - n_atoms_compose_map,
                  bond.loss =1-bond.score,
                  loss.distance = atom.loss+bond.loss
                  )

  return(map.stat)
}

get_CFM_data_trans_igraph <- function(object){



  node.df <- object@fragment_define %>%
    dplyr::mutate(id = fragment_id,
                  label = id,
                  no = 1:n(),
                  color.border =case_when(
                    no == 1~ "rgba(100, 100, 100, 0.8)",
                      T ~ "rgba(100, 100, 100, 0.8)"
                    ),
                  color.background = case_when(
                    T ~ "rgba(255, 255, 255, 0.8)"
                  ),
                  #color.highlight.background = "rgba(43, 124, 233, 0.1)",
                  #color.highlight.border = "rgba(43, 124, 233, 1)",
                  font.size = 30,
                  borderWidth = 5,
                  size = case_when(no == 1~ 100,
                                   T~50))
  edge.df <- object@fragment_transition%>%
    dplyr::mutate(match.str = paste0(from,to),
                  arrows.to.scaleFactor = 2,
                  atom.loss = n_atoms - n_atoms_compose_map,
                  bond.loss =1-bond.score,
                  color =  "rgba(100, 100, 100, 0.2)",
                  width = 10,
                  smooth = F,
                  loss.distance = atom.loss+bond.loss,
                  length = normalize_max_min(loss.distance)*1500,
                  length = 700,
                  no = 1:n())

  frag.trans.graph <- igraph::graph_from_data_frame(edge.df,
                                            vertices = node.df )
  return(frag.trans.graph)
}



#' cfm_data_get_FG_map
#'
#' @describeIn CFM_data cfm_data_get_FG_map
#'
#' @param cfm_data cfmd
#' @param ppm 10
#' @return cfmd
cfm_data_get_FG_map <- function(cfm_data,iso_ele = "[13]C",ppm = 5){


  ### Fragment group
  {

    fg <- groupMz(cfm_data@fragment_define$fragment_mz,ppm)
    cfm_data@fragment_define$fragment_group <- paste0("FG",num2str(fg))
    cfm_data@peak_assignment$fragment_group <-
      cfm_data@fragment_define$fragment_group[match(
        cfm_data@peak_assignment$fragment_id,
        cfm_data@fragment_define$fragment_id)]
    fg.count <- table(cfm_data@fragment_define$fragment_group)
    fg.df <- data.frame(
      fragment_group = paste0("FG",num2str(sort(unique(fg))))
    )%>%
      dplyr::mutate(
        fragment_count = as.numeric(fg.count[fragment_group]),
        fragment_mz = cfm_data@fragment_define$fragment_mz[match(fragment_group,
                                                   cfm_data@fragment_define$fragment_group)]
      )

  }


  ###  FG map
  {

    target_atoms <- get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfm_data),get_ele_uniso(iso_ele))
    frag.atom.matrix <- matrix(ncol = length(target_atoms),
                               nrow = nrow(fg.df),
                               dimnames = list(fg.df$fragment_group,
                                               target_atoms))
    for (i.fg in seq_len(nrow(fg.df))) {

      this.frag.group <- fg.df$fragment_group[i.fg]
      this.frags <- cfm_data@fragment_define[cfm_data@fragment_define$fragment_group==this.frag.group,]
      this.frag.atom <- get_cfm_data_fragment_group_atom_map(cfm_data,this.frag.group)
      this.frag.c <- this.frag.atom[target_atoms]
      #this.iso.expectation <- sum(str_extract_num(names(this.frag.ratio))*this.frag.ratio)
      #this.frag.c <- this.frag.c*this.iso.expectation/sum(this.frag.c)
      #this.frag.c <- this.frag.c[this.frag.c!=0]
      frag.atom.matrix[this.frag.group,names(this.frag.c)] <- this.frag.c
    }



  }

  ### stat
  {
    frag.certainty <- apply(frag.atom.matrix,1,function(x){
      sum(x==1)/sum(x)
    })
    frag.certainty[rowSums(frag.atom.matrix)==0] <- 0
    fg.df$certainty <- frag.certainty
  }
  cfm_data@fragment_group <- fg.df
  cfm_data@fragment_group_map <- frag.atom.matrix


  return(cfm_data)


}


cfm_data_add_seed <- function(cfm_data, smiles ){

  frag.str <- stringr::str_sub(cfm_data@fragment_define$fragment_id[1],1,-2)

  fragment_define <- cfm_data@fragment_define%>%
    bind_rows(data.frame(
      row.names = paste0(frag.str,0),
      fragment_id = paste0(frag.str,0),
      fragment_mz = 0,
      smiles = smiles,
      fragment_group =NA,
      formula = get_smile_formula(smiles)
    ),.)

  fragment_transition <- cfm_data@fragment_transition%>%
    bind_rows(data.frame(
      from = paste0(frag.str,0),
      to = paste0(frag.str,1)
    ),.)
  fragment_define -> cfm_data@fragment_define
  fragment_transition -> cfm_data@fragment_transition

  return(cfm_data)
}


get_cfm_data_fragment_group_atom_map <- function(cfm_data,frag.group){

  frag.idx <- which(cfm_data@fragment_define$fragment_group == frag.group)

  if (1%in%frag.idx) {
    ele <- get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfm_data))
    frag.atoms.prob <- rep(1,length(ele))
    names(frag.atoms.prob) <- ele
  }else{
    frag.def <- cfm_data@fragment_define[frag.idx,]
    #frag.score <- frag.def$ratio*0.7 + frag.def$bond.score*0.3
    #frag.idx<- frag.idx[frag.score == max(frag.score,na.rm = T)]
    frag.maps <- cfm_data@fragment_atom_map[frag.idx]
    frag.maps <- frag.maps[!sapply(frag.maps,is.null)]
    if (!length(frag.maps)) {
      ele <- get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfm_data))
      frag.atoms.prob <- rep(1,length(ele))
      names(frag.atoms.prob) <- ele
    }else{
      frag.atoms.prob <- sapply(frag.maps,rowSums)%>%
        rowMeans()
    }

  }


  return(frag.atoms.prob)
}


get_cfm_data_sdf_igraph <- function(cfm_data,fragment_id = 1 ){

  cfm_data@fragment_igraph[[fragment_id]]

}




heatmap_atom_iso_prob <- function(x){

  ComplexHeatmap::Heatmap(x,
                          na_col  ="#999999",
                          name = "isotope labeled\nprobability",
                          col = circlize::colorRamp2(breaks = c(0,0.5,1),
                                                     c("white","#F7844F","#B20C26")),
                          cluster_columns = F,
                          row_names_side  = "left",
                          rect_gp =  grid::gpar(lwd=2,col = "white"),
                          cluster_rows = F)

}


CFM_spectra_data_int_weight <- function(sp.data,iso_count){



  ### weighted intensity matrix
  {
    fg.idx <- split(1:nrow(sp.data),sp.data$fragment_group)
    if (!length(fg.idx))  return(sp.data)
    frag.iso.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso_count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso_count)))
    frag.int.sum <- c()
    for (i.fg in seq_along(fg.idx)) {
      x.df <- sp.data[fg.idx[[i.fg]],]
      x.int <- x.df%>%
        dplyr::select(-mz,-collisionEnergy)%>%
        tidyr::pivot_wider(names_from ="iso",
                           id_cols = "sp.id",
                           values_from = "intensity",
                           values_fn = sum)%>%
        tibble::column_to_rownames("sp.id")%>%
        dplyr::select(dplyr::starts_with("M"))%>%
        as.matrix()
      to.add <- setdiff(paste0("M",0:iso_count),colnames(x.int))
      x.int <- cbind(matrix(0,nrow(x.int),length(to.add),
                            dimnames = list(NULL,to.add)),x.int)
      x.int <- x.int[,paste0("M",0:iso_count),drop =F]
      x.int[is.na(x.int)] <- 0
      if (ncol(x.int)==1) {
        x.ratio <- x.int
        x.ratio[,1] <- 1
      }else{
        x.ratio <- t(apply(x.int,1,function(z) z/sum(z)))
      }
      x.weight <- rowSums(x.int)
      x.int.weighted <- apply(x.int,2,weighted.mean,w = x.weight)
      x.ratio.weighted <- apply(x.ratio,2,weighted.mean,w = x.weight)
      frag.iso.matrix[i.fg,] <- x.int.weighted
      frag.int.sum[i.fg] <- sum(x.int.weighted)
    }
    names(frag.int.sum) <- names(fg.idx)

  }


  ### merge
  {
    sp.mz <- sp.data%>%
      dplyr::filter(!is.na(fragment_group))%>%
      dplyr::mutate(fg.iso = paste0(fragment_group,iso))

    iso.peak.data <- frag.iso.matrix%>%
      as.data.frame()%>%
      rownames_to_column("fragment_group")%>%
      tidyr::pivot_longer(starts_with("M"),
                          names_to = "iso",
                          values_to = "intensity")%>%
      dplyr::mutate(sp.id= "combined_sp",
                    fg.iso = paste0(fragment_group,iso),
                    mz = sp.mz$mz[match(fg.iso,sp.mz$fg.iso)])%>%
      dplyr::filter(!is.na(mz))
    sp.data.weighted <- bind_rows(
      sp.data,
      iso.peak.data
    )%>%
      dplyr::select(-fg.iso)


  }
  return(sp.data.weighted)
}

CFM_spectra_data_merge <- function(sp.data,iso_count){


  sp.data$merged <- F
  ### weighted intensity matrix
  {
    fg.idx <- split(1:nrow(sp.data),sp.data$fragment_group)
    if (!length(fg.idx))  return(sp.data)
    frag.ratio.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso_count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso_count)))
    frag.df <- data.frame(
      fragment_group = names(fg.idx),
      int_sum = NA,
      peaks_count = NA,
      icc = NA,
      cos = NA
    )
    for (i.fg in seq_along(fg.idx)) {

      x.df <- sp.data[fg.idx[[i.fg]],]
      x.ratio <- x.df%>%
        tidyr::pivot_wider(names_from ="iso_count",
                           id_cols = "sp.id",
                           values_from = "ratio",
                           values_fn = mean
                           )%>%
        tibble::column_to_rownames("sp.id")%>%
        dplyr::select(dplyr::starts_with("M"))%>%
        as.matrix()
      x.ratio <- get_matrix_value_fill_with_NA(
        x.ratio,
        rownames_vec = row.names(x.ratio),
        colnames_vec = paste0("M",0:iso_count),drop = F)
      x.ratio[is.na(x.ratio)] <- 0
      if (ncol(x.ratio)==1) {
        x.ratio[,1] <- 1
      }
      x.int.sum <- x.df%>%
        dplyr::distinct(sp.id,.keep_all = T)%>%
        dplyr::pull(int_sum,name = sp.id)
      x.int.sum <- x.int.sum[rownames(x.ratio)]
      x.ratio.weighted <- apply(x.ratio,2,weighted.mean,w = x.int.sum)
      x.int.sum.weighted <- weighted.mean(x.int.sum,x.int.sum)
      x.icc <- irr::icc(t(x.ratio), model = "twoway",
          type = "consistency", unit = "single")$value
      x.cos <- lsa::cosine(x.ratio.weighted,t(x.ratio))
      x.cos.weight <- weighted.mean(x.cos,w = log10(x.int.sum))
      frag.ratio.matrix[i.fg,] <- x.ratio.weighted
      frag.df$int_sum[i.fg] <- x.int.sum.weighted
      frag.df$icc[i.fg] <- x.icc
      frag.df$cos[i.fg] <- x.cos.weight
      frag.df$peaks_count[i.fg] <- nrow(x.ratio)

      #Heatmap(x.ratio,cluster_columns = F,cluster_rows = F)

    }

  }


  ### merge
  {
    sp.mz <- sp.data%>%
      dplyr::filter(!is.na(fragment_group))%>%
      dplyr::mutate(fg.iso = paste0(fragment_group,iso_count))

    iso.peak.data <- frag.ratio.matrix%>%
      as.data.frame()%>%
      rownames_to_column("fragment_group")%>%
      tidyr::pivot_longer(starts_with("M"),
                          names_to = "iso_count",
                          values_to = "ratio")%>%
      dplyr::mutate(sp.id= "combined_sp",
                    fg.iso = paste0(fragment_group,iso_count),
                    mz = sp.mz$mz[match(fg.iso,sp.mz$fg.iso)],
                    frag.df[match(fragment_group,frag.df$fragment_group),],
                    intensity = int_sum*ratio ,
                    merged = T )%>%
      dplyr::filter(!is.na(mz))

    sp.data$merged <- F
    sp.data.weighted <- bind_rows(
        sp.data,
      iso.peak.data )%>%
      dplyr::select(-fg.iso)


  }

  return(sp.data.weighted)
}

CFM_spectra_data_remove_natural <-function(sp.data,
                                           natural.ratio,
                                           if.map){

  ###
  {
    if (sum(sp.data$sp.id=="combined_sp")==0) return(sp.data)
  }


  ### calculate intensity from natrual
  {
    sp.data.merged <- sp.data%>%
      dplyr::filter(sp.id== "combined_sp")%>%
      dplyr::mutate(fg.iso = paste0(fragment_group,"_",iso))
    int_sum <- sp.data.merged%>%
      dplyr::group_by(fragment_group)%>%
      dplyr::summarise(int_sum = sum(intensity))%>%
      #dplyr::mutate(int_sum = sum(intensity),
      #              fgi = paste0(fragment_group,"_",iso))%>%
      #column_to_rownames("fragment_group")%>%
      dplyr::pull(int_sum,name  = fragment_group)
    max.iso <- max(str_extract_num(sp.data.merged$iso))
    suffix <- paste0("_M",rep(0:max.iso,length(int_sum)))
    int_sum <- rep(int_sum,each = max.iso+1)
    names(int_sum) <- paste0(names(int_sum) , suffix)
    natural.distribution <- apply(if.map@isoform.map,1,
                                sum)/ncol(if.map@isoform.map)
    natural.int <- int_sum[names(natural.distribution)]*natural.distribution*natural.ratio

  }

  {
    sp.data.removed <-  sp.data.merged %>%
      dplyr::mutate(intensity = intensity - natural.int[fg.iso],
                    intensity = case_when(intensity <0 ~0,
                                          is.na(intensity)~0,
                                          T~intensity))%>%
      dplyr::group_by(fragment_group)%>%
      dplyr::filter(!all(intensity==0))%>%
      dplyr::select(-fg.iso)
    sp.data <- sp.data%>%
      dplyr::filter(!sp.id== "combined_sp")%>%
      rbind(sp.data.removed)
  }

  return(sp.data)


}


get_CFM_data_MSIPFragmentMap<- function(cfmd){

  cfmd.sp <- get_CFM_data_Spectra(cfmd)
  cfmd.msip.core <- get_MSIPCoreData(cfmd.sp,cfmd,0)
  cfmd.fg.map <- cfmd.msip.core@FG_map

  cfmd.fg.map
}



get_CFM_data_from_smiles <- function(smiles = "NCC(O)=O",
                                     compound_id = "temp_id",
                                     ppm = 5,
                                     adduct = "[M+H]+",
                                     check_temp = T,
                                     iso_ele = "[13]C",
                                     temp_dir = tempdir(),
                                     ...){
  if(check_temp){
    if (!dir.exists(temp_dir)) dir.create(temp_dir,recursive = T,showWarnings = F)
    temp_file <- paste0(temp_dir,"/",compound_id,"_",adduct,".rds")
    if(file.exists(temp_file)){
      message_with_time("loading from temp:",temp_file)
      cfmd <- readRDS(temp_file)
      return(cfmd)
    }
  }
  log.info <- c()
  log.info["smiles"] <- smiles
  message_with_time("CFM_annotate_by_predict")
  start.time <- Sys.time()
  cfmd <- CFM_annotate_by_predict(smiles_or_inchi = smiles,
                                  id = compound_id,
                                  ppm_mass_tol = ppm,
                                  abs_mass_tol = 0.005,
                                  param_adduct = adduct )
  cfmd.temp.file <- paste0(temp_dir,"/cfmd.temp.",compound_id,".",adduct,".rds")
  saveRDS(cfmd,cfmd.temp.file)
  log.info["cfm.time"] <- (Sys.time()-start.time)%>%
    as.numeric(units = "mins")
  start.time <- Sys.time()
  cfmd <- cfm_data_add_seed(cfmd,smiles)
  cfmd <- CFM_data_get_igraph(cfmd)

  message_with_time("CFM_data_get_atom_map")
  cfmd <- CFM_data_get_atom_map(cfmd,iso_ele = iso_ele)
  cfmd <- cfm_data_get_FG_map(cfmd,iso_ele = iso_ele)
  message_with_time("Done")
  log.info["map.time"] <- (Sys.time()-start.time)%>%as.numeric(units = "mins")
  file.remove(cfmd.temp.file)

  if(check_temp){
    dir.create(temp_dir,showWarnings = F,recursive = T)
    temp_file <- paste0(temp_dir,"/",compound_id,"_",adduct,".rds")
    saveRDS(cfmd,file = temp_file)

    ### log
    log.info["atom.count"] <- length(get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfmd)))
    cat(paste0(paste0(log.info,collapse = ","),"\n"),
        file = paste0(temp_dir,"/atm.log"),append = T)


  }
  return(cfmd)




}




shiny_vis_cfmd_FG_map <- function(cfmd){


  .ui <- function(){
    fluidPage(
      column(width = 6,
             shiny::plotOutput(outputId = "heatmap_atom_map",height  = "800px")),
      column(width = 6,
             selectInput(inputId = "fg_id",choices = sort(unique(cfmd@fragment_define$fragment_group)),
                         label = "Fragment group",selected = cfmd@fragment_define$fragment_group[1]),
             selectInput(inputId = "fragment",label = "fragment",
                         choices = sort(unique(cfmd@fragment_define$fragment_id)),
                         selected = cfmd@fragment_define$fragment_id[1] ),
             visNetwork::visNetworkOutput(outputId = "atom_map",height  = "800px"))
    )

  }

  .server <- function(){
    function(input, output, session) {

      output$heatmap_atom_map <- renderPlot({

        message_with_time("heatmap_atom_map")
        get_CFM_data_MSIPFragmentMap(cfmd)%>%
          heatmap_MSIPFragmentMap()
      })

      output$atom_map <- visNetwork::renderVisNetwork({
        message_with_time("atom_map")
        vis_cfm_data_fragment_atom_map(cfmd ,input$fragment,show_id = F)
      })

      observeEvent(input$fg_id,{

        message_with_time("fg_id")

        x <- cfmd@fragment_define%>%
          dplyr::filter(fragment_group == input$fg_id)
        updateSelectInput(inputId = "fragment",choices = x$fragment_id)
      })

    }
  }

  ### Start Shiny APP
  {
    shinyApp(ui = .ui(),
             server = .server(),
             options = list(host = "0.0.0.0",
                            #port = 6548,
                            launch.browser = T))
  }


}


shiny_vis_cfmd_trans <- function(cfmd){


  .ui <- function(){
    fluidPage(
      column(width = 6,
             visNetwork::visNetworkOutput(outputId = "vis_trans",height = "800px"),
             style = "border: 1px solid #aaa; padding: 6px;"),
      column(width = 6,
             plotlyOutput(outputId = "cfm_sp"),
             verbatimTextOutput("frag_info"),
             visNetwork::visNetworkOutput(outputId = "atom_map" ))
    )

  }

  .server <- function(){
    function(input, output, session) {

      output$vis_trans <-  visNetwork::renderVisNetwork({

        message_with_time("get_CFM_data_trans_igraph")
        vis_igraph(get_CFM_data_trans_igraph(cfmd)) %>%
          visOptions(nodesIdSelection =
                       list(enabled  = T,
                            selected  = cfmd@fragment_define$fragment_id[1]))


      })


      output$frag_info <- renderText({

        message_with_time("frag_info")
        x <- input$vis_trans_selected
        if(is.null(x)) return(NULL)
        if(x==""){      return(NULL)  }
        y <- cfmd@fragment_define
        paste0("Fragment: ",x,"\n",
               "Formula: ",y[x,"formula"],"\n",
               "mz: ",y[x,"fragment_mz"],"\n"  )

      })
      output$atom_map <-  visNetwork::renderVisNetwork({

        message_with_time("get_CFM_data_trans_igraph")
        x <- input$vis_trans_selected
        if(is.null(x)) return(NULL)
        if(x==""){      return(NULL)  }

        vis_cfm_data_fragment_atom_map(cfmd,input$vis_trans_selected,show_id = F)


      })


      output$cfm_sp <- renderPlotly({

        plotly_Spectra(get_CFM_data_Spectra(cfmd))

      })


    }
  }

  ### Start Shiny APP
  {
    shinyApp(ui = .ui(),
             server = .server(),
             options = list(host = "0.0.0.0",
                            #port = 6548,
                            launch.browser = T))
  }


}




get_cfmd_FG_map_check <- function(cfmd){

  fg.map <- get_CFM_data_MSIPFragmentMap(cfmd)
  fg.fragment.count <- table(cfmd@fragment_define$fragment_group)
  fg.certainty <- get_MSIPFragmentMap_certainty(fg.map)
  fg.df <- data.frame(
    fg = names(fg.certainty),
    certainty =fg.certainty,
    fragment.count = as.numeric(fg.fragment.count[names(fg.certainty)])
  )

  return(fg.df)

}


get_cfmd_map_evaluation <- function(cfmd){



}

