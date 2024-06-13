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
#' @description
#' see CFM doc
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
#' @import crayon stringr

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
    cfm.result <- read_CFM_predict_result(out.file)
    return(invisible(cfm.result))
  }else{
    message("Result export to: \n",crayon::red(output_file_or_dir))
    return(invisible(out.file))
  }
}

#' CFM_annotate
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
#' @import crayon stringr magrittr
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


#' CFM_annotate_by_fraggen
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
  cfm.anno <- CFM_annotate(smiles_or_inchi = smiles_or_inchi,
               spectrum_file = cfm.pred.sp,
               ppm_mass_tol=ppm_mass_tol,
               abs_mass_tol=abs_mass_tol,
               param_adduct=param_adduct   )
  return(cfm.anno)


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
      dplyr::select(from,to,smiles)

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
#'
#' @param result_path path
#'
#' @return list of cfm data
#' @import dplyr stringr
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
    plot_Spectra()

}



#' CFM_annotate_isotopologues
#'
#' @param sp Spectra
#' @param cfmd cfmd_data
#' @param isotope "[13]C"
#' @param iso.count num
#' @param ppm 20
#'
#' @return null
#' @export
#' @import MSCC
#'
CFM_annotate_isotopologues <- function(sp,
                                 cfmd,
                                 isotope = "[13]C",
                                 iso.count = 0 ,
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


  iso.mz.diff <- (0:iso.count)*chemform_mz("[13]CC-1")
  mz.labeled.m <- matrixSub(cfm.peaks.data$mz,-iso.mz.diff)%>%
    `colnames<-`(paste0("M",0:iso.count))%>%
    as.data.frame()%>%
    dplyr::mutate(fragment_group=cfm.peaks.data$fragment_group)%>%
    dplyr::distinct(M0,.keep_all = T)%>%
    tidyr::pivot_longer(paste0("M",0:iso.count),names_to = "iso",values_to = "mz")



  sp.data <- get_Spectra_data(sp,var = "collisionEnergy")%>%
    dplyr::mutate(fragment_group = NA,iso=NA)
  if (nrow(sp.data)) {


  sp.data <- sp.data%>%
    dplyr::mutate(idx =match_mz(mz1=mz,
                                mz2 = mz.labeled.m$mz,
                                mz.ppm = ppm),
                  fragment_group = mz.labeled.m$fragment_group[idx],
                   iso = mz.labeled.m$iso[idx])%>%
    dplyr::select(-idx)
  }


  return(sp.data)
}





setClass("CFM_data",
         slots = list(
           peak_assignment = "data.frame",
           fragment_define = "data.frame",
           fragment_transition = "data.frame",
           fragment_igraph = "list",
           fragment_sdf = "SDFset",
           fragment_atom_map = "list"
         ))


setMethod("show",signature ="CFM_data",definition =
            function(object){
              message("CFM_data with ",nrow(object@fragment_define)," fragment")
            } )


#' CFM_data_get_igraph
#' get cfm igraph and atom map
#'
#' @param object cfmd
#'
#' @return cfmd
#' @export
#' @import magrittr
CFM_data_get_igraph <- function(object){

  ### fragment def
  {
    fragment.data <- object@fragment_define
    fragment.sdf <- suppressWarnings(smiles2sdf(fragment.data$smiles))
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
    fragment.sdf <- suppressWarnings(smiles2sdf(fragment.data$smiles))
    cid(fragment.sdf) <- fragment.data$fragment_id
    fragment.igraph <- get_sdf_igraph(fragment.sdf)
    names(fragment.igraph) <- fragment.data$fragment_id
  }


  ### atom map
  {
    fragment.trans <- object@fragment_transition
    ### map for trans
    trans.atom.map <- as.list(rep(NA,nrow(fragment.trans)))
    fragment.atom.map <-list()
    .trace.atom.map <- function(i.trans){
      ig.parent <- fragment.igraph[[fragment.trans$from[i.trans]]]
      ig.product <- fragment.igraph[[fragment.trans$to[i.trans]]]
      sdf.parent <- fragment.sdf[[fragment.trans$from[i.trans]]]
      sdf.product <- fragment.sdf[[fragment.trans$to[i.trans]]]
      atom.map <- get_atom_map(sdf.parent ,sdf.product ,ig.parent ,ig.product )
      atom.map
    }
    if (nrow(fragment.trans)) {


      ### map for frag
      ig.trans <- get_CFM_data_trans_igraph(object)
      for (i in 1:nrow(fragment.data)) {
        this.frag <- fragment.data$fragment_id[i]
        if(i==1){
          ele <- get_atom_from_igraph(fragment.igraph[[1]])
          maps <- diag(nrow = length(ele))
          rownames(maps)<-colnames(maps)<-ele
          fragment.atom.map[[i]] <- maps
          next
        }
        this.path <- shortest_paths(ig.trans,
                                    1,this.frag,
                                    output = "epath")
        ep <- this.path$epath[[1]]
        if (!length(ep)==0) {
          trans.idx <- match(ep,E(ig.trans))
          maps <- lapply(trans.idx,function(x){
            if (all(is.na(trans.atom.map[[x]]))) {
              return( .trace.atom.map(x))
            }
            return(trans.atom.map[[x]])
          })
          trans.atom.map[trans.idx] <- maps

          while(length(maps)>1){
            maps[[2]] <- maps[[1]]%*% maps[[2]]
            maps[[1]] <-NULL
          }
          fragment.atom.map[[i]] <- maps[[1]]
        }

      }
      names(fragment.atom.map) <-fragment.data$fragment_id
    }


  }


  ### save to object
  {

    fragment.data -> object@fragment_define
    object@fragment_igraph <- fragment.igraph
    object@fragment_sdf <- fragment.sdf
    object@fragment_atom_map <- fragment.atom.map
  }

  return(object)
}



get_CFM_data_trans_igraph <- function(object){


  frag.trans.graph <- graph_from_data_frame(object@fragment_transition,
                                            vertices =object@fragment_define )
  return(frag.trans.graph)
}



#' Title
#'
#' @param cfm_data cfmd
#' @param ppm 10
#' @import magrittr
#' @return cfmd
cfm_data_get_fragment_group <- function(cfm_data,ppm = 10){

  fg <- groupMz(cfm_data@fragment_define$fragment_mz,ppm)
  cfm_data@fragment_define$fragment_group <- paste0("FG",num2str(fg))
  cfm_data@peak_assignment$fragment_group <-
    cfm_data@fragment_define$fragment_group[match(
      cfm_data@peak_assignment$fragment_id,
      cfm_data@fragment_define$fragment_id)]

  return(cfm_data)
}


get_cfm_data_fg_atom_map <- function(cfm_data,frag.group){

  frag.idx <- which(cfm_data@fragment_define$fragment_group == frag.group)
  if (1%in%frag.idx) {
    ele <- get_atom_from_igraph(get_cfm_data_sdf_igraph(cfm_data))
    frag.atoms.prob <- rep(1,length(ele))
    names(frag.atoms.prob) <- ele
  }else{
    frag.def <- cfm_data@fragment_define%>%
      dplyr::filter(fragment_group == frag.group)
    frag.maps <- cfm_data@fragment_atom_map[frag.idx]
    frag.maps <- frag.maps[!sapply(frag.maps,is.null)]
    if (!length(frag.maps)) {
      ele <- get_atom_from_igraph(get_cfm_data_sdf_igraph(cfm_data))
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


CFM_spectra_data_int_weight <- function(sp.data,iso.count){



  ### weighted intensity matrix
  {
    fg.idx <- split(1:nrow(sp.data),sp.data$fragment_group)
    if (!length(fg.idx))  return(sp.data)
    frag.iso.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso.count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso.count)))
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
      to.add <- setdiff(paste0("M",0:iso.count),colnames(x.int))
      x.int <- cbind(matrix(0,nrow(x.int),length(to.add),
                            dimnames = list(NULL,to.add)),x.int)
      x.int <- x.int[,paste0("M",0:iso.count),drop =F]
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
    natural.distribution <- apply(if.map$iso.form.map,1,
                                sum)/ncol(if.map$iso.form.map)
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
