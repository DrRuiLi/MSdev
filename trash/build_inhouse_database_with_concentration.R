#' @param compound.record.file
#'
#' @title build inhouse database
#'
#' @return
#' @export
#' @import xcms
#' @examples
build_inhouse_database_with_concentraion <-
  function(compound.record.file ) {

    ### test
    compound.record.file<- "C:\\Users\\91879\\OneDrive\\Documents\\Code\\R\\Projecct\\2022.7.22.MS.inhouse.database.building\\Standard.record.2021.12.18.STD_01.xlsx"
    compound.record.file<- "d:/2022.6.27.STD/Internal.standard.2022.6.27.xlsx"

    ###check files
    {
      if (!file.exists(compound.record.file)) {
        stop.info <- paste0(compound.record.file , " not exist ") %>%
          crayon::red()
        stop(stop.info)
      }


      compound.record <- Retrieve_compounds_data_from_pubchem(compound.record.file)
      compound.sample.info <- readxl::read_excel(path = compound.record.file ,
                                                 sheet = "Concentration_info")

    }


    ###pos
    {

      xcms.pos <- xcms_processing(compound.sample.info,"positive")
      plot_xcms_features_distribution(xcms.pos)
      xcms.pos.features <-featureDefinitions(xcms.pos)
      MS.network.pos <- expand_adduct_from_compounds(compound.record,"positive")
      MS.network.pos <- match_adduct_to_features(MS.network.pos,xcms.pos,ppm.thresh = 10)
      plot_adduct_distribution(MS.network.pos,3)


    }
    ###neg
    {

      xcms.neg <- xcms_processing(compound.sample.info,"negative")
      plot_xcms_features_distribution(xcms.neg)
      xcms.neg.features <-featureDefinitions(xcms.neg)
      MS.network.neg <- expand_adduct_from_compounds(compound.record,"negative")
      MS.network.neg <- match_adduct_to_features(MS.network.neg,xcms.neg,ppm.thresh = 10)
      plot_adduct_distribution(MS.network.neg,3)


    }




  }
