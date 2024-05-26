#' @title build inhouse database
#' @description Detect ion peaks from compounds standard. This method just read in single injection (pos and neg) and detect possible adduct form mz.
#' Usually can not determine whether peaks generated from compounds.
#' @param compound.record.file paht of xlsx file, contain column pubchem.cid
#'
#' @return
#' @export
#' @import xcms
#' @examples
build_inhouse_database_single_injection <-
  function(compound.record.file ) {

    ### test
    compound.record.file<- "C:\\Users\\91879\\OneDrive\\Documents\\Code\\R\\Projecct\\2022.7.22.MS.inhouse.database.building\\Standard.record.2021.12.18.STD_01.xlsx"
    compound.record.file<- "d:/2022.6.27.STD/Internal.standard.2022.6.27.xlsx"
    compound.record.file<- "C:\\Users\\91879\\OneDrive\\Documents\\Code\\R\\Default_Workspace\\Internal.standard.2022.3.7.xlsx"
    ###check files
    {
      if (!file.exists(compound.record.file)) {
        stop.info <- paste0(compound.record.file , " not exist ") %>%
          crayon::red()
        stop(stop.info)
      }


      compound.record <- Retrieve_compounds_data_from_pubchem(compound.record.file)
      if (is.null(compound.record$mzML.positive)) {
        message("Select positive mzML file")
        compound.record$mzML.positive <-
          choose.files(default = getwd(), caption = "Select positive")
        openxlsx::write.xlsx(compound.record , file = compound.record.file)
      }
      if (is.null(compound.record$mzML.negative)) {
        message("Select negative mzML file")
        compound.record$mzML.negative <-
          choose.files(default = getwd(), caption = "Select negative")
        openxlsx::write.xlsx(compound.record , file = compound.record.file)
      }
    }


    ###pos
    {

      xcms.pos <- readMSData(unique(compound.record$mzML.positive),
                             mode = "onDisk")
      centwave.param <- CentWaveParam(peakwidth = c(5,30),
                                      prefilter = c(3,100),
                                      snthresh = 10,
                                      ppm = 10)
      xcms.pos<-findChromPeaks(xcms.pos,
                               param = centwave.param,
                               BPPARAM = SerialParam())
      xcms.pos <- groupChromPeaks(xcms.pos , param = PeakDensityParam(sampleGroups = "sample",
                                                                        binSize = 0.015))

      plot_xcms_peaks_distribution(xcms.pos)
      plot_xcms_features_distribution(xcms.pos)
      xcms.pos.peaks <-chromPeaks(xcms.pos)
      MS.network.pos <- expand_adduct_from_compounds(compound.record,"positive")
      MS.network.pos <- match_adduct_to_features(MS.network.pos,
                                                 xcms.pos,
                                                 ppm.thresh = 10,
                                                 rt.tol = 10)
      MS.network.pos <-  match_adduct_with_eicSimilarity(MS.network.pos,xcms.pos)
      plot_adduct_distribution(MS.network.pos,3,rt.filter = T)
      plot_adduct_chromatogram(MS.network.pos,3,
                               rt.filter = T,cor.thresh = 0.3,
                               norm =F,
                               move = T)



    }
    ###neg
    {

      xcms.neg <- readMSData(unique(compound.record$mzML.negative),
                             mode = "onDisk")
      centwave.param <- CentWaveParam(peakwidth = c(5,30),
                                      prefilter = c(3,100),
                                      snthresh = 10,
                                      ppm = 10)
      xcms.neg<-findChromPeaks(xcms.neg,
                               param = centwave.param,
                               BPPARAM = SerialParam())
      xcms.neg <- groupChromPeaks(xcms.neg , param = PeakDensityParam(sampleGroups = "sample",
                                                                      binSize = 0.015))

      plot_xcms_peaks_distribution(xcms.neg)
      plot_xcms_features_distribution(xcms.neg)
      xcms.neg.peaks <-chromPeaks(xcms.neg)
      MS.network.neg <- expand_adduct_from_compounds(compound.record,"negative")
      MS.network.neg <- match_adduct_to_features(MS.network.neg,
                                                 xcms.neg,
                                                 ppm.thresh = 10,
                                                 rt.tol = 10)
      MS.network.neg <-  match_adduct_with_eicSimilarity(MS.network.neg,xcms.neg)
      plot_adduct_distribution(MS.network.neg,3,rt.filter = F)
      plot_adduct_chromatogram(MS.network.neg,3,
                               rt.filter = F,cor.thresh = -1,
                               norm = T,
                               move = T)
    }




  }

