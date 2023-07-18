#' @param compound.record.file
#'
#' @title build inhouse database
#'
#' @return
#' @export
#' @import xcms
#' @examples
build_inhouse_database <-
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
                                      snthresh = 100,
                                      ppm = 10)
      xcms.pos<-findChromPeaks(xcms.pos,
                                param = centwave.param,
                               BPPARAM = SerialParam())
      plot_xcms_peaks_distribution(xcms.pos)
      xcms.pos.peaks <-chromPeaks(xcms.pos)
      MS.network.pos <- expand_adduct_from_compounds(compound.record,"positive")
      MS.network.pos <- match_adduct_to_peaks(MS.network.pos,xcms.pos.peaks,ppm.thresh = 10)
      plot_adduct_distribution(MS.network.pos,5)


    }
    ###neg
    {

      xcms.neg <- readMSData(unique(compound.record$mzML.negative),
                             mode = "onDisk")
      centwave.param <- CentWaveParam(peakwidth = c(5,30),
                                      prefilter = c(3,100),
                                      snthresh = 100,
                                      ppm = 10)
      xcms.neg<-findChromPeaks(xcms.neg,
                               param = centwave.param,
                               BPPARAM = SerialParam())
      plot_xcms_peaks_distribution(xcms.neg)
      xcms.neg.peaks <-chromPeaks(xcms.neg)
      MS.network.neg <- expand_adduct_from_compounds(compound.record,"negative")
      MS.network.neg <- match_adduct_to_peaks(MS.network.neg,xcms.neg.peaks,ppm.thresh = 10)
      plot_adduct_distribution(MS.network.neg,5)


    }




  }
