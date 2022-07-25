ms_convert <- function(ms.ana){

  message(Sys.time()," MSconverting...")

  ### check
  {
    if (isTRUE(ms.ana$processing.info$raw.data.convert$done) ) {
      #message("Raw data convert done")
      return(ms.ana)
    }

  }

  ### convert
  {
    msconvert_wiff2mzML(wiff.files = ms.ana$sample.info$raw.file.positive,
                        mzML.files = ms.ana$sample.info$mzML.file.positive,
                        BPPARAM = SnowParam(workers = thread - 1))

    msconvert_wiff2mzML(wiff.files = ms.ana$sample.info$raw.file.negative,
                        mzML.files = ms.ana$sample.info$mzML.file.negative,
                        BPPARAM = SnowParam(workers = thread - 1))

  }

  ### return
  {
    ms.ana$processing.info$raw.data.convert <- list(
      done = T,
      time = Sys.time()

    )
    save(ms.ana , file = ms.ana$processing.info$project.info$ms.ana.file)
    return(ms.ana)

  }

}
