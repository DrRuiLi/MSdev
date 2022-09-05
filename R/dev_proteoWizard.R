msconvert_raw2mzML <- function(raw.files ,
                                mzML.files,
                                BPPARAM = SnowParam(workers = 10)){


  raw.files <- gsub(pattern = "\\",x = raw.files,replacement = "/",fixed = T)%>%
    na.omit()
  mzML.files <- gsub(pattern = "\\",x = mzML.files,replacement = "/",fixed = T)%>%
    na.omit()

  ###check msconvert
  {
    msconvert_return <- try(system("msconvert",
                                   intern = T))
    if(!any(grepl(pattern = "Usage: msconvert", x = msconvert_return))){
      stop("Command msconvert error, please check environment variables")
    }

    }
  ###check file and directory
  {
    if(!any(file.exists(raw.files))){
      stop(paste0("File not found : ",sum(!file.exists(raw.files)),"/", length(raw.files)))
    }
    if(length(raw.files) != length(mzML.files)){
      stop("raw files and mzml files not match")
    }
    sapply(unique(dirname(mzML.files)),dir.create,recursive =T,showWarnings =F)

  }

  ###msconvert
  {

    shell.commomd <- paste0("msconvert  --filter \"peakPicking true 1-\" --mzML ",
                            raw.files,
                            " -o ",
                            dirname(mzML.files),
                            " --outfile ",
                            mzML.files)
    #system(shell.commomd,intern = T)


    bplapply(shell.commomd,
             FUN = function(x){ system(x,intern = T)},
             BPPARAM = BPPARAM)
    return(0)

  }


}
